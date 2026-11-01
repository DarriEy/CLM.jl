# FUN-path crop N parity: soybean N fixation and the FUN uptake / allocation fluxes
# of a coupled ctsm5.3.012 crop column with use_fun = .true.
#
# Reference: the instrumented crop_fun_macos case
# (scripts/validation/crop_parity/build_crop_fun_case.sh), cold-started at the
# USplains crop point with use_cn + use_crop + use_fertilizer + use_flexiblecn +
# use_fun (docs/CROP_PARITY.md § "FUN crop N"). Per dumped step it writes
#   fundump_before_cnfun_n<n>.nc / fundump_after_cnfun_n<n>.nc   (every CNFUN input
#       and output, around the CNFUN call inside SoilBiogeochemCompetition)
#   cropdump_{after_soyfix,after_cropalloc,after_cnalloc,end_noleach}_n<n>.nc
#
# Every check INJECTS the Fortran inputs of step n into the Julia routine and
# compares with the Fortran outputs of the same step (the #338 method):
#   (1) cnfun!            before_cnfun(n)  -> after_cnfun(n)   (+ the p2c to column)
#   (2) FlexibleCN N demand under FUN (calc_plant_nutrient_demand_flexiblecn!, the
#       prognostic-crop call)   after_cropalloc(n) -> before_cnfun(n)
#   (3) FlexibleCN allocation under FUN (calc_plant_cn_alloc_flexiblecn!: sminn_to_npool,
#       plant_nalloc, plant_calloc = npp_growth, cpool_to_*, npool_to_*)
#                          after_cnfun(n) + after_cnalloc(n) state -> after_cnalloc(n)
#   (4) the harvest step of CNOffsetLitterfall with use_fun (FUN leaf-N split)
#   (5) SOYFIXN / FPG semantics under FUN.
#
#   julia --project=. scripts/fortran_parity_crop_fun.jl [REF]
#   REF default: $SYMFLUENCE_DATA/runs/crop_fun_ref
# Also `include`d by test/test_crop_fun_parity.jl (set CROP_FUN_NO_MAIN=1).
using CLM
using NCDatasets
using Printf

if !isdefined(@__MODULE__, :ct_load)
    ENV["CROP_TAILS_NO_MAIN"] = "1"
    include(joinpath(@__DIR__, "fortran_parity_crop_tails.jl"))
end

const CF_REF = get(ENV, "CROP_FUN_REF", joinpath(_SD, "runs", "crop_fun_ref"))
const CF_DT = 3600.0
const CF_SOY = (23, 24)   # temperate soybean, rainfed / irrigated (0-based itype)

function cf_load(dir, label, n)
    f = joinpath(dir, "fundump_$(label)_n$(n).nc")
    isfile(f) || return nothing
    d = Dict{String,Any}()
    NCDataset(f) do ds
        for k in keys(ds)
            v = ds[k]
            d[k] = ndims(v) == 0 ? v[] : Array(v)
        end
    end
    return d
end

function cf_steps(dir, label = "before_cnfun")
    rx = Regex("^fundump_$(label)_n(\\d+)\\.nc\$")
    s = Int[]
    for f in readdir(dir)
        m = match(rx, f)
        m === nothing || push!(s, parse(Int, m[1]))
    end
    return sort!(s)
end

# CLM5 "20SL_8.5m" decomposition-layer thicknesses (initVerticalMod.F90; same
# formula as CLM.varcon_init!), nlevdecomp = 20.
function cf_dzsoi_decomp(nlev)
    dz = zeros(nlev)
    for j in 1:nlev
        dz[j] = j <= 4 ? j * 0.02 : (j <= 13 ? 4 * 0.02 + (j - 4) * 0.04 : (4 * 0.02 + 9 * 0.04) + (j - 13) * 0.10)
    end
    return dz
end

# Accumulate one comparison. Besides CTStat's max abs / rel error, track the field's
# scale (max |Fortran|) so the report can give the error relative to the flux scale
# (a 1e-19 residual on a near-zero flux is a large relative but a negligible absolute
# error), and count SPVAL branch flips (one side 1e36, the other finite) separately.
const CF_SCALE = Dict{UInt,Float64}()
const CF_FLIPS = Dict{UInt,Int}()
function _cfacc!(s::CTStat, jl, f, tag; floor = 1e-12)
    if (abs(jl) >= 1e30) != (abs(f) >= 1e30)
        CF_FLIPS[objectid(s)] = get(CF_FLIPS, objectid(s), 0) + 1
        return
    end
    _acc!(s, jl, f, tag; floor)
    (isfinite(f) && abs(f) < 1e30) && (CF_SCALE[objectid(s)] = max(get(CF_SCALE, objectid(s), 0.0), abs(f)))
    return
end
cf_scale(s::CTStat) = get(CF_SCALE, objectid(s), 0.0)
cf_flips(s::CTStat) = get(CF_FLIPS, objectid(s), 0)
cf_normerr(s::CTStat) = (sc = cf_scale(s); sc > 0 ? s.maxabs / sc : s.maxabs)

# ---------------------------------------------------------------------------
# (1) cnfun!
# ---------------------------------------------------------------------------
const CF_FUN_OUT_P = ("NPP_NACTIVE", "NPP_NNONMYC", "NPP_NAM", "NPP_NECM",
    "NPP_NACTIVE_NO3", "NPP_NNONMYC_NO3", "NPP_NACTIVE_NH4", "NPP_NNONMYC_NH4",
    "NPP_NFIX", "NPP_NRETRANS", "NPP_NUPTAKE", "NPP_GROWTH", "NPP_BURNEDOFF", "SOILC_CHANGE",
    "NACTIVE", "NNONMYC", "NAM", "NECM", "NACTIVE_NO3", "NNONMYC_NO3", "NAM_NO3", "NECM_NO3",
    "NACTIVE_NH4", "NNONMYC_NH4", "NAM_NH4", "NECM_NH4", "NPASSIVE", "NFIX",
    "COST_NFIX", "COST_NACTIVE", "COST_NRETRANS", "NUPTAKE_NPP_FRACTION", "NRETRANS", "NUPTAKE",
    "RETRANSN_TO_NPOOL", "FREE_RETRANSN_TO_NPOOL", "SMINN_TO_PLANT_FUN", "PLANT_NDEMAND_RETRANS",
    "PLANTCN", "STORAGE_CDEMAND", "STORAGE_NDEMAND", "LEAFC_STORAGE_XFER_ACC", "LEAFN_STORAGE_XFER_ACC")

function _cf_jl_out(k, vs, cs, cf, ns, nf)
    k == "NPP_NACTIVE" && return cf.npp_Nactive_patch
    k == "NPP_NNONMYC" && return cf.npp_Nnonmyc_patch
    k == "NPP_NAM" && return cf.npp_Nam_patch
    k == "NPP_NECM" && return cf.npp_Necm_patch
    k == "NPP_NACTIVE_NO3" && return cf.npp_Nactive_no3_patch
    k == "NPP_NNONMYC_NO3" && return cf.npp_Nnonmyc_no3_patch
    k == "NPP_NACTIVE_NH4" && return cf.npp_Nactive_nh4_patch
    k == "NPP_NNONMYC_NH4" && return cf.npp_Nnonmyc_nh4_patch
    k == "NPP_NFIX" && return cf.npp_Nfix_patch
    k == "NPP_NRETRANS" && return cf.npp_Nretrans_patch
    k == "NPP_NUPTAKE" && return cf.npp_Nuptake_patch
    k == "NPP_GROWTH" && return cf.npp_growth_patch
    k == "NPP_BURNEDOFF" && return cf.npp_burnedoff_patch
    k == "SOILC_CHANGE" && return cf.soilc_change_patch
    k == "NACTIVE" && return nf.Nactive_patch
    k == "NNONMYC" && return nf.Nnonmyc_patch
    k == "NAM" && return nf.Nam_patch
    k == "NECM" && return nf.Necm_patch
    k == "NACTIVE_NO3" && return nf.Nactive_no3_patch
    k == "NNONMYC_NO3" && return nf.Nnonmyc_no3_patch
    k == "NAM_NO3" && return nf.Nam_no3_patch
    k == "NECM_NO3" && return nf.Necm_no3_patch
    k == "NACTIVE_NH4" && return nf.Nactive_nh4_patch
    k == "NNONMYC_NH4" && return nf.Nnonmyc_nh4_patch
    k == "NAM_NH4" && return nf.Nam_nh4_patch
    k == "NECM_NH4" && return nf.Necm_nh4_patch
    k == "NPASSIVE" && return nf.Npassive_patch
    k == "NFIX" && return nf.Nfix_patch
    k == "COST_NFIX" && return nf.cost_nfix_patch
    k == "COST_NACTIVE" && return nf.cost_nactive_patch
    k == "COST_NRETRANS" && return nf.cost_nretrans_patch
    k == "NUPTAKE_NPP_FRACTION" && return nf.nuptake_npp_fraction_patch
    k == "NRETRANS" && return nf.Nretrans_patch
    k == "NUPTAKE" && return nf.Nuptake_patch
    k == "RETRANSN_TO_NPOOL" && return nf.retransn_to_npool_patch
    k == "FREE_RETRANSN_TO_NPOOL" && return nf.free_retransn_to_npool_patch
    k == "SMINN_TO_PLANT_FUN" && return nf.sminn_to_plant_fun_patch
    k == "PLANT_NDEMAND_RETRANS" && return nf.plant_ndemand_retrans_patch
    k == "PLANTCN" && return vs.plantCN_patch
    k == "STORAGE_CDEMAND" && return cs.storage_cdemand_patch
    k == "STORAGE_NDEMAND" && return ns.storage_ndemand_patch
    k == "LEAFC_STORAGE_XFER_ACC" && return cs.leafc_storage_xfer_acc_patch
    k == "LEAFN_STORAGE_XFER_ACC" && return ns.leafn_storage_xfer_acc_patch
    error("no Julia field for $k")
end

"""
    cf_run_cnfun(b; pc, fun_params) -> (; vs, cs, cf, ns, nf, snf, scf, mask_p, mask_c)

Build the CNFUN inputs from a `before_cnfun` dump `b`, run `CLM.cnfun!` and
the two p2c calls of SoilBiogeochemCompetition, return the Julia state.
"""
function cf_run_cnfun(b; pc = ct_pftcon(), ndays_off = 15.0)
    np = length(b["ITYPE"]); nc = length(b["FILTER_C"]); nl = size(b["CROOTFR"], 2)
    mask_p = b["FILTER_P"] .== 1; mask_c = b["FILTER_C"] .== 1
    patch = CLM.PatchData(); patch.itype = copy(b["ITYPE"]); patch.column = copy(b["PCOLUMN"])
    patch.wtcol = copy(b["WTCOL"]); patch.gridcell = ones(Int, np)
    ws = CLM.WaterStateData(); ws.h2osoi_liq_col = copy(b["H2OSOI_LIQ"])
    temp = CLM.TemperatureData(); temp.t_soisno_col = copy(b["T_SOISNO"])
    ss = CLM.SoilStateData(); ss.crootfr_patch = copy(b["CROOTFR"])
    can = CLM.CanopyStateData{Float64}(); CLM.canopystate_init!(can, np); can.tlai_patch .= b["TLAI"]
    vs = CLM.CNVegStateData(); CLM.cnveg_state_init!(vs, np, nc)
    cs = CLM.CNVegCarbonStateData(); CLM.cnveg_carbon_state_init!(cs, np, nc, 1)
    cf = CLM.CNVegCarbonFluxData(); CLM.cnveg_carbon_flux_init!(cf, np, nc, 1)
    ns = CLM.CNVegNitrogenStateData(); CLM.cnveg_nitrogen_state_init!(ns, np, nc, 1)
    nf = CLM.CNVegNitrogenFluxData(); CLM.cnveg_nitrogen_flux_init!(nf, np, nc, 1; nlevdecomp_full = nl)
    snf = CLM.SoilBiogeochemNitrogenFluxData(); CLM.soil_bgc_nitrogen_flux_init!(snf, nc, nl, 7, 10)
    scf = CLM.SoilBiogeochemCarbonFluxData(); CLM.soil_bgc_carbon_flux_init!(scf, nc, nl, 7, 10)
    sns = CLM.SoilBiogeochemNitrogenStateData(); CLM.soil_bgc_nitrogen_state_init!(sns, nc, 1, nl, 7)
    CLM._zero_cnveg_flux_arrays!(cf); CLM._zero_cnveg_flux_arrays!(nf)   # SetValues (CNZero)
    vs.leafcn_offset_patch .= b["LEAFCN_OFFSET"]; vs.plantCN_patch .= b["PLANTCN"]
    vs.onset_flag_patch .= b["ONSET_FLAG"]; vs.offset_flag_patch .= b["OFFSET_FLAG"]
    vs.c_allometry_patch .= b["C_ALLOMETRY"]; vs.n_allometry_patch .= b["N_ALLOMETRY"]
    cf.availc_patch .= b["AVAILC"]; cf.leafc_to_litter_patch .= b["LEAFC_TO_LITTER"]
    cf.leafc_to_litter_fun_patch .= b["LEAFC_TO_LITTER_FUN"]
    cf.prev_leafc_to_litter_patch .= b["PREV_LEAFC_TO_LITTER"]
    cf.leafc_storage_to_xfer_patch .= b["LEAFC_STORAGE_TO_XFER"]
    cs.leafc_patch .= b["LEAFC"]; cs.leafc_storage_patch .= b["LEAFC_STORAGE"]
    cs.frootc_patch .= b["FROOTC"]; cs.frootc_storage_patch .= b["FROOTC_STORAGE"]
    cs.livestemc_patch .= b["LIVESTEMC"]; cs.livecrootc_patch .= b["LIVECROOTC"]
    cs.leafc_storage_xfer_acc_patch .= b["LEAFC_STORAGE_XFER_ACC"]
    cs.storage_cdemand_patch .= b["STORAGE_CDEMAND"]
    ns.leafn_patch .= b["LEAFN"]; ns.leafn_storage_patch .= b["LEAFN_STORAGE"]
    ns.frootn_patch .= b["FROOTN"]; ns.livestemn_patch .= b["LIVESTEMN"]
    ns.livecrootn_patch .= b["LIVECROOTN"]; ns.retransn_patch .= b["RETRANSN"]
    ns.leafn_storage_xfer_acc_patch .= b["LEAFN_STORAGE_XFER_ACC"]
    ns.storage_ndemand_patch .= b["STORAGE_NDEMAND"]
    nf.leafn_storage_to_xfer_patch .= b["LEAFN_STORAGE_TO_XFER"]
    nf.plant_ndemand_patch .= b["PLANT_NDEMAND"]
    snf.smin_no3_to_plant_vr_col[:, 1:nl] .= b["SMIN_NO3_TO_PLANT_VR"]
    snf.smin_nh4_to_plant_vr_col[:, 1:nl] .= b["SMIN_NH4_TO_PLANT_VR"]
    sns.smin_no3_vr_col[:, 1:nl] .= b["SMIN_NO3_VR"]; sns.smin_nh4_vr_col[:, 1:nl] .= b["SMIN_NH4_VR"]
    CLM.cnfun!(mask_p, mask_c, 1:np, 1:nc, CLM.FUNParams(ndays_off = ndays_off),
               CLM.pftcon_fun_from(pc), patch, ws, temp, ss, vs, cs, cf, ns, nf, snf, scf, can, sns;
               dt = CF_DT, nlevdecomp = nl, dzsoi_decomp_vals = cf_dzsoi_decomp(nl),
               use_flexiblecn = true, use_matrixcn = false)
    CLM._fun_p2c!(snf.sminn_to_plant_fun_no3_vr_col, nf.sminn_to_plant_fun_no3_vr_patch,
                  patch, mask_c, 1:nc, 1:np, nl)
    CLM._fun_p2c!(snf.sminn_to_plant_fun_nh4_vr_col, nf.sminn_to_plant_fun_nh4_vr_patch,
                  patch, mask_c, 1:nc, 1:np, nl)
    return (; vs, cs, cf, ns, nf, snf, scf, mask_p, mask_c, np, nc, nl)
end

function cf_check_cnfun(dir; pc = ct_pftcon())
    stats = Dict(k => CTStat() for k in (CF_FUN_OUT_P..., "SMINN_TO_PLANT_FUN_NO3_VR",
        "SMINN_TO_PLANT_FUN_NH4_VR", "SMINN_TO_PLANT_FUN_NO3_VR_COL", "SMINN_TO_PLANT_FUN_NH4_VR_COL",
        "NFIX_TO_SMINN", "SOILC_CHANGE_COL"))
    soy = (nsteps = Ref(0), nfix_pos = Ref(0), nfix_max = Ref(0.0), soyfixn_nonzero = Ref(0))
    nsteps = 0
    for n in cf_steps(dir)
        b = cf_load(dir, "before_cnfun", n); a = cf_load(dir, "after_cnfun", n)
        a === nothing && continue
        nsteps += 1
        r = cf_run_cnfun(b; pc)
        it = b["ITYPE"]
        for p in 1:r.np
            r.mask_p[p] || continue
            tag = "n=$n p=$p itype=$(it[p])"
            for k in CF_FUN_OUT_P
                _cfacc!(stats[k], _cf_jl_out(k, r.vs, r.cs, r.cf, r.ns, r.nf)[p], a[k][p], tag; floor = 1e-20)
            end
            for j in 1:r.nl
                _cfacc!(stats["SMINN_TO_PLANT_FUN_NO3_VR"], r.nf.sminn_to_plant_fun_no3_vr_patch[p, j],
                      a["SMINN_TO_PLANT_FUN_NO3_VR"][p, j], "$tag j=$j"; floor = 1e-20)
                _cfacc!(stats["SMINN_TO_PLANT_FUN_NH4_VR"], r.nf.sminn_to_plant_fun_nh4_vr_patch[p, j],
                      a["SMINN_TO_PLANT_FUN_NH4_VR"][p, j], "$tag j=$j"; floor = 1e-20)
            end
            if it[p] in CF_SOY && b["WTCOL"][p] > 0
                soy.nsteps[] += 1
                a["NFIX"][p] > 0 && (soy.nfix_pos[] += 1)
                soy.nfix_max[] = max(soy.nfix_max[], a["NFIX"][p])
                (a["SOYFIXN"][p] != 0 || b["SOYFIXN"][p] != 0) && (soy.soyfixn_nonzero[] += 1)
            end
        end
        for c in 1:r.nc
            r.mask_c[c] || continue
            for j in 1:r.nl
                _cfacc!(stats["SMINN_TO_PLANT_FUN_NO3_VR_COL"], r.snf.sminn_to_plant_fun_no3_vr_col[c, j],
                      a["SMINN_TO_PLANT_FUN_NO3_VR_COL"][c, j], "n=$n c=$c j=$j"; floor = 1e-20)
                _cfacc!(stats["SMINN_TO_PLANT_FUN_NH4_VR_COL"], r.snf.sminn_to_plant_fun_nh4_vr_col[c, j],
                      a["SMINN_TO_PLANT_FUN_NH4_VR_COL"][c, j], "n=$n c=$c j=$j"; floor = 1e-20)
            end
            _cfacc!(stats["NFIX_TO_SMINN"], r.snf.nfix_to_sminn_col[c], a["NFIX_TO_SMINN"][c], "n=$n c=$c"; floor = 1e-20)
            _cfacc!(stats["SOILC_CHANGE_COL"], r.scf.soilc_change_col[c], a["SOILC_CHANGE_COL"][c], "n=$n c=$c"; floor = 1e-20)
        end
    end
    return (; stats, nsteps, soy = map(x -> x[], soy))
end

ct_pfnc(pc) = CLM.PftConNutrientCompetition{Float64}(
    woody = Float64.(pc.woody), froot_leaf = Float64.(pc.froot_leaf),
    croot_stem = Float64.(pc.croot_stem), stem_leaf = Float64.(pc.stem_leaf),
    flivewd = Float64.(pc.flivewd), leafcn = Float64.(pc.leafcn),
    frootcn = Float64.(pc.frootcn), livewdcn = Float64.(pc.livewdcn),
    deadwdcn = Float64.(pc.deadwdcn), fcur = Float64.(pc.fcur),
    graincn = Float64.(pc.graincn), grperc = Float64.(pc.grperc),
    grpnow = Float64.(pc.grpnow), fleafcn = Float64.(pc.fleafcn),
    ffrootcn = Float64.(pc.ffrootcn), fstemcn = Float64.(pc.fstemcn),
    astemf = Float64.(pc.astemf), season_decid = Float64.(pc.season_decid),
    stress_decid = Float64.(pc.stress_decid), evergreen = Float64.(pc.evergreen))

function _cf_fun_sp()
    sp = CLM.CNSharedParamsData{Float64}(); sp.use_fun = true; sp.use_flexiblecn = true
    return sp
end

_cf_patch(d) = (p = CLM.PatchData(); p.itype = copy(d["ITYPE"]); p.column = copy(d["PCOLUMN"]);
                p.wtcol = copy(d["WTCOL"]); p.gridcell = ones(Int, length(d["ITYPE"])); p)

# ---------------------------------------------------------------------------
# (2) FlexibleCN plant N demand, prognostic-crop call, under FUN
#     (NutrientCompetitionFlexibleCNMod.F90:1393-1680, call_is_for_pcrop = .true.):
#     plant_ndemand = availc * n_allometry/c_allometry, the grain-fill
#     retranslocation, avail_retransn and retransn_to_npool (NOT subtracted from
#     the demand under FUN).
#   inputs : after_cropalloc(n)  (allocation/allometry done, demand not yet run)
#   outputs: plant_ndemand / retransn_to_npool at before_cnfun(n) (FUN not yet run);
#            avail_retransn, *_to_retransn, grain_flag at after_cnalloc(n)
# ---------------------------------------------------------------------------
function cf_check_ndemand(dir; pc = ct_pftcon())
    pfnc = ct_pfnc(pc); sp = _cf_fun_sp()
    keys_d = ("plant_ndemand", "retransn_to_npool", "avail_retransn", "leafn_to_retransn",
              "livestemn_to_retransn", "frootn_to_retransn", "grain_flag")
    stats = Dict(k => CTStat() for k in keys_d)
    ngrainfill_onset = 0
    for n in cf_steps(dir)
        pre = ct_load(dir, "after_cropalloc", n); post = ct_load(dir, "after_cnalloc", n)
        b = cf_load(dir, "before_cnfun", n)
        (pre === nothing || post === nothing || b === nothing) && continue
        haskey(pre, "GRAIN_FLAG") || return (; stats, ngrainfill_onset, available = false)
        it = pre["ITYPE"]; np = length(it); nc = length(pre["FPG"]); nl = 20
        crops = [p for p in 1:np if it[p] >= _CROP_MIN_ITYPE && pre["WTCOL"][p] > 0]
        isempty(crops) && continue
        mask = falses(np); mask[crops] .= true
        patch = _cf_patch(pre)
        cr = CLM.CropData(); CLM.crop_init!(cr, np); cr.croplive_patch .= pre["CROPLIVE"] .== 1
        cr.hui_patch .= pre["HUI"]; cr.gddtsoi_patch .= pre["GDDTSOI"]
        vs = CLM.CNVegStateData(); CLM.cnveg_state_init!(vs, np, nc)
        vs.huileaf_patch .= pre["HUILEAF"]; vs.huigrain_patch .= pre["HUIGRAIN"]
        vs.c_allometry_patch .= pre["C_ALLOMETRY"]; vs.n_allometry_patch .= pre["N_ALLOMETRY"]
        vs.astem_patch .= pre["ASTEM"]; vs.grain_flag_patch .= pre["GRAIN_FLAG"]
        vs.tempsum_potential_gpp_patch .= 0.0; vs.tempmax_retransn_patch .= pre["TEMPMAX_RETRANSN"]
        vs.annsum_potential_gpp_patch .= pre["ANNSUM_POTENTIAL_GPP"]; vs.annmax_retransn_patch .= pre["ANNMAX_RETRANSN"]
        cs = CLM.CNVegCarbonStateData(); CLM.cnveg_carbon_state_init!(cs, np, nc, 1)
        cs.leafc_patch .= pre["LEAFC"]; cs.livestemc_patch .= pre["LIVESTEMC"]; cs.frootc_patch .= pre["FROOTC"]
        ns = CLM.CNVegNitrogenStateData(); CLM.cnveg_nitrogen_state_init!(ns, np, nc, 1)
        ns.leafn_patch .= pre["LEAFN"]; ns.livestemn_patch .= pre["LIVESTEMN"]; ns.frootn_patch .= pre["FROOTN"]
        ns.retransn_patch .= pre["RETRANSN"]
        cf = CLM.CNVegCarbonFluxData(); CLM.cnveg_carbon_flux_init!(cf, np, nc, 1)
        cf.availc_patch .= pre["AVAILC"]; cf.gpp_before_downreg_patch .= 0.0
        nf = CLM.CNVegNitrogenFluxData(); CLM.cnveg_nitrogen_flux_init!(nf, np, nc, 1)
        CLM._zero_cnveg_flux_arrays!(nf)
        can = CLM.CanopyStateData{Float64}(); CLM.canopystate_init!(can, np)
        can.laisun_patch .= pre["LAISUN"]; can.laisha_patch .= pre["LAISHA"]
        sns = CLM.SoilBiogeochemNitrogenStateData(); CLM.soil_bgc_nitrogen_state_init!(sns, nc, 1, nl, 7)
        sns.sminn_vr_col .= 0.0
        scf = CLM.SoilBiogeochemCarbonFluxData(); CLM.soil_bgc_carbon_flux_init!(scf, nc, nl, 7, 10)
        scf.t_scalar_col .= 1.0
        CLM.calc_plant_nutrient_demand_flexiblecn!(mask, 1:np, true, pfnc, sp, patch, cr, can, vs,
            cs, cf, ns, nf, sns, scf; dt = CF_DT, dzsoi_decomp = cf_dzsoi_decomp(nl), nlevdecomp = nl)
        for p in crops
            tag = "n=$n p=$p itype=$(it[p])"
            _cfacc!(stats["plant_ndemand"], nf.plant_ndemand_patch[p], b["PLANT_NDEMAND"][p], tag; floor = 1e-20)
            _cfacc!(stats["retransn_to_npool"], nf.retransn_to_npool_patch[p], b["RETRANSN_TO_NPOOL"][p], tag; floor = 1e-20)
            _cfacc!(stats["avail_retransn"], nf.avail_retransn_patch[p], post["AVAIL_RETRANSN"][p], tag; floor = 1e-20)
            _cfacc!(stats["leafn_to_retransn"], nf.leafn_to_retransn_patch[p], post["LEAFN_TO_RETRANSN"][p], tag; floor = 1e-20)
            _cfacc!(stats["livestemn_to_retransn"], nf.livestemn_to_retransn_patch[p], post["LIVESTEMN_TO_RETRANSN"][p], tag; floor = 1e-20)
            _cfacc!(stats["frootn_to_retransn"], nf.frootn_to_retransn_patch[p], post["FROOTN_TO_RETRANSN"][p], tag; floor = 1e-20)
            _cfacc!(stats["grain_flag"], vs.grain_flag_patch[p], post["GRAIN_FLAG"][p], tag)
            post["LEAFN_TO_RETRANSN"][p] > 0 && (ngrainfill_onset += 1)
        end
    end
    return (; stats, ngrainfill_onset, available = true)
end

# ---------------------------------------------------------------------------
# (3) FlexibleCN allocation under FUN (NutrientCompetitionFlexibleCNMod.F90:439-460):
#     sminn_to_npool = sminn_to_plant_fun, plant_calloc = npp_growth, then the
#     C fluxes and the npool draw.
#   inputs : CNFUN outputs at after_cnfun(n); allocation state at after_cnalloc(n)
#   outputs: after_cnalloc(n)
# ---------------------------------------------------------------------------
function cf_check_alloc(dir; pc = ct_pftcon())
    pfnc = ct_pfnc(pc); sp = _cf_fun_sp()
    keys_a = ("SMINN_TO_NPOOL", "PLANT_NALLOC", "PLANT_CALLOC", "CPOOL_TO_LEAFC", "CPOOL_TO_LEAFC_STORAGE",
              "CPOOL_TO_FROOTC", "CPOOL_TO_FROOTC_STORAGE", "CPOOL_TO_LIVESTEMC", "CPOOL_TO_LIVESTEMC_STORAGE",
              "CPOOL_TO_REPRC", "NPOOL_TO_LEAFN", "NPOOL_TO_LEAFN_STORAGE", "NPOOL_TO_FROOTN",
              "NPOOL_TO_FROOTN_STORAGE", "NPOOL_TO_LIVESTEMN", "NPOOL_TO_LIVESTEMN_STORAGE", "NPOOL_TO_REPRN")
    stats = Dict(k => CTStat() for k in keys_a)
    nsoy = 0
    for n in cf_steps(dir)
        a = cf_load(dir, "after_cnfun", n); d = ct_load(dir, "after_cnalloc", n)
        (a === nothing || d === nothing) && continue
        haskey(d, "SMINN_TO_NPOOL") || return (; stats, nsoy, available = false)
        it = d["ITYPE"]; np = length(it); nc = length(d["FPG"])
        crops = [p for p in 1:np if it[p] >= _CROP_MIN_ITYPE && d["WTCOL"][p] > 0 && d["CROPLIVE"][p] == 1]
        isempty(crops) && continue
        mask = falses(np); mask[crops] .= true
        patch = _cf_patch(d)
        cr = CLM.CropData(); CLM.crop_init!(cr, np); cr.croplive_patch .= d["CROPLIVE"] .== 1
        vs = CLM.CNVegStateData(); CLM.cnveg_state_init!(vs, np, nc)
        vs.aleaf_patch .= d["ALEAF"]; vs.astem_patch .= d["ASTEM"]; vs.aroot_patch .= d["AROOT"]
        vs.arepr_patch[:, end] .= d["AREPR"]
        vs.c_allometry_patch .= d["C_ALLOMETRY"]; vs.n_allometry_patch .= d["N_ALLOMETRY"]
        cs = CLM.CNVegCarbonStateData(); CLM.cnveg_carbon_state_init!(cs, np, nc, 1)
        cs.leafc_storage_patch .= d["LEAFC_STORAGE"]; cs.frootc_storage_patch .= d["FROOTC_STORAGE"]
        cs.livestemc_storage_patch .= d["LIVESTEMC_STORAGE"]; cs.livecrootc_storage_patch .= 0.0
        ns = CLM.CNVegNitrogenStateData(); CLM.cnveg_nitrogen_state_init!(ns, np, nc, 1)
        ns.npool_patch .= d["NPOOL"]
        ns.leafn_storage_patch .= d["LEAFN_STORAGE"]; ns.frootn_storage_patch .= d["FROOTN_STORAGE"]
        ns.livestemn_storage_patch .= d["LIVESTEMN_STORAGE"]; ns.livecrootn_storage_patch .= 0.0
        cf = CLM.CNVegCarbonFluxData(); CLM.cnveg_carbon_flux_init!(cf, np, nc, 1)
        cf.availc_patch .= d["AVAILC"]; cf.npp_growth_patch .= a["NPP_GROWTH"]
        cf.annsum_npp_patch .= 0.0; cf.gpp_before_downreg_patch .= 1.0
        nf = CLM.CNVegNitrogenFluxData(); CLM.cnveg_nitrogen_flux_init!(nf, np, nc, 1)
        nf.sminn_to_plant_fun_patch .= a["SMINN_TO_PLANT_FUN"]; nf.retransn_to_npool_patch .= a["RETRANSN_TO_NPOOL"]
        nf.plant_ndemand_patch .= a["PLANT_NDEMAND"]
        can = CLM.CanopyStateData{Float64}(); CLM.canopystate_init!(can, np)
        can.laisun_patch .= d["LAISUN"]; can.laisha_patch .= d["LAISHA"]
        CLM.calc_plant_cn_alloc_flexiblecn!(mask, 1:np, pfnc, sp, patch, cr, can, vs, cs, cf, ns, nf;
                                 fpg_col = copy(d["FPG"]), dt = CF_DT, carbon_resp_opt = 0)
        jl = Dict("SMINN_TO_NPOOL" => nf.sminn_to_npool_patch, "PLANT_NALLOC" => nf.plant_nalloc_patch,
                  "PLANT_CALLOC" => cf.plant_calloc_patch, "CPOOL_TO_LEAFC" => cf.cpool_to_leafc_patch,
                  "CPOOL_TO_LEAFC_STORAGE" => cf.cpool_to_leafc_storage_patch,
                  "CPOOL_TO_FROOTC" => cf.cpool_to_frootc_patch, "CPOOL_TO_FROOTC_STORAGE" => cf.cpool_to_frootc_storage_patch,
                  "CPOOL_TO_LIVESTEMC" => cf.cpool_to_livestemc_patch,
                  "CPOOL_TO_LIVESTEMC_STORAGE" => cf.cpool_to_livestemc_storage_patch,
                  "CPOOL_TO_REPRC" => cf.cpool_to_reproductivec_patch[:, end],
                  "NPOOL_TO_LEAFN" => nf.npool_to_leafn_patch, "NPOOL_TO_LEAFN_STORAGE" => nf.npool_to_leafn_storage_patch,
                  "NPOOL_TO_FROOTN" => nf.npool_to_frootn_patch, "NPOOL_TO_FROOTN_STORAGE" => nf.npool_to_frootn_storage_patch,
                  "NPOOL_TO_LIVESTEMN" => nf.npool_to_livestemn_patch,
                  "NPOOL_TO_LIVESTEMN_STORAGE" => nf.npool_to_livestemn_storage_patch,
                  "NPOOL_TO_REPRN" => nf.npool_to_reproductiven_patch[:, end])
        for p in crops
            it[p] in CF_SOY && (nsoy += 1)
            for k in keys_a
                _cfacc!(stats[k], jl[k][p], d[k][p], "n=$n p=$p itype=$(it[p])"; floor = 1e-20)
            end
        end
    end
    return (; stats, nsoy, available = true)
end

# ---------------------------------------------------------------------------
# (4) npool update (CNNStateUpdate1Mod.F90:317-386, crop branch):
#     NPOOL end_noleach(n) = NPOOL after_cnalloc(n) + (sminn_to_npool +
#     retransn_to_npool + free_retransn_to_npool - npool_to_*) * dt,
#     driven through CLM.n_state_update1! with every other flux zero.
# ---------------------------------------------------------------------------
function cf_check_npool(dir)
    stats = Dict("NPOOL" => CTStat(), "RETRANSN" => CTStat())
    for n in cf_steps(dir)
        d = ct_load(dir, "after_cnalloc", n); e = ct_load(dir, "end_noleach", n)
        (d === nothing || e === nothing) && continue
        haskey(d, "SMINN_TO_NPOOL") || return (; stats, available = false)
        it = d["ITYPE"]; np = length(it); nc = length(d["FPG"]); nl = 20
        crops = [p for p in 1:np if it[p] >= _CROP_MIN_ITYPE && d["WTCOL"][p] > 0]
        isempty(crops) && continue
        mask = falses(np); mask[crops] .= true
        ns = CLM.CNVegNitrogenStateData(); CLM.cnveg_nitrogen_state_init!(ns, np, nc, 1)
        for f in fieldnames(typeof(ns))
            x = getfield(ns, f); (x isa AbstractArray{<:AbstractFloat} && !isempty(x)) && fill!(x, 0.0)
        end
        ns.npool_patch .= d["NPOOL"]; ns.retransn_patch .= d["RETRANSN"]
        nf = CLM.CNVegNitrogenFluxData(); CLM.cnveg_nitrogen_flux_init!(nf, np, nc, 1)
        CLM._zero_cnveg_flux_arrays!(nf)
        for (k, f) in (("SMINN_TO_NPOOL", :sminn_to_npool_patch), ("RETRANSN_TO_NPOOL", :retransn_to_npool_patch),
                       ("FREE_RETRANSN_TO_NPOOL", :free_retransn_to_npool_patch),
                       ("NPOOL_TO_LEAFN", :npool_to_leafn_patch), ("NPOOL_TO_LEAFN_STORAGE", :npool_to_leafn_storage_patch),
                       ("NPOOL_TO_FROOTN", :npool_to_frootn_patch), ("NPOOL_TO_FROOTN_STORAGE", :npool_to_frootn_storage_patch),
                       ("NPOOL_TO_LIVESTEMN", :npool_to_livestemn_patch),
                       ("NPOOL_TO_LIVESTEMN_STORAGE", :npool_to_livestemn_storage_patch),
                       ("LEAFN_TO_RETRANSN", :leafn_to_retransn_patch), ("LIVESTEMN_TO_RETRANSN", :livestemn_to_retransn_patch),
                       ("FROOTN_TO_RETRANSN", :frootn_to_retransn_patch))
            getfield(nf, f) .= e[k]
        end
        nf.npool_to_reproductiven_patch[:, end] .= e["NPOOL_TO_REPRN"]
        nf.npool_to_reproductiven_storage_patch[:, end] .= e["NPOOL_TO_REPRN_STOR"]
        snf = CLM.SoilBiogeochemNitrogenFluxData(); CLM.soil_bgc_nitrogen_flux_init!(snf, nc, nl, 7, 10)
        for f in fieldnames(typeof(snf))
            x = getfield(snf, f); (x isa AbstractArray{<:AbstractFloat} && !isempty(x)) && fill!(x, 0.0)
        end
        CLM.n_state_update1!(ns, nf, snf; mask_soilc = falses(nc), mask_soilp = mask, bounds_col = 1:nc,
            bounds_patch = 1:np, ivt = it, woody = zeros(80), col_is_fates = falses(nc), nlevdecomp = nl,
            i_litr_min = 1, i_litr_max = 3, i_cwd = 4, npcropmin = _CROP_MIN_ITYPE, nrepr = 1,
            repr_grain_min = 1, repr_grain_max = 1, use_fun = true, dt = CF_DT)
        for p in crops
            tag = "n=$n p=$p itype=$(it[p])"
            _cfacc!(stats["NPOOL"], ns.npool_patch[p], e["NPOOL"][p], tag; floor = 1e-12)
            _cfacc!(stats["RETRANSN"], ns.retransn_patch[p], e["RETRANSN"][p], tag; floor = 1e-12)
        end
    end
    return (; stats, available = true)
end

# ---------------------------------------------------------------------------
# (5) SOYFIXN and FPG under FUN. CTSM never calls CNSoyfix with FUN
#     (CNDriverMod.F90:324) and SetValues zeroes soyfixn every step, so SOYFIXN is
#     a computed 0 on every soybean patch-step (soybean fixation is CNFUN's Nfix);
#     FUN never writes fpg ("FUN has no concept of FPG",
#     SoilBiogeochemCompetitionMod.F90:978) so FPG keeps its allocation value.
# ---------------------------------------------------------------------------
function cf_check_soyfixn_fpg(dir)
    nsoy = 0; soyfixn_nonzero = 0; fpg_vals = Set{String}()
    for n in cf_steps(dir)
        for lab in ("after_soyfix", "after_cnalloc", "end_noleach")
            d = ct_load(dir, lab, n); d === nothing && continue
            for p in eachindex(d["ITYPE"])
                d["ITYPE"][p] in CF_SOY && d["WTCOL"][p] > 0 || continue
                nsoy += 1; d["SOYFIXN"][p] != 0 && (soyfixn_nonzero += 1)
            end
            for c in eachindex(d["FPG"]); push!(fpg_vals, string(d["FPG"][c])); end
        end
    end
    return (; nsoy, soyfixn_nonzero, fpg_vals = collect(fpg_vals))
end

# ---------------------------------------------------------------------------
# (6) harvest step of CNOffsetLitterfall WITH FUN: crop harvest branch + the FUN
#     clamp and leaf-N split (CNPhenologyMod.F90:3720-3858), which reads CNFUN's
#     retransn_to_npool / free_retransn_to_npool of the same step.
# ---------------------------------------------------------------------------
function cf_check_harvest(dir; pc = ct_pftcon(), dt = CF_DT)
    keys_h = ("LEAFC_TO_LITTER", "FROOTC_TO_LITTER", "LEAFN_TO_LITTER", "LEAFN_TO_RETRANSN",
              "FROOTN_TO_LITTER", "LEAFCN_OFFSET", "LEAFC_TO_LITTER_FUN", "REPR_GRAINC_TO_FOOD",
              "REPR_GRAINN_TO_FOOD", "LIVESTEMN_TO_LITTER")
    stats = Dict(k => CTStat() for k in keys_h)
    events = String[]
    pf = CLM.PftConPhenology{Float64}(
        lflitcn = Float64.(pc.lflitcn), leafcn = Float64.(pc.leafcn), frootcn = Float64.(pc.frootcn),
        biofuel_harvfrac = Float64.(pc.biofuel_harvfrac),
        repr_structure_harvfrac = Float64.(pc.repr_structure_harvfrac))
    pstate = CLM.PhenologyState(); pstate.dt = dt
    for n in cf_steps(dir)
        post = ct_load(dir, "end_noleach", n); post === nothing && continue
        it = post["ITYPE"]; np = length(it); nc = length(post["FPG"])
        hv = [p for p in 1:np if it[p] >= _CROP_MIN_ITYPE && post["WTCOL"][p] > 0 &&
              post["OFFSET_FLAG"][p] == 1.0 && abs(post["OFFSET_COUNTER"][p] - dt) <= dt / 2]
        isempty(hv) && continue
        pre = ct_load(dir, "after_cnalloc", n); a = cf_load(dir, "after_cnfun", n)
        (pre === nothing || a === nothing) && continue
        patch = _cf_patch(post)
        vs = CLM.CNVegStateData(); CLM.cnveg_state_init!(vs, np, nc)
        vs.offset_flag_patch .= post["OFFSET_FLAG"]; vs.offset_counter_patch .= post["OFFSET_COUNTER"]
        vs.leafcn_offset_patch .= pre["LEAFCN_OFFSET"]
        cr = CLM.CropData(); CLM.crop_init!(cr, np); cr.harvest_count .= post["HARVEST_COUNT"]
        cs = CLM.CNVegCarbonStateData(); CLM.cnveg_carbon_state_init!(cs, np, nc, 1)
        ns = CLM.CNVegNitrogenStateData(); CLM.cnveg_nitrogen_state_init!(ns, np, nc, 1)
        cf = CLM.CNVegCarbonFluxData(); CLM.cnveg_carbon_flux_init!(cf, np, nc, 1)
        nf = CLM.CNVegNitrogenFluxData(); CLM.cnveg_nitrogen_flux_init!(nf, np, nc, 1)
        CLM._zero_cnveg_flux_arrays!(cf); CLM._zero_cnveg_flux_arrays!(nf)
        cf.leafc_to_litter_fun_patch .= pre["LEAFC_TO_LITTER_FUN"]
        cf.prev_leafc_to_litter_patch .= pre["PREV_LEAFC_TO_LITTER"]; cf.prev_frootc_to_litter_patch .= pre["PREV_FROOTC_TO_LITTER"]
        cs.leafc_patch .= pre["LEAFC"]; cs.frootc_patch .= pre["FROOTC"]; cs.livestemc_patch .= pre["LIVESTEMC"]
        cs.reproductivec_patch[:, end] .= pre["REPRC"]; cs.cropseedc_deficit_patch .= pre["CROPSEEDC_DEFICIT"]
        ns.leafn_patch .= pre["LEAFN"]; ns.frootn_patch .= pre["FROOTN"]; ns.livestemn_patch .= pre["LIVESTEMN"]
        ns.reproductiven_patch[:, end] .= pre["REPRN"]; ns.cropseedn_deficit_patch .= pre["CROPSEEDN_DEFICIT"]
        cf.cpool_to_leafc_patch .= pre["CPOOL_TO_LEAFC"]; cf.cpool_to_livestemc_patch .= pre["CPOOL_TO_LIVESTEMC"]
        cf.cpool_to_frootc_patch .= pre["CPOOL_TO_FROOTC"]; cf.cpool_to_reproductivec_patch[:, end] .= pre["CPOOL_TO_REPRC"]
        nf.npool_to_reproductiven_patch[:, end] .= pre["NPOOL_TO_REPRN"]
        nf.retransn_to_npool_patch .= a["RETRANSN_TO_NPOOL"]; nf.free_retransn_to_npool_patch .= a["FREE_RETRANSN_TO_NPOOL"]
        mask = falses(np); mask[hv] .= true
        CLM.cn_offset_litterfall!(pstate, mask, pf, vs, cs, ns, cf, nf, cr, patch;
                                  mask_pcropp = mask, use_fun = true, CNratio_floating = true)
        jl = Dict("LEAFC_TO_LITTER" => cf.leafc_to_litter_patch, "FROOTC_TO_LITTER" => cf.frootc_to_litter_patch,
                  "LEAFN_TO_LITTER" => nf.leafn_to_litter_patch, "LEAFN_TO_RETRANSN" => nf.leafn_to_retransn_patch,
                  "FROOTN_TO_LITTER" => nf.frootn_to_litter_patch, "LEAFCN_OFFSET" => vs.leafcn_offset_patch,
                  "LEAFC_TO_LITTER_FUN" => cf.leafc_to_litter_fun_patch,
                  "REPR_GRAINC_TO_FOOD" => cf.repr_grainc_to_food_patch[:, end],
                  "REPR_GRAINN_TO_FOOD" => nf.repr_grainn_to_food_patch[:, end],
                  "LIVESTEMN_TO_LITTER" => nf.livestemn_to_litter_patch)
        for p in hv
            push!(events, "n=$n ymd=$(post["ymd"]) itype=$(it[p])")
            for k in keys_h
                _cfacc!(stats[k], jl[k][p], post[k][p], "n=$n p=$p itype=$(it[p])"; floor = 1e-15)
            end
        end
    end
    return (; stats, events)
end

# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
function cf_report(name, stats)
    println("== $name ==")
    for k in sort!(collect(keys(stats)))
        s = stats[k]
        @printf("  %-30s n=%6d  scale=%.3e  maxabs=%.3e  maxabs/scale=%.2e  spval-flips=%d  %s\n",
                k, s.n, cf_scale(s), s.maxabs, cf_normerr(s), cf_flips(s), s.maxabs > 0 ? s.worst : "")
    end
end

if get(ENV, "CROP_FUN_NO_MAIN", "0") != "1" && abspath(PROGRAM_FILE) == @__FILE__
    ref = length(ARGS) >= 1 ? ARGS[1] : CF_REF
    r1 = cf_check_cnfun(ref)
    println("cnfun! steps: ", r1.nsteps, "   soybean: ", r1.soy)
    cf_report("cnfun! (before_cnfun -> after_cnfun)", r1.stats)
    r2 = cf_check_ndemand(ref)
    r2.available ? (println("grain-fill retranslocation onsets: ", r2.ngrainfill_onset);
                    cf_report("FlexibleCN N demand, pcrop call, FUN", r2.stats)) :
                   println("(2) skipped: dumps predate GRAIN_FLAG")
    r3 = cf_check_alloc(ref)
    r3.available ? (println("soybean alloc patch-steps: ", r3.nsoy); cf_report("FlexibleCN allocation, FUN", r3.stats)) :
                   println("(3) skipped: dumps predate the FUN alloc fields")
    r4 = cf_check_npool(ref)
    r4.available ? cf_report("npool / retransn update (n_state_update1!)", r4.stats) : println("(4) skipped")
    r5 = cf_check_soyfixn_fpg(ref)
    println("SOYFIXN under FUN: ", r5.nsoy, " soybean patch-dumps, ", r5.soyfixn_nonzero, " non-zero; FPG values seen: ", r5.fpg_vals)
    r6 = cf_check_harvest(ref)
    println("harvests: ", r6.events)
    cf_report("CNOffsetLitterfall harvest step, FUN", r6.stats)
end
