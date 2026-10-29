# Crop-tail Fortran parity: grain-fill allocation + grain C flux, CNSoyfix,
# vernalization -- every step of the instrumented ctsm5.3.012 dumps
# (scripts/validation/crop_parity/cropdumpMod.F90), with Fortran inputs INJECTED
# into the Julia routines and the Julia outputs compared with the Fortran outputs
# of the same step (the restart-injection method of PR #314, at every step).
#
#   julia --project=. scripts/fortran_parity_crop_tails.jl [REF_A] [REF_SOY]
#
#   REF_A   cold crop run, dumps across grain fill / harvest / vernalization
#           (default $SYMFLUENCE_DATA/runs/crop_tails_refA)
#   REF_SOY continue run with the in-model CNSoyfix sweep
#           (default $SYMFLUENCE_DATA/runs/crop_tails_refB_soyfix)
#
# Reproduction recipe: docs/CROP_PARITY.md § "Crop tails (2026-10-28)".
# Also `include`d by test/test_crop_tails_parity.jl (functions only; set
# CROP_TAILS_NO_MAIN=1 before including to skip the report).
using CLM
using NCDatasets
using Printf

const _SD = get(ENV, "SYMFLUENCE_DATA", "/Users/darri/repos/SYMFLUENCE_data")
const CT_PARAMFILE = get(ENV, "CROP_TAILS_PARAMFILE", joinpath(_SD,
    "domain_Bow_at_Banff_lumped/optimization/CLM/dds_run_1/final_evaluation/settings/CLM/parameters/clm5_params.nc"))

# ---------------------------------------------------------------------------
# dump I/O
# ---------------------------------------------------------------------------
function ct_load(dir::AbstractString, label::AbstractString, n::Integer)
    f = joinpath(dir, "cropdump_$(label)_n$(n).nc")
    isfile(f) || return nothing
    d = Dict{String,Any}()
    NCDataset(f) do ds
        for k in keys(ds)
            v = ds[k]
            d[k] = ndims(v) == 0 ? v[] : Array(v[:])
        end
    end
    return d
end

function ct_steps(dir::AbstractString, label::AbstractString)
    rx = Regex("^cropdump_$(label)_n(\\d+)\\.nc\$")
    s = Int[]
    for f in readdir(dir)
        m = match(rx, f)
        m === nothing || push!(s, parse(Int, m[1]))
    end
    return sort!(s)
end

function ct_pftcon(paramfile::AbstractString = CT_PARAMFILE)
    p = CLM.PftconType()
    NCDataset(paramfile) do ds
        CLM.pftcon_read!(p, ds)
    end
    return p
end

ct_pfta(pc) = CLM.PftConAllocation{Float64}(
    woody = Float64.(pc.woody), froot_leaf = Float64.(pc.froot_leaf),
    croot_stem = Float64.(pc.croot_stem), stem_leaf = Float64.(pc.stem_leaf),
    flivewd = Float64.(pc.flivewd), leafcn = Float64.(pc.leafcn),
    frootcn = Float64.(pc.frootcn), livewdcn = Float64.(pc.livewdcn),
    deadwdcn = Float64.(pc.deadwdcn), graincn = Float64.(pc.graincn),
    grperc = Float64.(pc.grperc),
    arooti = Float64.(pc.arooti), arootf = Float64.(pc.arootf),
    fleafi = Float64.(pc.fleafi), bfact = Float64.(pc.bfact),
    astemf = Float64.(pc.astemf), declfact = Float64.(pc.declfact),
    allconss = Float64.(pc.allconss), aleaff = Float64.(pc.aleaff),
    allconsl = Float64.(pc.allconsl))

# relative error with an absolute floor (allocation fractions are O(1); fluxes are
# compared relative to the largest magnitude present)
_relerr(a, b; floor = 1e-30) = abs(a - b) / max(abs(b), floor)

mutable struct CTStat
    n::Int
    maxabs::Float64
    maxrel::Float64
    worst::String
end
CTStat() = CTStat(0, 0.0, 0.0, "")
function _acc!(s::CTStat, jl, f, tag; floor = 1e-12)
    (isfinite(f) || isfinite(jl)) || return
    s.n += 1
    ae = abs(jl - f)
    re = ae / max(abs(f), floor)
    s.maxabs = max(s.maxabs, ae)
    if re > s.maxrel
        s.maxrel = re
        s.worst = tag
    end
end

const _CROP_MIN_ITYPE = 17   # npcropmin (0-based Fortran itype)

# ---------------------------------------------------------------------------
# (1) calc_crop_allocation_fractions!  (all phases, incl. the grain-fill branch)
#
# inputs  : hui/gddtsoi/huileaf/huigrain/gddmaturity/peaklai/croplive at
#           'after_cropalloc'(n) (what allocation read -- nothing between the
#           allocation call and the dump touches them); the PRIOR aleaf/astem/
#           aleafi/astemi at 'end_noleach'(n-1) (CropPhenology of step n-1, which
#           runs after allocation in the non-FUN driver, may re-plant and reset them)
# outputs : aleaf/astem/aroot/arepr/aleafi/astemi at 'after_cropalloc'(n)
# chained : additionally, carry Julia's OWN aleaf/astem through each contiguous
#           grain-fill run (only hui etc. injected) -- tests the recursion itself.
# ---------------------------------------------------------------------------
function ct_check_allocation(dir; pfta = ct_pfta(ct_pftcon()), verbose = false)
    steps = ct_steps(dir, "after_cropalloc")
    stats = Dict(k => CTStat() for k in ("aleaf", "astem", "aroot", "arepr", "aleafi", "astemi"))
    chained = Dict(k => CTStat() for k in ("aleaf", "astem", "arepr"))
    phase_count = Dict(1.0 => 0, 2.0 => 0, 3.0 => 0)
    chain_state = Dict{Int,Tuple{Float64,Float64}}()   # patch => (aleaf, astem) Julia-carried
    last_n = -10
    for n in steps
        cur = ct_load(dir, "after_cropalloc", n)
        prev = ct_load(dir, "end_noleach", n - 1)
        prev === nothing && (chain_state = Dict{Int,Tuple{Float64,Float64}}(); last_n = n; continue)
        n == last_n + 1 || (chain_state = Dict{Int,Tuple{Float64,Float64}}())
        last_n = n
        np = length(cur["ITYPE"])
        it = cur["ITYPE"]
        crops = [p for p in 1:np if it[p] >= _CROP_MIN_ITYPE && cur["WTCOL"][p] > 0]

        patch = CLM.PatchData(); patch.itype = copy(it); patch.column = copy(cur["PCOLUMN"])
        patch.wtcol = copy(cur["WTCOL"]); patch.gridcell = ones(Int, np)
        cr = CLM.CropData(); CLM.crop_init!(cr, np)
        vs = CLM.CNVegStateData(); CLM.cnveg_state_init!(vs, np, length(cur["FPG"]))
        cr.croplive_patch .= cur["CROPLIVE"] .== 1
        cr.hui_patch .= cur["HUI"]; cr.gddtsoi_patch .= cur["GDDTSOI"]
        vs.huileaf_patch .= cur["HUILEAF"]; vs.huigrain_patch .= cur["HUIGRAIN"]
        vs.gddmaturity_patch .= cur["GDDMATURITY"]; vs.peaklai_patch .= cur["PEAKLAI"]
        vs.aleaf_patch .= prev["ALEAF"]; vs.astem_patch .= prev["ASTEM"]
        vs.aleafi_patch .= prev["ALEAFI"]; vs.astemi_patch .= prev["ASTEMI"]
        vs.aroot_patch .= prev["AROOT"]; vs.arepr_patch[:, end] .= prev["AREPR"]
        mask = falses(np); mask[crops] .= true
        # the injected-prior pass
        CLM.calc_crop_allocation_fractions!(mask, 1:np, pfta, patch, cr, vs)
        for p in crops
            cr.croplive_patch[p] || continue
            ph = cur["HUI"][p] >= cur["HUIGRAIN"][p] ? 3.0 :
                 cur["GDDTSOI"][p] >= cur["HUILEAF"][p] ? 2.0 : 1.0
            phase_count[ph] += 1
            tag = "n=$n p=$p itype=$(it[p]) phase=$(Int(ph))"
            _acc!(stats["aleaf"], vs.aleaf_patch[p], cur["ALEAF"][p], tag)
            _acc!(stats["astem"], vs.astem_patch[p], cur["ASTEM"][p], tag)
            _acc!(stats["aroot"], vs.aroot_patch[p], cur["AROOT"][p], tag)
            _acc!(stats["arepr"], vs.arepr_patch[p, end], cur["AREPR"][p], tag)
            _acc!(stats["aleafi"], vs.aleafi_patch[p], cur["ALEAFI"][p], tag)
            _acc!(stats["astemi"], vs.astemi_patch[p], cur["ASTEMI"][p], tag)
        end
        # the chained pass: Julia-carried prior aleaf/astem through grain fill
        for p in crops
            (cr.croplive_patch[p] && cur["HUI"][p] >= cur["HUIGRAIN"][p]) || (delete!(chain_state, p); continue)
            if haskey(chain_state, p)
                vs2 = deepcopy(vs)
                vs2.aleaf_patch[p], vs2.astem_patch[p] = chain_state[p]
                vs2.aleafi_patch[p] = prev["ALEAFI"][p]; vs2.astemi_patch[p] = prev["ASTEMI"][p]
                m1 = falses(np); m1[p] = true
                CLM.calc_crop_allocation_fractions!(m1, 1:np, pfta, patch, cr, vs2)
                tag = "n=$n p=$p itype=$(it[p]) chained"
                _acc!(chained["aleaf"], vs2.aleaf_patch[p], cur["ALEAF"][p], tag)
                _acc!(chained["astem"], vs2.astem_patch[p], cur["ASTEM"][p], tag)
                _acc!(chained["arepr"], vs2.arepr_patch[p, end], cur["AREPR"][p], tag)
                chain_state[p] = (vs2.aleaf_patch[p], vs2.astem_patch[p])
            else
                chain_state[p] = (cur["ALEAF"][p], cur["ASTEM"][p])   # seed from Fortran once
            end
        end
    end
    return (; stats, chained, phase_count)
end

# ---------------------------------------------------------------------------
# (2) calc_allometry! (c_allometry / n_allometry for crops, non-FUN g1 = 0.25)
#     and (3) the grain C flux cpool_to_reproductivec(+_storage) and the other
#     crop cpool_to_* fluxes, via calc_plant_cn_alloc! (plant_calloc = availc in
#     the non-FUN path).  inputs/outputs at 'after_cnalloc'(n).
# ---------------------------------------------------------------------------
function ct_check_cnalloc(dir; pc = ct_pftcon())
    pfta = ct_pfta(pc)
    pfnc = CLM.PftConNutrientCompetition{Float64}(
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
    sp = CLM.CNSharedParamsData{Float64}(); sp.use_fun = false
    sp.use_flexiblecn = true      # the reference runs use_flexiblecn=.true., use_fun=.false.
    keys_c = ("c_allometry", "n_allometry", "plant_calloc", "cpool_to_reprc", "cpool_to_reprc_stor",
              "cpool_to_leafc", "cpool_to_livestemc", "cpool_to_frootc")
    stats = Dict(k => CTStat() for k in keys_c)
    ngrain = 0
    for n in ct_steps(dir, "after_cnalloc")
        d = ct_load(dir, "after_cnalloc", n)
        np = length(d["ITYPE"]); nc = length(d["FPG"]); it = d["ITYPE"]
        crops = [p for p in 1:np if it[p] >= _CROP_MIN_ITYPE && d["WTCOL"][p] > 0 && d["CROPLIVE"][p] == 1]
        isempty(crops) && continue
        mask = falses(np); mask[crops] .= true
        patch = CLM.PatchData(); patch.itype = copy(it); patch.column = copy(d["PCOLUMN"])
        patch.wtcol = copy(d["WTCOL"]); patch.gridcell = ones(Int, np)
        cr = CLM.CropData(); CLM.crop_init!(cr, np); cr.croplive_patch .= d["CROPLIVE"] .== 1
        vs = CLM.CNVegStateData(); CLM.cnveg_state_init!(vs, np, nc)
        vs.aleaf_patch .= d["ALEAF"]; vs.astem_patch .= d["ASTEM"]; vs.aroot_patch .= d["AROOT"]
        vs.arepr_patch[:, end] .= d["AREPR"]
        cf = CLM.CNVegCarbonFluxData(); CLM.cnveg_carbon_flux_init!(cf, np, nc, 1)
        cf.annsum_npp_patch .= 0.0
        CLM.calc_allometry!(mask, 1:np, pfta, sp, patch, cf, vs)
        for p in crops
            tag = "n=$n p=$p itype=$(it[p])"
            _acc!(stats["c_allometry"], vs.c_allometry_patch[p], d["C_ALLOMETRY"][p], tag)
            _acc!(stats["n_allometry"], vs.n_allometry_patch[p], d["N_ALLOMETRY"][p], tag)
        end
        # allocation fluxes from the Fortran allometry (isolates calc_plant_cn_alloc)
        vs.c_allometry_patch .= d["C_ALLOMETRY"]; vs.n_allometry_patch .= d["N_ALLOMETRY"]
        cs = CLM.CNVegCarbonStateData(); CLM.cnveg_carbon_state_init!(cs, np, nc, 1)
        nf = CLM.CNVegNitrogenFluxData(); CLM.cnveg_nitrogen_flux_init!(nf, np, nc, 1)
        cf.availc_patch .= d["AVAILC"]
        nf.plant_ndemand_patch .= d["PLANT_NDEMAND"]; nf.retransn_to_npool_patch .= 0.0
        nf.sminn_to_plant_fun_patch .= 0.0
        cf.gpp_before_downreg_patch .= 1.0; cf.npp_growth_patch .= 0.0
        ns = CLM.CNVegNitrogenStateData(); CLM.cnveg_nitrogen_state_init!(ns, np, nc, 1)
        ns.npool_patch .= d["NPOOL"]
        for f in (:leafc_storage_patch, :frootc_storage_patch, :livestemc_storage_patch, :livecrootc_storage_patch)
            getfield(cs, f) .= 0.0
        end
        for f in (:leafn_storage_patch, :frootn_storage_patch, :livestemn_storage_patch, :livecrootn_storage_patch)
            getfield(ns, f) .= 0.0
        end
        can = CLM.CanopyStateData{Float64}(); CLM.canopystate_init!(can, np)
        can.laisun_patch .= 1.0; can.laisha_patch .= 1.0
        # FlexibleCN path (the driver dispatches on use_flexiblecn): non-FUN
        # plant_calloc = availc, no GPP downregulation.
        CLM.calc_plant_cn_alloc_flexiblecn!(mask, 1:np, pfnc, sp, patch, cr, can, vs, cs, cf, ns, nf;
                                 fpg_col = copy(d["FPG"]), dt = 3600.0, carbon_resp_opt = 0)
        for p in crops
            tag = "n=$n p=$p itype=$(it[p])"
            _acc!(stats["plant_calloc"], cf.plant_calloc_patch[p], d["PLANT_CALLOC"][p], tag)
            _acc!(stats["cpool_to_reprc"], cf.cpool_to_reproductivec_patch[p, end], d["CPOOL_TO_REPRC"][p], tag; floor = 1e-15)
            _acc!(stats["cpool_to_reprc_stor"], cf.cpool_to_reproductivec_storage_patch[p, end], d["CPOOL_TO_REPRC_STOR"][p], tag; floor = 1e-15)
            _acc!(stats["cpool_to_leafc"], cf.cpool_to_leafc_patch[p], d["CPOOL_TO_LEAFC"][p], tag; floor = 1e-15)
            _acc!(stats["cpool_to_livestemc"], cf.cpool_to_livestemc_patch[p], d["CPOOL_TO_LIVESTEMC"][p], tag; floor = 1e-15)
            _acc!(stats["cpool_to_frootc"], cf.cpool_to_frootc_patch[p], d["CPOOL_TO_FROOTC"][p], tag; floor = 1e-15)
            d["CPOOL_TO_REPRC"][p] > 0 && (ngrain += 1)
        end
    end
    return (; stats, ngrain)
end

# ---------------------------------------------------------------------------
# (4) cn_offset_litterfall! on the harvest (final offset) step -- the Julia
#      crop-harvest branch (grain -> food/seed C and N, leaf/livestem/froot litter,
#      biofuel) against CNOffsetLitterfall.
#   inputs : pools and this step's allocation fluxes at 'after_cnalloc'(n) (state
#            is not updated until CStateUpdate1, after phenology); offset flag/
#            counter and harvest_count at 'end_noleach'(n) (CropPhenology sets them
#            earlier in the same CNPhenology call)
#   outputs: the harvest fluxes at 'end_noleach'(n)
# ---------------------------------------------------------------------------
function ct_check_harvest(dir; pc = ct_pftcon(), dt = 3600.0)
    keys_h = ("repr_grainc_to_food", "repr_grainc_to_seed", "repr_grainn_to_food", "repr_grainn_to_seed",
              "leafc_to_litter", "frootc_to_litter", "livestemc_to_litter", "livestemn_to_litter",
              "leafn_to_litter", "leafn_to_retransn", "frootn_to_litter", "leafc_to_biofuelc")
    stats = Dict(k => CTStat() for k in keys_h)
    events = String[]
    pf = CLM.PftConPhenology{Float64}(
        lflitcn = Float64.(pc.lflitcn), leafcn = Float64.(pc.leafcn), frootcn = Float64.(pc.frootcn),
        biofuel_harvfrac = Float64.(pc.biofuel_harvfrac),
        repr_structure_harvfrac = Float64.(pc.repr_structure_harvfrac))
    pstate = CLM.PhenologyState(); pstate.dt = dt
    for n in ct_steps(dir, "end_noleach")
        post = ct_load(dir, "end_noleach", n)
        haskey(post, "LEAFC_TO_LITTER") || return (; stats, events, available = false)
        it = post["ITYPE"]; np = length(it); nc = length(post["FPG"])
        hv = [p for p in 1:np if it[p] >= _CROP_MIN_ITYPE && post["WTCOL"][p] > 0 &&
              post["OFFSET_FLAG"][p] == 1.0 && abs(post["OFFSET_COUNTER"][p] - dt) <= dt / 2]
        isempty(hv) && continue
        pre = ct_load(dir, "after_cnalloc", n); pre === nothing && continue
        patch = CLM.PatchData(); patch.itype = copy(it); patch.column = copy(post["PCOLUMN"])
        patch.wtcol = copy(post["WTCOL"]); patch.gridcell = ones(Int, np)
        vs = CLM.CNVegStateData(); CLM.cnveg_state_init!(vs, np, nc)
        vs.offset_flag_patch .= post["OFFSET_FLAG"]; vs.offset_counter_patch .= post["OFFSET_COUNTER"]
        cr = CLM.CropData(); CLM.crop_init!(cr, np); cr.harvest_count .= post["HARVEST_COUNT"]
        cs = CLM.CNVegCarbonStateData(); CLM.cnveg_carbon_state_init!(cs, np, nc, 1)
        ns = CLM.CNVegNitrogenStateData(); CLM.cnveg_nitrogen_state_init!(ns, np, nc, 1)
        cf = CLM.CNVegCarbonFluxData(); CLM.cnveg_carbon_flux_init!(cf, np, nc, 1)
        nf = CLM.CNVegNitrogenFluxData(); CLM.cnveg_nitrogen_flux_init!(nf, np, nc, 1)
        for a in (cf.leafc_to_litter_patch, cf.frootc_to_litter_patch, cf.livestemc_to_litter_patch,
                  cf.prev_leafc_to_litter_patch, cf.prev_frootc_to_litter_patch, cf.leafc_to_litter_fun_patch,
                  nf.leafn_to_litter_patch, nf.leafn_to_retransn_patch, nf.frootn_to_litter_patch,
                  nf.livestemn_to_litter_patch)
            a .= 0.0
        end
        for a in (cf.repr_grainc_to_food_thisyr_patch, cf.repr_grainc_to_seed_thisyr_patch,
                  nf.repr_grainn_to_food_thisyr_patch, nf.repr_grainn_to_seed_thisyr_patch)
            a .= 0.0
        end
        cs.leafc_patch .= pre["LEAFC"]; cs.frootc_patch .= pre["FROOTC"]; cs.livestemc_patch .= pre["LIVESTEMC"]
        cs.reproductivec_patch[:, end] .= pre["REPRC"]; cs.cropseedc_deficit_patch .= pre["CROPSEEDC_DEFICIT"]
        ns.leafn_patch .= pre["LEAFN"]; ns.frootn_patch .= pre["FROOTN"]; ns.livestemn_patch .= pre["LIVESTEMN"]
        ns.reproductiven_patch[:, end] .= pre["REPRN"]; ns.cropseedn_deficit_patch .= pre["CROPSEEDN_DEFICIT"]
        cf.cpool_to_leafc_patch .= pre["CPOOL_TO_LEAFC"]; cf.cpool_to_livestemc_patch .= pre["CPOOL_TO_LIVESTEMC"]
        cf.cpool_to_frootc_patch .= pre["CPOOL_TO_FROOTC"]; cf.cpool_to_reproductivec_patch[:, end] .= pre["CPOOL_TO_REPRC"]
        nf.npool_to_reproductiven_patch[:, end] .= pre["NPOOL_TO_REPRN"]
        mask = falses(np); mask[hv] .= true
        CLM.cn_offset_litterfall!(pstate, mask, pf, vs, cs, ns, cf, nf, cr, patch;
                                  mask_pcropp = mask, use_fun = false, CNratio_floating = true)
        jl = Dict("repr_grainc_to_food" => cf.repr_grainc_to_food_patch[:, end],
                  "repr_grainc_to_seed" => cf.repr_grainc_to_seed_patch[:, end],
                  "repr_grainn_to_food" => nf.repr_grainn_to_food_patch[:, end],
                  "repr_grainn_to_seed" => nf.repr_grainn_to_seed_patch[:, end],
                  "leafc_to_litter" => cf.leafc_to_litter_patch, "frootc_to_litter" => cf.frootc_to_litter_patch,
                  "livestemc_to_litter" => cf.livestemc_to_litter_patch, "livestemn_to_litter" => nf.livestemn_to_litter_patch,
                  "leafn_to_litter" => nf.leafn_to_litter_patch, "leafn_to_retransn" => nf.leafn_to_retransn_patch,
                  "frootn_to_litter" => nf.frootn_to_litter_patch, "leafc_to_biofuelc" => cf.leafc_to_biofuelc_patch)
        for p in hv
            push!(events, "n=$n ymd=$(post["ymd"]) itype=$(it[p]) grainc_to_food=$(post["REPR_GRAINC_TO_FOOD"][p])")
            for k in keys_h
                _acc!(stats[k], jl[k][p], post[uppercase(k)][p], "n=$n p=$p itype=$(it[p])"; floor = 1e-15)
            end
        end
    end
    return (; stats, events, available = true)
end

# ---------------------------------------------------------------------------
# (5) n_soyfix! against the in-model CNSoyfix sweep (every case of every step)
# ---------------------------------------------------------------------------
function ct_check_soyfix(dir)
    stat_p = CTStat(); stat_c = CTStat(); nz = 0; ncase_tot = 0
    branches = Set{String}()
    for f in readdir(dir)
        m = match(r"^cropdump_soyfix_sweep_n(\d+)\.nc$", f); m === nothing && continue
        ds = NCDataset(joinpath(dir, f))
        it = Array(ds["ITYPE"][:]); pcol = Array(ds["PCOLUMN"][:]); wt = Array(ds["WTCOL"][:])
        live = Array(ds["CROPLIVE"][:]) .== 1
        pnd = Array(ds["PLANT_NDEMAND"][:]); gddm = Array(ds["GDDMATURITY"][:])
        FPG = Array(ds["FPG"][:, :]); SMN = Array(ds["SMINN"][:, :]); WF = Array(ds["WF"][:, :])
        HUI = Array(ds["HUI"][:, :]); SFX = Array(ds["SOYFIXN"][:, :]); S2S = Array(ds["SOYFIXN_TO_SMINN"][:, :])
        close(ds)
        np = length(it); nc = size(FPG, 1)
        for k in axes(FPG, 2)
            ncase_tot += 1
            patch = CLM.PatchData(); patch.itype = copy(it); patch.column = copy(pcol)
            patch.wtcol = copy(wt); patch.gridcell = ones(Int, np)
            cr = CLM.CropData(); CLM.crop_init!(cr, np); cr.croplive_patch .= live
            cr.hui_patch .= HUI[:, k]
            vs = CLM.CNVegStateData(); CLM.cnveg_state_init!(vs, np, nc); vs.gddmaturity_patch .= gddm
            nf = CLM.CNVegNitrogenFluxData(); CLM.cnveg_nitrogen_flux_init!(nf, np, nc, 1)
            nf.plant_ndemand_patch .= pnd
            sbs = CLM.SoilBiogeochemStateData(); CLM.soil_bgc_state_init!(sbs, nc, np, 1, 1); sbs.fpg_col .= FPG[:, k]
            sbns = CLM.SoilBiogeochemNitrogenStateData(); CLM.soil_bgc_nitrogen_state_init!(sbns, nc, 1, 1, 1)
            sbns.sminn_col .= SMN[:, k]
            sbnf = CLM.SoilBiogeochemNitrogenFluxData(); CLM.soil_bgc_nitrogen_flux_init!(sbnf, nc, 1, 1, 1)
            wd = CLM.WaterDiagnosticBulkData{Float64}(); wd.wf_col = copy(WF[:, k])
            # CTSM passes filter_bgc_soilc / filter_bgc_vegp: every soil column and veg patch
            CLM.n_soyfix!(sbnf, nf, sbs, sbns, vs, cr, wd;
                mask_soilc = trues(nc), bounds = 1:nc, mask_soilp = trues(np), bounds_p = 1:np,
                patch = patch, ntmp_soybean = CLM.ntmp_soybean, nirrig_tmp_soybean = CLM.nirrig_tmp_soybean,
                ntrp_soybean = CLM.ntrp_soybean, nirrig_trp_soybean = CLM.nirrig_trp_soybean)
            for p in 1:np
                _acc!(stat_p, nf.soyfixn_patch[p], SFX[p, k], "$(f) case=$k p=$p"; floor = 1e-15)
                if it[p] in (23, 24) && live[p]
                    SFX[p, k] > 0 && (nz += 1)
                    gf = HUI[p, k] / gddm[p]; c = pcol[p]
                    push!(branches, string(FPG[c, k] < 1 ? "fpg<1" : "fpg>=1", "/",
                        SMN[c, k] > 30 ? "smn>30" : SMN[c, k] > 10 ? "10<smn<=30" : "smn<=10", "/",
                        gf <= 0.15 ? "g<=.15" : gf <= 0.30 ? "g<=.30" : gf <= 0.55 ? "g<=.55" : gf <= 0.75 ? "g<=.75" : "g>.75", "/",
                        WF[c, k] / 0.85 > 1 ? "fxw>1" : "fxw<1"))
                end
            end
            for c in 1:nc
                _acc!(stat_c, sbnf.soyfixn_to_sminn_col[c], S2S[c, k], "$(f) case=$k c=$c"; floor = 1e-15)
            end
        end
    end
    return (; stat_p, stat_c, nz, ncase_tot, branches)
end

# ---------------------------------------------------------------------------
# (6) vernalization!  inputs: 'after_cnalloc'(n) (vf/cumvd/hdidx before
#     CNPhenology, which in the non-FUN driver runs after allocation), outputs:
#     'end_noleach'(n). Gate (CropPhenology): croplive, winter wheat (21/22),
#     t_ref2m_min < 1e30, vf /= 1 -- checked on both sides: where the gate is
#     closed the Fortran state must not move.
# ---------------------------------------------------------------------------
function ct_check_vernalization(dir)
    stats = Dict(k => CTStat() for k in ("vf", "cumvd", "hdidx"))
    ncalls = 0; nclosed = 0; nharvest_forced = 0
    pstate = CLM.PhenologyState()
    for n in ct_steps(dir, "end_noleach")
        pre = ct_load(dir, "after_cnalloc", n); post = ct_load(dir, "end_noleach", n)
        pre === nothing && continue
        it = pre["ITYPE"]; np = length(it); nc = length(pre["FPG"])
        for p in 1:np
            it[p] in (21, 22) || continue
            pre["CROPLIVE"][p] == 1 || continue
            gate = pre["T_REF2M_MIN"][p] < 1e30 && pre["VF"][p] != 1.0
            if !gate
                nclosed += 1
                _acc!(stats["vf"], pre["VF"][p], post["VF"][p], "n=$n closed")
                continue
            end
            ncalls += 1
            patch = CLM.PatchData(); patch.itype = copy(it); patch.column = copy(pre["PCOLUMN"])
            cs_ = CLM.CanopyStateData{Float64}(); CLM.canopystate_init!(cs_, np); cs_.tlai_patch .= pre["TLAI"]
            tm = CLM.TemperatureData{Float64}()
            tm.t_ref2m_patch = copy(pre["T_REF2M"]); tm.t_ref2m_min_patch = copy(pre["T_REF2M_MIN"])
            tm.t_ref2m_max_patch = copy(pre["T_REF2M_MAX"])
            wd = CLM.WaterDiagnosticBulkData{Float64}(); wd.snow_depth_col = copy(pre["SNOW_DEPTH"])
            vs = CLM.CNVegStateData(); CLM.cnveg_state_init!(vs, np, nc)
            vs.cumvd_patch .= pre["CUMVD"]; vs.hdidx_patch .= pre["HDIDX"]
            cr = CLM.CropData(); CLM.crop_init!(cr, np); cr.vf_patch .= pre["VF"]
            fh = CLM.vernalization!(p, pstate, cs_, tm, wd, vs, cr, patch)
            fh && (nharvest_forced += 1)
            tag = "n=$n ymd=$(pre["ymd"]) tod=$(pre["tod"]) itype=$(it[p])"
            _acc!(stats["vf"], cr.vf_patch[p], post["VF"][p], tag)
            _acc!(stats["cumvd"], vs.cumvd_patch[p], post["CUMVD"][p], tag)
            _acc!(stats["hdidx"], vs.hdidx_patch[p], post["HDIDX"][p], tag)
        end
    end
    return (; stats, ncalls, nclosed, nharvest_forced)
end

function _report(name, s::CTStat)
    @printf("  %-22s n=%-6d max|abs|=%-11.3e max rel=%-11.3e %s\n", name, s.n, s.maxabs, s.maxrel,
            s.maxrel > 1e-10 ? "worst: " * s.worst : "")
end

function ct_main(refA, refSoy)
    println("crop-tail parity  refA=$refA  refSoy=$refSoy")
    pc = ct_pftcon()
    a = ct_check_allocation(refA; pfta = ct_pfta(pc))
    println("(1) calc_crop_allocation_fractions!  steps by phase: ", a.phase_count)
    for k in ("aleaf", "astem", "aroot", "arepr", "aleafi", "astemi"); _report(k, a.stats[k]); end
    println("    chained (Julia-carried aleaf/astem through grain fill):")
    for k in ("aleaf", "astem", "arepr"); _report(k, a.chained[k]); end
    b = ct_check_cnalloc(refA; pc = pc)
    println("(2/3) calc_allometry! + calc_plant_cn_alloc!  (grain-filling patch-steps: $(b.ngrain))")
    for k in ("c_allometry", "n_allometry", "plant_calloc", "cpool_to_reprc", "cpool_to_reprc_stor",
              "cpool_to_leafc", "cpool_to_livestemc", "cpool_to_frootc"); _report(k, b.stats[k]); end
    hv = ct_check_harvest(refA; pc = pc)
    if hv.available
        println("(4) cn_offset_litterfall! crop-harvest branch  ($(length(hv.events)) harvested patches)")
        for k in sort(collect(keys(hv.stats))); _report(k, hv.stats[k]); end
    else
        println("(4) skipped: dumps predate the harvest-flux fields")
    end
    if isdir(refSoy)
        s = ct_check_soyfix(refSoy)
        println("(5) n_soyfix! vs in-model CNSoyfix sweep: $(s.ncase_tot) cases, $(s.nz) soybean patch-cases with SOYFIXN>0, $(length(s.branches)) branch combos")
        _report("soyfixn", s.stat_p); _report("soyfixn_to_sminn", s.stat_c)
    end
    v = ct_check_vernalization(refA)
    println("(6) vernalization!  calls=$(v.ncalls)  gate-closed steps=$(v.nclosed)  forced harvests=$(v.nharvest_forced)")
    for k in ("vf", "cumvd", "hdidx"); _report(k, v.stats[k]); end
end

if get(ENV, "CROP_TAILS_NO_MAIN", "0") != "1" && abspath(PROGRAM_FILE) == @__FILE__
    refA = length(ARGS) >= 1 ? ARGS[1] : joinpath(_SD, "runs", "crop_tails_refA")
    refS = length(ARGS) >= 2 ? ARGS[2] : joinpath(_SD, "runs", "crop_tails_refB_soyfix")
    ct_main(refA, refS)
end
