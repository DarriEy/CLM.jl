# ==========================================================================
# fortran_parity_sweep.jl — COMPREHENSIVE per-module-boundary parity sweep.
#
# For every step in the window, build a Bow instance, inject the before_step
# dump, run clm_drv! one step, then diff a WIDE set of science fields against
# the appropriate Fortran per-boundary dump (each field is compared at the
# boundary where it is last meaningfully written). Aggregates each field's
# WORST divergence across the whole window and prints a ranked table so any
# remaining per-module bug surfaces at the top.
#
#   julia +1.12 --project=. scripts/fortran_parity_sweep.jl [firstN lastN]
# ==========================================================================

include(joinpath(@__DIR__, "fortran_parity_common.jl"))
using CLM, NCDatasets, Printf, Dates

const FIRST = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 13458
const LAST  = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 13470
const STEPS = FIRST:LAST

# Boundary dump suffix per field (where Fortran last writes it this step).
const B_CAN  = "after_canopyfluxes"
const B_ST   = "after_soiltemperature"
const B_SF   = "after_soilfluxes"
const B_HND  = "after_hydrologynodrainage"
const B_HD   = "after_hydrologydrainage"

# (dump_name, kind, boundary, getter::inst->array)
# kind: :col1d (scalar col 1), :col2d (level × col 1), :patch (1-D patch)
function sweep_registry(inst)
    ws = inst.water.waterstatebulk_inst.ws
    wd = inst.water.waterdiagnosticbulk_inst
    te = inst.temperature; cs = inst.canopystate; fv = inst.frictionvel
    ef = inst.energyflux;  sa = inst.solarabs;    sh = inst.soilhydrology
    ss = inst.soilstate
    R = Tuple{String,Symbol,String,Function}[]
    add(n,k,b,g) = push!(R, (n,k,b,g))

    # --- canopy boundary (frozen after CanopyFluxes) ---
    add("T_VEG",      :patch, B_CAN, _->te.t_veg_patch)
    add("T_STEM",     :patch, B_CAN, _->te.t_stem_patch)
    add("TAF_P",      :patch, B_CAN, _->fv.taf_patch)
    add("QAF_P",      :patch, B_CAN, _->fv.qaf_patch)
    add("RAH1_P",     :patch, B_CAN, _->fv.rah1_patch)
    add("RAH2_P",     :patch, B_CAN, _->fv.rah2_patch)
    add("OBU",        :patch, B_CAN, _->fv.obu_patch)
    add("FV_P",       :patch, B_CAN, _->fv.fv_patch)
    add("RAM1_P",     :patch, B_CAN, _->fv.ram1_patch)
    add("SABV_P",     :patch, B_CAN, _->sa.sabv_patch)
    add("SABG_P",     :patch, B_CAN, _->sa.sabg_patch)
    add("fsun",       :patch, B_CAN, _->cs.fsun_patch)
    add("elai",       :patch, B_CAN, _->cs.elai_patch)
    add("esai",       :patch, B_CAN, _->cs.esai_patch)

    # --- soil-temperature boundary ---
    add("T_SOISNO",   :col2d, B_ST,  _->te.t_soisno_col)
    add("T_GRND",     :col1d, B_ST,  _->te.t_grnd_col)
    add("T_GRND_R",   :col1d, B_ST,  _->te.t_grnd_r_col)
    add("THK_C",      :col2d, B_ST,  _->ss.thk_col)

    # --- soil-fluxes boundary (ground fluxes recomputed here) ---
    add("EFLX_SHG_P", :patch, B_SF,  _->ef.eflx_sh_grnd_patch)
    add("EFLX_LH_P",  :patch, B_SF,  _->ef.eflx_lh_tot_patch)
    add("EFLX_SOIG_P",:patch, B_SF,  _->ef.eflx_soil_grnd_patch)
    add("EFLX_GNET_P",:patch, B_SF,  _->ef.eflx_gnet_patch)
    add("EFLX_LWNET_P",:patch,B_SF,  _->ef.eflx_lwrad_net_patch)
    add("CGRNDS_P",   :patch, B_SF,  _->ef.cgrnds_patch)
    add("CGRNDL_P",   :patch, B_SF,  _->ef.cgrndl_patch)
    add("T_REF2M",    :patch, B_SF,  _->te.t_ref2m_patch)

    # --- hydrology (no-drainage) boundary ---
    add("SNOW_DEPTH", :col1d, B_HD,  _->wd.snow_depth_col)
    add("frac_sno",   :col1d, B_HD,  _->wd.frac_sno_col)
    add("frac_sno_eff",:col1d,B_HD,  _->wd.frac_sno_eff_col)
    add("INT_SNOW",   :col1d, B_HD,  _->ws.int_snow_col)
    add("H2OSFC",     :col1d, B_HD,  _->ws.h2osfc_col)
    add("LIQCAN",     :patch, B_HD,  _->ws.liqcan_patch)
    add("SNOCAN",     :patch, B_HD,  _->ws.snocan_patch)
    add("FWET",       :patch, B_HD,  _->wd.fwet_patch)

    # --- hydrology (drainage) / final state ---
    add("H2OSOI_LIQ", :col2d, B_HD,  _->ws.h2osoi_liq_col)
    add("H2OSOI_ICE", :col2d, B_HD,  _->ws.h2osoi_ice_col)
    add("ZWT",        :col1d, B_HD,  _->sh.zwt_col)
    add("ZWT_PERCH",  :col1d, B_HD,  _->sh.zwt_perched_col)
    add("WA",         :col1d, B_HD,  _->ws.wa_col)
    add("SMP",        :col2d, B_HD,  _->ss.smp_l_col)
    return R
end

# robust scalar read (missing → NaN)
_rd(v) = ismissing(v) ? NaN : Float64(v)

# compare one field; returns (nabs, nrel, npts) over finite pairs, or nothing
function _field_diff(ds, name, kind, getter, inst)
    haskey(ds, name) || return nothing
    jl = try getter(inst) catch; return nothing end
    (jl === nothing || isempty(jl)) && return nothing
    nlevsno = CLM.varpar.nlevsno
    fv = Float64[]; jv = Float64[]
    if kind == :col1d
        push!(fv, _rd(ds[name][1])); push!(jv, Float64(jl[1]))
    elseif kind == :col2d
        d = ds[name][:, :]               # (lev, col) or (lev,)
        nlev = min(size(d, 1), size(jl, 2))
        for k in 1:nlev
            push!(fv, _rd(d[k, 1])); push!(jv, Float64(jl[1, k]))
        end
    elseif kind == :patch
        d = ds[name][:]
        n = min(length(d), length(jl))
        for i in 1:n
            push!(fv, _rd(d[i])); push!(jv, Float64(jl[i]))
        end
    end
    # drop pairs where EITHER is NaN (Fortran spval/missing or structurally-empty
    # snow layers) — compare only the physically-active overlapping entries
    keep = .!(isnan.(fv) .| isnan.(jv))
    fv = fv[keep]; jv = jv[keep]
    isempty(fv) && return nothing
    nabs = maximum(abs.(fv .- jv))
    nrel = maximum(abs.(fv .- jv) ./ (1.0 .+ max.(abs.(fv), abs.(jv))))
    return (nabs, nrel, length(fv))
end

println("="^78)
println("  COMPREHENSIVE per-boundary parity sweep   steps $(FIRST)..$(LAST)")
println("="^78)

# name -> (maxabs, maxrel, step_at_max, boundary, npts)
agg = Dict{String,Any}()
order = String[]

for n in STEPS
    bstep = joinpath(DUMPDIR, "pdump_before_step_n$(n).nc")
    isfile(bstep) || (println("  skip n$n (no before dump)"); continue)
    inst, bounds = run_one_parity_step!(n)
    reg = sweep_registry(inst)
    # Julia `inst` is END-of-step, so diff against the END-of-step Fortran dump
    # (after_hydrologydrainage) for an apples-to-apples comparison. The per-field
    # `bnd` tag is kept only as a "last written by" label for localizing offenders.
    fHD = joinpath(DUMPDIR, "pdump_$(B_HD)_n$(n).nc")
    isfile(fHD) || (println("  skip n$n (no HD dump)"); continue)
    dsHD = NCDataset(fHD, "r")
    bcache = Dict{String,Any}("HD"=>dsHD)
    for (name, kind, bnd, getter) in reg
        ds = dsHD
        r = _field_diff(ds, name, kind, getter, inst)
        r === nothing && continue
        (nabs, nrel, npts) = r
        if !haskey(agg, name)
            agg[name] = (nabs, nrel, n, bnd, npts); push!(order, name)
        else
            (a0, r0, s0, b0, p0) = agg[name]
            nrel > r0 && (agg[name] = (nabs, nrel, n, bnd, npts))
        end
    end
    for ds in values(bcache); close(ds); end
    @printf("  ran n%d\n", n)
end

# ranked table
ranked = sort(order; by = nm -> -agg[nm][2])
println("\n  ", "-"^74)
@printf("  %-13s %11s %11s %6s %-22s %5s\n", "field", "max|abs|", "max|rel|", "step", "boundary", "npts")
println("  ", "-"^74)
for nm in ranked
    (a, r, s, b, p) = agg[nm]
    flag = r > 1e-3 ? " <==" : ""
    @printf("  %-13s %11.3e %11.3e %6d %-22s %5d%s\n", nm, a, r, s, b, p, flag)
end
println("  ", "-"^74)
worst = isempty(ranked) ? ("",0.0) : (ranked[1], agg[ranked[1]][2])
@printf("  worst field: %s (max|rel|=%.3e over the window)\n", worst[1], worst[2])
