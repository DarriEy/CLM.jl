# ==========================================================================
# fortran_parity_common.jl — shared infrastructure for Julia ↔ Fortran
# per-module / per-timestep parity.
#
# The Fortran side (instrumented cesm.exe via SourceMods restFile_write_dump)
# emits per-boundary snapshots `pdump_<boundary>_n<nstep>.nc` in CLM5
# restart format. This module:
#   - sets up a Bow-at-Banff single-column Julia CLM instance,
#   - injects a Fortran dump as the shared initial condition
#     (reusing CLM.read_fortran_restart!),
#   - compares the live Julia `inst` state to a Fortran dump field-by-field,
#     handling the snow-pad / levtot / patch-remap / fill→NaN conventions.
#
# Tolerances (tiered, per the parity plan): CPU vs Fortran from identical
# inputs ~1e-10 relative is the Float64 cross-compiler floor.
# ==========================================================================

using CLM
using NCDatasets
using Printf
using Dates

const _read_fortran_restart! = getfield(CLM, Symbol("read_fortran_restart!"))

# ---- Bow-at-Banff paths (the instrumented Fortran run's config) ----
const BOW_BASE = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped"
const BOW_CAL  = joinpath(BOW_BASE, "optimization/CLM/dds_run_1/final_evaluation/settings/CLM/parameters")
const FSURDAT  = joinpath(BOW_CAL, "surfdata_clm.nc")
const FPARAM   = joinpath(BOW_CAL, "clm5_params.nc")
const FFORCING = joinpath(BOW_BASE, "data/forcing/CLM_input/clmforc.2003.nc")
const FSNOWOPT = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_c013122.nc"
const FSNOWAGE = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc"
const DUMPDIR  = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_parity_run"

# Calibrated namelist values from the Fortran run (user_nl_clm / verify_vs_fortran.jl)
const BASEFLOW_SCALAR = 0.0022119554
const INT_SNOW_MAX    = 3113.2227

# ---- NaN-aware diffs (reused pattern from gpu_validate_*) ----
reldiff(a, b) = begin
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(A[i] - B[i]) / (1.0 + max(abs(A[i]), abs(B[i]))))
    end
    m
end
absdiff(a, b) = begin
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(A[i] - B[i]))
    end
    m
end

# Read a Fortran restart variable, converting missing/fill to NaN.
function _dumpvar(ds, name)
    haskey(ds, name) || return nothing
    raw = ds[name][:]
    return [ismissing(v) ? NaN : Float64(v) for v in raw]
end

"""
    build_bow_inst(; dtime=3600, use_aquifer_layer=false)

Initialize a single-column Bow-at-Banff Julia CLM instance matching the
instrumented Fortran run's configuration, with the runtime params wired in.
Returns (inst, bounds, filt, tm).
"""
function build_bow_inst(; dtime::Int=3600, use_aquifer_layer::Bool=false,
                          start_date::DateTime=DateTime(2003,1,1),
                          use_cn::Bool=false)
    (inst, bounds, filt, tm) = CLM.clm_initialize!(;
        fsurdat=FSURDAT, paramfile=FPARAM,
        start_date=start_date, dtime=dtime, use_cn=use_cn,
        use_bedrock=true, use_aquifer_layer=use_aquifer_layer,
        h2osfcflag=0, fsnowoptics=FSNOWOPT, fsnowaging=FSNOWAGE,
        int_snow_max=INT_SNOW_MAX)

    # The Bow run's lnd_in has use_fun=.true., use_flexiblecn=.true. (CLM5 default).
    # Enable them here — the parity harness runs from a warm Fortran restart, so the
    # canopy/availc are finite and FUN's carbon-cost N uptake is well-posed.
    if use_cn
        inst.bgc_vegetation.config.use_fun = true
        inst.bgc_vegetation.config.use_flexiblecn = true
        inst.bgc_vegetation.driver_config.use_fun = true          # set both: the
        inst.bgc_vegetation.driver_config.use_flexiblecn = true   # fc→dc sync already ran in init
    end

    # atm2lnd downscaling to match Fortran lnd_in defaults
    CLM.atm2lnd_read_namelist!(inst.atm2lnd;
        repartition_rain_snow=true, lapse_rate=0.006, lapse_rate_longwave=0.032,
        precip_repartition_nonglc_all_snow_t=0.0, precip_repartition_nonglc_all_rain_t=2.0,
        precip_repartition_glc_all_snow_t=-2.0, precip_repartition_glc_all_rain_t=0.0)
    CLM.init_soil_hydrology_config(baseflow_scalar=BASEFLOW_SCALAR)

    # Replicate clm_run.jl's runtime param wiring (n_melt, accum_factor, snow
    # density, grain radius, roughness, aging) — without this the snow-cover
    # n_melt is at its init default and int_snow diverges badly.
    nc = bounds.endc
    ds_p = NCDataset(FPARAM, "r")
    scf = inst.scf_method
    if haskey(ds_p, "n_melt_coef")
        nmc = Float64(ds_p["n_melt_coef"][1])
        for c in 1:nc
            topo_std = max(10.0, inst.column.topo_std[c])
            scf.n_melt[c] = nmc / topo_std
        end
    end
    haskey(ds_p, "accum_factor") && (scf.accum_factor = Float64(ds_p["accum_factor"][1]))
    haskey(ds_p, "SNOW_DENSITY_MAX") && (CLM.snowhydrology_params.rho_max = Float64(ds_p["SNOW_DENSITY_MAX"][1]))
    haskey(ds_p, "SNOW_DENSITY_MIN") && (CLM.snowhydrology_params.rho_min = Float64(ds_p["SNOW_DENSITY_MIN"][1]))
    haskey(ds_p, "fresh_snw_rds_max") && (CLM.snowhydrology_params.snw_rds_min = Float64(ds_p["fresh_snw_rds_max"][1]))
    haskey(ds_p, "SNO_Z0MV") && (inst.frictionvel.zsno = Float64(ds_p["SNO_Z0MV"][1]))
    if haskey(ds_p, "snw_aging_bst")
        CLM.snicar_params.xdrdt = Float64(ds_p["snw_aging_bst"][1])
    elseif haskey(ds_p, "xdrdt")
        CLM.snicar_params.xdrdt = Float64(ds_p["xdrdt"][1])
    end
    if haskey(ds_p, "pc")
        pc_val = Float64(ds_p["pc"][1])
        pc_val > 0 && (for c in 1:nc; inst.soilhydrology.hkdepth_col[c] = 1.0 / pc_val; end)
    end
    close(ds_p)

    return (inst, bounds, filt, tm)
end

# Field registry: (fortran_name, kind, julia_getter)
#   kind = :col1d  -> scalar column var, julia_getter(inst) returns inst-array[1]
#   kind = :col2d  -> levtot column var, julia_getter(inst) returns the (col,lev) matrix
#   kind = :patch  -> patch var, julia_getter(inst) returns the patch vector
function _parity_registry(inst)
    ws = inst.water.waterstatebulk_inst.ws
    wd = inst.water.waterdiagnosticbulk_inst
    [
        ("T_SOISNO",   :col2d, () -> inst.temperature.t_soisno_col),
        ("H2OSOI_LIQ", :col2d, () -> ws.h2osoi_liq_col),
        ("H2OSOI_ICE", :col2d, () -> ws.h2osoi_ice_col),
        ("T_GRND",     :col1d, () -> inst.temperature.t_grnd_col),
        ("WA",         :col1d, () -> ws.wa_col),
        ("H2OSFC",     :col1d, () -> ws.h2osfc_col),
        ("ZWT",        :col1d, () -> inst.soilhydrology.zwt_col),
        ("ZWT_PERCH",  :col1d, () -> inst.soilhydrology.zwt_perched_col),
        ("INT_SNOW",   :col1d, () -> inst.water.waterstatebulk_inst.int_snow_col),
        ("SNOW_DEPTH", :col1d, () -> wd.snow_depth_col),
        ("frac_sno",   :col1d, () -> wd.frac_sno_col),
        ("T_VEG",      :patch, () -> inst.temperature.t_veg_patch),
        ("elai",       :patch, () -> inst.canopystate.elai_patch),
        ("tlai",       :patch, () -> inst.canopystate.tlai_patch),
        ("LIQCAN",     :patch, () -> ws.liqcan_patch),
        ("SNOCAN",     :patch, () -> ws.snocan_patch),
    ]
end

"""
    compare_inst_to_dump(inst, dumpfile; label="", tol=1e-9)

Diff the live Julia `inst` prognostic state against a Fortran pdump file.
Prints per-field max-abs / max-rel diffs and returns a Vector of
(name, nabs, nrel, npts) plus the global max relative diff.
"""
function compare_inst_to_dump(inst, dumpfile; label::String="", tol::Float64=1e-9)
    ds = NCDataset(dumpfile, "r")
    nlevsno = CLM.varpar.nlevsno
    reg = _parity_registry(inst)
    results = Tuple{String,Float64,Float64,Int}[]
    gmax = 0.0
    @printf("  %-14s %12s %12s  %s\n", "field", "max|abs|", "max|rel|", "status")
    @printf("  %s\n", "-"^52)
    for (fname, kind, getter) in reg
        data = _dumpvar(ds, fname)
        data === nothing && continue
        jl = getter()
        nabs = NaN; nrel = NaN; npts = 0
        if kind == :col1d
            f = data[1]; j = Float64(jl[1])
            (isnan(f) && isnan(j)) && continue
            nabs = abs(f - j); nrel = abs(f - j) / (1.0 + max(abs(f), abs(j))); npts = 1
        elseif kind == :col2d
            # Fortran levtot (snow 1:nlevsno then soil); compare overlapping levels of col 1
            nlev = min(length(data), size(jl, 2))
            fv = Float64[]; jv = Float64[]
            for k in 1:nlev
                push!(fv, data[k]); push!(jv, Float64(jl[1, k]))
            end
            nabs = absdiff(fv, jv); nrel = reldiff(fv, jv); npts = nlev
        elseif kind == :patch
            # compare overlapping patches (Fortran np ≤ Julia np); naive index map
            npf = length(data); npj = length(jl)
            n = min(npf, npj)
            fv = Float64[data[i] for i in 1:n]; jv = Float64[Float64(jl[i]) for i in 1:n]
            nabs = absdiff(fv, jv); nrel = reldiff(fv, jv); npts = n
        end
        isnan(nrel) && continue
        gmax = max(gmax, nrel)
        status = nrel <= tol ? "PASS" : "DIFF"
        @printf("  %-14s %12.3e %12.3e  %s (%d pts)\n", fname, nabs, nrel, status, npts)
        push!(results, (fname, nabs, nrel, npts))
    end
    close(ds)
    @printf("  %s\n", "-"^52)
    @printf("  %s global max|rel| = %.3e  (tol %.0e)\n", label, gmax, tol)
    return results, gmax
end

"""
    inject_dump!(inst, bounds, dumpfile)

Inject a Fortran dump as the Julia initial condition (reuses
CLM.read_fortran_restart!).
"""
inject_dump!(inst, bounds, dumpfile) = _read_fortran_restart!(dumpfile, inst, bounds)

"""
    run_one_parity_step!(nstep) -> (inst, bounds)

Build a Bow instance, inject the `before_step_n<nstep>` dump as the shared IC,
read/downscale forcing at the step START, seed daylength, and run clm_drv! for
exactly one timestep. Mirrors the clm_run.jl loop body (and the single-step
setup in fortran_parity_validate.jl). Returns the post-step instance for diffing
against the Fortran per-boundary dumps.
"""
function run_one_parity_step!(nstep::Int; use_cn::Bool=false, dumpdir::String=DUMPDIR,
                              step_date::DateTime=DateTime(2003, 1, 1) + Hour(nstep - 8761),
                              forcing_file::String=FFORCING)
    (inst, bounds, filt, tm) = build_bow_inst(; dtime=3600, start_date=step_date, use_cn=use_cn)
    inject_dump!(inst, bounds, joinpath(dumpdir, "pdump_before_step_n$(nstep).nc"))

    config  = CLM.CLMDriverConfig(use_cn=use_cn, use_aquifer_layer=false)
    filt_ia = CLM.clump_filter_inactive_and_active
    ng, nc, np = bounds.endg, bounds.endc, bounds.endp

    fr = CLM.ForcingReader()
    CLM.forcing_reader_init!(fr, forcing_file)

    topo_file = replace(forcing_file, r"clmforc\.[^/]*\.nc$" => "topo_forcing.nc")
    if isfile(topo_file)
        ds_topo = NCDataset(topo_file, "r")
        if haskey(ds_topo, "TOPO")
            ft = Float64(ds_topo["TOPO"][1])
            for g in 1:ng; inst.atm2lnd.forc_topo_grc[g] = ft; end
            for c in 1:nc; inst.topo.topo_col[c] = ft; end
        end
        close(ds_topo)
    end

    step_start  = tm.current_date
    calday      = CLM.get_curr_calday(tm)
    (declin, _) = CLM.compute_orbital(calday)
    obliqr      = CLM.ORB_OBLIQR_DEFAULT
    nextsw_cday = calday + 3600.0 / CLM.SECSPDAY
    (declinm1, _) = CLM.compute_orbital(calday - 3600.0 / CLM.SECSPDAY)
    CLM.init_daylength!(inst.gridcell, declin, declinm1, obliqr, 1:bounds.endg)

    CLM.advance_timestep!(tm)
    CLM.read_forcing_step!(fr, inst.atm2lnd, step_start, ng, nc)
    CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
    (yr, mon, d, tod) = CLM.get_curr_date(tm)

    CLM.clm_drv!(config, inst, filt, filt_ia, bounds,
                 true, nextsw_cday, declin, declin, obliqr,
                 false, false, "", false;
                 nstep=tm.nstep, is_first_step=(tm.nstep==1),
                 is_beg_curr_day=CLM.is_beg_curr_day(tm),
                 is_end_curr_day=CLM.is_end_curr_day(tm),
                 is_beg_curr_year=CLM.is_beg_curr_year(tm),
                 dtime=3600.0, mon=mon, day=d, photosyns=inst.photosyns)
    CLM.forcing_reader_close!(fr)
    return inst, bounds
end
