# =============================================================================
# initcold_nan_probe.jl — sweep EVERY numeric field on EVERY state struct in
# CLMInstances after clm_initialize! (+ optionally N clm_drv! steps) and report
# which are NaN on ACTIVE columns / patches / landunits.
#
# Motivation: Fortran's `Init = InitAllocate + InitHistory + InitCold`, where
# InitAllocate NaN-fills and InitCold zeroes/seeds. Several `*_init_cold!` ports
# in CLM.jl were translated but never CALLED, so the NaN-fill survives into the
# run. Two such gaps were live bugs (waterflux -> qflx_evap_tot NaN on urban
# roof/wall columns; aerosol -> SNICAR albedo NaN -> glacier band solve NaN).
# This probe finds the rest, and stays in-tree as the class-level guard.
#
# The sweep is REFLECTIVE (fieldnames over the instance tree), so a newly added
# state field is covered automatically — no list to maintain.
#
# Usage:
#   julia +1.12 --project=. scripts/initcold_nan_probe.jl            # all domains
#   julia +1.12 --project=. scripts/initcold_nan_probe.jl bow_cn     # one domain
# =============================================================================
using CLM, NCDatasets, Dates, Printf

const BOW_CAL = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/" *
                "domain_Bow_at_Banff_lumped/settings/CLM/parameters"
const BOW_FS  = joinpath(BOW_CAL, "surfdata_clm.nc")
const BOW_FP  = joinpath(BOW_CAL, "clm5_params.nc")
const GLAC_FS = joinpath(@__DIR__, "..", "test_inputs", "glacier", "surfdata_glacier100.nc")
const LAKE_FS = joinpath(@__DIR__, "..", "test_inputs", "lake", "surfdata_lake100.nc")
const URB_FS  = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/installs/clm/" *
                "python/ctsm/test/testinputs/" *
                "surfdata_1x1_mexicocityMEX_hist_16pfts_CMIP6_2000_c231103.nc"
const SNOWOPT = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_c013122.nc"
const SNOWAGE = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc"

# -----------------------------------------------------------------------------
# Fields Fortran DELIBERATELY leaves at the NaN/spval missing-value flag on
# columns outside the routine's filter — zeroing these would be WRONG (a real
# 0.0 flux is not the same as "this column has no glacier-ice mass balance").
# Keyed by field name; matched on any struct.
# -----------------------------------------------------------------------------
const DELIBERATE_NAN = Set([
    # qflx_glcice* are only set inside the do_smb filter (glacier surface mass
    # balance); outside it Fortran leaves them at spval as a missing-value flag
    # and the history buffer masks them out. See WaterFluxType.F90::InitCold.
    :qflx_glcice_col, :qflx_glcice_frz_col, :qflx_glcice_melt_col,
    :qflx_glcice_dyn_water_flux_col,
    # CropType::InitCold sets latbaset = NaN unless baset_mapping == LATVARY.
    # It is READ only under LATVARY, so the NaN is the "not applicable" flag.
    :latbaset_patch,
])

# -----------------------------------------------------------------------------
# Reflective sweep
# -----------------------------------------------------------------------------
"""
    active_masks(inst, bounds) -> (act_c, act_p, act_l)

BitVectors of ACTIVE columns / patches / landunits (set_active! in init).
"""
function active_masks(inst, bounds)
    act_c = BitVector(undef, bounds.endc)
    act_p = BitVector(undef, bounds.endp)
    act_l = BitVector(undef, bounds.endl)
    for c in 1:bounds.endc; act_c[c] = inst.column.active[c];   end
    for p in 1:bounds.endp; act_p[p] = inst.patch.active[p];    end
    for l in 1:bounds.endl; act_l[l] = inst.landunit.active[l]; end
    return (act_c, act_p, act_l)
end

# Decide which subgrid level a field's FIRST dimension indexes, from its size.
# Ambiguity (e.g. nc == np) is resolved by the field-name suffix, which the port
# preserves from Fortran (_col / _patch / _lun / _grc).
function level_of(name::Symbol, n::Int, bounds)
    s = String(name)
    endswith(s, "_col")   && return (:col, n == bounds.endc)
    endswith(s, "_patch") && return (:patch, n == bounds.endp)
    endswith(s, "_lun")   && return (:lun, n == bounds.endl)
    endswith(s, "_grc")   && return (:grc, n == bounds.endg)
    # Unsuffixed (e.g. WaterStateData's liqcan_patch is suffixed, but some
    # container fields are not) — fall back to a unique size match.
    matches = Symbol[]
    n == bounds.endc && push!(matches, :col)
    n == bounds.endp && push!(matches, :patch)
    n == bounds.endl && push!(matches, :lun)
    length(matches) == 1 && return (matches[1], true)
    return (:unknown, false)
end

"""
    sweep_struct(obj, structname, bounds, masks) -> Vector{NamedTuple}

Every numeric array field of `obj` whose leading dim is a subgrid level and that
holds a non-finite value on an ACTIVE index of that level.
"""
function sweep_struct(obj, structname::String, bounds, masks)
    (act_c, act_p, act_l) = masks
    out = NamedTuple[]
    obj === nothing && return out
    for fname in fieldnames(typeof(obj))
        v = getfield(obj, fname)
        v isa AbstractArray || continue
        eltype(v) <: AbstractFloat || continue
        isempty(v) && continue
        n = size(v, 1)
        (lvl, ok) = level_of(fname, n, bounds)
        ok || continue
        mask = lvl === :col ? act_c : lvl === :patch ? act_p : lvl === :lun ? act_l : nothing
        mask === nothing && continue
        nbad = 0; firstidx = 0
        for i in 1:n
            mask[i] || continue
            row = ndims(v) == 1 ? (isfinite(v[i]) ? () : (v[i],)) :
                  (any(!isfinite, @view v[i, ntuple(_ -> :, ndims(v) - 1)...]) ? (NaN,) : ())
            if !isempty(row)
                nbad += 1
                firstidx == 0 && (firstidx = i)
            end
        end
        nbad > 0 && push!(out, (struct_name=structname, field=fname, level=lvl,
                                nbad=nbad, nactive=count(mask), firstidx=firstidx,
                                deliberate=(fname in DELIBERATE_NAN)))
    end
    return out
end

"""
    sweep_instances(inst, bounds) -> Vector{NamedTuple}

Recursively sweep the whole CLMInstances tree (one level of nested containers,
which is what `water` / `bgc_vegetation` are).
"""
function sweep_instances(inst, bounds)
    masks = active_masks(inst, bounds)
    findings = NamedTuple[]
    # skip the subgrid hierarchy itself (active/wt arrays are Int/Bool anyway) and
    # the non-numeric attachments.
    skip = Set([:surfdata, :fates, :scf_method, :overrides, :decomp_cascade])
    for f in fieldnames(typeof(inst))
        f in skip && continue
        obj = getfield(inst, f)
        obj === nothing && continue
        append!(findings, sweep_struct(obj, String(f), bounds, masks))
        # one level of nesting (water.waterstatebulk_inst, bgc_vegetation.*, …)
        for g in fieldnames(typeof(obj))
            sub = getfield(obj, g)
            (sub isa AbstractArray || sub isa Number || sub isa AbstractString ||
             sub isa Bool || sub === nothing) && continue
            isstructtype(typeof(sub)) || continue
            append!(findings, sweep_struct(sub, "$(f).$(g)", bounds, masks))
            for h in fieldnames(typeof(sub))
                sub2 = getfield(sub, h)
                (sub2 isa AbstractArray || sub2 isa Number || sub2 isa AbstractString ||
                 sub2 isa Bool || sub2 === nothing) && continue
                isstructtype(typeof(sub2)) || continue
                append!(findings, sweep_struct(sub2, "$(f).$(g).$(h)", bounds, masks))
            end
        end
    end
    return findings
end

function report(tag, findings)
    real_ = filter(x -> !x.deliberate, findings)
    delib = filter(x ->  x.deliberate, findings)
    @printf("\n=== %s: %d NaN-on-active fields (%d deliberate-missing-value) ===\n",
            tag, length(real_), length(delib))
    for x in sort(real_, by = y -> (y.struct_name, String(y.field)))
        @printf("  %-46s %-38s %-6s %d/%d active (first=%d)\n",
                x.struct_name, String(x.field), String(x.level),
                x.nbad, x.nactive, x.firstidx)
    end
    isempty(delib) || println("  [deliberate NaN, preserved]: ",
                              join(sort(unique(String.(getfield.(delib, :field)))), ", "))
    return real_
end

# -----------------------------------------------------------------------------
# Domain drivers
# -----------------------------------------------------------------------------
function set_generic_forcing!(inst, ng, nc; T=283.0, snow=0.0, rain=1.0e-5, topo=1400.0)
    a = inst.atm2lnd
    pbot = 90000.0; q = 0.005
    th = T * (1.0e5 / pbot)^0.286
    rho = pbot / (287.058 * T * (1.0 + 0.61 * q))
    vp = q * pbot / (0.622 + 0.378 * q)
    for g in 1:ng
        a.forc_t_not_downscaled_grc[g]    = T
        a.forc_th_not_downscaled_grc[g]   = th
        a.forc_pbot_not_downscaled_grc[g] = pbot
        a.forc_q_not_downscaled_grc[g]    = q
        a.forc_rho_not_downscaled_grc[g]  = rho
        a.forc_lwrad_not_downscaled_grc[g]= 300.0
        a.forc_rain_not_downscaled_grc[g] = rain
        a.forc_snow_not_downscaled_grc[g] = snow
        a.forc_u_grc[g] = 3.0; a.forc_v_grc[g] = 0.0
        isempty(a.forc_wind_grc) || (a.forc_wind_grc[g] = 3.0)
        a.forc_hgt_grc[g] = 30.0
        a.forc_hgt_u_grc[g] = 30.0; a.forc_hgt_t_grc[g] = 30.0; a.forc_hgt_q_grc[g] = 30.0
        a.forc_vp_grc[g] = vp
        a.forc_pco2_grc[g] = 367.0e-6 * pbot
        a.forc_po2_grc[g]  = 0.209 * pbot
        a.forc_topo_grc[g] = topo
        for b in 1:size(a.forc_solad_not_downscaled_grc, 2)
            a.forc_solad_not_downscaled_grc[g, b] = b == 1 ? 300.0 : 250.0
            a.forc_solai_grc[g, b]                = b == 1 ? 100.0 : 80.0
        end
    end
    for c in 1:nc; inst.topo.topo_col[c] = topo; end
    return nothing
end

function run_domain(tag; fsurdat, use_cn=false, nsteps=1, T=283.0, snow=0.0,
                    use_bedrock=true, use_aquifer_layer=false, h2osfcflag=0,
                    start_date=DateTime(2006, 6, 15, 20), dtime=3600)
    isfile(fsurdat) || (@info "$tag: surfdata absent, skipping" fsurdat; return nothing)
    isfile(BOW_FP)  || (@info "$tag: param file absent, skipping"; return nothing)
    CLM.rooting_profile_config.rooting_profile_method_water  = CLM.JACKSON_1996_ROOT
    CLM.rooting_profile_config.rooting_profile_method_carbon = CLM.JACKSON_1996_ROOT
    (inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat=fsurdat, paramfile=BOW_FP,
        start_date=start_date, dtime=dtime, use_cn=use_cn, use_luna=false,
        use_bedrock=use_bedrock, use_aquifer_layer=use_aquifer_layer,
        h2osfcflag=h2osfcflag, fsnowoptics=SNOWOPT, fsnowaging=SNOWAGE,
        int_snow_max=2000.0)
    ng, nc, np = bounds.endg, bounds.endc, bounds.endp

    post_init = sweep_instances(inst, bounds)
    report("$tag  [after clm_initialize!]", post_init)

    config  = CLM.CLMDriverConfig(use_cn=use_cn, use_aquifer_layer=use_aquifer_layer,
                                  use_hydrstress=false, use_luna=false)
    filt_ia = CLM.clump_filter_inactive_and_active
    for i in 1:nsteps
        set_generic_forcing!(inst, ng, nc; T=T, snow=snow)
        calday = CLM.get_curr_calday(tm); (declin, _) = CLM.compute_orbital(calday)
        CLM.init_daylength!(inst.gridcell, declin, declin, CLM.ORB_OBLIQR_DEFAULT, 1:ng)
        CLM.advance_timestep!(tm)
        CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
        (yr, mon, d, tod) = CLM.get_curr_date(tm)
        nextsw = calday + dtime / CLM.SECSPDAY
        CLM.clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw, declin, declin,
            CLM.ORB_OBLIQR_DEFAULT, false, false, "", false;
            nstep=tm.nstep, is_first_step=(i == 1),
            is_beg_curr_day=CLM.is_beg_curr_day(tm), is_end_curr_day=CLM.is_end_curr_day(tm),
            is_beg_curr_year=CLM.is_beg_curr_year(tm), dtime=Float64(dtime), mon=mon, day=d,
            photosyns=inst.photosyns)
    end
    post_step = sweep_instances(inst, bounds)
    report("$tag  [after clm_initialize! + $nsteps clm_drv! step(s)]", post_step)
    return (inst=inst, bounds=bounds, post_init=post_init, post_step=post_step)
end

const DOMAINS = Dict(
    "bow"     => () -> run_domain("BOW (soil+lake+…)"; fsurdat=BOW_FS),
    "bow_cn"  => () -> run_domain("BOW + use_cn"; fsurdat=BOW_FS, use_cn=true),
    "glacier" => () -> run_domain("GLACIER (istice)"; fsurdat=GLAC_FS, T=263.0, snow=5.0e-5),
    "lake"    => () -> run_domain("LAKE (istdlak)"; fsurdat=LAKE_FS),
    "urban"   => () -> run_domain("URBAN (mexicocity)"; fsurdat=URB_FS),
)

function main(which = collect(keys(DOMAINS)))
    for k in sort(collect(which))
        haskey(DOMAINS, k) || continue
        try
            DOMAINS[k]()
        catch e
            @error "domain $k threw" exception=(e, catch_backtrace())
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(isempty(ARGS) ? collect(keys(DOMAINS)) : ARGS)
end
