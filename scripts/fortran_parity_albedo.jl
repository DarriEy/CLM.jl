# ==========================================================================
# fortran_parity_albedo.jl — SurfaceAlbedo band parity, WITHOUT injection.
#
# The audit (docs/HARNESS_COVERAGE_AUDIT.md, Table B) flags that the dumps
# carry the surface-albedo bands albd/albi/fabd/fabi (and fsun, albgrd/albgri,
# albsod/albsoi) but NO harness diffs them — only the derived SABV/SABG are
# checked, AND the parity harness INJECTS the RT outputs (read_fortran_restart!
# sets albd/albi/fabd/fabi/fsun/fsun_z directly), which can MASK an albedo-RT
# residual.
#
# This probe lets Julia COMPUTE the albedo bands and diffs them against the
# dump WITHOUT injecting them. Procedure:
#   1. build the Bow inst, inject the before_step dump (sets canopy/soil/snow
#      state used by SurfaceAlbedo — h2osoi, snw_rds, t_grnd, elai/esai, etc.),
#   2. call surface_albedo! standalone with THIS step's declination + calday
#      (so the coszen matches the dump's `coszen` field, not next-step),
#   3. surface_albedo! OVERWRITES albd/albi/fabd/fabi/fsun → the injected
#      values are discarded, exposing any RT residual,
#   4. diff Julia's recomputed bands against the dump.
#
# Dump layout: albd/albi/fabd/fabi are (numrad=2, pft=3); Julia is (patch,
# numrad). fsun is (pft). albgrd/albgri/albsod/albsoi are (numrad, col).
#
#   julia +1.12 --project=. scripts/fortran_parity_albedo.jl
# ==========================================================================
include(joinpath(@__DIR__, "fortran_parity_common.jl"))

const BGC_DUMPDIR = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_summer"
const NSTEP = 1757852

function main()
    bdump = joinpath(BGC_DUMPDIR, "pdump_before_step_n$(NSTEP).nc")
    isfile(bdump) || (println("missing dump: $bdump"); return 1)

    # BGC run year-aligns 2202→2002; same step-date convention as cn_summer.
    step_date = DateTime(2002, 1, 1) + Hour(1757852 - 1753153)
    forcing   = replace(FFORCING, "clmforc.2003.nc" => "clmforc.2002.nc")

    (inst, bounds, filt, tm) = build_bow_inst(; dtime=3600, start_date=step_date, use_cn=true)

    # Allocate LUNA fields so inject_dump! is well-posed (matches run_one_parity_step!).
    if isempty(inst.photosyns.vcmx25_z_patch)
        inst.photosyns.vcmx25_z_patch = fill(30.0, bounds.endp, CLM.NLEVCAN)
        inst.photosyns.jmx25_z_patch  = fill(60.0, bounds.endp, CLM.NLEVCAN)
    end
    inject_dump!(inst, bounds, bdump)

    # --- record the INJECTED albedo bands so we can prove they were overwritten
    inj_albd = copy(inst.surfalb.albd_patch)

    # --- this step's orbital geometry (NOT nextsw_cday) so coszen matches dump
    calday = CLM.get_curr_calday(tm)
    (declin, _) = CLM.compute_orbital(calday)
    nextsw_cday = calday            # use THIS calday → albedo for the dumped step

    # local bindings mirroring clm_driver.jl
    alb  = inst.surfalb; alb_con = inst.surfalb_con
    grc  = inst.gridcell; col = inst.column; lun = inst.landunit; pch = inst.patch
    cs   = inst.canopystate; temp = inst.temperature
    wsb  = inst.water.waterstatebulk_inst; wdb = inst.water.waterdiagnosticbulk_inst
    ls   = inst.lakestate; aer = inst.aerosol
    bc_grc = 1:bounds.endg; bc_col = 1:bounds.endc; bc_patch = 1:bounds.endp
    pftc = CLM.pftcon

    coszen_cday = (cday, lat, lon, decl) -> begin
        hour_angle = 2.0 * π * mod(cday, 1.0) + lon - π
        max(sin(lat) * sin(decl) + cos(lat) * cos(decl) * cos(hour_angle), 0.0)
    end

    try
        CLM.surface_albedo!(alb, alb_con, grc, col, lun, pch, cs, temp, wsb, wdb, ls, aer,
                            filt.nourbanc, filt.nourbanp,
                            nextsw_cday, declin,
                            bc_grc, bc_col, bc_patch,
                            pftc.rhol, pftc.rhos, pftc.taul, pftc.taus, pftc.xl,
                            coszen_cday)
    catch e
        println("\n!!! surface_albedo! CRASHED:")
        Base.showerror(stdout, e, catch_backtrace()); println()
        return 2
    end

    overwritten = !isequal(inj_albd, inst.surfalb.albd_patch)
    @printf("Julia recomputed albd differs from injected? %s  (must be true → no masking)\n\n",
            overwritten ? "YES" : "NO — WARNING, injection NOT overwritten")

    ds = NCDataset(bdump, "r")
    fpfts = Int.(ds["pfts1d_itypveg"][:]); jpfts = Int.(inst.patch.itype)
    jl2f = Dict{Int,Int}()
    for pj in eachindex(jpfts), pf in eachindex(fpfts)
        if fpfts[pf] == jpfts[pj]; jl2f[pj] = pf; break; end
    end
    cz = haskey(ds, "coszen") ? Float64(ds["coszen"][1]) : NaN
    @printf("dump coszen[1] = %.4f   (Julia uses this step's declination → coszen match)\n\n", cz)

    gmax = 0.0; rows = Tuple{String,Float64,Float64,Int}[]
    # patch×numrad bands: dump (numrad, pft) → Julia (patch, numrad)
    function band2(name, jlmat)
        haskey(ds, name) || return
        f = ds[name][:, :]                    # (numrad, pft)
        nrad = size(f, 1)
        mabs = 0.0; mrel = 0.0; npts = 0
        for pj in eachindex(jpfts)
            haskey(jl2f, pj) || continue
            pf = jl2f[pj]
            for b in 1:nrad
                fv = ismissing(f[b, pf]) ? NaN : Float64(f[b, pf])
                jv = Float64(jlmat[pj, b])
                (isnan(fv) && isnan(jv)) && continue
                a = abs(fv - jv); mabs = max(mabs, a)
                mrel = max(mrel, a / (1 + max(abs(fv), abs(jv)))); npts += 1
            end
        end
        gmax = max(gmax, mrel); push!(rows, (name, mabs, mrel, npts))
    end
    # patch scalar
    function pscal(name, jlvec)
        haskey(ds, name) || return
        f = ds[name][:]
        mabs = 0.0; mrel = 0.0; npts = 0
        for pj in eachindex(jpfts)
            haskey(jl2f, pj) || continue
            fv = ismissing(f[jl2f[pj]]) ? NaN : Float64(f[jl2f[pj]])
            jv = Float64(jlvec[pj])
            (isnan(fv) && isnan(jv)) && continue
            a = abs(fv - jv); mabs = max(mabs, a)
            mrel = max(mrel, a / (1 + max(abs(fv), abs(jv)))); npts += 1
        end
        gmax = max(gmax, mrel); push!(rows, (name, mabs, mrel, npts))
    end
    # col×numrad bands: dump (numrad, col) → Julia (col, numrad)
    function colband(name, jlmat)
        haskey(ds, name) || return
        f = ds[name][:, :]; nrad = size(f, 1)
        mabs = 0.0; mrel = 0.0; npts = 0
        for b in 1:nrad
            fv = ismissing(f[b, 1]) ? NaN : Float64(f[b, 1])
            jv = Float64(jlmat[1, b])
            (isnan(fv) && isnan(jv)) && continue
            a = abs(fv - jv); mabs = max(mabs, a)
            mrel = max(mrel, a / (1 + max(abs(fv), abs(jv)))); npts += 1
        end
        gmax = max(gmax, mrel); push!(rows, (name, mabs, mrel, npts))
    end

    band2("albd", inst.surfalb.albd_patch)
    band2("albi", inst.surfalb.albi_patch)
    band2("fabd", inst.surfalb.fabd_patch)
    band2("fabi", inst.surfalb.fabi_patch)
    pscal("fsun", inst.canopystate.fsun_patch)
    colband("albgrd", inst.surfalb.albgrd_col)
    colband("albgri", inst.surfalb.albgri_col)
    colband("albsod", inst.surfalb.albsod_col)
    colband("albsoi", inst.surfalb.albsoi_col)
    close(ds)

    @printf("  %-10s %12s %12s %6s  %s\n", "band", "max|abs|", "max|rel|", "npts", "status")
    @printf("  %s\n", "-"^54)
    for (n, a, r, np) in rows
        @printf("  %-10s %12.3e %12.3e %6d  %s\n", n, a, r, np, r > 1e-4 ? "DIFF" : "ok")
    end
    @printf("  %s\n  surface-albedo band global max|rel| = %.3e  (computed, NOT injected)\n",
            "-"^54, gmax)
    return 0
end

exit(main())
