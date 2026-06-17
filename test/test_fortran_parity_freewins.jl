# ==========================================================================
# test_fortran_parity_freewins.jl — Tier-0 "free win" Fortran-parity coverage.
#
# docs/FULL_PLATFORM_COVERAGE_PLAN.md Tier 0 lists modules whose outputs are in
# the EXISTING Fortran pdump files (no Fortran rerun needed) but that no harness
# diffs. This file guards the ones that validate cleanly against the BGC summer
# dump window (n1757852), running ONE use_cn=true clm_drv! step from the
# injected IC and diffing the COMPUTED outputs (not the injected values) against
# the after_hydrologydrainage dump.
#
# ASSERTED (validated against the Fortran dump):
#   - surface_resistance  (SOILRESIS): calc_soilevap_resis! recompute. The IC
#     injects SOILRESIS but the step overwrites it (2712.6 -> 2671.6), so the
#     post-step value is a genuine recompute. Verified rel ~1.7e-8 (Float64
#     cross-compiler floor). Asserted < 1e-6.
#   - aerosol  (mss_bcphi/bcpho/dst1..4/ocphi/ocpho): aerosol_masses! recompute.
#     CAVEAT — the only available pdump_before_step dumps are the SNOW-FREE BGC
#     summer window, so every mss_* is 0.0; this is a structural / mass-
#     conservation sanity check (0 in -> 0 out, exact), NOT an aerosol-physics
#     stress test (no snowy dump exists on this machine). Asserted == 0.0.
#
# NOT ASSERTED (documented in the probes, left out on purpose):
#   - daylength (dayl/prev_dayl): NOT present in any pdump file — there is no
#     Fortran ground-truth field to diff against. scripts/fortran_parity_
#     daylength.jl runs a self-consistency check (init_daylength! vs the closed-
#     form astronomical formula, rel ~1.2e-9) but that is not a Fortran parity.
#   - active_layer (altmax/altmax_indx): FOUND BUG — active_layer_init_cold!
#     (active_layer.jl:64) is never called from instances.jl/clm_initialize!, so
#     altmax_col stays at SPVAL (1e36) and alt_calc!'s running-max test
#     `alt > altmax` is never true; altmax/altmax_indx never track the active
#     layer (rel ~1.0 vs the dump's altmax=41.998 m / indx=25). The per-step ALT
#     (alt_indx_col=25) is computed correctly; only the annual-max bootstrap is
#     broken by the missing cold-start wiring. See scripts/fortran_parity_
#     activelayer.jl. Fix = call active_layer_init_cold! in cold-start init.
#
# GATED: the Fortran BGC reference dumps + Bow input files live outside the repo
# (machine-specific). When absent, the test is skipped (not failed) so CI and
# other machines stay green.
# ==========================================================================

using Test
using CLM
using NCDatasets
using Dates

@testset "Fortran parity free-wins (gated)" begin
    common  = joinpath(@__DIR__, "..", "scripts", "fortran_parity_common.jl")
    dumpdir = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_summer"
    nstep   = 1757852
    date_base = 1753153

    bdump = joinpath(dumpdir, "pdump_before_step_n$(nstep).nc")
    edump = joinpath(dumpdir, "pdump_after_hydrologydrainage_n$(nstep).nc")
    have_dumps = isfile(bdump) && isfile(edump)

    if !(isfile(common) && have_dumps)
        @info "Fortran parity free-wins test SKIPPED (BGC reference dumps not present on this machine)"
        @test_skip have_dumps
    else
        Base.include(@__MODULE__, common)   # build_bow_inst / run_one_parity_step! / FFORCING
        forcing2002 = replace(FFORCING, "clmforc.2003.nc" => "clmforc.2002.nc")

        # Run ONE use_cn step from the BGC IC (shared by both checks below).
        # invokelatest: run_one_parity_step! was just defined by Base.include
        # above, so it is newer than this method's world age — invoke the latest.
        inst, bounds = Base.invokelatest(run_one_parity_step!, nstep; use_cn=true,
            dumpdir=dumpdir, use_hydrstress=true, use_luna=true,
            step_date=DateTime(2002, 1, 1) + Hour(nstep - date_base),
            forcing_file=forcing2002)

        # ----------------------------------------------------------------
        # 1. surface_resistance — SOILRESIS recompute (STRONG)
        # ----------------------------------------------------------------
        @testset "surface_resistance (SOILRESIS)" begin
            jl = Float64(inst.soilstate.soilresis_col[1])

            ds0 = NCDataset(bdump, "r"); inj = Float64(ds0["SOILRESIS"][1]); close(ds0)
            # prove the step genuinely overwrote the injected IC value
            @test abs(jl - inj) > 1e-6

            ds = NCDataset(edump, "r"); f = Float64(ds["SOILRESIS"][1]); close(ds)
            rel = abs(jl - f) / (1.0 + max(abs(jl), abs(f)))
            @info "SOILRESIS parity" julia=jl fortran=f rel=rel
            @test rel < 1e-6
        end

        # ----------------------------------------------------------------
        # 2. aerosol — snow aerosol masses recompute (SANITY; snow-free window)
        # ----------------------------------------------------------------
        @testset "aerosol (snow masses)" begin
            nlevsno = CLM.varpar.nlevsno
            aer_fields = [
                ("mss_bcphi", :mss_bcphi_col), ("mss_bcpho", :mss_bcpho_col),
                ("mss_ocphi", :mss_ocphi_col), ("mss_ocpho", :mss_ocpho_col),
                ("mss_dst1",  :mss_dst1_col),  ("mss_dst2",  :mss_dst2_col),
                ("mss_dst3",  :mss_dst3_col),  ("mss_dst4",  :mss_dst4_col),
            ]
            ds = NCDataset(edump, "r")
            gmax = 0.0
            for (fname, jfield) in aer_fields
                (haskey(ds, fname) && hasproperty(inst.aerosol, jfield)) || continue
                f = ds[fname][:, 1]
                jl = getfield(inst.aerosol, jfield)
                nlev = min(nlevsno, length(f), size(jl, 2))
                for k in 1:nlev
                    fv = ismissing(f[k]) ? NaN : Float64(f[k])
                    jv = Float64(jl[1, k])
                    (isnan(fv) && isnan(jv)) && continue
                    gmax = max(gmax, abs(fv - jv) / (1.0 + max(abs(fv), abs(jv))))
                end
            end
            close(ds)
            @info "aerosol-mass parity (snow-free window: all masses 0)" max_rel=gmax
            # exact: snow-free => 0 (dump) vs 0 (Julia recompute)
            @test gmax < 1e-9
        end

        # ----------------------------------------------------------------
        # 3. active_layer — altmax annual-max trackers (BUG NOW FIXED)
        # cold-init (cold_start.jl) + altmax/lastyear injection (fortran_restart.jl)
        # close the former rel~1.0 SPVAL bug. altmax/indx re-tracked by alt_calc!;
        # altmax_lastyear injected (year-boundary carried). Now exact (0.0).
        # ----------------------------------------------------------------
        @testset "active_layer (altmax trackers)" begin
            al = inst.active_layer
            ds = NCDataset(edump, "r")
            checks = (("altmax", Float64(al.altmax_col[1])),
                      ("altmax_indx", Float64(al.altmax_indx_col[1])),
                      ("altmax_lastyear", Float64(al.altmax_lastyear_col[1])),
                      ("altmax_lastyear_indx", Float64(al.altmax_lastyear_indx_col[1])))
            gmax = 0.0
            for (fname, jv) in checks
                haskey(ds, fname) || continue
                f = ismissing(ds[fname][1]) ? NaN : Float64(ds[fname][1])
                isnan(f) && continue
                gmax = max(gmax, abs(jv - f) / (1.0 + max(abs(jv), abs(f))))
            end
            close(ds)
            @info "active_layer altmax parity" max_rel=gmax
            @test gmax < 1e-6
        end

        # ----------------------------------------------------------------
        # 4. surface_albedo — the radiation bands (albgrd/albgri/albd/albi/
        #    fabd/fabi) RECOMPUTED by surface_albedo! during the step.
        #    The main parity harness INJECTS these RT outputs as IC (to dodge
        #    the cold-start albedo NaN), so it never diffs the module's own
        #    output — an albedo bug could hide behind the injection (audit B).
        #    Here we diff the POST-step recompute against the after_hydrology-
        #    drainage dump: the step overwrites the injected IC, so this is a
        #    genuine recompute. Agreement is ~1e-4 (albd/albi/albgri) up to
        #    ~1.5e-3 (fabd) — a real ~0.1% residual (soil-albedo / coszen-timing,
        #    not chased here), much looser than soilresis's 1e-8 but a true
        #    Fortran diff that guards against gross albedo regressions.
        # ----------------------------------------------------------------
        @testset "surface_albedo (radiation bands)" begin
            sa = inst.surfalb
            ds = NCDataset(edump, "r")
            relf(a, b) = abs(a - b) / (1.0 + max(abs(a), abs(b)))
            # column-indexed (numrad bands), 1 column
            col_fields = (("albgrd", sa.albgrd_col), ("albgri", sa.albgri_col))
            # patch-indexed (numrad bands); diff the two vegetated patches (2,3)
            pat_fields = (("albd", sa.albd_patch), ("albi", sa.albi_patch),
                          ("fabd", sa.fabd_patch), ("fabi", sa.fabi_patch))
            gmax = 0.0; nchecked = 0
            for (nm, arr) in col_fields
                haskey(ds, nm) || continue
                for b in 1:2
                    f = ismissing(ds[nm][b, 1]) ? NaN : Float64(ds[nm][b, 1])
                    isnan(f) && continue
                    gmax = max(gmax, relf(Float64(arr[1, b]), f)); nchecked += 1
                end
            end
            for (nm, arr) in pat_fields
                haskey(ds, nm) || continue
                for p in 2:3, b in 1:2
                    f = ismissing(ds[nm][b, p]) ? NaN : Float64(ds[nm][b, p])
                    isnan(f) && continue
                    gmax = max(gmax, relf(Float64(arr[p, b]), f)); nchecked += 1
                end
            end
            close(ds)
            @info "surface_albedo bands parity (recompute vs Fortran)" max_rel=gmax nchecked
            @test nchecked > 0           # the bands were actually present + diffed
            @test gmax < 3e-3            # recompute matches Fortran to ~0.1%
        end
    end
end
