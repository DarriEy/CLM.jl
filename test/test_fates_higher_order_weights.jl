# test_fates_higher_order_weights.jl
# After FATES re-weights its patches on the column (dyn_ED: patch%wtcol = wt_ed),
# CTSM recomputes the higher-order weights before anything reads them:
#
#   dynSubgrid_driver (dynSubgridDriverMod.F90:313)  call dyn_ED(bounds_clump)
#     -> dynSubgrid_wrapup_weight_changes (:333)
#        -> compute_higher_order_weights (:389)   patch%wtlunit/wtgcell, col%wtgcell
#        -> reweight_wrapup (:395)                set_active + setFilters
#
# The port ran dyn_ED (fates_apply_deferred_activation!) but never the recompute,
# so FATES patches kept their SURFDATA wtgcell, and the history writer averaged
# patch fields to the gridcell with those stale weights and without the
# patch%active test of p2g_1d (subgridAveMod.F90:531-538).
#
#   PART A (synthetic): the history writer's gridcell aggregation follows p2g_1d /
#     c2g_1d: inactive, zero-weight and spval units are skipped; l2g_scale_type
#     'veg' drops non-vegetated landunits.
#   PART B (data-gated, Bow FATES-SP site through the real clm_drv!): at every step
#     patch%wtgcell == wtcol * col%wtgcell and wtlunit == wtcol * col%wtlunit; an
#     inactive slot carries zero gridcell weight; and the history-writer gridcell
#     BTRAN / TV / FSH match the Fortran h0 (hist_dov2xy default) of the same run.
using Test
using CLM

include("testdata.jl")

@testset "FATES higher-order weights + history p2g" begin
    @testset "A: history-writer p2g/c2g skips inactive / zero-weight / spval" begin
        inst = CLM.CLMInstances()
        lun = inst.landunit; col = inst.column; pch = inst.patch
        lun.itype = [CLM.ISTSOIL, CLM.ISTDLAK]; lun.gridcell = [1, 1]; lun.wtgcell = [0.8, 0.2]
        lun.active = [true, true]
        col.landunit = [1, 2]; col.gridcell = [1, 1]; col.wtgcell = [0.8, 0.2]
        col.active = [true, true]
        pch.column = [1, 1, 1, 2]; pch.landunit = [1, 1, 1, 2]; pch.gridcell = [1, 1, 1, 1]
        # p3: inactive with a stale nonzero weight (the pre-fix FATES state)
        pch.wtgcell = [0.5, 0.3, 0.1, 0.2]
        pch.active = [true, true, false, true]
        vals = [1.0, 2.0, 100.0, 4.0]
        f = CLM.HistFieldDef("X", "x", "1", "patch", _ -> vals)
        @test CLM._aggregate_to_gridcell(f, inst, 1)[1] ≈ (0.5*1 + 0.3*2 + 0.2*4) / 1.0
        # spval is "no value", as in p2g_1d
        vals[2] = CLM.SPVAL
        @test CLM._aggregate_to_gridcell(f, inst, 1)[1] ≈ (0.5*1 + 0.2*4) / 0.7
        vals[2] = 2.0
        # l2g_scale_type 'veg' drops the lake landunit
        fveg = CLM.HistFieldDef("X", "x", "1", "patch", _ -> vals, :avg, "veg")
        @test CLM._aggregate_to_gridcell(fveg, inst, 1)[1] ≈ (0.5*1 + 0.3*2) / 0.8
        # c2g: inactive column skipped
        col.active[2] = false
        fc = CLM.HistFieldDef("Y", "y", "1", "column", _ -> [3.0, 50.0])
        @test CLM._aggregate_to_gridcell(fc, inst, 1)[1] ≈ 3.0
    end

    @testset "B: clm_drv! recomputes wtgcell after dyn_ED; history vs Fortran h0" begin
        sitedir = symfluence_path("domain_Bow_at_Banff_lumped_era5", "settings", "CLM", "parameters")
        fsurdat = joinpath(sitedir, "surfdata_clm.nc")
        forcing = symfluence_path("domain_Bow_at_Banff_lumped_era5", "data", "forcing",
                                  "CLM_input", "clmforc.2002.nc")
        if !(isfile(fsurdat) && isfile(forcing))
            testdata_missing("FATES higher-order weights run", fsurdat, forcing)
        else
            harness = joinpath(@__DIR__, "..", "scripts", "validation", "fates_fortran_parity",
                               "fates_modes_parity.jl")
            ENV["SYMFLUENCE_DATA"] = symfluence_data_root()
            m = Module(:_FatesHiWeightsHarness)
            Core.eval(m, :(include(p::AbstractString) = Base.include($m, p)))
            Base.include(m, harness)
            G(sym) = getglobal(m, sym)
            nsteps = 36   # 1.5 days, crossing the nstep-25 daily FATES step
            bad_hi = Ref(0); bad_inactive = Ref(0); bad_sum = Ref(0); nveg = Ref(0)
            hook = (inst, i) -> begin
                i >= 1 || return
                col = inst.column; pch = inst.patch
                c = findfirst(col.is_fates)
                for p in eachindex(pch.column)
                    cc = pch.column[p]
                    pch.wtgcell[p] == pch.wtcol[p] * col.wtgcell[cc] || (bad_hi[] += 1)
                    pch.wtlunit[p] == pch.wtcol[p] * col.wtlunit[cc] || (bad_hi[] += 1)
                end
                s = 0.0
                for p in col.patchi[c]:col.patchf[c]
                    if pch.active[p]
                        s += pch.wtgcell[p]
                        pch.is_veg[p] && (nveg[] += 1)
                    else
                        pch.wtgcell[p] == 0.0 || (bad_inactive[] += 1)
                    end
                end
                isapprox(s, col.wtgcell[c]; rtol = 1e-12) || (bad_sum[] += 1)
            end
            jl = Base.invokelatest(G(:run_julia), "sp", nsteps; hook = hook)
            @test bad_hi[] == 0
            @test bad_inactive[] == 0
            @test bad_sum[] == 0
            @test nveg[] > 0

            ft = Base.invokelatest(G(:read_fortran), G(:fortran_csv)("sp"))
            # Per-step relative error over nstep 2..nsteps of the history-writer
            # gridcell value against the Fortran h0 (floor 1e-3 x peak |F|). nstep 1
            # is left out: there the per-patch btran itself differs (0.135 vs 0.036;
            # TV/FSH agree), which is not a weighting question.
            function relerr(key)
                F = ft[key][3:nsteps+1]; J = [jl[k][key] for k in 3:nsteps+1]
                peak = maximum(abs, F)
                return maximum(abs(J[k] - F[k]) / max(abs(F[k]), 1e-3peak) for k in eachindex(F))
            end
            errs = Dict(k => relerr("HIST_" * k) for k in ("BTRAN", "TV", "FSH", "FSA", "FIRE"))
            @info "FATES-SP history vs Fortran h0, max per-step rel err (nstep 2-$nsteps)" errs
            # Pre-fix (stale surfdata wtgcell, inactive slots averaged in) the 240-step
            # window-mean errors were BTRAN 0.82, TV 4.5e-3 (1.3 K), FSH 0.53, FSA 0.52,
            # FIRE 0.53. Measured after: BTRAN ~1e-5, TV ~1e-6, FSH 1.7e-3, FSA/FIRE ~2e-5.
            @test errs["BTRAN"] < 1e-4
            @test errs["TV"] < 5e-5
            @test errs["FSH"] < 5e-3
            @test errs["FSA"] < 1e-4
            @test errs["FIRE"] < 1e-4
        end
    end
end
