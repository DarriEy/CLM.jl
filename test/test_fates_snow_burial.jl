# ==========================================================================
# test_fates_snow_burial.jl — the host snow pack must reach FATES.
#
# Fortran wrap_update_hlmfates_dyn (clmfates_interfaceMod.F90:1494-1504), which runs
# after the daily dynamics (and at the cold start), opens with
#     bc_in(s)%snow_depth_si   = snow_depth(c)
#     bc_in(s)%frac_sno_eff_si = frac_sno_eff(c)
#     call UpdateFatesAvgSnowDepth(sites, bc_in)    ! site%snow_depth = depth * fsno
# and the leaf-area profile buries the canopy layers below site%snow_depth, so the
# HLM's elai/esai fall below tlai/tsai in winter. The port never packed the two
# fields, so FATES saw no snow: at Bow the exposed LAI equalled the total LAI all
# winter, the canopy absorbed far more SW, and the Feb-May leaf ran 0.3-1.1 K warm
# (FATES_TVEG24), which also sped up the spring grass flush.
#
# A. fates_case_snow_namelist! sets the Fortran FATES case snow namelist
#    (Vionnet2012 overburden compaction, wind-dependent new-snow density).
# B. Live Bow FATES site (data-gated): the daily step packs the host snow depth,
#    site%snow_depth = snow_depth * frac_sno_eff, and a 1.2 m snowpack buries part
#    of the ~1.9 m canopy (elai < tlai), while no snow leaves it exposed.
# ==========================================================================
using Test
using CLM
include(joinpath(@__DIR__, "testdata.jl"))

@testset "FATES snow burial + case snow namelist" begin
    ENV["SYMFLUENCE_DATA"] = symfluence_data_root()
    m = Module(:_FatesSnowBurialHarness)
    Core.eval(m, :(include(p::AbstractString) = Base.include($m, p)))
    Base.include(m, joinpath(@__DIR__, "..", "scripts", "fates_fortran_parity.jl"))
    saved = (CLM.WIND_DEPENDENT_SNOW_DENSITY[], CLM.OVERBURDEN_COMPACTION_METHOD[],
             CLM.NEW_SNOW_DENSITY[])
    try
        @testset "A. Fortran case snow namelist" begin
            CLM.snow_hydrology_set_control_for_testing!(; wind_dep_snow_density = false,
                overburden_compaction_method = CLM.OVERBURDEN_COMPACTION_ANDERSON1976)
            old = Base.invokelatest(getglobal(m, :fates_case_snow_namelist!))
            @test old.wind_dep_snow_density == false
            @test CLM.WIND_DEPENDENT_SNOW_DENSITY[] == true
            @test CLM.OVERBURDEN_COMPACTION_METHOD[] == CLM.OVERBURDEN_COMPACTION_VIONNET2012
            @test CLM.NEW_SNOW_DENSITY[] == CLM.LO_TMP_DNS_SLATER2017
        end

        @testset "B. live FATES daily step buries the canopy under the host snow" begin
            sitedir = symfluence_path("domain_Bow_at_Banff_lumped_era5", "settings", "CLM", "parameters")
            fsurdat = joinpath(sitedir, "surfdata_clm.nc")
            fparam  = joinpath(sitedir, "clm5_params.nc")
            if !(isfile(fsurdat) && isfile(fparam))
                testdata_missing("FATES snow burial (live site)", fsurdat, fparam)
                return
            end
            inst, bounds, filt, filt_ia, config, _tm = Base.invokelatest(getglobal(m, :build))
            c = findfirst(inst.column.is_fates)
            wdb = inst.water.waterdiagnosticbulk_inst
            cs = inst.canopystate
            site = inst.fates.sites[1]
            nlevsoil = inst.fates.bc_in[1].nlevsoil
            nlevdecomp = length(inst.fates.bc_in[1].w_scalar_sisl)
            veg = [p for (ifp, p) in CLM.fates_veg_patches(site, c, inst.column)]
            @test !isempty(veg)
            daily!() = CLM.fates_daily_dynamics_step!(inst; nlevsoil = nlevsoil,
                nlevdecomp = nlevdecomp, use_fates_bgc = config.use_fates_bgc)

            # no snow: the whole canopy is exposed
            wdb.snow_depth_col[c] = 0.0; wdb.frac_sno_eff_col[c] = 0.0
            daily!()
            @test site.snow_depth == 0.0
            @test all(isapprox(cs.elai_patch[p], cs.tlai_patch[p]; rtol = 1e-12) for p in veg)

            # 1.2 m of snow over 90 % of the column: site%snow_depth = 1.08 m
            wdb.snow_depth_col[c] = 1.2; wdb.frac_sno_eff_col[c] = 0.9
            daily!()
            @test inst.fates.bc_in[1].snow_depth_si == 1.2
            @test inst.fates.bc_in[1].frac_sno_eff_si == 0.9
            @test site.snow_depth ≈ 1.08
            p = argmax(p -> cs.htop_patch[p], veg)
            @test cs.htop_patch[p] > site.snow_depth           # canopy top above the snow
            @test cs.elai_patch[p] < cs.tlai_patch[p] - 1e-3    # but its lower part buried
            @test all(cs.elai_patch[q] <= cs.tlai_patch[q] + 1e-12 for q in veg)
        end
    finally
        CLM.WIND_DEPENDENT_SNOW_DENSITY[], CLM.OVERBURDEN_COMPACTION_METHOD[],
            CLM.NEW_SNOW_DENSITY[] = saved
    end
end
