# ==========================================================================
# test_coldstart_radiation.jl — the cold-start radiation that drives the first-steps
# soil thaw (the "cold-start thaw" host residual of #336 / #345).
#
#   * init_surface_albedo_cold! keeps Fortran surfalb_type%InitCold's 0.2 in every band
#     on non-glacier columns (SurfaceAlbedoType.F90:370-371). It used to write vis 0.20 /
#     NIR 0.40 (broadband 0.3), which cut SABG by 48 W/m2 on the first steps at a summer
#     Bow cold start and delayed the top-soil thaw (SOILLIQ_1..4 9-16 % low at nstep 2).
#   * nuopc_doalb mirrors the NUOPC cap: no albedo solve on nsteps 0 and 1
#     (lnd_comp_nuopc.F90:830-841), so nsteps 0-2 radiate with the InitCold albedos.
# ==========================================================================
using Test, CLM

@testset "cold-start radiation (InitCold albedo + NUOPC doalb)" begin
    @testset "nuopc_doalb: no albedo solve on nsteps 0 and 1" begin
        @test CLM.nuopc_doalb(0) == false
        @test CLM.nuopc_doalb(1) == false
        @test CLM.nuopc_doalb(2) == true
        @test all(CLM.nuopc_doalb, 3:100)
    end

    @testset "init_surface_albedo_cold!: Fortran 0.2 on soil, ALBICE on glacier" begin
        nc, np, ng = 3, 3, 1
        inst = CLM.CLMInstances()
        inst.landunit = CLM.LandunitData{Float64}()
        inst.landunit.itype = [CLM.ISTSOIL, CLM.ISTICE, CLM.ISTDLAK]
        inst.column = CLM.ColumnData{Float64}()
        inst.column.landunit = [1, 2, 3]
        CLM.surfalb_init!(inst.surfalb, np, nc, ng)
        CLM.surfalb_init_cold!(inst.surfalb, 1:nc, 1:np)   # Fortran InitCold (0.2)
        bounds = CLM.BoundsType(begg=1, endg=1, begl=1, endl=3, begc=1, endc=nc,
                                begp=1, endp=np, level=CLM.BOUNDS_LEVEL_PROC)
        CLM.init_surface_albedo_cold!(inst, bounds)
        sa = inst.surfalb
        for c in (1, 3), ib in 1:CLM.NUMRAD     # soil + lake: untouched InitCold value
            @test sa.albgrd_col[c, ib] == 0.2
            @test sa.albgri_col[c, ib] == 0.2
        end
        @test sa.albgrd_col[2, :] == [0.80, 0.55]   # glacier-only deviation
        @test sa.albgri_col[2, :] == [0.80, 0.55]
        # InitCold canopy transmission: everything reaches the ground on the first steps.
        @test all(==(1.0), sa.ftdd_patch) && all(==(0.0), sa.fabd_patch)
    end
end
