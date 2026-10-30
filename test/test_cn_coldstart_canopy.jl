# ==========================================================================
# test_cn_coldstart_canopy.jl — a use_cn cold start begins with CTSM's zero canopy.
#
# CTSM fills tlai/tsai/elai/esai/htop/hbot from the surface-data monthly LAI/SAI
# (SatellitePhenologyInit / interpMonthlyVeg) only for SP, FATES-SP or dry
# deposition (clm_initializeMod.F90:461-486, 697-707). Under use_cn they keep
# CanopyStateType::InitCold's zeros (CanopyStateType.F90:529-536) and
# CNVegStructUpdate builds the canopy from leaf/stem C.
#
# The port seeded the SP values on every cold start. Under use_cn that handed
# CNVegStructUpdate the surface SAI (2.77 for Bow's needleleaf, vs CN's 0.5) as
# tsai_old, which its tsai = max(tsai_alpha*tsai_old + ..., tsai_min) relaxes only
# over months: the cold-start tree carried 5x Fortran's stem area through the
# first winter (root conductance x2.9 via the root-area index, xylem conductance x5.5).
#
#   julia +1.12 --project=. test/test_cn_coldstart_canopy.jl
# ==========================================================================
using Test, CLM, Dates

include(joinpath(@__DIR__, "testdata.jl"))

@testset "use_cn cold start: zero canopy state (CanopyStateType InitCold)" begin
    common = joinpath(@__DIR__, "..", "scripts", "fortran_parity_common.jl")
    Base.include(@__MODULE__, common)
    if !isfile(FSURDAT) || !isfile(FPARAM)
        testdata_missing("use_cn cold-start canopy", FSURDAT, FPARAM) && return
    end
    init(use_cn) = CLM.clm_initialize!(;
        fsurdat = FSURDAT, paramfile = FPARAM,
        start_date = DateTime(2002, 1, 1), dtime = 3600,
        use_cn = use_cn, use_luna = true, use_bedrock = true, use_aquifer_layer = false,
        h2osfcflag = 0, fsnowoptics = FSNOWOPT, fsnowaging = FSNOWAGE,
        int_snow_max = INT_SNOW_MAX)

    (inst, bounds, _, _) = init(true)
    cs = inst.canopystate
    veg = [p for p in bounds.begp:bounds.endp if inst.patch.active[p] && inst.patch.itype[p] > 0]
    @test !isempty(veg)
    for f in (:tlai_patch, :tsai_patch, :elai_patch, :esai_patch, :htop_patch, :hbot_patch)
        @test all(getfield(cs, f)[p] == 0.0 for p in bounds.begp:bounds.endp)
    end
    @test all(cs.frac_veg_nosno_alb_patch[p] == 0 for p in veg)
    # The SP interpolation arrays are still read (harmless, and what SP needs).
    @test any(inst.satellite_phenology.msai2t[p, 1] > 0 for p in veg)

    # Non-vacuity: the SP cold start still seeds the canopy from the surface data.
    (inst_sp, bounds_sp, _, _) = init(false)
    @test any(inst_sp.canopystate.tsai_patch[p] > 0 for p in veg)
    @test any(inst_sp.canopystate.tlai_patch[p] > 0 for p in veg)
end
