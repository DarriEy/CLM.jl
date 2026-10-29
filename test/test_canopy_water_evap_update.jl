# Canopy water phase bookkeeping after leaf-surface evaporation / dew.
#
# Regression test for the FATES-SP canopy-snow residual (SNOCAN about -9 % on the Bow
# tree patches vs ctsm5.3.012): the port used to re-split the canopy store in proportion
# to the existing phase amounts, while Fortran (CanopyFluxesMod.F90:1602-1629) works on
# the leaf-temperature phase first and only takes the remainder from the other phase.
using Test
using CLM

# Literal transcription of CanopyFluxesMod.F90:1602-1629 (+ truncate_small_values,
# NumericsMod.F90:92-97, custom_rel_epsilon = 1e-10).
function fortran_dew_update(liqcan, snocan, tran, evap, t_veg, dtime)
    tfrz = CLM.TFRZ
    snocan_baseline = snocan
    if t_veg > tfrz
        if (evap - tran) * dtime > liqcan
            snocan = max(0.0, snocan + liqcan + (tran - evap) * dtime)
        end
        liqcan = max(0.0, liqcan + (tran - evap) * dtime)
    elseif t_veg <= tfrz
        if (evap - tran) * dtime > snocan
            liqcan = liqcan + snocan + (tran - evap) * dtime
        end
        snocan = max(0.0, snocan + (tran - evap) * dtime)
    end
    if abs(snocan) < 1.0e-10 * abs(snocan_baseline)
        snocan = 0.0
    end
    return liqcan, snocan
end

@testset "canopy_water_evap_update matches CanopyFluxesMod.F90:1602-1629" begin
    dt = 3600.0
    tf = CLM.TFRZ
    # (liqcan, snocan, qflx_tran_veg, qflx_evap_veg, t_veg)
    cases = [
        (0.40, 0.30, 0.0, 5.0e-5, tf + 5.0),    # warm, evap < liqcan: snow untouched
        (0.10, 0.30, 0.0, 1.0e-4, tf + 5.0),    # warm, evap > liqcan: remainder from snow
        (0.10, 0.30, 0.0, 1.1111e-4, tf + 5.0), # warm, evap = whole store
        (0.40, 0.30, 0.0, 5.0e-5, tf - 3.0),    # cold, evap < snocan: liquid untouched
        (0.40, 0.10, 0.0, 1.0e-4, tf - 3.0),    # cold, evap > snocan: remainder from liquid
        (0.20, 0.10, 0.0, -2.0e-5, tf + 1.0),   # warm dew -> liqcan
        (0.20, 0.10, 0.0, -2.0e-5, tf - 1.0),   # frost   -> snocan
        (0.20, 0.10, 0.0, 2.0e-5, tf),          # exactly tfrz: Fortran's cold branch
        (0.00, 0.00, 1.0e-5, 1.0e-5, tf + 2.0), # transpiration only: no canopy change
        (0.30, 0.30, 0.0, 0.0, tf + 2.0),
    ]
    for (liq, sno, tran, evap, tv) in cases
        lf, sf = fortran_dew_update(liq, sno, tran, evap, tv, dt)
        lj, sj = CLM.canopy_water_evap_update(liq, sno, (tran - evap) * dt, tv)
        @test lj ≈ lf atol = 1e-15
        @test sj ≈ sf atol = 1e-15
    end

    # The defect: a warm leaf under a snow-laden canopy. Fortran evaporates the
    # liquid and keeps the snow; the old proportional re-split took 3/7 of the
    # loss from the snow.
    l, s = CLM.canopy_water_evap_update(0.4, 0.3, -0.2, tf + 4.0)
    @test l ≈ 0.2
    @test s == 0.3

    # truncate_small_values: a snocan driven to a relative 1e-10 residue is zeroed.
    _, s = CLM.canopy_water_evap_update(0.0, 0.3, -0.3 * (1 - 1e-12), tf - 1.0)
    @test s == 0.0

    # Float32 stays Float32 (GPU element type).
    l32, s32 = CLM.canopy_water_evap_update(0.4f0, 0.3f0, -0.2f0, Float32(tf + 4))
    @test l32 isa Float32 && s32 isa Float32
end
