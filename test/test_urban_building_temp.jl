using ForwardDiff

@testset "UrbanBuildingTemperature (prognostic, Oleson2015)" begin
    # ---------------------------------------------------------------------
    # Minimal single-urban-landunit setup. One urban landunit (l=1) with
    # three columns: roof (c=1), sunwall (c=2), shadewall (c=3). We build only
    # the fields building_temperature! reads, then exercise the prognostic
    # 5x5 building energy-balance solve directly.
    # ---------------------------------------------------------------------
    vp = CLM.varpar
    saved_nlevsno        = vp.nlevsno
    saved_nlevgrnd       = vp.nlevgrnd
    saved_nlevurb        = vp.nlevurb
    saved_nlevmaxurbgrnd = vp.nlevmaxurbgrnd
    saved_nlevsoi        = vp.nlevsoi

    vp.nlevsno        = 5
    vp.nlevgrnd       = 10
    vp.nlevurb        = 5
    vp.nlevmaxurbgrnd = 10
    vp.nlevsoi        = 8

    nlevsno = vp.nlevsno
    nlevurb = vp.nlevurb
    joff    = nlevsno
    nlevurb_jj = nlevurb + joff

    # Urban control: prognostic method, HAC on so the HAC clamp branch runs.
    saved_read_nml             = CLM.urban_ctrl.read_namelist
    saved_method               = CLM.urban_ctrl.building_temp_method
    saved_hac                  = CLM.urban_ctrl.urban_hac
    saved_explicit             = CLM.urban_ctrl.urban_explicit_ac
    CLM.urban_ctrl.read_namelist       = true
    CLM.urban_ctrl.building_temp_method = CLM.BUILDING_TEMP_METHOD_PROG
    CLM.urban_ctrl.urban_hac            = CLM.URBAN_HAC_ON
    CLM.urban_ctrl.urban_explicit_ac    = true

    nc = 3   # roof, sunwall, shadewall
    nl = 1
    ng = 1
    np = 1

    function setup(::Type{T}=Float64) where {T<:Real}
        col = CLM.ColumnData{T}();       CLM.column_init!(col, nc)
        lun = CLM.LandunitData{T}();     CLM.landunit_init!(lun, nl)
        temp = CLM.TemperatureData{T}(); CLM.temperature_init!(temp, np, nc, nl, ng)
        ef  = CLM.EnergyFluxData{T}();   CLM.energyflux_init!(ef, np, nc, nl, ng)
        up  = CLM.UrbanParamsData{T}();  CLM.urbanparams_init!(up, nl; nlevurb=nlevurb)

        # Landunit l=1 is urban
        lun.urbpoi[1]       = true
        lun.itype[1]        = CLM.ISTURB_MD
        lun.canyon_hwr[1]   = 1.0
        lun.wtlunit_roof[1] = 0.5
        lun.ht_roof[1]      = 10.0

        up.t_building_min[1] = 290.0

        # Columns -> landunit 1
        col.landunit .= 1
        col.itype[1] = CLM.ICOL_ROOF
        col.itype[2] = CLM.ICOL_SUNWALL
        col.itype[3] = CLM.ICOL_SHADEWALL
        col.snl     .= 0

        # Geometry at inner (nlevurb) node: z node depth, zi lower interface.
        # building_temperature! reads col.z[c, nlevurb_jj] and col.zi[c, nlevurb_jj+1].
        for c in 1:nc
            col.z[c, nlevurb_jj]       = 0.20   # node depth (m)
            col.zi[c, nlevurb_jj + 1]  = 0.25   # interface below node (m)
        end

        # Temperatures: inner-node soil temps (just-solved) + previous-step.
        for c in 1:nc
            temp.t_soisno_col[c, nlevurb_jj] = 291.0   # t_*_innerl (this step)
            temp.t_ssbef_col[c, nlevurb_jj]  = 290.5   # t_*_innerl_bef (prev step)
        end

        # Landunit interior temperatures (previous step) + canopy air temp.
        temp.t_roof_inner_lun[1] = 292.0
        temp.t_sunw_inner_lun[1] = 292.5
        temp.t_shdw_inner_lun[1] = 291.5
        temp.t_floor_lun[1]      = 293.0
        temp.t_building_lun[1]   = 294.0
        temp.taf_lun[1]          = 300.0   # warm canyon air

        # Interface thermal conductivity at the inner node (W m-1 K-1).
        tk = zeros(T, nc, nlevsno + vp.nlevmaxurbgrnd)
        for c in 1:nc
            tk[c, nlevurb_jj] = 1.5
        end

        urbantv_t_building_max = fill(T(297.0), nl)  # AC setpoint (K)
        urbantv_p_ac           = fill(T(1.0), nl)     # full AC adoption

        return col, lun, temp, ef, up, tk, urbantv_t_building_max, urbantv_p_ac
    end

    mask_urbanl  = trues(nl)
    mask_nolakec = trues(nc)
    dtime = 1800.0

    # --- 1. Finite / sanity ---------------------------------------------
    col, lun, temp, ef, up, tk, tbmax, p_ac = setup()
    CLM.building_temperature!(col, lun, temp, ef, up, tbmax, p_ac, tk,
                              mask_urbanl, mask_nolakec, 1:nl, dtime)

    for x in (temp.t_roof_inner_lun[1], temp.t_sunw_inner_lun[1],
              temp.t_shdw_inner_lun[1], temp.t_floor_lun[1], temp.t_building_lun[1])
        @test isfinite(x)
    end
    @test isfinite(ef.eflx_building_lun[1])
    @test isfinite(ef.eflx_ventilation_lun[1])
    @test isfinite(ef.eflx_urban_ac_lun[1])
    @test isfinite(ef.eflx_urban_heat_lun[1])

    # Inner-surface temperatures should land in a physical range (K).
    for x in (temp.t_roof_inner_lun[1], temp.t_sunw_inner_lun[1],
              temp.t_shdw_inner_lun[1], temp.t_floor_lun[1], temp.t_building_lun[1])
        @test 250.0 < x < 350.0
    end

    # Warm canyon air (300K) > AC setpoint (297K): building air should be cooled
    # toward/at the max setpoint and AC flux should be non-negative.
    @test temp.t_building_lun[1] <= 297.0 + 1.0
    @test ef.eflx_urban_ac_lun[1] >= 0.0

    # --- 2. Energy-balance check of the building-air equation ------------
    # Re-solve and verify the building-air budget closes BEFORE the HAC clamp
    # by turning HAC off, so t_building is the raw 5x5 solution.
    CLM.urban_ctrl.urban_hac = CLM.URBAN_HAC_OFF
    col, lun, temp, ef, up, tk, tbmax, p_ac = setup()
    # Snapshot previous-step values used in the balance.
    trib = temp.t_roof_inner_lun[1]; tsib = temp.t_sunw_inner_lun[1]
    thib = temp.t_shdw_inner_lun[1]; tfb  = temp.t_floor_lun[1]
    tbb  = temp.t_building_lun[1];   taf  = temp.taf_lun[1]
    rho_dair = CLM.PSTD / (CLM.RAIR * tbb)
    bhwr = lun.canyon_hwr[1]*(1.0 - lun.wtlunit_roof[1])/lun.wtlunit_roof[1]
    hcv_roofi  = trib <= tbb ? CLM.HCV_ROOF_ENHANCED  : CLM.HCV_ROOF
    hcv_floori = tfb  >= tbb ? CLM.HCV_FLOOR_ENHANCED : CLM.HCV_FLOOR

    CLM.building_temperature!(col, lun, temp, ef, up, tbmax, p_ac, tk,
                              mask_urbanl, mask_nolakec, 1:nl, dtime)
    tri = temp.t_roof_inner_lun[1]; tsi = temp.t_sunw_inner_lun[1]
    thi = temp.t_shdw_inner_lun[1]; tf  = temp.t_floor_lun[1]
    tbi = temp.t_building_lun[1]

    # Building-air energy balance residual (UrbBuildTempOleson2015Mod, enrgy_bal_buildair).
    bal = (lun.ht_roof[1]*rho_dair*CLM.CPAIR/dtime)*(tbi - tbb) -
          lun.ht_roof[1]*(CLM.VENT_ACH/3600.0)*rho_dair*CLM.CPAIR*(taf - tbi) -
          0.5*hcv_roofi*(tri - tbi) -
          0.5*hcv_roofi*(trib - tbb) -
          0.5*CLM.HCV_SUNW*(tsi - tbi)*bhwr -
          0.5*CLM.HCV_SUNW*(tsib - tbb)*bhwr -
          0.5*CLM.HCV_SHDW*(thi - tbi)*bhwr -
          0.5*CLM.HCV_SHDW*(thib - tbb)*bhwr -
          0.5*hcv_floori*(tf - tbi) -
          0.5*hcv_floori*(tfb - tbb)
    @test abs(bal) < 1.0e-6   # Fortran's hard tolerance is 0.10 W/m2

    # --- 3. FD-derivative check: d t_building / d taf --------------------
    # The 5x5 system is linear in taf (only res[5] depends on it via the
    # ventilation term), so the prognostic t_building is an analytic function
    # of taf. Verify ForwardDiff matches a finite difference.
    CLM.urban_ctrl.urban_hac = CLM.URBAN_HAC_OFF  # avoid the clamp kink
    function tbuilding_of_taf(taf_val::T) where {T<:Real}
        col, lun, temp, ef, up, tk, tbmax, p_ac = setup(T)
        temp.taf_lun[1] = taf_val
        CLM.building_temperature!(col, lun, temp, ef, up, tbmax, p_ac, tk,
                                  mask_urbanl, mask_nolakec, 1:nl, dtime)
        return temp.t_building_lun[1]
    end

    taf0 = 300.0
    d_ad = ForwardDiff.derivative(tbuilding_of_taf, taf0)
    h = 1.0e-4
    d_fd = (tbuilding_of_taf(taf0 + h) - tbuilding_of_taf(taf0 - h)) / (2h)
    @test isfinite(d_ad)
    @test d_ad > 0.0                         # warmer canyon air -> warmer building
    @test isapprox(d_ad, d_fd; rtol=1e-5, atol=1e-8)

    # Restore globals
    CLM.urban_ctrl.read_namelist        = saved_read_nml
    CLM.urban_ctrl.building_temp_method = saved_method
    CLM.urban_ctrl.urban_hac            = saved_hac
    CLM.urban_ctrl.urban_explicit_ac    = saved_explicit
    vp.nlevsno        = saved_nlevsno
    vp.nlevgrnd       = saved_nlevgrnd
    vp.nlevurb        = saved_nlevurb
    vp.nlevmaxurbgrnd = saved_nlevmaxurbgrnd
    vp.nlevsoi        = saved_nlevsoi
end
