# Guards that the decomposed canopy forward (canopy_rev_forward!, the basis of the
# compositional reverse pass in src/biogeophys/canopy_fluxes_reverse.jl) stays
# byte-faithful to canopy_fluxes_core! on a converged single patch. NO Enzyme here —
# this is a fast, suite-safe primal-parity guard against the decomposition drifting
# from the production forward. The reverse-AD gradient itself is validated in
# scripts/enzyme_canopy_psn.jl (Enzyme; run on Julia 1.10), and the full canopy
# reverse + canopy→soil_temperature two-module reverse CHAIN (both through the one
# CLM.compositional_reverse! engine) in scripts/enzyme_driver_reverse.jl.

@testset "canopy_fluxes compositional-reverse forward parity" begin
    FT = Float64
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    NP = 1; NC = 1; NG = 1; NL = 1; NLEVCAN = CLM.NLEVCAN; MP = CLM.MXPFT + 1

    CLM.soil_resistance_read_nl!(soil_resis_method = CLM.SOIL_RESIS_LEEPIELKE_1992)
    CLM.canopy_fluxes_read_nml!(use_undercanopy_stability = false,
        use_biomass_heat_storage = false, itmax_canopy_fluxes = 40)
    CLM.canopy_fluxes_read_params!()

    # photosynthesis scalar params (else photosynthesis! is NaN)
    CLM.photo_params_init!(CLM.params_inst)
    let p = CLM.params_inst
        p.theta_cj = fill(0.98, MP)
        p.theta_ip = 0.95; p.act25 = 72.0; p.fnr = 7.16; p.cp25_yr2000 = 42.75e-6
        p.kc25_coef = 404.9e-6; p.ko25_coef = 278.4e-3; p.fnps = 0.15; p.theta_psii = 0.7
        p.vcmaxha = 65330.0; p.jmaxha = 43540.0; p.tpuha = 53100.0; p.lmrha = 46390.0
        p.kcha = 79430.0; p.koha = 36380.0; p.cpha = 37830.0
        p.vcmaxhd = 149250.0; p.jmaxhd = 152040.0; p.tpuhd = 150650.0; p.lmrhd = 150650.0
        p.lmrse = 490.0; p.tpu25ratio = 0.167; p.kp25ratio = 20160.0
        p.vcmaxse_sf = 1.0; p.jmaxse_sf = 1.0; p.tpuse_sf = 1.0; p.jmax25top_sf = 1.0
    end

    psn = (; c3psn = fill(FT(1.0), MP), leafcn = fill(FT(25.0), MP), flnr = fill(FT(0.1), MP),
        fnitr = fill(FT(1.0), MP), slatop = fill(FT(0.01), MP), mbbopt = fill(FT(9.0), MP),
        medlynintercept = fill(FT(100.0), MP), medlynslope = fill(FT(6.0), MP),
        nrad = fill(1, NP), tlai_z = fill(FT(1.0), NP, NLEVCAN),
        parsun_z = fill(FT(250.0), NP, NLEVCAN), parsha_z = fill(FT(120.0), NP, NLEVCAN),
        laisun_z = fill(FT(1.0), NP, NLEVCAN), laisha_z = fill(FT(1.0), NP, NLEVCAN),
        vcmaxcint_sun = fill(FT(1.0), NP), vcmaxcint_sha = fill(FT(0.6), NP),
        o3coefv = fill(FT(1.0), NP), o3coefg = fill(FT(1.0), NP), t10 = fill(FT(290.0), NP))

    function build()
        canopystate = CLM.CanopyStateData{FT}(); CLM.canopystate_init!(canopystate, NP)
        canopystate.elai_patch[1] = 2.0;  canopystate.esai_patch[1] = 0.5
        canopystate.laisun_patch[1] = 1.2; canopystate.laisha_patch[1] = 0.8
        canopystate.displa_patch[1] = 5.0; canopystate.htop_patch[1] = 10.0
        # Bare-canopy momentum z0 (z0mr*htop): canopy_fluxes! re-seeds z0mv from this
        # each step, so set it when bypassing the preflux z0m/displa kernel.
        canopystate.z0m_patch[1] = 0.5
        canopystate.frac_veg_nosno_patch[1] = 1; canopystate.dleaf_patch[1] = 0.04
        canopystate.stem_biomass_patch[1] = 0.0; canopystate.leaf_biomass_patch[1] = 0.0
        canopystate.tlai_patch[1] = 2.0
        canopystate.laisun_z_patch = copy(psn.laisun_z); canopystate.laisha_z_patch = copy(psn.laisha_z)

        energyflux = CLM.EnergyFluxData{FT}(); CLM.energyflux_init!(energyflux, NP, NC, NL, NG)
        energyflux.btran_patch[1] = 0.8; energyflux.bsun_patch[1] = 0.8; energyflux.bsha_patch[1] = 0.8
        energyflux.htvp_col[1] = CLM.HVAP

        frictionvel = CLM.FrictionVelocityData{FT}(); CLM.frictionvel_init!(frictionvel, NP, NC)
        frictionvel.zetamaxstable = 0.5; frictionvel.zsno = 0.00085; frictionvel.zlnd = 0.000775
        frictionvel.z0mv_patch[1] = 0.5; frictionvel.z0hv_patch[1] = 0.5; frictionvel.z0qv_patch[1] = 0.5
        frictionvel.z0mg_col[1] = 0.01; frictionvel.z0hg_col[1] = 0.01; frictionvel.z0qg_col[1] = 0.01
        frictionvel.forc_hgt_u_patch[1] = 30.0; frictionvel.forc_hgt_t_patch[1] = 30.0
        frictionvel.forc_hgt_q_patch[1] = 30.0
        frictionvel.ustar_patch[1] = 0.5; frictionvel.um_patch[1] = 5.0; frictionvel.uaf_patch[1] = 3.0
        frictionvel.taf_patch[1] = 290.0; frictionvel.qaf_patch[1] = 0.008
        frictionvel.obu_patch[1] = -100.0; frictionvel.zeta_patch[1] = -0.1
        frictionvel.vpd_patch[1] = 1.0; frictionvel.rb1_patch[1] = 50.0
        frictionvel.ram1_patch[1] = 50.0; frictionvel.num_iter_patch[1] = 0.0

        temperature = CLM.TemperatureData{FT}(); CLM.temperature_init!(temperature, NP, NC, NL, NG)
        temperature.t_veg_patch[1] = 290.0; temperature.t_stem_patch[1] = 289.0
        temperature.t_skin_patch[1] = 290.0; temperature.thm_patch[1] = 290.0
        temperature.t_grnd_col[1] = 288.0; temperature.t_h2osfc_col[1] = 288.0
        temperature.thv_col[1] = 291.0; temperature.emv_patch[1] = 0.97; temperature.emg_col[1] = 0.96
        temperature.t_ref2m_patch[1] = 290.0; temperature.t_ref2m_r_patch[1] = 290.0
        temperature.t_a10_patch[1] = 290.0
        for j in 1:size(temperature.t_soisno_col, 2); temperature.t_soisno_col[1, j] = 288.0; end

        solarabs = CLM.SolarAbsorbedData{FT}(); CLM.solarabs_init!(solarabs, NP, NL)
        solarabs.sabv_patch[1] = 150.0
        solarabs.parsun_z_patch = copy(psn.parsun_z); solarabs.parsha_z_patch = copy(psn.parsha_z)

        soilstate = CLM.SoilStateData{FT}(); CLM.soilstate_init!(soilstate, NC, NL)
        soilstate.soilbeta_col[1] = 0.8; soilstate.soilresis_col[1] = 100.0
        for j in 1:CLM.varpar.nlevgrnd; soilstate.watsat_col[1, j] = 0.45; end

        waterfluxbulk = CLM.WaterFluxBulkData{FT}(); CLM.waterfluxbulk_init!(waterfluxbulk, NC, NP, NL, NG)
        waterstatebulk = CLM.WaterStateBulkData{FT}(); CLM.waterstatebulk_init!(waterstatebulk, NC, NP, NL, NG)
        waterstatebulk.ws.liqcan_patch[1] = 0.1; waterstatebulk.ws.snocan_patch[1] = 0.0

        waterdiagbulk = CLM.WaterDiagnosticBulkData{FT}(); CLM.waterdiagnosticbulk_init!(waterdiagbulk, NC, NP, NL, NG)
        waterdiagbulk.fwet_patch[1] = 0.1; waterdiagbulk.fdry_patch[1] = 0.8
        waterdiagbulk.frac_sno_eff_col[1] = 0.0; waterdiagbulk.frac_h2osfc_col[1] = 0.0
        waterdiagbulk.snow_depth_col[1] = 0.0; waterdiagbulk.qg_col[1] = 0.005
        waterdiagbulk.qg_snow_col[1] = 0.005; waterdiagbulk.qg_soil_col[1] = 0.005
        waterdiagbulk.qg_h2osfc_col[1] = 0.005; waterdiagbulk.dqgdT_col[1] = 0.0003
        waterdiagbulk.rh_af_patch[1] = 0.6

        photosyns = CLM.PhotosynthesisData{FT}(); CLM.photosynthesis_data_init!(photosyns, NP)
        photosyns.rssun_patch[1] = 200.0; photosyns.rssha_patch[1] = 300.0
        photosyns.stomatalcond_mtd = CLM.STOMATALCOND_MTD_BB1987
        photosyns.leafresp_method  = CLM.LEAFRESP_MTD_RYAN1991
        photosyns.light_inhibit    = false

        patch_data = CLM.PatchData{FT}(); CLM.patch_init!(patch_data, NP)
        patch_data.column[1] = 1; patch_data.gridcell[1] = 1; patch_data.landunit[1] = 1
        patch_data.itype[1] = 1; patch_data.active[1] = true
        col_data = CLM.ColumnData{FT}(); CLM.column_init!(col_data, NC); col_data.snl[1] = 0
        gridcell_data = CLM.GridcellData{FT}(); CLM.gridcell_init!(gridcell_data, NG)

        S = (; canopystate, energyflux, frictionvel, temperature, solarabs, soilstate,
               waterfluxbulk, waterstatebulk, waterdiagbulk, photosyns, patch_data,
               col_data, gridcell_data)
        forc = (; lwrad = FT[350.0], q = FT[0.008], pbot = FT[101325.0], th = FT[290.0],
                  rho = FT[1.2], t = FT[288.0], u_grc = FT[3.0], v_grc = FT[1.0],
                  pco2 = FT[40.0], po2 = FT[21000.0], hgt_t = FT[10.0], hgt_u = FT[10.0],
                  hgt_q = FT[10.0], dayl = FT[43200.0], max_dayl = FT[50000.0],
                  downreg = FT[1.0], leafn = FT[1.0])
        pft = (; dleaf = fill(FT(0.04), MP), z0v_Cr = fill(FT(0.35), MP),
                 z0v_Cs = fill(FT(0.003), MP), z0v_c = fill(FT(0.25), MP),
                 z0v_cw = fill(FT(2.0), MP), z0v_LAImax = fill(FT(8.0), MP),
                 grnd_ch4 = fill(FT(0.0), NP))
        return (; S, forc, pft, dtime = FT(1800.0))
    end

    # reference: real canopy_fluxes! with photosynthesis enabled
    Bref = build()
    let S = Bref.S, forc = Bref.forc, pft = Bref.pft
        CLM.canopy_fluxes!(S.canopystate, S.energyflux, S.frictionvel, S.temperature,
            S.solarabs, S.soilstate, S.waterfluxbulk, S.waterstatebulk, S.waterdiagbulk,
            S.photosyns, S.patch_data, S.col_data, S.gridcell_data, trues(NP), 1:NP, 1:NC,
            forc.lwrad, forc.q, forc.pbot, forc.th, forc.rho, forc.t, forc.u_grc, forc.v_grc,
            forc.pco2, forc.po2, forc.hgt_t, forc.hgt_u, forc.hgt_q, forc.dayl, forc.max_dayl,
            forc.downreg, forc.leafn, Bref.dtime;
            dleaf_pft = pft.dleaf, z0v_Cr_pft = pft.z0v_Cr, z0v_Cs_pft = pft.z0v_Cs,
            z0v_c_pft = pft.z0v_c, z0v_cw_pft = pft.z0v_cw, z0v_LAImax_pft = pft.z0v_LAImax,
            grnd_ch4_cond_patch = pft.grnd_ch4,
            c3psn_pft = psn.c3psn, leafcn_pft = psn.leafcn, flnr_pft = psn.flnr,
            fnitr_pft = psn.fnitr, slatop_pft = psn.slatop, mbbopt_pft = psn.mbbopt,
            medlynintercept_pft = psn.medlynintercept, medlynslope_pft = psn.medlynslope,
            t10_patch = psn.t10, nrad_patch = psn.nrad, tlai_z_patch = psn.tlai_z,
            vcmaxcint_sun_patch = psn.vcmaxcint_sun, vcmaxcint_sha_patch = psn.vcmaxcint_sha,
            parsun_z_patch = psn.parsun_z, parsha_z_patch = psn.parsha_z,
            laisun_z_patch = psn.laisun_z, laisha_z_patch = psn.laisha_z,
            o3coefv_patch = psn.o3coefv, o3coefg_patch = psn.o3coefg)
    end
    nref = Int(Bref.S.frictionvel.num_iter_patch[1])

    # decomposed forward via the src module
    B = build()
    aux = (; patch = B.S.patch_data, col = B.S.col_data, grid = B.S.gridcell_data,
        forc = B.forc, pft = B.pft, psn = psn, filterp = Int[1], fn = 1, active = trues(NP),
        mask = trues(NP), ivt = B.S.patch_data.itype .+ 1,
        forc_pbot_patch = FT[B.forc.pbot[B.S.patch_data.column[p]] for p in 1:NP],
        soilevap_beta = CLM.do_soilevap_beta(), soil_resis_sl14 = CLM.do_soil_resistance_sl14(),
        nlevsno = CLM.varpar.nlevsno, dtime = B.dtime, use_psn = true)
    b = CLM.cf_rev_bundle(B.S.canopystate, B.S.energyflux, B.S.frictionvel, B.S.temperature,
        B.S.solarabs, B.S.soilstate, B.S.waterfluxbulk, B.S.waterstatebulk, B.S.waterdiagbulk,
        B.S.photosyns, CLM.cf_rev_scratch(FT, NP))
    CLM.canopy_rev_forward!(b, aux, max(nref, 1))

    @test isfinite(b.temperature.t_veg_patch[1])
    @test b.temperature.t_veg_patch[1] ≈ Bref.S.temperature.t_veg_patch[1] atol = 1e-9
    @test b.photosyns.rssun_patch[1] ≈ Bref.S.photosyns.rssun_patch[1] rtol = 1e-9
end
