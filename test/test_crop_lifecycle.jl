using Test
using CLM

# ============================================================================
# Crop lifecycle — validated against the CTSM CN-crop Fortran reference.
#
# Reference: SYMFLUENCE_data/clm_bgc_spinup/crop_ref_usplains/
#   refs_crop_may/bgcdump_after_fire_n38665.nc   (fertilization window OPEN)
#   refs_crop_july/bgcdump_after_fire_n39865.nc  (window CLOSED, FERT == 0)
#   Crop_USplains.clm2.r.2020-05-30-14400.nc     (restart: carries idop,
#                                                 gddmaturity, huileaf, huigrain,
#                                                 fert_counter, ...)
# See docs/CROP_LIFECYCLE_MAP.md for the dependency ordering.
#
# WHAT THIS VALIDATES: the derivation chain that the lifecycle produces —
#   gddmaturity -> huileaf / huigrain -> fertilizer flux
# against the reference's own numbers, plus the fact that `croplive` actually
# flips (it never could before: plant_crop! had no caller in the entire port).
#
# WHAT THIS DOES *NOT* VALIDATE, stated explicitly so a green run is not
# over-read:
#   * The planting DATE selection (reference idop = 136/127/145/106). Reproducing
#     it needs the run's 20-year GDD climatology (gdd020/gdd820/gdd1020 are
#     20-yr running means), which is not reconstructible from the dumps. Here
#     the planting conditions are synthesized; the date is an input, not an
#     output under test.
#   * vernalization!. The reference has NO winter-wheat CFT (no itype 21/22) and
#     its cumvd/hdidx are NaN on every patch. That path is wired but UNVALIDATED.
#   * n_soyfix!. SOYFIXN == 0 in BOTH reference windows, so there is no non-zero
#     ground truth; it is deliberately left unwired.
# ============================================================================

@testset "Crop lifecycle vs Fortran crop reference" begin

    # --- Reference: the 8 live crop patches (0-based dump index 5..12) -------
    REF_ITYPE       = [17, 18, 19, 20, 23, 24, 75, 76]
    # itype 17/18 temperate corn, 19/20 spring wheat, 23/24 temperate soybean,
    # 75/76 TROPICAL corn (confirmed from clm5_params.nc pftname) — 75/76 are in
    # the corn family, which is what selects the crmcorn huigrain branch.
    # FULL PRECISION, straight out of the reference files (not rounded — rounding
    # these is what forces a loose tolerance and hides real error).
    REF_GDDMATURITY = [1477.2873853187507, 1476.7794835074303, 1700.0, 1700.0,
                       1259.2610846357402, 1257.359842574543,
                       1334.9881712188765, 1334.2627384206983]
    REF_HUILEAF     = [44.318621559562516, 44.30338450522291, 85.0, 85.0,
                       37.77783253907221, 37.72079527723629,
                       40.04964513656629, 40.02788215262095]
    REF_HUIGRAIN    = [850.0767538833009, 849.8925535972276, 1020.0, 1020.0,
                       629.6305423178701, 628.6799212872716,
                       595.3141022625068, 595.1300578305277]
    REF_FERTNITRO   = [14.748827666355744, 14.748827666355744,
                       7.896248057624567, 7.98679573715093,
                       0.7075451110211234, 0.7075451110211234,
                       14.748827666355742, 14.748827666355744]
    REF_FERT        = [9.69260860321513e-6, 9.69260860321513e-6,
                       5.726995403717921e-6, 5.779395681221604e-6,
                       1.5668663836927799e-6, 1.5668663836927799e-6,
                       9.69260860321513e-6, 9.69260860321513e-6]

    # --- Real parameter values, read from the paramfile the reference ran with
    # (clm5_params.nc, pft = 79). NOT invented: these reproduce the reference's
    # huileaf/huigrain/FERT exactly, which is the point of the test.
    P_LFEMERG   = [0.03, 0.03, 0.05, 0.05, 0.03, 0.03, 0.03, 0.03]
    P_GRNFILL   = [0.65, 0.65, 0.60, 0.60, 0.50, 0.50, 0.50, 0.50]
    P_HYBGDD    = [1700.0, 1700.0, 1700.0, 1700.0, 1900.0, 1900.0, 1800.0, 1800.0]
    P_MXMAT     = [165, 165, 150, 150, 150, 150, 160, 160]
    P_LEAFCN    = [25.0, 25.0, 20.0, 20.0, 20.0, 20.0, 25.0, 25.0]
    P_MANUNITRO = fill(0.002, 8)
    P_GDDMIN    = fill(50.0, 8)

    NPFT = 79
    np, nc, ng = 8, 1, 1
    DT = 1800.0
    JDAY_PLANT = 136

    # gdd20 inputs chosen to INVERT to the reference gddmaturity through CTSM's
    # own formulas — so gddmaturity is a genuine output under test, not an input.
    #   corn family : gddmaturity = max(950, min(gdd820*0.85, hybgdd)) (+150 if normal)
    #   soybean     : gddmaturity = min(gdd1020, hybgdd)
    #   spring wheat: gddmaturity = min(gdd020,  hybgdd)
    # Inverted from REF_GDDMATURITY rather than hardcoded, so the inversion is
    # explicit and the forward direction stays a real test.
    # 100.0 entries are "irrelevant but >= gddmin(50)" — the planting-decision GDD
    # gate uses gdd820 for every non-winter-cereal crop even when gddmaturity
    # keys off gdd020/gdd1020 (CNPhenologyMod.F90:2279).
    GDD820  = [(REF_GDDMATURITY[1] - 150.0) / 0.85, (REF_GDDMATURITY[2] - 150.0) / 0.85,
               100.0, 100.0, 100.0, 100.0,
               (REF_GDDMATURITY[7] - 150.0) / 0.85, (REF_GDDMATURITY[8] - 150.0) / 0.85]
    GDD1020 = [0.0, 0.0, 0.0, 0.0, REF_GDDMATURITY[5], REF_GDDMATURITY[6], 0.0, 0.0]
    GDD020  = [0.0, 0.0, 1800.0, 1800.0, 0.0, 0.0, 0.0, 0.0]  # min(1800,1700)=1700

    function make_crop_state()
        pstate = CLM.PhenologyState()
        params = CLM.PhenologyParams()
        CLM.cn_phenology_init!(pstate, params, DT)

        pft = CLM.PftConPhenology(
            evergreen = zeros(NPFT), season_decid = zeros(NPFT),
            season_decid_temperate = zeros(NPFT), stress_decid = zeros(NPFT),
            woody = zeros(NPFT), leaf_long = fill(1.0, NPFT),
            leafcn = fill(25.0, NPFT), frootcn = fill(42.0, NPFT),
            lflitcn = fill(50.0, NPFT), livewdcn = fill(50.0, NPFT),
            deadwdcn = fill(500.0, NPFT), ndays_on = fill(30.0, NPFT),
            crit_onset_gdd_sf = fill(1.0, NPFT),
            lf_f = fill(1.0/3, NPFT, 3), fr_f = fill(1.0/3, NPFT, 3),
            biofuel_harvfrac = zeros(NPFT), repr_structure_harvfrac = zeros(NPFT, 1),
            minplanttemp = fill(280.0, NPFT), planttemp = fill(283.0, NPFT),
            gddmin = fill(50.0, NPFT), lfemerg = zeros(NPFT), grnfill = zeros(NPFT),
            hybgdd = zeros(NPFT), mxmat = fill(200, NPFT), manunitro = zeros(NPFT),
            is_pft_known_to_model = fill(true, NPFT),
            mnNHplantdate = fill(100.0, NPFT), mxNHplantdate = fill(170.0, NPFT),
            mnSHplantdate = fill(280.0, NPFT), mxSHplantdate = fill(350.0, NPFT),
        )
        # Install the REAL per-crop params at their 1-based Julia slots.
        for (k, t) in enumerate(REF_ITYPE)
            j = t + 1
            pft.lfemerg[j]   = P_LFEMERG[k]
            pft.grnfill[j]   = P_GRNFILL[k]
            pft.hybgdd[j]    = P_HYBGDD[k]
            pft.mxmat[j]     = P_MXMAT[k]
            pft.leafcn[j]    = P_LEAFCN[k]
            pft.manunitro[j] = P_MANUNITRO[k]
            pft.gddmin[j]    = P_GDDMIN[k]
        end

        patch = CLM.PatchData()
        patch.itype    = copy(REF_ITYPE)
        patch.column   = fill(1, np)
        patch.gridcell = fill(1, np)
        patch.wtcol    = fill(1.0, np)

        grc = CLM.GridcellData()
        grc.latdeg = fill(44.8, ng)          # the reference point (US Great Plains)
        grc.dayl = fill(45000.0, ng); grc.prev_dayl = fill(44000.0, ng)

        temp = CLM.TemperatureData()
        temp.t_ref2m_patch     = fill(290.0, np)
        temp.t_ref2m_min_patch = fill(285.0, np)
        temp.t_ref2m_max_patch = fill(295.0, np)
        temp.t_a10_patch       = fill(290.0, np)   # > planttemp(283)
        temp.t_a10min_patch    = fill(285.0, np)   # > minplanttemp(280)
        temp.t_a5min_patch     = fill(285.0, np)
        temp.gdd020_patch      = copy(GDD020)
        temp.gdd820_patch      = copy(GDD820)
        temp.gdd1020_patch     = copy(GDD1020)

        wdb = CLM.WaterDiagnosticBulkData()
        wdb.snow_depth_col = fill(0.0, nc); wdb.snow_5day_col = fill(0.0, nc)

        can = CLM.CanopyStateData(); can.tlai_patch = fill(1.0, np)

        vs = CLM.CNVegStateData(); CLM.cnveg_state_init!(vs, np, nc)
        cs = CLM.CNVegCarbonStateData(); CLM.cnveg_carbon_state_init!(cs, np, nc, ng)
        ns = CLM.CNVegNitrogenStateData(); CLM.cnveg_nitrogen_state_init!(ns, np, nc, ng)
        cf = CLM.CNVegCarbonFluxData(); CLM.cnveg_carbon_flux_init!(cf, np, nc, ng)
        nf = CLM.CNVegNitrogenFluxData(); CLM.cnveg_nitrogen_flux_init!(nf, np, nc, ng)

        cr = CLM.CropData(); CLM.crop_init!(cr, np)
        cr.fertnitro_patch .= REF_FERTNITRO
        cr.croplive_patch  .= false
        cr.sowing_count    .= 0
        cr.harvest_count   .= 0
        cr.hui_patch       .= 0.0
        cr.gddtsoi_patch   .= 0.0
        cr.gddaccum_patch  .= 0.0
        cr.vf_patch        .= 0.0
        cr.sown_in_this_window .= false
        cr.rx_swindow_starts_thisyr_patch .= -1
        cr.rx_swindow_ends_thisyr_patch   .= -1

        vs.idop_patch .= 0
        vs.onset_counter_patch .= 0.0
        vs.peaklai_patch .= 0
        nf.fert_patch .= 0.0
        nf.fert_counter_patch .= 0.0

        # THE call that had no live call site anywhere in the port.
        CLM.crop_phenology_init!(pstate, pft, patch, grc, 17, NPFT - 1, NPFT - 1)

        return (; pstate, params, pft, patch, grc, temp, wdb, can, vs, cs, ns, cf, nf, cr)
    end

    step!(d, jday) = CLM.crop_phenology!(
        d.pstate, d.params, trues(np), d.pft, d.wdb, d.temp, d.cr, d.can,
        d.vs, d.cs, d.ns, d.cf, d.nf, d.patch, d.grc;
        jday = jday, kyr = 2020, dayspyr = 366.0, use_fertilizer = true,
        avg_dayspyr = 365.0)

    # =====================================================================
    # S0 — the planting-window tables are real (they were EMPTY: no caller)
    # =====================================================================
    @testset "S0 sowing window tables are populated" begin
        d = make_crop_state()
        @test !isempty(d.pstate.minplantjday)
        @test !isempty(d.pstate.inhemi)
        @test all(d.pstate.inhemi .== CLM.inNH)          # lat 44.8 > 0
        for t in REF_ITYPE
            @test d.pstate.minplantjday[t + 1, CLM.inNH] == 100
            @test d.pstate.maxplantjday[t + 1, CLM.inNH] == 170
        end
        # A non-crop PFT must NOT have been given a planting window.
        @test d.pstate.minplantjday[5, CLM.inNH] == typemax(Int)
    end

    # =====================================================================
    # S1 — planting fires and gddmaturity matches the reference
    # =====================================================================
    @testset "S1 planting + gddmaturity" begin
        d = make_crop_state()
        @test !any(d.cr.croplive_patch)     # dead before
        step!(d, JDAY_PLANT)

        # croplive flipping AT ALL is the structural fact that was impossible
        # before this change (plant_crop! had zero callers).
        @test all(d.cr.croplive_patch)
        @test all(d.vs.idop_patch .== JDAY_PLANT)
        @test all(d.vs.iyop_patch .== 2020)
        @test all(d.cr.harvdate_patch .== CLM.NOT_Harvested)
        @test all(d.cr.sowing_count .== 1)

        for k in 1:np
            @test d.vs.gddmaturity_patch[k] ≈ REF_GDDMATURITY[k] rtol=1e-12
        end

        # The seed transfer really happened (non-no-op body, not just a flag).
        seed = CLM._initial_seed_at_planting[]
        @test all(d.cs.leafc_xfer_patch .≈ seed)
        for k in 1:np
            @test d.ns.leafn_xfer_patch[k] ≈ seed / P_LEAFCN[k]
        end
    end

    # =====================================================================
    # S2 — phase thresholds, incl. the corn-family vs plain huigrain BRANCH
    # =====================================================================
    @testset "S2 huileaf / huigrain vs reference" begin
        d = make_crop_state()
        step!(d, JDAY_PLANT)

        # Full-precision reference values -> demand round-off agreement, not a
        # loose band. Anything looser would not distinguish the huigrain branches
        # from a mis-tuned coefficient.
        for k in 1:np
            @test d.vs.huileaf_patch[k]  ≈ REF_HUILEAF[k]  rtol=1e-12
            @test d.vs.huigrain_patch[k] ≈ REF_HUIGRAIN[k] rtol=1e-12
        end

        # Branch discrimination, asserted directly: corn family (17/18/75/76)
        # must NOT equal grnfill*gddmaturity, everything else MUST.
        for k in 1:np
            plain = P_GRNFILL[k] * d.vs.gddmaturity_patch[k]
            if REF_ITYPE[k] in (17, 18, 75, 76)
                @test !isapprox(d.vs.huigrain_patch[k], plain; rtol=1e-3)
            else
                @test d.vs.huigrain_patch[k] ≈ plain rtol=1e-12
            end
        end
    end

    # =====================================================================
    # S3/S4 — leaf emergence onset opens the fertilizer window; FERT value
    # =====================================================================
    @testset "S3/S4 leaf emergence + FERT vs May reference" begin
        d = make_crop_state()
        step!(d, JDAY_PLANT)
        @test all(d.nf.fert_patch .== 0.0)      # not yet emerged -> no fertilizer

        # Emerge: gddtsoi >= huileaf, hui still below huigrain.
        d.cr.gddtsoi_patch .= d.vs.huileaf_patch .+ 1.0
        step!(d, JDAY_PLANT + 1)

        @test all(d.cr.cphase_patch .== CLM.cphase_leafemerge)
        @test all(d.vs.onset_flag_patch .== 1.0)
        # fert_counter is set to 20 days then decremented once in the same step.
        @test all(d.nf.fert_counter_patch .≈ (20.0 * CLM.SECSPDAY - DT))

        # THE value check. fert = (manunitro*1000 + fertnitro)/(20*86400).
        # With full-precision FERTNITRO this reproduces the reference FERT to
        # relative error 0.0 on all 8 patches, so it is asserted at round-off.
        for k in 1:np
            @test d.nf.fert_patch[k] ≈ REF_FERT[k] rtol=1e-13
        end
        @test all(d.nf.fert_patch .> 0.0)
    end

    # =====================================================================
    # S4b — the July signature: the window CLOSES and FERT returns to 0
    # =====================================================================
    @testset "S4b fertilizer window closes (July reference: FERT == 0)" begin
        d = make_crop_state()
        step!(d, JDAY_PLANT)
        d.cr.gddtsoi_patch .= d.vs.huileaf_patch .+ 1.0
        step!(d, JDAY_PLANT + 1)
        @test all(d.nf.fert_patch .> 0.0)

        # Run the counter out — 20 days of 1800 s steps.
        d.nf.fert_counter_patch .= 0.0
        step!(d, JDAY_PLANT + 2)
        @test all(d.nf.fert_patch .== 0.0)   # matches July window: FERT identically 0
        # ...and the crop is still alive: a zero FERT here is a CLOSED window,
        # not a dead path. This is the distinction the July-only dump hides.
        @test all(d.cr.croplive_patch)
    end

    # =====================================================================
    # S6 — harvest at maturity
    # =====================================================================
    @testset "S6 harvest at gddmaturity" begin
        d = make_crop_state()
        step!(d, JDAY_PLANT)
        @test all(d.cr.croplive_patch)

        # Push hui past gddmaturity -> harvest.
        d.cr.hui_patch .= d.vs.gddmaturity_patch .+ 1.0
        step!(d, JDAY_PLANT + 10)

        @test !any(d.cr.croplive_patch)
        @test all(d.cr.cphase_patch .== CLM.cphase_harvest)
        @test all(d.cr.harvdate_patch .== JDAY_PLANT + 10)
        @test all(d.cr.harvest_count .== 1)
        # tlai > 0 so the offset (litterfall) path must be taken.
        @test all(d.vs.offset_flag_patch .== 1.0)
        for k in 1:np
            @test d.cr.harvest_reason_thisyr_patch[k, 1] == CLM.HARVEST_REASON_MATURE
            @test d.cr.hui_thisyr_patch[k, 1] ≈ d.vs.gddmaturity_patch[k] + 1.0
        end
    end

    # =====================================================================
    # Guard: the init assertion actually fires (dead-wiring cannot return)
    # =====================================================================
    @testset "crop_phenology! refuses to run un-initialized" begin
        d = make_crop_state()
        d.pstate.minplantjday = Matrix{Int}(undef, 0, 0)
        @test_throws ErrorException step!(d, JDAY_PLANT)
    end
end
