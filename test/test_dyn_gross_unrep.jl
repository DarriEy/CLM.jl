# Tests for src/biogeochem/dyn_gross_unrep.jl — the gross unrepresented landcover
# change reader (UNREPRESENTED_PFT_LULCC) and the CN disturbance-flux kernel.
#
# (1) File-reading test: write a temp NetCDF with a YEAR axis and an
#     UNREPRESENTED_PFT_LULCC variable shaped (gridcell, natpft, time), run
#     dynGrossUnrep_init + dynGrossUnrep_interp!, and assert the current-year
#     fraction (and the before-time-series do_grossunrep gating).
#
# (2) CNGrossUnrep flux test: synthetic veg C/N pools + a known gross-unrep
#     fraction, run cn_gross_unrep!, and assert
#       - C/N conservation (every pool loss == its litter + atmosphere +
#         wood-product gain), and
#       - correct patch -> column aggregation (litter / CWD / product columns).

using Test
using NCDatasets

@testset "dyn_gross_unrep" begin

    # ==================================================================
    # (1) File reading: init + interp picks the correct current-year frac
    # ==================================================================
    @testset "dynGrossUnrep_init + interp" begin
        file_years = [2000, 2001, 2002]
        ngrid = 2
        natpft = 4   # natpft_size

        # frac[g, k, t] distinct & checkable: 0.01*g + 0.001*k + 0.0001*(t)
        frac = Array{Float64}(undef, ngrid, natpft, length(file_years))
        for t in 1:length(file_years), k in 1:natpft, g in 1:ngrid
            frac[g, k, t] = 0.01 * g + 0.001 * k + 0.0001 * t
        end

        mktempdir() do dir
            fn = joinpath(dir, "grossunrep_test.nc")
            NCDataset(fn, "c") do ds
                defDim(ds, "gridcell", ngrid)
                defDim(ds, "natpft", natpft)
                defDim(ds, "time", length(file_years))
                defVar(ds, "YEAR", Int, ("time",))[:] = file_years
                defVar(ds, "UNREPRESENTED_PFT_LULCC", Float64,
                       ("gridcell", "natpft", "time"))[:, :, :] = frac
            end

            # --- Init positioned before the axis (1999): do_grossunrep OFF ---
            state = CLM.dynGrossUnrep_init(1, ngrid, natpft, fn; current_year = 1999)
            @test size(state.grossunrepfrac) == (ngrid, natpft)

            CLM.dynGrossUnrep_interp!(state, 1999)
            @test state.do_grossunrep == false
            # Before the series the fraction array is zeroed.
            @test all(state.grossunrepfrac .== 0.0)

            # --- Step to the first file year (2000): ON, year-1 slice ---
            CLM.dynGrossUnrep_interp!(state, 2000)
            @test state.do_grossunrep == true
            @test state.grossunrepfrac ≈ frac[:, :, 1]

            # --- Step mid-axis (2001): year-2 slice ---
            CLM.dynGrossUnrep_interp!(state, 2001)
            @test state.do_grossunrep == true
            @test state.grossunrepfrac ≈ frac[:, :, 2]

            # --- Step past the end (2050): last-year rates maintained, still ON ---
            CLM.dynGrossUnrep_interp!(state, 2050)
            @test state.do_grossunrep == true
            @test state.grossunrepfrac ≈ frac[:, :, 3]

            CLM.dyn_file_close!(state.dyngrossunrep_file)
        end
    end

    # ==================================================================
    # (2) CNGrossUnrep flux: conservation + patch -> column aggregation
    # ==================================================================
    @testset "cn_gross_unrep! conservation + aggregation" begin
        np = 2            # two patches
        nc = 1            # one column (both patches in column 1)
        ngrid = 1
        nlevdecomp = 1
        n_litr = 3
        i_litr_min, i_litr_max, i_met_lit = 1, 3, 1
        natpft = 15       # enough to cover PFT indices up to nc4_grass

        dt = 1800.0
        noveg = 0
        nc4_grass = 14

        # Known gross-unrep fraction per (gridcell, PFT). Both patches are natural
        # PFTs in (noveg, nc4_grass]. PFT axis stored 1-based (PFT k -> col k+1).
        am1 = 0.10   # annual fraction for patch 1 (PFT 2)
        am2 = 0.20   # annual fraction for patch 2 (PFT 7)
        grossunrepfrac = zeros(ngrid, natpft)
        grossunrepfrac[1, 2 + 1] = am1   # PFT 2
        grossunrepfrac[1, 7 + 1] = am2   # PFT 7

        state = CLM.DynGrossUnrepState(
            grossunrepfrac = grossunrepfrac,
            do_grossunrep  = true,
        )

        # convfrac (pconv): fraction of deadstem routed to atmosphere vs product.
        # Distinct per PFT so the wood-product split is exercised.
        pconv = fill(0.0, natpft)
        pconv[2 + 1] = 0.6
        pconv[7 + 1] = 0.25

        pftcon = CLM.PftConGrossUnrep(
            pconv = pconv,
            lf_f  = fill(1.0 / n_litr, natpft, n_litr),
            fr_f  = fill(1.0 / n_litr, natpft, n_litr),
        )

        patch = CLM.PatchData()
        patch.itype    = [2, 7]
        patch.gridcell = [1, 1]
        patch.column   = [1, 1]
        patch.wtcol    = [0.4, 0.6]

        # --- Carbon state pools (distinct values per patch) ---
        cs = CLM.CNVegCarbonStateData()
        cs.leafc_patch              = [10.0, 12.0]
        cs.frootc_patch             = [ 4.0,  5.0]
        cs.livestemc_patch          = [20.0, 22.0]
        cs.deadstemc_patch          = [50.0, 60.0]
        cs.livecrootc_patch         = [ 8.0,  9.0]
        cs.deadcrootc_patch         = [25.0, 30.0]
        cs.xsmrpool_patch           = [ 1.0,  2.0]
        cs.leafc_storage_patch      = [ 1.0,  1.2]
        cs.frootc_storage_patch     = [ 0.5,  0.6]
        cs.livestemc_storage_patch  = [ 2.0,  2.2]
        cs.deadstemc_storage_patch  = [ 3.0,  3.3]
        cs.livecrootc_storage_patch = [ 1.0,  1.1]
        cs.deadcrootc_storage_patch = [ 1.5,  1.6]
        cs.gresp_storage_patch      = [ 0.2,  0.3]
        cs.leafc_xfer_patch         = [ 0.5,  0.6]
        cs.frootc_xfer_patch        = [ 0.2,  0.25]
        cs.livestemc_xfer_patch     = [ 1.0,  1.1]
        cs.deadstemc_xfer_patch     = [ 1.5,  1.6]
        cs.livecrootc_xfer_patch    = [ 0.5,  0.55]
        cs.deadcrootc_xfer_patch    = [ 0.8,  0.9]
        cs.gresp_xfer_patch         = [ 0.1,  0.15]

        # --- Nitrogen state pools ---
        ns = CLM.CNVegNitrogenStateData()
        ns.leafn_patch              = [ 1.0,  1.2]
        ns.frootn_patch             = [ 0.4,  0.5]
        ns.livestemn_patch          = [ 0.5,  0.6]
        ns.deadstemn_patch          = [ 1.2,  1.4]
        ns.livecrootn_patch         = [ 0.3,  0.35]
        ns.deadcrootn_patch         = [ 0.6,  0.7]
        ns.retransn_patch           = [ 0.2,  0.25]
        ns.leafn_storage_patch      = [ 0.1,  0.12]
        ns.frootn_storage_patch     = [ 0.05, 0.06]
        ns.livestemn_storage_patch  = [ 0.2,  0.22]
        ns.deadstemn_storage_patch  = [ 0.3,  0.33]
        ns.livecrootn_storage_patch = [ 0.1,  0.11]
        ns.deadcrootn_storage_patch = [ 0.15, 0.16]
        ns.leafn_xfer_patch         = [ 0.05, 0.06]
        ns.frootn_xfer_patch        = [ 0.02, 0.025]
        ns.livestemn_xfer_patch     = [ 0.1,  0.11]
        ns.deadstemn_xfer_patch     = [ 0.15, 0.16]
        ns.livecrootn_xfer_patch    = [ 0.05, 0.055]
        ns.deadcrootn_xfer_patch    = [ 0.08, 0.09]

        # --- Flux outputs (zero-initialized) ---
        cf = CLM.CNVegCarbonFluxData()
        for fld in (:gru_leafc_to_litter_patch, :gru_frootc_to_litter_patch,
                    :gru_livestemc_to_atm_patch, :gru_deadstemc_to_atm_patch,
                    :gru_wood_productc_gain_patch, :gru_livecrootc_to_litter_patch,
                    :gru_deadcrootc_to_litter_patch, :gru_xsmrpool_to_atm_patch,
                    :gru_leafc_storage_to_atm_patch, :gru_frootc_storage_to_atm_patch,
                    :gru_livestemc_storage_to_atm_patch, :gru_deadstemc_storage_to_atm_patch,
                    :gru_livecrootc_storage_to_atm_patch, :gru_deadcrootc_storage_to_atm_patch,
                    :gru_gresp_storage_to_atm_patch, :gru_leafc_xfer_to_atm_patch,
                    :gru_frootc_xfer_to_atm_patch, :gru_livestemc_xfer_to_atm_patch,
                    :gru_deadstemc_xfer_to_atm_patch, :gru_livecrootc_xfer_to_atm_patch,
                    :gru_deadcrootc_xfer_to_atm_patch, :gru_gresp_xfer_to_atm_patch)
            setfield!(cf, fld, zeros(np))
        end
        cf.gru_c_to_litr_c_col       = zeros(nc, nlevdecomp, n_litr)
        cf.gru_c_to_cwdc_col         = zeros(nc, nlevdecomp)
        cf.gru_wood_productc_gain_col = zeros(nc)

        nf = CLM.CNVegNitrogenFluxData()
        for fld in (:gru_leafn_to_litter_patch, :gru_frootn_to_litter_patch,
                    :gru_livestemn_to_atm_patch, :gru_deadstemn_to_atm_patch,
                    :gru_wood_productn_gain_patch, :gru_livecrootn_to_litter_patch,
                    :gru_deadcrootn_to_litter_patch, :gru_retransn_to_litter_patch,
                    :gru_leafn_storage_to_atm_patch, :gru_frootn_storage_to_atm_patch,
                    :gru_livestemn_storage_to_atm_patch, :gru_deadstemn_storage_to_atm_patch,
                    :gru_livecrootn_storage_to_atm_patch, :gru_deadcrootn_storage_to_atm_patch,
                    :gru_leafn_xfer_to_atm_patch, :gru_frootn_xfer_to_atm_patch,
                    :gru_livestemn_xfer_to_atm_patch, :gru_deadstemn_xfer_to_atm_patch,
                    :gru_livecrootn_xfer_to_atm_patch, :gru_deadcrootn_xfer_to_atm_patch)
            setfield!(nf, fld, zeros(np))
        end
        nf.gru_n_to_litr_n_col        = zeros(nc, nlevdecomp, n_litr)
        nf.gru_n_to_cwdn_col          = zeros(nc, nlevdecomp)
        nf.gru_wood_productn_gain_col = zeros(nc)

        # --- soilbiogeochem profiles (1/m); uniform = 1.0 so column sums are clean ---
        sbs = CLM.SoilBiogeochemStateData()
        sbs.leaf_prof_patch  = fill(1.0, np, nlevdecomp)
        sbs.froot_prof_patch = fill(1.0, np, nlevdecomp)
        sbs.croot_prof_patch = fill(1.0, np, nlevdecomp)
        sbs.stem_prof_patch  = fill(1.0, np, nlevdecomp)

        mask = trues(np)

        CLM.cn_gross_unrep!(mask, state, pftcon, patch, cs, ns, cf, nf, sbs;
                            dt = dt, noveg = noveg, nc4_grass = nc4_grass,
                            is_beg_curr_year = true,
                            nlevdecomp = nlevdecomp,
                            i_litr_min = i_litr_min, i_litr_max = i_litr_max,
                            i_met_lit = i_met_lit)

        # per-patch rate (1/s): am applied at start of year
        m = [am1 / dt, am2 / dt]

        # ----- (2a) Per-patch CARBON conservation -----
        # For each pool, the sum of its gru fluxes must equal pool * m.
        for p in 1:np
            @test cf.gru_leafc_to_litter_patch[p]      ≈ cs.leafc_patch[p]      * m[p]
            @test cf.gru_frootc_to_litter_patch[p]     ≈ cs.frootc_patch[p]     * m[p]
            @test cf.gru_livestemc_to_atm_patch[p]     ≈ cs.livestemc_patch[p]  * m[p]
            @test cf.gru_livecrootc_to_litter_patch[p] ≈ cs.livecrootc_patch[p] * m[p]
            @test cf.gru_deadcrootc_to_litter_patch[p] ≈ cs.deadcrootc_patch[p] * m[p]
            @test cf.gru_xsmrpool_to_atm_patch[p]      ≈ cs.xsmrpool_patch[p]   * m[p]
            # deadstem splits between atm (convfrac) and wood product (1-convfrac);
            # together they must reconstruct the full pool loss.
            @test cf.gru_deadstemc_to_atm_patch[p] + cf.gru_wood_productc_gain_patch[p] ≈
                  cs.deadstemc_patch[p] * m[p]
        end

        # ----- (2b) Per-patch NITROGEN conservation -----
        for p in 1:np
            @test nf.gru_leafn_to_litter_patch[p]      ≈ ns.leafn_patch[p]      * m[p]
            @test nf.gru_frootn_to_litter_patch[p]     ≈ ns.frootn_patch[p]     * m[p]
            @test nf.gru_livestemn_to_atm_patch[p]     ≈ ns.livestemn_patch[p]  * m[p]
            @test nf.gru_livecrootn_to_litter_patch[p] ≈ ns.livecrootn_patch[p] * m[p]
            @test nf.gru_deadcrootn_to_litter_patch[p] ≈ ns.deadcrootn_patch[p] * m[p]
            @test nf.gru_retransn_to_litter_patch[p]   ≈ ns.retransn_patch[p]   * m[p]
            @test nf.gru_deadstemn_to_atm_patch[p] + nf.gru_wood_productn_gain_patch[p] ≈
                  ns.deadstemn_patch[p] * m[p]
        end

        # ----- (2c) Patch -> column aggregation: leaf+froot C to litter -----
        # With profiles = 1.0 and equal litter fractions 1/n_litr, the total over
        # the litter pools for each (c,j) must equal sum_p (leafc + frootc) flux * wt.
        c, j = 1, 1
        expected_leaffroot_c = sum(
            (cf.gru_leafc_to_litter_patch[p] + cf.gru_frootc_to_litter_patch[p]) *
            patch.wtcol[p] for p in 1:np)
        @test sum(cf.gru_c_to_litr_c_col[c, j, :]) ≈ expected_leaffroot_c

        expected_leaffroot_n = sum(
            (nf.gru_leafn_to_litter_patch[p] + nf.gru_frootn_to_litter_patch[p]) *
            patch.wtcol[p] for p in 1:np)
        # N litter also includes the retransn -> i_met_lit contribution.
        expected_retransn = sum(nf.gru_retransn_to_litter_patch[p] * patch.wtcol[p]
                                for p in 1:np)
        @test sum(nf.gru_n_to_litr_n_col[c, j, :]) ≈ expected_leaffroot_n + expected_retransn

        # ----- (2d) Coarse root C/N to CWD -----
        expected_cwdc = sum(
            (cf.gru_livecrootc_to_litter_patch[p] + cf.gru_deadcrootc_to_litter_patch[p]) *
            patch.wtcol[p] for p in 1:np)
        @test cf.gru_c_to_cwdc_col[c, j] ≈ expected_cwdc

        expected_cwdn = sum(
            (nf.gru_livecrootn_to_litter_patch[p] + nf.gru_deadcrootn_to_litter_patch[p]) *
            patch.wtcol[p] for p in 1:np)
        @test nf.gru_n_to_cwdn_col[c, j] ≈ expected_cwdn

        # ----- (2e) Wood product gain aggregated to column -----
        expected_prodc = sum(cf.gru_wood_productc_gain_patch[p] * patch.wtcol[p]
                             for p in 1:np)
        @test cf.gru_wood_productc_gain_col[c] ≈ expected_prodc
        expected_prodn = sum(nf.gru_wood_productn_gain_patch[p] * patch.wtcol[p]
                             for p in 1:np)
        @test nf.gru_wood_productn_gain_col[c] ≈ expected_prodn

        # ----- (2f) Gating: not start-of-year -> all fluxes zero -----
        for fld in (:gru_leafc_to_litter_patch, :gru_deadstemc_to_atm_patch)
            setfield!(cf, fld, zeros(np))
        end
        cf.gru_c_to_litr_c_col       .= 0.0
        cf.gru_c_to_cwdc_col         .= 0.0
        cf.gru_wood_productc_gain_col .= 0.0
        CLM.cn_gross_unrep!(mask, state, pftcon, patch, cs, ns, cf, nf, sbs;
                            dt = dt, noveg = noveg, nc4_grass = nc4_grass,
                            is_beg_curr_year = false,
                            nlevdecomp = nlevdecomp,
                            i_litr_min = i_litr_min, i_litr_max = i_litr_max,
                            i_met_lit = i_met_lit)
        @test all(cf.gru_leafc_to_litter_patch .== 0.0)
        @test all(cf.gru_deadstemc_to_atm_patch .== 0.0)
    end
end
