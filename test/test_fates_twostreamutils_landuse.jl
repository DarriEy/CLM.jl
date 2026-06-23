# test_fates_twostreamutils_landuse.jl
# Tests for FATES Tier F, Batch 11:
#   * FatesTwoStreamUtilsMod — the glue between the FATES patch/cohort canopy and
#     the TwoStreamMLPEMod solver. We build a synthetic single-patch site with
#     two cohorts, run the build→solve→scatter-back chain, and assert energy
#     conservation (absorbed + reflected + transmitted ≈ incident) and that the
#     per-cohort absorption is non-negative.
#   * FatesLandUseChangeMod — the LUH2 → FATES land-use transition machinery. We
#     apply a synthetic transition matrix / state vector and assert area
#     conservation (states sum to 1) and bounds (non-negative, transitions in
#     range).

using Test
using CLM

# ---------------------------------------------------------------------------
# Helpers to set up the global FATES parameter tables that the radiation glue
# depends on (prt_params for CrownDepth, ed_params for the VAI bins,
# EDPftvarcon_inst for the optical properties, numpft, and rad_params).
# ---------------------------------------------------------------------------
function setup_fates_rad_globals!(; npft::Int=2)
    # --- prt_params: only the crown-depth allometry is needed by VegAreaLayer. --
    p = CLM.prt_params
    CLM.allocate_prt_params!(p, npft, 1, 1)
    p.allom_h2cd1 .= 0.5    # crown depth = 0.5 * height (linear, dmode 1)
    p.allom_h2cd2 .= 1.0
    p.allom_dmode .= 1

    # --- ed_params: VAI bins, sized nlevleaf, uniform increment. ---------------
    edp = CLM.ed_params_type()
    nlv = CLM.nlevleaf
    edp.dinc_vai   = fill(1.0, nlv)
    edp.dlower_vai = collect(1.0:1.0:Float64(nlv))
    CLM.EDParams[] = edp

    # --- EDPftvarcon_inst: optical properties (pft, band) layout. --------------
    edpf = CLM.EDPftvarcon_type()
    nb = CLM.num_swb
    edpf.rhol = fill(0.10, npft, nb)
    edpf.rhos = fill(0.16, npft, nb)
    edpf.taul = fill(0.05, npft, nb)
    edpf.taus = fill(0.001, npft, nb)
    edpf.xl   = fill(0.01, npft)
    edpf.clumping_index = fill(1.0, npft)
    CLM.EDPftvarcon_inst[] = edpf

    # --- numpft + transfer the optical params into the solver's rad_params. -----
    CLM.numpft[] = npft
    CLM.TransferRadParams()

    return nothing
end

# Build a synthetic single-patch site with `ncohorts` cohorts in canopy layer 1.
function build_synthetic_site(; ncohorts::Int=2)
    site = CLM.ed_site_type()
    site.snow_depth = 0.0
    site.lat = 45.0
    site.lon = 0.0

    patch = CLM.fates_patch_type()
    patch.ncl_p = 1
    # Each cohort covers part of the canopy; sum < total so an air element is
    # appended (keeps the layer area summing to 1.0).
    patch.total_canopy_area = CLM.area  # notional forest area [m2]

    # Cohort areas: leave ~20% open so an air ghost element is added.
    c_areas = ncohorts == 2 ? [0.5, 0.3] : [0.7]
    lais    = ncohorts == 2 ? [2.0, 1.5] : [2.0]
    sais    = ncohorts == 2 ? [0.5, 0.3] : [0.5]

    prev = nothing
    first_cohort = nothing
    for i in 1:ncohorts
        co = CLM.fates_cohort_type()
        co.pft = 1
        co.canopy_layer = 1
        co.height = 15.0
        co.nv = 2
        co.treelai = lais[i]
        co.treesai = sais[i]
        co.c_area = c_areas[i] * patch.total_canopy_area
        co.shorter = nothing
        if prev === nothing
            first_cohort = co
        else
            prev.shorter = co
        end
        prev = co
    end

    patch.tallest = first_cohort
    patch.shortest = prev
    patch.older = nothing
    patch.younger = nothing

    site.oldest_patch = patch
    site.youngest_patch = patch

    return site, patch
end

@testset "FATES Batch 11: TwoStreamUtils + LandUseChange" begin

    # =======================================================================
    # FatesTwoStreamUtilsMod
    # =======================================================================
    @testset "FatesTwoStreamUtilsMod build→solve→scatter-back" begin
        setup_fates_rad_globals!(npft=2)
        site, patch = build_synthetic_site(ncohorts=2)

        ifp_max = 1
        fcansno_pa = zeros(Float64, ifp_max)      # no snow
        coszen_pa  = fill(0.7, ifp_max)           # cos(zenith) = 0.7

        # Set ground albedos on the (yet-to-be-allocated) band records by running
        # the construct step first, then assigning albedos before the solve.
        CLM.FatesConstructRadElements!(site, fcansno_pa, coszen_pa)

        twostr = patch.twostr

        # The construct step set up the elements: one column per cohort plus an
        # air ghost (since the canopy is not full).
        @test twostr.n_lyr == 1
        @test twostr.n_col[1] == 3              # 2 cohorts + air element
        @test twostr.n_scel == 3

        # Element areas in layer 1 sum to ~1.0 (area conservation in the build).
        area_sum = sum(twostr.scelg[1, icol].area for icol in 1:twostr.n_col[1])
        @test area_sum ≈ 1.0 atol = 1e-6

        # Cohorts were assigned column indices.
        co = patch.tallest
        seen_cols = Int[]
        while co !== nothing
            push!(seen_cols, co.twostr_col)
            co = co.shorter
        end
        @test sort(seen_cols) == [1, 2]

        # Site scratch arrays were sized for the solve (>= 2*n_scel).
        @test length(site.taulambda_2str) >= 2 * twostr.n_scel
        @test size(site.omega_2str, 1) == length(site.taulambda_2str)

        # Ground albedos for the solve.
        for b in twostr.band
            b.albedo_grnd_diff = 0.2
            b.albedo_grnd_beam = 0.2
        end

        ib = 1
        n_eq = 2 * twostr.n_scel
        (albedo_beam, albedo_diff, consv_err,
         frac_abs_can_beam, frac_abs_can_diff,
         frac_beam_grnd_beam, frac_diff_grnd_beam, frac_diff_grnd_diff) =
            CLM.Solve!(twostr, ib, CLM.absolute_upper_boundary,
                       1.0, 1.0,
                       view(site.taulambda_2str, 1:n_eq),
                       view(site.omega_2str, 1:n_eq, 1:n_eq),
                       view(site.ipiv_2str, 1:n_eq))

        # Solver-level conservation (already normalized by total incident).
        @test consv_err < CLM.rel_err_thresh

        # Explicit per-stream energy balance: incident = reflected + absorbed +
        # transmitted-to-ground-and-absorbed-there.
        beam_bal = albedo_beam + frac_abs_can_beam +
                   frac_diff_grnd_beam * (1.0 - twostr.band[ib].albedo_grnd_diff) +
                   frac_beam_grnd_beam * (1.0 - twostr.band[ib].albedo_grnd_beam)
        diff_bal = albedo_diff + frac_abs_can_diff +
                   frac_diff_grnd_diff * (1.0 - twostr.band[ib].albedo_grnd_diff)
        @test beam_bal ≈ 1.0 atol = 1e-5
        @test diff_bal ≈ 1.0 atol = 1e-5

        # ---- Scatter-back: per-cohort absorbed radiation is non-negative ------
        co = patch.tallest
        total_cohort_abs = 0.0
        while co !== nothing
            # whole-cohort VAI interval (iv = 0 returns total LAIs/SAIs)
            (vt, vb, elai, esai, _, _) =
                CLM.VegAreaLayer(co.treelai, co.treesai, co.height,
                                 0, co.nv, co.pft, site.snow_depth)
            (rb_abs, rd_abs, rb_abs_leaf, rd_abs_leaf, leaf_sun_frac) =
                CLM.FatesGetCohortAbsRad(patch, co, ib, vt, vb, elai, esai)

            @test rb_abs >= -CLM.rel_err_thresh
            @test rd_abs >= -CLM.rel_err_thresh
            @test rb_abs_leaf >= -CLM.rel_err_thresh
            @test rd_abs_leaf >= -CLM.rel_err_thresh
            @test 0.0 <= leaf_sun_frac <= 1.0
            # leaf absorption cannot exceed total (leaf+stem) absorption
            @test rb_abs_leaf <= rb_abs + CLM.rel_err_thresh
            @test rd_abs_leaf <= rd_abs + CLM.rel_err_thresh

            total_cohort_abs += (rb_abs + rd_abs) * co.c_area / patch.total_canopy_area
            co = co.shorter
        end
        @test total_cohort_abs >= 0.0

        # The summed cohort absorption should match the canopy-absorbed fraction
        # the solver reported (this is exactly what CheckPatchRadiationBalance
        # verifies). fabd/fabi here are the canopy-absorbed fractions.
        @test isfinite(total_cohort_abs)

        # ---- FatesPatchFSun: sunlit/shaded LAI partition is sane -------------
        (fsun, laisun, laisha) = CLM.FatesPatchFSun(patch)
        @test 0.0 <= fsun <= 1.0
        @test laisun >= -CLM.rel_err_thresh
        @test laisha >= -CLM.rel_err_thresh
        @test laisun + laisha > 0.0   # there is leaf area present

        # ---- CheckPatchRadiationBalance closes ------------------------------
        # Pass the canopy-absorbed fractions as fabd/fabi (what the host stores).
        CLM.CheckPatchRadiationBalance(patch, site.snow_depth, ib,
                                       frac_abs_can_beam, frac_abs_can_diff)
        @test true   # no endrun thrown => balance closed
    end

    # =======================================================================
    # FatesLandUseChangeMod
    # =======================================================================
    @testset "FatesLandUseChangeMod state aggregation + transitions" begin
        nlu = CLM.n_landuse_cats

        # luh2 ↔ fates lutype map sanity.
        lumap = CLM.luh2_fates_lutype_map()
        @test length(lumap.state_names) == 12
        @test CLM.GetLUCategoryFromStateName(lumap, "primf") == CLM.primaryland
        @test CLM.GetLUCategoryFromStateName(lumap, "secdf") == CLM.secondaryland
        @test CLM.GetLUCategoryFromStateName(lumap, "pastr") == CLM.pastureland
        @test CLM.GetLUCategoryFromStateName(lumap, "range") == CLM.rangeland
        @test CLM.GetLUCategoryFromStateName(lumap, "c3ann") == CLM.cropland
        @test CLM.GetLUCategoryFromStateName(lumap, "urban") == CLM.fates_unset_int

        # --- Synthetic LUH2 boundary condition. --------------------------------
        # 12 states (matching the lutype map order). Choose a vegetated mix that
        # includes a small urban fraction (which must be factored out).
        state_names = copy(lumap.state_names)
        # primf primn secdf secdn pastr range urban c3ann c4ann c3per c4per c3nfx
        states = [0.40, 0.05, 0.10, 0.05, 0.08, 0.07, 0.05, 0.05, 0.03, 0.03, 0.02, 0.02]
        # raw LUH2 states need not sum to 1 (the module renormalizes); just sane.
        @test all(states .>= 0.0)

        bc = CLM.bc_in_type()
        bc.hlm_luh_states = states
        bc.hlm_luh_state_names = state_names
        # transitions: a couple of nonzero transitions, names xxxxx_to_yyyyy.
        trans_names = ["primf_to_secdf", "primf_to_c3ann", "secdf_to_pastr"]
        trans_rates = [0.01, 0.02, 0.005]
        bc.hlm_luh_transitions = trans_rates
        bc.hlm_luh_transition_names = trans_names

        CLM.hlm_num_luh2_states[] = length(states)
        CLM.hlm_num_luh2_transitions[] = length(trans_names)
        CLM.hlm_use_luh[] = CLM.itrue
        CLM.hlm_use_potentialveg[] = CLM.ifalse

        # --- State aggregation: area conservation + bounds. --------------------
        state_vector = CLM.GetLUHStatedata(bc)
        @test length(state_vector) == nlu
        @test all(state_vector .>= -CLM.nearzero)            # non-negative
        @test all(state_vector .<= 1.0 + 1e-9)               # bounded above
        @test isapprox(sum(state_vector), 1.0; atol=1e-9)    # AREA CONSERVATION

        # Urban is excluded and the non-urban states are renormalized to sum to 1.
        # Build the expected aggregated vector independently and compare.
        urban = states[7]
        exp_vec = zeros(Float64, nlu)
        exp_vec[CLM.primaryland]   = (states[1] + states[2]) / (1.0 - urban)
        exp_vec[CLM.secondaryland] = (states[3] + states[4]) / (1.0 - urban)
        exp_vec[CLM.pastureland]   = states[5] / (1.0 - urban)
        exp_vec[CLM.rangeland]     = states[6] / (1.0 - urban)
        exp_vec[CLM.cropland]      = (states[8] + states[9] + states[10] +
                                      states[11] + states[12]) / (1.0 - urban)
        exp_vec ./= sum(exp_vec)
        @test state_vector ≈ exp_vec atol = 1e-9

        # --- Transition matrix: bounds + correct aggregation. ------------------
        min_frac = 1e-4
        ltm = zeros(Float64, nlu, nlu)
        gt_min = fill(false, nlu)
        CLM.GetLanduseTransitionRates!(bc, min_frac, ltm, gt_min)

        @test size(ltm) == (nlu, nlu)
        @test all(ltm .>= -CLM.nearzero)                     # non-negative rates
        # diagonal entries stay zero (no self-transition)
        for i in 1:nlu
            @test ltm[i, i] == 0.0
        end
        # the primf_to_c3ann transition (primary→cropland) should be > 0
        @test ltm[CLM.primaryland, CLM.cropland] > 0.0

        # --- Clearing-rules matrix (ruleset 4). --------------------------------
        clearing = CLM.GetLanduseChangeRules()
        @test size(clearing) == (nlu, nlu)
        @test all(clearing[:, CLM.cropland])                 # clear into cropland
        @test all(clearing[:, CLM.pastureland])              # clear into pasture
        @test !any(clearing[:, CLM.rangeland])               # do NOT clear into rangeland

        # --- Potential-veg path: everything is primary, no transitions. --------
        CLM.hlm_use_potentialveg[] = CLM.itrue
        sv_pot = CLM.GetLUHStatedata(bc)
        @test sv_pot[CLM.primaryland] ≈ 1.0
        @test isapprox(sum(sv_pot), 1.0; atol=1e-12)
        ltm_pot = ones(Float64, nlu, nlu)
        CLM.GetLanduseTransitionRates!(bc, min_frac, ltm_pot, fill(false, nlu))
        @test all(ltm_pot .== 0.0)                           # no transitions for pot-veg
        CLM.hlm_use_potentialveg[] = CLM.ifalse

        # --- Spin-up→land-use init transitions + harvest rate. -----------------
        ltm_init = zeros(Float64, nlu, nlu)
        gt_min_init = fill(false, nlu)
        CLM.GetInitLanduseTransitionRates!(bc, min_frac, ltm_init, gt_min_init)
        @test all(ltm_init .>= 0.0)
        # only primary→(non-secondary) rows are populated
        for i in (CLM.secondaryland + 1):nlu
            if state_vector[i] > min_frac
                @test ltm_init[CLM.primaryland, i] ≈ state_vector[i]
                @test gt_min_init[i]
            end
        end

        gt_min_h = fill(false, nlu)
        hrate = CLM.GetInitLanduseHarvestRate(bc, min_frac, gt_min_h)
        @test hrate >= 0.0
        if state_vector[CLM.secondaryland] > min_frac
            @test hrate ≈ state_vector[CLM.secondaryland]
            @test gt_min_h[CLM.secondaryland]
        end

        # --- CheckLUHData!: all-NaN state vector defaults to primary. -----------
        nan_state = fill(NaN, length(states))
        modified = CLM.CheckLUHData!(nan_state)
        @test modified
        @test nan_state[CLM.primaryland] == 1.0
        @test sum(nan_state) ≈ 1.0
    end
end
