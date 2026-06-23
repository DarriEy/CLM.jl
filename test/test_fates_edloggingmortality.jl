# test_fates_edloggingmortality.jl
# Tests for FATES Tier F, Batch 12: EDLoggingMortalityMod — logging / harvest
# disturbance mortality.
#
# Strategy (mirrors test_fates_soilbgcflux.jl):
#   * Define a lightweight stub PRT object (subtype of CLM.AbstractPRTVartypes)
#     with known per-(organ,element) biomass so we can exercise the mortality-
#     fraction and litter-flux math without standing up the full PRT global
#     descriptor.
#   * Set up the global FATES parameter tables the logging code reads
#     (ed_params logging_*, prt_params woody/allom_agb_frac, sf_params CWD frac,
#     EDPftvarcon_inst decomposability + product fractions, numpft, num_elements).
#   * Verify:
#       - LoggingMortality_frac : direct/collateral/mechanical fractions are in
#         range and consistent (direct+infra+collateral+l_degrad == harvest_rate
#         in the canopy), per harvest scenario and PFT (woody vs grass) / size.
#       - get_harvest_rate_area / get_harvest_rate_carbon : sensible bounded rates
#         and harvest-tag semantics.
#       - logging_litter_fluxes! : MASS CONSERVATION — every kg of killed biomass
#         appears in the CWD/fine-litter destinations + the exported trunk product.

using Test
using CLM

# ---------------------------------------------------------------------------
# A minimal stub PRT object. Logging only needs GetState (organ,element)->mass.
# ---------------------------------------------------------------------------
mutable struct _LogStubPRT <: CLM.AbstractPRTVartypes
    state::Dict{Tuple{Int,Int},Float64}     # (organ, element) -> mass [kg]
end
_LogStubPRT() = _LogStubPRT(Dict{Tuple{Int,Int},Float64}())

function CLM.GetState(this::_LogStubPRT, organ_id::Integer, element_id::Integer,
                      position_id::Union{Nothing,Integer}=nothing)
    return get(this.state, (Int(organ_id), Int(element_id)), 0.0)
end

@testset "FATES Batch 12: EDLoggingMortalityMod" begin

    # --- Preserve and restore mutated module-global state --------------------
    old_numpft   = CLM.numpft[]
    old_numel    = CLM.num_elements[]
    old_ellist   = copy(CLM.element_list)
    old_elpos    = copy(CLM.element_pos)
    old_edparams = CLM.EDParams[]
    old_edpft    = CLM.EDPftvarcon_inst[]
    old_sfparams = CLM.SFParams[]
    old_logtime  = CLM.logging_time[]

    # PRT param table is a live const instance; remember the woody/agb vectors.
    old_woody    = copy(CLM.prt_params.woody)
    old_agb      = copy(CLM.prt_params.allom_agb_frac)
    old_fnrtmode = copy(CLM.prt_params.fnrt_prof_mode)
    old_fnrta    = copy(CLM.prt_params.fnrt_prof_a)
    old_fnrtb    = copy(CLM.prt_params.fnrt_prof_b)

    # hlm flags
    old_use_log    = CLM.hlm_use_logging[]
    old_use_luh    = CLM.hlm_use_luh[]
    old_use_luharv = CLM.hlm_use_lu_harvest[]
    old_nharv      = CLM.hlm_num_lu_harvest_cats[]
    old_dpy        = CLM.hlm_days_per_year[]
    old_planthydro = CLM.hlm_use_planthydro[]
    old_cday       = CLM.hlm_current_day[]

    try
        npft     = 2   # pft 1 = woody tree, pft 2 = grass
        nlevsoil = 3
        ncat     = 5   # five LUH harvest categories (VH1 VH2 SH1 SH2 SH3)

        CLM.numpft[]       = npft
        CLM.num_elements[] = 1                  # carbon-only
        empty!(CLM.element_list); push!(CLM.element_list, CLM.carbon12_element)
        fill!(CLM.element_pos, 0)
        CLM.element_pos[CLM.carbon12_element] = 1

        CLM.hlm_use_logging[]        = CLM.itrue
        CLM.hlm_use_luh[]            = CLM.ifalse
        CLM.hlm_use_lu_harvest[]     = CLM.ifalse  # fates-parameter logging
        CLM.hlm_num_lu_harvest_cats[] = ncat
        CLM.hlm_days_per_year[]      = 365
        CLM.hlm_use_planthydro[]     = CLM.ifalse
        CLM.hlm_current_day[]        = 1

        # --- ed_params: logging fractions + event code (3 = every day). --------
        ep = CLM.ed_params_type()
        ep.logging_direct_frac     = 0.7
        ep.logging_collateral_frac = 0.1
        ep.logging_mechanical_frac = 0.1
        ep.logging_coll_under_frac = 0.5
        ep.logging_dbhmin          = 30.0
        ep.logging_dbhmax          = CLM.fates_check_param_set   # dbhmax disabled
        ep.logging_dbhmax_infra    = 35.0
        ep.logging_export_frac     = 0.8
        ep.logging_event_code      = 2.0    # apply annual rate once (no per-day divisor)
        CLM.EDParams[] = ep

        # --- EDPftvarcon: leaf/fineroot decomposability + product fractions. ---
        pft_con = CLM.EDPftvarcon_type()
        pft_con.lf_flab = [0.30, 0.30]
        pft_con.lf_fcel = [0.45, 0.45]
        pft_con.lf_flig = [0.25, 0.25]      # leaf dcmpy fracs sum to 1
        pft_con.fr_flab = [0.30, 0.30]
        pft_con.fr_fcel = [0.45, 0.45]
        pft_con.fr_flig = [0.25, 0.25]      # fineroot dcmpy fracs sum to 1
        pft_con.harvest_pprod10       = [0.6, 0.6]
        pft_con.landusechange_pprod10 = [0.4, 0.4]
        CLM.EDPftvarcon_inst[] = pft_con

        # --- sf_params: CWD partition fractions (twig/sb/lb/trunk) sum to 1. ---
        sfp = CLM.sf_params_type()
        sfp.SF_val_CWD_frac = [0.045, 0.075, 0.21, 0.67]
        CLM.SFParams[] = sfp

        # --- prt_params: woody flag + AGB fraction + fine-root profile. --------
        CLM.prt_params.woody          = [CLM.itrue, CLM.ifalse]    # tree, grass
        CLM.prt_params.allom_agb_frac = [0.6, 0.6]
        # Fine-root profile: 1-parameter exponential (type 2), scale a (used by
        # set_root_fraction for the below-ground litter split).
        CLM.prt_params.fnrt_prof_mode = [2.0, 2.0]
        CLM.prt_params.fnrt_prof_a    = [4.0, 4.0]
        CLM.prt_params.fnrt_prof_b    = [0.0, 0.0]

        # =====================================================================
        # Helper: build a single-site / single-patch / single-cohort system.
        # =====================================================================
        function build_logging_site(; pft::Int, dbh::Float64, canopy_layer::Int,
                                       nplant::Float64,
                                       lmort_direct::Float64=0.0,
                                       lmort_collateral::Float64=0.0,
                                       lmort_infra::Float64=0.0,
                                       sapw=2.0, struct_c=3.0, leaf=0.5,
                                       fnrt=0.4, store=0.3, repro=0.1)
            site = CLM.ed_site_type()
            site.nlevsoil = nlevsoil
            site.zi_soil  = [0.0, 0.1, 0.3, 0.6]    # length nlevsoil+1 (zero index)
            site.rootfrac_scr = zeros(Float64, nlevsoil)
            site.spread = 1.0
            site.transition_landuse_from_off_to_on = false
            site.min_allowed_landuse_fraction = 1e-4
            site.landuse_vector_gt_min = fill(false, CLM.n_landuse_cats)

            # per-element mass balance + flux diagnostics
            mb = CLM.site_massbal_type()
            site.mass_balance = [mb]
            ed = CLM.elem_diag_type()
            ed.surf_fine_litter_input = zeros(Float64, npft)
            ed.root_litter_input      = zeros(Float64, npft)
            site.flux_diags.elem = [ed]

            # PRT stub with known biomass
            prt = _LogStubPRT()
            el = CLM.carbon12_element
            prt.state[(CLM.sapw_organ,   el)] = sapw
            prt.state[(CLM.struct_organ, el)] = struct_c
            prt.state[(CLM.leaf_organ,   el)] = leaf
            prt.state[(CLM.fnrt_organ,   el)] = fnrt
            prt.state[(CLM.store_organ,  el)] = store
            prt.state[(CLM.repro_organ,  el)] = repro

            co = CLM.fates_cohort_type()
            co.pft = pft
            co.dbh = dbh
            co.n = nplant
            co.canopy_layer = canopy_layer
            co.crowndamage = 1
            co.prt = prt
            co.lmort_direct     = lmort_direct
            co.lmort_collateral = lmort_collateral
            co.lmort_infra      = lmort_infra

            patch = CLM.fates_patch_type()
            patch.area = CLM.area
            patch.land_use_label = CLM.primaryland
            patch.age_since_anthro_disturbance = 0.0
            patch.fract_ldist_not_harvested = 0.0
            patch.nocomp_pft_label = pft
            patch.tallest  = co
            patch.shortest = co
            patch.older = nothing
            patch.younger = nothing

            litt = CLM.litter_type()
            CLM.InitAllocate!(litt, npft, nlevsoil, el)
            CLM.ZeroFlux!(litt)
            CLM.InitConditions!(litt, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
            patch.litter = [litt]

            site.oldest_patch = patch
            site.youngest_patch = patch
            return site, patch, co
        end

        # Dummy harvest inputs (unused on the fates-parameter path).
        harv_rates = zeros(Float64, ncat)
        harv_names = ["HARVEST_VH1", "HARVEST_VH2", "HARVEST_SH1",
                      "HARVEST_SH2", "HARVEST_SH3"]
        harvestable_c = zeros(Float64, ncat)

        # =====================================================================
        # 1. IsItLoggingTime — event code 2 on the first model day.
        # =====================================================================
        @testset "IsItLoggingTime" begin
            site, _, _ = build_logging_site(pft=1, dbh=50.0, canopy_layer=1,
                                            nplant=0.1)
            CLM.hlm_model_day[] = 1
            CLM.IsItLoggingTime(1, site)
            @test CLM.logging_time[] == true

            CLM.hlm_model_day[] = 5
            CLM.IsItLoggingTime(1, site)
            @test CLM.logging_time[] == false

            # logging off entirely
            CLM.hlm_use_logging[] = CLM.ifalse
            CLM.IsItLoggingTime(1, site)
            @test CLM.logging_time[] == false
            CLM.hlm_use_logging[] = CLM.itrue
        end

        # =====================================================================
        # 2. LoggingMortality_frac — fraction logic per scenario / PFT / size.
        # =====================================================================
        @testset "LoggingMortality_frac fractions" begin
            CLM.logging_time[] = true   # logging is happening

            # --- (a) Big woody tree in the canopy: harvest_rate = 1 (fates param)
            site, _, _ = build_logging_site(pft=1, dbh=50.0, canopy_layer=1, nplant=0.1)
            (ld, lc, li, lg, htag) = CLM.LoggingMortality_frac(
                site, CLM.bc_in_type(), 1, 50.0, 1,
                harv_rates, harv_names, CLM.hlm_harvest_area_fraction,
                CLM.primaryland, 0.0, 1.0, -1.0, harvestable_c)

            @test ld ≈ 0.7   # harvest_rate(1) * logging_direct_frac(0.7)
            @test lc ≈ 0.1   # canopy collateral
            # dbh 50 >= dbhmax_infra 35 -> lmort_infra = 0
            @test li ≈ 0.0
            # l_degrad = harvest_rate - (direct + infra + collateral)
            @test lg ≈ 1.0 - (ld + li + lc)
            # all fractions in [0,1]
            for f in (ld, lc, li, lg)
                @test 0.0 <= f <= 1.0 + 1e-12
            end
            # killed + degraded + retained == harvest_rate (mass/area closure)
            @test (ld + lc + li + lg) ≈ 1.0

            # --- (b) Small woody tree below dbhmin: no direct logging.
            site2, _, _ = build_logging_site(pft=1, dbh=10.0, canopy_layer=1, nplant=0.1)
            (ld2, lc2, li2, lg2, _) = CLM.LoggingMortality_frac(
                site2, CLM.bc_in_type(), 1, 10.0, 1,
                harv_rates, harv_names, CLM.hlm_harvest_area_fraction,
                CLM.primaryland, 0.0, 1.0, -1.0, harvestable_c)
            @test ld2 ≈ 0.0           # below dbhmin -> no direct
            @test li2 ≈ 0.1           # dbh 10 < dbhmax_infra 35 -> mechanical applies
            @test lc2 ≈ 0.1           # canopy collateral
            @test lg2 ≈ 1.0 - (ld2 + li2 + lc2)
            @test (ld2 + lc2 + li2 + lg2) ≈ 1.0

            # --- (c) Grass (non-woody): only mechanical mortality.
            site3, _, _ = build_logging_site(pft=2, dbh=2.0, canopy_layer=1, nplant=0.1)
            (ld3, lc3, li3, lg3, _) = CLM.LoggingMortality_frac(
                site3, CLM.bc_in_type(), 2, 2.0, 1,
                harv_rates, harv_names, CLM.hlm_harvest_area_fraction,
                CLM.primaryland, 0.0, 1.0, -1.0, harvestable_c)
            @test ld3 ≈ 0.0
            @test lc3 ≈ 0.0
            @test li3 ≈ 0.1
            @test lg3 ≈ 1.0 - li3

            # --- (d) Understory cohort (canopy_layer=2): no collateral/degrad.
            site4, _, _ = build_logging_site(pft=1, dbh=50.0, canopy_layer=2, nplant=0.1)
            (ld4, lc4, li4, lg4, _) = CLM.LoggingMortality_frac(
                site4, CLM.bc_in_type(), 1, 50.0, 2,
                harv_rates, harv_names, CLM.hlm_harvest_area_fraction,
                CLM.primaryland, 0.0, 1.0, -1.0, harvestable_c)
            @test lc4 ≈ 0.0   # collateral only applies in canopy layer 1
            @test lg4 ≈ 0.0   # l_degrad only in canopy layer 1
        end

        # =====================================================================
        # 3. get_harvest_rate_area — bounded rate + category aggregation.
        # =====================================================================
        @testset "get_harvest_rate_area" begin
            rates = [0.3, 0.1, 0.0, 0.0, 0.0]   # VH1=0.3, VH2=0.1 (primary)
            # event code 2 -> apply annual rate once (no divisor)
            hr = CLM.get_harvest_rate_area(CLM.primaryland, harv_names, rates,
                                           1.0, -1.0, -1.0, 0.0)
            @test hr ≈ 0.4                       # VH1 + VH2, frac_site_primary=1
            @test 0.0 <= hr <= 1.0

            # capping at 1: huge rate / small primary fraction
            hr_cap = CLM.get_harvest_rate_area(CLM.primaryland, harv_names,
                                               [0.9, 0.9, 0.0, 0.0, 0.0],
                                               0.5, -1.0, -1.0, 0.0)
            @test hr_cap ≈ 1.0                   # min(1.8/0.5, 1) = 1

            # no primary area -> zero
            hr_zero = CLM.get_harvest_rate_area(CLM.primaryland, harv_names, rates,
                                                0.0, -1.0, -1.0, 0.0)
            @test hr_zero ≈ 0.0
        end

        # =====================================================================
        # 4. get_harvest_rate_carbon — tag semantics + carbon->area conversion.
        # =====================================================================
        @testset "get_harvest_rate_carbon" begin
            # demand for primary (VH1 + VH2) = 0.2 + 0.0; supply (harvestable_c) ample
            rates = [2.0, 0.0, 0.0, 0.0, 0.0]   # 2 kgC demand on VH1
            supply = [10.0, 0.0, 0.0, 0.0, 0.0] # 10 kgC available on VH1
            htag = fill(2, ncat)
            (hr, cur) = CLM.get_harvest_rate_carbon(CLM.primaryland, harv_names,
                                                    rates, 0.0, supply, htag)
            @test htag[1] == 0                   # enough carbon -> success
            @test cur == 0
            @test hr ≈ 2.0 / 10.0                # demand / supply
            @test 0.0 <= hr <= 1.0

            # not enough carbon -> tag 1, rate 0
            htag2 = fill(2, ncat)
            supply2 = [1.0, 0.0, 0.0, 0.0, 0.0]  # 1 kgC < 2 demand
            (hr2, cur2) = CLM.get_harvest_rate_carbon(CLM.primaryland, harv_names,
                                                      rates, 0.0, supply2, htag2)
            @test htag2[1] == 1                  # insufficient carbon
            @test cur2 == 1
            @test hr2 ≈ 0.0                      # no supply -> rate 0
        end

        # =====================================================================
        # 5. logging_litter_fluxes! — MASS CONSERVATION of killed biomass.
        # =====================================================================
        @testset "logging_litter_fluxes! mass conservation" begin
            CLM.logging_time[] = true

            # Canopy tree with prescribed mortality fractions. Use a new-patch
            # that takes a fraction of the area (disturbance generation).
            sapw = 2.0; struct_c = 3.0; leaf = 0.5; fnrt = 0.4; store = 0.3; repro = 0.1
            nplant = 0.1
            ld = 0.7; lc = 0.1; li = 0.1

            site, cpatch, co = build_logging_site(pft=1, dbh=50.0, canopy_layer=1,
                nplant=nplant, lmort_direct=ld, lmort_collateral=lc, lmort_infra=li,
                sapw=sapw, struct_c=struct_c, leaf=leaf, fnrt=fnrt, store=store,
                repro=repro)

            # Build a new (disturbed) patch with its own litter pool + cohort.
            patch_site_areadis = 0.25 * CLM.area
            npatch = CLM.fates_patch_type()
            npatch.area = patch_site_areadis
            npatch.land_use_label = CLM.secondaryland
            npatch.age_since_anthro_disturbance = 0.0
            npatch.fract_ldist_not_harvested = 0.0
            nlitt = CLM.litter_type()
            CLM.InitAllocate!(nlitt, npft, nlevsoil, CLM.carbon12_element)
            CLM.ZeroFlux!(nlitt)
            CLM.InitConditions!(nlitt, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
            npatch.litter = [nlitt]
            npatch.tallest = nothing       # no cohorts moved in (skip carea recompute)
            npatch.shortest = nothing

            bc = CLM.bc_in_type()
            bc.max_rooting_depth_index_col = nlevsoil
            CLM.logging_litter_fluxes!(site, cpatch, npatch, patch_site_areadis, bc)

            # --- expected total killed biomass (number of individuals) ----------
            direct_dead   = nplant * ld
            indirect_dead = nplant * (lc + li)
            total_dead    = direct_dead + indirect_dead

            # Total killed biomass per organ (kgC/site).
            killed_wood  = total_dead * (struct_c + sapw)
            killed_leaf  = total_dead * (leaf + repro)
            killed_root  = total_dead * (fnrt + store)
            killed_total = killed_wood + killed_leaf + killed_root

            # --- where did the killed mass go? ---------------------------------
            # CWD/litter pools are area-normalized [kg/m2]; multiply each patch's
            # pool by that patch's area to recover [kg/site], then sum donor+new.
            cl = cpatch.litter[1]
            nl = npatch.litter[1]
            remainder_area = cpatch.area - patch_site_areadis

            litter_in_cur = (sum(cl.ag_cwd) + sum(cl.bg_cwd) +
                             sum(cl.leaf_fines) + sum(cl.root_fines)) * remainder_area
            litter_in_new = (sum(nl.ag_cwd) + sum(nl.bg_cwd) +
                             sum(nl.leaf_fines) + sum(nl.root_fines)) * npatch.area

            # exported trunk product [kg/site]
            exported = site.mass_balance[1].wood_product_harvest[1]

            recovered = litter_in_cur + litter_in_new + exported

            # MASS CONSERVATION: every kg of killed biomass is accounted for in
            # litter destinations + the exported product (within rounding).
            @test isapprox(recovered, killed_total; rtol=1e-10)

            # The exported product is a fraction of the directly-logged trunk only,
            # so it must be strictly positive and below the total killed mass.
            @test exported > 0.0
            @test exported < killed_total

            # Resources-management diagnostics are populated and non-negative.
            @test site.resources_management.delta_biomass_stock > 0.0
            @test site.resources_management.delta_individual ≈ total_dead
            @test site.resources_management.trunk_product_site ≈ exported

            # Diagnostics: exported_harvest [kg/m2] = exported * area_inv.
            @test site.flux_diags.elem[1].exported_harvest ≈ exported * CLM.area_inv
        end

        # =====================================================================
        # 6. UpdateHarvestC! — product-pool split (10yr vs 100yr).
        # =====================================================================
        @testset "UpdateHarvestC!" begin
            site, _, _ = build_logging_site(pft=1, dbh=50.0, canopy_layer=1, nplant=0.1)
            # Seed the wood-product pools [kg/site/day].
            site.mass_balance[1].wood_product_harvest[1] = 100.0
            site.mass_balance[1].wood_product_landusechange[1] = 50.0

            bc_out = CLM.bc_out_type()
            CLM.UpdateHarvestC!(site, bc_out)

            # 10yr + 100yr must sum to the total flux (converted kgC/day -> gC/s).
            uf = CLM.g_per_kg * CLM.days_per_sec * CLM.area_inv
            tot10  = bc_out.hrv_deadstemc_to_prod10c
            tot100 = bc_out.hrv_deadstemc_to_prod100c
            expected_total = (100.0 + 50.0) * uf
            @test isapprox(tot10 + tot100, expected_total; rtol=1e-10)
            @test tot10 > 0.0
            @test tot100 > 0.0
        end

        # =====================================================================
        # 7. get_harvest_debt! — accumulate debt where tag == 1.
        # =====================================================================
        @testset "get_harvest_debt!" begin
            CLM.logging_time[] = true
            site, _, _ = build_logging_site(pft=1, dbh=50.0, canopy_layer=1, nplant=0.1)
            site.resources_management.harvest_debt = 0.0
            site.resources_management.harvest_debt_sec = 0.0

            bc = CLM.bc_in_type()
            bc.hlm_harvest_rates = [0.3, 0.1, 0.0, 0.0, 0.0]   # primary demand 0.4
            bc.hlm_harvest_catnames = harv_names

            htag = [1, 2, 2, 2, 2]   # primary (VH1) unsuccessful
            CLM.get_harvest_debt!(site, bc, htag)
            @test site.resources_management.harvest_debt ≈ 0.4   # VH1 + VH2 demand

            # When not logging time, nothing accumulates.
            CLM.logging_time[] = false
            site.resources_management.harvest_debt = 0.0
            CLM.get_harvest_debt!(site, bc, htag)
            @test site.resources_management.harvest_debt ≈ 0.0
        end

    finally
        # --- restore mutated module-global state -----------------------------
        CLM.numpft[]       = old_numpft
        CLM.num_elements[] = old_numel
        empty!(CLM.element_list); append!(CLM.element_list, old_ellist)
        copyto!(CLM.element_pos, old_elpos)
        CLM.EDParams[]         = old_edparams
        CLM.EDPftvarcon_inst[] = old_edpft
        CLM.SFParams[]         = old_sfparams
        CLM.logging_time[]     = old_logtime
        CLM.prt_params.woody          = old_woody
        CLM.prt_params.allom_agb_frac = old_agb
        CLM.prt_params.fnrt_prof_mode = old_fnrtmode
        CLM.prt_params.fnrt_prof_a    = old_fnrta
        CLM.prt_params.fnrt_prof_b    = old_fnrtb
        CLM.hlm_use_logging[]         = old_use_log
        CLM.hlm_use_luh[]             = old_use_luh
        CLM.hlm_use_lu_harvest[]      = old_use_luharv
        CLM.hlm_num_lu_harvest_cats[] = old_nharv
        CLM.hlm_days_per_year[]       = old_dpy
        CLM.hlm_use_planthydro[]      = old_planthydro
        CLM.hlm_current_day[]         = old_cday
    end
end
