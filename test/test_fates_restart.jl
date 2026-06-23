# test_fates_restart.jl
# Tests for FATES Batch 18 (Tier F): the RESTART I/O pair —
# FatesRestartVariableType (one restart variable: metadata + int/real payload)
# and FatesRestartInterfaceMod (the registry of restart variables + the
# site -> patch -> cohort linked-list pack/unpack serialization).
#
# Coverage:
#   * fates_restart_variable_type Init!/Flush! for BOTH int and real payloads.
#   * The interface type Init! + assemble_restart_output_types! + the registry
#     build (non-empty; a sample of registered names match expected dim-kind).
#   * The load-bearing ROUND-TRIP: build a site with 2 patches + cohorts, PACK
#     via set_restart_vectors!, build a FRESH skeleton via
#     create_patchcohort_structure! (from the stored counts), UNPACK via
#     get_restart_vectors!, and assert the demographic state (patch area/age/
#     land-use/nocomp + cohort n/dbh/pft/height/coage/canopy_layer/crowndamage/
#     status) is reconstructed identically, including the count consistency
#     checks.

using Test
using CLM

# --- single woody evergreen carbon PFT (mirrors the cohort/patch suites) -----
function _setup_restart_pft!()
    npft = 1
    p = CLM.prt_params
    CLM.allocate_prt_params!(p, npft, CLM.num_organ_types, 1)

    p.c2b          .= 2.0
    p.wood_density .= 0.6
    p.slatop       .= 0.012
    p.slamax       .= 0.020
    p.allom_agb_frac .= 0.6
    p.allom_dbh_maxheight .= 90.0
    p.allom_la_per_sa_int .= 0.8e3
    p.allom_la_per_sa_slp .= 0.0
    p.allom_sai_scaler    .= 0.1
    p.allom_l2fr          .= 1.0
    p.cushion             .= 1.0
    p.allom_blca_expnt_diff      .= 0.0
    p.allom_d2ca_coefficient_min .= 0.3
    p.allom_d2ca_coefficient_max .= 0.6
    p.allom_h2cd1 .= 0.5
    p.allom_h2cd2 .= 1.0
    p.allom_dmode .= 1
    p.woody       .= 1
    p.allom_cmode .= 1
    p.allom_hmode .= 1
    p.allom_lmode .= 1
    p.allom_fmode .= 1
    p.allom_amode .= 1
    p.allom_smode .= 1
    p.allom_stmode .= 1
    p.allom_d2h1 .= 0.64
    p.allom_d2h2 .= 0.37
    p.allom_d2h3 .= -999.9
    p.allom_d2bl1 .= 0.06
    p.allom_d2bl2 .= 1.3
    p.allom_d2bl3 .= 0.55
    p.allom_agb1 .= 0.0673
    p.allom_agb2 .= 0.976
    p.allom_agb3 .= 1.0
    p.allom_agb4 .= 1.0
    p.allom_zroot_max_dbh .= 100.0
    p.allom_zroot_max_z   .= 100.0
    p.allom_zroot_min_dbh .= 1.0
    p.allom_zroot_min_z   .= 100.0
    p.allom_zroot_k       .= 2.0
    p.fnrt_prof_a .= 7.0
    p.fnrt_prof_b .= 1.0
    p.fnrt_prof_mode .= 1

    # param_derived: branch_frac (used by bsap) + canopy-top derived rates.
    pd = CLM.param_derived_type()
    pd.branch_frac = fill(0.25, npft)
    pd.jmax25top = fill(85.0, npft, 1)
    pd.tpu25top  = fill(8.0, npft, 1)
    pd.kp25top   = fill(0.6, npft, 1)
    CLM.ParamDerived[] = pd

    ev = CLM.EDPftvarcon_type()
    ev.vcmax25top = fill(50.0, npft, 1)
    ev.damage_recovery_scalar = fill(0.0, npft)
    ev.lf_flab = fill(0.25, npft); ev.lf_fcel = fill(0.50, npft); ev.lf_flig = fill(0.25, npft)
    ev.fr_flab = fill(0.25, npft); ev.fr_fcel = fill(0.50, npft); ev.fr_flig = fill(0.25, npft)
    CLM.EDPftvarcon_inst[] = ev

    edp = CLM.ed_params_type()
    edp.regeneration_model = CLM.default_regeneration
    edp.ED_val_history_damage_bin_edges = [0.0, 20.0, 50.0, 80.0]
    edp.ED_val_history_sizeclass_bin_edges  = [0.0, 5.0, 10.0, 20.0, 50.0]
    edp.ED_val_history_coageclass_bin_edges = [0.0, 3.0, 6.0, 9.0, 12.0]
    edp.ED_val_history_ageclass_bin_edges   = [0.0, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]
    edp.ED_val_cohort_size_fusion_tol = 0.06
    edp.ED_val_cohort_age_fusion_tol  = 0.08
    edp.max_cohort_per_patch          = 100
    edp.ED_val_patch_fusion_tol       = 0.05
    edp.maxpatches_by_landuse         = fill(10, CLM.n_landuse_cats)
    nlv = CLM.nlevleaf
    edp.dinc_vai   = fill(1.0, nlv)
    edp.dlower_vai = collect(1.0:1.0:Float64(nlv))
    CLM.EDParams[] = edp

    CLM.nlevsclass[] = length(edp.ED_val_history_sizeclass_bin_edges)
    CLM.nlevcoage[]  = length(edp.ED_val_history_coageclass_bin_edges)

    sfp = CLM.sf_params_type()
    sfp.SF_val_CWD_frac = [0.045, 0.075, 0.21, 0.67]
    CLM.SFParams[] = sfp

    return npft
end

function _seed_prt_r(ipft::Int, dbh0::Float64)
    prt = CLM.callom_prt_vartypes()
    CLM.InitPRTVartype!(prt)
    canopy_trim = 1.0; crowndamage = 1; ef = 1.0
    l2fr = CLM.prt_params.allom_l2fr[ipft]
    tgt_leaf, _    = CLM.bleaf(dbh0, ipft, crowndamage, canopy_trim, ef)
    tgt_fnrt, _    = CLM.bfineroot(dbh0, ipft, canopy_trim, l2fr, ef)
    _, tgt_sapw, _ = CLM.bsap_allom(dbh0, ipft, crowndamage, canopy_trim, ef)
    tgt_store, _   = CLM.bstore_allom(dbh0, ipft, crowndamage, canopy_trim)
    tgt_agw, _     = CLM.bagw_allom(dbh0, ipft, crowndamage, ef)
    tgt_bgw, _     = CLM.bbgw_allom(dbh0, ipft, ef)
    tgt_struct, _  = CLM.bdead_allom(tgt_agw, tgt_bgw, tgt_sapw, ipft)
    CLM.SetState!(prt, CLM.leaf_organ,   CLM.carbon12_element, tgt_leaf)
    CLM.SetState!(prt, CLM.fnrt_organ,   CLM.carbon12_element, tgt_fnrt)
    CLM.SetState!(prt, CLM.sapw_organ,   CLM.carbon12_element, tgt_sapw)
    CLM.SetState!(prt, CLM.store_organ,  CLM.carbon12_element, tgt_store)
    CLM.SetState!(prt, CLM.struct_organ, CLM.carbon12_element, tgt_struct)
    CLM.SetState!(prt, CLM.repro_organ,  CLM.carbon12_element, 0.0)
    return prt
end

# Allocate + deterministically fill the R5 site-level demographic / mortality /
# flux-diagnostic arrays so the round-trip exercises them. Values are seeded by
# index so the round-trip can assert exact recovery.
function _alloc_fill_site_diagnostics!(site, npft::Int)
    nsc = CLM.nlevsclass[]
    nmt = CLM.n_term_mort_types
    nel = CLM.num_elements[]
    nlu = CLM.n_landuse_cats

    site.term_nindivs_canopy = [Float64(100*t + 10*c + p) for t in 1:nmt, c in 1:nsc, p in 1:npft]
    site.term_nindivs_ustory = [Float64(200*t + 10*c + p) for t in 1:nmt, c in 1:nsc, p in 1:npft]
    site.demotion_rate  = [Float64(3*c)    for c in 1:nsc]
    site.promotion_rate = [Float64(4*c)    for c in 1:nsc]
    site.imort_rate         = [Float64(5*c + p)  for c in 1:nsc, p in 1:npft]
    site.fmort_rate_canopy  = [Float64(6*c + p)  for c in 1:nsc, p in 1:npft]
    site.fmort_rate_ustory  = [Float64(7*c + p)  for c in 1:nsc, p in 1:npft]
    site.fmort_rate_cambial = [Float64(8*c + p)  for c in 1:nsc, p in 1:npft]
    site.fmort_rate_crown   = [Float64(9*c + p)  for c in 1:nsc, p in 1:npft]
    site.growthflux_fusion  = [Float64(11*c + p) for c in 1:nsc, p in 1:npft]
    site.term_carbonflux_canopy = [Float64(12*t + p) for t in 1:nmt, p in 1:npft]
    site.term_carbonflux_ustory = [Float64(13*t + p) for t in 1:nmt, p in 1:npft]
    site.imort_carbonflux        = [Float64(14 + p) for p in 1:npft]
    site.fmort_carbonflux_canopy = [Float64(15 + p) for p in 1:npft]
    site.fmort_carbonflux_ustory = [Float64(16 + p) for p in 1:npft]
    site.term_abg_flux  = [Float64(17*c + p) for c in 1:nsc, p in 1:npft]
    site.imort_abg_flux = [Float64(18*c + p) for c in 1:nsc, p in 1:npft]
    site.fmort_abg_flux = [Float64(19*c + p) for c in 1:nsc, p in 1:npft]

    site.rec_l2fr = [Float64(20*p + cl) for p in 1:npft, cl in 1:CLM.nclmax]
    site.area_PFT = [Float64(p) + 0.1*lu for p in 1:npft, lu in 1:nlu]
    site.use_this_pft = [p for p in 1:npft]
    site.recruitment_rate = [Float64(21 + p) for p in 1:npft]
    site.seed_in  = [Float64(22 + p) for p in 1:npft]
    site.seed_out = [Float64(23 + p) for p in 1:npft]
    site.landuse_vector_gt_min = [isodd(lu) for lu in 1:nlu]
    site.min_allowed_landuse_fraction = 0.001
    site.area_bareground = 0.05

    site.dstatus       = [10 + p for p in 1:npft]
    site.dleafondate   = [50 + p for p in 1:npft]
    site.dleafoffdate  = [60 + p for p in 1:npft]
    site.dndaysleafon  = [70 + p for p in 1:npft]
    site.dndaysleafoff = [80 + p for p in 1:npft]
    site.elong_factor  = [Float64(p) / 10 for p in 1:npft]

    site.liqvol_memory = [Float64(24 + i + p) for i in 1:CLM.numWaterMem, p in 1:npft]
    site.smp_memory    = [Float64(25 + i + p) for i in 1:CLM.numWaterMem, p in 1:npft]
    site.vegtemp_memory = [Float64(26 + i) for i in 1:CLM.num_vegtemp_mem]

    site.term_crownarea_canopy = 31.0; site.term_crownarea_ustory = 32.0
    site.ema_npp = 33.0
    site.demotion_carbonflux = 34.0; site.promotion_carbonflux = 35.0
    site.imort_crownarea = 36.0
    site.fmort_crownarea_canopy = 37.0; site.fmort_crownarea_ustory = 38.0

    # flux diagnostics + mass balance (one per element)
    site.flux_diags = CLM.site_fluxdiags_type()
    site.flux_diags.elem = [CLM.elem_diag_type() for _ in 1:nel]
    site.mass_balance  = [CLM.site_massbal_type()  for _ in 1:nel]
    site.iflux_balance = [CLM.site_ifluxbal_type() for _ in 1:nel]
    for el in 1:nel
        e = site.flux_diags.elem[el]
        e.cwd_ag_input = [Float64(40 + el + c) for c in 1:CLM.ncwd]
        e.cwd_bg_input = [Float64(41 + el + c) for c in 1:CLM.ncwd]
        e.surf_fine_litter_input = [Float64(42 + el + p) for p in 1:npft]
        e.root_litter_input      = [Float64(43 + el + p) for p in 1:npft]
        e.err_liveveg = Float64(44 + el); e.err_litter = Float64(45 + el)
        m = site.mass_balance[el]
        m.old_stock = Float64(46 + el); m.err_fates = Float64(47 + el)
        m.wood_product_harvest       = [Float64(48 + el + p) for p in 1:npft]
        m.wood_product_landusechange = [Float64(49 + el + p) for p in 1:npft]
        ib = site.iflux_balance[el]
        ib.iflux_liveveg = Float64(50 + el); ib.iflux_litter = Float64(51 + el)
    end
    return nothing
end

function _build_restart_site(npft::Int, nlevsoil::Int)
    site = CLM.ed_site_type()
    site.nlevsoil     = nlevsoil
    site.zi_soil      = collect(0.0:0.1:0.1 * nlevsoil)
    site.dz_soil      = fill(0.1, nlevsoil)
    site.rootfrac_scr = zeros(nlevsoil)
    site.spread       = 1.0
    # site-level phenology state we round-trip
    site.cstatus = 2; site.nchilldays = 5; site.ncolddays = 3
    site.cleafondate = 100; site.cleafoffdate = 280
    site.cndaysleafon = 30; site.cndaysleafoff = 10
    site.phen_model_date = 12345
    site.grow_deg_days = 42.5; site.snow_depth = 0.13
    site.resources_management.trunk_product_site = 0.7
    _alloc_fill_site_diagnostics!(site, npft)
    return site
end

function _make_restart_patch(area::Float64, age::Float64, npft::Int, nlevsoil::Int)
    patch = CLM.fates_patch_type()
    CLM.Create!(patch, age, area, CLM.primaryland, CLM.fates_unset_int, CLM.num_swb,
                npft, nlevsoil, 0, CLM.ed_params().regeneration_model)
    patch.canopy_layer_tlai = zeros(Float64, CLM.nclmax)
    patch.livegrass = 0.0
    patch.fcansno = 0.0
    patch.solar_zenith_flag = false
    patch.solar_zenith_angle = 0.0
    patch.zstar = 0.0
    return patch
end

function _add_restart_cohort!(site, patch, ipft, dbh0, nn; clayer=1)
    prt = _seed_prt_r(ipft, dbh0)
    h, _ = CLM.h_allom(dbh0, ipft)
    CLM.create_cohort(site, patch, ipft, nn, h, 0.0, dbh0, prt, 1.0, 1.0, 1.0,
                      CLM.leaves_on, 0, 1.0, 0.0, clayer, 1, site.spread)
    c = patch.tallest
    while c !== nothing; c.isnew = false; c = c.shorter; end
    return patch.tallest
end

# Build the restart interface (Init + dim setup + registry).
function _build_interface()
    ri = CLM.fates_restart_interface_type()
    fb = CLM.fates_bounds_type()
    # one site, two patches worth of cohort slots.
    fb.cohort_begin = 1
    fb.cohort_end   = 5 * CLM.fates_maxElementsPerPatch[]   # plenty of room
    fb.column_begin = 1
    fb.column_end   = 1
    CLM.Init!(ri, 1, fb)
    CLM.SetThreadBoundsEach!(ri, 1, fb)
    CLM.assemble_restart_output_types!(ri)
    CLM.initialize_restart_vars!(ri)
    # restart_map: site 1 starts at cohort slot 1.
    ri.restart_map[1].site_index    = [1]
    ri.restart_map[1].cohort1_index = [1]
    return ri
end

# walk a patch shortest->taller collecting (pft, dbh, n, height, canopy_layer).
function _cohort_tuples(patch)
    out = Tuple{Int,Float64,Float64,Float64,Int}[]
    c = patch.shortest
    while c !== nothing
        push!(out, (c.pft, c.dbh, c.n, c.height, c.canopy_layer))
        c = c.taller
    end
    return out
end

function _patch_list(site)
    out = CLM.fates_patch_type[]
    p = site.oldest_patch
    while p !== nothing; push!(out, p); p = p.younger; end
    return out
end

@testset "FATES Batch 18: Restart I/O" begin

    # ---- save/restore global FATES config ------------------------------------
    old_global   = CLM.prt_global[]
    old_numel    = CLM.num_elements[]
    old_ellist   = copy(CLM.element_list)
    old_parteh   = CLM.hlm_parteh_mode[]
    old_use_sp   = CLM.hlm_use_sp[]
    old_hydro    = CLM.hlm_use_planthydro[]
    old_nocomp   = CLM.hlm_use_nocomp[]
    old_fixbio   = CLM.hlm_use_fixed_biogeog[]
    old_potveg   = CLM.hlm_use_potentialveg[]
    old_freqday  = CLM.hlm_freq_day[]
    old_tod      = CLM.hlm_current_tod[]
    old_numpft   = CLM.numpft[]
    old_name     = CLM.hlm_name[]
    old_maxpp    = CLM.fates_maxElementsPerPatch[]
    old_maxps    = CLM.fates_maxElementsPerSite[]
    old_nlevsc   = CLM.nlevsclass[]
    old_nlevca   = CLM.nlevcoage[]
    old_sfp      = CLM.SFParams[]
    old_edp      = CLM.EDParams[]
    old_pd       = CLM.ParamDerived[]

    try
        npft = _setup_restart_pft!()
        ipft = 1
        CLM.InitPRTGlobalAllometricCarbon!()

        CLM.num_elements[] = 1
        empty!(CLM.element_list)
        push!(CLM.element_list, CLM.carbon12_element)

        CLM.hlm_parteh_mode[]       = CLM.prt_carbon_allom_hyp
        CLM.hlm_use_sp[]            = CLM.ifalse
        CLM.hlm_use_planthydro[]    = CLM.ifalse
        CLM.hlm_use_nocomp[]        = CLM.ifalse
        CLM.hlm_use_fixed_biogeog[] = CLM.ifalse
        CLM.hlm_use_potentialveg[]  = CLM.ifalse
        CLM.hlm_freq_day[]          = 1.0 / 365.0
        CLM.hlm_current_tod[]       = 0
        CLM.numpft[]                = npft
        CLM.hlm_name[]              = "CLM"
        # generous element capacity per patch/site (cohorts + counts)
        CLM.fates_maxElementsPerPatch[] = 100
        CLM.fates_maxElementsPerSite[]  = 100

        nlevsoil = 3

        # ==================================================================
        # 1. fates_restart_variable_type Init!/Flush! — real + int payloads
        # ==================================================================
        @testset "restart_variable_type Init/Flush (int + real)" begin
            # Build the 4 dim-kinds + 2 dim-bounds standalone (as the interface does).
            dim_kinds = [CLM.fates_io_variable_kind_type() for _ in 1:CLM.fates_restart_num_dim_kinds]
            CLM.Init!(dim_kinds[1], CLM.cohort_r8, 1)
            CLM.Init!(dim_kinds[2], CLM.site_r8, 1)
            CLM.Init!(dim_kinds[3], CLM.cohort_int, 1)
            CLM.Init!(dim_kinds[4], CLM.site_int, 1)

            dim_bounds = [CLM.fates_io_dimension_type() for _ in 1:CLM.fates_restart_num_dimensions]
            CLM.Init!(dim_bounds[1], CLM.dimname_cohort, 1, 1, 8)  # cohort dim: 1..8
            CLM.Init!(dim_bounds[2], CLM.dimname_column, 1, 1, 4)  # column dim: 1..4
            # wire dim-kinds to dim indices
            dim_kinds[1].dim1_index = 1; dim_kinds[1].dimsize[1] = 8
            dim_kinds[3].dim1_index = 1; dim_kinds[3].dimsize[1] = 8
            dim_kinds[2].dim1_index = 2; dim_kinds[2].dimsize[1] = 4
            dim_kinds[4].dim1_index = 2; dim_kinds[4].dimsize[1] = 4

            # real cohort var
            rvr = CLM.fates_restart_variable_type()
            CLM.Init!(rvr, "fates_dbh", "cm", "diameter", CLM.cohort_r8, CLM.flushzero,
                      CLM.fates_restart_num_dim_kinds, dim_kinds, dim_bounds)
            @test !rvr.is_int
            @test length(rvr.r81d) == 8
            @test isempty(rvr.int1d)
            @test all(rvr.r81d .== 0.0)
            rvr.r81d[3] = 9.9
            CLM.Flush!(rvr, 0, dim_bounds, dim_kinds)
            @test all(rvr.r81d .== 0.0)   # flushed back to flushval (0)

            # int site var with non-zero flushval (-9999)
            rvi = CLM.fates_restart_variable_type()
            CLM.Init!(rvi, "fates_PatchesPerSite", "none", "npatch", CLM.site_int,
                      CLM.flushinvalid, CLM.fates_restart_num_dim_kinds, dim_kinds, dim_bounds)
            @test rvi.is_int
            @test length(rvi.int1d) == 4
            @test isempty(rvi.r81d)
            @test all(rvi.int1d .== -9999)
            rvi.int1d[2] = 7
            CLM.Flush!(rvi, 0, dim_bounds, dim_kinds)
            @test all(rvi.int1d .== -9999)
        end

        # ==================================================================
        # 2. Interface Init + registry build
        # ==================================================================
        ri = _build_interface()

        @testset "registry build" begin
            @test CLM.num_restart_vars(ri) > 100         # non-empty registry
            @test length(ri.rvars) == CLM.num_restart_vars(ri)
            # the load-bearing count handles are registered and non-zero
            @test CLM.ir_get(ri, :ir_npatch_si) > 0
            @test CLM.ir_get(ri, :ir_ncohort_pa) > 0
            @test CLM.ir_get(ri, :ir_dbh_co) > 0
            @test CLM.ir_get(ri, :ir_pft_co) > 0
            @test CLM.ir_get(ri, :ir_area_pa) > 0
            # dim-kinds: PatchesPerSite is a SITE INT; dbh is a COHORT R8.
            npatch_var = ri.rvars[CLM.ir_get(ri, :ir_npatch_si)]
            @test npatch_var.vtype == CLM.site_int
            @test npatch_var.is_int
            dbh_var = ri.rvars[CLM.ir_get(ri, :ir_dbh_co)]
            @test dbh_var.vtype == CLM.cohort_r8
            @test !dbh_var.is_int
            # a couple of expected names exist verbatim
            names = [v.vname for v in ri.rvars]
            @test "fates_PatchesPerSite" in names
            @test "fates_CohortsPerPatch" in names
            @test "fates_dbh" in names
            @test "fates_nplant" in names
        end

        # ==================================================================
        # 3. DEMOGRAPHIC ROUND-TRIP
        # ==================================================================
        @testset "pack/unpack demographic round-trip" begin
            # --- source site: 2 patches (oldest=A age 20 / youngest=B age 5),
            #     A has 2 cohorts, B has 1 cohort. ------------------------------
            site = _build_restart_site(npft, nlevsoil)
            pA = _make_restart_patch(6000.0, 20.0, npft, nlevsoil)
            pB = _make_restart_patch(3000.0, 5.0,  npft, nlevsoil)
            pA.patchno = 1; pB.patchno = 2
            # build oldest->youngest list manually (A oldest, B youngest)
            site.oldest_patch   = pA
            site.youngest_patch = pB
            pA.younger = pB; pA.older = nothing
            pB.older   = pA; pB.younger = nothing

            cA1 = _add_restart_cohort!(site, pA, ipft, 25.0, 0.20; clayer=1)
            cA2 = _add_restart_cohort!(site, pA, ipft, 12.0, 0.50; clayer=2)
            cB1 = _add_restart_cohort!(site, pB, ipft, 8.0,  1.00; clayer=1)
            # set a few extra demographic fields to non-default to test round-trip
            for (c, dam, st) in ((cA1,1,2),(cA2,1,1),(cB1,1,2))
                c.crowndamage = dam; c.status_coh = st
                c.canopy_layer_yesterday = Float64(c.canopy_layer)
                c.size_class_lasttimestep = 2
            end

            # --- R1: cohort diagnostic scalars + mortality fluxes (seed values) -
            # walk all cohorts (pA shortest->taller, then pB) to seed unique values.
            all_src_coh = CLM.fates_cohort_type[]
            let c = pA.shortest
                while c !== nothing; push!(all_src_coh, c); c = c.taller; end
            end
            let c = pB.shortest
                while c !== nothing; push!(all_src_coh, c); c = c.taller; end
            end
            for (i, c) in enumerate(all_src_coh)
                c.gpp_acc = 1.0 * i; c.npp_acc = 2.0 * i; c.resp_acc = 3.0 * i
                c.gpp_acc_hold = 4.0 * i; c.npp_acc_hold = 5.0 * i; c.resp_acc_hold = 6.0 * i
                c.resp_excess = 7.0 * i; c.bmort = 0.1 * i; c.hmort = 0.2 * i
                c.cmort = 0.3 * i; c.smort = 0.4 * i; c.asmort = 0.5 * i; c.dgmort = 0.6 * i
                c.frmort = 0.7 * i; c.lmort_direct = 0.8 * i; c.lmort_collateral = 0.9 * i
                c.lmort_infra = 1.1 * i; c.ddbhdt = 1.2 * i; c.resp_tstep = 1.3 * i
            end

            # --- R2: patch fuel/scorch + ground albedo (seed values) -----------
            for (k, p) in enumerate((pA, pB))
                for i in 1:npft; p.scorch_ht[i] = 2.0 * k + i; end
                for i in 1:CLM.num_fuel_classes; p.fuel.effective_moisture[i] = 0.01 * (k + i); end
                for i in 1:CLM.num_swb
                    p.gnd_alb_dir[i] = 0.10 * k + 0.01 * i
                    p.gnd_alb_dif[i] = 0.20 * k + 0.01 * i
                end
            end

            # --- R4: patch litter pools (seed values) --------------------------
            for (k, p) in enumerate((pA, pB))
                for el in 1:CLM.num_elements[]
                    L = p.litter[el]
                    for i in 1:npft
                        L.seed[i] = 10.0 * k + i; L.seed_germ[i] = 20.0 * k + i
                        L.seed_decay[i] = 0.1 * k + i; L.seed_germ_decay[i] = 0.2 * k + i
                    end
                    for i in 1:CLM.ndcmpy
                        L.leaf_fines[i] = 30.0 * k + i; L.leaf_fines_frag[i] = 0.3 * k + i
                        for j in 1:nlevsoil
                            L.root_fines[i, j] = 40.0 * k + i + 0.5 * j
                            L.root_fines_frag[i, j] = 0.4 * k + i + 0.5 * j
                        end
                    end
                    for i in 1:CLM.ncwd
                        L.ag_cwd[i] = 50.0 * k + i; L.ag_cwd_frag[i] = 0.5 * k + i
                        for j in 1:nlevsoil
                            L.bg_cwd[i, j] = 60.0 * k + i + 0.5 * j
                            L.bg_cwd_frag[i, j] = 0.6 * k + i + 0.5 * j
                        end
                    end
                end
            end

            # capture source demographic snapshot
            src_patches = _patch_list(site)
            @test length(src_patches) == 2
            src_area  = [p.area for p in src_patches]
            src_age   = [p.age  for p in src_patches]
            src_cohA  = _cohort_tuples(pA)
            src_cohB  = _cohort_tuples(pB)
            @test length(src_cohA) == 2
            @test length(src_cohB) == 1

            # --- PACK -------------------------------------------------------
            CLM.set_restart_vectors!(ri, 1, 1, [site])

            # the counts must be in the flat vectors
            @test ri.rvars[CLM.ir_get(ri, :ir_npatch_si)].int1d[1] == 2
            # patch A base = slot 1, patch B base = slot 1 + maxElementsPerPatch
            base_b = 1 + CLM.fates_maxElementsPerPatch[]
            @test ri.rvars[CLM.ir_get(ri, :ir_ncohort_pa)].int1d[1] == 2      # patch A
            @test ri.rvars[CLM.ir_get(ri, :ir_ncohort_pa)].int1d[base_b] == 1 # patch B

            # --- ALLOCATE FRESH SKELETON from the stored counts -------------
            site2 = _build_restart_site(npft, nlevsoil)
            CLM.create_patchcohort_structure!(ri, 1, 1, [site2]; current_tod=0)

            dst_patches = _patch_list(site2)
            @test length(dst_patches) == 2                  # 2 patches rebuilt
            # cohort counts on the empty skeleton (via the linked list)
            @test length(_cohort_tuples(dst_patches[1])) == 2
            @test length(_cohort_tuples(dst_patches[2])) == 1

            # --- UNPACK -----------------------------------------------------
            CLM.get_restart_vectors!(ri, 1, 1, [site2])

            # --- ASSERT demographic state reconstructed identically ----------
            dst_area = [p.area for p in dst_patches]
            dst_age  = [p.age  for p in dst_patches]
            @test dst_area ≈ src_area
            @test dst_age  ≈ src_age
            @test dst_patches[1].land_use_label == pA.land_use_label
            @test dst_patches[2].land_use_label == pB.land_use_label
            @test dst_patches[1].nocomp_pft_label == pA.nocomp_pft_label

            # cohorts (shortest->taller) match (pft, dbh, n, height, canopy_layer)
            dst_cohA = _cohort_tuples(dst_patches[1])
            dst_cohB = _cohort_tuples(dst_patches[2])
            for (s, d) in zip(src_cohA, dst_cohA)
                @test s[1] == d[1]          # pft
                @test s[2] ≈ d[2]           # dbh
                @test s[3] ≈ d[3]           # n
                @test s[4] ≈ d[4]           # height
                @test s[5] == d[5]          # canopy_layer
            end
            for (s, d) in zip(src_cohB, dst_cohB)
                @test s[1] == d[1] && s[2] ≈ d[2] && s[3] ≈ d[3] && s[4] ≈ d[4] && s[5] == d[5]
            end
            # extra fields round-trip (crowndamage / status) on patch A cohorts
            dcoh = CLM.fates_cohort_type[]
            c = dst_patches[1].shortest
            while c !== nothing; push!(dcoh, c); c = c.taller; end
            @test all(c.crowndamage == 1 for c in dcoh)
            @test dcoh[1].status_coh == 2 || dcoh[1].status_coh == 1   # statuses preserved
            @test count(c -> c.status_coh == 2, dcoh) == 1
            @test count(c -> c.status_coh == 1, dcoh) == 1

            # site-level phenology round-trips
            @test site2.cstatus == site.cstatus
            @test site2.nchilldays == site.nchilldays
            @test site2.phen_model_date == site.phen_model_date
            @test site2.grow_deg_days ≈ site.grow_deg_days
            @test site2.snow_depth ≈ site.snow_depth

            # ============================================================
            # B18-followup R1/R2/R4/R5 round-trip assertions
            # ============================================================
            # R1: cohort diagnostic scalars + mortality fluxes.
            dco = CLM.fates_cohort_type[]
            c = dst_patches[1].shortest
            while c !== nothing; push!(dco, c); c = c.taller; end
            c = dst_patches[2].shortest
            while c !== nothing; push!(dco, c); c = c.taller; end
            # the source-order is (cA1, cA2, cB1); each cohort carries a unique
            # gpp_acc = 1*i. find by gpp_acc identity (robust to list order).
            gpps = sort([cc.gpp_acc for cc in dco])
            @test gpps ≈ [1.0, 2.0, 3.0]
            # for cohort with gpp_acc==2 (i=2 → cA2): hmort should be 0.2*2=0.4
            c2 = dco[findfirst(cc -> cc.gpp_acc ≈ 2.0, dco)]
            @test c2.hmort ≈ 0.4
            @test c2.resp_acc ≈ 6.0
            @test c2.lmort_direct ≈ 1.6
            @test c2.ddbhdt ≈ 2.4

            # R2: patch ground albedo + scorch height (pA is k=1, pB is k=2).
            @test dst_patches[1].gnd_alb_dir[1] ≈ 0.11
            @test dst_patches[2].gnd_alb_dif[2] ≈ 0.42
            @test dst_patches[1].scorch_ht[1]   ≈ 3.0   # 2*1 + 1
            @test dst_patches[1].fuel.effective_moisture[1] ≈ 0.02  # 0.01*(1+1)

            # R4: patch litter pool round-trip (pA, element 1).
            L1 = dst_patches[1].litter[1]
            @test L1.seed[1]       ≈ 11.0           # 10*1 + 1
            @test L1.ag_cwd[1]     ≈ 51.0           # 50*1 + 1
            @test L1.bg_cwd[1, 1]  ≈ 60.0 + 1 + 0.5 # 60*1 + 1 + 0.5*1
            @test L1.root_fines[1, nlevsoil] ≈ 40.0 + 1 + 0.5 * nlevsoil

            # R5: site demographic / mortality / flux array round-trip.
            @test site2.fmort_rate_canopy ≈ site.fmort_rate_canopy
            @test site2.imort_rate ≈ site.imort_rate
            @test site2.term_nindivs_canopy ≈ site.term_nindivs_canopy
            @test site2.demotion_rate ≈ site.demotion_rate
            @test site2.recruitment_rate ≈ site.recruitment_rate
            @test site2.dstatus == site.dstatus
            @test site2.ema_npp ≈ site.ema_npp
            @test site2.area_PFT ≈ site.area_PFT
            @test site2.landuse_vector_gt_min == site.landuse_vector_gt_min
            @test site2.liqvol_memory ≈ site.liqvol_memory
            @test site2.vegtemp_memory ≈ site.vegtemp_memory
            @test site2.mass_balance[1].old_stock ≈ site.mass_balance[1].old_stock
            @test site2.flux_diags.elem[1].cwd_ag_input ≈ site.flux_diags.elem[1].cwd_ag_input
            @test site2.iflux_balance[1].iflux_liveveg ≈ site.iflux_balance[1].iflux_liveveg
            @test site2.term_carbonflux_canopy ≈ site.term_carbonflux_canopy
            @test site2.term_abg_flux ≈ site.term_abg_flux
        end

        # ==================================================================
        # 4. R3 — update_3dpatch_radiation! runs + finite albd_parb
        # ==================================================================
        @testset "R3 update_3dpatch_radiation!" begin
            site = _build_restart_site(npft, nlevsoil)
            # bare-ground patch (no leaf layers): nrad zeroed, daylight on.
            p = _make_restart_patch(10000.0, 1.0, npft, nlevsoil)
            fill!(p.nrad, 0)
            p.solar_zenith_flag = true
            p.solar_zenith_angle = 0.5
            for i in 1:CLM.num_swb
                p.gnd_alb_dir[i] = 0.15 + 0.01 * i
                p.gnd_alb_dif[i] = 0.25 + 0.01 * i
            end
            site.oldest_patch = p; site.youngest_patch = p
            p.older = nothing; p.younger = nothing

            # minimal bc_out with the 7 radiation arrays sized (npatch x num_swb).
            bc_out = CLM.bc_out_type()
            npa = 1
            bc_out.albd_parb = fill(NaN, npa, CLM.num_swb)
            bc_out.albi_parb = fill(NaN, npa, CLM.num_swb)
            bc_out.fabd_parb = fill(NaN, npa, CLM.num_swb)
            bc_out.fabi_parb = fill(NaN, npa, CLM.num_swb)
            bc_out.ftdd_parb = fill(NaN, npa, CLM.num_swb)
            bc_out.ftid_parb = fill(NaN, npa, CLM.num_swb)
            bc_out.ftii_parb = fill(NaN, npa, CLM.num_swb)

            CLM.update_3dpatch_radiation!(ri, 1, [site], [bc_out])

            @test all(isfinite, bc_out.albd_parb)
            @test all(isfinite, bc_out.albi_parb)
            # bare-ground branch copies gnd albedos into albd/albi_parb
            @test bc_out.albd_parb[1, 1] ≈ p.gnd_alb_dir[1]
            @test bc_out.albi_parb[1, 2] ≈ p.gnd_alb_dif[2]
        end

    finally
        CLM.prt_global[]                 = old_global
        CLM.num_elements[]               = old_numel
        empty!(CLM.element_list); append!(CLM.element_list, old_ellist)
        CLM.hlm_parteh_mode[]            = old_parteh
        CLM.hlm_use_sp[]                 = old_use_sp
        CLM.hlm_use_planthydro[]         = old_hydro
        CLM.hlm_use_nocomp[]             = old_nocomp
        CLM.hlm_use_fixed_biogeog[]      = old_fixbio
        CLM.hlm_use_potentialveg[]       = old_potveg
        CLM.hlm_freq_day[]               = old_freqday
        CLM.hlm_current_tod[]            = old_tod
        CLM.numpft[]                     = old_numpft
        CLM.hlm_name[]                   = old_name
        CLM.fates_maxElementsPerPatch[]  = old_maxpp
        CLM.fates_maxElementsPerSite[]   = old_maxps
        CLM.nlevsclass[]                 = old_nlevsc
        CLM.nlevcoage[]                  = old_nlevca
        CLM.SFParams[]                   = old_sfp
        CLM.EDParams[]                   = old_edp
        CLM.ParamDerived[]               = old_pd
    end
end
