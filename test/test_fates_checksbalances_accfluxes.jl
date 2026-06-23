# test_fates_checksbalances_accfluxes.jl
#
# Tests for the FATES Tier-F Batch-11 modules:
#   * ChecksBalancesMod     — SiteMassStock / PatchMassStock site C/mass stock
#     summation + CheckIntegratedMassPools flux-vs-state balance check.
#   * EDAccumulateFluxesMod — AccumulateFluxes_ED daily flux accumulation.

using Test

@testset "FATES ChecksBalances + AccumulateFluxes" begin

    # ----------------------------------------------------------------------
    # Build a carbon-only PRT global descriptor with all SIX organs so that
    # GetState(organ, carbon12_element) resolves for the mass-stock sum.
    # ----------------------------------------------------------------------
    el = CLM.carbon12_element

    g = CLM.prt_global_type()
    CLM.ZeroGlobal!(g)
    g.hyp_name = "test_carbon_balance"
    g.hyp_id   = CLM.prt_carbon_allom_hyp
    g.num_vars     = 6
    g.num_bc_in    = 0
    g.num_bc_out   = 0
    g.num_bc_inout = 0
    g.state_descriptor = [CLM.state_descriptor_type() for _ in 1:g.num_vars]

    CLM.RegisterVarInGlobal!(g, 1, "leaf carbon",      "leaf_c",  CLM.leaf_organ,   el, 1)
    CLM.RegisterVarInGlobal!(g, 2, "fineroot carbon",  "fnrt_c",  CLM.fnrt_organ,   el, 1)
    CLM.RegisterVarInGlobal!(g, 3, "sapwood carbon",   "sapw_c",  CLM.sapw_organ,   el, 1)
    CLM.RegisterVarInGlobal!(g, 4, "storage carbon",   "store_c", CLM.store_organ,  el, 1)
    CLM.RegisterVarInGlobal!(g, 5, "repro carbon",     "repro_c", CLM.repro_organ,  el, 1)
    CLM.RegisterVarInGlobal!(g, 6, "structure carbon", "str_c",   CLM.struct_organ, el, 1)

    # Install the global PRT descriptor + a single-element FATES element set.
    old_global       = CLM.prt_global[]
    old_num_elements = CLM.num_elements[]
    old_element_list = copy(CLM.element_list)
    old_numpft       = CLM.numpft[]

    CLM.prt_global[]     = g
    CLM.num_elements[]   = 1
    empty!(CLM.element_list)
    push!(CLM.element_list, CLM.carbon12_element)
    CLM.numpft[]         = 2

    try
        numpft     = CLM.numpft[]
        numlevsoil = 3

        # Helper: build a cohort whose 6 organs carry a known per-organ mass.
        function make_cohort(; leaf, fnrt, sapw, store, repro, struc, n)
            prt = CLM.prt_vartypes()
            CLM.InitPRTVartype!(prt)
            CLM.SetState!(prt, CLM.leaf_organ,   el, leaf)
            CLM.SetState!(prt, CLM.fnrt_organ,   el, fnrt)
            CLM.SetState!(prt, CLM.sapw_organ,   el, sapw)
            CLM.SetState!(prt, CLM.store_organ,  el, store)
            CLM.SetState!(prt, CLM.repro_organ,  el, repro)
            CLM.SetState!(prt, CLM.struct_organ, el, struc)
            coh = CLM.fates_cohort_type()
            coh.prt = prt
            coh.n   = n
            return coh
        end

        # Helper: build a litter object with known stocks.
        function make_litter(; ag_cwd, bg_cwd, leaf_fines, root_fines, seed, seed_germ)
            litt = CLM.litter_type()
            litt.element_id = CLM.carbon12_element
            litt.ag_cwd     = fill(ag_cwd, CLM.ncwd)
            litt.bg_cwd     = fill(bg_cwd, CLM.ncwd, numlevsoil)
            litt.leaf_fines = fill(leaf_fines, CLM.ndcmpy)
            litt.root_fines = fill(root_fines, CLM.ndcmpy, numlevsoil)
            litt.seed       = fill(seed, numpft)
            litt.seed_germ  = fill(seed_germ, numpft)
            return litt
        end

        @testset "PatchMassStock + SiteMassStock" begin
            # --- Patch 1: two cohorts ---
            c1 = make_cohort(leaf=1.0, fnrt=0.5, sapw=2.0, store=0.25, repro=0.1, struc=4.0, n=10.0)
            c2 = make_cohort(leaf=0.5, fnrt=0.2, sapw=1.0, store=0.10, repro=0.0, struc=2.0, n=5.0)
            # tallest -> shorter linked list: c1 (tallest) -> c2
            c1.shorter = c2

            p1 = CLM.fates_patch_type()
            p1.area    = 100.0
            p1.tallest = c1
            p1.litter  = [make_litter(ag_cwd=0.1, bg_cwd=0.05, leaf_fines=0.02,
                                      root_fines=0.01, seed=0.5, seed_germ=0.25)]

            # Hand-summed expectations for patch 1
            organ_sum_c1 = 1.0 + 0.5 + 2.0 + 0.25 + 0.1 + 4.0
            organ_sum_c2 = 0.5 + 0.2 + 1.0 + 0.10 + 0.0 + 2.0
            live1 = organ_sum_c1 * 10.0 + organ_sum_c2 * 5.0
            litter1 = p1.area * (CLM.ncwd*0.1 + CLM.ncwd*numlevsoil*0.05 +
                                 CLM.ndcmpy*0.02 + CLM.ndcmpy*numlevsoil*0.01)
            seed1 = p1.area * (numpft*0.5 + numpft*0.25)

            pb, ps, pl = CLM.PatchMassStock(p1, 1)
            @test pb ≈ live1
            @test ps ≈ seed1
            @test pl ≈ litter1

            # --- Patch 2: one cohort ---
            c3 = make_cohort(leaf=2.0, fnrt=1.0, sapw=3.0, store=0.5, repro=0.2, struc=6.0, n=2.0)
            p2 = CLM.fates_patch_type()
            p2.area    = 50.0
            p2.tallest = c3
            p2.litter  = [make_litter(ag_cwd=0.2, bg_cwd=0.0, leaf_fines=0.0,
                                      root_fines=0.0, seed=0.0, seed_germ=0.0)]

            organ_sum_c3 = 2.0 + 1.0 + 3.0 + 0.5 + 0.2 + 6.0
            live2 = organ_sum_c3 * 2.0
            litter2 = p2.area * (CLM.ncwd*0.2)
            seed2 = 0.0

            # --- Site: oldest_patch -> younger linked list (p1 oldest -> p2) ---
            site = CLM.ed_site_type()
            site.oldest_patch = p1
            p1.younger = p2

            tot, bio, lit, seed = CLM.SiteMassStock(site, 1)
            @test bio  ≈ live1 + live2
            @test lit  ≈ litter1 + litter2
            @test seed ≈ seed1 + seed2
            @test tot  ≈ bio + seed + lit
            # total equals the fully hand-summed grand total
            @test tot  ≈ (live1 + live2) + (seed1 + seed2) + (litter1 + litter2)
        end

        @testset "CheckIntegratedMassPools: balanced vs imbalanced" begin
            # Single patch, single cohort, simple known stocks.
            coh = make_cohort(leaf=1.0, fnrt=0.0, sapw=0.0, store=0.0, repro=0.0, struc=0.0, n=1.0)
            patch = CLM.fates_patch_type()
            patch.area    = CLM.area              # full site area => area*area_inv == 1
            patch.tallest = coh
            patch.litter  = [make_litter(ag_cwd=0.0, bg_cwd=0.0, leaf_fines=0.0,
                                         root_fines=0.0, seed=0.0, seed_germ=0.0)]

            site = CLM.ed_site_type()
            site.oldest_patch = patch

            # Per-element accounting structures (one element).
            site.iflux_balance = [CLM.site_ifluxbal_type()]
            site.mass_balance  = [CLM.site_massbal_type()]
            ediag = CLM.elem_diag_type()
            ediag.surf_fine_litter_input = zeros(CLM.ndcmpy)
            ediag.root_litter_input      = zeros(CLM.ndcmpy)
            site.flux_diags = CLM.site_fluxdiags_type(elem = [ediag])

            # --- Balanced case ---
            # state_liveveg = (biomass + seed)*area_inv. biomass = 1*1 = 1.0, seed = 0.
            # area = full site so area_inv*area = 1, biomass_stock = live*area? No:
            # live_stock is NOT multiplied by area in PatchMassStock. biomass=1.0.
            biomass = 1.0
            # Set the integrated live-veg flux to exactly match the state so err=0.
            # iflux_liveveg starts 0, gets += net_uptake (=npp) with all else 0.
            site.flux_diags.npp = biomass * CLM.area_inv   # net_uptake = npp (carbon)
            CLM.CheckIntegratedMassPools(site)

            @test site.iflux_balance[1].state_liveveg ≈ biomass * CLM.area_inv
            @test site.flux_diags.elem[1].err_liveveg ≈ 0.0 atol=1e-15
            @test site.flux_diags.elem[1].err_litter  ≈ 0.0 atol=1e-15

            # --- Imbalanced case ---
            # Reset the integrated balances + diagnostics, then feed a flux that
            # does NOT match the state, and assert the error term is flagged.
            site.iflux_balance[1] = CLM.site_ifluxbal_type()
            ediag2 = CLM.elem_diag_type()
            ediag2.surf_fine_litter_input = zeros(CLM.ndcmpy)
            ediag2.root_litter_input      = zeros(CLM.ndcmpy)
            site.flux_diags = CLM.site_fluxdiags_type(elem = [ediag2])
            site.flux_diags.npp = 0.0   # net uptake 0 but biomass stock is 1.0 -> imbalance
            CLM.CheckIntegratedMassPools(site)

            @test site.iflux_balance[1].state_liveveg ≈ biomass * CLM.area_inv
            @test site.iflux_balance[1].iflux_liveveg ≈ 0.0
            # err = iflux - state = 0 - biomass*area_inv  => nonzero, flagged
            @test site.flux_diags.elem[1].err_liveveg ≈ -biomass * CLM.area_inv
            @test abs(site.flux_diags.elem[1].err_liveveg) > 1e-12
        end

        @testset "AccumulateFluxes_ED daily accumulation" begin
            # One cohort accumulated over N steps with a known per-step flux.
            coh = CLM.fates_cohort_type()
            coh.nv = 2
            coh.npp_acc  = 0.0
            coh.gpp_acc  = 0.0
            coh.resp_acc = 0.0
            coh.sym_nfix_daily = 0.0
            coh.c13disc_acc = 0.0
            coh.year_net_uptake = fill(999.0, coh.nv)   # sentinel: leaves present
            coh.ts_net_uptake   = [0.1, 0.2]

            patch = CLM.fates_patch_type()
            patch.nocomp_pft_label = 1                  # not bareground
            patch.shortest = coh                        # single cohort (taller -> nothing)

            site = CLM.ed_site_type()
            site.oldest_patch = patch

            bc_in  = CLM.bc_in_type()
            bc_in.filter_photo_pa = [3]                 # ifp 1 -> run
            bc_out = CLM.bc_out_type()

            sites   = [site]
            bcin_v  = [bc_in]
            bcout_v = [bc_out]

            nstep = 4
            npp_step, gpp_step, resp_step, nfix_step = 0.3, 0.5, 0.2, 0.05
            c13_step = 10.0
            # Mirror the in-routine GPP-weighted-mean recurrence exactly. Note the
            # Fortran updates gpp_acc BEFORE the c13 weighting, so the weighting
            # uses the already-incremented gpp_acc + the same-step gpp_tstep.
            exp_gpp_acc = 0.0
            exp_c13 = 0.0
            for _ in 1:nstep
                exp_gpp_acc += gpp_step
                exp_c13 = ((exp_c13 * exp_gpp_acc) + (c13_step * gpp_step)) /
                          (exp_gpp_acc + gpp_step)

                coh.npp_tstep  = npp_step
                coh.gpp_tstep  = gpp_step
                coh.resp_tstep = resp_step
                coh.sym_nfix_tstep = nfix_step
                coh.c13disc_clm = c13_step
                CLM.AccumulateFluxes_ED(1, sites, bcin_v, bcout_v, 1800.0)
            end

            @test coh.npp_acc  ≈ nstep * npp_step
            @test coh.gpp_acc  ≈ nstep * gpp_step
            @test coh.resp_acc ≈ nstep * resp_step
            @test coh.sym_nfix_daily ≈ nstep * nfix_step
            # GPP-weighted running mean of a constant c13disc_clm (matches the
            # in-routine recurrence exactly).
            @test coh.c13disc_acc ≈ exp_c13
            # year_net_uptake: sentinel reset to 0 on first step, then += ts each step.
            @test coh.year_net_uptake[1] ≈ nstep * 0.1
            @test coh.year_net_uptake[2] ≈ nstep * 0.2

            # --- Bare-ground patch + filter != 3 are skipped ---
            coh2 = CLM.fates_cohort_type()
            coh2.nv = 1
            coh2.npp_acc = 0.0; coh2.gpp_acc = 0.0; coh2.resp_acc = 0.0
            coh2.sym_nfix_daily = 0.0; coh2.c13disc_acc = 0.0
            coh2.year_net_uptake = [999.0]; coh2.ts_net_uptake = [0.5]
            coh2.npp_tstep = 1.0; coh2.gpp_tstep = 1.0; coh2.resp_tstep = 1.0
            coh2.sym_nfix_tstep = 1.0; coh2.c13disc_clm = 1.0

            bare = CLM.fates_patch_type()
            bare.nocomp_pft_label = CLM.nocomp_bareground   # bareground -> skipped, ifp not incremented
            bare.shortest = coh2

            site2 = CLM.ed_site_type()
            site2.oldest_patch = bare
            bc_in2 = CLM.bc_in_type()
            bc_in2.filter_photo_pa = Int[]   # nothing to index (patch skipped)
            CLM.AccumulateFluxes_ED(1, [site2], [bc_in2], [CLM.bc_out_type()], 1800.0)
            @test coh2.npp_acc ≈ 0.0          # untouched (bareground skipped)
            @test coh2.year_net_uptake[1] == 999.0
        end

    finally
        CLM.prt_global[]   = old_global
        CLM.num_elements[] = old_num_elements
        empty!(CLM.element_list)
        append!(CLM.element_list, old_element_list)
        CLM.numpft[]       = old_numpft
    end
end
