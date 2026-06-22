# test_fates_litter_radmem.jl
# Tests for FatesLitterMod (litter_type) and FatesRadiationMemMod (Tier F, Batch 1).

using Test
using CLM

@testset "FATES litter + radiation-mem" begin

    @testset "FatesLitterMod indices" begin
        @test CLM.ncwd == 4
        @test CLM.ndcmpy == 3
        @test CLM.ilabile == 1
        @test CLM.icellulose == 2
        @test CLM.ilignin == 3
    end

    @testset "litter_type InitAllocate! shapes + unset-init" begin
        numpft = 5
        numlevsoil = 7
        element_id = 11

        litt = CLM.litter_type()
        CLM.InitAllocate!(litt, numpft, numlevsoil, element_id)

        @test litt.element_id == element_id

        # CWD pools
        @test size(litt.ag_cwd) == (CLM.ncwd,)
        @test size(litt.ag_cwd_in) == (CLM.ncwd,)
        @test size(litt.ag_cwd_frag) == (CLM.ncwd,)
        @test size(litt.bg_cwd) == (CLM.ncwd, numlevsoil)
        @test size(litt.bg_cwd_in) == (CLM.ncwd, numlevsoil)
        @test size(litt.bg_cwd_frag) == (CLM.ncwd, numlevsoil)

        # fine litter (dcmpy)
        @test size(litt.leaf_fines) == (CLM.ndcmpy,)
        @test size(litt.leaf_fines_in) == (CLM.ndcmpy,)
        @test size(litt.leaf_fines_frag) == (CLM.ndcmpy,)
        @test size(litt.root_fines) == (CLM.ndcmpy, numlevsoil)
        @test size(litt.root_fines_in) == (CLM.ndcmpy, numlevsoil)
        @test size(litt.root_fines_frag) == (CLM.ndcmpy, numlevsoil)

        # seed pools (pft)
        @test size(litt.seed) == (numpft,)
        @test size(litt.seed_germ) == (numpft,)
        @test size(litt.seed_in_local) == (numpft,)
        @test size(litt.seed_in_extern) == (numpft,)
        @test size(litt.seed_germ_in) == (numpft,)
        @test size(litt.seed_germ_decay) == (numpft,)
        @test size(litt.seed_decay) == (numpft,)

        # everything set to the unset sentinel
        @test all(litt.ag_cwd .== CLM.fates_unset_r8)
        @test all(litt.bg_cwd .== CLM.fates_unset_r8)
        @test all(litt.leaf_fines .== CLM.fates_unset_r8)
        @test all(litt.root_fines .== CLM.fates_unset_r8)
        @test all(litt.seed .== CLM.fates_unset_r8)
        @test all(litt.seed_germ_in .== CLM.fates_unset_r8)
        @test all(litt.seed_decay .== CLM.fates_unset_r8)
    end

    @testset "ZeroFlux! zeroes flux pools only" begin
        litt = CLM.litter_type()
        CLM.InitAllocate!(litt, 4, 6, 1)
        CLM.ZeroFlux!(litt)

        @test all(litt.ag_cwd_in .== 0.0)
        @test all(litt.bg_cwd_in .== 0.0)
        @test all(litt.leaf_fines_in .== 0.0)
        @test all(litt.root_fines_in .== 0.0)
        @test all(litt.seed_in_local .== 0.0)
        @test all(litt.seed_in_extern .== 0.0)
        @test all(litt.ag_cwd_frag .== 0.0)
        @test all(litt.bg_cwd_frag .== 0.0)
        @test all(litt.leaf_fines_frag .== 0.0)
        @test all(litt.root_fines_frag .== 0.0)
        @test all(litt.seed_germ_in .== 0.0)
        @test all(litt.seed_decay .== 0.0)
        @test all(litt.seed_germ_decay .== 0.0)

        # prognostic pools untouched by ZeroFlux! (still unset)
        @test all(litt.ag_cwd .== CLM.fates_unset_r8)
        @test all(litt.seed .== CLM.fates_unset_r8)
    end

    @testset "InitConditions! + GetTotalLitterMass" begin
        numpft = 3
        numlevsoil = 4
        litt = CLM.litter_type()
        CLM.InitAllocate!(litt, numpft, numlevsoil, 1)
        CLM.InitConditions!(litt, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0)

        @test all(litt.leaf_fines .== 1.0)
        @test all(litt.root_fines .== 2.0)
        @test all(litt.ag_cwd .== 3.0)
        @test all(litt.bg_cwd .== 4.0)
        @test all(litt.seed .== 5.0)
        @test all(litt.seed_germ .== 6.0)

        # total = ag_cwd + bg_cwd + root_fines + leaf_fines + seed + seed_germ
        expected = 3.0 * CLM.ncwd +
                   4.0 * CLM.ncwd * numlevsoil +
                   2.0 * CLM.ndcmpy * numlevsoil +
                   1.0 * CLM.ndcmpy +
                   5.0 * numpft +
                   6.0 * numpft
        @test CLM.GetTotalLitterMass(litt) ≈ expected
    end

    @testset "FuseLitter! area-weighted + CopyLitter!" begin
        numpft = 2
        numlevsoil = 3
        a = CLM.litter_type(); CLM.InitAllocate!(a, numpft, numlevsoil, 1)
        b = CLM.litter_type(); CLM.InitAllocate!(b, numpft, numlevsoil, 1)
        CLM.InitConditions!(a, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
        CLM.ZeroFlux!(a)
        CLM.InitConditions!(b, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0)
        CLM.ZeroFlux!(b)

        # equal areas -> simple average of prognostic pools
        CLM.FuseLitter!(a, 1.0, 1.0, b)
        @test all(litt_pool ≈ 2.0 for litt_pool in a.ag_cwd)
        @test all(a.bg_cwd .≈ 2.0)
        @test all(a.seed .≈ 2.0)
        @test all(a.leaf_fines .≈ 2.0)
        @test all(a.root_fines .≈ 2.0)

        # CopyLitter! makes a value-equal copy
        c = CLM.litter_type(); CLM.InitAllocate!(c, numpft, numlevsoil, 1)
        CLM.CopyLitter!(c, a)
        @test c.ag_cwd == a.ag_cwd
        @test c.bg_cwd == a.bg_cwd
        @test c.seed == a.seed
        @test c.root_fines == a.root_fines
    end

    @testset "adjust_SF_CWD_frac partitioning" begin
        ncwd = CLM.ncwd
        frac = [0.4, 0.3, 0.2, 0.1]   # twig..trunk, sums to 1.0
        adj = zeros(ncwd)

        # large dbh (>7.6) -> unchanged
        CLM.adjust_SF_CWD_frac(10.0, ncwd, frac, adj)
        @test adj == frac

        # mid dbh (2.5 < dbh <= 7.6) -> trunk class (ncwd) zeroed, rest rescaled
        CLM.adjust_SF_CWD_frac(5.0, ncwd, frac, adj)
        @test adj[ncwd] == 0.0
        @test adj[ncwd-1] ≈ frac[ncwd-1] / (1.0 - frac[ncwd])
        @test sum(adj) ≈ 1.0  # rescaled back to unity

        # tiny dbh (<=0.6) -> everything to twigs
        CLM.adjust_SF_CWD_frac(0.3, ncwd, frac, adj)
        @test adj[ncwd] == 0.0
        @test adj[ncwd-1] == 0.0
        @test adj[ncwd-2] == 0.0
        @test adj[ncwd-3] ≈ sum(frac)
    end

    @testset "FatesRadiationMemMod indices + params" begin
        @test CLM.norman_solver == 1
        @test CLM.twostr_solver == 2
        @test CLM.num_rad_stream_types == 2
        @test CLM.idirect == 1
        @test CLM.idiffuse == 2
        @test CLM.num_swb == 2
        @test CLM.ivis == 1
        @test CLM.inir == 2
        @test CLM.ipar == CLM.ivis

        @test length(CLM.alb_ice) == CLM.num_swb
        @test CLM.alb_ice == [0.80, 0.55]
        @test CLM.rho_snow == [0.80, 0.55]
        @test CLM.tau_snow == [0.01, 0.01]
    end

end
