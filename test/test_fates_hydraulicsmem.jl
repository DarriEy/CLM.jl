# test_fates_hydraulicsmem.jl
# Tests for FATES plant-hydraulics memory/state types (Tier F, Batch 2).
#
# Verifies the compartment-index constants, then constructs + allocates both the
# site- and cohort-level hydraulics state with a small synthetic dimension set,
# asserting array shapes, zero/NaN/sentinel initialization, and the static
# node/connection topology produced by SetConnections! for both solver families.

using Test
using CLM

@testset "FATES HydraulicsMem" begin

    # =======================================================================
    # Compartment-index constants
    # =======================================================================
    @testset "constants" begin
        @test CLM.hydr_solver_1DTaylor == 1
        @test CLM.hydr_solver_2DPicard == 2
        @test CLM.hydr_solver_2DNewton == 3

        @test CLM.nlevsoi_hyd_max == 40
        @test CLM.n_porous_media == 5
        @test CLM.n_plant_media  == 4
        @test CLM.n_hypool_leaf  == 1
        @test CLM.n_hypool_stem  == 1
        @test CLM.n_hypool_troot == 1
        @test CLM.n_hypool_aroot == 1
        @test CLM.nshell         == 1

        @test CLM.n_hypool_ag    == CLM.n_hypool_leaf + CLM.n_hypool_stem
        @test CLM.n_hypool_tot   == CLM.n_hypool_ag + CLM.n_hypool_troot +
                                    CLM.n_hypool_aroot + CLM.nshell
        @test CLM.n_hypool_plant == CLM.n_hypool_tot - CLM.nshell

        # Porous-medium indices
        @test CLM.stomata_p_media == 0
        @test CLM.leaf_p_media    == 1
        @test CLM.stem_p_media    == 2
        @test CLM.troot_p_media   == 3
        @test CLM.aroot_p_media   == 4
        @test CLM.rhiz_p_media    == 5

        @test CLM.rwcft  == (1.0, 0.958, 0.958, 0.958)
        @test CLM.rwccap == (1.0, 0.947, 0.947, 0.947)
        @test CLM.fine_root_radius_const == 0.0001
    end

    # =======================================================================
    # Cohort-level allocate / copy
    # =======================================================================
    @testset "cohort allocate" begin
        nlevrhiz = 4
        ch = CLM.ed_cohort_hydr_type()

        # Fixed-size aboveground arrays present at construction, zero-init.
        @test length(ch.z_node_ag) == CLM.n_hypool_ag
        @test length(ch.th_ag)     == CLM.n_hypool_ag
        @test all(ch.th_ag .== 0.0)
        @test ch.btran == 0.0
        @test ch.is_newly_recruited == false

        # Per-layer arrays are empty before allocation.
        @test isempty(ch.th_aroot)
        @test isempty(ch.v_aroot_layer)

        CLM.AllocateHydrCohortArrays!(ch, nlevrhiz)
        for f in (:kmax_troot_lower, :kmax_aroot_upper, :kmax_aroot_lower,
                  :kmax_aroot_radial_in, :kmax_aroot_radial_out,
                  :v_aroot_layer_init, :v_aroot_layer, :l_aroot_layer,
                  :th_aroot, :psi_aroot, :ftc_aroot)
            v = getfield(ch, f)
            @test length(v) == nlevrhiz
            @test all(v .== 0.0)
        end

        # Copy semantics: independent storage.
        ch.th_troot = 3.5
        ch.th_aroot[1] = 9.0
        ch2 = CLM.ed_cohort_hydr_type()
        CLM.AllocateHydrCohortArrays!(ch2, nlevrhiz)
        CLM.CopyCohortHydraulics!(ch2, ch)
        @test ch2.th_troot == 3.5
        @test ch2.th_aroot[1] == 9.0
        ch2.th_aroot[1] = -1.0
        @test ch.th_aroot[1] == 9.0   # source unchanged

        CLM.DeallocateHydrCohortArrays!(ch)
        @test isempty(ch.th_aroot)
        @test isempty(ch.kmax_aroot_radial_out)
    end

    # =======================================================================
    # Site-level init: 1D (Taylor) solver
    # =======================================================================
    @testset "site init 1D" begin
        nlevrhiz = 5
        nlevsoil = 8
        numpft = 3
        numlevsclass = 6

        site = CLM.ed_site_hydr_type(nlevrhiz = nlevrhiz)
        CLM.InitHydrSite!(site, numpft, numlevsclass, CLM.hydr_solver_1DTaylor, nlevsoil)

        # Rhizosphere geometry shapes
        @test length(site.zi_rhiz) == nlevrhiz
        @test size(site.map_r2s) == (nlevrhiz, 2)
        @test size(site.v_shell) == (nlevrhiz, CLM.nshell)
        @test size(site.kmax_upper_shell) == (nlevrhiz, CLM.nshell)
        @test size(site.h2osoi_liqvol_shell) == (nlevrhiz, CLM.nshell)

        # Sentinel / NaN / constant init
        @test all(isnan, site.zi_rhiz)
        @test all(site.map_s2r .== -999)
        @test all(site.supsub_flag .== -999)
        @test all(site.rs1 .== CLM.fine_root_radius_const)
        @test all(site.rootl_sl .== 0.0)
        @test length(site.rootuptake_sl) == nlevsoil
        @test all(isnan, site.rootuptake_sl)

        # Diagnostics by size x pft
        @test size(site.sapflow_scpf) == (numlevsclass, numpft)
        @test size(site.rootuptake100_scpf) == (numlevsclass, numpft)

        # Scalars
        @test site.h2oveg == 0.0
        @test isnan(site.errh2o_hyd)

        # WTF holder vectors
        @test length(site.wrf_soil) == nlevrhiz
        @test length(site.wkf_soil) == nlevrhiz
        @test eltype(site.wrf_soil) == CLM.WRFType
        @test eltype(site.wkf_soil) == CLM.WKFType

        # Scratch arrays fixed length
        @test length(site.cohort_recruit_water_layer) == CLM.nlevsoi_hyd_max
        @test length(site.recruit_water_avail_layer) == CLM.nlevsoi_hyd_max

        # 1D solver: node/connection counts and topology
        @test site.num_nodes == CLM.n_hypool_leaf + CLM.n_hypool_stem +
                                CLM.n_hypool_troot + CLM.n_hypool_aroot + CLM.nshell
        @test site.num_connections == site.num_nodes - 1
        @test length(site.conn_up) == site.num_connections
        @test length(site.pm_node) == site.num_nodes
        @test site.pm_node[1] == CLM.leaf_p_media
        @test site.pm_node[CLM.n_hypool_ag] == CLM.stem_p_media
        @test site.pm_node[CLM.n_hypool_ag + 1] == CLM.troot_p_media
        @test site.pm_node[CLM.n_hypool_ag + 2] == CLM.aroot_p_media
        @test site.pm_node[CLM.n_hypool_ag + 3] == CLM.rhiz_p_media
    end

    # =======================================================================
    # Site-level init: 2D (Newton) solver + flush
    # =======================================================================
    @testset "site init 2D" begin
        nlevrhiz = 5
        nlevsoil = 8
        numpft = 3
        numlevsclass = 6

        site = CLM.ed_site_hydr_type(nlevrhiz = nlevrhiz)
        CLM.InitHydrSite!(site, numpft, numlevsclass, CLM.hydr_solver_2DNewton, nlevsoil)

        expected_nodes = CLM.n_hypool_leaf + CLM.n_hypool_stem + CLM.n_hypool_troot +
                         (CLM.n_hypool_aroot + CLM.nshell) * nlevrhiz
        expected_cnxs  = CLM.n_hypool_leaf + CLM.n_hypool_stem + CLM.n_hypool_troot - 1 +
                         (CLM.n_hypool_aroot + CLM.nshell) * nlevrhiz
        @test site.num_nodes == expected_nodes
        @test site.num_connections == expected_cnxs

        # Matrix-solver arrays sized
        @test size(site.ajac) == (expected_nodes, expected_nodes)
        @test length(site.residual) == expected_nodes
        @test length(site.node_layer) == expected_nodes
        @test length(site.q_flux) == expected_cnxs
        @test length(site.kmax_up) == expected_cnxs

        # node_layer topology: aboveground = 0, transporting root = 1, then layers
        @test all(site.node_layer[1:CLM.n_hypool_ag] .== 0)
        @test site.node_layer[CLM.n_hypool_ag + CLM.n_hypool_troot] == 1
        @test maximum(site.node_layer) == nlevrhiz

        # pm_node: troot node first, then each rhiz layer's junction node is aroot.
        @test site.pm_node[CLM.n_hypool_ag + CLM.n_hypool_troot] == CLM.troot_p_media
        @test site.pm_node[CLM.n_hypool_ag + CLM.n_hypool_troot + 1] == CLM.aroot_p_media

        # Flush resets scratch to fates_unset_r8 (2D only).
        CLM.FlushSiteScratch!(site, CLM.hydr_solver_2DNewton)
        @test all(site.residual .== CLM.fates_unset_r8)
        @test all(site.th_node .== CLM.fates_unset_r8)
        @test all(site.q_flux .== CLM.fates_unset_r8)
    end

    # =======================================================================
    # AggBCToRhiz harmonic mean
    # =======================================================================
    @testset "AggBCToRhiz" begin
        site = CLM.ed_site_hydr_type(nlevrhiz = 2)
        CLM.InitHydrSite!(site, 1, 1, CLM.hydr_solver_1DTaylor, 4)
        # Map rhiz layer 1 to soil layers 1:2.
        site.map_r2s[1, 1] = 1
        site.map_r2s[1, 2] = 2
        var_in = [2.0, 4.0, 99.0, 99.0]
        weight = [1.0, 1.0, 1.0, 1.0]
        # Harmonic mean of 2 and 4 with equal weights = 2/(1/2 + 1/4) = 2.6666...
        got = CLM.AggBCToRhiz(site, var_in, 1, weight)
        @test isapprox(got, 2.0 / (1.0/2.0 + 1.0/4.0); atol = 1.0e-12)
    end
end
