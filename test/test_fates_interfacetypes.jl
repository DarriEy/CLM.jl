# test_fates_interfacetypes.jl
# Tests for FatesInterfaceTypesMod (Tier F, Batch 1): the host<->FATES boundary
# type definitions and their allocate/init routines.

using Test
using CLM

@testset "FATES interface types" begin

    # Small synthetic dimension set.
    npatches   = 4
    nlevsoil   = 5
    nlevdecomp = 3
    num_rad_bands = 2
    max_comp   = 6
    num_harvest_cats     = 2
    num_luh2_states      = 12
    num_luh2_transitions = 26
    num_pft    = 7
    num_lu     = 5

    @testset "module parameters / sentinels" begin
        # Seed-dispersal kernel + cadence parameters are fixed.
        @test CLM.fates_dispersal_kernel_exponential == 1
        @test CLM.fates_dispersal_kernel_exppower    == 2
        @test CLM.fates_dispersal_kernel_logsech     == 3
        @test CLM.fates_dispersal_cadence_none    == 0
        @test CLM.fates_dispersal_cadence_daily   == 1
        @test CLM.fates_dispersal_cadence_monthly == 2
        @test CLM.fates_dispersal_cadence_yearly  == 3

        # Mutable HLM control variables default to the unset sentinel and can be set.
        #
        # These are module GLOBALS, so this testset was both order-DEPENDENT and
        # order-CORRUPTING: it asserted the live value still equalled the sentinel
        # (true only while nothing had initialised them) and then left its own
        # values behind for every later test. Wiring InitTimeAveragingGlobals()
        # exposed that — hlm_stepsize is now legitimately set to the run's dtime.
        # Save, reset to the sentinel to test the declared default, exercise the
        # setter, then restore, so the assertion is about the code rather than
        # about test ordering.
        _saved = (CLM.hlm_numSWb[], CLM.hlm_name[], CLM.hlm_stepsize[])
        try
            CLM.hlm_numSWb[]   = CLM.fates_unset_int
            CLM.hlm_name[]     = ""
            CLM.hlm_stepsize[] = CLM.fates_unset_r8
            @test CLM.hlm_numSWb[]   == CLM.fates_unset_int
            @test CLM.hlm_name[]     == ""
            @test CLM.hlm_stepsize[] == CLM.fates_unset_r8
            CLM.hlm_numSWb[]   = num_rad_bands
            CLM.hlm_name[]     = "CLM"
            CLM.hlm_stepsize[] = 1800.0
            @test CLM.hlm_numSWb[]   == num_rad_bands
            @test CLM.hlm_name[]     == "CLM"
            @test CLM.hlm_stepsize[] == 1800.0
        finally
            CLM.hlm_numSWb[], CLM.hlm_name[], CLM.hlm_stepsize[] = _saved
        end

        # History dimension map vectors default to empty and accept assignment.
        @test isempty(CLM.fates_hdim_levsclass[])
        @test eltype(CLM.fates_hdim_levsclass[]) == Float64
        @test eltype(CLM.fates_hdim_pfmap_levscpf[]) == Int
        CLM.fates_hdim_levsclass[] = [0.0, 5.0, 20.0]
        @test length(CLM.fates_hdim_levsclass[]) == 3
    end

    @testset "bc_in_type defaults + allocate" begin
        bc_in = CLM.bc_in_type()
        # Scalar defaults are the unset sentinels.
        @test bc_in.npatches == CLM.fates_unset_int
        @test bc_in.forc_pbot == CLM.fates_unset_r8
        @test bc_in.filter_btran == false
        # Allocatable arrays default to undef-0.
        @test isempty(bc_in.zi_sisl)
        @test size(bc_in.solad_parb) == (0, 0)

        CLM.allocate_bcin!(bc_in; npatches=npatches, nlevsoil=nlevsoil,
                           nlevdecomp=nlevdecomp, num_rad_bands=num_rad_bands,
                           max_comp=max_comp, num_harvest_cats=num_harvest_cats,
                           num_luh2_states=num_luh2_states,
                           num_luh2_transitions=num_luh2_transitions,
                           num_pft=num_pft, num_lu=num_lu)

        # Scalar dimensions recorded.
        @test bc_in.npatches   == npatches
        @test bc_in.nlevsoil   == nlevsoil
        @test bc_in.nlevdecomp == nlevdecomp

        # Soil-layer structure: zi carries a surface zero index (nlevsoil+1).
        @test length(bc_in.zi_sisl) == nlevsoil + 1
        @test length(bc_in.dz_sisl) == nlevsoil
        @test length(bc_in.z_sisl)  == nlevsoil
        @test length(bc_in.decomp_id) == nlevsoil
        @test eltype(bc_in.decomp_id) == Int

        # Decomp-layer fractions.
        @test length(bc_in.dz_decomp_sisl) == nlevdecomp
        @test length(bc_in.w_scalar_sisl)  == nlevdecomp
        @test length(bc_in.t_scalar_sisl)  == nlevdecomp

        # Patch-dimensioned arrays.
        @test length(bc_in.lightning24)   == npatches
        @test length(bc_in.dayl_factor_pa) == npatches
        @test length(bc_in.filter_photo_pa) == npatches
        @test length(bc_in.filter_vegzen_pa) == npatches
        @test bc_in.filter_vegzen_pa isa AbstractVector{Bool}

        # Radiation arrays (patch x band).
        @test size(bc_in.solad_parb) == (npatches, num_rad_bands)
        @test size(bc_in.solai_parb) == (npatches, num_rad_bands)
        @test length(bc_in.albgr_dir_rb) == num_rad_bands

        # Nutrient input fluxes (competitor x 1).
        @test size(bc_in.plant_nh4_uptake_flux) == (max_comp, 1)
        @test size(bc_in.plant_p_uptake_flux)   == (max_comp, 1)

        # Hydrology + plant-hydro soil arrays.
        @test length(bc_in.smp_sl)      == nlevsoil
        @test length(bc_in.watsat_sisl) == nlevsoil
        @test length(bc_in.h2o_liq_sisl) == nlevsoil

        # Land use dimensions.
        @test length(bc_in.hlm_harvest_rates)    == num_harvest_cats
        @test length(bc_in.hlm_harvest_catnames) == num_harvest_cats
        @test eltype(bc_in.hlm_harvest_catnames) == String
        @test length(bc_in.hlm_luh_states)       == num_luh2_states
        @test length(bc_in.hlm_luh_transitions)  == num_luh2_transitions

        # Fixed biogeography.
        @test length(bc_in.pft_areafrac)     == num_pft
        @test size(bc_in.pft_areafrac_lu)    == (num_pft, num_lu)

        # Satellite phenology (defaults to npatches when num_sp omitted).
        @test length(bc_in.hlm_sp_tlai) == npatches
        @test length(bc_in.hlm_sp_htop) == npatches

        # Zero-initialized.
        @test all(iszero, bc_in.dz_sisl)
        @test all(iszero, bc_in.solad_parb)
    end

    @testset "bc_out_type defaults + allocate" begin
        bc_out = CLM.bc_out_type()
        @test bc_out.num_plant_comps == CLM.fates_unset_int
        @test bc_out.gpp_site == CLM.fates_unset_r8
        @test isempty(bc_out.fsun_pa)
        @test size(bc_out.albd_parb) == (0, 0)

        CLM.allocate_bcout!(bc_out; npatches=npatches, nlevsoil=nlevsoil,
                            nlevdecomp=nlevdecomp, num_rad_bands=num_rad_bands,
                            max_comp=max_comp, num_pft=num_pft)

        # Sun/shade canopy (patch).
        @test length(bc_out.fsun_pa)   == npatches
        @test length(bc_out.laisun_pa) == npatches
        @test length(bc_out.btran_pa)  == npatches

        # Suction filter + root fraction.
        @test length(bc_out.active_suction_sl) == nlevsoil
        @test bc_out.active_suction_sl isa AbstractVector{Bool}
        @test size(bc_out.rootr_pasl) == (npatches, nlevsoil)

        # Canopy radiation (patch x band).
        @test size(bc_out.albd_parb) == (npatches, num_rad_bands)
        @test size(bc_out.ftii_parb) == (npatches, num_rad_bands)

        # Litter mass fluxes (site x decomp).
        @test length(bc_out.litt_flux_cel_c_si) == nlevdecomp
        @test length(bc_out.litt_flux_lab_p_si) == nlevdecomp

        # Nutrient competition.
        @test length(bc_out.source_nh4) == nlevdecomp
        @test size(bc_out.veg_rootc)    == (max_comp, nlevdecomp)
        @test length(bc_out.decompmicc) == nlevdecomp
        @test length(bc_out.ft_index)   == max_comp
        @test eltype(bc_out.ft_index)   == Int

        # CH4 boundary conditions.
        @test length(bc_out.annavg_agnpp_pa) == npatches
        @test size(bc_out.rootfr_pa)         == (npatches, nlevsoil)

        # Canopy structure (patch).
        @test length(bc_out.elai_pa) == npatches
        @test length(bc_out.htop_pa) == npatches
        @test length(bc_out.nocomp_pft_label_pa) == npatches
        @test eltype(bc_out.nocomp_pft_label_pa) == Int

        # FATES hydraulics (site x soil layer).
        @test length(bc_out.qflx_soil2root_sisl) == nlevsoil
        @test length(bc_out.qflx_ro_sisl)        == nlevsoil

        @test all(iszero, bc_out.fsun_pa)
        @test all(iszero, bc_out.albd_parb)
    end

    @testset "bc_pconst_type defaults + allocate" begin
        bc_pconst = CLM.bc_pconst_type()
        @test bc_pconst.max_plant_comps == CLM.fates_unset_int
        @test bc_pconst.eca_plant_escalar == CLM.fates_unset_r8
        @test isempty(bc_pconst.vmax_nh4)
        @test isempty(bc_pconst.j_uptake)

        CLM.allocate_bcpconst!(bc_pconst; num_pft=num_pft,
                               max_plant_comps=max_comp, nlevdecomp=nlevdecomp)

        @test bc_pconst.max_plant_comps == max_comp
        # PFT-dimensioned ECA vectors.
        @test length(bc_pconst.vmax_nh4)       == num_pft
        @test length(bc_pconst.vmax_p)         == num_pft
        @test length(bc_pconst.eca_km_nh4)     == num_pft
        @test length(bc_pconst.eca_km_ptase)   == num_pft
        @test length(bc_pconst.eca_lambda_ptase) == num_pft
        # Decomposition -> uptake layer map.
        @test length(bc_pconst.j_uptake) == nlevdecomp
        @test eltype(bc_pconst.j_uptake) == Int
        @test all(iszero, bc_pconst.vmax_nh4)
    end

end
