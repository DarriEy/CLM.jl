# test_fates_interface.jl
# Tests for FatesInterfaceMod (Tier F, Batch 18): the CLM<->FATES coupling point.
#
# Covers:
#   * set_fates_ctrlparms — flush + round-trips the hlm_* control flags
#     (ival / rval / cval dispatch).
#   * SetFatesTime — sets the FATES clock module-globals.
#   * allocate_bcin / allocate_bcout / zero_bcs — correctly-sized + zeroed bc
#     structures, driven through the fates_interface_type container.
#   * allocate_bcpconst / set_bcpconst — PFT-dimensioned ECA param constants.
#   * SetFatesGlobalElements2 (use_fates=false) + fates_history_maps +
#     InitPARTEHGlobals — sane global dimension counts and history maps.
#   * DetermineGridCellNeighbors — sane neighbor list on a small grid.
#
# fates_hist: NOT referenced by this module (the Fortran `use ... only :
# fates_hist` import is vestigial — never called), so nothing here needs the
# absent history module.

using Test
using CLM

@testset "FATES Batch 18: FatesInterfaceMod" begin

    # =======================================================================
    # set_fates_ctrlparms — flush + round-trip the control flags.
    # =======================================================================
    @testset "set_fates_ctrlparms round-trip" begin
        CLM.set_fates_ctrlparms("flush_to_unset")
        # After a flush, the unset sentinel (-999) is in place.
        @test CLM.hlm_numSWb[]   == -999
        @test CLM.hlm_name[]     == "unset"
        @test CLM.hlm_parteh_mode[] == -999
        @test CLM.hlm_use_planthydro[] == -999

        # Integer tags.
        CLM.set_fates_ctrlparms("masterproc";   ival = 1)
        CLM.set_fates_ctrlparms("num_sw_bbands"; ival = 2)
        CLM.set_fates_ctrlparms("vis_sw_index"; ival = 1)
        CLM.set_fates_ctrlparms("nir_sw_index"; ival = 2)
        CLM.set_fates_ctrlparms("parteh_mode";  ival = 1)
        CLM.set_fates_ctrlparms("use_planthydro"; ival = 0)
        CLM.set_fates_ctrlparms("use_sp";       ival = 0)
        CLM.set_fates_ctrlparms("num_lev_soil"; ival = 10)
        @test CLM.hlm_masterproc[]     == 1
        @test CLM.hlm_numSWb[]         == 2
        @test CLM.hlm_ivis[]           == 1
        @test CLM.hlm_inir[]           == 2
        @test CLM.hlm_parteh_mode[]    == 1
        @test CLM.hlm_use_planthydro[] == 0
        @test CLM.hlm_use_sp[]         == 0
        @test CLM.hlm_maxlevsoil[]     == 10

        # Real tag.
        CLM.set_fates_ctrlparms("hio_ignore_val"; rval = -9.9e30)
        @test CLM.hlm_hio_ignore_val[] == -9.9e30

        # Character tags.
        CLM.set_fates_ctrlparms("hlm_name";      cval = "CLM")
        CLM.set_fates_ctrlparms("nu_com";        cval = "RD")
        CLM.set_fates_ctrlparms("decomp_method"; cval = "CENTURY")
        @test CLM.hlm_name[]   == "CLM"
        @test CLM.hlm_nu_com[] == "RD"
        @test CLM.hlm_decomp[] == "CENTURY"

        # An unrecognized tag is non-fatal (just a log line), leaves state alone.
        prev = CLM.hlm_numSWb[]
        CLM.set_fates_ctrlparms("a_bogus_tag"; ival = 7)
        @test CLM.hlm_numSWb[] == prev
    end

    # =======================================================================
    # SetFatesTime — set the clock.
    # =======================================================================
    @testset "SetFatesTime" begin
        CLM.SetFatesTime(2026, 6, 23, 43200, 20260623, 20000101,
                         9670.5, 174, 365, 1.0 / 365.0)
        @test CLM.hlm_current_year[]   == 2026
        @test CLM.hlm_current_month[]  == 6
        @test CLM.hlm_current_day[]    == 23
        @test CLM.hlm_current_tod[]    == 43200
        @test CLM.hlm_current_date[]   == 20260623
        @test CLM.hlm_reference_date[] == 20000101
        @test CLM.hlm_model_day[]      == 9670.5
        @test CLM.hlm_day_of_year[]    == 174
        @test CLM.hlm_days_per_year[]  == 365
        @test CLM.hlm_freq_day[]       ≈ 1.0 / 365.0
    end

    # =======================================================================
    # allocate_bcin / allocate_bcout / zero_bcs via the interface container.
    # =======================================================================
    @testset "allocate + zero bcs" begin
        # Configure the module globals the allocators read.
        numpft_test = 3
        nlevsoil    = 6
        nlevdecomp  = 1     # non-vertsoilc => 1
        CLM.numpft[]            = numpft_test
        CLM.hlm_use_vertsoilc[] = CLM.ifalse
        CLM.hlm_use_planthydro[] = CLM.ifalse
        CLM.hlm_parteh_mode[]   = CLM.prt_carbon_allom_hyp
        CLM.hlm_use_lu_harvest[] = 0
        CLM.hlm_use_luh[]       = CLM.ifalse
        CLM.hlm_use_fixed_biogeog[] = CLM.ifalse
        CLM.hlm_use_sp[]        = CLM.ifalse
        CLM.max_comp_per_site[] = 1
        CLM.fates_maxPatchesPerSite[] = 4

        # maxpatch_total lives on ed_params; size it for the test.
        edp = CLM.ed_params_type()
        edp.maxpatch_total = 3
        CLM.EDParams[] = edp
        @test CLM.maxpatch_total() == 3

        fates = CLM.fates_interface_type(
            nsites = 1,
            bc_in  = [CLM.bc_in_type()],
            bc_out = [CLM.bc_out_type()],
        )

        CLM.allocate_bcin(fates.bc_in[1], nlevsoil, nlevdecomp,
                          0, 0, 0, 1, numpft_test)
        CLM.allocate_bcout(fates.bc_out[1], nlevsoil, nlevdecomp)

        bc_in  = fates.bc_in[1]
        bc_out = fates.bc_out[1]

        # Soil-layer + decomp dims recorded + sized.
        @test bc_in.nlevsoil   == nlevsoil
        @test bc_in.nlevdecomp == nlevdecomp
        @test length(bc_in.zi_sisl)  == nlevsoil + 1   # surface zero index
        @test length(bc_in.smp_sl)   == nlevsoil
        @test length(bc_in.w_scalar_sisl) == nlevdecomp

        # Patch-dimensioned arrays use maxpatch_total (3).
        @test length(bc_in.lightning24) == 3
        @test size(bc_in.solad_parb)    == (3, CLM.num_swb)
        @test length(bc_out.elai_pa)    == 3
        @test size(bc_out.albd_parb)    == (3, CLM.num_swb)

        # Litter / radiation outputs sized over nlevdecomp / soil.
        @test length(bc_out.litt_flux_cel_c_si) == nlevdecomp
        @test length(bc_out.qflx_soil2root_sisl) == nlevsoil

        # Bare-ground rootfr row (row 1) gets nominal 1/nlevsoil.
        @test all(≈(1.0 / nlevsoil), bc_out.rootfr_pa[1, 1:nlevsoil])

        # Now stamp non-zero values, then zero_bcs should clear them.
        fill!(bc_in.lightning24, 5.0)
        fill!(bc_in.smp_sl,      -123.0)
        fill!(bc_out.elai_pa,    9.0)
        fill!(bc_out.litt_flux_cel_c_si, 7.0)
        bc_out.gpp_site = 42.0

        CLM.zero_bcs(fates, 1)

        @test all(iszero, bc_in.lightning24)
        @test all(iszero, bc_in.smp_sl)
        @test all(iszero, bc_out.elai_pa)
        @test all(iszero, bc_out.litt_flux_cel_c_si)
        @test bc_out.gpp_site == 0.0
        # MIMICS litter quality is set to the unset sentinel, not zero.
        @test bc_out.litt_flux_ligc_per_n == CLM.fates_unset_r8
    end

    # =======================================================================
    # allocate_bcpconst / set_bcpconst — ECA param constants.
    # =======================================================================
    @testset "allocate + set bcpconst" begin
        numpft_test = 3
        nlevdecomp  = 2
        CLM.numpft[] = numpft_test

        # Populate a synthetic EDPftvarcon_inst with the ECA vectors set_bcpconst reads.
        pft = CLM.EDPftvarcon_type()
        pft.vmax_nh4         = [1.0, 2.0, 3.0]
        pft.vmax_no3         = [4.0, 5.0, 6.0]
        pft.vmax_p           = [7.0, 8.0, 9.0]
        pft.eca_km_nh4       = [0.1, 0.2, 0.3]
        pft.eca_km_no3       = [0.4, 0.5, 0.6]
        pft.eca_km_p         = [0.7, 0.8, 0.9]
        pft.eca_km_ptase     = [1.1, 1.2, 1.3]
        pft.eca_vmax_ptase   = [1.4, 1.5, 1.6]
        pft.eca_alpha_ptase  = [1.7, 1.8, 1.9]
        pft.eca_lambda_ptase = [2.1, 2.2, 2.3]
        CLM.EDPftvarcon_inst[] = pft

        # eca_plant_escalar lives on ed_params.
        edp = CLM.ed_params_type()
        edp.eca_plant_escalar = 0.5
        CLM.EDParams[] = edp

        bc_pconst = CLM.bc_pconst_type()
        CLM.allocate_bcpconst(bc_pconst, nlevdecomp)

        @test length(bc_pconst.vmax_nh4)      == numpft_test
        @test length(bc_pconst.eca_lambda_ptase) == numpft_test
        @test length(bc_pconst.j_uptake)      == nlevdecomp

        CLM.set_bcpconst(bc_pconst, nlevdecomp)
        @test bc_pconst.vmax_nh4         == pft.vmax_nh4
        @test bc_pconst.vmax_p           == pft.vmax_p
        @test bc_pconst.eca_km_ptase     == pft.eca_km_ptase
        @test bc_pconst.eca_lambda_ptase == pft.eca_lambda_ptase
        @test bc_pconst.eca_plant_escalar == 0.5
    end

    # =======================================================================
    # Global dimension counts: SetFatesGlobalElements2 (use_fates=false) +
    # InitPARTEHGlobals + fates_history_maps.
    # =======================================================================
    @testset "global dimension counts + history maps" begin
        # use_fates = false collapses the per-patch/per-site dims to 1.
        CLM.SetFatesGlobalElements2(false)
        @test CLM.fates_maxElementsPerPatch[] == 1
        @test CLM.fates_maxElementsPerSite[]  == 1

        # InitPARTEHGlobals (carbon path) sets the element list/positions.
        # Needs prt_params.leaf_long (numpft x nleafage) for the carbon allom init.
        numpft_test = 3
        nleafage    = 1
        CLM.numpft[]   = numpft_test
        CLM.nleafage[] = nleafage
        CLM.prt_params.leaf_long    = fill(1.5, numpft_test, nleafage)
        CLM.prt_params.wood_density = fill(0.6, numpft_test)

        CLM.hlm_parteh_mode[] = CLM.prt_carbon_allom_hyp
        CLM.InitPARTEHGlobals()
        @test CLM.num_elements[] == 1
        @test CLM.element_list[1] == CLM.carbon12_element
        @test CLM.element_pos[CLM.carbon12_element] == 1

        # fates_history_maps: set the history dim counts + the bin-edge value
        # arrays it reads, then verify the multiplexed maps are sized + filled.
        nlevsclass = 4
        nlevage    = 3
        nlevheight = 2
        nlevcoage  = 2
        nlevdamage = 2
        CLM.nlevsclass[] = nlevsclass
        CLM.nlevage[]    = nlevage
        CLM.nlevheight[] = nlevheight
        CLM.nlevcoage[]  = nlevcoage
        CLM.nlevdamage[] = nlevdamage

        edp = CLM.ed_params_type()
        edp.ED_val_history_sizeclass_bin_edges  = collect(range(0.0, length=nlevsclass, step=10.0))
        edp.ED_val_history_ageclass_bin_edges   = collect(range(0.0, length=nlevage,    step=5.0))
        edp.ED_val_history_height_bin_edges     = collect(range(0.0, length=nlevheight, step=2.0))
        edp.ED_val_history_coageclass_bin_edges = collect(range(0.0, length=nlevcoage,  step=3.0))
        edp.ED_val_history_damage_bin_edges     = collect(range(0.0, length=nlevdamage, step=0.5))
        edp.dlower_vai = collect(range(0.1, length=CLM.nlevleaf, step=0.1))
        CLM.EDParams[] = edp

        CLM.fates_history_maps()

        # Single dimensions.
        @test length(CLM.fates_hdim_levsclass[]) == nlevsclass
        @test length(CLM.fates_hdim_levpft[])    == numpft_test
        @test CLM.fates_hdim_levpft[] == collect(1:numpft_test)
        @test length(CLM.fates_hdim_levfuel[])   == CLM.num_fuel_classes

        # Multiplexed size-class x pft map.
        @test length(CLM.fates_hdim_pfmap_levscpf[]) == nlevsclass * numpft_test
        @test length(CLM.fates_hdim_scmap_levscpf[]) == nlevsclass * numpft_test
        # First nlevsclass entries are pft 1, sizes 1..nlevsclass.
        @test CLM.fates_hdim_pfmap_levscpf[][1] == 1
        @test CLM.fates_hdim_scmap_levscpf[][1:nlevsclass] == collect(1:nlevsclass)

        # Bin-edge value arrays were copied through.
        @test CLM.fates_hdim_levsclass[] == edp.ED_val_history_sizeclass_bin_edges
        @test CLM.fates_hdim_levleaf[]   == edp.dlower_vai

        # Element x pft map (num_elements == 1).
        @test length(CLM.fates_hdim_elmap_levelpft[]) == 1 * numpft_test
        @test all(==(1), CLM.fates_hdim_elmap_levelpft[])
    end

    # =======================================================================
    # DetermineGridCellNeighbors — neighbor list on a small grid.
    # =======================================================================
    @testset "DetermineGridCellNeighbors" begin
        numpft_test = 2
        CLM.numpft[] = numpft_test

        # Three gridcells: two very close, one far away. With a large max
        # dispersal distance, all become mutual neighbors; with a tiny one, none.
        gclat = [0.0, 0.0, 80.0]
        gclon = [0.0, 0.1, 0.0]

        # Use the exponential kernel; needs seed_dispersal_pdf_scale (per pft).
        pft = CLM.EDPftvarcon_type()
        pft.seed_dispersal_max_dist  = [2.0e7, 2.0e7]   # >Earth span -> everyone neighbors
        pft.seed_dispersal_pdf_scale = [1.0e-5, 1.0e-5]
        CLM.EDPftvarcon_inst[] = pft
        CLM.fates_dispersal_kernel_mode[] = CLM.fates_dispersal_kernel_exponential

        seeds = CLM.dispersal_type()
        neighbors = CLM.DetermineGridCellNeighbors(seeds, 3, gclat, gclon)

        @test length(neighbors) == 3
        # Each cell should neighbor the other two.
        @test neighbors[1].neighbor_count == 2
        @test neighbors[2].neighbor_count == 2
        @test neighbors[3].neighbor_count == 2
        @test sort(neighbors[1].neighbor_indices) == [2, 3]
        @test sort(neighbors[3].neighbor_indices) == [1, 2]
        # Per-PFT density-prob filled on the first neighbor of cell 1.
        @test length(neighbors[1].first_neighbor.density_prob) == numpft_test
        @test all(p -> 0.0 <= p <= 1.0, neighbors[1].first_neighbor.density_prob)
        # Single-process book-keeping.
        @test seeds.ncells_array == [3]
        @test seeds.begg_array   == [0]

        # Now with a tiny max distance: no neighbors at all.
        pft2 = CLM.EDPftvarcon_type()
        pft2.seed_dispersal_max_dist  = [1.0, 1.0]   # 1 m -> nobody neighbors
        pft2.seed_dispersal_pdf_scale = [1.0e-5, 1.0e-5]
        CLM.EDPftvarcon_inst[] = pft2

        neighbors2 = CLM.DetermineGridCellNeighbors(seeds, 3, gclat, gclon)
        @test all(n -> n.neighbor_count == 0, neighbors2)
        @test all(n -> isempty(n.neighbor_indices), neighbors2)
    end

end
