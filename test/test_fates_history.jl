# test_fates_history.jl
# Tests for the FATES history I/O pair (Tier F, Batch 18):
#   * FatesHistoryVariableType — one history output variable (metadata + buffer)
#   * FatesHistoryInterfaceMod — the registry/dimension bookkeeping + define/update
#
# Verifies: var-type Init!/HFlush!/GetBounds; interface Init! dimension bookkeeping;
# initialize_history_vars! building a non-empty, correctly-counted registry; a
# sample of registered vars carrying the expected name/units/dimension-kind; and
# one update_history_dyn1! path aggregating a constructed site/patch/cohort into
# the right output slot.

using Test
using CLM

@testset "FATES history I/O (Batch 18)" begin

    # ------------------------------------------------------------------
    # 1. fates_history_variable_type: Init! / GetBounds / HFlush!
    # ------------------------------------------------------------------
    @testset "history variable type" begin
        # Build a minimal dim_kinds + dim_bounds machinery for a 1D (site) and a
        # 2D (site x pft) variable.
        dim_bounds = [CLM.fates_io_dimension_type() for _ in 1:CLM.fates_history_num_dimensions]
        dim_kinds  = [CLM.fates_io_variable_kind_type() for _ in 1:CLM.fates_history_num_dim_kinds]

        # column dim -> slot 1 (sites 1:3); levpft dim -> slot 2 (pfts 1:4)
        CLM.Init!(dim_bounds[1], CLM.dimname_column, 1, 1, 3)
        CLM.Init!(dim_bounds[2], CLM.dimname_levpft,  1, 1, 4)

        # site_r8 kind in slot 1 (1D) bound to dim 1
        CLM.Init!(dim_kinds[1], CLM.site_r8, 1)
        dim_kinds[1].dim1_index = 1
        # site_pft_r8 kind in slot 2 (2D) bound to dims 1,2
        CLM.Init!(dim_kinds[2], CLM.site_pft_r8, 2)
        dim_kinds[2].dim1_index = 1
        dim_kinds[2].dim2_index = 2

        # --- 1D variable ---
        hv1 = CLM.fates_history_variable_type()
        CLM.Init!(hv1, "FATES_TEST_SI", "kg", "test site var", "active",
                  CLM.site_r8, "A", -9999.0, CLM.group_dyna_simple, 2, dim_kinds, dim_bounds)
        @test hv1.vname == "FATES_TEST_SI"
        @test hv1.units == "kg"
        @test hv1.allocated == :r81d
        @test length(hv1.r81d) == 3
        @test all(hv1.r81d .== -9999.0)        # flushed to flushval on Init
        @test hv1.dim_kinds_index == 1
        @test CLM.is_active(dim_kinds[1])      # Init! marks the kind active

        # GetBounds returns the whole-proc bounds for thread 0
        (lb1, ub1, lb2, ub2) = CLM.GetBounds(hv1, 0, dim_bounds, dim_kinds)
        @test (lb1, ub1, lb2, ub2) == (1, 3, 0, 0)

        # write then flush
        hv1.r81d .= 5.0
        CLM.HFlush!(hv1, 0, dim_bounds, dim_kinds)
        @test all(hv1.r81d .== -9999.0)

        # --- 2D variable ---
        hv2 = CLM.fates_history_variable_type()
        CLM.Init!(hv2, "FATES_TEST_PF", "1", "test pft var", "active",
                  CLM.site_pft_r8, "A", 0.0, CLM.group_dyna_complx, 2, dim_kinds, dim_bounds)
        @test hv2.allocated == :r82d
        @test size(hv2.r82d) == (3, 4)
        @test all(hv2.r82d .== 0.0)
        (lb1, ub1, lb2, ub2) = CLM.GetBounds(hv2, 0, dim_bounds, dim_kinds)
        @test (lb1, ub1, lb2, ub2) == (1, 3, 1, 4)
    end

    # ------------------------------------------------------------------
    # 2. fates_history_interface_type: Init! dimension bookkeeping
    # ------------------------------------------------------------------
    @testset "interface dimension bookkeeping" begin
        # Build a synthetic bounds object covering the dims our sample uses.
        nsite = 2; npft = 5; nscls = 4; nscpf = npft * nscls
        fb = CLM.fates_bounds_type(
            column_begin = 1,            column_end = nsite,
            pft_class_begin = 1,         pft_class_end = npft,
            size_class_begin = 1,        size_class_end = nscls,
            sizepft_class_begin = 1,     sizepft_class_end = nscpf,
            soil_begin = 1,              soil_end = 3,
            age_class_begin = 1,         age_class_end = 2,
        )
        hist = CLM.fates_history_interface_type()
        CLM.Init!(hist, 1, fb)

        # 29 dimensions get index handles (1-based, in _DIM_SETUP order)
        @test CLM.column_index(hist) == 1
        @test CLM.levsoil_index(hist) == 2
        @test CLM.levpft_index(hist) == 7
        @test hist.dim_bounds[CLM.column_index(hist)].upper_bound == nsite
        @test hist.dim_bounds[CLM.levpft_index(hist)].upper_bound == npft

        # assemble dim-kind maps + indices
        CLM.assemble_history_output_types!(hist)
        ityp = CLM.iotype_index(CLM.site_pft_r8, CLM.fates_history_num_dim_kinds, hist.dim_kinds)
        dk = hist.dim_kinds[ityp]
        @test dk.ndims == 2
        @test dk.dim1_index == CLM.column_index(hist)
        @test dk.dim2_index == CLM.levpft_index(hist)
        @test dk.dimsize[2] == npft       # set from bounds upper-lower+1
        @test dk.dimsize[1] == nsite
    end

    # ------------------------------------------------------------------
    # 3. define_history_vars! / initialize_history_vars!: registry build
    # ------------------------------------------------------------------
    @testset "registry build + sample variables" begin
        # Save + set the HLM control flags so the full dynamics+hifrq registry
        # (minus hydro/treedamage/nocomp/element gating) is built.
        saved = (CLM.hlm_hist_level_dynam[], CLM.hlm_hist_level_hifrq[],
                 CLM.hlm_use_planthydro[], CLM.hlm_use_tree_damage[],
                 CLM.hlm_use_nocomp[], CLM.hlm_use_ed_st3[], CLM.hlm_use_sp[],
                 CLM.hlm_hio_ignore_val[], CLM.hlm_name[])
        try
            CLM.hlm_hist_level_dynam[] = 2   # dynam0 + dynam1 active
            CLM.hlm_hist_level_hifrq[] = 2   # hifrq0 + hifrq1 active
            CLM.hlm_use_planthydro[]   = CLM.ifalse
            CLM.hlm_use_tree_damage[]  = CLM.ifalse
            CLM.hlm_use_nocomp[]       = CLM.ifalse
            CLM.hlm_use_ed_st3[]       = CLM.ifalse
            CLM.hlm_use_sp[]           = CLM.ifalse
            CLM.hlm_hio_ignore_val[]   = -9999.0
            CLM.hlm_name[]             = ""   # build all vars listing CLM:ALM

            # element_list is empty in this standalone build -> has_n/has_p gates
            # are FALSE, so N/P vars (25 + 21) are excluded. Tree-damage (18),
            # nocomp (3), planthydro (3+7+26), and not-st3-not-sp (3) vars: tree
            # damage/nocomp/planthydro OFF here.
            nscls = 4; npft = 5; nscpf = npft * nscls
            fb = CLM.fates_bounds_type(
                column_begin = 1,            column_end = 2,
                pft_class_begin = 1,         pft_class_end = npft,
                size_class_begin = 1,        size_class_end = nscls,
                sizepft_class_begin = 1,     sizepft_class_end = nscpf,
                soil_begin = 1,              soil_end = 3,
                age_class_begin = 1,         age_class_end = 2,
                fuel_begin = 1,              fuel_end = 4,
                cwdsc_begin = 1,             cwdsc_end = 4,
                can_begin = 1,               can_end = 2,
                cnlf_begin = 1,              cnlf_end = 6,
                cnlfpft_begin = 1,           cnlfpft_end = 30,
                cdpf_begin = 1,              cdpf_end = nscpf,
                cdsc_begin = 1,              cdsc_end = nscls,
                cdam_begin = 1,              cdam_end = 3,
                sizeage_class_begin = 1,     sizeage_class_end = 8,
                sizeagepft_class_begin = 1,  sizeagepft_class_end = 40,
                agepft_class_begin = 1,      agepft_class_end = 10,
                coage_class_begin = 1,       coage_class_end = 2,
                coagepf_class_begin = 1,     coagepf_class_end = 10,
                height_begin = 1,            height_end = 4,
                elem_begin = 1,              elem_end = 1,
                elpft_begin = 1,             elpft_end = 5,
                elcwd_begin = 1,             elcwd_end = 4,
                elage_begin = 1,             elage_end = 2,
                agefuel_begin = 1,           agefuel_end = 8,
                clscpf_begin = 1,            clscpf_end = nscpf,
                landuse_begin = 1,           landuse_end = 5,
                lulu_begin = 1,              lulu_end = 25,
                lupft_begin = 1,             lupft_end = 25,
            )
            hist = CLM.fates_history_interface_type()
            CLM.Init!(hist, 1, fb)
            CLM.assemble_history_output_types!(hist)
            CLM.initialize_history_vars!(hist)

            @test hist.num_history_vars_ > 0
            @test length(hist.hvars) == hist.num_history_vars_

            # Count expected = all rows whose gates pass under these flags.
            expected = count(r -> CLM._all_gates_pass(r[1]), CLM._HISTORY_VAR_REGISTRY)
            @test hist.num_history_vars_ == expected
            # Sanity: should be a few hundred (dynam + hifrq, no N/P/hydro/damage)
            @test hist.num_history_vars_ > 250

            # Each active var got an Init! buffer
            @test all(hv -> hv.allocated != :none, hist.hvars)

            # The unconditional vars should all be active with a nonzero handle.
            @test hist.ih[:ih_npatches_si] > 0
            @test hist.ih[:ih_ncohorts_si] > 0
            @test hist.ih[:ih_gpp_si]      > 0   # hifrq var

            # The 3 declared-but-never-registered handles stay 0.
            @test hist.ih[:ih_l2fr_clscpf]   == 0
            @test hist.ih[:ih_supsub_scpf]   == 0
            @test hist.ih[:ih_crownarea_clll] == 0

            # Sample: FATES_NPATCHES metadata matches the registry/Fortran.
            hv = hist.hvars[hist.ih[:ih_npatches_si]]
            @test hv.vname == "FATES_NPATCHES"
            @test hv.units == ""
            @test hv.vtype == CLM.site_r8
            @test hv.upfreq == CLM.group_dyna_simple
            @test hv.allocated == :r81d

            # Sample: a 2D size-pft var (FATES_TRIMMING is site_r8; pick a scpf one).
            # FATES_VEGC_SZPF -> ih_totvegc_scpf (site_size_pft_r8). Look it up by name.
            idx_scpf = findfirst(hv -> hv.vtype == CLM.site_size_pft_r8, hist.hvars)
            @test idx_scpf !== nothing
            @test hist.hvars[idx_scpf].allocated == :r82d

            # FATES_FRACTION flushes to zero (flush_to_zero=true), not ignore_val.
            hvfrac = hist.hvars[hist.ih[:ih_fates_fraction_si]]
            @test hvfrac.flushval == 0.0

            # ------------------------------------------------------------------
            # 4. update_history_dyn1!: aggregate a site/patch/cohort
            # ------------------------------------------------------------------
            @testset "update_history_dyn1! aggregation" begin
                # Two patches: one primary, one secondary; the primary holds 2
                # cohorts, the secondary holds 1. site h_gid = 1.
                co1 = CLM.fates_cohort_type(n = 1.0, dbh = 10.0)
                co2 = CLM.fates_cohort_type(n = 2.0, dbh = 20.0)
                co3 = CLM.fates_cohort_type(n = 3.0, dbh = 5.0)
                # primary patch cohort list (doubly linked, as the real cold-start
                # builds): co1 (tallest) <-> co2 (shortest). update_history_dyn1!
                # walks shortest -> taller, matching Fortran.
                co1.shorter = co2
                co2.taller  = co1

                p_primary = CLM.fates_patch_type(
                    land_use_label = 1, area = 6000.0,
                    total_canopy_area = 6000.0, total_tree_area = 4000.0,
                    tallest = co1, shortest = co2)
                p_secondary = CLM.fates_patch_type(
                    land_use_label = CLM.secondaryland, area = 4000.0,
                    total_canopy_area = 4000.0, total_tree_area = 1000.0,
                    tallest = co3, shortest = co3)
                # age order: oldest -> younger. Make primary oldest.
                p_primary.younger = p_secondary

                site = CLM.ed_site_type(
                    h_gid = 1, oldest_patch = p_primary,
                    cstatus = 2, nchilldays = 7, ncolddays = 3,
                    grow_deg_days = 123.5,
                    phen_model_date = 100, cleafoffdate = 60, cleafondate = 90)

                CLM.update_history_dyn1!(hist, 1, 1, [site], CLM.bc_in_type[])

                io = site.h_gid
                @test hist.hvars[hist.ih[:ih_npatches_si]].r81d[io]     == 2.0
                @test hist.hvars[hist.ih[:ih_npatches_sec_si]].r81d[io] == 1.0
                @test hist.hvars[hist.ih[:ih_ncohorts_si]].r81d[io]     == 3.0
                @test hist.hvars[hist.ih[:ih_ncohorts_sec_si]].r81d[io] == 1.0

                # area_plant = (6000 + 4000)/AREA(=10000) = 1.0
                @test hist.hvars[hist.ih[:ih_area_plant_si]].r81d[io] ≈ 1.0
                # area_trees = (4000 + 1000)/10000 = 0.5
                @test hist.hvars[hist.ih[:ih_area_trees_si]].r81d[io] ≈ 0.5

                # site-level scalar state copies
                @test hist.hvars[hist.ih[:ih_site_cstatus_si]].r81d[io]    == 2.0
                @test hist.hvars[hist.ih[:ih_gdd_si]].r81d[io]             == 123.5
                @test hist.hvars[hist.ih[:ih_site_nchilldays_si]].r81d[io] == 7.0
                @test hist.hvars[hist.ih[:ih_site_ncolddays_si]].r81d[io]  == 3.0
                @test hist.hvars[hist.ih[:ih_cleafoff_si]].r81d[io]        == 40.0  # 100-60
                @test hist.hvars[hist.ih[:ih_cleafon_si]].r81d[io]         == 10.0  # 100-90
            end

            # ------------------------------------------------------------------
            # 5. flush_hvars! restores the flush value
            # ------------------------------------------------------------------
            @testset "flush_hvars!" begin
                hv = hist.hvars[hist.ih[:ih_npatches_si]]
                @test hv.r81d[1] != hv.flushval
                CLM.flush_hvars!(hist, 1, CLM.group_dyna_simple)
                @test hv.r81d[1] == hv.flushval
            end

        finally
            (CLM.hlm_hist_level_dynam[], CLM.hlm_hist_level_hifrq[],
             CLM.hlm_use_planthydro[], CLM.hlm_use_tree_damage[],
             CLM.hlm_use_nocomp[], CLM.hlm_use_ed_st3[], CLM.hlm_use_sp[],
             CLM.hlm_hio_ignore_val[], CLM.hlm_name[]) = saved
        end
    end

end
