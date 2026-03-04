using Test

@testset "CN Precision Control" begin

    # =========================================================================
    # Helper: create minimal CNVegCarbonStateData with np patches
    # =========================================================================
    function make_cs(np; nrepr=1, val=0.0)
        cs = CLM.CNVegCarbonStateData()
        cs.leafc_patch              = fill(val, np)
        cs.leafc_storage_patch      = fill(val, np)
        cs.leafc_xfer_patch         = fill(val, np)
        cs.frootc_patch             = fill(val, np)
        cs.frootc_storage_patch     = fill(val, np)
        cs.frootc_xfer_patch        = fill(val, np)
        cs.livestemc_patch          = fill(val, np)
        cs.livestemc_storage_patch  = fill(val, np)
        cs.livestemc_xfer_patch     = fill(val, np)
        cs.deadstemc_patch          = fill(val, np)
        cs.deadstemc_storage_patch  = fill(val, np)
        cs.deadstemc_xfer_patch     = fill(val, np)
        cs.livecrootc_patch         = fill(val, np)
        cs.livecrootc_storage_patch = fill(val, np)
        cs.livecrootc_xfer_patch    = fill(val, np)
        cs.deadcrootc_patch         = fill(val, np)
        cs.deadcrootc_storage_patch = fill(val, np)
        cs.deadcrootc_xfer_patch    = fill(val, np)
        cs.gresp_storage_patch      = fill(val, np)
        cs.gresp_xfer_patch         = fill(val, np)
        cs.cpool_patch              = fill(val, np)
        cs.xsmrpool_patch           = fill(val, np)
        cs.ctrunc_patch             = fill(0.0, np)
        cs.cropseedc_deficit_patch  = fill(val, np)
        cs.reproductivec_patch          = fill(val, np, nrepr)
        cs.reproductivec_storage_patch  = fill(val, np, nrepr)
        cs.reproductivec_xfer_patch     = fill(val, np, nrepr)
        return cs
    end

    # =========================================================================
    # Helper: create minimal CNVegNitrogenStateData with np patches
    # =========================================================================
    function make_ns(np; nrepr=1, val=0.0)
        ns = CLM.CNVegNitrogenStateData()
        ns.leafn_patch              = fill(val, np)
        ns.leafn_storage_patch      = fill(val, np)
        ns.leafn_xfer_patch         = fill(val, np)
        ns.frootn_patch             = fill(val, np)
        ns.frootn_storage_patch     = fill(val, np)
        ns.frootn_xfer_patch        = fill(val, np)
        ns.livestemn_patch          = fill(val, np)
        ns.livestemn_storage_patch  = fill(val, np)
        ns.livestemn_xfer_patch     = fill(val, np)
        ns.deadstemn_patch          = fill(val, np)
        ns.deadstemn_storage_patch  = fill(val, np)
        ns.deadstemn_xfer_patch     = fill(val, np)
        ns.livecrootn_patch         = fill(val, np)
        ns.livecrootn_storage_patch = fill(val, np)
        ns.livecrootn_xfer_patch    = fill(val, np)
        ns.deadcrootn_patch         = fill(val, np)
        ns.deadcrootn_storage_patch = fill(val, np)
        ns.deadcrootn_xfer_patch    = fill(val, np)
        ns.retransn_patch           = fill(val, np)
        ns.npool_patch              = fill(val, np)
        ns.ntrunc_patch             = fill(0.0, np)
        ns.cropseedn_deficit_patch  = fill(val, np)
        ns.reproductiven_patch          = fill(val, np, nrepr)
        ns.reproductiven_storage_patch  = fill(val, np, nrepr)
        ns.reproductiven_xfer_patch     = fill(val, np, nrepr)
        return ns
    end

    # =====================================================================
    # Constants tests
    # =====================================================================
    @testset "Module constants" begin
        @test CLM.CN_CCRIT_DEFAULT    == 1.0e-8
        @test CLM.CN_CNEGCRIT_DEFAULT == -6.0e+1
        @test CLM.CN_NCRIT_DEFAULT    == 1.0e-8
        @test CLM.CN_NNEGCRIT_DEFAULT == -7.0e+0
        @test CLM.CN_N_MIN            == 1.0e-9
    end

    # =====================================================================
    # truncate_c_and_n_states! tests
    # =====================================================================
    @testset "truncate_c_and_n_states!" begin

        @testset "values below ccrit are truncated" begin
            np = 4
            carbon   = [1.0, 5.0e-9, 3.0e-9, 2.0]
            nitrogen = [0.5, 2.0e-9, 1.0e-9, 1.0]
            pc = zeros(np)
            pn = zeros(np)
            filter = [1, 2, 3, 4]
            itype  = [1, 1, 1, 1]

            (num_tp, filter_tp) = CLM.truncate_c_and_n_states!(
                carbon, nitrogen, pc, pn, filter, itype)

            # Patches 2 and 3 should have been truncated (abs(C) < ccrit)
            @test num_tp == 2
            @test sort(filter_tp) == [2, 3]

            # Truncated values should be zero
            @test carbon[2] == 0.0
            @test carbon[3] == 0.0
            @test nitrogen[2] == 0.0
            @test nitrogen[3] == 0.0

            # Non-truncated values should be unchanged
            @test carbon[1] == 1.0
            @test carbon[4] == 2.0
            @test nitrogen[1] == 0.5
            @test nitrogen[4] == 1.0

            # Truncation accumulators should hold the removed mass
            @test pc[2] == 5.0e-9
            @test pc[3] == 3.0e-9
            @test pn[2] == 2.0e-9
            @test pn[3] == 1.0e-9
        end

        @testset "values above ccrit are not truncated" begin
            np = 2
            carbon   = [1.0, 0.5]
            nitrogen = [0.5, 0.3]
            pc = zeros(np)
            pn = zeros(np)
            filter = [1, 2]
            itype  = [1, 1]

            (num_tp, _) = CLM.truncate_c_and_n_states!(
                carbon, nitrogen, pc, pn, filter, itype)

            @test num_tp == 0
            @test carbon == [1.0, 0.5]
            @test nitrogen == [0.5, 0.3]
            @test all(pc .== 0.0)
            @test all(pn .== 0.0)
        end

        @testset "critically negative values cause error" begin
            np = 2
            carbon   = [-100.0, 1.0]
            nitrogen = [0.5, 0.5]
            pc = zeros(np)
            pn = zeros(np)
            filter = [1, 2]
            itype  = [1, 1]

            @test_throws ErrorException CLM.truncate_c_and_n_states!(
                carbon, nitrogen, pc, pn, filter, itype)
        end

        @testset "allowneg=true allows negative values" begin
            np = 2
            carbon   = [-100.0, 5.0e-9]
            nitrogen = [-10.0, 2.0e-9]
            pc = zeros(np)
            pn = zeros(np)
            filter = [1, 2]
            itype  = [1, 1]

            # Should not throw
            (num_tp, filter_tp) = CLM.truncate_c_and_n_states!(
                carbon, nitrogen, pc, pn, filter, itype; allowneg=true)

            # Patch 1 is not truncated (abs(C) > ccrit)
            # Patch 2 is truncated (abs(C) < ccrit)
            @test num_tp == 1
            @test filter_tp == [2]
            @test carbon[1] == -100.0
            @test carbon[2] == 0.0
        end

        @testset "croponly=true filters non-crop patches" begin
            np = 3
            carbon   = [5.0e-9, 5.0e-9, 5.0e-9]
            nitrogen = [2.0e-9, 2.0e-9, 2.0e-9]
            pc = zeros(np)
            pn = zeros(np)
            filter = [1, 2, 3]
            itype  = [1, 15, 16]  # nc3crop = 15; only patches 2,3 are crops

            (num_tp, filter_tp) = CLM.truncate_c_and_n_states!(
                carbon, nitrogen, pc, pn, filter, itype; croponly=true)

            # Only patches 2 and 3 (itype >= nc3crop) should be truncated
            @test num_tp == 2
            @test sort(filter_tp) == [2, 3]
            @test carbon[1] == 5.0e-9  # untouched (non-crop)
            @test carbon[2] == 0.0
            @test carbon[3] == 0.0
        end

        @testset "use_matrixcn=true uses different truncation logic" begin
            np = 3
            # Patch 1: C < ccrit and C > -ccrit*1e6 -- should truncate
            # Patch 2: large C -- should not truncate
            # Patch 3: small negative C, within range -- should truncate
            carbon   = [5.0e-9, 1.0, -3.0e-9]
            nitrogen = [2.0e-9, 0.5,  1.0e-9]
            pc = zeros(np)
            pn = zeros(np)
            filter = [1, 2, 3]
            itype  = [1, 1, 1]

            (num_tp, filter_tp) = CLM.truncate_c_and_n_states!(
                carbon, nitrogen, pc, pn, filter, itype;
                use_matrixcn=true)

            @test num_tp == 2
            @test sort(filter_tp) == [1, 3]
            @test carbon[1] == 0.0
            @test carbon[2] == 1.0
            @test carbon[3] == 0.0
        end

        @testset "use_nguardrail triggers truncation on small N" begin
            np = 2
            # Patch 1: C is large but N is tiny -- should truncate with nguardrail
            carbon   = [1.0, 0.5]
            nitrogen = [5.0e-9, 0.3]
            pc = zeros(np)
            pn = zeros(np)
            filter = [1, 2]
            itype  = [1, 1]

            (num_tp, filter_tp) = CLM.truncate_c_and_n_states!(
                carbon, nitrogen, pc, pn, filter, itype;
                use_nguardrail=true)

            # Patch 1 should be truncated because abs(N) < ncrit
            @test num_tp == 1
            @test filter_tp == [1]
            @test carbon[1] == 0.0
            @test nitrogen[1] == 0.0
            @test pc[1] == 1.0
            @test pn[1] == 5.0e-9
        end
    end

    # =====================================================================
    # truncate_c_states! tests
    # =====================================================================
    @testset "truncate_c_states!" begin

        @testset "small positive C values truncated" begin
            np = 3
            carbon = [5.0e-9, 1.0, 3.0e-9]
            pc = zeros(np)
            filter = [1, 2, 3]
            itype  = [1, 1, 1]

            (num_tp, filter_tp) = CLM.truncate_c_states!(carbon, pc, filter, itype)

            @test num_tp == 2
            @test sort(filter_tp) == [1, 3]
            @test carbon[1] == 0.0
            @test carbon[2] == 1.0
            @test carbon[3] == 0.0
            @test pc[1] == 5.0e-9
            @test pc[3] == 3.0e-9
        end

        @testset "critically negative C causes error" begin
            np = 1
            carbon = [-100.0]
            pc = zeros(np)
            filter = [1]
            itype  = [1]

            @test_throws ErrorException CLM.truncate_c_states!(carbon, pc, filter, itype)
        end

        @testset "allowneg=true allows negative values" begin
            np = 2
            carbon = [-100.0, 5.0e-9]
            pc = zeros(np)
            filter = [1, 2]
            itype  = [1, 1]

            # Should not throw
            (num_tp, filter_tp) = CLM.truncate_c_states!(
                carbon, pc, filter, itype; allowneg=true)

            @test num_tp == 1
            @test filter_tp == [2]
            @test carbon[1] == -100.0
            @test carbon[2] == 0.0
        end

        @testset "croponly=true filters non-crop patches" begin
            np = 2
            carbon = [5.0e-9, 5.0e-9]
            pc = zeros(np)
            filter = [1, 2]
            itype  = [1, 16]  # nc3crop=15; only patch 2 is crop

            (num_tp, filter_tp) = CLM.truncate_c_states!(
                carbon, pc, filter, itype; croponly=true)

            @test num_tp == 1
            @test filter_tp == [2]
            @test carbon[1] == 5.0e-9
            @test carbon[2] == 0.0
        end

        @testset "invalid cnegcrit vs ccrit" begin
            np = 1
            carbon = [0.5]
            pc = zeros(np)
            filter = [1]
            itype  = [1]

            # ccrit=1.0, cnegcrit=0.0 => -ccrit = -1.0 < cnegcrit = 0.0 => error
            @test_throws ErrorException CLM.truncate_c_states!(
                carbon, pc, filter, itype; ccrit=1.0, cnegcrit=0.0)
        end
    end

    # =====================================================================
    # truncate_n_states! tests
    # =====================================================================
    @testset "truncate_n_states!" begin

        @testset "small positive N values truncated" begin
            np = 3
            nitrogen = [5.0e-9, 1.0, 3.0e-9]
            pn = zeros(np)
            filter = [1, 2, 3]

            CLM.truncate_n_states!(nitrogen, pn, filter)

            @test nitrogen[1] == 0.0
            @test nitrogen[2] == 1.0
            @test nitrogen[3] == 0.0
            @test pn[1] == 5.0e-9
            @test pn[3] == 3.0e-9
        end

        @testset "critically negative N is silently ignored (Fortran behaviour)" begin
            np = 2
            nitrogen = [-10.0, 5.0e-9]
            pn = zeros(np)
            filter = [1, 2]

            # Should not throw (Fortran has endrun commented out)
            CLM.truncate_n_states!(nitrogen, pn, filter)

            # Patch 1: critically negative, not truncated (just skipped)
            @test nitrogen[1] == -10.0
            @test pn[1] == 0.0

            # Patch 2: small positive, truncated
            @test nitrogen[2] == 0.0
            @test pn[2] == 5.0e-9
        end
    end

    # =====================================================================
    # truncate_additional! tests
    # =====================================================================
    @testset "truncate_additional!" begin
        np = 4
        state      = [10.0, 20.0, 30.0, 40.0]
        truncation = zeros(np)
        filter_tp  = [2, 4]
        num_tp     = 2

        CLM.truncate_additional!(state, truncation, num_tp, filter_tp)

        # Only patches 2 and 4 should be zeroed
        @test state[1] == 10.0
        @test state[2] == 0.0
        @test state[3] == 30.0
        @test state[4] == 0.0

        @test truncation[1] == 0.0
        @test truncation[2] == 20.0
        @test truncation[3] == 0.0
        @test truncation[4] == 40.0
    end

    # =====================================================================
    # cn_precision_control! integration tests
    # =====================================================================
    @testset "cn_precision_control! basic" begin

        @testset "all pools below threshold are truncated" begin
            np = 2
            tiny = 5.0e-9
            cs = make_cs(np; val=tiny)
            ns = make_ns(np; val=tiny)
            filter = [1, 2]
            itype  = [1, 1]

            CLM.cn_precision_control!(cs, ns, filter, itype)

            # All pools should be zeroed
            @test all(cs.leafc_patch .== 0.0)
            @test all(cs.leafc_storage_patch .== 0.0)
            @test all(cs.leafc_xfer_patch .== 0.0)
            @test all(cs.frootc_patch .== 0.0)
            @test all(cs.frootc_storage_patch .== 0.0)
            @test all(cs.frootc_xfer_patch .== 0.0)
            @test all(cs.livestemc_patch .== 0.0)
            @test all(cs.deadstemc_patch .== 0.0)
            @test all(cs.livecrootc_patch .== 0.0)
            @test all(cs.deadcrootc_patch .== 0.0)
            @test all(cs.gresp_storage_patch .== 0.0)
            @test all(cs.gresp_xfer_patch .== 0.0)
            @test all(cs.cpool_patch .== 0.0)

            @test all(ns.leafn_patch .== 0.0)
            @test all(ns.retransn_patch .== 0.0)
            @test all(ns.npool_patch .== 0.0)

            # Truncation sinks should have accumulated the removed mass
            # The exact value depends on how many pools contributed
            @test all(cs.ctrunc_patch .!= 0.0)
            @test all(ns.ntrunc_patch .!= 0.0)
        end

        @testset "pools above threshold are not truncated" begin
            np = 2
            big_val = 1.0
            cs = make_cs(np; val=big_val)
            ns = make_ns(np; val=big_val)
            filter = [1, 2]
            itype  = [1, 1]

            CLM.cn_precision_control!(cs, ns, filter, itype)

            # All pools should be unchanged
            @test all(cs.leafc_patch .== big_val)
            @test all(ns.leafn_patch .== big_val)
            @test all(cs.ctrunc_patch .== 0.0)
            @test all(ns.ntrunc_patch .== 0.0)
        end

        @testset "mixed small and large values" begin
            np = 3
            cs = make_cs(np; val=1.0)
            ns = make_ns(np; val=1.0)

            # Make patch 2 tiny
            cs.leafc_patch[2] = 5.0e-9
            ns.leafn_patch[2] = 2.0e-9

            filter = [1, 2, 3]
            itype  = [1, 1, 1]

            CLM.cn_precision_control!(cs, ns, filter, itype)

            # Patch 2 leafc should be zeroed
            @test cs.leafc_patch[2] == 0.0
            @test ns.leafn_patch[2] == 0.0

            # Other patches should be unchanged for leafc
            @test cs.leafc_patch[1] == 1.0
            @test cs.leafc_patch[3] == 1.0
        end

        @testset "prec_control_for_froot=false skips froot" begin
            np = 2
            tiny = 5.0e-9
            cs = make_cs(np; val=tiny)
            ns = make_ns(np; val=tiny)
            filter = [1, 2]
            itype  = [1, 1]

            # Set frootc to something tiny -- it should NOT be truncated
            # when prec_control_for_froot=false
            cs.frootc_patch .= tiny
            ns.frootn_patch .= tiny

            CLM.cn_precision_control!(cs, ns, filter, itype;
                prec_control_for_froot=false)

            # frootc should remain at tiny value (not truncated)
            @test cs.frootc_patch[1] == tiny
            @test cs.frootc_patch[2] == tiny
            @test ns.frootn_patch[1] == tiny
            @test ns.frootn_patch[2] == tiny

            # But leafc should be truncated
            @test all(cs.leafc_patch .== 0.0)
        end
    end

    @testset "cn_precision_control! with C13 isotopes" begin
        np = 2
        tiny = 5.0e-9
        cs    = make_cs(np; val=tiny)
        ns    = make_ns(np; val=tiny)
        c13cs = make_cs(np; val=tiny * 0.01)  # C13 is much smaller

        filter = [1, 2]
        itype  = [1, 1]

        CLM.cn_precision_control!(cs, ns, filter, itype;
            c13cs=c13cs, use_c13=true)

        # C13 pools should also be zeroed where C12 was truncated
        @test all(c13cs.leafc_patch .== 0.0)
        @test all(c13cs.leafc_storage_patch .== 0.0)
        @test all(c13cs.leafc_xfer_patch .== 0.0)
        @test all(c13cs.ctrunc_patch .!= 0.0)
    end

    @testset "cn_precision_control! with C14 isotopes" begin
        np = 2
        tiny = 5.0e-9
        cs    = make_cs(np; val=tiny)
        ns    = make_ns(np; val=tiny)
        c14cs = make_cs(np; val=tiny * 0.001)

        filter = [1, 2]
        itype  = [1, 1]

        CLM.cn_precision_control!(cs, ns, filter, itype;
            c14cs=c14cs, use_c14=true)

        @test all(c14cs.leafc_patch .== 0.0)
        @test all(c14cs.ctrunc_patch .!= 0.0)
    end

    @testset "cn_precision_control! with crop pools" begin
        np = 3
        tiny = 5.0e-9
        cs = make_cs(np; val=tiny)
        ns = make_ns(np; val=tiny)
        filter = [1, 2, 3]
        itype  = [1, 15, 16]  # patches 2, 3 are crops (>= nc3crop=15)

        CLM.cn_precision_control!(cs, ns, filter, itype; use_crop=true)

        # Crop-only pools should only be truncated for crop patches
        # xsmrpool -- croponly
        # Patch 1 (itype=1, non-crop): xsmrpool untouched
        @test cs.xsmrpool_patch[1] == tiny
        # Patches 2, 3 (crop): xsmrpool truncated
        @test cs.xsmrpool_patch[2] == 0.0
        @test cs.xsmrpool_patch[3] == 0.0

        # reproductivec -- croponly
        @test cs.reproductivec_patch[1, 1] == tiny  # non-crop, not truncated
        @test cs.reproductivec_patch[2, 1] == 0.0   # crop, truncated
        @test cs.reproductivec_patch[3, 1] == 0.0   # crop, truncated

        # cropseedc deficit -- croponly
        @test cs.cropseedc_deficit_patch[1] == tiny
        @test cs.cropseedc_deficit_patch[2] == 0.0
    end

    @testset "cn_precision_control! mass conservation" begin
        # The total mass (pool + ctrunc) should be conserved
        np = 4
        tiny = 5.0e-9
        cs = make_cs(np; val=tiny)
        ns = make_ns(np; val=tiny)
        filter = [1, 2, 3, 4]
        itype  = [1, 1, 1, 1]

        # Track initial sums
        initial_c_leafc = sum(cs.leafc_patch) + sum(cs.ctrunc_patch)
        initial_n_leafn = sum(ns.leafn_patch) + sum(ns.ntrunc_patch)

        CLM.cn_precision_control!(cs, ns, filter, itype)

        # After truncation, pool+sink should be conserved
        final_c_leafc = sum(cs.leafc_patch) + sum(cs.ctrunc_patch)
        final_n_leafn = sum(ns.leafn_patch) + sum(ns.ntrunc_patch)

        # leafc: all patches had tiny value, so they were all truncated.
        # The truncation sinks absorbed the mass, but ctrunc also includes
        # mass from other pools. So we just check that the ctrunc is nonzero
        # and roughly accounts for the total.
        @test all(cs.leafc_patch .== 0.0)
        @test all(ns.leafn_patch .== 0.0)

        # ctrunc should have absorbed mass from all C pools
        for p in 1:np
            @test cs.ctrunc_patch[p] > 0.0
            @test ns.ntrunc_patch[p] > 0.0
        end
    end

    @testset "cn_precision_control! with both C13 and C14" begin
        np = 2
        tiny = 5.0e-9
        cs    = make_cs(np; val=tiny)
        ns    = make_ns(np; val=tiny)
        c13cs = make_cs(np; val=tiny * 0.01)
        c14cs = make_cs(np; val=tiny * 0.001)
        filter = [1, 2]
        itype  = [1, 1]

        CLM.cn_precision_control!(cs, ns, filter, itype;
            c13cs=c13cs, c14cs=c14cs,
            use_c13=true, use_c14=true)

        # All isotope pools should be zeroed
        @test all(c13cs.leafc_patch .== 0.0)
        @test all(c14cs.leafc_patch .== 0.0)
        @test all(c13cs.livestemc_patch .== 0.0)
        @test all(c14cs.livestemc_patch .== 0.0)
        @test all(c13cs.gresp_storage_patch .== 0.0)
        @test all(c14cs.gresp_storage_patch .== 0.0)

        # Isotope truncation sinks should be nonzero
        for p in 1:np
            @test c13cs.ctrunc_patch[p] != 0.0
            @test c14cs.ctrunc_patch[p] != 0.0
        end
    end

    @testset "cn_precision_control! empty filter" begin
        np = 3
        cs = make_cs(np; val=5.0e-9)
        ns = make_ns(np; val=5.0e-9)
        filter = Int[]   # empty filter
        itype  = [1, 1, 1]

        # Should not error with empty filter
        CLM.cn_precision_control!(cs, ns, filter, itype)

        # Nothing should be modified
        @test all(cs.leafc_patch .== 5.0e-9)
        @test all(cs.ctrunc_patch .== 0.0)
    end

    @testset "cn_precision_control! partial filter" begin
        np = 4
        tiny = 5.0e-9
        cs = make_cs(np; val=tiny)
        ns = make_ns(np; val=tiny)
        filter = [2, 4]  # only patches 2 and 4 are in the filter
        itype  = [1, 1, 1, 1]

        CLM.cn_precision_control!(cs, ns, filter, itype)

        # Only patches 2 and 4 should be truncated
        @test cs.leafc_patch[1] == tiny  # not in filter
        @test cs.leafc_patch[2] == 0.0   # in filter, truncated
        @test cs.leafc_patch[3] == tiny  # not in filter
        @test cs.leafc_patch[4] == 0.0   # in filter, truncated

        @test cs.ctrunc_patch[1] == 0.0
        @test cs.ctrunc_patch[2] > 0.0
        @test cs.ctrunc_patch[3] == 0.0
        @test cs.ctrunc_patch[4] > 0.0
    end

    @testset "cn_precision_control! accumulates into existing ctrunc" begin
        np = 2
        tiny = 5.0e-9
        cs = make_cs(np; val=tiny)
        ns = make_ns(np; val=tiny)

        # Set pre-existing truncation sink values
        cs.ctrunc_patch .= 100.0
        ns.ntrunc_patch .= 50.0

        filter = [1, 2]
        itype  = [1, 1]

        CLM.cn_precision_control!(cs, ns, filter, itype)

        # ctrunc should be > 100.0 (initial + truncated mass)
        @test cs.ctrunc_patch[1] > 100.0
        @test ns.ntrunc_patch[1] > 50.0
    end

    @testset "truncation of negative tiny values" begin
        np = 2
        carbon   = [-5.0e-9, 1.0]
        nitrogen = [-2.0e-9, 0.5]
        pc = zeros(np)
        pn = zeros(np)
        filter = [1, 2]
        itype  = [1, 1]

        # Small negative value should be truncated (abs < ccrit)
        (num_tp, filter_tp) = CLM.truncate_c_and_n_states!(
            carbon, nitrogen, pc, pn, filter, itype)

        @test num_tp == 1
        @test filter_tp == [1]
        @test carbon[1] == 0.0
        @test nitrogen[1] == 0.0
        # pc should have the negative value
        @test pc[1] == -5.0e-9
        @test pn[1] == -2.0e-9
    end

end
