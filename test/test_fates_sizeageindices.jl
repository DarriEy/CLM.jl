# test_fates_sizeageindices.jl
# Tests for FATES FatesSizeAgeTypeIndicesMod (Tier F, Batch 2): the pure
# index-mapping helpers for the history/diagnostic multiplexed dimensions.
#
# Strategy: install synthetic monotonically-increasing bin edges on the global
# ed_params() instance and synthetic nlev* dimension counts on the
# FatesInterfaceTypesMod Refs, then assert
#   (a) the class-index lookups bin correctly (below-min -> 0, on a boundary,
#       between edges, at/above max),
#   (b) the flat (size,pft)/(age,pft)/(coage,pft)/(layer,size,pft) etc. index
#       combiners invert back to their factors (round-trip).

using Test
using CLM

@testset "FATES SizeAgeTypeIndices" begin

    # ---- synthetic bin edges (first edge is the minimum) --------------------
    # size: 0,10,20,30,40 -> 5 classes; age: 0,1,5,10 -> 4 classes;
    # height: 0,2,4 -> 3; coage: 0,3,6,9,12 -> 5.
    p = CLM.ed_params()
    size_edges   = [0.0, 10.0, 20.0, 30.0, 40.0]
    age_edges    = [0.0, 1.0, 5.0, 10.0]
    height_edges = [0.0, 2.0, 4.0]
    coage_edges  = [0.0, 3.0, 6.0, 9.0, 12.0]

    # Preserve & restore the global state so we don't disturb other testsets.
    saved = (p.ED_val_history_sizeclass_bin_edges,
             p.ED_val_history_ageclass_bin_edges,
             p.ED_val_history_height_bin_edges,
             p.ED_val_history_coageclass_bin_edges,
             CLM.nlevsclass[], CLM.nlevage[], CLM.nlevheight[],
             CLM.nlevcoage[], CLM.nlevdamage[])

    p.ED_val_history_sizeclass_bin_edges  = size_edges
    p.ED_val_history_ageclass_bin_edges   = age_edges
    p.ED_val_history_height_bin_edges     = height_edges
    p.ED_val_history_coageclass_bin_edges = coage_edges

    nsize   = length(size_edges)    # 5
    nage    = length(age_edges)     # 4
    nheight = length(height_edges)  # 3
    ncoage  = length(coage_edges)   # 5
    ndamage = 3
    CLM.nlevsclass[] = nsize
    CLM.nlevage[]    = nage
    CLM.nlevheight[] = nheight
    CLM.nlevcoage[]  = ncoage
    CLM.nlevdamage[] = ndamage

    try
        @testset "get_size_class_index binning" begin
            @test CLM.get_size_class_index(-1.0) == 0      # below the minimum
            @test CLM.get_size_class_index(0.0)  == 1      # exactly on first edge
            @test CLM.get_size_class_index(5.0)  == 1      # in [0,10)
            @test CLM.get_size_class_index(10.0) == 2      # on the second edge
            @test CLM.get_size_class_index(25.0) == 3      # in [20,30)
            @test CLM.get_size_class_index(40.0) == nsize  # on the last edge
            @test CLM.get_size_class_index(99.0) == nsize  # above the max
        end

        @testset "get_age_class_index binning" begin
            @test CLM.get_age_class_index(-0.5) == 0
            @test CLM.get_age_class_index(0.0)  == 1
            @test CLM.get_age_class_index(0.5)  == 1
            @test CLM.get_age_class_index(1.0)  == 2
            @test CLM.get_age_class_index(7.0)  == 3      # in [5,10)
            @test CLM.get_age_class_index(10.0) == nage   # on last edge
            @test CLM.get_age_class_index(50.0) == nage   # above max
        end

        @testset "get_height_index / get_coage_class_index binning" begin
            @test CLM.get_height_index(-1.0) == 0
            @test CLM.get_height_index(0.0)  == 1
            @test CLM.get_height_index(3.0)  == 2          # in [2,4)
            @test CLM.get_height_index(4.0)  == nheight
            @test CLM.get_height_index(99.0) == nheight

            @test CLM.get_coage_class_index(-1.0) == 0
            @test CLM.get_coage_class_index(0.0)  == 1
            @test CLM.get_coage_class_index(7.0)  == 3     # in [6,9)
            @test CLM.get_coage_class_index(12.0) == ncoage
            @test CLM.get_coage_class_index(99.0) == ncoage
        end

        @testset "(size, pft) flat index round-trip" begin
            # size_by_pft_class = (pft-1)*nlevsclass + size_class
            for pft in 1:6, dbh in (3.0, 15.0, 35.0, 99.0)
                size_class, idx = CLM.sizetype_class_index(dbh, pft)
                @test size_class == CLM.get_size_class_index(dbh)
                @test idx == (pft - 1) * nsize + size_class
                # invert (1-based): size_class in 1..nsize, pft from the quotient
                pft_back  = div(idx - 1, nsize) + 1
                size_back = idx - (pft_back - 1) * nsize
                @test pft_back == pft
                @test size_back == size_class
            end
        end

        @testset "(coage, pft) flat index round-trip" begin
            for pft in 1:5, coage in (1.0, 7.0, 99.0)
                coage_class, idx = CLM.coagetype_class_index(coage, pft)
                @test coage_class == CLM.get_coage_class_index(coage)
                @test idx == (pft - 1) * ncoage + coage_class
                pft_back   = div(idx - 1, ncoage) + 1
                coage_back = idx - (pft_back - 1) * ncoage
                @test pft_back == pft
                @test coage_back == coage_class
            end
        end

        @testset "(size, age) flat index round-trip" begin
            # size_by_age = (age_class-1)*nlevsclass + size_class
            for dbh in (3.0, 25.0, 99.0), age in (0.5, 7.0, 50.0)
                idx = CLM.get_sizeage_class_index(dbh, age)
                sc = CLM.get_size_class_index(dbh)
                ac = CLM.get_age_class_index(age)
                @test idx == (ac - 1) * nsize + sc
                age_back  = div(idx - 1, nsize) + 1
                size_back = idx - (age_back - 1) * nsize
                @test age_back == ac
                @test size_back == sc
            end
        end

        @testset "(age, pft) and (age, fuel) flat indices" begin
            # age_by_pft = age_class + (pft-1)*nlevage
            for pft in 1:4, age in (0.5, 7.0, 50.0)
                idx = CLM.get_agepft_class_index(age, pft)
                ac  = CLM.get_age_class_index(age)
                @test idx == ac + (pft - 1) * nage
                pft_back = div(idx - 1, nage) + 1
                age_back = idx - (pft_back - 1) * nage
                @test pft_back == pft
                @test age_back == ac
            end
            # get_agefuel_class_index has identical arithmetic with `fuel`.
            for fuel in 1:3, age in (0.5, 7.0)
                @test CLM.get_agefuel_class_index(age, fuel) ==
                      CLM.get_agepft_class_index(age, fuel)
            end
        end

        @testset "(size, age, pft) flat index round-trip" begin
            # idx = (age-1)*nsize + size + (pft-1)*nsize*nage
            for pft in 1:4, dbh in (3.0, 25.0), age in (0.5, 7.0)
                idx = CLM.get_sizeagepft_class_index(dbh, age, pft)
                sc  = CLM.get_size_class_index(dbh)
                ac  = CLM.get_age_class_index(age)
                @test idx == (ac - 1) * nsize + sc + (pft - 1) * nsize * nage
                stride = nsize * nage
                pft_back = div(idx - 1, stride) + 1
                rem      = idx - (pft_back - 1) * stride          # (ac-1)*nsize + sc
                age_back = div(rem - 1, nsize) + 1
                size_back = rem - (age_back - 1) * nsize
                @test pft_back == pft
                @test age_back == ac
                @test size_back == sc
            end
        end

        @testset "(layer, size, pft) flat index round-trip" begin
            # idx = (pft-1)*nlevsclass*nclmax + (size-1)*nclmax + layer
            nclmax = CLM.nclmax
            for pft in 1:4, dbh in (3.0, 25.0), layer in 1:nclmax
                idx = CLM.get_layersizetype_class_index(layer, dbh, pft)
                sc  = CLM.get_size_class_index(dbh)
                @test idx == (pft - 1) * nsize * nclmax + (sc - 1) * nclmax + layer
                # reverse mapping from the preserved Fortran comment (1-based)
                pft_back   = cld(idx, nsize * nclmax)
                size_back  = cld(idx - (pft_back - 1) * nsize * nclmax, nclmax)
                layer_back = idx - ((pft_back - 1) * nsize * nclmax + (size_back - 1) * nclmax)
                @test pft_back == pft
                @test size_back == sc
                @test layer_back == layer
            end
        end

        @testset "crown-damage flat indices" begin
            # cdamage_by_size = (cdamage-1)*nlevsclass + size
            for cdamage in 1:ndamage, dbh in (3.0, 25.0, 99.0)
                idx = CLM.get_cdamagesize_class_index(dbh, cdamage)
                sc  = CLM.get_size_class_index(dbh)
                @test idx == (cdamage - 1) * nsize + sc
            end
            # cdamage_by_size_by_pft = (cdamage-1)*nsize + size + (pft-1)*nsize*ndamage
            for pft in 1:3, cdamage in 1:ndamage, dbh in (3.0, 25.0)
                idx = CLM.get_cdamagesizepft_class_index(dbh, cdamage, pft)
                sc  = CLM.get_size_class_index(dbh)
                @test idx == (cdamage - 1) * nsize + sc + (pft - 1) * nsize * ndamage
            end
        end

    finally
        # Restore the global state.
        p.ED_val_history_sizeclass_bin_edges  = saved[1]
        p.ED_val_history_ageclass_bin_edges   = saved[2]
        p.ED_val_history_height_bin_edges     = saved[3]
        p.ED_val_history_coageclass_bin_edges = saved[4]
        CLM.nlevsclass[] = saved[5]
        CLM.nlevage[]    = saved[6]
        CLM.nlevheight[] = saved[7]
        CLM.nlevcoage[]  = saved[8]
        CLM.nlevdamage[] = saved[9]
    end
end
