# test/test_dyn_harvest.jl
#
# Tests for src/biogeochem/dyn_harvest.jl (port of dyn_subgrid/dynHarvestMod.F90).
#
# (1) File-reading test: write a temp NetCDF with a YEAR axis + the 5 HARVEST_*
#     variables, run dynHarvest_init + dynHarvest_interp! / resolve, and assert the
#     resolved per-PFT (here per-gridcell) harvest rate for a couple of years.
# (2) CNHarvest flux test: synthetic veg C/N pools + a known harvest fraction; run
#     cn_harvest!, assert C/N conservation (sum of pool losses == product-pool +
#     litter gains) and correct patch->column area-weighted aggregation.

using Test
using NCDatasets

@testset "dyn_harvest" begin

    # ---------------------------------------------------------------------
    # (1) File reading: init + interp over a couple of years
    # ---------------------------------------------------------------------
    @testset "file reading / interp" begin
        file_years = [2000, 2001, 2002]
        nspace = 2                       # 2 gridcells (begg=1, endg=2)
        ny = length(file_years)

        # Per-type, per-gridcell, per-year data. Index: [space, time].
        # Build distinct values so the per-type and the summed rate are checkable.
        # data[v][s, t]
        data = Dict{String,Matrix{Float64}}()
        for (vi, vn) in enumerate(CLM.harvest_varnames)
            m = Matrix{Float64}(undef, nspace, ny)
            for s in 1:nspace, t in 1:ny
                m[s, t] = vi + 0.1 * s + 0.01 * t   # arbitrary distinct values
            end
            data[vn] = m
        end

        mktempdir() do dir
            fn = joinpath(dir, "harvest.nc")
            NCDataset(fn, "c") do ds
                defDim(ds, "space", nspace)
                defDim(ds, "time", ny)
                defVar(ds, "YEAR", Int, ("time",))[:] = file_years
                for vn in CLM.harvest_varnames
                    v = defVar(ds, vn, Float64, ("space", "time"))
                    v[:, :] = data[vn]
                    v.attrib["units"] = CLM.harvest_unitless_units
                end
            end

            begg, endg = 1, 2
            state = CLM.dynHarvest_init(begg, endg, fn;
                                        current_year = 2001, use_fates = false)

            @test length(state.harvest_inst) == CLM.num_harvest_inst
            @test state.harvest_units == CLM.harvest_unitless_units
            @test length(state.harvest) == nspace

            # --- Year 2001 (interior; lower index = year 2001 = t=2) ---
            CLM.dynHarvest_interp!(state, begg, endg, 2001)
            @test state.do_harvest
            for s in 1:nspace
                expected = sum(data[vn][s, 2] for vn in CLM.harvest_varnames)
                @test state.harvest[s] ≈ expected
            end

            # resolve-types keeps the 5 distinct, for 2001 (t=2)
            rates, after = CLM.dynHarvest_interp_resolve_harvesttypes(state, begg, endg, 2001)
            @test after
            @test size(rates) == (nspace, CLM.num_harvest_inst)
            for s in 1:nspace, (vi, vn) in enumerate(CLM.harvest_varnames)
                @test rates[s, vi] ≈ data[vn][s, 2]
            end

            # --- Before the time series (1999): harvest off, rates zero ---
            CLM.dynHarvest_interp!(state, begg, endg, 1999)
            @test !state.do_harvest
            @test all(state.harvest .== 0.0)
            rates0, after0 = CLM.dynHarvest_interp_resolve_harvesttypes(state, begg, endg, 1999)
            @test !after0
            @test all(rates0 .== 0.0)

            # --- Past the end (2050): rates held at the LAST year (t=ny=3) ---
            CLM.dynHarvest_interp!(state, begg, endg, 2050)
            @test state.do_harvest
            for s in 1:nspace
                expected = sum(data[vn][s, ny] for vn in CLM.harvest_varnames)
                @test state.harvest[s] ≈ expected
            end
        end
    end

    # ---------------------------------------------------------------------
    # (2) CNHarvest flux test: conservation + patch->column aggregation
    # ---------------------------------------------------------------------
    @testset "cn_harvest! conservation + aggregation" begin
        # Two tree patches on the SAME column/gridcell (so we can check the
        # area-weighted patch->column aggregation), plus one non-tree patch that
        # must be left untouched.
        np = 3
        nc = 1
        ng = 1
        nlevdecomp = 1
        ndecomp_pools = 3
        i_litr_min, i_litr_max, i_met_lit = 1, 3, 1

        # --- PatchData ---
        patch = CLM.PatchData{Float64}()
        patch.itype    = [1, 7, 0]            # two tree PFTs (1,7), one bare (0)
        patch.column   = [1, 1, 1]
        patch.gridcell = [1, 1, 1]
        patch.wtcol    = [0.4, 0.5, 0.1]

        # --- pftcon: only lf_f / fr_f needed (1-based; ivt+1 indexing) ---
        # PFT rows are 1-based Julia for Fortran 0-based pft; size must cover ivt+1.
        npft_rows = 8
        pftcon = (lf_f = zeros(npft_rows, ndecomp_pools),
                  fr_f = zeros(npft_rows, ndecomp_pools))
        # Litter fractions must sum to 1 across litter pools (per the model);
        # give each used PFT row (ivt 1 -> row 2, ivt 7 -> row 8) a partition.
        for r in (2, 8)
            pftcon.lf_f[r, :] .= [0.5, 0.3, 0.2]
            pftcon.fr_f[r, :] .= [0.4, 0.4, 0.2]
        end

        # --- SoilBiogeochemStateData: unit profiles (1 level => prof = 1) ---
        sbs = CLM.SoilBiogeochemStateData{Float64}()
        sbs.leaf_prof_patch  = ones(np, nlevdecomp)
        sbs.froot_prof_patch = ones(np, nlevdecomp)
        sbs.croot_prof_patch = ones(np, nlevdecomp)
        sbs.stem_prof_patch  = ones(np, nlevdecomp)

        # --- CN veg state: give each patch distinct pool values ---
        cs = CLM.CNVegCarbonStateData{Float64}()
        ns = CLM.CNVegNitrogenStateData{Float64}()

        cpools_c = (:leafc_patch, :frootc_patch, :livestemc_patch, :deadstemc_patch,
                    :livecrootc_patch, :deadcrootc_patch, :xsmrpool_patch,
                    :leafc_storage_patch, :frootc_storage_patch, :livestemc_storage_patch,
                    :deadstemc_storage_patch, :livecrootc_storage_patch,
                    :deadcrootc_storage_patch, :gresp_storage_patch,
                    :leafc_xfer_patch, :frootc_xfer_patch, :livestemc_xfer_patch,
                    :deadstemc_xfer_patch, :livecrootc_xfer_patch, :deadcrootc_xfer_patch,
                    :gresp_xfer_patch)
        for (k, f) in enumerate(cpools_c)
            setfield!(cs, f, Float64[10.0 * k + p for p in 1:np])
        end

        npools_n = (:leafn_patch, :frootn_patch, :livestemn_patch, :deadstemn_patch,
                    :livecrootn_patch, :deadcrootn_patch, :retransn_patch,
                    :leafn_storage_patch, :frootn_storage_patch, :livestemn_storage_patch,
                    :deadstemn_storage_patch, :livecrootn_storage_patch,
                    :deadcrootn_storage_patch,
                    :leafn_xfer_patch, :frootn_xfer_patch, :livestemn_xfer_patch,
                    :deadstemn_xfer_patch, :livecrootn_xfer_patch, :deadcrootn_xfer_patch)
        for (k, f) in enumerate(npools_n)
            setfield!(ns, f, Float64[1.0 * k + 0.1 * p for p in 1:np])
        end

        # --- CN veg flux outputs: allocate the patch + col arrays we touch ---
        cf = CLM.CNVegCarbonFluxData{Float64}()
        nf = CLM.CNVegNitrogenFluxData{Float64}()

        cflux_patch = (:hrv_leafc_to_litter_patch, :hrv_frootc_to_litter_patch,
                       :hrv_livestemc_to_litter_patch, :wood_harvestc_patch,
                       :hrv_livecrootc_to_litter_patch, :hrv_deadcrootc_to_litter_patch,
                       :hrv_xsmrpool_to_atm_patch,
                       :hrv_leafc_storage_to_litter_patch, :hrv_frootc_storage_to_litter_patch,
                       :hrv_livestemc_storage_to_litter_patch, :hrv_deadstemc_storage_to_litter_patch,
                       :hrv_livecrootc_storage_to_litter_patch, :hrv_deadcrootc_storage_to_litter_patch,
                       :hrv_gresp_storage_to_litter_patch,
                       :hrv_leafc_xfer_to_litter_patch, :hrv_frootc_xfer_to_litter_patch,
                       :hrv_livestemc_xfer_to_litter_patch, :hrv_deadstemc_xfer_to_litter_patch,
                       :hrv_livecrootc_xfer_to_litter_patch, :hrv_deadcrootc_xfer_to_litter_patch,
                       :hrv_gresp_xfer_to_litter_patch)
        for f in cflux_patch
            setfield!(cf, f, zeros(np))
        end
        nflux_patch = (:hrv_leafn_to_litter_patch, :hrv_frootn_to_litter_patch,
                       :hrv_livestemn_to_litter_patch, :wood_harvestn_patch,
                       :hrv_livecrootn_to_litter_patch, :hrv_deadcrootn_to_litter_patch,
                       :hrv_retransn_to_litter_patch,
                       :hrv_leafn_storage_to_litter_patch, :hrv_frootn_storage_to_litter_patch,
                       :hrv_livestemn_storage_to_litter_patch, :hrv_deadstemn_storage_to_litter_patch,
                       :hrv_livecrootn_storage_to_litter_patch, :hrv_deadcrootn_storage_to_litter_patch,
                       :hrv_leafn_xfer_to_litter_patch, :hrv_frootn_xfer_to_litter_patch,
                       :hrv_livestemn_xfer_to_litter_patch, :hrv_deadstemn_xfer_to_litter_patch,
                       :hrv_livecrootn_xfer_to_litter_patch, :hrv_deadcrootn_xfer_to_litter_patch)
        for f in nflux_patch
            setfield!(nf, f, zeros(np))
        end

        # column accumulators (start at zero)
        cf.harvest_c_to_litr_c_col = zeros(nc, nlevdecomp, ndecomp_pools)
        cf.harvest_c_to_cwdc_col   = zeros(nc, nlevdecomp)
        cf.wood_harvestc_col       = zeros(nc)
        nf.harvest_n_to_litr_n_col = zeros(nc, nlevdecomp, ndecomp_pools)
        nf.harvest_n_to_cwdn_col   = zeros(nc, nlevdecomp)
        nf.wood_harvestn_col       = zeros(nc)

        # --- harvest state: unitless rate so m = harvest[g] / dt at year start ---
        # Choose a fraction am = 0.25 (1/yr) directly via unitless units.
        state = CLM.DynHarvestState(harvest = [0.25],
                                    do_harvest = true,
                                    harvest_units = CLM.harvest_unitless_units)

        dt = 1800.0
        mask = [true, true, true]

        CLM.cn_harvest!(state, mask, patch, pftcon, sbs, cs, ns, cf, nf;
                        dt = dt, is_beg_curr_year = true,
                        nlevdecomp = nlevdecomp,
                        i_litr_min = i_litr_min, i_litr_max = i_litr_max,
                        i_met_lit = i_met_lit)

        m = 0.25 / dt   # fractional rate per second, unitless path

        # --- non-tree patch (p=3) must be untouched ---
        @test cf.hrv_leafc_to_litter_patch[3] == 0.0
        @test nf.hrv_leafn_to_litter_patch[3] == 0.0
        @test cf.wood_harvestc_patch[3] == 0.0

        # --- a couple of explicit patch-level flux checks (tree p=1) ---
        @test cf.hrv_leafc_to_litter_patch[1] ≈ cs.leafc_patch[1] * m
        @test cf.wood_harvestc_patch[2]       ≈ cs.deadstemc_patch[2] * m
        @test nf.hrv_retransn_to_litter_patch[1] ≈ ns.retransn_patch[1] * m

        # --- Carbon conservation: total C removed from pools (patch, area-wtd)
        #     == column litter + CWD + product-pool gains. ---
        # Total C harvested out of every C pool, per patch, weighted by wtcol.
        total_c_out = 0.0
        for p in (1, 2)
            wt = patch.wtcol[p]
            for f in cpools_c
                total_c_out += getfield(cs, f)[p] * m * wt
            end
        end
        # xsmrpool goes to atmosphere (not to litter/cwd/product) — but it IS a
        # pool loss, so it must be accounted for as a separate sink.
        xsmr_to_atm = sum(cf.hrv_xsmrpool_to_atm_patch[p] * patch.wtcol[p] for p in (1, 2))

        col_c_gain = sum(cf.harvest_c_to_litr_c_col) + sum(cf.harvest_c_to_cwdc_col) +
                     sum(cf.wood_harvestc_col)
        @test col_c_gain + xsmr_to_atm ≈ total_c_out atol = 1e-9

        # --- Nitrogen conservation: all N pool losses route to litter/cwd/product
        #     (no atmosphere sink for N harvest). ---
        total_n_out = 0.0
        for p in (1, 2)
            wt = patch.wtcol[p]
            for f in npools_n
                total_n_out += getfield(ns, f)[p] * m * wt
            end
        end
        col_n_gain = sum(nf.harvest_n_to_litr_n_col) + sum(nf.harvest_n_to_cwdn_col) +
                     sum(nf.wood_harvestn_col)
        @test col_n_gain ≈ total_n_out atol = 1e-9

        # --- patch->column product-pool aggregation explicit check ---
        @test cf.wood_harvestc_col[1] ≈
            cf.wood_harvestc_patch[1] * patch.wtcol[1] +
            cf.wood_harvestc_patch[2] * patch.wtcol[2]
        @test nf.wood_harvestn_col[1] ≈
            nf.wood_harvestn_patch[1] * patch.wtcol[1] +
            nf.wood_harvestn_patch[2] * patch.wtcol[2]

        # --- with do_harvest=false, all fluxes stay zero ---
        state2 = CLM.DynHarvestState(harvest = [0.25], do_harvest = false,
                                     harvest_units = CLM.harvest_unitless_units)
        for f in cflux_patch; setfield!(cf, f, zeros(np)); end
        for f in nflux_patch; setfield!(nf, f, zeros(np)); end
        cf.harvest_c_to_litr_c_col = zeros(nc, nlevdecomp, ndecomp_pools)
        cf.harvest_c_to_cwdc_col   = zeros(nc, nlevdecomp)
        cf.wood_harvestc_col       = zeros(nc)
        CLM.cn_harvest!(state2, mask, patch, pftcon, sbs, cs, ns, cf, nf;
                        dt = dt, is_beg_curr_year = true, nlevdecomp = nlevdecomp,
                        i_litr_min = i_litr_min, i_litr_max = i_litr_max, i_met_lit = i_met_lit)
        @test all(cf.hrv_leafc_to_litter_patch[1:2] .== 0.0)
        @test all(cf.harvest_c_to_litr_c_col .== 0.0)
        @test all(cf.wood_harvestc_col .== 0.0)
    end

end
