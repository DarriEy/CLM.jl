# Tests for src/driver/dyn_subgrid_driver.jl — the top-level transient-land-use
# orchestrator (dynSubgrid_init! / dynSubgrid_driver! /
# dynSubgrid_wrapup_weight_changes!).
#
# We build a small synthetic single-gridcell domain (1 soil landunit, 1 column,
# 3 natural-PFT patches) plus a temp NetCDF flanduse_timeseries with a PCT_NAT_PFT
# time series across a few years, then:
#   - run dynSubgrid_init! (do_transient_pfts on),
#   - step dynSubgrid_driver! across 2001 / 2002 (changes) and a no-change re-step,
# asserting that:
#   * patch weights snap to the per-year file values,
#   * landunit weights still sum to 1 after each step,
#   * prior weights are snapshotted before the change,
#   * patch-level state is conserved across the area change (reusing the
#     PatchStateUpdater conservation accounting that the driver wires up), and
#   * a no-change year produces ~no weight adjustment.
# A small harvest-on smoke case is added at the end.

using Test
using NCDatasets

@testset "dyn_subgrid_driver" begin

    # ==================================================================
    # Core case: do_transient_pfts on, all others off.
    # ==================================================================
    @testset "init + driver across years (transient PFTs)" begin
        file_years = [2000, 2001, 2002]
        natpft_size = 3                       # 3 natural PFTs (itype 0,1,2)
        ng = 1                                # single grid cell

        # PCT for each (pft, year): columns (over pft) must sum to 100.
        pct = Array{Float64}(undef, natpft_size, ng, length(file_years))
        pct[:, 1, 1] = [50.0, 30.0, 20.0]     # year 2000 -> [.5 .3 .2]
        pct[:, 1, 2] = [20.0, 50.0, 30.0]     # year 2001 -> [.2 .5 .3]
        pct[:, 1, 3] = [10.0, 10.0, 80.0]     # year 2002 -> [.1 .1 .8]
        wexp = pct ./ 100.0

        mktempdir() do dir
            fn = joinpath(dir, "flanduse_pft.nc")
            NCDataset(fn, "c") do ds
                defDim(ds, "natpft", natpft_size)
                defDim(ds, "lndgrid", ng)
                defDim(ds, "time", length(file_years))
                defVar(ds, "YEAR", Int, ("time",))[:] = file_years
                pv = defVar(ds, "PCT_NAT_PFT", Float64, ("natpft", "lndgrid", "time"))
                for t in 1:length(file_years)
                    pv[:, :, t] = pct[:, :, t]      # (natpft, ng)
                end
            end

            np  = natpft_size
            ncol = 1
            bounds = CLM.BoundsType(begg = 1, endg = 1, begl = 1, endl = 1,
                                    begc = 1, endc = ncol, begp = 1, endp = np,
                                    level = CLM.BOUNDS_LEVEL_PROC)

            # --- subgrid: 1 gridcell, 1 soil landunit, 1 column, 3 PFT patches ---
            grc = CLM.GridcellData{Float64}()
            grc.landunit_indices = fill(CLM.ISPVAL, CLM.MAX_LUNIT, ng)
            grc.landunit_indices[CLM.ISTSOIL, 1] = 1

            lun = CLM.LandunitData{Float64}()
            lun.itype   = [CLM.ISTSOIL]
            lun.wtgcell = [1.0]
            lun.active  = falses(1)
            lun.coli    = [1]               # column index range for the landunit
            lun.colf    = [ncol]

            col = CLM.ColumnData{Float64}()
            col.landunit = [1]
            col.gridcell = [1]
            col.wtlunit  = [1.0]
            col.wtgcell  = [1.0]
            col.active   = falses(ncol)

            pch = CLM.PatchData{Float64}()
            pch.gridcell = fill(1, np)
            pch.landunit = fill(1, np)
            pch.column   = fill(1, np)
            pch.itype    = collect(0:(np - 1))    # natpft itype 0..np-1
            pch.wtcol    = fill(NaN, np)
            pch.wtlunit  = fill(NaN, np)
            pch.wtgcell  = fill(NaN, np)
            pch.active   = falses(np)

            ctl = CLM.dyn_subgrid_control_init(
                flanduse_timeseries = fn, do_transient_pfts = true)

            # --- init at the first file year ---
            state = CLM.dynSubgrid_init!(bounds, ctl, grc, lun, col, pch;
                current_year = 2000, natpft_size = natpft_size,
                check_dynpft_consistency = false)

            @test state isa CLM.DynSubgridState
            @test state.dynpft_state !== nothing
            @test state.dyncrop_state === nothing
            # init set the year-2000 weights and wrapped up higher-order weights.
            @test pch.wtcol ≈ wexp[:, 1, 1]
            @test pch.wtgcell ≈ wexp[:, 1, 1]      # single column at wtgcell 1
            @test isapprox(sum(CLM.get_landunit_weight(grc, lun, 1, lt)
                               for lt in 1:CLM.MAX_LUNIT), 1.0; atol = 1e-12)
            @test all(pch.active)                  # all 3 PFTs have positive weight

            # ------------------------------------------------------------
            # Step to 2001 (weights change). Track a per-area patch state for
            # the conservation check using the driver's PatchStateUpdater.
            # ------------------------------------------------------------
            var0   = [10.0, 20.0, 5.0]             # per-area state before the change
            pwt_old = copy(pch.wtgcell)

            CLM.dynSubgrid_driver!(state, bounds, grc, lun, col, pch; year = 2001)

            # prior weights snapshotted BEFORE the change == pre-step weights.
            @test state.prior_weights.pwtgcell[1:np] ≈ pwt_old
            # patch weights snapped to year-2001 values.
            @test pch.wtcol ≈ wexp[:, 1, 2]
            @test pch.wtgcell ≈ wexp[:, 1, 2]
            @test isapprox(sum(CLM.get_landunit_weight(grc, lun, 1, lt)
                               for lt in 1:CLM.MAX_LUNIT), 1.0; atol = 1e-12)

            # The updater (set up inside the driver) carries old/new patch weights.
            @test state.patch_state_updater.pwtgcell_old[1:np] ≈ pwt_old
            @test state.patch_state_updater.pwtgcell_new[1:np] ≈ pch.wtgcell

            # --- conservation: update a per-area state and check mass balance ---
            var = copy(var0)
            flux_grc = zeros(Float64, np)
            filterp  = collect(1:np)
            CLM.update_patch_state!(state.patch_state_updater, bounds, pch, filterp,
                                    var; flux_out_grc_area = flux_grc)
            mass_before = sum(var0 .* pwt_old)
            mass_after  = sum(var  .* pch.wtgcell)
            total_flux  = sum(flux_grc)            # negative: mass leaving shrinking patches
            @test isapprox(mass_after - total_flux, mass_before; atol = 1e-12)

            # ------------------------------------------------------------
            # Step to 2002 (after the time series -> last year's weights).
            # ------------------------------------------------------------
            CLM.dynSubgrid_driver!(state, bounds, grc, lun, col, pch; year = 2002)
            @test pch.wtcol ≈ wexp[:, 1, 3]
            @test isapprox(sum(CLM.get_landunit_weight(grc, lun, 1, lt)
                               for lt in 1:CLM.MAX_LUNIT), 1.0; atol = 1e-12)

            # ------------------------------------------------------------
            # No-change re-step at 2002 -> no weight adjustment, updater flags none.
            # ------------------------------------------------------------
            wt_before = copy(pch.wtgcell)
            CLM.dynSubgrid_driver!(state, bounds, grc, lun, col, pch; year = 2002,
                                   clump_index = 1)
            @test pch.wtgcell ≈ wt_before
            @test all(state.patch_state_updater.dwt[1:np] .≈ 0.0)
            @test state.column_state_updater.any_changes[1] == false

            CLM.dynpft_close!(state.dynpft_state)
        end
    end

    # ==================================================================
    # Harvest-on smoke: do_harvest on (use_cn so the consistency check passes),
    # do_transient_pfts off. Confirms the driver wires the harvest reader.
    # ==================================================================
    @testset "driver with do_harvest" begin
        file_years = [2000, 2001, 2002]
        ng = 1
        mktempdir() do dir
            fn = joinpath(dir, "flanduse_harvest.nc")
            NCDataset(fn, "c") do ds
                defDim(ds, "lndgrid", ng)
                defDim(ds, "time", length(file_years))
                defVar(ds, "YEAR", Int, ("time",))[:] = file_years
                for (vi, vn) in enumerate(CLM.harvest_varnames)
                    v = defVar(ds, vn, Float64, ("lndgrid", "time"))
                    for t in 1:length(file_years)
                        v[:, t] = [Float64(vi) + 0.01 * t]
                    end
                    v.attrib["units"] = CLM.harvest_unitless_units
                end
            end

            bounds = CLM.BoundsType(begg = 1, endg = 1, begl = 1, endl = 1,
                                    begc = 1, endc = 1, begp = 1, endp = 1,
                                    level = CLM.BOUNDS_LEVEL_PROC)

            grc = CLM.GridcellData{Float64}()
            grc.landunit_indices = fill(CLM.ISPVAL, CLM.MAX_LUNIT, ng)
            grc.landunit_indices[CLM.ISTSOIL, 1] = 1
            lun = CLM.LandunitData{Float64}()
            lun.itype = [CLM.ISTSOIL]; lun.wtgcell = [1.0]; lun.active = falses(1)
            lun.coli = [1]; lun.colf = [1]
            col = CLM.ColumnData{Float64}()
            col.landunit = [1]; col.gridcell = [1]
            col.wtlunit = [1.0]; col.wtgcell = [1.0]; col.active = falses(1)
            pch = CLM.PatchData{Float64}()
            pch.gridcell = [1]; pch.landunit = [1]; pch.column = [1]
            pch.itype = [0]; pch.wtcol = [1.0]
            pch.wtlunit = [1.0]; pch.wtgcell = [1.0]; pch.active = falses(1)

            ctl = CLM.dyn_subgrid_control_init(
                flanduse_timeseries = fn, do_harvest = true, use_cn = true)

            state = CLM.dynSubgrid_init!(bounds, ctl, grc, lun, col, pch;
                current_year = 2000)
            @test state.dynharvest_state !== nothing
            # No interp during init (matching Fortran) -> harvest not yet active.

            # Driver step reads the harvest rates for the year.
            CLM.dynSubgrid_driver!(state, bounds, grc, lun, col, pch; year = 2001)
            @test state.dynharvest_state.do_harvest == true
            @test length(state.dynharvest_state.harvest) == 1
            # Landunit weights still sum to 1.
            @test isapprox(sum(CLM.get_landunit_weight(grc, lun, 1, lt)
                               for lt in 1:CLM.MAX_LUNIT), 1.0; atol = 1e-12)

            CLM.dyn_file_close!(state.dynharvest_state.dynHarvest_file)
        end
    end

end
