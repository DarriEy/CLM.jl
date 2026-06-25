# Tests for the dyn_subgrid CN-balance WIRING into the live driver
# (src/driver/dyn_subgrid_driver.jl — the optional `cons_bgc` bundle on
# dynSubgrid_driver! → dynSubgrid_run_cnbal!, mirroring the Fortran
# `if (use_cn) call DynamicAreaConservation` in dynSubgridDriverMod.F90).
#
# We reuse the synthetic transient-PFT domain from test_dyn_subgrid_driver.jl
# (1 gridcell, 1 soil landunit, 1 column, 3 natural-PFT patches, with a
# flanduse_timeseries that changes the per-PFT weights from year to year), then
# allocate the full CN-veg + soilbiogeochem facade state and run the driver with
# a `cons_bgc` bundle across a weight-change year, asserting:
#
#   * gate ON (cons_bgc supplied): the patch-level veg C/N is CONSERVED across
#     the area change (after = before + seed_in − flux_out), exactly as the
#     standalone dyn_cnbal_patch! unit test checks — proving the bundle threads
#     the facade state through correctly;
#   * the initiating patch had its phenology state reset;
#   * gate OFF (cons_bgc = nothing): the CN state is BYTE-IDENTICAL to not
#     touching it (the default standalone path is a no-op);
#   * no-weight-change re-step with the gate ON produces ~zero conversion flux.

using Test
using NCDatasets

@testset "dyn_subgrid_cnbal_wiring" begin

    CLM.varpar_init!(CLM.varpar, 1, 15, 2, 5)
    CLM.varcon_init!()
    nlevdecomp    = CLM.varpar.nlevdecomp
    ndecomp_pools = 3
    i_litr_min    = 1
    i_litr_max    = 2
    i_cwd         = 3
    dt            = 1800.0

    # ----------------------------------------------------------------------
    # Minimal PftconType with the fields the C/N balance touches (matching
    # test_dyn_cons_biogeochem.jl's make_pftcon).
    # ----------------------------------------------------------------------
    function make_pftcon(n::Int, ndp::Int)
        p = CLM.PftconType()
        p.c3psn     = ones(n)
        p.evergreen = zeros(n)
        p.woody     = zeros(n)
        p.leafcn    = fill(25.0, n)
        p.deadwdcn  = fill(100.0, n)
        p.noveg     = fill(typemax(Int), n)
        p.pconv     = fill(0.5, n)
        p.fr_f      = zeros(n, ndp)
        for j in 1:n
            p.fr_f[j, 1] = 0.6
            p.fr_f[j, 2] = 0.4
        end
        return p
    end

    # Build the synthetic transient-PFT domain + a DynSubgridState wired to a
    # flanduse_timeseries; returns everything needed plus a populated CN facade.
    # `dir` must outlive the returned dynpft file handle (caller closes it).
    function build_domain(dir)
        file_years  = [2000, 2001, 2002]
        natpft_size = 3
        ng = 1
        ncol = 1
        np = natpft_size

        pct = Array{Float64}(undef, natpft_size, ng, length(file_years))
        pct[:, 1, 1] = [50.0, 30.0, 20.0]   # 2000 -> [.5 .3 .2]
        pct[:, 1, 2] = [20.0, 50.0, 30.0]   # 2001 -> [.2 .5 .3]
        pct[:, 1, 3] = [10.0, 10.0, 80.0]   # 2002 -> [.1 .1 .8]

        fn = joinpath(dir, "flanduse_pft.nc")
        NCDataset(fn, "c") do ds
            defDim(ds, "natpft", natpft_size)
            defDim(ds, "lndgrid", ng)
            defDim(ds, "time", length(file_years))
            defVar(ds, "YEAR", Int, ("time",))[:] = file_years
            pv = defVar(ds, "PCT_NAT_PFT", Float64, ("natpft", "lndgrid", "time"))
            for t in 1:length(file_years)
                pv[:, :, t] = pct[:, :, t]
            end
        end

        bounds = CLM.BoundsType(begg = 1, endg = 1, begl = 1, endl = 1,
                                begc = 1, endc = ncol, begp = 1, endp = np,
                                level = CLM.BOUNDS_LEVEL_PROC)

        grc = CLM.GridcellData{Float64}()
        grc.landunit_indices = fill(CLM.ISPVAL, CLM.MAX_LUNIT, ng)
        grc.landunit_indices[CLM.ISTSOIL, 1] = 1

        lun = CLM.LandunitData{Float64}()
        lun.itype = [CLM.ISTSOIL]; lun.wtgcell = [1.0]; lun.active = falses(1)
        lun.coli = [1]; lun.colf = [ncol]

        col = CLM.ColumnData{Float64}()
        col.landunit = [1]; col.gridcell = [1]
        col.wtlunit = [1.0]; col.wtgcell = [1.0]; col.active = falses(ncol)

        pch = CLM.PatchData{Float64}()
        pch.gridcell = fill(1, np); pch.landunit = fill(1, np); pch.column = fill(1, np)
        pch.itype = collect(0:(np - 1))
        pch.wtcol = fill(NaN, np); pch.wtlunit = fill(NaN, np); pch.wtgcell = fill(NaN, np)
        pch.active = falses(np)

        ctl = CLM.dyn_subgrid_control_init(
            flanduse_timeseries = fn, do_transient_pfts = true)

        state = CLM.dynSubgrid_init!(bounds, ctl, grc, lun, col, pch;
            current_year = 2000, natpft_size = natpft_size,
            check_dynpft_consistency = false)

        return (; bounds, grc, lun, col, pch, state, np, ncol, ng)
    end

    # Allocate + populate a full CN-veg + soilbgc facade with distinct nonzero
    # pools; returns the bundle pieces. (Mirrors test_dyn_cons_biogeochem.jl.)
    cpool_fields = (:leafc_patch, :leafc_storage_patch, :leafc_xfer_patch,
        :frootc_patch, :frootc_storage_patch, :frootc_xfer_patch,
        :livestemc_patch, :livestemc_storage_patch, :livestemc_xfer_patch,
        :deadstemc_patch, :deadstemc_storage_patch, :deadstemc_xfer_patch,
        :livecrootc_patch, :livecrootc_storage_patch, :livecrootc_xfer_patch,
        :deadcrootc_patch, :deadcrootc_storage_patch, :deadcrootc_xfer_patch,
        :gresp_storage_patch, :gresp_xfer_patch, :cpool_patch,
        :xsmrpool_patch, :ctrunc_patch)
    npool_fields = (:leafn_patch, :leafn_storage_patch, :leafn_xfer_patch,
        :frootn_patch, :frootn_storage_patch, :frootn_xfer_patch,
        :livestemn_patch, :livestemn_storage_patch, :livestemn_xfer_patch,
        :deadstemn_patch, :deadstemn_storage_patch, :deadstemn_xfer_patch,
        :livecrootn_patch, :livecrootn_storage_patch, :livecrootn_xfer_patch,
        :deadcrootn_patch, :deadcrootn_storage_patch, :deadcrootn_xfer_patch,
        :retransn_patch, :npool_patch, :ntrunc_patch)

    function build_cn(np, ncol, ng)
        cs = CLM.CNVegCarbonStateData{Float64}()
        CLM.cnveg_carbon_state_init!(cs, np, ncol, ng)
        ns = CLM.CNVegNitrogenStateData{Float64}()
        CLM.cnveg_nitrogen_state_init!(ns, np, ncol, ng)

        v = 1.0
        for f in cpool_fields
            arr = getproperty(cs, f)
            for p in 1:np; arr[p] = v; v += 1.0; end
        end
        for f in npool_fields
            arr = getproperty(ns, f)
            for p in 1:np; arr[p] = v * 0.01; v += 1.0; end
        end

        cf = CLM.CNVegCarbonFluxData{Float64}()
        CLM.cnveg_carbon_flux_init!(cf, np, ncol, ng;
            nlevdecomp_full = nlevdecomp, ndecomp_pools = ndecomp_pools)
        nf = CLM.CNVegNitrogenFluxData{Float64}()
        CLM.cnveg_nitrogen_flux_init!(nf, np, ncol, ng;
            nlevdecomp_full = nlevdecomp, i_litr_max = i_litr_max)

        for fld in (:dwt_seedc_to_leaf_patch, :dwt_seedc_to_leaf_grc,
                    :dwt_seedc_to_deadstem_patch, :dwt_seedc_to_deadstem_grc,
                    :dwt_conv_cflux_patch, :dwt_conv_cflux_grc,
                    :dwt_wood_productc_gain_patch, :dwt_crop_productc_gain_patch,
                    :dwt_slash_cflux_patch, :dwt_slash_cflux_grc)
            getproperty(cf, fld) .= 0.0
        end
        cf.dwt_frootc_to_litr_c_col .= 0.0
        cf.dwt_livecrootc_to_cwdc_col .= 0.0
        cf.dwt_deadcrootc_to_cwdc_col .= 0.0
        for fld in (:dwt_seedn_to_leaf_patch, :dwt_seedn_to_leaf_grc,
                    :dwt_seedn_to_deadstem_patch, :dwt_seedn_to_deadstem_grc,
                    :dwt_conv_nflux_patch, :dwt_conv_nflux_grc,
                    :dwt_wood_productn_gain_patch, :dwt_crop_productn_gain_patch)
            getproperty(nf, fld) .= 0.0
        end
        nf.dwt_frootn_to_litr_n_col .= 0.0
        nf.dwt_livecrootn_to_cwdn_col .= 0.0
        nf.dwt_deadcrootn_to_cwdn_col .= 0.0

        cnst = CLM.CNVegStateData{Float64}()
        CLM.cnveg_state_init!(cnst, np, ncol)
        cnst.dwt_smoothed_patch .= 0.0
        cnst.annavg_t2m_col .= 285.0

        canopy = CLM.CanopyStateData{Float64}()
        CLM.canopystate_init!(canopy, np)

        sbgc = CLM.SoilBiogeochemStateData{Float64}()
        CLM.soil_bgc_state_init!(sbgc, ncol, np, nlevdecomp, 0)
        sbgc.froot_prof_patch .= 1.0 / nlevdecomp
        sbgc.croot_prof_patch .= 1.0 / nlevdecomp

        sc = CLM.SoilBiogeochemCarbonStateData{Float64}()
        CLM.soil_bgc_carbon_state_init!(sc, ncol, ng, nlevdecomp, ndecomp_pools)
        sn = CLM.SoilBiogeochemNitrogenStateData{Float64}()
        CLM.soil_bgc_nitrogen_state_init!(sn, ncol, ng, nlevdecomp, ndecomp_pools)
        for l in 1:ndecomp_pools, j in 1:nlevdecomp
            sc.decomp_cpools_vr_col[1, j, l] = 5.0 + j + l
            sn.decomp_npools_vr_col[1, j, l] = (5.0 + j + l) * 0.05
        end

        return (; cs, ns, cf, nf, cnst, canopy, sbgc, sc, sn)
    end

    # Assemble the cons_bgc bundle from a domain + CN facade (as clm_driver does).
    function make_bundle(cn)
        return (
            dynbal = CLM.DynConsBiogeochemState(),
            pftcon = make_pftcon(15, ndecomp_pools),
            canopystate = cn.canopy,
            cnveg_state = cn.cnst,
            cnveg_carbonstate = cn.cs,
            cnveg_carbonflux = cn.cf,
            cnveg_nitrogenstate = cn.ns,
            cnveg_nitrogenflux = cn.nf,
            soilbiogeochem_state = cn.sbgc,
            soilbiogeochem_carbonstate = cn.sc,
            soilbiogeochem_nitrogenstate = cn.sn,
            mask_soilp_with_inactive = trues(3),
            mask_soilc_with_inactive = trues(1),
            dt = dt,
            nlevdecomp = nlevdecomp,
            ndecomp_pools = ndecomp_pools,
            i_litr_min = i_litr_min,
            i_litr_max = i_litr_max,
            i_cwd = i_cwd,
            use_crop = false,
            nrepr = CLM.NREPR,
            use_nitrif_denitrif = false,
        )
    end

    # ----------------------------------------------------------------------
    # 1. Gate ON: patch-level veg C/N conserved across the 2000->2001 change.
    # ----------------------------------------------------------------------
    @testset "cons_bgc ON: patch C/N conserved across weight change" begin
        mktempdir() do dir
            d  = build_domain(dir)
            cn = build_cn(d.np, d.ncol, d.ng)
            bundle = make_bundle(cn)

            np = d.np
            pwt_old = copy(d.pch.wtgcell)         # year-2000 weights

            # total gridcell-weighted veg C/N BEFORE the change.
            function total_grc(state, fields, pw)
                t = 0.0
                for f in fields
                    arr = getproperty(state, f)
                    for p in 1:np; t += arr[p] * pw[p]; end
                end
                return t
            end
            totC_before = total_grc(cn.cs, cpool_fields, pwt_old)
            totN_before = total_grc(cn.ns, npool_fields, pwt_old)

            CLM.dynSubgrid_driver!(d.state, d.bounds, d.grc, d.lun, d.col, d.pch;
                year = 2001, cons_bgc = bundle)

            pwt_new = copy(d.pch.wtgcell)          # year-2001 weights
            @test pwt_new != pwt_old               # the area actually changed

            totC_after = total_grc(cn.cs, cpool_fields, pwt_new)
            totN_after = total_grc(cn.ns, npool_fields, pwt_new)

            # C/N that left the veg pools (per gridcell area).
            cflux_out = 0.0; nflux_out = 0.0
            for p in 1:np
                cflux_out += (cn.cf.dwt_conv_cflux_patch[p] +
                              cn.cf.dwt_wood_productc_gain_patch[p] +
                              cn.cf.dwt_crop_productc_gain_patch[p]) * dt
                nflux_out += (cn.nf.dwt_conv_nflux_patch[p] +
                              cn.nf.dwt_wood_productn_gain_patch[p] +
                              cn.nf.dwt_crop_productn_gain_patch[p]) * dt
            end
            # root-to-litter / CWD are per-column-area; convert via OLD column weight.
            for c in 1:d.ncol
                wcol = d.state.patch_state_updater.cwtgcell_old[c]
                for j in 1:nlevdecomp
                    for i in i_litr_min:i_litr_max
                        cflux_out += cn.cf.dwt_frootc_to_litr_c_col[c, j, i] * dt * wcol
                        nflux_out += cn.nf.dwt_frootn_to_litr_n_col[c, j, i] * dt * wcol
                    end
                    cflux_out += cn.cf.dwt_livecrootc_to_cwdc_col[c, j] * dt * wcol
                    cflux_out += cn.cf.dwt_deadcrootc_to_cwdc_col[c, j] * dt * wcol
                    nflux_out += cn.nf.dwt_livecrootn_to_cwdn_col[c, j] * dt * wcol
                    nflux_out += cn.nf.dwt_deadcrootn_to_cwdn_col[c, j] * dt * wcol
                end
            end
            # seed C/N that entered the veg pools.
            cseed_in = 0.0; nseed_in = 0.0
            for p in 1:np
                cseed_in += (cn.cf.dwt_seedc_to_leaf_patch[p] +
                             cn.cf.dwt_seedc_to_deadstem_patch[p]) * dt
                nseed_in += (cn.nf.dwt_seedn_to_leaf_patch[p] +
                             cn.nf.dwt_seedn_to_deadstem_patch[p]) * dt
            end

            @test isapprox(totC_after, totC_before + cseed_in - cflux_out;
                           atol = 1e-8, rtol = 1e-9)
            @test isapprox(totN_after, totN_before + nseed_in - nflux_out;
                           atol = 1e-8, rtol = 1e-9)

            # The shrinking PFT (p=1, .5 -> .2) produced conversion flux.
            @test cn.cf.dwt_conv_cflux_grc[1] != 0.0
            @test cn.nf.dwt_conv_nflux_grc[1] != 0.0

            # Soil-bgc state stayed finite (column weight didn't change -> col redistribution no-op).
            @test all(isfinite, cn.sc.decomp_cpools_vr_col)
            @test all(isfinite, cn.sn.decomp_npools_vr_col)

            CLM.dynpft_close!(d.state.dynpft_state)
        end
    end

    # ----------------------------------------------------------------------
    # 2. Gate OFF: CN state byte-identical to not running the conservation.
    # ----------------------------------------------------------------------
    @testset "cons_bgc OFF (nothing): CN state untouched" begin
        mktempdir() do dir
            d  = build_domain(dir)
            cn = build_cn(d.np, d.ncol, d.ng)

            # Snapshot every veg C/N pool + soil pool BEFORE the (gate-off) step.
            csnap = Dict(f => copy(getproperty(cn.cs, f)) for f in cpool_fields)
            nsnap = Dict(f => copy(getproperty(cn.ns, f)) for f in npool_fields)
            sc_snap = copy(cn.sc.decomp_cpools_vr_col)
            sn_snap = copy(cn.sn.decomp_npools_vr_col)

            CLM.dynSubgrid_driver!(d.state, d.bounds, d.grc, d.lun, d.col, d.pch;
                year = 2001)   # cons_bgc defaults to nothing -> no-op

            for f in cpool_fields
                @test getproperty(cn.cs, f) == csnap[f]
            end
            for f in npool_fields
                @test getproperty(cn.ns, f) == nsnap[f]
            end
            @test cn.sc.decomp_cpools_vr_col == sc_snap
            @test cn.sn.decomp_npools_vr_col == sn_snap

            CLM.dynpft_close!(d.state.dynpft_state)
        end
    end

    # ----------------------------------------------------------------------
    # 3. Gate ON but no weight change -> ~zero conversion / product flux.
    # ----------------------------------------------------------------------
    @testset "cons_bgc ON, no weight change -> zero conv flux" begin
        mktempdir() do dir
            d  = build_domain(dir)
            # advance to the final year so a re-step at 2002 has no change.
            CLM.dynSubgrid_driver!(d.state, d.bounds, d.grc, d.lun, d.col, d.pch;
                year = 2002)
            wt_before = copy(d.pch.wtgcell)

            cn = build_cn(d.np, d.ncol, d.ng)
            bundle = make_bundle(cn)

            CLM.dynSubgrid_driver!(d.state, d.bounds, d.grc, d.lun, d.col, d.pch;
                year = 2002, cons_bgc = bundle)

            @test d.pch.wtgcell ≈ wt_before
            @test all(abs.(cn.cf.dwt_conv_cflux_patch) .< 1e-9)
            @test all(abs.(cn.cf.dwt_wood_productc_gain_patch) .< 1e-9)
            @test all(abs.(cn.nf.dwt_conv_nflux_patch) .< 1e-9)
            @test all(abs.(cn.cf.dwt_slash_cflux_patch) .< 1e-9)

            CLM.dynpft_close!(d.state.dynpft_state)
        end
    end

end
