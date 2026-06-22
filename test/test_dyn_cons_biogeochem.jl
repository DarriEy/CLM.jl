@testset "dyn_cons_biogeochem" begin

    tol = 1.0e-9

    # ----------------------------------------------------------------------
    # Helper: minimal PftconType with the fields the C/N balance touches
    # (fr_f, pconv, woody, evergreen, leafcn, deadwdcn, c3psn, noveg).
    # ndecomp_pools = 3 -> i_litr_min=1, i_litr_max=2 (litter pools 1..2).
    # ----------------------------------------------------------------------
    function make_pftcon(n::Int, ndecomp_pools::Int)
        p = CLM.PftconType()
        p.c3psn    = ones(n)
        p.evergreen = zeros(n)
        p.woody    = zeros(n)
        p.leafcn   = fill(25.0, n)
        p.deadwdcn = fill(100.0, n)
        p.noveg    = fill(typemax(Int), n)
        p.pconv    = fill(0.5, n)                 # half of deadstem -> conv, half -> wood product
        p.fr_f     = zeros(n, ndecomp_pools)
        # fine root litter fractions: split across the 2 litter pools, summing to 1
        for j in 1:n
            p.fr_f[j, 1] = 0.6
            p.fr_f[j, 2] = 0.4
        end
        return p
    end

    # ----------------------------------------------------------------------
    # Patch-level test: 3 soil patches on 3 columns in one gridcell.
    #   p=1 shrinking   : 0.40 -> 0.20
    #   p=2 growing     : 0.10 -> 0.30
    #   p=3 initiating  : 0.00 -> 0.25
    # One patch per column; column weight mirrors patch weight.
    # ----------------------------------------------------------------------
    @testset "dyn_cnbal_patch! C & N conservation" begin
        CLM.varpar_init!(CLM.varpar, 1, 15, 2, 5)
        CLM.varcon_init!()

        ndecomp_pools = 3
        nlevdecomp    = CLM.varpar.nlevdecomp
        i_litr_min = 1
        i_litr_max = 2

        npatch = 3
        ncol   = 3
        ng     = 1

        bounds = CLM.BoundsType(begg=1, endg=ng, begl=1, endl=1,
                                begc=1, endc=ncol, begp=1, endp=npatch,
                                level=CLM.BOUNDS_LEVEL_PROC)

        pwt_old = [0.40, 0.10, 0.00]
        pwt_new = [0.20, 0.30, 0.25]

        pch = CLM.PatchData{Float64}()
        pch.column   = collect(1:npatch)
        pch.landunit = fill(1, npatch)
        pch.gridcell = fill(1, npatch)
        pch.itype    = [7, 13, 8]          # all vegetated (non-bareground) PFTs
        pch.wtgcell  = copy(pwt_old)

        col = CLM.ColumnData{Float64}()
        col.gridcell = fill(1, ncol)
        col.wtgcell  = copy(pwt_old)

        lun = CLM.LandunitData()
        lun.itype = [CLM.ISTSOIL]

        pftcon = make_pftcon(15, ndecomp_pools)

        # All patches are soil patches (in the filter, incl. inactive).
        mask_soilp = trues(npatch)

        # --- updater snapshots ---
        updater = CLM.PatchStateUpdater(bounds)
        CLM.set_old_weights!(updater, bounds, pch, col)
        pch.wtgcell .= pwt_new
        col.wtgcell .= pwt_new
        CLM.set_new_weights!(updater, bounds, pch)

        # --- allocate + fill C/N veg state with nonzero pools ---
        cs = CLM.CNVegCarbonStateData{Float64}()
        CLM.cnveg_carbon_state_init!(cs, npatch, ncol, ng)
        ns = CLM.CNVegNitrogenStateData{Float64}()
        CLM.cnveg_nitrogen_state_init!(ns, npatch, ncol, ng)

        # Helper to fill all the patch C pools the routine reads/writes.
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

        # Distinct nonzero values so accidental aliasing would show up.
        v = 1.0
        for f in cpool_fields
            arr = getproperty(cs, f)
            for p in 1:npatch
                arr[p] = v; v += 1.0
            end
        end
        for f in npool_fields
            arr = getproperty(ns, f)
            for p in 1:npatch
                arr[p] = v * 0.01; v += 1.0
            end
        end

        # --- C/N flux structs: zero out the dwt accumulator fields used ---
        cf = CLM.CNVegCarbonFluxData{Float64}()
        CLM.cnveg_carbon_flux_init!(cf, npatch, ncol, ng; nlevdecomp_full=nlevdecomp, ndecomp_pools=ndecomp_pools)
        nf = CLM.CNVegNitrogenFluxData{Float64}()
        CLM.cnveg_nitrogen_flux_init!(nf, npatch, ncol, ng; nlevdecomp_full=nlevdecomp, i_litr_max=i_litr_max)

        # zero the gridcell/patch/col dwt accumulators (init fills with NaN)
        for fld in (:dwt_seedc_to_leaf_patch, :dwt_seedc_to_leaf_grc,
                    :dwt_seedc_to_deadstem_patch, :dwt_seedc_to_deadstem_grc,
                    :dwt_conv_cflux_patch, :dwt_conv_cflux_grc,
                    :dwt_wood_productc_gain_patch, :dwt_crop_productc_gain_patch,
                    :dwt_slash_cflux_patch, :dwt_slash_cflux_grc)
            getproperty(cf, fld) .= 0.0
        end
        cf.dwt_frootc_to_litr_c_col   .= 0.0
        cf.dwt_livecrootc_to_cwdc_col .= 0.0
        cf.dwt_deadcrootc_to_cwdc_col .= 0.0
        for fld in (:dwt_seedn_to_leaf_patch, :dwt_seedn_to_leaf_grc,
                    :dwt_seedn_to_deadstem_patch, :dwt_seedn_to_deadstem_grc,
                    :dwt_conv_nflux_patch, :dwt_conv_nflux_grc,
                    :dwt_wood_productn_gain_patch, :dwt_crop_productn_gain_patch)
            getproperty(nf, fld) .= 0.0
        end
        nf.dwt_frootn_to_litr_n_col   .= 0.0
        nf.dwt_livecrootn_to_cwdn_col .= 0.0
        nf.dwt_deadcrootn_to_cwdn_col .= 0.0

        # --- cnveg_state + canopy + soilbgc state ---
        cnst = CLM.CNVegStateData{Float64}()
        CLM.cnveg_state_init!(cnst, npatch, ncol)
        cnst.dwt_smoothed_patch .= 0.0
        cnst.annavg_t2m_col     .= 285.0

        canopy = CLM.CanopyStateData{Float64}()
        CLM.canopystate_init!(canopy, npatch)

        sbgc = CLM.SoilBiogeochemStateData{Float64}()
        CLM.soil_bgc_state_init!(sbgc, ncol, npatch, nlevdecomp, 0)
        # uniform root profiles that integrate to ~1 over the decomp levels
        sbgc.froot_prof_patch .= 1.0 / nlevdecomp
        sbgc.croot_prof_patch .= 1.0 / nlevdecomp

        dynbal = CLM.DynConsBiogeochemState()
        dt = 1800.0

        # --- total gridcell-weighted veg C/N BEFORE ---
        function total_grc(state, fields, pw)
            t = 0.0
            for f in fields
                arr = getproperty(state, f)
                for p in 1:npatch
                    t += arr[p] * pw[p]
                end
            end
            return t
        end
        totC_before = total_grc(cs, cpool_fields, pwt_old)
        totN_before = total_grc(ns, npool_fields, pwt_old)

        CLM.dyn_cnbal_patch!(dynbal, bounds, mask_soilp, pwt_old, updater,
            pch, lun, col, pftcon, canopy, cnst, cs, cf, ns, nf, sbgc;
            dt=dt, i_litr_min=i_litr_min, i_litr_max=i_litr_max)

        totC_after = total_grc(cs, cpool_fields, pwt_new)
        totN_after = total_grc(ns, npool_fields, pwt_new)

        # --- C/N that LEFT the veg pools (per gridcell area), via the flux
        #     fields the module produced (all *dt to undo the /dt). ---
        # conversion + wood-product + crop-product (per-gridcell-area patch flux)
        cflux_out = 0.0
        nflux_out = 0.0
        for p in 1:npatch
            cflux_out += (cf.dwt_conv_cflux_patch[p] +
                          cf.dwt_wood_productc_gain_patch[p] +
                          cf.dwt_crop_productc_gain_patch[p]) * dt
            nflux_out += (nf.dwt_conv_nflux_patch[p] +
                          nf.dwt_wood_productn_gain_patch[p] +
                          nf.dwt_crop_productn_gain_patch[p]) * dt
        end
        # root-to-litter / CWD: these are per-COLUMN-area fluxes computed using
        # the OLD column weight (flux_out_col_area divides by cwtgcell_old), so
        # convert to gridcell area with the OLD column weight.
        for c in 1:ncol
            wcol = updater.cwtgcell_old[c]
            for j in 1:nlevdecomp
                for i in i_litr_min:i_litr_max
                    cflux_out += cf.dwt_frootc_to_litr_c_col[c, j, i] * dt * wcol
                    nflux_out += nf.dwt_frootn_to_litr_n_col[c, j, i] * dt * wcol
                end
                cflux_out += cf.dwt_livecrootc_to_cwdc_col[c, j] * dt * wcol
                cflux_out += cf.dwt_deadcrootc_to_cwdc_col[c, j] * dt * wcol
                nflux_out += nf.dwt_livecrootn_to_cwdn_col[c, j] * dt * wcol
                nflux_out += nf.dwt_deadcrootn_to_cwdn_col[c, j] * dt * wcol
            end
        end
        # seed C/N that ENTERED the veg pools (per gridcell area).
        cseed_in = 0.0
        nseed_in = 0.0
        for p in 1:npatch
            cseed_in += (cf.dwt_seedc_to_leaf_patch[p] +
                         cf.dwt_seedc_to_deadstem_patch[p]) * dt
            nseed_in += (nf.dwt_seedn_to_leaf_patch[p] +
                         nf.dwt_seedn_to_deadstem_patch[p]) * dt
        end

        # Conservation: after = before + seed_in - flux_out
        @test isapprox(totC_after, totC_before + cseed_in - cflux_out; atol=tol,
                       rtol=1e-9)
        @test isapprox(totN_after, totN_before + nseed_in - nflux_out; atol=tol,
                       rtol=1e-9)

        # Initiating patch (p=3) had its phenology/flux state reset.
        @test cnst.dormant_flag_patch[3] == 1.0
        @test cnst.days_active_patch[3] == 0.0
        @test cnst.annavg_t2m_patch[3] == cnst.annavg_t2m_col[1]
        @test cf.availc_patch[3] == 0.0
        @test canopy.laisun_patch[3] == 0.0

        # Some conversion flux WAS produced (shrinking p=1 lost C/N).
        @test cf.dwt_conv_cflux_grc[1] != 0.0
        @test nf.dwt_conv_nflux_grc[1] != 0.0
    end

    # ----------------------------------------------------------------------
    # No-weight-change -> ~0 conversion / product fluxes.
    # ----------------------------------------------------------------------
    @testset "dyn_cnbal_patch! no weight change -> zero conv flux" begin
        CLM.varpar_init!(CLM.varpar, 1, 15, 2, 5)
        CLM.varcon_init!()
        ndecomp_pools = 3
        nlevdecomp = CLM.varpar.nlevdecomp
        i_litr_min = 1; i_litr_max = 2

        npatch = 2; ncol = 2; ng = 1
        bounds = CLM.BoundsType(begg=1, endg=ng, begl=1, endl=1,
                                begc=1, endc=ncol, begp=1, endp=npatch,
                                level=CLM.BOUNDS_LEVEL_PROC)
        pwt = [0.3, 0.3]
        pch = CLM.PatchData{Float64}()
        pch.column = collect(1:npatch); pch.landunit = fill(1, npatch)
        pch.gridcell = fill(1, npatch); pch.itype = [7, 13]; pch.wtgcell = copy(pwt)
        col = CLM.ColumnData{Float64}(); col.gridcell = fill(1, ncol); col.wtgcell = copy(pwt)
        lun = CLM.LandunitData(); lun.itype = [CLM.ISTSOIL]
        pftcon = make_pftcon(15, ndecomp_pools)
        mask_soilp = trues(npatch)

        updater = CLM.PatchStateUpdater(bounds)
        CLM.set_old_weights!(updater, bounds, pch, col)
        CLM.set_new_weights!(updater, bounds, pch)   # new == old

        cs = CLM.CNVegCarbonStateData{Float64}(); CLM.cnveg_carbon_state_init!(cs, npatch, ncol, ng)
        ns = CLM.CNVegNitrogenStateData{Float64}(); CLM.cnveg_nitrogen_state_init!(ns, npatch, ncol, ng)
        for f in (:leafc_patch,:leafc_storage_patch,:leafc_xfer_patch,:frootc_patch,
                  :frootc_storage_patch,:frootc_xfer_patch,:livestemc_patch,
                  :livestemc_storage_patch,:livestemc_xfer_patch,:deadstemc_patch,
                  :deadstemc_storage_patch,:deadstemc_xfer_patch,:livecrootc_patch,
                  :livecrootc_storage_patch,:livecrootc_xfer_patch,:deadcrootc_patch,
                  :deadcrootc_storage_patch,:deadcrootc_xfer_patch,:gresp_storage_patch,
                  :gresp_xfer_patch,:cpool_patch,:xsmrpool_patch,:ctrunc_patch)
            getproperty(cs, f) .= 5.0
        end
        for f in (:leafn_patch,:leafn_storage_patch,:leafn_xfer_patch,:frootn_patch,
                  :frootn_storage_patch,:frootn_xfer_patch,:livestemn_patch,
                  :livestemn_storage_patch,:livestemn_xfer_patch,:deadstemn_patch,
                  :deadstemn_storage_patch,:deadstemn_xfer_patch,:livecrootn_patch,
                  :livecrootn_storage_patch,:livecrootn_xfer_patch,:deadcrootn_patch,
                  :deadcrootn_storage_patch,:deadcrootn_xfer_patch,:retransn_patch,
                  :npool_patch,:ntrunc_patch)
            getproperty(ns, f) .= 0.2
        end

        cf = CLM.CNVegCarbonFluxData{Float64}(); CLM.cnveg_carbon_flux_init!(cf, npatch, ncol, ng; nlevdecomp_full=nlevdecomp, ndecomp_pools=ndecomp_pools)
        nf = CLM.CNVegNitrogenFluxData{Float64}(); CLM.cnveg_nitrogen_flux_init!(nf, npatch, ncol, ng; nlevdecomp_full=nlevdecomp, i_litr_max=i_litr_max)
        for fld in (:dwt_seedc_to_leaf_patch,:dwt_seedc_to_leaf_grc,:dwt_seedc_to_deadstem_patch,
                    :dwt_seedc_to_deadstem_grc,:dwt_conv_cflux_patch,:dwt_conv_cflux_grc,
                    :dwt_wood_productc_gain_patch,:dwt_crop_productc_gain_patch,
                    :dwt_slash_cflux_patch,:dwt_slash_cflux_grc)
            getproperty(cf, fld) .= 0.0
        end
        cf.dwt_frootc_to_litr_c_col .= 0.0; cf.dwt_livecrootc_to_cwdc_col .= 0.0
        cf.dwt_deadcrootc_to_cwdc_col .= 0.0
        for fld in (:dwt_seedn_to_leaf_patch,:dwt_seedn_to_leaf_grc,:dwt_seedn_to_deadstem_patch,
                    :dwt_seedn_to_deadstem_grc,:dwt_conv_nflux_patch,:dwt_conv_nflux_grc,
                    :dwt_wood_productn_gain_patch,:dwt_crop_productn_gain_patch)
            getproperty(nf, fld) .= 0.0
        end
        nf.dwt_frootn_to_litr_n_col .= 0.0; nf.dwt_livecrootn_to_cwdn_col .= 0.0
        nf.dwt_deadcrootn_to_cwdn_col .= 0.0

        cnst = CLM.CNVegStateData{Float64}(); CLM.cnveg_state_init!(cnst, npatch, ncol)
        cnst.dwt_smoothed_patch .= 0.0; cnst.annavg_t2m_col .= 285.0
        canopy = CLM.CanopyStateData{Float64}(); CLM.canopystate_init!(canopy, npatch)
        sbgc = CLM.SoilBiogeochemStateData{Float64}(); CLM.soil_bgc_state_init!(sbgc, ncol, npatch, nlevdecomp, 0)
        sbgc.froot_prof_patch .= 1.0 / nlevdecomp; sbgc.croot_prof_patch .= 1.0 / nlevdecomp

        dynbal = CLM.DynConsBiogeochemState()
        CLM.dyn_cnbal_patch!(dynbal, bounds, mask_soilp, pwt, updater,
            pch, lun, col, pftcon, canopy, cnst, cs, cf, ns, nf, sbgc;
            dt=1800.0, i_litr_min=1, i_litr_max=2)

        @test all(abs.(cf.dwt_conv_cflux_patch) .< tol)
        @test all(abs.(cf.dwt_wood_productc_gain_patch) .< tol)
        @test all(abs.(nf.dwt_conv_nflux_patch) .< tol)
        @test all(abs.(cf.dwt_slash_cflux_patch) .< tol)
    end

    # ----------------------------------------------------------------------
    # Column-level test: 2 natveg columns on one gridcell, area swap
    # (c1 shrinks, c2 grows, equal magnitude) -> soil C/N conserved.
    # ----------------------------------------------------------------------
    @testset "dyn_cnbal_col! soil C & N conservation" begin
        CLM.varpar_init!(CLM.varpar, 1, 15, 2, 5)
        CLM.varcon_init!()
        nlevdecomp = CLM.varpar.nlevdecomp
        ndecomp_pools = 3

        ng, nl, nc, np = 1, 1, 2, 2
        grc = CLM.GridcellData(); lun = CLM.LandunitData()
        col = CLM.ColumnData(); pch = CLM.PatchData()
        CLM.gridcell_init!(grc, ng); CLM.landunit_init!(lun, nl)
        CLM.column_init!(col, nc); CLM.patch_init!(pch, np)

        li = Ref(0); ci = Ref(0); pi = Ref(0)
        CLM.add_landunit!(lun, li, 1, CLM.ISTSOIL, 1.0)
        CLM.add_column!(col, lun, ci, li[], CLM.ISTSOIL, 0.5)
        CLM.add_patch!(pch, col, lun, pi, ci[], 1, 1.0)
        CLM.add_column!(col, lun, ci, li[], CLM.ISTSOIL, 0.5)
        CLM.add_patch!(pch, col, lun, pi, ci[], 1, 1.0)

        bounds = CLM.BoundsType(begg=1, endg=ng, begl=1, endl=nl,
                                begc=1, endc=nc, begp=1, endp=np,
                                level=CLM.BOUNDS_LEVEL_PROC)
        CLM.clm_ptrs_compdown!(bounds, grc, lun, col, pch)

        csu = CLM.ColumnStateUpdater(bounds, 1)
        col.active .= [true, true]
        col.wtgcell .= [0.5, 0.5]
        CLM.set_old_weights!(csu, bounds, grc, lun, col)
        col.wtgcell .= [0.3, 0.7]      # c1 shrinks 0.2, c2 grows 0.2
        CLM.set_new_weights!(csu, bounds, 1, col)

        sc = CLM.SoilBiogeochemCarbonStateData{Float64}()
        CLM.soil_bgc_carbon_state_init!(sc, nc, ng, nlevdecomp, ndecomp_pools)
        sn = CLM.SoilBiogeochemNitrogenStateData{Float64}()
        CLM.soil_bgc_nitrogen_state_init!(sn, nc, ng, nlevdecomp, ndecomp_pools)

        # Fill the vr pools with distinct nonzero values per column.
        for l in 1:ndecomp_pools, j in 1:nlevdecomp
            sc.decomp_cpools_vr_col[1, j, l] = 10.0 + j + l
            sc.decomp_cpools_vr_col[2, j, l] = 3.0  + j + l
            sn.decomp_npools_vr_col[1, j, l] = (10.0 + j + l) * 0.05
            sn.decomp_npools_vr_col[2, j, l] = (3.0  + j + l) * 0.05
        end
        for j in 1:nlevdecomp
            sc.ctrunc_vr_col[1, j] = 0.5; sc.ctrunc_vr_col[2, j] = 0.2
            sn.ntrunc_vr_col[1, j] = 0.01; sn.ntrunc_vr_col[2, j] = 0.005
            sn.sminn_vr_col[1, j]  = 1.0; sn.sminn_vr_col[2, j]  = 0.4
        end

        dzs = CLM.dzsoi_decomp[]

        # Gridcell-weighted, depth-integrated soil C BEFORE.
        function soilC_total()
            t = 0.0
            for c in 1:nc, l in 1:ndecomp_pools, j in 1:nlevdecomp
                t += sc.decomp_cpools_vr_col[c, j, l] * dzs[j] * col.wtgcell[c]
            end
            for c in 1:nc, j in 1:nlevdecomp
                t += sc.ctrunc_vr_col[c, j] * dzs[j] * col.wtgcell[c]
            end
            return t
        end
        function soilN_total()
            t = 0.0
            for c in 1:nc, l in 1:ndecomp_pools, j in 1:nlevdecomp
                t += sn.decomp_npools_vr_col[c, j, l] * dzs[j] * col.wtgcell[c]
            end
            for c in 1:nc, j in 1:nlevdecomp
                t += (sn.ntrunc_vr_col[c, j] + sn.sminn_vr_col[c, j]) * dzs[j] * col.wtgcell[c]
            end
            return t
        end

        # BEFORE totals use the NEW weights for the post-redistribution sum and
        # OLD weights for the pre sum. We snapshot pre with the OLD weights.
        col.wtgcell .= [0.5, 0.5]
        C_before = soilC_total()
        N_before = soilN_total()
        col.wtgcell .= [0.3, 0.7]

        CLM.dyn_cnbal_col!(bounds, 1, csu, col, sc, sn)

        C_after = soilC_total()
        N_after = soilN_total()

        @test isapprox(C_after, C_before; atol=1e-8, rtol=1e-9)
        @test isapprox(N_after, N_before; atol=1e-8, rtol=1e-9)

        # The per-column depth-integrated apparent state change is tracked
        # (nonzero on the growing column, finite everywhere).
        @test all(isfinite, sc.dyn_cbal_adjustments_col)
        @test all(isfinite, sn.dyn_nbal_adjustments_col)
        @test sc.dyn_cbal_adjustments_col[2] != 0.0   # c2 grew
        @test sn.dyn_nbal_adjustments_col[2] != 0.0
    end

end
