@testset "NutrientCompetitionFlexibleCN matrixcn + AgSys" begin

    dt = 1800.0
    NCPOOL = CLM.NVEGPOOL_NATVEG + CLM.NVEGPOOL_CROP   # 21, covers igrain_st=20

    # =====================================================================
    # Builder: one woody non-crop patch (matrixcn) OR one crop patch (agsys).
    # The matrix_* flux arrays are allocated (the raw struct constructors leave
    # them empty); everything else mirrors test_flexiblecn.jl's make_flexcn_data.
    # =====================================================================
    function make_data(; is_crop=false, use_matrixcn=false, use_crop_agsys=false,
            use_c13=false, use_c14=false,
            availc=8.0, annsum_npp=400.0, c_allometry=1.5, n_allometry=0.06,
            fpg=0.8, npool=0.4)

        np = 1; nc = 1; nlevdecomp = 1
        nrepr = CLM.NREPR

        # itype 0 => generic non-crop; itype >= npcropmin-1 (0-based) => crop.
        # npcropmin (1-based Julia) default; a crop ivt is npcropmin in 1-based,
        # i.e. itype = npcropmin-1 in 0-based.
        npcropmin = CLM.NPCROPMIN
        itype0 = is_crop ? (npcropmin - 1) : 0
        ivt = itype0 + 1

        # pftcon must be indexable at ivt
        n_pft = ivt
        woody = fill(is_crop ? 0.0 : 1.0, n_pft)
        pftcon = CLM.PftConNutrientCompetition(
            woody = woody,
            froot_leaf = fill(1.0, n_pft), croot_stem = fill(0.3, n_pft),
            stem_leaf = fill(-1.0, n_pft), flivewd = fill(0.1, n_pft),
            leafcn = fill(25.0, n_pft), frootcn = fill(42.0, n_pft),
            livewdcn = fill(50.0, n_pft), deadwdcn = fill(500.0, n_pft),
            fcur = fill(0.5, n_pft), graincn = fill(50.0, n_pft),
            grperc = fill(0.3, n_pft), grpnow = fill(0.5, n_pft),
            fleafcn = fill(65.0, n_pft), ffrootcn = fill(100.0, n_pft),
            fstemcn = fill(130.0, n_pft), astemf = fill(0.0, n_pft),
            season_decid = fill(0.0, n_pft), stress_decid = fill(0.0, n_pft),
            evergreen = fill(0.0, n_pft))

        cn_shared_params = CLM.CNSharedParamsData(use_matrixcn=use_matrixcn,
            use_fun=false, use_flexiblecn=true)

        patch = CLM.PatchData(); patch.column = [1]; patch.itype = [itype0]
        crop = CLM.CropData(); crop.croplive_patch = [is_crop]

        vs = CLM.CNVegStateData()
        vs.c_allometry_patch          = [c_allometry]
        vs.n_allometry_patch          = [n_allometry]
        vs.downreg_patch              = [0.0]
        if is_crop
            vs.aleaf_patch = [0.4]; vs.astem_patch = [0.3]; vs.aroot_patch = [0.2]
            vs.arepr_patch = fill(0.1, np, nrepr)
        else
            vs.aleaf_patch = [NaN]; vs.astem_patch = [NaN]; vs.aroot_patch = [NaN]
            vs.arepr_patch = fill(NaN, np, nrepr)
        end
        # AgSys N allocation coefficients
        vs.aleaf_n_patch = [0.35]; vs.astem_n_patch = [0.25]; vs.aroot_n_patch = [0.30]
        vs.arepr_n_patch = fill(0.10, np, nrepr)
        vs.peaklai_patch              = [0]
        vs.tempsum_potential_gpp_patch = [0.0]
        vs.annsum_potential_gpp_patch  = [5000.0]
        vs.tempmax_retransn_patch      = [0.0]
        vs.annmax_retransn_patch       = [1.0]
        vs.grain_flag_patch            = [0.0]

        cs = CLM.CNVegCarbonStateData()
        cs.leafc_patch = [10.0]; cs.frootc_patch = [5.0]
        cs.livestemc_patch = [20.0]; cs.livecrootc_patch = [15.0]
        cs.leafc_storage_patch = [8.0]; cs.frootc_storage_patch = [4.0]
        cs.livestemc_storage_patch = [10.0]; cs.livecrootc_storage_patch = [7.0]

        cf = CLM.CNVegCarbonFluxData()
        cf.gpp_before_downreg_patch = [10.0]
        cf.availc_patch = [availc]; cf.npp_growth_patch = [availc]
        cf.excess_cflux_patch = [0.0]; cf.plant_calloc_patch = [0.0]
        cf.psnsun_to_cpool_patch = [6.0]; cf.psnshade_to_cpool_patch = [4.0]
        cf.annsum_npp_patch = [annsum_npp]
        for fld in (:cpool_to_leafc_patch, :cpool_to_leafc_storage_patch,
                    :cpool_to_frootc_patch, :cpool_to_frootc_storage_patch,
                    :cpool_to_livestemc_patch, :cpool_to_livestemc_storage_patch,
                    :cpool_to_deadstemc_patch, :cpool_to_deadstemc_storage_patch,
                    :cpool_to_livecrootc_patch, :cpool_to_livecrootc_storage_patch,
                    :cpool_to_deadcrootc_patch, :cpool_to_deadcrootc_storage_patch,
                    :cpool_to_gresp_storage_patch,
                    :cpool_to_resp_patch, :cpool_to_leafc_resp_patch,
                    :cpool_to_leafc_storage_resp_patch, :cpool_to_frootc_resp_patch,
                    :cpool_to_frootc_storage_resp_patch, :cpool_to_livecrootc_resp_patch,
                    :cpool_to_livecrootc_storage_resp_patch, :cpool_to_livestemc_resp_patch,
                    :cpool_to_livestemc_storage_resp_patch)
            setfield!(cf, fld, [0.0])
        end
        cf.cpool_to_reproductivec_patch         = fill(0.0, np, nrepr)
        cf.cpool_to_reproductivec_storage_patch = fill(0.0, np, nrepr)
        # matrix flux arrays
        cf.matrix_Cinput_patch   = [0.0]
        cf.matrix_C13input_patch = [0.0]
        cf.matrix_C14input_patch = [0.0]
        cf.matrix_alloc_patch    = fill(0.0, np, NCPOOL)

        # isotope carbon flux instances (only psn fields read)
        c13 = CLM.CNVegCarbonFluxData()
        c13.psnsun_to_cpool_patch = [3.0]; c13.psnshade_to_cpool_patch = [2.0]
        c14 = CLM.CNVegCarbonFluxData()
        c14.psnsun_to_cpool_patch = [1.0]; c14.psnshade_to_cpool_patch = [0.5]

        ns = CLM.CNVegNitrogenStateData()
        ns.retransn_patch = [0.5]; ns.npool_patch = [npool]
        ns.leafn_patch = [0.4]; ns.frootn_patch = [0.12]; ns.livestemn_patch = [0.4]
        ns.leafn_storage_patch = [0.32]; ns.frootn_storage_patch = [0.1]
        ns.livestemn_storage_patch = [0.2]; ns.livecrootn_storage_patch = [0.14]

        nf = CLM.CNVegNitrogenFluxData()
        nf.plant_ndemand_patch = [0.2]; nf.avail_retransn_patch = [0.0]
        nf.retransn_to_npool_patch = [0.05]
        nf.sminn_to_npool_patch = [0.0]; nf.plant_nalloc_patch = [0.0]
        for fld in (:npool_to_leafn_patch, :npool_to_leafn_storage_patch,
                    :npool_to_frootn_patch, :npool_to_frootn_storage_patch,
                    :npool_to_livestemn_patch, :npool_to_livestemn_storage_patch,
                    :npool_to_deadstemn_patch, :npool_to_deadstemn_storage_patch,
                    :npool_to_livecrootn_patch, :npool_to_livecrootn_storage_patch,
                    :npool_to_deadcrootn_patch, :npool_to_deadcrootn_storage_patch,
                    :sminn_to_plant_fun_patch, :leafn_to_retransn_patch,
                    :frootn_to_retransn_patch, :livestemn_to_retransn_patch)
            setfield!(nf, fld, [0.0])
        end
        nf.npool_to_reproductiven_patch         = fill(0.0, np, nrepr)
        nf.npool_to_reproductiven_storage_patch = fill(0.0, np, nrepr)
        nf.matrix_Ninput_patch = [0.0]
        nf.matrix_nalloc_patch = fill(0.0, np, NCPOOL)

        canopystate = CLM.CanopyStateData()
        canopystate.laisun_patch = [1.0]; canopystate.laisha_patch = [0.5]

        return (pftcon=pftcon, cn_shared_params=cn_shared_params, patch=patch,
                crop=crop, canopystate=canopystate, cnveg_state=vs, cnveg_cs=cs,
                cnveg_cf=cf, cnveg_ns=ns, cnveg_nf=nf, c13=c13, c14=c14,
                mask_soilp=trues(np), fpg_col=[fpg], bounds=1:np,
                use_crop_agsys=use_crop_agsys, ivt=ivt, nrepr=nrepr,
                use_c13=use_c13, use_c14=use_c14)
    end

    run_alloc!(d) = CLM.calc_plant_nutrient_competition_flexiblecn!(d.mask_soilp,
        d.bounds, d.pftcon, d.cn_shared_params, d.patch, d.crop, d.canopystate,
        d.cnveg_state, d.cnveg_cs, d.cnveg_cf, d.cnveg_ns, d.cnveg_nf;
        fpg_col=d.fpg_col, dt=dt,
        c13_cnveg_cf=d.c13, c14_cnveg_cf=d.c14,
        use_c13=d.use_c13, use_c14=d.use_c14,
        use_crop_agsys=d.use_crop_agsys)

    # =====================================================================
    # use_matrixcn writes the matrix inputs/B-matrix WITHOUT changing the
    # sequential cpool_to_* / npool_to_* fluxes (matrix is a reformulation).
    # =====================================================================
    @testset "matrixcn writes B-matrix; sequential fluxes unchanged" begin
        dm = make_data(is_crop=false, use_matrixcn=true)
        ds = make_data(is_crop=false, use_matrixcn=false)
        run_alloc!(dm); run_alloc!(ds)
        p = 1

        # sequential C/N fluxes identical between matrix-on and matrix-off
        for fld in (:cpool_to_leafc_patch, :cpool_to_leafc_storage_patch,
                    :cpool_to_frootc_patch, :cpool_to_livestemc_patch,
                    :cpool_to_deadstemc_patch, :cpool_to_livecrootc_patch,
                    :cpool_to_gresp_storage_patch)
            @test getfield(dm.cnveg_cf, fld)[p] == getfield(ds.cnveg_cf, fld)[p]
        end
        for fld in (:npool_to_leafn_patch, :npool_to_frootn_patch,
                    :npool_to_livestemn_patch, :npool_to_deadcrootn_storage_patch)
            @test getfield(dm.cnveg_nf, fld)[p] == getfield(ds.cnveg_nf, fld)[p]
        end

        # matrix C input = total C allocated to veg (cpool_to_veg, Fortran 531);
        # matrix_alloc fractions sum to 1
        cpool_to_veg = dm.cnveg_cf.cpool_to_leafc_patch[p] + dm.cnveg_cf.cpool_to_leafc_storage_patch[p] +
            dm.cnveg_cf.cpool_to_frootc_patch[p] + dm.cnveg_cf.cpool_to_frootc_storage_patch[p] +
            dm.cnveg_cf.cpool_to_livestemc_patch[p] + dm.cnveg_cf.cpool_to_livestemc_storage_patch[p] +
            dm.cnveg_cf.cpool_to_deadstemc_patch[p] + dm.cnveg_cf.cpool_to_deadstemc_storage_patch[p] +
            dm.cnveg_cf.cpool_to_livecrootc_patch[p] + dm.cnveg_cf.cpool_to_livecrootc_storage_patch[p] +
            dm.cnveg_cf.cpool_to_deadcrootc_patch[p] + dm.cnveg_cf.cpool_to_deadcrootc_storage_patch[p]
        @test dm.cnveg_cf.matrix_Cinput_patch[p] ≈ cpool_to_veg
        @test sum(dm.cnveg_cf.matrix_alloc_patch[p, :]) ≈ 1.0 rtol=1e-12

        # matrix_alloc * Cinput reconstructs the per-pool C flux (B-matrix identity)
        @test dm.cnveg_cf.matrix_alloc_patch[p, CLM.FLEXCN_ILEAF] *
              dm.cnveg_cf.matrix_Cinput_patch[p] ≈ dm.cnveg_cf.cpool_to_leafc_patch[p]
        @test dm.cnveg_cf.matrix_alloc_patch[p, CLM.FLEXCN_IFROOT] *
              dm.cnveg_cf.matrix_Cinput_patch[p] ≈ dm.cnveg_cf.cpool_to_frootc_patch[p]

        # matrix N B-matrix: nalloc fractions sum to 1; Ninput consistent
        @test sum(dm.cnveg_nf.matrix_nalloc_patch[p, :]) ≈ 1.0 rtol=1e-12
        npool_to_veg = dm.cnveg_nf.npool_to_leafn_patch[p] + dm.cnveg_nf.npool_to_leafn_storage_patch[p] +
            dm.cnveg_nf.npool_to_frootn_patch[p] + dm.cnveg_nf.npool_to_frootn_storage_patch[p] +
            dm.cnveg_nf.npool_to_livestemn_patch[p] + dm.cnveg_nf.npool_to_livestemn_storage_patch[p] +
            dm.cnveg_nf.npool_to_deadstemn_patch[p] + dm.cnveg_nf.npool_to_deadstemn_storage_patch[p] +
            dm.cnveg_nf.npool_to_livecrootn_patch[p] + dm.cnveg_nf.npool_to_livecrootn_storage_patch[p] +
            dm.cnveg_nf.npool_to_deadcrootn_patch[p] + dm.cnveg_nf.npool_to_deadcrootn_storage_patch[p]
        @test dm.cnveg_nf.matrix_Ninput_patch[p] ≈ npool_to_veg - dm.cnveg_nf.retransn_to_npool_patch[p]

        # matrix-off run leaves the matrix arrays untouched (all zero)
        @test ds.cnveg_cf.matrix_Cinput_patch[p] == 0.0
        @test all(ds.cnveg_cf.matrix_alloc_patch[p, :] .== 0.0)
        @test all(ds.cnveg_nf.matrix_nalloc_patch[p, :] .== 0.0)
    end

    # =====================================================================
    # matrixcn isotope C inputs (use_c13 / use_c14)
    # =====================================================================
    @testset "matrixcn C13/C14 inputs" begin
        d = make_data(is_crop=false, use_matrixcn=true, use_c13=true, use_c14=true)
        run_alloc!(d)
        p = 1
        psn = d.cnveg_cf.psnsun_to_cpool_patch[p] + d.cnveg_cf.psnshade_to_cpool_patch[p]
        c13ratio = (d.c13.psnsun_to_cpool_patch[p] + d.c13.psnshade_to_cpool_patch[p]) / psn
        c14ratio = (d.c14.psnsun_to_cpool_patch[p] + d.c14.psnshade_to_cpool_patch[p]) / psn
        @test d.cnveg_cf.matrix_C13input_patch[p] ≈ d.cnveg_cf.plant_calloc_patch[p] * c13ratio
        @test d.cnveg_cf.matrix_C14input_patch[p] ≈ d.cnveg_cf.plant_calloc_patch[p] * c14ratio
    end

    # =====================================================================
    # matrixcn carbon_resp_opt: cpool_to_resp removed from matrix C input
    # =====================================================================
    @testset "matrixcn carbon_resp_opt removes resp from Cinput" begin
        d = make_data(is_crop=false, use_matrixcn=true)
        # drive high leaf storage C:N so the resp turnover branch fires
        d.cnveg_cs.leafc_storage_patch[1] = 80.0
        d.cnveg_ns.leafn_storage_patch[1] = 0.32
        CLM.calc_plant_nutrient_competition_flexiblecn!(d.mask_soilp, d.bounds,
            d.pftcon, d.cn_shared_params, d.patch, d.crop, d.canopystate,
            d.cnveg_state, d.cnveg_cs, d.cnveg_cf, d.cnveg_ns, d.cnveg_nf;
            fpg_col=d.fpg_col, dt=dt, carbon_resp_opt=1)
        p = 1
        @test d.cnveg_cf.cpool_to_resp_patch[p] > 0.0
        # Cinput = cpool_to_veg - cpool_to_resp
        cpool_to_veg = d.cnveg_cf.cpool_to_leafc_patch[p] + d.cnveg_cf.cpool_to_leafc_storage_patch[p] +
            d.cnveg_cf.cpool_to_frootc_patch[p] + d.cnveg_cf.cpool_to_frootc_storage_patch[p] +
            d.cnveg_cf.cpool_to_livestemc_patch[p] + d.cnveg_cf.cpool_to_livestemc_storage_patch[p] +
            d.cnveg_cf.cpool_to_deadstemc_patch[p] + d.cnveg_cf.cpool_to_deadstemc_storage_patch[p] +
            d.cnveg_cf.cpool_to_livecrootc_patch[p] + d.cnveg_cf.cpool_to_livecrootc_storage_patch[p] +
            d.cnveg_cf.cpool_to_deadcrootc_patch[p] + d.cnveg_cf.cpool_to_deadcrootc_storage_patch[p]
        @test d.cnveg_cf.matrix_Cinput_patch[p] ≈ cpool_to_veg - d.cnveg_cf.cpool_to_resp_patch[p]
    end

    # =====================================================================
    # AgSys: crop N allocation drawn from AgSys coefficients (aleaf_n etc.),
    # differs from the FlexibleCN fractional-demand draw.
    # =====================================================================
    @testset "AgSys crop N allocation from coefficients" begin
        da = make_data(is_crop=true, use_crop_agsys=true)
        df = make_data(is_crop=true, use_crop_agsys=false)
        run_alloc!(da); run_alloc!(df)
        p = 1; npool = 0.4
        fcur = 0.5; f4 = 0.1  # flivewd

        # AgSys formula: npool_to_leafn = aleaf_n * fcur * npool/dt
        @test da.cnveg_nf.npool_to_leafn_patch[p] ≈ 0.35 * fcur * npool / dt
        @test da.cnveg_nf.npool_to_leafn_storage_patch[p] ≈ 0.35 * (1 - fcur) * npool / dt
        @test da.cnveg_nf.npool_to_frootn_patch[p] ≈ 0.30 * fcur * npool / dt
        @test da.cnveg_nf.npool_to_livestemn_patch[p] ≈ 0.25 * f4 * fcur * npool / dt
        @test da.cnveg_nf.npool_to_deadstemn_patch[p] ≈ 0.25 * (1 - f4) * fcur * npool / dt
        @test da.cnveg_nf.npool_to_reproductiven_patch[p, 1] ≈ 0.10 * fcur * npool / dt
        # AgSys: no coarse-root allocation for crops
        @test da.cnveg_nf.npool_to_livecrootn_patch[p] == 0.0
        @test da.cnveg_nf.npool_to_deadcrootn_patch[p] == 0.0

        # AgSys differs from the FlexibleCN fractional-demand draw
        @test da.cnveg_nf.npool_to_leafn_patch[p] != df.cnveg_nf.npool_to_leafn_patch[p]

        # the C allocation fluxes are NOT changed by use_crop_agsys (only N draw)
        @test da.cnveg_cf.cpool_to_leafc_patch[p] == df.cnveg_cf.cpool_to_leafc_patch[p]
    end

    # =====================================================================
    # AgSys default OFF: a crop patch without use_crop_agsys uses FlexibleCN.
    # =====================================================================
    @testset "AgSys off => FlexibleCN crop draw" begin
        df = make_data(is_crop=true, use_crop_agsys=false)
        run_alloc!(df)
        p = 1
        # FlexibleCN crop draw is finite and conserves npool when demand>0
        total = df.cnveg_nf.npool_to_leafn_patch[p] + df.cnveg_nf.npool_to_leafn_storage_patch[p] +
            df.cnveg_nf.npool_to_frootn_patch[p] + df.cnveg_nf.npool_to_frootn_storage_patch[p] +
            df.cnveg_nf.npool_to_livestemn_patch[p] + df.cnveg_nf.npool_to_livestemn_storage_patch[p] +
            df.cnveg_nf.npool_to_deadstemn_patch[p] + df.cnveg_nf.npool_to_deadstemn_storage_patch[p] +
            df.cnveg_nf.npool_to_livecrootn_patch[p] + df.cnveg_nf.npool_to_livecrootn_storage_patch[p] +
            df.cnveg_nf.npool_to_deadcrootn_patch[p] + df.cnveg_nf.npool_to_deadcrootn_storage_patch[p] +
            df.cnveg_nf.npool_to_reproductiven_patch[p, 1] + df.cnveg_nf.npool_to_reproductiven_storage_patch[p, 1]
        @test total ≈ 0.4 / dt rtol=1e-12
        @test all(isfinite, (df.cnveg_nf.npool_to_leafn_patch[p],
                             df.cnveg_nf.npool_to_reproductiven_patch[p, 1]))
    end

end
