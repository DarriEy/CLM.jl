# =====================================================================
# collapse_crop_types! on the STATIC surfdata path
#
# Fortran calls collapse_crop_types TWICE — surfrdMod.F90:991 (static surfdata
# read) and dyncropFileMod.F90:176 (transient landuse). CLM.jl only ever had the
# transient call, so a static crop run never collapsed the 64-CFT set: it built
# 22 crop patches where Fortran builds 8, with every itype shifted +1.
#
# The reference numbers below come from the CN-crop Fortran run documented in
# docs/CROP_PARITY.md (US Great Plains point, irrigate=.true.).
# =====================================================================

@testset "collapse_crop_types! — static surfdata path" begin

    # ------------------------------------------------------------------
    # Part 1: self-contained. Exercises the merge arithmetic and the
    # ordering guard without needing any external data file.
    # ------------------------------------------------------------------
    @testset "merge arithmetic + ordering guard" begin
        old_merge = copy(CLM.pftcon.mergetoclmpft)
        old_known = copy(CLM.pftcon.is_pft_known_to_model)
        old_irrigate = CLM.varctl.irrigate
        old_maxveg = CLM.varpar.maxveg
        old_cft_lb = CLM.varpar.cft_lb
        try
            # 5 PFTs (0..4); CFTs live at 3..4. Map PFT 4 -> PFT 3.
            CLM.varpar.cft_lb = 3
            CLM.varpar.maxveg = 4
            CLM.varctl.irrigate = true          # skip the rainfed/irrigated merge
            CLM.pftcon.mergetoclmpft = [0, 1, 2, 3, 3]
            CLM.pftcon.is_pft_known_to_model = fill(true, 5)

            #                       cft 3   cft 4
            wt_cft   = reshape([0.25, 0.75], 1, 2)
            fert_cft = reshape([10.0, 20.0], 1, 2)

            CLM.collapse_crop_types!(wt_cft, fert_cft, 1, 1)

            # All weight collapses onto CFT 3; CFT 4 is emptied.
            @test wt_cft[1, 1] ≈ 1.0
            @test wt_cft[1, 2] == 0.0
            # Fertilizer is the weight-weighted mean of the two sources:
            # (0.25*10 + 0.75*20) / 1.0 = 17.5
            @test fert_cft[1, 1] ≈ 17.5
            # The merged-away type is no longer modelled.
            @test CLM.pftcon.is_pft_known_to_model[5] == false
            @test CLM.pftcon.is_pft_known_to_model[4] == true

            # A missing merge map must FAIL, not silently pass the uncollapsed
            # set through — that silent no-op is the bug this routine prevents.
            CLM.pftcon.mergetoclmpft = Int[]
            @test_throws ErrorException CLM.collapse_crop_types!(
                reshape([0.5, 0.5], 1, 2), reshape([0.0, 0.0], 1, 2), 1, 1)
        finally
            CLM.pftcon.mergetoclmpft = old_merge
            CLM.pftcon.is_pft_known_to_model = old_known
            CLM.varctl.irrigate = old_irrigate
            CLM.varpar.maxveg = old_maxveg
            CLM.varpar.cft_lb = old_cft_lb
        end
    end

    # ------------------------------------------------------------------
    # Part 2: against the real crop-CFT surfdata and the Fortran reference.
    # Gated on the data files, which are not in the repo.
    # ------------------------------------------------------------------
    _dataroot = joinpath(homedir(), "Library/CloudStorage",
                         "GoogleDrive-dareyt@gmail.com", "My Drive", "data",
                         "SYMFLUENCE_data")
    _surfdata = joinpath(_dataroot, "crop_cft_surfdata",
                         "surfdata_cropCFT_USplains_1pt.nc")
    _paramfile = joinpath(_dataroot, "domain_Aripuana_Amazon", "settings",
                          "CLM", "parameters", "clm5_params.nc")

    if !(isfile(_surfdata) && isfile(_paramfile))
        @info "skipping crop-CFT reference test (data files absent)" _surfdata
    else
        @testset "US-plains crop point matches the Fortran patch vector" begin
            # Fortran's 8 crop CFT patches. itypes are the crops CLM has
            # parameterizations for: temperate corn (17/18), spring wheat
            # (19/20), soybean (23/24) and the generic crop (75/76); each is a
            # rainfed/irrigated pair, both of which survive because the
            # reference run sets irrigate=.true.
            fortran_itypes = [17, 18, 19, 20, 23, 24, 75, 76]
            # Collapsed weights: every one of the 22 nonzero input CFTs lands on
            # one of these 8 via mergetoclmpft.
            fortran_weights = [0.311670, 0.004981, 0.319289, 0.002998,
                               0.353584, 0.003512, 0.003966, 0.0]

            old_crop = CLM.varctl.use_crop
            old_clu = CLM.varctl.create_crop_landunit
            old_irrigate = CLM.varctl.irrigate
            try
                CLM.varctl.use_crop = true
                CLM.varctl.create_crop_landunit = true
                CLM.varctl.irrigate = true

                (numpft, numcft) = CLM.surfrd_get_num_patches(_surfdata)
                nlevurb = CLM.surfrd_get_nlevurb(_surfdata)

                # numpft is the highest natural PFT INDEX (0-based), so the
                # 15-long natpft dimension gives 14 — matching Fortran's
                # natpft_ub. Returning the raw dimension here is what pushed
                # cft_lb to 16 and maxveg to 79.
                @test numpft == 14
                @test numcft == 64

                CLM.varpar_init!(CLM.varpar, 1, numpft, numcft, nlevurb)

                # Exactly Fortran's clm_varpar values.
                @test CLM.varpar.natpft_lb == 0
                @test CLM.varpar.natpft_ub == 14
                @test CLM.varpar.cft_lb == 15
                @test CLM.varpar.cft_ub == 78
                @test CLM.varpar.maxveg == 78

                CLM.varcon_init!()
                CLM.readParameters!(_paramfile)

                # The merge map must be in bounds for the m in 1:maxveg loop,
                # which indexes mergetoclmpft[m+1]. With maxveg=79 this read
                # index 80 of a 79-element array.
                @test length(CLM.pftcon.mergetoclmpft) == 79
                @test CLM.varpar.maxveg + 1 <= length(CLM.pftcon.mergetoclmpft)

                surf = CLM.SurfaceInputData()
                CLM.surfrd_get_data!(surf, 1, 1, _surfdata)

                # The collapse must conserve total crop weight exactly.
                @test sum(surf.wt_cft[1, :]) ≈ 1.0 atol = 1e-12

                # Build the real subgrid and read the actual patch vector — this
                # is the deliverable: the crop patches CLM.jl constructs, not a
                # derived quantity.
                ng = 1
                (nl, nc, np) = CLM.count_subgrid_elements(surf, ng)
                bounds = CLM.BoundsType(begg=1, endg=ng, begl=1, endl=nl,
                                        begc=1, endc=nc, begp=1, endp=np,
                                        level=CLM.BOUNDS_LEVEL_CLUMP, clump_index=1)
                grc = CLM.GridcellData(); CLM.gridcell_init!(grc, ng)
                lun = CLM.LandunitData(); CLM.landunit_init!(lun, nl)
                col = CLM.ColumnData();   CLM.column_init!(col, nc)
                pch = CLM.PatchData();    CLM.patch_init!(pch, np)
                CLM.initGridCells!(bounds, surf, grc, lun, col, pch)

                crop_p = [p for p in 1:np
                          if lun.itype[col.landunit[pch.column[p]]] == CLM.ISTCROP]

                # 8 crop patches, not 22.
                @test length(crop_p) == length(fortran_itypes)
                @test [pch.itype[p] for p in crop_p] == fortran_itypes

                # Weights live on the crop columns (each crop patch is the whole
                # of its column), and must match Fortran cft-for-cft.
                got_w = [col.wtlunit[pch.column[p]] for p in crop_p]
                for (i, w) in enumerate(fortran_weights)
                    @test got_w[i] ≈ w atol = 1e-6
                end
                @test sum(got_w) ≈ 1.0 atol = 1e-12
            finally
                CLM.varctl.use_crop = old_crop
                CLM.varctl.create_crop_landunit = old_clu
                CLM.varctl.irrigate = old_irrigate
            end
        end
    end
end
