@testset "Initialization Pipeline" begin
    using Dates

    fsurdat  = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc"
    paramfile = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc"

    # Skip if input files are not available
    if !isfile(fsurdat) || !isfile(paramfile)
        @warn "Skipping initialization tests: input files not found"
        @test true  # placeholder
    else
        @testset "surfrd_get_num_patches" begin
            (numpft, numcft) = CLM.surfrd_get_num_patches(fsurdat)
            @test numpft > 0
            @test numcft >= 0
        end

        @testset "surfrd_get_nlevurb" begin
            nlevurb = CLM.surfrd_get_nlevurb(fsurdat)
            @test nlevurb > 0
        end

        @testset "Full clm_initialize!" begin
            (inst, bounds, filt, tm) = CLM.clm_initialize!(
                fsurdat=fsurdat,
                paramfile=paramfile,
                start_date=DateTime(2000, 1, 1),
                dtime=1800,
                use_cn=false,
                use_crop=false)

            # ---- Bounds checks ----
            @test bounds.endg == 1  # single gridcell
            @test bounds.endl >= 2  # at least soil + lake landunits
            @test bounds.endc >= 2
            @test bounds.endp >= 1

            # ---- Gridcell checks ----
            @test !isnan(inst.gridcell.latdeg[1])
            @test !isnan(inst.gridcell.londeg[1])
            @test inst.gridcell.area[1] > 0.0

            # ---- Landunit checks ----
            @test any(inst.landunit.itype[1:bounds.endl] .== CLM.ISTSOIL)

            # ---- Vertical structure checks ----
            nlevsno = CLM.varpar.nlevsno
            nlevgrnd = CLM.varpar.nlevgrnd
            joff = nlevsno
            # Soil layer 1 depth should be > 0
            c_test = 1
            @test inst.column.z[c_test, 1 + joff] > 0.0
            @test !isnan(inst.column.z[c_test, 1 + joff])
            @test inst.column.dz[c_test, 1 + joff] > 0.0

            # Bedrock should be set
            @test inst.column.nbedrock[c_test] > 0
            @test inst.column.nbedrock[c_test] <= CLM.varpar.nlevsoi

            # ---- Temperature checks ----
            @test !isnan(inst.temperature.t_grnd_col[c_test])
            @test inst.temperature.t_grnd_col[c_test] > 200.0
            @test inst.temperature.t_grnd_col[c_test] < 400.0

            # ---- Soil properties checks ----
            @test inst.soilstate.watsat_col[c_test, 1] > 0.0
            @test inst.soilstate.watsat_col[c_test, 1] < 1.0

            # ---- PFT parameters checks ----
            @test length(CLM.pftcon.slatop) > 0
            # Check that at least some PFTs have non-zero slatop
            @test any(CLM.pftcon.slatop .> 0.0)

            # ---- Filter checks ----
            @test length(filt.soilc) == bounds.endc
            @test any(filt.soilc)  # at least one soil column

            # ---- Time manager checks ----
            @test CLM.get_nstep(tm) == 0
            @test CLM.get_curr_calday(tm) ≈ 1.0
            @test inst.soilhydrology.h2osfcflag == 0

            # ---- Active checks ----
            @test any(inst.column.active[1:bounds.endc])
            @test any(inst.patch.active[1:bounds.endp])

            # ---- Topography checks ----
            @test !isnan(inst.column.topo_slope[c_test])
            @test inst.column.topo_slope[c_test] > 0.0
            @test !isnan(inst.column.micro_sigma[c_test])
        end

        @testset "clm_initialize! h2osfcflag override" begin
            (inst1, _, _, _) = CLM.clm_initialize!(
                fsurdat=fsurdat,
                paramfile=paramfile,
                h2osfcflag=1)
            @test inst1.soilhydrology.h2osfcflag == 1
        end

        @testset "clm_initialize! use_aquifer_layer=false initializes ZWT at bedrock" begin
            (inst2, bounds2, _, _) = CLM.clm_initialize!(
                fsurdat=fsurdat,
                paramfile=paramfile,
                use_aquifer_layer=false)

            joff_zi = CLM.varpar.nlevsno + 1
            found_soil = false
            for c in bounds2.begc:bounds2.endc
                l = inst2.column.landunit[c]
                itype_l = inst2.landunit.itype[l]
                if itype_l == CLM.ISTSOIL || itype_l == CLM.ISTCROP
                    found_soil = true
                    nb = inst2.column.nbedrock[c]
                    expected = inst2.column.zi[c, joff_zi + nb]
                    @test isfinite(inst2.soilhydrology.zwt_col[c])
                    @test inst2.soilhydrology.zwt_col[c] ≈ expected atol=1e-10
                end
            end
            @test found_soil
        end
    end
end
