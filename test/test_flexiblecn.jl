@testset "NutrientCompetitionFlexibleCN" begin

    # =====================================================================
    # Helper: build state for a single non-crop woody patch, with the extra
    # fields the FlexibleCN method reads (canopystate LAI, soil mineral N,
    # t_scalar, npool, leaf/froot storage C & N). The same builder feeds both
    # the FlexibleCN and the CLM4.5-default methods so they can be compared.
    # =====================================================================
    function make_flexcn_data(; gpp=10.0, availc=8.0, annsum_npp=400.0,
            plant_ndemand=0.2, retransn_to_npool=0.05,
            c_allometry=1.5, n_allometry=0.06, fpg=0.8, retransn=0.5,
            annsum_potential_gpp=5000.0, annmax_retransn=1.0,
            use_fun=false, npool=0.4, leafc=10.0, leafn=0.4,
            leafc_storage=8.0, leafn_storage=0.32, frootc=5.0,
            sminn_vr=2.0, t_scalar=0.7, laisun=1.0, laisha=0.5)

        np = 1; nc = 1; nlevdecomp = 1
        nrepr = CLM.NREPR

        pftcon = CLM.PftConNutrientCompetition(
            woody = [1.0], froot_leaf = [1.0], croot_stem = [0.3],
            stem_leaf = [-1.0], flivewd = [0.1], leafcn = [25.0],
            frootcn = [42.0], livewdcn = [50.0], deadwdcn = [500.0],
            fcur = [0.5], graincn = [50.0], grperc = [0.3], grpnow = [0.5],
            fleafcn = [65.0], ffrootcn = [100.0], fstemcn = [130.0],
            astemf = [0.0], season_decid = [0.0], stress_decid = [0.0],
            evergreen = [0.0])

        cn_shared_params = CLM.CNSharedParamsData(use_matrixcn=false,
            use_fun=use_fun, use_flexiblecn=true)

        patch = CLM.PatchData(); patch.column = [1]; patch.itype = [0]
        crop = CLM.CropData(); crop.croplive_patch = [false]

        vs = CLM.CNVegStateData()
        vs.c_allometry_patch          = [c_allometry]
        vs.n_allometry_patch          = [n_allometry]
        vs.downreg_patch              = [0.0]
        vs.aleaf_patch = [NaN]; vs.astem_patch = [NaN]; vs.aroot_patch = [NaN]
        vs.arepr_patch = fill(NaN, np, nrepr)
        vs.peaklai_patch              = [0]
        vs.tempsum_potential_gpp_patch = [0.0]
        vs.annsum_potential_gpp_patch  = [annsum_potential_gpp]
        vs.tempmax_retransn_patch      = [0.0]
        vs.annmax_retransn_patch       = [annmax_retransn]
        vs.grain_flag_patch            = [0.0]

        cs = CLM.CNVegCarbonStateData()
        cs.leafc_patch = [leafc]; cs.frootc_patch = [frootc]
        cs.livestemc_patch = [20.0]; cs.livecrootc_patch = [15.0]
        cs.leafc_storage_patch = [leafc_storage]; cs.frootc_storage_patch = [4.0]
        cs.livestemc_storage_patch = [10.0]; cs.livecrootc_storage_patch = [7.0]

        cf = CLM.CNVegCarbonFluxData()
        cf.gpp_before_downreg_patch = [gpp]
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

        ns = CLM.CNVegNitrogenStateData()
        ns.retransn_patch = [retransn]; ns.npool_patch = [npool]
        ns.leafn_patch = [leafn]; ns.frootn_patch = [0.12]; ns.livestemn_patch = [0.4]
        ns.leafn_storage_patch = [leafn_storage]; ns.frootn_storage_patch = [0.1]
        ns.livestemn_storage_patch = [0.2]; ns.livecrootn_storage_patch = [0.14]

        nf = CLM.CNVegNitrogenFluxData()
        nf.plant_ndemand_patch = [plant_ndemand]; nf.avail_retransn_patch = [0.0]
        nf.retransn_to_npool_patch = [retransn_to_npool]
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

        canopystate = CLM.CanopyStateData()
        canopystate.laisun_patch = [laisun]; canopystate.laisha_patch = [laisha]

        soilbgc_ns = CLM.SoilBiogeochemNitrogenStateData()
        soilbgc_ns.sminn_vr_col = fill(sminn_vr, nc, nlevdecomp)
        soilbgc_cf = CLM.SoilBiogeochemCarbonFluxData()
        soilbgc_cf.t_scalar_col = fill(t_scalar, nc, nlevdecomp)

        return (pftcon=pftcon, cn_shared_params=cn_shared_params, patch=patch,
                crop=crop, canopystate=canopystate, cnveg_state=vs, cnveg_cs=cs,
                cnveg_cf=cf, cnveg_ns=ns, cnveg_nf=nf, soilbgc_ns=soilbgc_ns,
                soilbgc_cf=soilbgc_cf, mask_soilp=trues(np), fpg_col=[fpg],
                bounds=1:np, dzsoi_decomp=[1.0], nlevdecomp=nlevdecomp)
    end

    dt = 1800.0

    # =====================================================================
    # FlexibleCN plant N demand: finite, non-negative, Michaelis-Menten form
    # =====================================================================
    @testset "FlexibleCN N demand finite + physical" begin
        d = make_flexcn_data()
        CLM.calc_plant_nutrient_demand_flexiblecn!(d.mask_soilp, d.bounds, false,
            d.pftcon, d.cn_shared_params, d.patch, d.crop, d.canopystate,
            d.cnveg_state, d.cnveg_cs, d.cnveg_cf, d.cnveg_ns, d.cnveg_nf,
            d.soilbgc_ns, d.soilbgc_cf;
            dt=dt, dzsoi_decomp=d.dzsoi_decomp, nlevdecomp=d.nlevdecomp)

        p = 1
        @test isfinite(d.cnveg_nf.plant_ndemand_patch[p])
        @test d.cnveg_nf.plant_ndemand_patch[p] >= 0.0
        # tempsum/tempmax accumulators advanced
        @test d.cnveg_state.tempsum_potential_gpp_patch[p] == 10.0
        @test d.cnveg_state.tempmax_retransn_patch[p] == 0.5
        @test isfinite(d.cnveg_nf.avail_retransn_patch[p])
        @test d.cnveg_nf.retransn_to_npool_patch[p] >= 0.0
        # retrans drawn does not exceed storage rate
        @test d.cnveg_nf.avail_retransn_patch[p] <= d.cnveg_ns.retransn_patch[p] / dt + 1e-15
    end

    # =====================================================================
    # FlexibleCN N demand differs from CLM4.5 default (different physics)
    # =====================================================================
    @testset "FlexibleCN N demand differs from default" begin
        # FlexibleCN (Vmax_N kinetics) vs default (availc * n/c allometry)
        df = make_flexcn_data(); dd = make_flexcn_data()
        CLM.calc_plant_nutrient_demand_flexiblecn!(df.mask_soilp, df.bounds, false,
            df.pftcon, df.cn_shared_params, df.patch, df.crop, df.canopystate,
            df.cnveg_state, df.cnveg_cs, df.cnveg_cf, df.cnveg_ns, df.cnveg_nf,
            df.soilbgc_ns, df.soilbgc_cf;
            dt=dt, dzsoi_decomp=df.dzsoi_decomp, nlevdecomp=df.nlevdecomp)
        CLM.calc_plant_nitrogen_demand!(dd.mask_soilp, dd.bounds, false,
            dd.pftcon, dd.cn_shared_params, dd.patch, dd.crop,
            dd.cnveg_state, dd.cnveg_cs, dd.cnveg_cf, dd.cnveg_ns, dd.cnveg_nf;
            dt=dt)
        # the methods compute different plant_ndemand
        @test df.cnveg_nf.plant_ndemand_patch[1] != dd.cnveg_nf.plant_ndemand_patch[1]
    end

    # =====================================================================
    # FlexibleCN C/N allocation: finite, npool-conserving N allocation
    # =====================================================================
    @testset "FlexibleCN allocation finite + npool draw" begin
        d = make_flexcn_data()
        # need plant_ndemand & retransn_to_npool set (demand call first)
        CLM.calc_plant_nutrient_demand_flexiblecn!(d.mask_soilp, d.bounds, false,
            d.pftcon, d.cn_shared_params, d.patch, d.crop, d.canopystate,
            d.cnveg_state, d.cnveg_cs, d.cnveg_cf, d.cnveg_ns, d.cnveg_nf,
            d.soilbgc_ns, d.soilbgc_cf;
            dt=dt, dzsoi_decomp=d.dzsoi_decomp, nlevdecomp=d.nlevdecomp)

        CLM.calc_plant_nutrient_competition_flexiblecn!(d.mask_soilp, d.bounds,
            d.pftcon, d.cn_shared_params, d.patch, d.crop, d.canopystate,
            d.cnveg_state, d.cnveg_cs, d.cnveg_cf, d.cnveg_ns, d.cnveg_nf;
            fpg_col=d.fpg_col, dt=dt)

        p = 1
        @test isfinite(d.cnveg_cf.plant_calloc_patch[p])
        @test isfinite(d.cnveg_nf.npool_to_leafn_patch[p])
        @test d.cnveg_cf.cpool_to_leafc_patch[p] >= 0.0
        # plant_calloc = availc (FlexibleCN, non-FUN)
        @test d.cnveg_cf.plant_calloc_patch[p] ≈ d.cnveg_cf.availc_patch[p]

        # N-conservation: total npool_to_* drawn equals npool/dt (all N in the
        # pool is allocated when total demand > 0), since fractions sum to 1.
        npool = 0.4
        total_alloc = d.cnveg_nf.npool_to_leafn_patch[p] +
            d.cnveg_nf.npool_to_leafn_storage_patch[p] +
            d.cnveg_nf.npool_to_frootn_patch[p] +
            d.cnveg_nf.npool_to_frootn_storage_patch[p] +
            d.cnveg_nf.npool_to_livestemn_patch[p] +
            d.cnveg_nf.npool_to_livestemn_storage_patch[p] +
            d.cnveg_nf.npool_to_deadstemn_patch[p] +
            d.cnveg_nf.npool_to_deadstemn_storage_patch[p] +
            d.cnveg_nf.npool_to_livecrootn_patch[p] +
            d.cnveg_nf.npool_to_livecrootn_storage_patch[p] +
            d.cnveg_nf.npool_to_deadcrootn_patch[p] +
            d.cnveg_nf.npool_to_deadcrootn_storage_patch[p]
        @test total_alloc ≈ npool / dt rtol=1e-12
    end

    # =====================================================================
    # carbon_resp_opt: high actual storage-C:N triggers extra-C turnover
    # =====================================================================
    @testset "carbon_resp_opt high C:N turnover" begin
        # leafc_storage/leafn_storage = 40/0.32 ... make storage C:N very high
        # leafc_storage=80, leafn_storage=0.32 -> 250 >> leafcn(25)+15
        d = make_flexcn_data(leafc_storage=80.0, leafn_storage=0.32)
        CLM.calc_plant_nutrient_demand_flexiblecn!(d.mask_soilp, d.bounds, false,
            d.pftcon, d.cn_shared_params, d.patch, d.crop, d.canopystate,
            d.cnveg_state, d.cnveg_cs, d.cnveg_cf, d.cnveg_ns, d.cnveg_nf,
            d.soilbgc_ns, d.soilbgc_cf;
            dt=dt, dzsoi_decomp=d.dzsoi_decomp, nlevdecomp=d.nlevdecomp)
        CLM.calc_plant_nutrient_competition_flexiblecn!(d.mask_soilp, d.bounds,
            d.pftcon, d.cn_shared_params, d.patch, d.crop, d.canopystate,
            d.cnveg_state, d.cnveg_cs, d.cnveg_cf, d.cnveg_ns, d.cnveg_nf;
            fpg_col=d.fpg_col, dt=dt, carbon_resp_opt=1)
        p = 1
        # high leaf storage C:N → cpool_to_leafc_resp > 0
        @test d.cnveg_cf.cpool_to_leafc_resp_patch[p] > 0.0
        @test d.cnveg_cf.cpool_to_resp_patch[p] >= d.cnveg_cf.cpool_to_leafc_resp_patch[p]
        @test isfinite(d.cnveg_cf.cpool_to_resp_patch[p])

        # carbon_resp_opt=0 leaves all resp terms at zero
        d0 = make_flexcn_data(leafc_storage=80.0, leafn_storage=0.32)
        CLM.calc_plant_nutrient_demand_flexiblecn!(d0.mask_soilp, d0.bounds, false,
            d0.pftcon, d0.cn_shared_params, d0.patch, d0.crop, d0.canopystate,
            d0.cnveg_state, d0.cnveg_cs, d0.cnveg_cf, d0.cnveg_ns, d0.cnveg_nf,
            d0.soilbgc_ns, d0.soilbgc_cf;
            dt=dt, dzsoi_decomp=d0.dzsoi_decomp, nlevdecomp=d0.nlevdecomp)
        CLM.calc_plant_nutrient_competition_flexiblecn!(d0.mask_soilp, d0.bounds,
            d0.pftcon, d0.cn_shared_params, d0.patch, d0.crop, d0.canopystate,
            d0.cnveg_state, d0.cnveg_cs, d0.cnveg_cf, d0.cnveg_ns, d0.cnveg_nf;
            fpg_col=d0.fpg_col, dt=dt, carbon_resp_opt=0)
        @test d0.cnveg_cf.cpool_to_leafc_resp_patch[1] == 0.0
        @test d0.cnveg_cf.cpool_to_resp_patch[1] == 0.0
    end

    # =====================================================================
    # FlexibleCN allocation differs from CLM4.5 default where stoichiometry
    # flexes (plant_calloc form + npool-draw N allocation).
    # =====================================================================
    @testset "FlexibleCN allocation differs from default" begin
        df = make_flexcn_data(); dd = make_flexcn_data()
        # Same demand inputs for both
        for d in (df, dd)
            d.cnveg_nf.plant_ndemand_patch[1] = 0.2
            d.cnveg_nf.retransn_to_npool_patch[1] = 0.05
        end
        CLM.calc_plant_nutrient_competition_flexiblecn!(df.mask_soilp, df.bounds,
            df.pftcon, df.cn_shared_params, df.patch, df.crop, df.canopystate,
            df.cnveg_state, df.cnveg_cs, df.cnveg_cf, df.cnveg_ns, df.cnveg_nf;
            fpg_col=df.fpg_col, dt=dt)
        # default method (no flexiblecn) on the same inputs
        CLM.calc_plant_nutrient_competition!(dd.mask_soilp, dd.bounds,
            dd.pftcon, dd.cn_shared_params, dd.patch, dd.crop, dd.cnveg_state,
            dd.cnveg_cs, dd.cnveg_cf, dd.cnveg_nf;
            fpg_col=dd.fpg_col, cnveg_ns=nothing, dt=dt)

        # plant_calloc differs: FlexibleCN = availc; default = nalloc*c/n
        @test df.cnveg_cf.plant_calloc_patch[1] != dd.cnveg_cf.plant_calloc_patch[1]
        # leaf N allocation differs (npool-draw vs C-paired demand)
        @test df.cnveg_nf.npool_to_leafn_patch[1] != dd.cnveg_nf.npool_to_leafn_patch[1]
    end

end
