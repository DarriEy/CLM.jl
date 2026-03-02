@testset "FUN (Fixation and Uptake of Nitrogen)" begin

    # =====================================================================
    # Test: PftConFUN defaults
    # =====================================================================
    @testset "PftConFUN defaults" begin
        pfc = CLM.PftConFUN()
        @test isempty(pfc.leafcn)
        @test isempty(pfc.a_fix)
        @test isempty(pfc.perecm)
        @test isempty(pfc.c3psn)
    end

    # =====================================================================
    # Test: FUNParams defaults
    # =====================================================================
    @testset "FUNParams defaults" begin
        fp = CLM.FUNParams()
        @test fp.ndays_off == 21.0
    end

    # =====================================================================
    # Test: fun_cost_fix — fixers vs non-fixers
    # =====================================================================
    @testset "fun_cost_fix" begin
        big_cost = CLM.BIG_COST

        # Non-fixer should return big_cost
        @test CLM.fun_cost_fix(0, -3.62, 0.27, 25.15, big_cost, 0.5, -6.0, 20.0) == big_cost

        # Fixer with negligible root fraction → big_cost
        @test CLM.fun_cost_fix(1, -3.62, 0.27, 25.15, big_cost, 1e-7, -6.0, 20.0) == big_cost

        # Fixer with real roots → finite cost
        cost = CLM.fun_cost_fix(1, -3.62, 0.27, 25.15, big_cost, 0.5, -6.0, 20.0)
        @test isfinite(cost)
        @test cost > 0.0

        # Cost should increase at low temperatures
        cost_warm = CLM.fun_cost_fix(1, -3.62, 0.27, 25.15, big_cost, 0.5, -6.0, 25.0)
        cost_cold = CLM.fun_cost_fix(1, -3.62, 0.27, 25.15, big_cost, 0.5, -6.0, 5.0)
        @test cost_cold > cost_warm   # colder = more expensive
    end

    # =====================================================================
    # Test: fun_cost_active
    # =====================================================================
    @testset "fun_cost_active" begin
        big_cost = CLM.BIG_COST
        sv = CLM.SMALLVALUE

        # Low root density → big_cost
        @test CLM.fun_cost_active(1.0, big_cost, 0.5, 0.5, 1e-7, 0.5, sv) == big_cost

        # Low N → big_cost
        @test CLM.fun_cost_active(1e-13, big_cost, 0.5, 0.5, 100.0, 0.5, sv) == big_cost

        # Normal conditions → finite cost
        cost = CLM.fun_cost_active(0.5, big_cost, 10.0, 5.0, 100.0, 0.5, sv)
        @test isfinite(cost)
        @test cost > 0.0
        # Check formula: kn/sminn + kc/rootc
        @test cost ≈ 5.0/0.5 + 10.0/100.0

        # More N = cheaper
        cost_high_n = CLM.fun_cost_active(2.0, big_cost, 10.0, 5.0, 100.0, 0.5, sv)
        @test cost_high_n < cost
    end

    # =====================================================================
    # Test: fun_cost_nonmyc
    # =====================================================================
    @testset "fun_cost_nonmyc" begin
        big_cost = CLM.BIG_COST
        sv = CLM.SMALLVALUE

        @test CLM.fun_cost_nonmyc(1.0, big_cost, 0.5, 0.5, 1e-7, 0.5, sv) == big_cost
        @test CLM.fun_cost_nonmyc(1e-13, big_cost, 0.5, 0.5, 100.0, 0.5, sv) == big_cost

        cost = CLM.fun_cost_nonmyc(0.5, big_cost, 10.0, 5.0, 100.0, 0.5, sv)
        @test cost ≈ 5.0/0.5 + 10.0/100.0
    end

    # =====================================================================
    # Test: fun_retranslocation
    # =====================================================================
    @testset "fun_retranslocation" begin
        dt = 1800.0
        npp_to_spend = 10.0          # gC/m2
        total_falling_leaf_c = 5.0   # gC/m2
        total_falling_leaf_n = 0.2   # gN/m2
        total_n_resistance = 50.0    # gC/gN
        target_leafcn = 25.0         # gC/gN
        grperc = 0.3
        plantCN = 30.0

        result = CLM.fun_retranslocation(
            dt, npp_to_spend, total_falling_leaf_c, total_falling_leaf_n,
            total_n_resistance, target_leafcn, grperc, plantCN)

        # Should return named tuple
        @test haskey(result, :total_c_spent_retrans)
        @test haskey(result, :total_c_accounted_retrans)
        @test haskey(result, :free_n_retrans)
        @test haskey(result, :paid_for_n_retrans)

        # Free N should be non-negative
        @test result.free_n_retrans >= 0.0

        # Total C spent should be non-negative (or very close to zero)
        @test result.total_c_spent_retrans >= -1e-10

        # Paid-for N should be non-negative
        @test result.paid_for_n_retrans >= 0.0

        # Free retrans: falling_leaf_n - falling_leaf_c/min_falling_leaf_cn
        min_cn = target_leafcn * 1.5
        expected_free = max(total_falling_leaf_n - total_falling_leaf_c / min_cn, 0.0)
        @test result.free_n_retrans ≈ expected_free
    end

    # =====================================================================
    # Test: fun_retranslocation with zero NPP
    # =====================================================================
    @testset "fun_retranslocation zero NPP" begin
        result = CLM.fun_retranslocation(
            1800.0, 0.0, 5.0, 0.2, 50.0, 25.0, 0.3, 30.0)

        # With zero NPP, only free retrans should occur
        @test result.total_c_spent_retrans ≈ 0.0
        @test result.paid_for_n_retrans ≈ 0.0
    end

    # =====================================================================
    # Helper: create minimal test data for CNFUN
    # =====================================================================
    function make_fun_test_data(; np=1, nc=1, nlevdecomp=2,
                                  dt=1800.0,
                                  availc=0.001,     # gC/m2/s
                                  leafc=100.0,
                                  frootc=50.0,
                                  leafn=4.0,
                                  frootn=2.0,
                                  livestemn=1.0,
                                  livecrootn=0.5,
                                  plantCN=30.0,
                                  c_allometry=1.0,
                                  n_allometry=0.033,
                                  smin_no3=0.5,
                                  smin_nh4=0.3,
                                  t_soil=290.0,
                                  h2osoi_liq=10.0,
                                  crootfr_vals=[0.7, 0.3])

        nlevdecomp_full = nlevdecomp

        # Ensure varpar is initialized
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

        # PFT constants (1 PFT)
        pftcon = CLM.PftConFUN(
            leafcn        = [25.0],
            season_decid  = [0.0],
            stress_decid  = [0.0],
            a_fix         = [-3.62],
            b_fix         = [0.27],
            c_fix         = [25.15],
            s_fix         = [-6.0],
            akc_active    = [1.0],
            akn_active    = [1.0],
            ekc_active    = [1.5],
            ekn_active    = [1.5],
            kc_nonmyc     = [2.0],
            kn_nonmyc     = [2.0],
            perecm        = [0.5],
            grperc        = [0.3],
            fun_cn_flex_a = [5.0],
            fun_cn_flex_b = [200.0],
            fun_cn_flex_c = [80.0],
            FUN_fracfixers = [0.3],
            c3psn         = [1.0],
        )

        fun_params = CLM.FUNParams(ndays_off=21.0)

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.itype[1] = 1
        patch.column[1] = 1
        patch.wtcol[1] = 1.0

        nl = 1  # landunits
        ng = 1  # gridcells

        # Water state — directly set fields we need
        waterstate = CLM.WaterStateData()
        CLM.waterstate_init!(waterstate, nc, np, nl, ng)
        for j in 1:nlevdecomp
            waterstate.h2osoi_liq_col[1, j] = h2osoi_liq
        end

        # Temperature
        temperature = CLM.TemperatureData()
        CLM.temperature_init!(temperature, np, nc, nl, ng)
        for j in 1:nlevdecomp
            temperature.t_soisno_col[1, j] = t_soil
        end

        # Soil state
        soilstate = CLM.SoilStateData()
        CLM.soilstate_init!(soilstate, np, nc)
        for j in 1:nlevdecomp
            soilstate.crootfr_patch[1, j] = crootfr_vals[j]
        end

        # Canopy state
        canopystate = CLM.CanopyStateData()
        CLM.canopystate_init!(canopystate, np)
        canopystate.tlai_patch[1] = 3.0

        # CN Veg state
        cnveg_state = CLM.CNVegStateData()
        CLM.cnveg_state_init!(cnveg_state, np, nc)
        cnveg_state.leafcn_offset_patch[1] = 25.0
        cnveg_state.plantCN_patch[1] = plantCN
        cnveg_state.onset_flag_patch[1] = 0.0
        cnveg_state.offset_flag_patch[1] = 0.0
        cnveg_state.c_allometry_patch[1] = c_allometry
        cnveg_state.n_allometry_patch[1] = n_allometry

        # CN Veg carbon state
        cnveg_cs = CLM.CNVegCarbonStateData()
        CLM.cnveg_carbon_state_init!(cnveg_cs, np, nc, ng)
        cnveg_cs.leafc_patch[1] = leafc
        cnveg_cs.frootc_patch[1] = frootc
        cnveg_cs.livestemc_patch[1] = 200.0
        cnveg_cs.livecrootc_patch[1] = 100.0
        cnveg_cs.leafc_storage_patch[1] = 10.0
        cnveg_cs.leafc_storage_xfer_acc_patch[1] = 0.0
        cnveg_cs.storage_cdemand_patch[1] = 0.0

        # CN Veg carbon flux
        cnveg_cf = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cnveg_cf, np, nc, ng)
        cnveg_cf.availc_patch[1] = availc
        cnveg_cf.leafc_to_litter_fun_patch[1] = 0.001
        cnveg_cf.leafc_storage_to_xfer_patch[1] = 0.0

        # CN Veg nitrogen state
        cnveg_ns = CLM.CNVegNitrogenStateData()
        CLM.cnveg_nitrogen_state_init!(cnveg_ns, np, nc, ng)
        cnveg_ns.leafn_patch[1] = leafn
        cnveg_ns.frootn_patch[1] = frootn
        cnveg_ns.livestemn_patch[1] = livestemn
        cnveg_ns.livecrootn_patch[1] = livecrootn
        cnveg_ns.retransn_patch[1] = 1.0
        cnveg_ns.leafn_storage_patch[1] = 0.5
        cnveg_ns.leafn_storage_xfer_acc_patch[1] = 0.0
        cnveg_ns.storage_ndemand_patch[1] = 0.0

        # CN Veg nitrogen flux
        cnveg_nf = CLM.CNVegNitrogenFluxData()
        CLM.cnveg_nitrogen_flux_init!(cnveg_nf, np, nc, ng;
                                      nlevdecomp_full=nlevdecomp_full)
        cnveg_nf.plant_ndemand_patch[1] = 0.0001
        cnveg_nf.leafn_storage_to_xfer_patch[1] = 0.0
        cnveg_nf.plant_ndemand_retrans_patch[1] = 0.0

        # Soil BGC nitrogen flux
        ndecomp_pools = 1
        ndecomp_cascade_transitions = 1
        soilbgc_nf = CLM.SoilBiogeochemNitrogenFluxData()
        CLM.soil_bgc_nitrogen_flux_init!(soilbgc_nf, nc, nlevdecomp_full,
                                         ndecomp_pools, ndecomp_cascade_transitions)
        for j in 1:nlevdecomp
            soilbgc_nf.smin_no3_to_plant_vr_col[1, j] = smin_no3
            soilbgc_nf.smin_nh4_to_plant_vr_col[1, j] = smin_nh4
        end

        # Soil BGC carbon flux
        soilbgc_cf = CLM.SoilBiogeochemCarbonFluxData()
        CLM.soil_bgc_carbon_flux_init!(soilbgc_cf, nc, nlevdecomp_full,
                                       ndecomp_pools, ndecomp_cascade_transitions)

        # Soil BGC nitrogen state
        soilbgc_ns = CLM.SoilBiogeochemNitrogenStateData()
        CLM.soil_bgc_nitrogen_state_init!(soilbgc_ns, nc, ng, nlevdecomp_full,
                                          ndecomp_pools)
        for j in 1:nlevdecomp
            soilbgc_ns.smin_no3_vr_col[1, j] = smin_no3
            soilbgc_ns.smin_nh4_vr_col[1, j] = smin_nh4
        end

        dzsoi_decomp_vals = fill(0.1, nlevdecomp)

        mask_soilp = BitVector([true])
        mask_soilc = BitVector([true])
        bounds_p = 1:np
        bounds_c = 1:nc

        return (pftcon=pftcon, fun_params=fun_params, patch=patch,
                waterstate=waterstate, temperature=temperature,
                soilstate=soilstate, canopystate=canopystate,
                cnveg_state=cnveg_state, cnveg_cs=cnveg_cs, cnveg_cf=cnveg_cf,
                cnveg_ns=cnveg_ns, cnveg_nf=cnveg_nf,
                soilbgc_nf=soilbgc_nf, soilbgc_cf=soilbgc_cf, soilbgc_ns=soilbgc_ns,
                dzsoi_decomp_vals=dzsoi_decomp_vals,
                mask_soilp=mask_soilp, mask_soilc=mask_soilc,
                bounds_p=bounds_p, bounds_c=bounds_c,
                dt=dt, nlevdecomp=nlevdecomp)
    end

    # =====================================================================
    # Test: cnfun! produces valid output fluxes
    # =====================================================================
    @testset "cnfun! basic functionality" begin
        d = make_fun_test_data()

        CLM.cnfun!(d.mask_soilp, d.mask_soilc, d.bounds_p, d.bounds_c,
                   d.fun_params, d.pftcon, d.patch,
                   d.waterstate, d.temperature, d.soilstate,
                   d.cnveg_state, d.cnveg_cs, d.cnveg_cf,
                   d.cnveg_ns, d.cnveg_nf, d.soilbgc_nf,
                   d.soilbgc_cf, d.canopystate, d.soilbgc_ns;
                   dt=d.dt, nlevdecomp=d.nlevdecomp,
                   dzsoi_decomp_vals=d.dzsoi_decomp_vals,
                   use_flexiblecn=false, use_matrixcn=false)

        # All output N fluxes should be finite
        @test isfinite(d.cnveg_nf.Nuptake_patch[1])
        @test isfinite(d.cnveg_nf.Nfix_patch[1])
        @test isfinite(d.cnveg_nf.Nactive_patch[1])
        @test isfinite(d.cnveg_nf.Nnonmyc_patch[1])
        @test isfinite(d.cnveg_nf.Npassive_patch[1])
        @test isfinite(d.cnveg_nf.sminn_to_plant_fun_patch[1])
        @test isfinite(d.cnveg_nf.retransn_to_npool_patch[1])
        @test isfinite(d.cnveg_nf.free_retransn_to_npool_patch[1])

        # C fluxes should be finite
        @test isfinite(d.cnveg_cf.npp_Nuptake_patch[1])
        @test isfinite(d.cnveg_cf.npp_Nfix_patch[1])
        @test isfinite(d.cnveg_cf.npp_Nactive_patch[1])
        @test isfinite(d.cnveg_cf.npp_Nretrans_patch[1])
        @test isfinite(d.cnveg_cf.soilc_change_patch[1])
        @test isfinite(d.cnveg_cf.npp_growth_patch[1])

        # N uptake should be non-negative
        @test d.cnveg_nf.Nuptake_patch[1] >= 0.0
        @test d.cnveg_nf.Nfix_patch[1] >= 0.0

        # Soil N extraction should be non-negative
        @test d.cnveg_nf.sminn_to_plant_fun_patch[1] >= 0.0

        # Column-level aggregation should be consistent
        @test isfinite(d.soilbgc_cf.soilc_change_col[1])
        @test d.soilbgc_cf.soilc_change_col[1] ≈ d.cnveg_cf.soilc_change_patch[1] * d.patch.wtcol[1]
    end

    # =====================================================================
    # Test: cnfun! with zero available C → zero uptake
    # =====================================================================
    @testset "cnfun! zero available C" begin
        d = make_fun_test_data(availc=0.0)

        CLM.cnfun!(d.mask_soilp, d.mask_soilc, d.bounds_p, d.bounds_c,
                   d.fun_params, d.pftcon, d.patch,
                   d.waterstate, d.temperature, d.soilstate,
                   d.cnveg_state, d.cnveg_cs, d.cnveg_cf,
                   d.cnveg_ns, d.cnveg_nf, d.soilbgc_nf,
                   d.soilbgc_cf, d.canopystate, d.soilbgc_ns;
                   dt=d.dt, nlevdecomp=d.nlevdecomp,
                   dzsoi_decomp_vals=d.dzsoi_decomp_vals)

        # With zero available C, active uptake should be zero or very small
        @test d.cnveg_cf.npp_Nfix_patch[1] ≈ 0.0 atol=1e-15
        @test d.cnveg_nf.Nfix_patch[1] ≈ 0.0 atol=1e-15
    end

    # =====================================================================
    # Test: cnfun! with zero soil N → zero soil extraction
    # =====================================================================
    @testset "cnfun! zero soil N" begin
        d = make_fun_test_data(smin_no3=0.0, smin_nh4=0.0)

        CLM.cnfun!(d.mask_soilp, d.mask_soilc, d.bounds_p, d.bounds_c,
                   d.fun_params, d.pftcon, d.patch,
                   d.waterstate, d.temperature, d.soilstate,
                   d.cnveg_state, d.cnveg_cs, d.cnveg_cf,
                   d.cnveg_ns, d.cnveg_nf, d.soilbgc_nf,
                   d.soilbgc_cf, d.canopystate, d.soilbgc_ns;
                   dt=d.dt, nlevdecomp=d.nlevdecomp,
                   dzsoi_decomp_vals=d.dzsoi_decomp_vals)

        # Fluxes should still be finite
        @test isfinite(d.cnveg_nf.Nuptake_patch[1])
        @test isfinite(d.cnveg_cf.npp_Nuptake_patch[1])
    end

    # =====================================================================
    # Test: cnfun! with masked-out patch
    # =====================================================================
    @testset "cnfun! masked-out patch" begin
        d = make_fun_test_data()

        # Mask out the only patch
        d.mask_soilp[1] = false

        CLM.cnfun!(d.mask_soilp, d.mask_soilc, d.bounds_p, d.bounds_c,
                   d.fun_params, d.pftcon, d.patch,
                   d.waterstate, d.temperature, d.soilstate,
                   d.cnveg_state, d.cnveg_cs, d.cnveg_cf,
                   d.cnveg_ns, d.cnveg_nf, d.soilbgc_nf,
                   d.soilbgc_cf, d.canopystate, d.soilbgc_ns;
                   dt=d.dt, nlevdecomp=d.nlevdecomp,
                   dzsoi_decomp_vals=d.dzsoi_decomp_vals)

        # Column aggregation should be zero (no active patches)
        @test d.soilbgc_cf.soilc_change_col[1] ≈ 0.0
        @test d.soilbgc_nf.nfix_to_sminn_col[1] ≈ 0.0
    end

    # =====================================================================
    # Test: cnfun! with flexibleCN enabled
    # =====================================================================
    @testset "cnfun! with flexibleCN" begin
        d = make_fun_test_data()

        CLM.cnfun!(d.mask_soilp, d.mask_soilc, d.bounds_p, d.bounds_c,
                   d.fun_params, d.pftcon, d.patch,
                   d.waterstate, d.temperature, d.soilstate,
                   d.cnveg_state, d.cnveg_cs, d.cnveg_cf,
                   d.cnveg_ns, d.cnveg_nf, d.soilbgc_nf,
                   d.soilbgc_cf, d.canopystate, d.soilbgc_ns;
                   dt=d.dt, nlevdecomp=d.nlevdecomp,
                   dzsoi_decomp_vals=d.dzsoi_decomp_vals,
                   use_flexiblecn=true)

        # All fluxes should still be finite
        @test isfinite(d.cnveg_nf.Nuptake_patch[1])
        @test isfinite(d.cnveg_cf.npp_Nuptake_patch[1])
        @test isfinite(d.cnveg_cf.npp_growth_patch[1])
    end

    # =====================================================================
    # Test: cnfun_init!
    # =====================================================================
    @testset "cnfun_init!" begin
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
        np = 1; nc = 1; ng = 1
        pftcon = CLM.PftConFUN(leafcn=[25.0])

        fun_params = CLM.FUNParams(ndays_off=21.0)

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.itype[1] = 1

        cnveg_state = CLM.CNVegStateData()
        CLM.cnveg_state_init!(cnveg_state, np, nc)
        cnveg_state.leafcn_offset_patch[1] = 0.0

        cnveg_cs = CLM.CNVegCarbonStateData()
        CLM.cnveg_carbon_state_init!(cnveg_cs, np, nc, ng)
        cnveg_cs.storage_cdemand_patch[1] = 99.0
        cnveg_cs.leafc_storage_xfer_acc_patch[1] = 99.0

        cnveg_ns = CLM.CNVegNitrogenStateData()
        CLM.cnveg_nitrogen_state_init!(cnveg_ns, np, nc, ng)
        cnveg_ns.storage_ndemand_patch[1] = 99.0
        cnveg_ns.leafn_storage_xfer_acc_patch[1] = 99.0

        mask = BitVector([true])
        bounds = 1:np

        # Call with nstep that triggers reset (nstep == nstep_fun)
        dt = 1800.0
        dayspyr = 365.0
        nstep_fun = round(Int, CLM.SECSPDAY * dayspyr / dt)

        CLM.cnfun_init!(mask, bounds, fun_params, pftcon, patch,
                        cnveg_state, cnveg_cs, cnveg_ns;
                        dt=dt, nstep=nstep_fun, dayspyr=dayspyr)

        @test cnveg_state.leafcn_offset_patch[1] ≈ 25.0
        @test cnveg_cs.storage_cdemand_patch[1] ≈ 0.0
        @test cnveg_ns.storage_ndemand_patch[1] ≈ 0.0
        @test cnveg_ns.leafn_storage_xfer_acc_patch[1] ≈ 0.0
        @test cnveg_cs.leafc_storage_xfer_acc_patch[1] ≈ 0.0
    end

    # =====================================================================
    # Test: Module-level constants are correct
    # =====================================================================
    @testset "Module constants" begin
        @test CLM.NSTP == 2
        @test CLM.NCOST6 == 6
        @test CLM.BIG_COST == 1.0e9
        @test CLM.ECM_STEP == 1
        @test CLM.AM_STEP == 2
        @test CLM.PLANTS_ARE_FIXING == 1
        @test CLM.PLANTS_NOT_FIXING == 2
    end

end
