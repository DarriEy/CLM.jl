@testset "matrix-CN veg-C wiring — CROP (grain pools)" begin
    # Crop extension of the veg-C matrix: the 18 natveg C pools + 3 grain (reproductive)
    # pools (IGRAIN/IGRAIN_ST/IGRAIN_XF = 19/20/21). Crop adds ONE phenology transfer
    # (igrain_to_iout = grain harvest) + grain allocation (cpool_to_reproductivec[_storage]).
    # gap/fire stay natveg-only (grain does not gap/burn). The matrix advance must
    # reproduce the sequential update EXACTLY, including the grain pools.
    np = 2; dt = 1800.0
    nveg = CLM.NVEGPOOL_NATVEG + CLM.NVEGPOOL_CROP    # 21
    counts = CLM.veg_matrix_transfer_counts(true)
    ivt = [15, 16]; woody = ones(Float64, 20); npcropmin = 15

    donerpool    = CLM._PHC_DONER_NATVEG
    receiverpool = CLM._PHC_RECEIVER_NATVEG
    fldlist = [:leafc_storage_to_xfer_patch, :leafc_xfer_to_leafc_patch,
        :frootc_storage_to_xfer_patch, :frootc_xfer_to_frootc_patch,
        :livestemc_storage_to_xfer_patch, :livestemc_xfer_to_livestemc_patch,
        :deadstemc_storage_to_xfer_patch, :deadstemc_xfer_to_deadstemc_patch,
        :livecrootc_storage_to_xfer_patch, :livecrootc_xfer_to_livecrootc_patch,
        :deadcrootc_storage_to_xfer_patch, :deadcrootc_xfer_to_deadcrootc_patch,
        :livestemc_to_deadstemc_patch, :livecrootc_to_deadcrootc_patch,
        :leafc_to_litter_patch, :frootc_to_litter_patch, :livestemc_to_litter_patch]
    rates = [3e-6 2e-6 3e-6 2e-6 1e-6 8e-7 6e-7 4e-7 3e-7 2e-7 1e-7 9e-8 5e-7 4e-7 3e-7 2e-7 1e-7;
             2e-6 6e-6 1e-6 5e-6 7e-7 9e-7 5e-7 3e-7 2e-7 4e-7 8e-8 6e-8 6e-7 3e-7 2e-7 3e-7 2e-7]

    poolget(cs, i, p) = CLM._vegc_pool_val(cs, i, p)
    poolfield = Dict(
        CLM.ILEAF=>:leafc_patch, CLM.ILEAF_ST=>:leafc_storage_patch, CLM.ILEAF_XF=>:leafc_xfer_patch,
        CLM.IFROOT=>:frootc_patch, CLM.IFROOT_ST=>:frootc_storage_patch, CLM.IFROOT_XF=>:frootc_xfer_patch,
        CLM.ILIVESTEM=>:livestemc_patch, CLM.ILIVESTEM_ST=>:livestemc_storage_patch, CLM.ILIVESTEM_XF=>:livestemc_xfer_patch,
        CLM.IDEADSTEM=>:deadstemc_patch, CLM.IDEADSTEM_ST=>:deadstemc_storage_patch, CLM.IDEADSTEM_XF=>:deadstemc_xfer_patch,
        CLM.ILIVECROOT=>:livecrootc_patch, CLM.ILIVECROOT_ST=>:livecrootc_storage_patch, CLM.ILIVECROOT_XF=>:livecrootc_xfer_patch,
        CLM.IDEADCROOT=>:deadcrootc_patch, CLM.IDEADCROOT_ST=>:deadcrootc_storage_patch, CLM.IDEADCROOT_XF=>:deadcrootc_xfer_patch)
    # grain pools are 2D [p,1]
    graincol = Dict(CLM.IGRAIN=>:reproductivec_patch, CLM.IGRAIN_ST=>:reproductivec_storage_patch,
                    CLM.IGRAIN_XF=>:reproductivec_xfer_patch)
    pooladd!(cs, i, p, v) = begin
        if haskey(poolfield, i); getfield(cs, poolfield[i])[p] += v
        elseif haskey(graincol, i); getfield(cs, graincol[i])[p, 1] += v; end
    end

    function fresh()
        cs = CLM.CNVegCarbonStateData(); CLM.cnveg_carbon_state_init!(cs, np, 1, 1; use_matrixcn=false, nrepr=1)
        cs.leafc_patch .= [50.,30.];  cs.leafc_storage_patch .= [20.,15.];  cs.leafc_xfer_patch .= [10.,8.]
        cs.frootc_patch .= [40.,25.]; cs.frootc_storage_patch .= [18.,12.]; cs.frootc_xfer_patch .= [9.,7.]
        cs.livestemc_patch .= [60.,45.]; cs.livestemc_storage_patch .= [22.,16.]; cs.livestemc_xfer_patch .= [11.,9.]
        cs.deadstemc_patch .= [100.,80.]; cs.deadstemc_storage_patch .= [5.,4.]; cs.deadstemc_xfer_patch .= [3.,2.5]
        cs.livecrootc_patch .= [35.,28.]; cs.livecrootc_storage_patch .= [12.,9.]; cs.livecrootc_xfer_patch .= [6.,5.]
        cs.deadcrootc_patch .= [70.,55.]; cs.deadcrootc_storage_patch .= [4.,3.]; cs.deadcrootc_xfer_patch .= [2.,1.8]
        cs.reproductivec_patch = reshape([25.,18.], np, 1)
        cs.reproductivec_storage_patch = reshape([6.,4.], np, 1)
        cs.reproductivec_xfer_patch = reshape([3.,2.], np, 1)
        return cs
    end
    function set_fluxes!(cf, cs)
        for (k, f) in enumerate(fldlist)
            length(getfield(cf, f)) == 0 && setfield!(cf, f, zeros(np))
            arr = getfield(cf, f); for p in 1:np; arr[p] = poolget(cs, donerpool[k], p) * rates[p, k]; end
        end
    end
    gm_rate = [1e-6,2e-6,3e-6,1.5e-6,2.5e-6,1e-6,8e-7,6e-7,4e-7,3e-7,2e-7,1e-7,7e-7,5e-7,3e-7,2e-7,1e-7,9e-8]
    set_gm!(cf, cs) = for i in 1:CLM.NVEGPOOL_NATVEG
        f = CLM._GMC_FLUX[i]; length(getfield(cf, f)) == 0 && setfield!(cf, f, zeros(np))
        arr = getfield(cf, f); for p in 1:np; arr[p] = poolget(cs, i, p) * gm_rate[i]; end
    end
    fi_out_rate = [1e-6,1.5e-6,2e-6,8e-7,1e-6,1.2e-6,7e-7,5e-7,3e-7,2e-7,1e-7,9e-8,6e-7,4e-7,2e-7,1e-7,8e-8,6e-8]
    fi_ls_rate = 3e-7; fi_lc_rate = 2e-7
    set_fi!(cf, cs) = begin
        for i in 1:CLM.NVEGPOOL_NATVEG
            f = CLM._FIC_TO_FIRE[i]; length(getfield(cf, f)) == 0 && setfield!(cf, f, zeros(np))
            fl = CLM._FIC_TO_LITTER[i]; length(getfield(cf, fl)) == 0 && setfield!(cf, fl, zeros(np)); fill!(getfield(cf, fl), 0.0)
            arr = getfield(cf, f); for p in 1:np; arr[p] = poolget(cs, i, p)*fi_out_rate[i]; end
        end
        length(cf.m_livestemc_to_deadstemc_fire_patch)==0 && (cf.m_livestemc_to_deadstemc_fire_patch=zeros(np))
        length(cf.m_livecrootc_to_deadcrootc_fire_patch)==0 && (cf.m_livecrootc_to_deadcrootc_fire_patch=zeros(np))
        for p in 1:np
            cf.m_livestemc_to_deadstemc_fire_patch[p]  = poolget(cs, CLM.ILIVESTEM, p)*fi_ls_rate
            cf.m_livecrootc_to_deadcrootc_fire_patch[p] = poolget(cs, CLM.ILIVECROOT, p)*fi_lc_rate
        end
    end
    alloc_flux = Dict(CLM.ILEAF=>3e-5, CLM.ILEAF_ST=>1e-5, CLM.IFROOT=>2.5e-5, CLM.IFROOT_ST=>8e-6,
        CLM.ILIVESTEM=>2e-5, CLM.ILIVESTEM_ST=>6e-6, CLM.IDEADSTEM=>1e-5, CLM.IDEADSTEM_ST=>4e-6,
        CLM.ILIVECROOT=>1.5e-5, CLM.ILIVECROOT_ST=>5e-6, CLM.IDEADCROOT=>8e-6, CLM.IDEADCROOT_ST=>3e-6)
    grain_alloc = Dict(CLM.IGRAIN=>2e-5, CLM.IGRAIN_ST=>7e-6)   # cpool -> grain
    grain_harv_food = 4e-7; grain_harv_seed = 1e-7             # grain -> out (per gC/gC/s)
    set_alloc!(cf, cs) = begin
        for (pool, fld) in CLM._ALLOC_TARGET
            length(getfield(cf, fld)) == 0 && setfield!(cf, fld, zeros(np))
            arr = getfield(cf, fld); for p in 1:np; arr[p] = alloc_flux[pool]; end
        end
        cf.cpool_to_reproductivec_patch = reshape(fill(grain_alloc[CLM.IGRAIN], np), np, 1)
        cf.cpool_to_reproductivec_storage_patch = reshape(fill(grain_alloc[CLM.IGRAIN_ST], np), np, 1)
        cf.repr_grainc_to_food_patch = reshape([poolget(cs, CLM.IGRAIN, p)*grain_harv_food for p in 1:np], np, 1)
        cf.repr_grainc_to_seed_patch = reshape([poolget(cs, CLM.IGRAIN, p)*grain_harv_seed for p in 1:np], np, 1)
    end

    # (a) sequential reference
    cs_seq = fresh()
    cf_seq = CLM.CNVegCarbonFluxData(); CLM.cnveg_carbon_flux_init!(cf_seq, np, 1, 1; use_matrixcn=true)
    set_fluxes!(cf_seq, cs_seq); set_gm!(cf_seq, cs_seq); set_fi!(cf_seq, cs_seq); set_alloc!(cf_seq, cs_seq)
    for p in 1:np, (pool, fld) in CLM._ALLOC_TARGET
        pooladd!(cs_seq, pool, p, getfield(cf_seq, fld)[p] * dt)
    end
    for p in 1:np                                  # grain allocation
        pooladd!(cs_seq, CLM.IGRAIN,    p, cf_seq.cpool_to_reproductivec_patch[p, 1] * dt)
        pooladd!(cs_seq, CLM.IGRAIN_ST, p, cf_seq.cpool_to_reproductivec_storage_patch[p, 1] * dt)
    end
    for p in 1:np, k in 1:17
        flux = getfield(cf_seq, fldlist[k])[p]
        # skip the virtual "outside veg" receiver (== NVEGPOOL_NATVEG+1, which collides
        # with IGRAIN=19); litterfall out-transfers only remove from the donor.
        receiverpool[k] != CLM.NVEGPOOL_NATVEG + 1 && pooladd!(cs_seq, receiverpool[k], p, flux*dt)
        pooladd!(cs_seq, donerpool[k], p, -flux*dt)
    end
    for p in 1:np                                  # grain harvest (grain -> out)
        pooladd!(cs_seq, CLM.IGRAIN, p, -(cf_seq.repr_grainc_to_food_patch[p, 1] + cf_seq.repr_grainc_to_seed_patch[p, 1]) * dt)
    end
    for p in 1:np, i in 1:CLM.NVEGPOOL_NATVEG
        pooladd!(cs_seq, i, p, -getfield(cf_seq, CLM._GMC_FLUX[i])[p] * dt)
    end
    for p in 1:np
        fls = cf_seq.m_livestemc_to_deadstemc_fire_patch[p]; flc = cf_seq.m_livecrootc_to_deadcrootc_fire_patch[p]
        pooladd!(cs_seq, CLM.ILIVESTEM, p, -fls*dt); pooladd!(cs_seq, CLM.IDEADSTEM, p, fls*dt)
        pooladd!(cs_seq, CLM.ILIVECROOT, p, -flc*dt); pooladd!(cs_seq, CLM.IDEADCROOT, p, flc*dt)
        for i in 1:CLM.NVEGPOOL_NATVEG
            pooladd!(cs_seq, i, p, -(getfield(cf_seq, CLM._FIC_TO_FIRE[i])[p]+getfield(cf_seq, CLM._FIC_TO_LITTER[i])[p])*dt)
        end
    end

    # (b) matrix advance
    cs_mat = fresh()
    cf = CLM.CNVegCarbonFluxData(); CLM.cnveg_carbon_flux_init!(cf, np, 1, 1; use_matrixcn=true)
    set_fluxes!(cf, cs_mat); set_gm!(cf, cs_mat); set_fi!(cf, cs_mat); set_alloc!(cf, cs_mat)
    CLM.cn_veg_matrix_c_topology!(cf; use_crop=true, nvegcpool=nveg)
    CLM.cn_veg_matrix_alloc_c!(cf, trues(np), 1:np; use_crop=true)
    CLM.cn_veg_matrix_accumulate_ph_c!(cf, cs_mat, trues(np), 1:np; dt=dt, matrixcheck=false, use_crop=true)
    CLM.cn_veg_matrix_accumulate_gm_c!(cf, cs_mat, trues(np), 1:np; dt=dt, matrixcheck=false)
    CLM.cn_veg_matrix_accumulate_fi_c!(cf, cs_mat, trues(np), 1:np; dt=dt, matrixcheck=false)
    CLM.cn_veg_matrix_solve_c!(cs_mat, cf; mask_soilp=trues(np), bounds_patch=1:np,
        ivt=ivt, woody=woody, npcropmin=npcropmin, nvegcpool=nveg, counts=counts, dt=dt, num_actfirep=np)

    for p in 1:np, i in 1:nveg
        @test isapprox(poolget(cs_mat, i, p), poolget(cs_seq, i, p); atol=1e-9, rtol=1e-10)
    end
    @test cf.igrain_to_iout_ph == 18
    @test length(cf.matrix_phtransfer_doner_patch) == counts.ncphtrans   # 18
    # grain actually changed (allocation + harvest, nonzero)
    @test !(cs_mat.reproductivec_patch[1, 1] ≈ 25.0)
end
