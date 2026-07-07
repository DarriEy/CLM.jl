@testset "matrix-CN veg-C wiring (accumulate from live flux fields)" begin
    # Phase-1 matrix-CN wiring: cn_veg_matrix_c_topology! + accumulate_ph_c! read
    # the per-flux phenology transfer fields (as computed by phenology) and record
    # them into the matrix; cn_veg_matrix_solve_c! then advances the pools. At
    # matrixcheck=false this must reproduce the sequential c_state_update EXACTLY.
    np = 3; dt = 1800.0
    nveg = CLM.NVEGPOOL_NATVEG
    counts = CLM.veg_matrix_transfer_counts(false)
    ivt = [1, 2, 3]; woody = ones(Float64, 5); npcropmin = 15
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

    rates = [
     5e-6 4e-6 3e-6 2e-6 1e-6 8e-7 6e-7 4e-7 3e-7 2e-7 1e-7 9e-8 5e-7 4e-7 3e-7 2e-7 0.0;
     2e-6 6e-6 1e-6 5e-6 7e-7 9e-7 5e-7 3e-7 2e-7 4e-7 8e-8 6e-8 6e-7 3e-7 2e-7 3e-7 0.0;
     1e-6 1e-6 1e-6 1e-6 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 0.0]

    poolfield = Dict(
        CLM.ILEAF=>:leafc_patch, CLM.ILEAF_ST=>:leafc_storage_patch, CLM.ILEAF_XF=>:leafc_xfer_patch,
        CLM.IFROOT=>:frootc_patch, CLM.IFROOT_ST=>:frootc_storage_patch, CLM.IFROOT_XF=>:frootc_xfer_patch,
        CLM.ILIVESTEM=>:livestemc_patch, CLM.ILIVESTEM_ST=>:livestemc_storage_patch, CLM.ILIVESTEM_XF=>:livestemc_xfer_patch,
        CLM.IDEADSTEM=>:deadstemc_patch, CLM.IDEADSTEM_ST=>:deadstemc_storage_patch, CLM.IDEADSTEM_XF=>:deadstemc_xfer_patch,
        CLM.ILIVECROOT=>:livecrootc_patch, CLM.ILIVECROOT_ST=>:livecrootc_storage_patch, CLM.ILIVECROOT_XF=>:livecrootc_xfer_patch,
        CLM.IDEADCROOT=>:deadcrootc_patch, CLM.IDEADCROOT_ST=>:deadcrootc_storage_patch, CLM.IDEADCROOT_XF=>:deadcrootc_xfer_patch)
    poolget(cs, i, p) = CLM._vegc_pool_val(cs, i, p)
    pooladd!(cs, i, p, v) = (haskey(poolfield, i) && (getfield(cs, poolfield[i])[p] += v))

    function fresh()
        cs = CLM.CNVegCarbonStateData(); CLM.cnveg_carbon_state_init!(cs, np, 1, 1; use_matrixcn=false, nrepr=1)
        cs.leafc_patch .= [50.,30.,10.];  cs.leafc_storage_patch .= [20.,15.,5.];  cs.leafc_xfer_patch .= [10.,8.,3.]
        cs.frootc_patch .= [40.,25.,12.]; cs.frootc_storage_patch .= [18.,12.,6.]; cs.frootc_xfer_patch .= [9.,7.,4.]
        cs.livestemc_patch .= [60.,45.,20.]; cs.livestemc_storage_patch .= [22.,16.,8.]; cs.livestemc_xfer_patch .= [11.,9.,5.]
        cs.deadstemc_patch .= [100.,80.,40.]; cs.deadstemc_storage_patch .= [5.,4.,2.]; cs.deadstemc_xfer_patch .= [3.,2.5,1.]
        cs.livecrootc_patch .= [35.,28.,14.]; cs.livecrootc_storage_patch .= [12.,9.,5.]; cs.livecrootc_xfer_patch .= [6.,5.,2.5]
        cs.deadcrootc_patch .= [70.,55.,28.]; cs.deadcrootc_storage_patch .= [4.,3.,1.5]; cs.deadcrootc_xfer_patch .= [2.,1.8,.8]
        return cs
    end
    function set_fluxes!(cf, cs)
        for (k, f) in enumerate(fldlist)
            length(getfield(cf, f)) == 0 && setfield!(cf, f, zeros(np))
            arr = getfield(cf, f)
            for p in 1:np; arr[p] = poolget(cs, donerpool[k], p) * rates[p, k]; end
        end
    end

    # (a) sequential reference
    cs_seq = fresh()
    cf_seq = CLM.CNVegCarbonFluxData(); CLM.cnveg_carbon_flux_init!(cf_seq, np, 1, 1; use_matrixcn=true)
    set_fluxes!(cf_seq, cs_seq)
    for p in 1:np, k in 1:17
        flux = getfield(cf_seq, fldlist[k])[p]
        pooladd!(cs_seq, receiverpool[k], p,  flux*dt)
        pooladd!(cs_seq, donerpool[k],    p, -flux*dt)
    end

    # (b) matrix advance via the wiring under test
    cs_mat = fresh()
    cf = CLM.CNVegCarbonFluxData(); CLM.cnveg_carbon_flux_init!(cf, np, 1, 1; use_matrixcn=true)
    set_fluxes!(cf, cs_mat)
    CLM.cn_veg_matrix_c_topology!(cf; use_crop=false, nvegcpool=nveg)
    cf.matrix_alloc_patch = zeros(np, nveg); cf.matrix_Cinput_patch = zeros(np)
    cf.matrix_gmtransfer_doner_patch = fill(CLM.ILEAF, counts.ncgmtrans)
    cf.matrix_gmtransfer_receiver_patch = fill(CLM.ILEAF, counts.ncgmtrans)
    cf.matrix_gmtransfer_patch = zeros(np, counts.ncgmtrans); cf.matrix_gmturnover_patch = zeros(np, nveg)
    cf.matrix_fitransfer_doner_patch = fill(CLM.ILEAF, counts.ncfitrans)
    cf.matrix_fitransfer_receiver_patch = fill(CLM.ILEAF, counts.ncfitrans)
    cf.matrix_fitransfer_patch = zeros(np, counts.ncfitrans); cf.matrix_fiturnover_patch = zeros(np, nveg)

    CLM.cn_veg_matrix_accumulate_ph_c!(cf, cs_mat, trues(np), 1:np; dt=dt, matrixcheck=false)
    CLM.cn_veg_matrix_solve_c!(cs_mat, cf; mask_soilp=trues(np), bounds_patch=1:np,
        ivt=ivt, woody=woody, npcropmin=npcropmin, nvegcpool=nveg, counts=counts, dt=dt, num_actfirep=0)

    # matrix advance must reproduce the sequential update exactly
    for p in 1:np, i in 1:nveg
        @test isapprox(poolget(cs_mat, i, p), poolget(cs_seq, i, p); atol=1e-9, rtol=1e-10)
    end

    # topology sanity
    @test cf.ileafst_to_ileafxf_ph == 1
    @test cf.ilivestem_to_iout_ph == 17
    @test length(cf.matrix_phtransfer_doner_patch) == counts.ncphtrans
end
