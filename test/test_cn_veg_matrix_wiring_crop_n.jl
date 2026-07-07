@testset "matrix-CN veg-N wiring — CROP (grain + retransn=22)" begin
    # Crop veg-N matrix: 18 natveg N pools + 3 grain (19..21) + retransn shifted to pool
    # 22 (nvegnpool). Crop adds 3 phenology transfers (2 retransn→grain supply + grain
    # harvest) and 1 more out-transfer. gap/fire stay natveg+retransn (grain excluded).
    # The matrix advance must reproduce an independent sequential N update EXACTLY.
    np = 2; dt = 1800.0
    nveg = CLM.NVEGPOOL_NATVEG + CLM.NVEGPOOL_CROP + 1   # 22 (incl grain + retransn)
    IRET = nveg                                          # retransn = pool 22
    counts = CLM.veg_matrix_transfer_counts(true)
    ivt = [15, 16]; npcropmin = 15

    npget(ns, i, p) = CLM._vegn_pool_val(ns, i, p, IRET)
    gmpool(i) = i <= CLM.NVEGPOOL_NATVEG ? i : IRET      # gm/fin index 19 → retransn
    nfield = Dict(
        CLM.ILEAF=>:leafn_patch, CLM.ILEAF_ST=>:leafn_storage_patch, CLM.ILEAF_XF=>:leafn_xfer_patch,
        CLM.IFROOT=>:frootn_patch, CLM.IFROOT_ST=>:frootn_storage_patch, CLM.IFROOT_XF=>:frootn_xfer_patch,
        CLM.ILIVESTEM=>:livestemn_patch, CLM.ILIVESTEM_ST=>:livestemn_storage_patch, CLM.ILIVESTEM_XF=>:livestemn_xfer_patch,
        CLM.IDEADSTEM=>:deadstemn_patch, CLM.IDEADSTEM_ST=>:deadstemn_storage_patch, CLM.IDEADSTEM_XF=>:deadstemn_xfer_patch,
        CLM.ILIVECROOT=>:livecrootn_patch, CLM.ILIVECROOT_ST=>:livecrootn_storage_patch, CLM.ILIVECROOT_XF=>:livecrootn_xfer_patch,
        CLM.IDEADCROOT=>:deadcrootn_patch, CLM.IDEADCROOT_ST=>:deadcrootn_storage_patch, CLM.IDEADCROOT_XF=>:deadcrootn_xfer_patch,
        IRET=>:retransn_patch)
    graincol = Dict(CLM.IGRAIN=>:reproductiven_patch, CLM.IGRAIN_ST=>:reproductiven_storage_patch,
                    CLM.IGRAIN_XF=>:reproductiven_xfer_patch)
    nadd!(ns, i, p, v) = begin
        if haskey(nfield, i); getfield(ns, nfield[i])[p] += v
        elseif haskey(graincol, i); getfield(ns, graincol[i])[p, 1] += v; end
    end

    function fresh()
        ns = CLM.CNVegNitrogenStateData(); CLM.cnveg_nitrogen_state_init!(ns, np, 1, 1; use_matrixcn=false, nrepr=1)
        ns.leafn_patch .= [5.,3.];   ns.leafn_storage_patch .= [2.,1.5];  ns.leafn_xfer_patch .= [1.,.8]
        ns.frootn_patch .= [4.,2.5]; ns.frootn_storage_patch .= [1.8,1.2]; ns.frootn_xfer_patch .= [.9,.7]
        ns.livestemn_patch .= [6.,4.5]; ns.livestemn_storage_patch .= [2.2,1.6]; ns.livestemn_xfer_patch .= [1.1,.9]
        ns.deadstemn_patch .= [10.,8.]; ns.deadstemn_storage_patch .= [.5,.4]; ns.deadstemn_xfer_patch .= [.3,.25]
        ns.livecrootn_patch .= [3.5,2.8]; ns.livecrootn_storage_patch .= [1.2,.9]; ns.livecrootn_xfer_patch .= [.6,.5]
        ns.deadcrootn_patch .= [7.,5.5]; ns.deadcrootn_storage_patch .= [.4,.3]; ns.deadcrootn_xfer_patch .= [.2,.18]
        ns.retransn_patch .= [2.,1.5]
        ns.reproductiven_patch = reshape([1.6,1.1], np, 1)
        ns.reproductiven_storage_patch = reshape([.5,.35], np, 1)
        ns.reproductiven_xfer_patch = reshape([.2,.15], np, 1)
        return ns
    end

    ph_specs = [
        (:leafn_storage_to_xfer_patch, CLM.ILEAF_ST, CLM.ILEAF_XF, 5e-6),
        (:leafn_xfer_to_leafn_patch, CLM.ILEAF_XF, CLM.ILEAF, 4e-6),
        (:frootn_storage_to_xfer_patch, CLM.IFROOT_ST, CLM.IFROOT_XF, 3e-6),
        (:frootn_xfer_to_frootn_patch, CLM.IFROOT_XF, CLM.IFROOT, 2e-6),
        (:livestemn_storage_to_xfer_patch, CLM.ILIVESTEM_ST, CLM.ILIVESTEM_XF, 1e-6),
        (:livestemn_xfer_to_livestemn_patch, CLM.ILIVESTEM_XF, CLM.ILIVESTEM, 8e-7),
        (:deadstemn_storage_to_xfer_patch, CLM.IDEADSTEM_ST, CLM.IDEADSTEM_XF, 6e-7),
        (:deadstemn_xfer_to_deadstemn_patch, CLM.IDEADSTEM_XF, CLM.IDEADSTEM, 4e-7),
        (:livecrootn_storage_to_xfer_patch, CLM.ILIVECROOT_ST, CLM.ILIVECROOT_XF, 3e-7),
        (:livecrootn_xfer_to_livecrootn_patch, CLM.ILIVECROOT_XF, CLM.ILIVECROOT, 2e-7),
        (:deadcrootn_storage_to_xfer_patch, CLM.IDEADCROOT_ST, CLM.IDEADCROOT_XF, 1e-7),
        (:deadcrootn_xfer_to_deadcrootn_patch, CLM.IDEADCROOT_XF, CLM.IDEADCROOT, 9e-8),
        (:livestemn_to_deadstemn_patch, CLM.ILIVESTEM, CLM.IDEADSTEM, 5e-7),
        (:livecrootn_to_deadcrootn_patch, CLM.ILIVECROOT, CLM.IDEADCROOT, 4e-7),
        (:leafn_to_retransn_patch, CLM.ILEAF, IRET, 7e-7),
        (:frootn_to_retransn_patch, CLM.IFROOT, IRET, 6e-7),
        (:livestemn_to_retransn_patch, CLM.ILIVESTEM, IRET, 3e-7),
        (:livecrootn_to_retransn_patch, CLM.ILIVECROOT, IRET, 2e-7),
        (:leafn_to_litter_patch, CLM.ILEAF, 0, 3e-7),
        (:frootn_to_litter_patch, CLM.IFROOT, 0, 2e-7),
        (:livestemn_to_litter_patch, CLM.ILIVESTEM, 0, 1e-7)]
    gm_rate = [1e-6,2e-6,3e-6,1.5e-6,2.5e-6,1e-6,8e-7,6e-7,4e-7,3e-7,2e-7,1e-7,7e-7,5e-7,3e-7,2e-7,1e-7,9e-8,5e-7]
    fi_out_rate = [1e-6,1.5e-6,2e-6,8e-7,1e-6,1.2e-6,7e-7,5e-7,3e-7,2e-7,1e-7,9e-8,6e-7,4e-7,2e-7,1e-7,8e-8,6e-8,4e-7]
    fi_ls_rate = 3e-7; fi_lc_rate = 2e-7
    alloc_flux = Dict(CLM.ILEAF=>3e-5, CLM.ILEAF_ST=>1e-5, CLM.IFROOT=>2.5e-5, CLM.IFROOT_ST=>8e-6,
        CLM.ILIVESTEM=>2e-5, CLM.ILIVESTEM_ST=>6e-6, CLM.IDEADSTEM=>1e-5, CLM.IDEADSTEM_ST=>4e-6,
        CLM.ILIVECROOT=>1.5e-5, CLM.ILIVECROOT_ST=>5e-6, CLM.IDEADCROOT=>8e-6, CLM.IDEADCROOT_ST=>3e-6)
    grain_alloc = Dict(CLM.IGRAIN=>1.2e-5, CLM.IGRAIN_ST=>4e-6)
    grain_harv_food = 3e-7; grain_harv_seed = 1e-7
    r2n_val = [2e-5, 1.5e-5]

    setv!(nf, f, v) = (length(getfield(nf, f))==0 && setfield!(nf, f, zeros(np)); getfield(nf, f) .= v)
    function set_all!(nf, ns)
        for (f, d, _, r) in ph_specs
            setv!(nf, f, 0.0); arr=getfield(nf,f); for p in 1:np; arr[p]=npget(ns,d,p)*r; end
        end
        for i in 1:19
            f = CLM._GMN_FLUX[i]; setv!(nf, f, 0.0); arr=getfield(nf,f); for p in 1:np; arr[p]=npget(ns,gmpool(i),p)*gm_rate[i]; end
            ff = CLM._FIN_TO_FIRE[i]; setv!(nf, ff, 0.0); arr2=getfield(nf,ff); for p in 1:np; arr2[p]=npget(ns,gmpool(i),p)*fi_out_rate[i]; end
            setv!(nf, CLM._FIN_TO_LITTER[i], 0.0)
        end
        setv!(nf, :m_livestemn_to_deadstemn_fire_patch, 0.0); setv!(nf, :m_livecrootn_to_deadcrootn_fire_patch, 0.0)
        for p in 1:np
            nf.m_livestemn_to_deadstemn_fire_patch[p] = npget(ns, CLM.ILIVESTEM, p)*fi_ls_rate
            nf.m_livecrootn_to_deadcrootn_fire_patch[p] = npget(ns, CLM.ILIVECROOT, p)*fi_lc_rate
        end
        for (pool, fld) in CLM._NALLOC_TARGET; setv!(nf, fld, alloc_flux[pool]); end
        nf.npool_to_reproductiven_patch = reshape(fill(grain_alloc[CLM.IGRAIN], np), np, 1)
        nf.npool_to_reproductiven_storage_patch = reshape(fill(grain_alloc[CLM.IGRAIN_ST], np), np, 1)
        nf.repr_grainn_to_food_patch = reshape([npget(ns, CLM.IGRAIN, p)*grain_harv_food for p in 1:np], np, 1)
        nf.repr_grainn_to_seed_patch = reshape([npget(ns, CLM.IGRAIN, p)*grain_harv_seed for p in 1:np], np, 1)
        setv!(nf, :retransn_to_npool_patch, 0.0); nf.retransn_to_npool_patch .= r2n_val
    end

    # ============ (a) independent sequential reference ============
    ns_seq = fresh()
    nf_seq = CLM.CNVegNitrogenFluxData(); CLM.cnveg_nitrogen_flux_init!(nf_seq, np, 1, 1; use_matrixcn=true, nvegnpool=nveg)
    set_all!(nf_seq, ns_seq)
    for p in 1:np
        for (pool, fld) in CLM._NALLOC_TARGET; nadd!(ns_seq, pool, p, getfield(nf_seq, fld)[p]*dt); end
        nadd!(ns_seq, CLM.IGRAIN, p, nf_seq.npool_to_reproductiven_patch[p, 1]*dt)          # grain alloc net
        nadd!(ns_seq, CLM.IGRAIN_ST, p, nf_seq.npool_to_reproductiven_storage_patch[p, 1]*dt)
        nadd!(ns_seq, IRET, p, -nf_seq.retransn_to_npool_patch[p]*dt)                        # retransn depletion
        nadd!(ns_seq, CLM.IGRAIN, p, -(nf_seq.repr_grainn_to_food_patch[p, 1]+nf_seq.repr_grainn_to_seed_patch[p, 1])*dt)  # grain harvest
        for (f, d, r, _) in ph_specs
            flux = getfield(nf_seq, f)[p]; nadd!(ns_seq, d, p, -flux*dt); r == 0 || nadd!(ns_seq, r, p, flux*dt)
        end
        for i in 1:19; nadd!(ns_seq, gmpool(i), p, -getfield(nf_seq, CLM._GMN_FLUX[i])[p]*dt); end
        fls = nf_seq.m_livestemn_to_deadstemn_fire_patch[p]; flc = nf_seq.m_livecrootn_to_deadcrootn_fire_patch[p]
        nadd!(ns_seq, CLM.ILIVESTEM, p, -fls*dt); nadd!(ns_seq, CLM.IDEADSTEM, p, fls*dt)
        nadd!(ns_seq, CLM.ILIVECROOT, p, -flc*dt); nadd!(ns_seq, CLM.IDEADCROOT, p, flc*dt)
        for i in 1:19; nadd!(ns_seq, gmpool(i), p, -(getfield(nf_seq, CLM._FIN_TO_FIRE[i])[p]+getfield(nf_seq, CLM._FIN_TO_LITTER[i])[p])*dt); end
    end

    # ============ (b) matrix advance ============
    ns_mat = fresh()
    nf = CLM.CNVegNitrogenFluxData(); CLM.cnveg_nitrogen_flux_init!(nf, np, 1, 1; use_matrixcn=true, nvegnpool=nveg)
    set_all!(nf, ns_mat)
    CLM.cn_veg_matrix_n_topology!(nf; use_crop=true, nvegnpool=nveg)
    CLM.cn_veg_matrix_alloc_n!(nf, trues(np), 1:np; use_crop=true)
    CLM.cn_veg_matrix_accumulate_phn!(nf, ns_mat, trues(np), 1:np; dt=dt, matrixcheck=false, use_crop=true, nvegnpool=nveg)
    CLM.cn_veg_matrix_accumulate_gmn!(nf, ns_mat, trues(np), 1:np; dt=dt, matrixcheck=false, nvegnpool=nveg)
    CLM.cn_veg_matrix_accumulate_fin!(nf, ns_mat, trues(np), 1:np; dt=dt, matrixcheck=false, nvegnpool=nveg)
    CLM.cn_veg_matrix_solve_n!(ns_mat, nf; mask_soilp=trues(np), bounds_patch=1:np,
        ivt=ivt, npcropmin=npcropmin, nvegnpool=nveg, counts=counts, dt=dt, num_actfirep=np)

    for p in 1:np, i in 1:nveg
        @test isapprox(npget(ns_mat, i, p), npget(ns_seq, i, p); atol=1e-11, rtol=1e-9)
    end
    @test nf.iretransn_to_igrain_ph == 31
    @test nf.iretransn_to_igrainst_ph == 32
    @test nf.igrain_to_iout_ph == 36
    @test nf.iretransn_to_iout_ph == 37
    @test length(nf.matrix_nphtransfer_doner_patch) == counts.nnphtrans   # 37
    @test !(ns_mat.reproductiven_patch[1, 1] ≈ 1.6)   # grain changed
end
