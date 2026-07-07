@testset "matrix-CN veg-N wiring (accumulate from live flux fields)" begin
    # Phase-2 matrix-CN veg-N wiring: cn_veg_matrix_n_topology! + accumulate_{phn,
    # gmn,fin}! + alloc_n! read the live N flux fields and record them into the
    # matrix; cn_veg_matrix_solve_n! then advances all 19 pools (18 C-parallel +
    # retransn). At matrixcheck=false this must reproduce a physically-independent
    # sequential N update EXACTLY (incl. the retranslocation accounting).
    np = 3; dt = 1800.0
    nveg = CLM.IRETRANSN_NATVEG          # 19 (INCLUDES retransn) — solve convention
    counts = CLM.veg_matrix_transfer_counts(false)
    ivt = [1, 2, 3]; npcropmin = 15
    IRET = CLM.IRETRANSN_NATVEG

    npget(ns, i, p) = CLM._vegn_pool_val(ns, i, p)
    # pool index -> N state field for mutation (incl. retransn)
    nfield = Dict(
        CLM.ILEAF=>:leafn_patch, CLM.ILEAF_ST=>:leafn_storage_patch, CLM.ILEAF_XF=>:leafn_xfer_patch,
        CLM.IFROOT=>:frootn_patch, CLM.IFROOT_ST=>:frootn_storage_patch, CLM.IFROOT_XF=>:frootn_xfer_patch,
        CLM.ILIVESTEM=>:livestemn_patch, CLM.ILIVESTEM_ST=>:livestemn_storage_patch, CLM.ILIVESTEM_XF=>:livestemn_xfer_patch,
        CLM.IDEADSTEM=>:deadstemn_patch, CLM.IDEADSTEM_ST=>:deadstemn_storage_patch, CLM.IDEADSTEM_XF=>:deadstemn_xfer_patch,
        CLM.ILIVECROOT=>:livecrootn_patch, CLM.ILIVECROOT_ST=>:livecrootn_storage_patch, CLM.ILIVECROOT_XF=>:livecrootn_xfer_patch,
        CLM.IDEADCROOT=>:deadcrootn_patch, CLM.IDEADCROOT_ST=>:deadcrootn_storage_patch, CLM.IDEADCROOT_XF=>:deadcrootn_xfer_patch,
        IRET=>:retransn_patch)
    nadd!(ns, i, p, v) = (getfield(ns, nfield[i])[p] += v)

    function fresh()
        ns = CLM.CNVegNitrogenStateData(); CLM.cnveg_nitrogen_state_init!(ns, np, 1, 1; use_matrixcn=false, nrepr=1)
        ns.leafn_patch .= [5.,3.,1.];   ns.leafn_storage_patch .= [2.,1.5,.5];  ns.leafn_xfer_patch .= [1.,.8,.3]
        ns.frootn_patch .= [4.,2.5,1.2]; ns.frootn_storage_patch .= [1.8,1.2,.6]; ns.frootn_xfer_patch .= [.9,.7,.4]
        ns.livestemn_patch .= [6.,4.5,2.]; ns.livestemn_storage_patch .= [2.2,1.6,.8]; ns.livestemn_xfer_patch .= [1.1,.9,.5]
        ns.deadstemn_patch .= [10.,8.,4.]; ns.deadstemn_storage_patch .= [.5,.4,.2]; ns.deadstemn_xfer_patch .= [.3,.25,.1]
        ns.livecrootn_patch .= [3.5,2.8,1.4]; ns.livecrootn_storage_patch .= [1.2,.9,.5]; ns.livecrootn_xfer_patch .= [.6,.5,.25]
        ns.deadcrootn_patch .= [7.,5.5,2.8]; ns.deadcrootn_storage_patch .= [.4,.3,.15]; ns.deadcrootn_xfer_patch .= [.2,.18,.08]
        ns.retransn_patch .= [2.,1.5,.8]
        return ns
    end

    # ---- phenology fluxes: 21 (12 storage/xfer + 2 livewood + 4 resorption + 3 litterfall)
    # each as (field, donor pool, receiver pool, rate); flux = donor*rate.
    ph_specs = [
        (:leafn_storage_to_xfer_patch,      CLM.ILEAF_ST,     CLM.ILEAF_XF,     5e-6),
        (:leafn_xfer_to_leafn_patch,        CLM.ILEAF_XF,     CLM.ILEAF,        4e-6),
        (:frootn_storage_to_xfer_patch,     CLM.IFROOT_ST,    CLM.IFROOT_XF,    3e-6),
        (:frootn_xfer_to_frootn_patch,      CLM.IFROOT_XF,    CLM.IFROOT,       2e-6),
        (:livestemn_storage_to_xfer_patch,  CLM.ILIVESTEM_ST, CLM.ILIVESTEM_XF, 1e-6),
        (:livestemn_xfer_to_livestemn_patch, CLM.ILIVESTEM_XF, CLM.ILIVESTEM,   8e-7),
        (:deadstemn_storage_to_xfer_patch,  CLM.IDEADSTEM_ST, CLM.IDEADSTEM_XF, 6e-7),
        (:deadstemn_xfer_to_deadstemn_patch, CLM.IDEADSTEM_XF, CLM.IDEADSTEM,   4e-7),
        (:livecrootn_storage_to_xfer_patch, CLM.ILIVECROOT_ST, CLM.ILIVECROOT_XF, 3e-7),
        (:livecrootn_xfer_to_livecrootn_patch, CLM.ILIVECROOT_XF, CLM.ILIVECROOT, 2e-7),
        (:deadcrootn_storage_to_xfer_patch, CLM.IDEADCROOT_ST, CLM.IDEADCROOT_XF, 1e-7),
        (:deadcrootn_xfer_to_deadcrootn_patch, CLM.IDEADCROOT_XF, CLM.IDEADCROOT, 9e-8),
        (:livestemn_to_deadstemn_patch,     CLM.ILIVESTEM,    CLM.IDEADSTEM,    5e-7),
        (:livecrootn_to_deadcrootn_patch,   CLM.ILIVECROOT,   CLM.IDEADCROOT,   4e-7),
        (:leafn_to_retransn_patch,          CLM.ILEAF,        IRET,             7e-7),
        (:frootn_to_retransn_patch,         CLM.IFROOT,       IRET,             6e-7),
        (:livestemn_to_retransn_patch,      CLM.ILIVESTEM,    IRET,             3e-7),
        (:livecrootn_to_retransn_patch,     CLM.ILIVECROOT,   IRET,             2e-7),
        (:leafn_to_litter_patch,            CLM.ILEAF,        0,                3e-7),
        (:frootn_to_litter_patch,           CLM.IFROOT,       0,                2e-7),
        (:livestemn_to_litter_patch,        CLM.ILIVESTEM,    0,                1e-7)]

    # gap: 19 pools -> litter
    gm_rate = [1e-6,2e-6,3e-6,1.5e-6,2.5e-6,1e-6,8e-7,6e-7,4e-7,3e-7,2e-7,1e-7,7e-7,5e-7,3e-7,2e-7,1e-7,9e-8,5e-7]
    # fire: 19 out (to_fire) + 2 live->dead
    fi_out_rate = [1e-6,1.5e-6,2e-6,8e-7,1e-6,1.2e-6,7e-7,5e-7,3e-7,2e-7,1e-7,9e-8,6e-7,4e-7,2e-7,1e-7,8e-8,6e-8,4e-7]
    fi_ls_rate = 3e-7; fi_lc_rate = 2e-7
    # allocation: npool_to_Xn (fixed fluxes) + retransn_to_npool
    alloc_flux = Dict(CLM.ILEAF=>3e-5, CLM.ILEAF_ST=>1e-5, CLM.IFROOT=>2.5e-5, CLM.IFROOT_ST=>8e-6,
        CLM.ILIVESTEM=>2e-5, CLM.ILIVESTEM_ST=>6e-6, CLM.IDEADSTEM=>1e-5, CLM.IDEADSTEM_ST=>4e-6,
        CLM.ILIVECROOT=>1.5e-5, CLM.ILIVECROOT_ST=>5e-6, CLM.IDEADCROOT=>8e-6, CLM.IDEADCROOT_ST=>3e-6)
    r2n_val = [2e-5, 1.5e-5, 8e-6]     # retransn_to_npool per patch

    setv!(nf, f, v) = (length(getfield(nf, f))==0 && setfield!(nf, f, zeros(np)); getfield(nf, f) .= v)
    function set_all!(nf, ns)
        for (f, d, _, r) in ph_specs
            setv!(nf, f, 0.0); arr=getfield(nf,f); for p in 1:np; arr[p]=npget(ns,d,p)*r; end
        end
        for i in 1:19
            f = CLM._GMN_FLUX[i]; setv!(nf, f, 0.0); arr=getfield(nf,f); for p in 1:np; arr[p]=npget(ns,i,p)*gm_rate[i]; end
            ff = CLM._FIN_TO_FIRE[i]; setv!(nf, ff, 0.0); arr2=getfield(nf,ff); for p in 1:np; arr2[p]=npget(ns,i,p)*fi_out_rate[i]; end
            setv!(nf, CLM._FIN_TO_LITTER[i], 0.0)   # keep to_litter_fire zero (init NaN-fills)
        end
        setv!(nf, :m_livestemn_to_deadstemn_fire_patch, 0.0); setv!(nf, :m_livecrootn_to_deadcrootn_fire_patch, 0.0)
        for p in 1:np
            nf.m_livestemn_to_deadstemn_fire_patch[p]  = npget(ns, CLM.ILIVESTEM, p)*fi_ls_rate
            nf.m_livecrootn_to_deadcrootn_fire_patch[p] = npget(ns, CLM.ILIVECROOT, p)*fi_lc_rate
        end
        for (pool, fld) in CLM._NALLOC_TARGET; setv!(nf, fld, alloc_flux[pool]); end
        setv!(nf, :retransn_to_npool_patch, 0.0); nf.retransn_to_npool_patch .= r2n_val
    end

    # ============ (a) independent sequential reference ============
    ns_seq = fresh()
    nf_seq = CLM.CNVegNitrogenFluxData(); CLM.cnveg_nitrogen_flux_init!(nf_seq, np, 1, 1; use_matrixcn=true, nvegnpool=nveg)
    set_all!(nf_seq, ns_seq)
    for p in 1:np
        # allocation net: pools gain npool_to_Xn*dt; retransn loses retransn_to_npool*dt
        for (pool, fld) in CLM._NALLOC_TARGET; nadd!(ns_seq, pool, p, getfield(nf_seq, fld)[p]*dt); end
        nadd!(ns_seq, IRET, p, -nf_seq.retransn_to_npool_patch[p]*dt)
        # phenology transfers 1-18 + litterfall 31-33 (skip retransn->pool 19-30 = allocation)
        for (f, d, r, _) in ph_specs
            flux = getfield(nf_seq, f)[p]
            nadd!(ns_seq, d, p, -flux*dt)
            r == 0 || nadd!(ns_seq, r, p, flux*dt)
        end
        # gap losses (19 incl retransn)
        for i in 1:19; nadd!(ns_seq, i, p, -getfield(nf_seq, CLM._GMN_FLUX[i])[p]*dt); end
        # fire: 2 live->dead + 19 out losses
        fls = nf_seq.m_livestemn_to_deadstemn_fire_patch[p]; flc = nf_seq.m_livecrootn_to_deadcrootn_fire_patch[p]
        nadd!(ns_seq, CLM.ILIVESTEM, p, -fls*dt); nadd!(ns_seq, CLM.IDEADSTEM, p, fls*dt)
        nadd!(ns_seq, CLM.ILIVECROOT, p, -flc*dt); nadd!(ns_seq, CLM.IDEADCROOT, p, flc*dt)
        for i in 1:19
            nadd!(ns_seq, i, p, -(getfield(nf_seq, CLM._FIN_TO_FIRE[i])[p]+getfield(nf_seq, CLM._FIN_TO_LITTER[i])[p])*dt)
        end
    end

    # ============ (b) matrix advance via the wiring under test ============
    ns_mat = fresh()
    nf = CLM.CNVegNitrogenFluxData(); CLM.cnveg_nitrogen_flux_init!(nf, np, 1, 1; use_matrixcn=true, nvegnpool=nveg)
    set_all!(nf, ns_mat)
    CLM.cn_veg_matrix_n_topology!(nf; use_crop=false, nvegnpool=nveg)
    CLM.cn_veg_matrix_alloc_n!(nf, trues(np), 1:np)
    CLM.cn_veg_matrix_accumulate_phn!(nf, ns_mat, trues(np), 1:np; dt=dt, matrixcheck=false)
    CLM.cn_veg_matrix_accumulate_gmn!(nf, ns_mat, trues(np), 1:np; dt=dt, matrixcheck=false)
    CLM.cn_veg_matrix_accumulate_fin!(nf, ns_mat, trues(np), 1:np; dt=dt, matrixcheck=false)
    CLM.cn_veg_matrix_solve_n!(ns_mat, nf; mask_soilp=trues(np), bounds_patch=1:np,
        ivt=ivt, npcropmin=npcropmin, nvegnpool=nveg, counts=counts, dt=dt, num_actfirep=np)

    # matrix advance must reproduce the sequential update exactly (all 19 pools)
    for p in 1:np, i in 1:19
        @test isapprox(npget(ns_mat, i, p), npget(ns_seq, i, p); atol=1e-11, rtol=1e-10)
    end

    # topology sanity (exact Fortran indices)
    @test nf.ileaf_to_iretransn_ph == 1
    @test nf.iretransn_to_ideadcrootst_ph == 30
    @test nf.iretransn_to_iout_ph == 34
    @test nf.iretransn_to_iout_gm == 19
    @test nf.iretransn_to_iout_fi == 21
    @test length(nf.matrix_nphtransfer_doner_patch) == counts.nnphtrans
    @test nf.matrix_nphtransfer_receiver_patch[1] == IRET   # ileaf -> retransn
end
