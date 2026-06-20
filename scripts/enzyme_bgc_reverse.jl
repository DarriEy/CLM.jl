# =============================================================================
# BGC (use_cn) carbon/nitrogen phase reverses through the production
# CLM.compositional_reverse! engine. Minimal SYNTHETIC CN states (mirroring the
# test_*_resp.jl / test_c_state_update*.jl patterns) — no full use_cn inst needed:
# the differentiated bundle is the relevant CN struct(s); PFT/shared params + patch +
# masks + (Const) state structs are the aux. Each section perturbs an input and seeds an
# output, FD-validating the reverse gradient.
#
#   [P1] cn_gresp!  (growth respiration:   cpool_*_gr = grperc·grpnow·allocation)
#   [P2] cn_mresp!  (maintenance resp:     livestem_mr = livestemn·br·tc(Q10), the `br` param)
#   [P3] c_state_update1! (C pool integration: leafc += cpool_to_leafc·dt — the C-cycle update)
#   [P4] n_state_update1! (N pool integration: leafn += npool_to_leafn·dt — the N-cycle update)
#   [P5] decomp_rate_constants_bgc! (soil-C turnover: t_scalar=Q10^((Tsoi-Tref)/10) — the `Q10` param)
#   [P6] soil_bgc_potential! (potential decomp + mineral-N: p_decomp_cpool_loss=Cpool·decomp_k·pathfrac)
#   [P7] soilbiogeochem_n_state_update1! (mineral-N update: smin_nh4 += (ndep/fix/nmin−immob−...)·dt)
#
#   julia +1.10 --project=/tmp/clm_jl10_<id> scripts/enzyme_bgc_reverse.jl
# =============================================================================
using CLM, Enzyme, Printf
const C = CLM; const FT = Float64
C.varpar_init!(C.varpar, 1, 14, 2, 5)   # nlevsoi/nlevgrnd/nlevsno needed by the soil-layer phases

richardson(L_pert) = (cfd(h) = (L_pert(h) - L_pert(-h)) / (2h); (4*cfd(5e-3) - cfd(1e-2)) / 3)
function verdict(tag, g_fd, g_rev)
    rl = abs(g_rev - g_fd) / max(abs(g_fd), 1e-12)
    @printf("[%s] FD = % .8e  rev = % .8e  rel = %.3e  %s\n", tag, g_fd, g_rev, rl,
            rl < 1e-5 ? "PASS ✓" : "FAIL ✗")
    return rl
end

# =============================================================================
# [P1] cn_gresp! — growth respiration. Bundle = CNVegCarbonFluxData (input cpool_to_*
#      and output cpool_*_gr both live in it). Perturb cpool_to_leafc, seed cpool_*_gr.
# =============================================================================
function section_gresp()
    println("\n", "#"^70, "\n# [P1] cn_gresp! (growth respiration) REVERSE\n", "#"^70)
    np = 1; nrepr = C.NREPR
    pftcon = C.PftConGrowthResp(woody=[0.0,1.0], grperc=[0.0,0.3], grpnow=[0.0,0.5])
    patch = C.PatchData(); C.patch_init!(patch, np); patch.itype[1] = 1
    function fresh()
        cf = C.CNVegCarbonFluxData(); C.cnveg_carbon_flux_init!(cf, np, 1, 1)
        cf.cpool_to_leafc_patch[1]=1.0; cf.cpool_to_leafc_storage_patch[1]=0.5; cf.leafc_xfer_to_leafc_patch[1]=0.2
        cf.cpool_to_frootc_patch[1]=0.8; cf.cpool_to_frootc_storage_patch[1]=0.4; cf.frootc_xfer_to_frootc_patch[1]=0.3
        cf.cpool_to_livestemc_patch[1]=0.6; cf.cpool_to_livestemc_storage_patch[1]=0.3; cf.livestemc_xfer_to_livestemc_patch[1]=0.25
        for k in 1:nrepr; cf.cpool_to_reproductivec_patch[1,k]=0.0; cf.cpool_to_reproductivec_storage_patch[1,k]=0.0; cf.reproductivec_xfer_to_reproductivec_patch[1,k]=0.0; end
        return cf
    end
    aux = (; mask=BitVector([true]), bounds=1:np, pftcon=pftcon, patch=patch, npcropmin=17, nrepr=nrepr)
    phase!(cf, a) = (C.cn_gresp!(a.mask, a.bounds, a.pftcon, a.patch, cf; npcropmin=a.npcropmin, nrepr=a.nrepr); nothing)
    grf = (:cpool_leaf_gr_patch,:cpool_froot_gr_patch,:cpool_livestem_gr_patch)
    L(cf) = sum(sum(abs2, getfield(cf,f)) for f in grf)
    g_fd = richardson(δ -> (cf=fresh(); cf.cpool_to_leafc_patch[1]+=δ; phase!(cf,aux); L(cf)))
    db = C.compositional_reverse!(Any[(phase!,(aux,))], fresh(),
        (db,b)->(for f in grf; getfield(db,f) .= 2 .* getfield(b,f); end))
    return verdict("P1 cn_gresp", g_fd, db.cpool_to_leafc_patch[1])
end

# =============================================================================
# [P2] cn_mresp! — maintenance respiration. Bundle = (cnveg_cf out, cnveg_ns in); the
#      MaintResp/shared/PFT params + patch + canopy/soil/temperature/photosyns structs are
#      Const aux. Perturb livestemn → livestem_mr = livestemn·br·tc (t_soisno=293.15 → tc=1).
# =============================================================================
function section_mresp()
    println("\n", "#"^70, "\n# [P2] cn_mresp! (maintenance respiration) REVERSE\n", "#"^70)
    np=1; nc=1; nlevgrnd=C.varpar.nlevgrnd; nlevsno=C.varpar.nlevsno; nrepr=C.NREPR; npcropmin=17
    br=2.525e-6; br_root=0.83e-6; Q10=1.5
    params = C.MaintRespParams(br=br, br_root=br_root)
    cn_params = C.CNSharedParamsData(Q10=Q10)
    pftcon = C.PftConMaintResp(woody=[0.0, 1.0])
    patch = C.PatchData(); C.patch_init!(patch, np); patch.itype[1]=1
    cs = C.CanopyStateData(); C.canopystate_init!(cs, np)
    ss = C.SoilStateData(); C.soilstate_init!(ss, np, nc)
    temp = C.TemperatureData(); C.temperature_init!(temp, np, nc, 1, 1)
    temp.t_ref2m_patch[1] = 293.15; temp.t_a10_patch[1] = 293.15  # above-ground MR temps (tc=1)
    for j in 1:nlevgrnd; temp.t_soisno_col[1, j + nlevsno] = 293.15; end
    ps = C.PhotosynthesisData(); C.photosynthesis_data_init!(ps, np)
    function fresh()
        cf = C.CNVegCarbonFluxData(); C.cnveg_carbon_flux_init!(cf, np, nc, 1)
        ns = C.CNVegNitrogenStateData(); C.cnveg_nitrogen_state_init!(ns, np, nc, 1)
        ns.livestemn_patch[1]=0.05; ns.livecrootn_patch[1]=0.03; ns.frootn_patch[1]=0.01
        return (; cnveg_cf=cf, cnveg_ns=ns)
    end
    aux = (; mask_c=BitVector([true]), mask_p=BitVector([true]), bounds_c=1:nc, bounds_p=1:np,
             params=params, cn_params=cn_params, pftcon=pftcon, patch=patch, cs=cs, ss=ss, temp=temp, ps=ps,
             nlevgrnd=nlevgrnd, nlevsno=nlevsno, npcropmin=npcropmin, nrepr=nrepr)
    phase!(b, a) = (C.cn_mresp!(a.mask_c, a.mask_p, a.bounds_c, a.bounds_p, a.params, a.cn_params,
        a.pftcon, a.patch, a.cs, a.ss, a.temp, a.ps, b.cnveg_cf, b.cnveg_ns;
        nlevgrnd=a.nlevgrnd, nlevsno=a.nlevsno, npcropmin=a.npcropmin, nrepr=a.nrepr); nothing)
    let b=fresh(); phase!(b,aux)
        @printf("primal: livestem_mr=%.4e (expect livestemn·br=%.4e) finite=%s\n",
            b.cnveg_cf.livestem_mr_patch[1], 0.05*br, string(isfinite(b.cnveg_cf.livestem_mr_patch[1])))
    end
    L(b) = sum(abs2, b.cnveg_cf.livestem_mr_patch) + sum(abs2, b.cnveg_cf.livecroot_mr_patch)
    g_fd = richardson(δ -> (b=fresh(); b.cnveg_ns.livestemn_patch[1]+=δ; phase!(b,aux); L(b)))
    db = C.compositional_reverse!(Any[(phase!,(aux,))], fresh(),
        (db,b)->(db.cnveg_cf.livestem_mr_patch .= 2 .* b.cnveg_cf.livestem_mr_patch;
                 db.cnveg_cf.livecroot_mr_patch .= 2 .* b.cnveg_cf.livecroot_mr_patch))
    return verdict("P2 cn_mresp", g_fd, db.cnveg_ns.livestemn_patch[1])
end

# =============================================================================
# [P3] c_state_update1! — carbon pool integration (the C-cycle state update). Bundle =
#      (cs_veg pools out, cf_veg fluxes in, cf_soil litter sink). All index/cascade args
#      (patch_column, cascade_donor/receiver_pool, woody, …) are Const aux. Perturb
#      cpool_to_leafc → leafc_patch (leafc += cpool_to_leafc·dt, c_state_update1.jl:489).
# =============================================================================
function section_cstate1()
    println("\n", "#"^70, "\n# [P3] c_state_update1! (carbon pool integration) REVERSE\n", "#"^70)
    np=1; nc=1; ng=1; nlevdecomp=1; ndecomp_pools=7; ndecomp_cascade_transitions=5; nrepr=C.NREPR
    i_litr_min=1; i_litr_max=3; i_cwd=4; dt=1800.0
    patch_column=[1]; ivt=[1]; woody=zeros(Float64,80); harvdate=fill(999,np); col_is_fates=fill(false,nc)
    cascade_donor_pool=[1,2,3,1,2]; cascade_receiver_pool=[2,3,0,4,4]
    function fresh()
        cs = C.CNVegCarbonStateData(); C.cnveg_carbon_state_init!(cs, np, nc, ng; nrepr=nrepr)
        cs.cpool_patch .= 100.0; cs.leafc_patch .= 50.0
        cf = C.CNVegCarbonFluxData()
        C.cnveg_carbon_flux_init!(cf, np, nc, ng; nrepr=nrepr, nlevdecomp_full=nlevdecomp, ndecomp_pools=ndecomp_pools)
        for f in fieldnames(typeof(cf))          # init NaN-fills (surfaces unset fluxes) → zero for a clean probe
            v = getfield(cf, f); v isa AbstractArray{<:AbstractFloat} && fill!(v, 0.0)
        end
        cf.cpool_to_leafc_patch[1]=5.0e-7; cf.cpool_to_leafc_storage_patch[1]=2.0e-7
        cf.leafc_xfer_to_leafc_patch[1]=2.0e-7
        cf.psnsun_to_cpool_patch[1]=1.0e-6; cf.psnshade_to_cpool_patch[1]=0.5e-6
        cfs = C.SoilBiogeochemCarbonFluxData()
        C.soil_bgc_carbon_flux_init!(cfs, nc, nlevdecomp, ndecomp_pools, ndecomp_cascade_transitions)
        cfs.decomp_cpools_sourcesink_col .= 0.0
        return (; cs_veg=cs, cf_veg=cf, cf_soil=cfs)
    end
    aux = (; mask_c=BitVector([true]), mask_p=BitVector([true]), bounds_c=1:nc, bounds_p=1:np,
             patch_column=patch_column, ivt=ivt, woody=woody,
             cascade_donor_pool=cascade_donor_pool, cascade_receiver_pool=cascade_receiver_pool,
             harvdate=harvdate, col_is_fates=col_is_fates, nlevdecomp=nlevdecomp,
             ndecomp_cascade_transitions=ndecomp_cascade_transitions,
             i_litr_min=i_litr_min, i_litr_max=i_litr_max, i_cwd=i_cwd, nrepr=nrepr, dt=dt)
    phase!(b, a) = (C.c_state_update1!(b.cs_veg, b.cf_veg, b.cf_soil;
        mask_soilc=a.mask_c, mask_soilp=a.mask_p, bounds_col=a.bounds_c, bounds_patch=a.bounds_p,
        patch_column=a.patch_column, ivt=a.ivt, woody=a.woody,
        cascade_donor_pool=a.cascade_donor_pool, cascade_receiver_pool=a.cascade_receiver_pool,
        harvdate=a.harvdate, col_is_fates=a.col_is_fates, nlevdecomp=a.nlevdecomp,
        ndecomp_cascade_transitions=a.ndecomp_cascade_transitions,
        i_litr_min=a.i_litr_min, i_litr_max=a.i_litr_max, i_cwd=a.i_cwd, nrepr=a.nrepr, dt=a.dt); nothing)
    let b=fresh(); l0=b.cs_veg.leafc_patch[1]; phase!(b,aux)
        @printf("primal: leafc %.4f→%.8f (Δ=%.4e expect cpool_to_leafc·dt=%.4e)\n",
            l0, b.cs_veg.leafc_patch[1], b.cs_veg.leafc_patch[1]-l0, 5.0e-7*dt)
    end
    L(b) = sum(abs2, b.cs_veg.leafc_patch)
    g_fd = richardson(δ -> (b=fresh(); b.cf_veg.cpool_to_leafc_patch[1]+=δ; phase!(b,aux); L(b)))
    db = C.compositional_reverse!(Any[(phase!,(aux,))], fresh(),
        (db,b)->(db.cs_veg.leafc_patch .= 2 .* b.cs_veg.leafc_patch))
    return verdict("P3 c_state_update1", g_fd, db.cf_veg.cpool_to_leafc_patch[1])
end

# =============================================================================
# [P4] n_state_update1! — nitrogen pool integration (the N-cycle state update; mirror of
#      [P3] with N structs, fewer consts — no cascade pools/harvdate). Bundle = (ns_veg N
#      pools out, nf_veg N fluxes in, nf_soil litter-N sink). Perturb npool_to_leafn →
#      leafn_patch (leafn += npool_to_leafn·dt, n_state_update1.jl:362).
# =============================================================================
function section_nstate1()
    println("\n", "#"^70, "\n# [P4] n_state_update1! (nitrogen pool integration) REVERSE\n", "#"^70)
    np=1; nc=1; ng=1; nlevdecomp=1; ndecomp_pools=7; ndecomp_cascade_transitions=5; nrepr=C.NREPR
    i_litr_min=1; i_litr_max=3; i_cwd=4; dt=1800.0
    ivt=[1]; woody=zeros(Float64,80); col_is_fates=fill(false,nc)
    function fresh()
        ns = C.CNVegNitrogenStateData(); C.cnveg_nitrogen_state_init!(ns, np, nc, ng; nrepr=nrepr)
        ns.leafn_patch .= 5.0; ns.leafn_xfer_patch .= 0.5; ns.npool_patch .= 10.0
        nf = C.CNVegNitrogenFluxData()
        C.cnveg_nitrogen_flux_init!(nf, np, nc, ng; nrepr=nrepr, nlevdecomp_full=nlevdecomp,
                                    ndecomp_pools=ndecomp_pools, i_litr_max=i_litr_max)
        for f in fieldnames(typeof(nf))          # init NaN-fills → zero for a clean probe
            v = getfield(nf, f); v isa AbstractArray{<:AbstractFloat} && fill!(v, 0.0)
        end
        nf.npool_to_leafn_patch[1]=5.0e-7; nf.npool_to_leafn_storage_patch[1]=2.0e-7
        nf.leafn_xfer_to_leafn_patch[1]=2.0e-7
        nfs = C.SoilBiogeochemNitrogenFluxData()
        C.soil_bgc_nitrogen_flux_init!(nfs, nc, nlevdecomp, ndecomp_pools, ndecomp_cascade_transitions)
        nfs.decomp_npools_sourcesink_col .= 0.0
        return (; ns_veg=ns, nf_veg=nf, nf_soil=nfs)
    end
    aux = (; mask_c=BitVector([true]), mask_p=BitVector([true]), bounds_c=1:nc, bounds_p=1:np,
             ivt=ivt, woody=woody, col_is_fates=col_is_fates, nlevdecomp=nlevdecomp,
             i_litr_min=i_litr_min, i_litr_max=i_litr_max, i_cwd=i_cwd, nrepr=nrepr, dt=dt)
    phase!(b, a) = (C.n_state_update1!(b.ns_veg, b.nf_veg, b.nf_soil;
        mask_soilc=a.mask_c, mask_soilp=a.mask_p, bounds_col=a.bounds_c, bounds_patch=a.bounds_p,
        ivt=a.ivt, woody=a.woody, col_is_fates=a.col_is_fates, nlevdecomp=a.nlevdecomp,
        i_litr_min=a.i_litr_min, i_litr_max=a.i_litr_max, i_cwd=a.i_cwd, nrepr=a.nrepr, dt=a.dt); nothing)
    let b=fresh(); l0=b.ns_veg.leafn_patch[1]; phase!(b,aux)
        @printf("primal: leafn %.4f→%.8f (Δ=%.4e expect npool_to_leafn·dt=%.4e)\n",
            l0, b.ns_veg.leafn_patch[1], b.ns_veg.leafn_patch[1]-l0, 5.0e-7*dt)
    end
    L(b) = sum(abs2, b.ns_veg.leafn_patch)
    g_fd = richardson(δ -> (b=fresh(); b.nf_veg.npool_to_leafn_patch[1]+=δ; phase!(b,aux); L(b)))
    db = C.compositional_reverse!(Any[(phase!,(aux,))], fresh(),
        (db,b)->(db.ns_veg.leafn_patch .= 2 .* b.ns_veg.leafn_patch))
    return verdict("P4 n_state_update1", g_fd, db.nf_veg.npool_to_leafn_patch[1])
end

# =============================================================================
# [P5] decomp_rate_constants_bgc! — soil-carbon decomposition rate constants (the Q10/tau
#      turnover calibration). Bundle = (cf out [t_scalar/w_scalar/decomp_k], t_soisno in —
#      the kwarg temperature array); params/cn_params(Q10)/bgc_state/cascade_con/soilpsi/
#      zsoi/col_dz are Const aux. Perturb t_soisno, seed t_scalar = Q10^((Tsoi-Tref)/10).
# =============================================================================
function section_decomprate()
    println("\n", "#"^70, "\n# [P5] decomp_rate_constants_bgc! (soil-C turnover rates) REVERSE\n", "#"^70)
    nc=1; nlevdecomp=1; ndecomp_pools=7; ndecomp_cascade_transitions=10; nlev=max(nlevdecomp,5)
    params = C.DecompBGCParams(cn_s1_bgc=12.0, cn_s2_bgc=12.0, cn_s3_bgc=10.0,
        rf_l1s1_bgc=0.39, rf_l2s1_bgc=0.55, rf_l3s2_bgc=0.29, rf_s2s1_bgc=0.55, rf_s2s3_bgc=0.55,
        rf_s3s1_bgc=0.55, rf_cwdl3_bgc=0.0, tau_l1_bgc=1.0/18.5, tau_l2_l3_bgc=1.0/4.9,
        tau_s1_bgc=1.0/7.3, tau_s2_bgc=1.0/0.2, tau_s3_bgc=1.0/0.0045, cwd_fcel_bgc=0.45,
        bgc_initial_Cstocks=fill(200.0,7), bgc_initial_Cstocks_depth=0.3)
    cn_params = C.CNSharedParamsData(Q10=1.5, minpsi=-10.0, maxpsi=-0.1, rf_cwdl2=0.0, tau_cwd=10.0,
        cwd_flig=0.24, froz_q10=1.5, decomp_depth_efolding=0.5, mino2lim=0.0)
    bgc_state = C.DecompBGCState(); cascade_con = C.DecompCascadeConData()
    soilpsi = fill(-1.0, nc, nlev)
    zsoi_vals = [0.01,0.04,0.09,0.16,0.26,0.40,0.58,0.80,1.06,1.36]
    col_dz = fill(0.1, nc, nlev)
    function fresh()
        cf = C.SoilBiogeochemCarbonFluxData()
        C.soil_bgc_carbon_flux_init!(cf, nc, nlev, ndecomp_pools, ndecomp_cascade_transitions)
        t_soisno = fill(C.TFRZ + 15.0, nc, nlev)
        return (; cf=cf, t_soisno=t_soisno)
    end
    aux = (; bgc_state=bgc_state, params=params, cn_params=cn_params, cascade_con=cascade_con,
             mask=BitVector([true]), bounds=1:nc, nlevdecomp=nlevdecomp,
             soilpsi=soilpsi, zsoi_vals=zsoi_vals, col_dz=col_dz)
    phase!(b, a) = (C.decomp_rate_constants_bgc!(b.cf, a.bgc_state, a.params, a.cn_params, a.cascade_con;
        mask_bgc_soilc=a.mask, bounds=a.bounds, nlevdecomp=a.nlevdecomp, t_soisno=b.t_soisno,
        soilpsi=a.soilpsi, days_per_year=365.0, dt=1800.0, zsoi_vals=a.zsoi_vals, col_dz=a.col_dz); nothing)
    let b=fresh(); phase!(b,aux)
        @printf("primal: t_scalar=%.6f w_scalar=%.6f finite=%s\n",
            b.cf.t_scalar_col[1,1], b.cf.w_scalar_col[1,1], string(isfinite(b.cf.t_scalar_col[1,1])))
    end
    L(b) = sum(abs2, @view b.cf.t_scalar_col[:, 1])
    g_fd = richardson(δ -> (b=fresh(); b.t_soisno[1,1]+=δ; phase!(b,aux); L(b)))
    db = C.compositional_reverse!(Any[(phase!,(aux,))], fresh(),
        (db,b)->(db.cf.t_scalar_col[:,1] .= 2 .* @view b.cf.t_scalar_col[:,1]))
    return verdict("P5 decomp_rate_bgc", g_fd, db.t_soisno[1,1])
end

# =============================================================================
# [P6] soil_bgc_potential! — potential decomposition + mineral-N flux (consumes decomp_k
#      from [P5]): p_decomp_cpool_loss[tr] = Cpool[donor]·decomp_k[donor]·pathfrac. Bundle =
#      (cf in [decomp_k], ploss out [the p_decomp_cpool_loss scratch]); cs/ns/st/cascade_con
#      + the other potential-flux scratch arrays are Const aux. Perturb decomp_k, seed ploss.
# =============================================================================
function section_potential()
    println("\n", "#"^70, "\n# [P6] soil_bgc_potential! (potential decomp + mineral-N) REVERSE\n", "#"^70)
    nc=1; nlevdecomp=1; ndecomp_pools=7; nct=5
    cascade_con = C.DecompCascadeConData()
    cascade_con.cascade_donor_pool    = [1, 2, 3, 3, 4]
    cascade_con.cascade_receiver_pool = [3, 3, 4, C.I_ATM, 3]
    cascade_con.floating_cn_ratio_decomp_pools = BitVector([true, true, false, false, false, false, false])
    cascade_con.initial_cn_ratio = [20.0, 25.0, 12.0, 10.0, 10.0, 10.0, 10.0]
    cs = C.SoilBiogeochemCarbonStateData(); cs.decomp_cpools_vr_col = zeros(nc, nlevdecomp, ndecomp_pools)
    cs.decomp_cpools_vr_col[1,1,1]=110.0; cs.decomp_cpools_vr_col[1,1,2]=85.0
    cs.decomp_cpools_vr_col[1,1,3]=53.0; cs.decomp_cpools_vr_col[1,1,4]=32.0
    ns = C.SoilBiogeochemNitrogenStateData(); ns.decomp_npools_vr_col = zeros(nc, nlevdecomp, ndecomp_pools)
    ns.decomp_npools_vr_col[1,1,1]=110/20; ns.decomp_npools_vr_col[1,1,2]=85/16
    ns.decomp_npools_vr_col[1,1,3]=53/12; ns.decomp_npools_vr_col[1,1,4]=32/10
    st = C.SoilBiogeochemStateData(); st.nue_decomp_cascade_col = ones(nct) .* 0.5
    nf = C.SoilBiogeochemNitrogenFluxData()
    nf.potential_immob_vr_col = zeros(nc, nlevdecomp); nf.gross_nmin_vr_col = zeros(nc, nlevdecomp)
    function fresh()
        cf = C.SoilBiogeochemCarbonFluxData()
        cf.rf_decomp_cascade_col = zeros(nc, nlevdecomp, nct); cf.pathfrac_decomp_cascade_col = ones(nc, nlevdecomp, nct)
        cf.decomp_k_col = zeros(nc, nlevdecomp, ndecomp_pools); cf.cn_col = zeros(nc, ndecomp_pools); cf.phr_vr_col = zeros(nc, nlevdecomp)
        cf.rf_decomp_cascade_col[1,1,:] .= [0.39,0.55,0.28,0.55,0.55]
        cf.decomp_k_col[1,1,1]=0.01; cf.decomp_k_col[1,1,2]=0.005; cf.decomp_k_col[1,1,3]=0.002; cf.decomp_k_col[1,1,4]=0.001
        ploss = zeros(nc, nlevdecomp, nct)
        return (; cf=cf, ploss=ploss)
    end
    aux = (; cs=cs, nf=nf, ns=ns, st=st, cascade_con=cascade_con, mask=BitVector([true]), bounds=1:nc,
             nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools, nct=nct,
             cn_decomp_pools=zeros(nc,nlevdecomp,ndecomp_pools), p_decomp_cn_gain=zeros(nc,nlevdecomp,ndecomp_pools),
             pmnf=zeros(nc,nlevdecomp,nct), p_npool_to_din=zeros(nc,nlevdecomp,nct))
    phase!(b, a) = (C.soil_bgc_potential!(b.cf, a.cs, a.nf, a.ns, a.st, a.cascade_con;
        mask_bgc_soilc=a.mask, bounds=a.bounds, nlevdecomp=a.nlevdecomp, ndecomp_pools=a.ndecomp_pools,
        ndecomp_cascade_transitions=a.nct, cn_decomp_pools=a.cn_decomp_pools,
        p_decomp_cpool_loss=b.ploss, p_decomp_cn_gain=a.p_decomp_cn_gain,
        pmnf_decomp_cascade=a.pmnf, p_decomp_npool_to_din=a.p_npool_to_din); nothing)
    let b=fresh(); phase!(b,aux)
        @printf("primal: p_decomp_cpool_loss[1]=%.6e (expect Cpool1·k1=%.6e) finite=%s\n",
            b.ploss[1,1,1], 110.0*0.01, string(isfinite(b.ploss[1,1,1])))
    end
    L(b) = sum(abs2, b.ploss)
    g_fd = richardson(δ -> (b=fresh(); b.cf.decomp_k_col[1,1,1]+=δ; phase!(b,aux); L(b)))
    db = C.compositional_reverse!(Any[(phase!,(aux,))], fresh(),
        (db,b)->(db.ploss .= 2 .* b.ploss))
    return verdict("P6 soil_bgc_potential", g_fd, db.cf.decomp_k_col[1,1,1])
end

# =============================================================================
# [P7] soilbiogeochem_n_state_update1! — mineral-N pool update (applies dep/fix, gross
#      mineralization, immobilization, plant uptake, nitrif/denitrif, supplement to
#      smin_nh4_vr/smin_no3_vr). Bundle = (ns out [smin_nh4/no3], nf in [the mineral-N
#      fluxes]); st (ndep/nfixation profiles) is Const aux. Perturb gross_nmin → smin_nh4.
# =============================================================================
function section_sminn()
    println("\n", "#"^70, "\n# [P7] soilbiogeochem_n_state_update1! (mineral-N update) REVERSE\n", "#"^70)
    nc=1; nlevdecomp=1; dt=1800.0
    st = C.SoilBiogeochemStateData()
    st.ndep_prof_col = zeros(nc, nlevdecomp); st.nfixation_prof_col = zeros(nc, nlevdecomp)
    function fresh()
        ns = C.SoilBiogeochemNitrogenStateData()
        ns.smin_nh4_vr_col = fill(2.0, nc, nlevdecomp); ns.smin_no3_vr_col = fill(1.0, nc, nlevdecomp)
        ns.sminn_vr_col = zeros(nc, nlevdecomp)
        nf = C.SoilBiogeochemNitrogenFluxData()
        for f in (:gross_nmin_vr_col, :actual_immob_nh4_vr_col, :actual_immob_no3_vr_col,
                  :smin_nh4_to_plant_vr_col, :smin_no3_to_plant_vr_col, :f_nit_vr_col,
                  :f_denit_vr_col, :supplement_to_sminn_vr_col)
            setfield!(nf, f, zeros(nc, nlevdecomp))
        end
        nf.ndep_to_sminn_col = zeros(nc); nf.nfix_to_sminn_col = zeros(nc)
        nf.gross_nmin_vr_col[1,1] = 1.0e-6
        return (; ns=ns, nf=nf)
    end
    aux = (; st=st, mask=BitVector([true]), bounds=1:nc, nlevdecomp=nlevdecomp, dt=dt)
    phase!(b, a) = (C.soilbiogeochem_n_state_update1!(b.ns, b.nf, a.st;
        mask_bgc_soilc=a.mask, bounds_col=a.bounds, nlevdecomp=a.nlevdecomp, dt=a.dt, use_fun=false); nothing)
    let b=fresh(); s0=b.ns.smin_nh4_vr_col[1,1]; phase!(b,aux)
        @printf("primal: smin_nh4 %.6f→%.8f (Δ=%.4e expect gross_nmin·dt=%.4e)\n",
            s0, b.ns.smin_nh4_vr_col[1,1], b.ns.smin_nh4_vr_col[1,1]-s0, 1.0e-6*dt)
    end
    L(b) = sum(abs2, b.ns.smin_nh4_vr_col)
    g_fd = richardson(δ -> (b=fresh(); b.nf.gross_nmin_vr_col[1,1]+=δ; phase!(b,aux); L(b)))
    db = C.compositional_reverse!(Any[(phase!,(aux,))], fresh(),
        (db,b)->(db.ns.smin_nh4_vr_col .= 2 .* b.ns.smin_nh4_vr_col))
    return verdict("P7 sminn_update", g_fd, db.nf.gross_nmin_vr_col[1,1])
end

rP1 = try section_gresp() catch e; @printf("[P1] ERRORED: %s\n", sprint(showerror,e)); NaN end
rP2 = try section_mresp() catch e; @printf("[P2] ERRORED: %s\n", sprint(showerror,e)); NaN end
rP3 = try section_cstate1() catch e; @printf("[P3] ERRORED: %s\n", sprint(showerror,e)); NaN end
rP4 = try section_nstate1() catch e; @printf("[P4] ERRORED: %s\n", sprint(showerror,e)); NaN end
rP5 = try section_decomprate() catch e; @printf("[P5] ERRORED: %s\n", sprint(showerror,e)); NaN end
rP6 = try section_potential() catch e; @printf("[P6] ERRORED: %s\n", sprint(showerror,e)); NaN end
rP7 = try section_sminn() catch e; @printf("[P7] ERRORED: %s\n", sprint(showerror,e)); NaN end
println("\n", "="^70)
@printf("BGC REVERSE SUMMARY  [P1] cn_gresp=%.3e  [P2] cn_mresp=%.3e  [P3] c_state_update1=%.3e  [P4] n_state_update1=%.3e  [P5] decomp_rate=%.3e  [P6] potential=%.3e  [P7] sminn=%.3e\n", rP1, rP2, rP3, rP4, rP5, rP6, rP7)
println("="^70)
