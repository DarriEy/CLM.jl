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

rP1 = try section_gresp() catch e; @printf("[P1] ERRORED: %s\n", sprint(showerror,e)); NaN end
rP2 = try section_mresp() catch e; @printf("[P2] ERRORED: %s\n", sprint(showerror,e)); NaN end
println("\n", "="^70)
@printf("BGC REVERSE SUMMARY  [P1] cn_gresp=%.3e  [P2] cn_mresp=%.3e\n", rP1, rP2)
println("="^70)
