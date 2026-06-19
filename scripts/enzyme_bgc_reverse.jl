# =============================================================================
# FIRST BGC (use_cn) phase reverse — proves the carbon/nitrogen biogeochemistry
# domain reverse-differentiates through the production CLM.compositional_reverse!
# engine, on the cleanest central CN phase: cn_gresp! (growth respiration,
# cpool_*_gr = grperc · grpnow · allocation fluxes — smooth/linear).
#
# Minimal SYNTHETIC CN state (mirrors test/test_growth_resp.jl) — no full use_cn
# inst needed: the differentiated bundle is the CNVegCarbonFluxData; the PFT growth
# params + patch + masks are Const aux. Perturb an allocation flux (cpool_to_leafc),
# seed the growth-respiration outputs (cpool_*_gr), FD-validate the reverse gradient.
#   julia +1.10 --project=/tmp/clm_jl10_<id> scripts/enzyme_bgc_reverse.jl
# =============================================================================
using CLM, Enzyme, Printf
const C = CLM; const FT = Float64

# ---- minimal synthetic CN carbon-flux state (test_growth_resp.jl pattern) ----
np = 1; nrepr = C.NREPR
grperc_val = 0.3; grpnow_val = 0.5
pftcon = C.PftConGrowthResp(woody = [0.0, 1.0], grperc = [0.0, grperc_val], grpnow = [0.0, grpnow_val])
patch = C.PatchData(); C.patch_init!(patch, np); patch.itype[1] = 1

function fresh_cf()
    cf = C.CNVegCarbonFluxData(); C.cnveg_carbon_flux_init!(cf, np, 1, 1)
    cf.cpool_to_leafc_patch[1] = 1.0;   cf.cpool_to_leafc_storage_patch[1] = 0.5
    cf.leafc_xfer_to_leafc_patch[1] = 0.2
    cf.cpool_to_frootc_patch[1] = 0.8;  cf.cpool_to_frootc_storage_patch[1] = 0.4
    cf.frootc_xfer_to_frootc_patch[1] = 0.3
    cf.cpool_to_livestemc_patch[1] = 0.6; cf.cpool_to_livestemc_storage_patch[1] = 0.3
    cf.livestemc_xfer_to_livestemc_patch[1] = 0.25
    for k in 1:nrepr
        cf.cpool_to_reproductivec_patch[1, k] = 0.0
        cf.cpool_to_reproductivec_storage_patch[1, k] = 0.0
        cf.reproductivec_xfer_to_reproductivec_patch[1, k] = 0.0
    end
    return cf
end

aux = (; mask = BitVector([true]), bounds = 1:np, pftcon = pftcon, patch = patch,
         npcropmin = 17, nrepr = nrepr)
function gresp_phase!(cf, aux)
    C.cn_gresp!(aux.mask, aux.bounds, aux.pftcon, aux.patch, cf;
                npcropmin = aux.npcropmin, nrepr = aux.nrepr)
    return nothing
end

# growth-respiration outputs (the live state the phase writes)
gr_fields = (:cpool_leaf_gr_patch, :cpool_leaf_storage_gr_patch, :transfer_leaf_gr_patch,
             :cpool_froot_gr_patch, :cpool_froot_storage_gr_patch,
             :cpool_livestem_gr_patch, :cpool_livestem_storage_gr_patch)
L(cf) = sum(sum(abs2, getfield(cf, f)) for f in gr_fields)

# primal sanity
let cf = fresh_cf(); gresp_phase!(cf, aux)
    @printf("primal: cpool_leaf_gr=%.6f (expect grperc·cpool_to_leafc=%.6f) finite=%s\n",
        cf.cpool_leaf_gr_patch[1], grperc_val * 1.0, string(isfinite(cf.cpool_leaf_gr_patch[1])))
end

# FD: perturb cpool_to_leafc → growth-resp outputs
function L_pert(δ)
    cf = fresh_cf(); cf.cpool_to_leafc_patch[1] += δ
    gresp_phase!(cf, aux); return L(cf)
end
hs = (1e-2, 5e-3); cfd(h) = (L_pert(h) - L_pert(-h)) / (2h)
g_fd = (4*cfd(hs[2]) - cfd(hs[1])) / 3

# reverse through the production engine
seed!(db, b) = (for f in gr_fields; getfield(db, f) .= 2 .* getfield(b, f); end)
db = C.compositional_reverse!(Any[(gresp_phase!, (aux,))], fresh_cf(), seed!)
g_rev = db.cpool_to_leafc_patch[1]
rl = abs(g_rev - g_fd) / max(abs(g_fd), 1e-12)

println("\n", "="^70)
@printf("BGC cn_gresp! REVERSE: perturb cpool_to_leafc, L=sum(cpool_*_gr^2)\n")
@printf("  FD = % .8e   rev = % .8e   rel = %.3e   %s\n", g_fd, g_rev, rl, rl < 1e-5 ? "PASS ✓" : "FAIL ✗")
println("="^70)
