# ==========================================================================
# gpu_validate_growth_resp_e2e.jl — end-to-end Metal parity for cn_gresp!
# (growth respiration; one per-patch kernel grouped into _GrespOut/_GrespIn
# device-view bundles incl. reproductive [p,k] matrix fields).
#
#   julia --project=scripts scripts/gpu_validate_growth_resp_e2e.jl
# ==========================================================================
using CLM
using Printf
import Metal
include(joinpath(@__DIR__, "gpu_backends.jl"))

struct _F32 end
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:AbstractFloat}) = Float32.(x)
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:Integer}) = x
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{Bool}) = x

function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + max(abs(Float64(A[i])), abs(Float64(B[i])))))
    end
    return m
end
allfinite(a) = all(isfinite, Array(a))

function build(; np=2, woody=true)
    nrepr = CLM.NREPR
    pftcon = CLM.PftConGrowthResp(woody=[0.0, woody ? 1.0 : 0.0, 1.0],
        grperc=[0.0, 0.3, 0.3], grpnow=[0.0, 0.5, 0.5])
    patch = CLM.PatchData(); CLM.patch_init!(patch, np)
    patch.itype .= 1   # ivt+1=2 (woody)
    cnveg_cf = CLM.CNVegCarbonFluxData(); CLM.cnveg_carbon_flux_init!(cnveg_cf, np, 1, 1)
    for (f,v) in ((:cpool_to_leafc_patch,1.0),(:cpool_to_leafc_storage_patch,0.5),(:leafc_xfer_to_leafc_patch,0.2),
                  (:cpool_to_frootc_patch,0.8),(:cpool_to_frootc_storage_patch,0.4),(:frootc_xfer_to_frootc_patch,0.15),
                  (:cpool_to_livestemc_patch,0.6),(:cpool_to_livestemc_storage_patch,0.3),(:livestemc_xfer_to_livestemc_patch,0.1),
                  (:cpool_to_deadstemc_patch,0.5),(:cpool_to_deadstemc_storage_patch,0.25),(:deadstemc_xfer_to_deadstemc_patch,0.08),
                  (:cpool_to_livecrootc_patch,0.4),(:cpool_to_livecrootc_storage_patch,0.2),(:livecrootc_xfer_to_livecrootc_patch,0.06),
                  (:cpool_to_deadcrootc_patch,0.3),(:cpool_to_deadcrootc_storage_patch,0.15),(:deadcrootc_xfer_to_deadcrootc_patch,0.05))
        getfield(cnveg_cf,f) .= v
    end
    return (; pftcon, patch, cnveg_cf, mask=BitVector(fill(true,np)), bounds=1:np, nrepr)
end

run_g!(d) = CLM.cn_gresp!(d.mask, d.bounds, d.pftcon, d.patch, d.cnveg_cf; npcropmin=17, nrepr=d.nrepr)

const OUT = (:cpool_leaf_gr_patch, :transfer_leaf_gr_patch, :cpool_froot_gr_patch,
             :cpool_livestem_gr_patch, :cpool_deadstem_gr_patch, :cpool_livecroot_gr_patch,
             :cpool_deadcroot_gr_patch, :cpool_livestem_storage_gr_patch)

function main(backend)
    println("=" ^ 64); println("END-TO-END Metal parity for cn_gresp!"); println("=" ^ 64)
    if backend === nothing; println("  No GPU backend."); return 0; end
    name, _, FT = backend
    @printf("  Backend: %s   (precision: %s)\n", name, FT)
    H = build(); B = build()
    run_g!(H)
    mf(x) = CLM.Adapt.adapt(Metal.MtlArray, CLM.Adapt.adapt(_F32(), x))
    D = (; pftcon=mf(B.pftcon), patch=mf(B.patch), cnveg_cf=mf(B.cnveg_cf),
         mask=Metal.MtlArray(collect(B.mask)), bounds=B.bounds, nrepr=B.nrepr)
    if !(D.cnveg_cf.cpool_leaf_gr_patch isa Metal.MtlArray); println("  BLOCKED."); return 2; end
    run_g!(D)
    nfail = 0
    for f in OUT
        a = getfield(H.cnveg_cf, f); b = getfield(D.cnveg_cf, f)
        allfinite(a) || (@printf("  [WARN] %-30s non-finite\n", f); continue)
        dd = reldiff(a, b); ok = dd < 1f-3
        @printf("  [%s] %-30s rel = %.3e\n", ok ? "PASS" : "FAIL", f, dd); ok || (nfail += 1)
    end
    println(); println(nfail == 0 ? "  cn_gresp! MATCHES CPU ON $name ($FT)" : "  DIVERGENCE ($nfail).")
    return nfail == 0 ? 0 : 1
end
const BACKEND = detect_backend()
exit(main(BACKEND))
