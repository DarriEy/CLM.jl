# ==========================================================================
# gpu_validate_cnveg_summaries.jl — host-vs-Metal parity for the 3 cnveg per-patch
# summary reductions (cnveg_carbon_state_summary! / cnveg_nitrogen_state_summary! /
# cnveg_carbon_flux_summary!), kernelized as one-thread-per-patch reductions.
# Exercises the crop reproductive branch (use_crop + nrepr=2 → the `for k` path) so
# a dropped reproductive term or the KA-CPU accumulation trap would show as a mismatch.
#
#   julia --project=scripts scripts/gpu_validate_cnveg_summaries.jl
# ==========================================================================
using CLM
using Printf
import Metal
include(joinpath(@__DIR__, "gpu_backends.jl"))
include(joinpath(@__DIR__, "gpu_adapt.jl"))
mf(x) = mf(Metal.MtlArray, x)

function reldiff(H, D)
    A = Array(H); B = Array(D); m = 0.0; n = 0
    for i in eachindex(A, B)
        (isfinite(A[i]) && isfinite(B[i])) || continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + abs(Float64(A[i])))); n += 1
    end
    return m, n
end

# Fill every patch-length vector and (np×nrepr) matrix field with field-dependent,
# non-uniform values so a swapped/dropped field surfaces as a mismatch.
function fill_patch_state!(x, np, nrepr)
    for (i, f) in enumerate(fieldnames(typeof(x)))
        v = getfield(x, f)
        if v isa Vector{Float64} && length(v) == np
            for p in 1:np
                v[p] = 0.5 + 0.13 * p + 0.017 * i
            end
        elseif v isa Matrix{Float64} && size(v) == (np, nrepr)
            for p in 1:np, k in 1:nrepr
                v[p, k] = 0.2 + 0.07 * p + 0.05 * k + 0.011 * i
            end
        end
    end
    return x
end

function main(backend)
    println("="^64); println("cnveg per-patch summaries — host vs device parity"); println("="^64)
    if backend === nothing; println("  No GPU backend."); return 0; end
    name, _, FT = backend
    @printf("  Backend: %s (%s)\n", name, FT)

    np = 6; nc = 4; ng = 2; nrepr = 2
    mask = trues(np)
    # patches 2,4,6 are crops (itype >= npcropmin) → exercise the reproductive branch
    patch_itype = [1, 17, 1, 18, 1, 19]
    npcropmin = 17
    nfail = 0; ncmp = 0

    # ---- carbon state ----
    csH = CLM.CNVegCarbonStateData{Float64}(); CLM.cnveg_carbon_state_init!(csH, np, nc, ng; nrepr=nrepr)
    fill_patch_state!(csH, np, nrepr)
    csD = mf(deepcopy(csH))
    csD.dispvegc_patch isa Metal.MtlArray || (println("  BLOCKED: cs adapt"); return 2)
    kw = (; use_crop=true, patch_itype=patch_itype, npcropmin=npcropmin, nrepr=nrepr)
    CLM.cnveg_carbon_state_summary!(csH, mask, 1:np; kw...)
    CLM.cnveg_carbon_state_summary!(csD, mf(mask), 1:np; kw...); Metal.synchronize()
    for f in (:dispvegc_patch, :storvegc_patch, :totvegc_patch, :totc_patch, :woodc_patch)
        r, n = reldiff(getfield(csH, f), getfield(csD, f)); ncmp += n
        ok = r < 1f-3; ok || (nfail += 1)
        @printf("  [%s] C.%-18s rel=%.2e over %d\n", ok ? "PASS" : "FAIL", f, r, n)
    end

    # ---- nitrogen state ----
    nsH = CLM.CNVegNitrogenStateData{Float64}(); CLM.cnveg_nitrogen_state_init!(nsH, np, nc, ng; nrepr=nrepr)
    fill_patch_state!(nsH, np, nrepr)
    nsD = mf(deepcopy(nsH))
    CLM.cnveg_nitrogen_state_summary!(nsH, mask, 1:np; kw...)
    CLM.cnveg_nitrogen_state_summary!(nsD, mf(mask), 1:np; kw...); Metal.synchronize()
    for f in (:dispvegn_patch, :storvegn_patch, :totvegn_patch, :totn_patch)
        r, n = reldiff(getfield(nsH, f), getfield(nsD, f)); ncmp += n
        ok = r < 1f-3; ok || (nfail += 1)
        @printf("  [%s] N.%-18s rel=%.2e over %d\n", ok ? "PASS" : "FAIL", f, r, n)
    end

    # ---- carbon flux ----
    cfH = CLM.CNVegCarbonFluxData{Float64}()
    CLM.cnveg_carbon_flux_init!(cfH, np, nc, ng; nrepr=nrepr)
    fill_patch_state!(cfH, np, nrepr)
    cfD = mf(deepcopy(cfH))
    kwf = (; use_crop=true, use_fun=true, carbon_resp_opt=1,
           patch_itype=patch_itype, npcropmin=npcropmin, nrepr=nrepr)
    CLM.cnveg_carbon_flux_summary!(cfH, mask, 1:np; kwf...)
    CLM.cnveg_carbon_flux_summary!(cfD, mf(mask), 1:np; kwf...); Metal.synchronize()
    for f in (:mr_patch, :current_gr_patch, :transfer_gr_patch, :storage_gr_patch, :gr_patch,
              :ar_patch, :gpp_patch, :npp_patch, :rr_patch, :agnpp_patch, :bgnpp_patch,
              :litfall_patch, :fire_closs_patch, :frootc_alloc_patch, :frootc_loss_patch,
              :leafc_alloc_patch, :leafc_loss_patch, :woodc_alloc_patch, :woodc_loss_patch)
        r, n = reldiff(getfield(cfH, f), getfield(cfD, f)); ncmp += n
        ok = r < 1f-3; ok || (nfail += 1)
        @printf("  [%s] F.%-20s rel=%.2e over %d\n", ok ? "PASS" : "FAIL", f, r, n)
    end

    println()
    println(nfail == 0 ? "  cnveg SUMMARIES MATCH host on $name over $ncmp finite outputs" : "  DIVERGENCE ($nfail).")
    return nfail == 0 ? 0 : 1
end
exit(main(detect_backend()))
