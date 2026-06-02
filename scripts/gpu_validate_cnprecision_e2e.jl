# ==========================================================================
# gpu_validate_cnprecision_e2e.jl — end-to-end Metal parity for the WHOLE
# cn_precision_control! driver (the per-patch precision-truncation cascade over
# ~40 C/N pool pairs + crop reproductive pools + optional c13/c14 isotopes).
#
# cn_precision_control! is a sequence of truncate_*! helpers, each a fully
# independent per-filter-entry loop (thread fp owns patch p = filter[fp]). Each
# helper's arithmetic is a KernelAbstractions kernel (CPU + Metal). The host-side
# error() pre-scan and the dynamic filter_truncatep build are byte-identical
# bookkeeping kept off the device.
#
# This harness mirrors test/test_cn_precision_control.jl's multi-patch setup,
# runs cn_precision_control! on the CPU (Float64 reference), adapts the state to
# Metal at Float32 (the MetalF32 downcast adaptor below), runs the SAME call on
# the device, and compares every mutated output with reldiff.
#
# Exercises three configurations: base (no isotopes), crop pools (use_crop), and
# c13+c14 isotopes (use_c13/use_c14).
#
#   julia --project=scripts scripts/gpu_validate_cnprecision_e2e.jl
# ==========================================================================

using CLM
using Printf
import Metal
include(joinpath(@__DIR__, "gpu_backends.jl"))

function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(A[i] - B[i]) / (1.0 + max(abs(A[i]), abs(B[i]))))
    end
    return m
end

# --------------------------------------------------------------------------
# MetalF32 adaptor: the CN state structs are parametric ({FT,V,M}); we build the
# host reference at Float32 then Adapt.adapt(MtlArray, ·) to move every array
# field onto the device (Metal has no Float64, so the host struct is Float32).
# --------------------------------------------------------------------------
to_metal(x) = CLM.Adapt.adapt(Metal.MtlArray, x)

# Build a minimal CNVegCarbonStateData with np patches at host precision FT.
function make_cs(::Type{FT}, np; nrepr=1, val=FT(0.0)) where {FT}
    cs = CLM.CNVegCarbonStateData{FT}()
    cs.leafc_patch              = fill(val, np)
    cs.leafc_storage_patch      = fill(val, np)
    cs.leafc_xfer_patch         = fill(val, np)
    cs.frootc_patch             = fill(val, np)
    cs.frootc_storage_patch     = fill(val, np)
    cs.frootc_xfer_patch        = fill(val, np)
    cs.livestemc_patch          = fill(val, np)
    cs.livestemc_storage_patch  = fill(val, np)
    cs.livestemc_xfer_patch     = fill(val, np)
    cs.deadstemc_patch          = fill(val, np)
    cs.deadstemc_storage_patch  = fill(val, np)
    cs.deadstemc_xfer_patch     = fill(val, np)
    cs.livecrootc_patch         = fill(val, np)
    cs.livecrootc_storage_patch = fill(val, np)
    cs.livecrootc_xfer_patch    = fill(val, np)
    cs.deadcrootc_patch         = fill(val, np)
    cs.deadcrootc_storage_patch = fill(val, np)
    cs.deadcrootc_xfer_patch    = fill(val, np)
    cs.gresp_storage_patch      = fill(val, np)
    cs.gresp_xfer_patch         = fill(val, np)
    cs.cpool_patch              = fill(val, np)
    cs.xsmrpool_patch           = fill(val, np)
    cs.ctrunc_patch             = fill(FT(0.0), np)
    cs.cropseedc_deficit_patch  = fill(val, np)
    cs.reproductivec_patch          = fill(val, np, nrepr)
    cs.reproductivec_storage_patch  = fill(val, np, nrepr)
    cs.reproductivec_xfer_patch     = fill(val, np, nrepr)
    return cs
end

function make_ns(::Type{FT}, np; nrepr=1, val=FT(0.0)) where {FT}
    ns = CLM.CNVegNitrogenStateData{FT}()
    ns.leafn_patch              = fill(val, np)
    ns.leafn_storage_patch      = fill(val, np)
    ns.leafn_xfer_patch         = fill(val, np)
    ns.frootn_patch             = fill(val, np)
    ns.frootn_storage_patch     = fill(val, np)
    ns.frootn_xfer_patch        = fill(val, np)
    ns.livestemn_patch          = fill(val, np)
    ns.livestemn_storage_patch  = fill(val, np)
    ns.livestemn_xfer_patch     = fill(val, np)
    ns.deadstemn_patch          = fill(val, np)
    ns.deadstemn_storage_patch  = fill(val, np)
    ns.deadstemn_xfer_patch     = fill(val, np)
    ns.livecrootn_patch         = fill(val, np)
    ns.livecrootn_storage_patch = fill(val, np)
    ns.livecrootn_xfer_patch    = fill(val, np)
    ns.deadcrootn_patch         = fill(val, np)
    ns.deadcrootn_storage_patch = fill(val, np)
    ns.deadcrootn_xfer_patch    = fill(val, np)
    ns.retransn_patch           = fill(val, np)
    ns.npool_patch              = fill(val, np)
    ns.ntrunc_patch             = fill(FT(0.0), np)
    ns.cropseedn_deficit_patch  = fill(val, np)
    ns.reproductiven_patch          = fill(val, np, nrepr)
    ns.reproductiven_storage_patch  = fill(val, np, nrepr)
    ns.reproductiven_xfer_patch     = fill(val, np, nrepr)
    return ns
end

# Fields of the carbon state that cn_precision_control! mutates (for comparison).
const CS_FIELDS = (:leafc_patch, :leafc_storage_patch, :leafc_xfer_patch,
    :frootc_patch, :frootc_storage_patch, :frootc_xfer_patch,
    :livestemc_patch, :livestemc_storage_patch, :livestemc_xfer_patch,
    :deadstemc_patch, :deadstemc_storage_patch, :deadstemc_xfer_patch,
    :livecrootc_patch, :livecrootc_storage_patch, :livecrootc_xfer_patch,
    :deadcrootc_patch, :deadcrootc_storage_patch, :deadcrootc_xfer_patch,
    :gresp_storage_patch, :gresp_xfer_patch, :cpool_patch, :xsmrpool_patch,
    :ctrunc_patch, :cropseedc_deficit_patch,
    :reproductivec_patch, :reproductivec_storage_patch, :reproductivec_xfer_patch)

const NS_FIELDS = (:leafn_patch, :leafn_storage_patch, :leafn_xfer_patch,
    :frootn_patch, :frootn_storage_patch, :frootn_xfer_patch,
    :livestemn_patch, :livestemn_storage_patch, :livestemn_xfer_patch,
    :deadstemn_patch, :deadstemn_storage_patch, :deadstemn_xfer_patch,
    :livecrootn_patch, :livecrootn_storage_patch, :livecrootn_xfer_patch,
    :deadcrootn_patch, :deadcrootn_storage_patch, :deadcrootn_xfer_patch,
    :retransn_patch, :npool_patch, :ntrunc_patch, :cropseedn_deficit_patch,
    :reproductiven_patch, :reproductiven_storage_patch, :reproductiven_xfer_patch)

# Compare every mutated field of a CPU reference cs/ns against the device cs/ns.
function compare(label, csH, nsH, csD, nsD; c13H=nothing, c13D=nothing,
                 c14H=nothing, c14D=nothing)
    nfail = 0; worst = 0.0
    function chk(nm, a, b)
        (isempty(a) && isempty(b)) && return
        d = reldiff(a, b); ok = d < 1f-2
        worst = max(worst, d)
        ok || (nfail += 1)
        @printf("    [%s] %-32s rel = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
    end
    println("  -- $label --")
    for f in CS_FIELDS
        chk("cs.$f", getfield(csH, f), getfield(csD, f))
    end
    for f in NS_FIELDS
        chk("ns.$f", getfield(nsH, f), getfield(nsD, f))
    end
    if c13H !== nothing
        for f in CS_FIELDS
            chk("c13.$f", getfield(c13H, f), getfield(c13D, f))
        end
    end
    if c14H !== nothing
        for f in CS_FIELDS
            chk("c14.$f", getfield(c14H, f), getfield(c14D, f))
        end
    end
    @printf("    worst reldiff = %.3e\n\n", worst)
    return nfail
end

function main(backend)
    println("=" ^ 70)
    println("END-TO-END Metal parity for cn_precision_control! (whole driver)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU driver exercised by the suite).")
        return 0
    end
    name, _, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    tiny = 5.0e-9
    nfail = 0

    # ----------------------------------------------------------------------
    # Case 1: base path (mixed tiny/large values, no isotopes, no crop).
    # ----------------------------------------------------------------------
    let np = 4
        csH = make_cs(Float32, np; val=1.0f0)
        nsH = make_ns(Float32, np; val=1.0f0)
        # Make some pools tiny so truncation actually fires.
        csH.leafc_patch  .= Float32[1.0, tiny, -tiny, 1.0]
        nsH.leafn_patch  .= Float32[0.5, tiny,  tiny, 0.5]
        csH.cpool_patch  .= Float32[tiny, 1.0, tiny, 2.0]
        csH.gresp_xfer_patch .= Float32[tiny, tiny, 1.0, 1.0]
        nsH.npool_patch  .= Float32[tiny, 1.0, tiny, 1.0]
        filter = [1, 2, 3, 4]; itype = [1, 1, 1, 1]

        csD = to_metal(deepcopy(csH)); nsD = to_metal(deepcopy(nsH))
        fD = Metal.MtlArray(filter); iD = Metal.MtlArray(itype)

        CLM.cn_precision_control!(csH, nsH, filter, itype)
        CLM.cn_precision_control!(csD, nsD, fD, iD)
        nfail += compare("Case 1: base path", csH, nsH, csD, nsD)
    end

    # ----------------------------------------------------------------------
    # Case 2: crop pools (use_crop=true; crop-only pools on crop patches).
    # ----------------------------------------------------------------------
    let np = 3
        csH = make_cs(Float32, np; val=Float32(tiny))
        nsH = make_ns(Float32, np; val=Float32(tiny))
        filter = [1, 2, 3]; itype = [1, 15, 16]  # patches 2,3 are crops

        csD = to_metal(deepcopy(csH)); nsD = to_metal(deepcopy(nsH))
        fD = Metal.MtlArray(filter); iD = Metal.MtlArray(itype)

        CLM.cn_precision_control!(csH, nsH, filter, itype; use_crop=true)
        CLM.cn_precision_control!(csD, nsD, fD, iD; use_crop=true)
        nfail += compare("Case 2: crop pools", csH, nsH, csD, nsD)
    end

    # ----------------------------------------------------------------------
    # Case 3: c13 + c14 isotopes (use_c13/use_c14=true).
    # ----------------------------------------------------------------------
    let np = 2
        csH   = make_cs(Float32, np; val=Float32(tiny))
        nsH   = make_ns(Float32, np; val=Float32(tiny))
        c13H  = make_cs(Float32, np; val=Float32(tiny * 0.01))
        c14H  = make_cs(Float32, np; val=Float32(tiny * 0.001))
        filter = [1, 2]; itype = [1, 1]

        csD  = to_metal(deepcopy(csH));  nsD  = to_metal(deepcopy(nsH))
        c13D = to_metal(deepcopy(c13H)); c14D = to_metal(deepcopy(c14H))
        fD = Metal.MtlArray(filter); iD = Metal.MtlArray(itype)

        CLM.cn_precision_control!(csH, nsH, filter, itype;
            c13cs=c13H, c14cs=c14H, use_c13=true, use_c14=true)
        CLM.cn_precision_control!(csD, nsD, fD, iD;
            c13cs=c13D, c14cs=c14D, use_c13=true, use_c14=true)
        nfail += compare("Case 3: c13+c14 isotopes", csH, nsH, csD, nsD;
            c13H=c13H, c13D=c13D, c14H=c14H, c14D=c14D)
    end

    # ----------------------------------------------------------------------
    # Case 4: matrix-CN guardrail path (use_matrixcn=true + use_nguardrail=true).
    # This is the ONLY kernel branch carrying numeric literals — the
    # `-ccrit*T(1.0e+6)` / `-ncrit*T(1.0e+6)` window check — so it must be
    # exercised on Metal. Values are placed inside the (-crit*1e6, crit) window
    # so truncation actually fires under the matrixcn predicate.
    # ----------------------------------------------------------------------
    let np = 4
        csH = make_cs(Float32, np; val=1.0f0)
        nsH = make_ns(Float32, np; val=1.0f0)
        # mix: tiny+, tiny-, mid-negative (inside the 1e6 window), large (kept)
        csH.leafc_patch  .= Float32[tiny, -tiny, -1.0e-4, 1.0]
        nsH.leafn_patch  .= Float32[tiny, -tiny, -1.0e-4, 0.5]
        csH.cpool_patch  .= Float32[-tiny, tiny, -1.0e-5, 2.0]
        nsH.npool_patch  .= Float32[-tiny, tiny, -1.0e-5, 1.0]
        filter = [1, 2, 3, 4]; itype = [1, 1, 1, 1]

        csD = to_metal(deepcopy(csH)); nsD = to_metal(deepcopy(nsH))
        fD = Metal.MtlArray(filter); iD = Metal.MtlArray(itype)

        CLM.cn_precision_control!(csH, nsH, filter, itype;
            use_matrixcn=true, use_nguardrail=true)
        CLM.cn_precision_control!(csD, nsD, fD, iD;
            use_matrixcn=true, use_nguardrail=true)
        nfail += compare("Case 4: matrixcn guardrail", csH, nsH, csD, nsD)
    end

    println(nfail == 0 ?
        "  WHOLE cn_precision_control! MATCHES CPU ON $name ($FT)" :
        "  DIVERGENCE — investigate ($nfail field(s) failed).")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
