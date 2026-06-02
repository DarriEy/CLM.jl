# ==========================================================================
# gpu_validate_cstateupdate0_dyn_e2e.jl — end-to-end Metal parity for the two
# small sibling C-state updates: c_state_update0! (photosynthesis cpool update,
# one thread per patch) and c_state_update_dyn_patch! (dyn_cnbal_patch litter/CWD
# input, one thread per column + seed carbon, one thread per gridcell).
#
# Builds a small multi-column / multi-patch / multi-gridcell instance mirroring
# test/test_c_state_update1.jl, runs each function on the CPU, converts every state
# struct to Float32 + adapts to Metal, runs the SAME calls on the device, and
# compares the mutated outputs field-by-field. The dyn scenario exercises the
# per-column j/i litter loops (litter pools + CWD) and the per-gridcell seed loop;
# update0 exercises the per-patch cpool accumulation.
#
#   julia --project=scripts scripts/gpu_validate_cstateupdate0_dyn_e2e.jl
#
# The CN *_init! routines allocate Float64 (fill(NaN, …)) regardless of the struct
# type param, so we build at Float64 and down-convert array fields to Float32 for
# the device snapshot (Metal has no Float64).
# ==========================================================================

using CLM
using Printf
import Metal   # MtlArray is the Adapt adaptor type for the device-view structs
include(joinpath(@__DIR__, "gpu_backends.jl"))

# NaN-aware mixed abs/rel diff; asserts the CPU reference is finite so a both-NaN
# false PASS can't slip through (fields *_init! leaves as NaN stay NaN on both).
function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(A[i] - B[i]) / (1.0 + max(abs(A[i]), abs(B[i]))))
    end
    return m
end
cpu_has_finite(a) = any(isfinite, Array(a))

# Float32-down-converting Metal adaptor. The CN state structs are concretely
# Float64-typed (default ctor pins ::Vector{Float64}/::Matrix{Float64}) and the
# *_init! routines fill Float64, so we can't setfield! a Float32 array. Instead we
# adapt with a custom storage rule that down-converts float arrays to Float32 as it
# reconstructs the struct (integer/bool arrays move as-is). @adapt_structure rebuilds
# each struct positionally, inferring the new {Float32,…} params from the adapted fields.
struct MetalF32 end
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{<:AbstractFloat}) = Metal.MtlArray(Float32.(x))
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{<:Integer})       = Metal.MtlArray(x)
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{Bool})            = Metal.MtlArray(x)

const NC = 4
const NP = 6
const NG = 2
const NLEVDECOMP = 1
const NDECOMP_POOLS = 7
const NCASCADE = 5
const NREPR = 1
const I_LITR_MIN = 1
const I_LITR_MAX = 3
const I_CWD = 4

# Build the CPU reference state (Float64). Returns a NamedTuple of the structs +
# the masks update0! / dyn_patch! need.
function build()
    cs_veg = CLM.CNVegCarbonStateData()
    CLM.cnveg_carbon_state_init!(cs_veg, NP, NC, NG; nrepr=NREPR)
    cf_veg = CLM.CNVegCarbonFluxData()
    CLM.cnveg_carbon_flux_init!(cf_veg, NP, NC, NG; nrepr=NREPR,
                                nlevdecomp_full=NLEVDECOMP, ndecomp_pools=NDECOMP_POOLS)
    cs_soil = CLM.SoilBiogeochemCarbonStateData()
    CLM.soil_bgc_carbon_state_init!(cs_soil, NC, NG, NLEVDECOMP, NDECOMP_POOLS)

    # --- update0! inputs: cpool + photosynthesis fluxes ---
    cs_veg.cpool_patch             .= 100.0
    cf_veg.psnsun_to_cpool_patch   .= 1.0e-6
    cf_veg.psnshade_to_cpool_patch .= 0.5e-6

    # --- dyn_patch! inputs: seed carbon + decomp pools + dyn fluxes ---
    cs_veg.seedc_grc               .= 1000.0
    cs_soil.decomp_cpools_vr_col   .= 100.0
    for i in I_LITR_MIN:I_LITR_MAX
        cf_veg.dwt_frootc_to_litr_c_col[:, :, i] .= 1.0e-6
    end
    cf_veg.dwt_livecrootc_to_cwdc_col .= 2.0e-6
    cf_veg.dwt_deadcrootc_to_cwdc_col .= 3.0e-6
    cf_veg.dwt_seedc_to_leaf_grc      .= 0.5e-6
    cf_veg.dwt_seedc_to_deadstem_grc  .= 0.3e-6

    mask_soilp                = trues(NP)
    mask_soilc_with_inactive  = trues(NC)

    return (; cs_veg, cf_veg, cs_soil, mask_soilp, mask_soilc_with_inactive)
end

run_csu0!(S, m_p, dt) =
    CLM.c_state_update0!(S.cs_veg, S.cf_veg;
        mask_soilp=m_p, bounds_patch=1:NP, dt=dt)

run_dyn!(S, m_c, dt) =
    CLM.c_state_update_dyn_patch!(S.cs_veg, S.cf_veg, S.cs_soil;
        mask_soilc_with_inactive=m_c, bounds_col=1:NC, bounds_grc=1:NG,
        nlevdecomp=NLEVDECOMP, i_litr_min=I_LITR_MIN, i_litr_max=I_LITR_MAX,
        i_cwd=I_CWD, dt=dt)

function main(backend)
    println("=" ^ 70)
    println("END-TO-END Metal parity for c_state_update0! + c_state_update_dyn_patch!")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU path exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    H = build()    # CPU reference (Float64)
    B = build()    # source for the device snapshot

    # Adapt the device snapshot to Metal, down-converting float arrays to Float32.
    ad(x) = CLM.Adapt.adapt(MetalF32(), x)
    cs_d  = ad(B.cs_veg); cf_d = ad(B.cf_veg); css_d = ad(B.cs_soil)
    Sd = (; cs_veg=cs_d, cf_veg=cf_d, cs_soil=css_d)

    if !(cs_d.cpool_patch isa Metal.MtlArray)
        println("  BLOCKED: a state struct did not move to the device under adapt.")
        return 2
    end

    dt = FT(1800.0)
    dmask(m) = to(collect(Bool, m))

    # CPU
    run_csu0!(H, H.mask_soilp, 1800.0)
    run_dyn!(H, H.mask_soilc_with_inactive, 1800.0)
    # Device
    run_csu0!(Sd, dmask(B.mask_soilp), dt)
    run_dyn!(Sd, dmask(B.mask_soilc_with_inactive), dt)

    checks = [
        ("cpool (update0)",        H.cs_veg.cpool_patch,           Sd.cs_veg.cpool_patch),
        ("seedc_grc (dyn)",        H.cs_veg.seedc_grc,             Sd.cs_veg.seedc_grc),
        ("decomp_cpools_vr (dyn)", H.cs_soil.decomp_cpools_vr_col, Sd.cs_soil.decomp_cpools_vr_col),
    ]
    nfail = 0
    for (nm, a, b) in checks
        if !cpu_has_finite(a)
            @printf("  [WARN] %-24s CPU reference is all-NaN/Inf — skipping (no parity signal)\n", nm)
            continue
        end
        d = reldiff(a, b); ok = d < 1f-3
        @printf("  [%s] %-24s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end
    println()
    println(nfail == 0 ? "  update0! + dyn_patch! MATCH CPU ON $name ($FT) ✓" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
