# ==========================================================================
# gpu_validate_cstateupdate3_e2e.jl — end-to-end Metal parity for the WHOLE
# c_state_update3! BGC driver (fire C fluxes; part of the Phase B C/N
# state-update cascade).
#
# Builds a small multi-column / multi-patch instance mirroring
# test/test_c_state_update3.jl, runs c_state_update3! on the CPU, converts every
# state struct to Float32 + adapts to Metal, runs the SAME call on the device, and
# compares the mutated outputs field-by-field. The scenario exercises the fused
# column kernel (fire mortality -> cwd/litter additions, then per-pool fire losses)
# and the patch kernel (gresp + displayed/storage/transfer pool fire fluxes,
# including the live-to-dead stem/croot transfers).
#
#   julia --project=scripts scripts/gpu_validate_cstateupdate3_e2e.jl
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
# false PASS can't slip through.
function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(A[i] - B[i]) / (1.0 + max(abs(A[i]), abs(B[i]))))
    end
    return m
end
cpu_has_finite(a) = any(isfinite, Array(a))

# Float32-down-converting Metal adaptor (see gpu_validate_cstateupdate1_e2e.jl).
struct MetalF32 end
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{<:AbstractFloat}) = Metal.MtlArray(Float32.(x))
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{<:Integer})       = Metal.MtlArray(x)
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{Bool})            = Metal.MtlArray(x)

const NC = 4
const NP = 6
const NG = 2
const NLEVDECOMP = 1
const NDECOMP_POOLS = 7
const NREPR = 1
const I_LITR_MIN = 1
const I_LITR_MAX = 3
const I_CWD = 4

# Build the CPU reference state (Float64). Mirrors test/test_c_state_update3.jl.
function build()
    cs_veg = CLM.CNVegCarbonStateData()
    CLM.cnveg_carbon_state_init!(cs_veg, NP, NC, NG; nrepr=NREPR)
    cs_veg.leafc_patch             .= 50.0
    cs_veg.leafc_storage_patch     .= 10.0
    cs_veg.leafc_xfer_patch        .= 5.0
    cs_veg.frootc_patch            .= 40.0
    cs_veg.frootc_storage_patch    .= 8.0
    cs_veg.frootc_xfer_patch       .= 4.0
    cs_veg.livestemc_patch         .= 30.0
    cs_veg.livestemc_storage_patch .= 6.0
    cs_veg.livestemc_xfer_patch    .= 3.0
    cs_veg.deadstemc_patch         .= 200.0
    cs_veg.deadstemc_storage_patch .= 4.0
    cs_veg.deadstemc_xfer_patch    .= 2.0
    cs_veg.livecrootc_patch        .= 20.0
    cs_veg.livecrootc_storage_patch .= 4.0
    cs_veg.livecrootc_xfer_patch   .= 2.0
    cs_veg.deadcrootc_patch        .= 100.0
    cs_veg.deadcrootc_storage_patch .= 2.0
    cs_veg.deadcrootc_xfer_patch   .= 1.0
    cs_veg.gresp_storage_patch     .= 5.0
    cs_veg.gresp_xfer_patch        .= 2.0

    cf_veg = CLM.CNVegCarbonFluxData()
    CLM.cnveg_carbon_flux_init!(cf_veg, NP, NC, NG; nrepr=NREPR,
                                nlevdecomp_full=NLEVDECOMP, ndecomp_pools=NDECOMP_POOLS)

    # Fire mortality fluxes — patch level to combustion
    cf_veg.m_leafc_to_fire_patch              .= 1.0e-6
    cf_veg.m_frootc_to_fire_patch             .= 0.8e-6
    cf_veg.m_livestemc_to_fire_patch          .= 0.5e-6
    cf_veg.m_deadstemc_to_fire_patch          .= 0.3e-6
    cf_veg.m_livecrootc_to_fire_patch         .= 0.4e-6
    cf_veg.m_deadcrootc_to_fire_patch         .= 0.2e-6
    cf_veg.m_leafc_storage_to_fire_patch      .= 0.1e-6
    cf_veg.m_frootc_storage_to_fire_patch     .= 0.08e-6
    cf_veg.m_livestemc_storage_to_fire_patch  .= 0.05e-6
    cf_veg.m_deadstemc_storage_to_fire_patch  .= 0.03e-6
    cf_veg.m_livecrootc_storage_to_fire_patch .= 0.04e-6
    cf_veg.m_deadcrootc_storage_to_fire_patch .= 0.02e-6
    cf_veg.m_leafc_xfer_to_fire_patch         .= 0.05e-6
    cf_veg.m_frootc_xfer_to_fire_patch        .= 0.04e-6
    cf_veg.m_livestemc_xfer_to_fire_patch     .= 0.03e-6
    cf_veg.m_deadstemc_xfer_to_fire_patch     .= 0.02e-6
    cf_veg.m_livecrootc_xfer_to_fire_patch    .= 0.015e-6
    cf_veg.m_deadcrootc_xfer_to_fire_patch    .= 0.01e-6
    cf_veg.m_gresp_storage_to_fire_patch      .= 0.01e-6
    cf_veg.m_gresp_xfer_to_fire_patch         .= 0.005e-6

    # Fire mortality fluxes — patch level to litter
    cf_veg.m_leafc_to_litter_fire_patch              .= 0.5e-6
    cf_veg.m_frootc_to_litter_fire_patch             .= 0.4e-6
    cf_veg.m_livestemc_to_litter_fire_patch          .= 0.3e-6
    cf_veg.m_deadstemc_to_litter_fire_patch          .= 0.15e-6
    cf_veg.m_livecrootc_to_litter_fire_patch         .= 0.2e-6
    cf_veg.m_deadcrootc_to_litter_fire_patch         .= 0.1e-6
    cf_veg.m_leafc_storage_to_litter_fire_patch      .= 0.05e-6
    cf_veg.m_frootc_storage_to_litter_fire_patch     .= 0.04e-6
    cf_veg.m_livestemc_storage_to_litter_fire_patch  .= 0.03e-6
    cf_veg.m_deadstemc_storage_to_litter_fire_patch  .= 0.015e-6
    cf_veg.m_livecrootc_storage_to_litter_fire_patch .= 0.02e-6
    cf_veg.m_deadcrootc_storage_to_litter_fire_patch .= 0.01e-6
    cf_veg.m_leafc_xfer_to_litter_fire_patch         .= 0.025e-6
    cf_veg.m_frootc_xfer_to_litter_fire_patch        .= 0.02e-6
    cf_veg.m_livestemc_xfer_to_litter_fire_patch     .= 0.015e-6
    cf_veg.m_deadstemc_xfer_to_litter_fire_patch     .= 0.01e-6
    cf_veg.m_livecrootc_xfer_to_litter_fire_patch    .= 0.008e-6
    cf_veg.m_deadcrootc_xfer_to_litter_fire_patch    .= 0.005e-6
    cf_veg.m_gresp_storage_to_litter_fire_patch      .= 0.005e-6
    cf_veg.m_gresp_xfer_to_litter_fire_patch         .= 0.003e-6

    # Live-to-dead fire transfers
    cf_veg.m_livestemc_to_deadstemc_fire_patch       .= 0.2e-6
    cf_veg.m_livecrootc_to_deadcrootc_fire_patch     .= 0.15e-6

    # Column-level fire fluxes
    cf_veg.fire_mortality_c_to_cwdc_col   .= 1.0e-6
    cf_veg.m_c_to_litr_fire_col           .= 2.0e-6
    cf_veg.m_decomp_cpools_to_fire_vr_col .= 0.5e-6

    cs_soil = CLM.SoilBiogeochemCarbonStateData()
    CLM.soil_bgc_carbon_state_init!(cs_soil, NC, NG, NLEVDECOMP, NDECOMP_POOLS)
    cs_soil.decomp_cpools_vr_col .= 100.0

    mask_soilc = trues(NC)
    mask_soilp = trues(NP)

    return (; cs_veg, cf_veg, cs_soil, mask_soilc, mask_soilp)
end

run_csu3!(S, m_c, m_p, dt) =
    CLM.c_state_update3!(S.cs_veg, S.cf_veg, S.cs_soil;
        mask_soilc=m_c, mask_soilp=m_p, bounds_col=1:NC, bounds_patch=1:NP,
        nlevdecomp=NLEVDECOMP, ndecomp_pools=NDECOMP_POOLS,
        i_litr_min=I_LITR_MIN, i_litr_max=I_LITR_MAX, i_cwd=I_CWD,
        use_matrixcn=false, use_soil_matrixcn=false, dt=dt)

function main(backend)
    println("=" ^ 70)
    println("END-TO-END Metal parity for c_state_update3! (BGC fire C-state)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU path exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    H = build()    # CPU reference (Float64)
    B = build()    # source for the device snapshot

    ad(x) = CLM.Adapt.adapt(MetalF32(), x)
    cs_d  = ad(B.cs_veg); cf_d = ad(B.cf_veg); css_d = ad(B.cs_soil)
    Sd = (; cs_veg=cs_d, cf_veg=cf_d, cs_soil=css_d)

    if !(cs_d.leafc_patch isa Metal.MtlArray)
        println("  BLOCKED: a state struct did not move to the device under adapt.")
        return 2
    end

    dt = FT(1800.0)
    dmask(m) = to(collect(Bool, m))

    # CPU
    run_csu3!(H, H.mask_soilc, H.mask_soilp, 1800.0)
    # Device
    run_csu3!(Sd, dmask(B.mask_soilc), dmask(B.mask_soilp), dt)

    checks = [
        ("gresp_storage",    H.cs_veg.gresp_storage_patch,    Sd.cs_veg.gresp_storage_patch),
        ("gresp_xfer",       H.cs_veg.gresp_xfer_patch,       Sd.cs_veg.gresp_xfer_patch),
        ("leafc",            H.cs_veg.leafc_patch,            Sd.cs_veg.leafc_patch),
        ("frootc",           H.cs_veg.frootc_patch,           Sd.cs_veg.frootc_patch),
        ("livestemc",        H.cs_veg.livestemc_patch,        Sd.cs_veg.livestemc_patch),
        ("deadstemc",        H.cs_veg.deadstemc_patch,        Sd.cs_veg.deadstemc_patch),
        ("livecrootc",       H.cs_veg.livecrootc_patch,       Sd.cs_veg.livecrootc_patch),
        ("deadcrootc",       H.cs_veg.deadcrootc_patch,       Sd.cs_veg.deadcrootc_patch),
        ("leafc_storage",    H.cs_veg.leafc_storage_patch,    Sd.cs_veg.leafc_storage_patch),
        ("frootc_storage",   H.cs_veg.frootc_storage_patch,   Sd.cs_veg.frootc_storage_patch),
        ("livestemc_storage",H.cs_veg.livestemc_storage_patch,Sd.cs_veg.livestemc_storage_patch),
        ("deadstemc_storage",H.cs_veg.deadstemc_storage_patch,Sd.cs_veg.deadstemc_storage_patch),
        ("livecrootc_storage",H.cs_veg.livecrootc_storage_patch,Sd.cs_veg.livecrootc_storage_patch),
        ("deadcrootc_storage",H.cs_veg.deadcrootc_storage_patch,Sd.cs_veg.deadcrootc_storage_patch),
        ("leafc_xfer",       H.cs_veg.leafc_xfer_patch,       Sd.cs_veg.leafc_xfer_patch),
        ("frootc_xfer",      H.cs_veg.frootc_xfer_patch,      Sd.cs_veg.frootc_xfer_patch),
        ("livestemc_xfer",   H.cs_veg.livestemc_xfer_patch,   Sd.cs_veg.livestemc_xfer_patch),
        ("deadstemc_xfer",   H.cs_veg.deadstemc_xfer_patch,   Sd.cs_veg.deadstemc_xfer_patch),
        ("livecrootc_xfer",  H.cs_veg.livecrootc_xfer_patch,  Sd.cs_veg.livecrootc_xfer_patch),
        ("deadcrootc_xfer",  H.cs_veg.deadcrootc_xfer_patch,  Sd.cs_veg.deadcrootc_xfer_patch),
        ("decomp_cpools_vr", H.cs_soil.decomp_cpools_vr_col,  Sd.cs_soil.decomp_cpools_vr_col),
    ]
    nfail = 0
    for (nm, a, b) in checks
        if !cpu_has_finite(a)
            @printf("  [WARN] %-20s CPU reference is all-NaN/Inf — skipping (no parity signal)\n", nm)
            continue
        end
        d = reldiff(a, b); ok = d < 1f-3
        @printf("  [%s] %-20s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end
    println()
    println(nfail == 0 ? "  WHOLE c_state_update3! MATCHES CPU ON $name ($FT) ✓" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
