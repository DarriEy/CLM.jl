# ==========================================================================
# gpu_validate_cstateupdate2_e2e.jl — end-to-end GPU parity for the WHOLE
# c_state_update2! / c_state_update2h! / c_state_update2g! BGC mortality drivers
# (gap-phase mortality, harvest, gross unrepresented landcover change C fluxes).
#
# Builds a small multi-column / multi-patch instance mirroring
# test/test_c_state_update2.jl, runs each of the three functions on the CPU,
# converts every state struct to Float32 + adapts to the GPU, runs the SAME calls
# on the device, and compares the mutated outputs field-by-field. The scenario
# exercises woody + non-woody + crop patches in one build (these mortality
# functions are not ivt/woody-branchy, so all patches take the same path, but the
# multi-patch / multi-column topology stresses the per-patch + per-column kernels).
#
#   julia --project=scripts scripts/gpu_validate_cstateupdate2_e2e.jl
#
# The CN *_init! routines allocate Float64 (fill(NaN, …)) regardless of the struct
# type param, so we build at Float64 and down-convert array fields to Float32 for
# the device snapshot (Metal has no Float64).
# ==========================================================================

using CLM
using Printf
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

# Float32-down-converting Metal adaptor (see c_state_update1 harness for rationale).
struct MetalF32 end
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{<:AbstractFloat}) = device_array_type()(Float32.(x))
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{<:Integer})       = device_array_type()(x)
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{Bool})            = device_array_type()(x)

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
const NPCROPMIN = 17

# Build the CPU reference state (Float64). Mirrors test/test_c_state_update2.jl.
function build()
    cs_veg = CLM.CNVegCarbonStateData()
    CLM.cnveg_carbon_state_init!(cs_veg, NP, NC, NG; nrepr=NREPR)
    cf_veg = CLM.CNVegCarbonFluxData()
    CLM.cnveg_carbon_flux_init!(cf_veg, NP, NC, NG; nrepr=NREPR,
                                nlevdecomp_full=NLEVDECOMP, ndecomp_pools=NDECOMP_POOLS)
    cs_soil = CLM.SoilBiogeochemCarbonStateData()
    CLM.soil_bgc_carbon_state_init!(cs_soil, NC, NG, NLEVDECOMP, NDECOMP_POOLS)

    # --- Carbon-state pools (known nonzero so every flux moves something) ---
    cs_veg.leafc_patch              .= 50.0
    cs_veg.leafc_storage_patch      .= 10.0
    cs_veg.leafc_xfer_patch         .= 5.0
    cs_veg.frootc_patch             .= 40.0
    cs_veg.frootc_storage_patch     .= 8.0
    cs_veg.frootc_xfer_patch        .= 4.0
    cs_veg.livestemc_patch          .= 30.0
    cs_veg.livestemc_storage_patch  .= 6.0
    cs_veg.livestemc_xfer_patch     .= 3.0
    cs_veg.deadstemc_patch          .= 200.0
    cs_veg.deadstemc_storage_patch  .= 4.0
    cs_veg.deadstemc_xfer_patch     .= 2.0
    cs_veg.livecrootc_patch         .= 20.0
    cs_veg.livecrootc_storage_patch .= 4.0
    cs_veg.livecrootc_xfer_patch    .= 2.0
    cs_veg.deadcrootc_patch         .= 100.0
    cs_veg.deadcrootc_storage_patch .= 2.0
    cs_veg.deadcrootc_xfer_patch    .= 1.0
    cs_veg.xsmrpool_patch           .= 15.0
    cs_veg.gresp_storage_patch      .= 5.0
    cs_veg.gresp_xfer_patch         .= 2.0

    # --- Gap mortality fluxes (CStateUpdate2) ---
    cf_veg.m_leafc_to_litter_patch              .= 1.0e-6
    cf_veg.m_frootc_to_litter_patch             .= 0.8e-6
    cf_veg.m_livestemc_to_litter_patch          .= 0.5e-6
    cf_veg.m_deadstemc_to_litter_patch          .= 0.3e-6
    cf_veg.m_livecrootc_to_litter_patch         .= 0.4e-6
    cf_veg.m_deadcrootc_to_litter_patch         .= 0.2e-6
    cf_veg.m_leafc_storage_to_litter_patch      .= 0.1e-6
    cf_veg.m_frootc_storage_to_litter_patch     .= 0.08e-6
    cf_veg.m_livestemc_storage_to_litter_patch  .= 0.05e-6
    cf_veg.m_deadstemc_storage_to_litter_patch  .= 0.03e-6
    cf_veg.m_livecrootc_storage_to_litter_patch .= 0.04e-6
    cf_veg.m_deadcrootc_storage_to_litter_patch .= 0.02e-6
    cf_veg.m_leafc_xfer_to_litter_patch         .= 0.05e-6
    cf_veg.m_frootc_xfer_to_litter_patch        .= 0.04e-6
    cf_veg.m_livestemc_xfer_to_litter_patch     .= 0.03e-6
    cf_veg.m_deadstemc_xfer_to_litter_patch     .= 0.02e-6
    cf_veg.m_livecrootc_xfer_to_litter_patch    .= 0.015e-6
    cf_veg.m_deadcrootc_xfer_to_litter_patch    .= 0.01e-6
    cf_veg.m_gresp_storage_to_litter_patch      .= 0.01e-6
    cf_veg.m_gresp_xfer_to_litter_patch         .= 0.005e-6
    cf_veg.gap_mortality_c_to_litr_c_col        .= 2.0e-6
    cf_veg.gap_mortality_c_to_cwdc_col          .= 1.0e-6

    # --- Harvest mortality fluxes (CStateUpdate2h) ---
    cf_veg.hrv_leafc_to_litter_patch            .= 0.5e-6
    cf_veg.hrv_frootc_to_litter_patch           .= 0.4e-6
    cf_veg.hrv_livestemc_to_litter_patch        .= 0.3e-6
    cf_veg.hrv_livecrootc_to_litter_patch       .= 0.2e-6
    cf_veg.hrv_deadcrootc_to_litter_patch       .= 0.1e-6
    cf_veg.wood_harvestc_patch                  .= 0.6e-6
    cf_veg.hrv_leafc_storage_to_litter_patch    .= 0.05e-6
    cf_veg.hrv_frootc_storage_to_litter_patch   .= 0.04e-6
    cf_veg.hrv_livestemc_storage_to_litter_patch .= 0.03e-6
    cf_veg.hrv_deadstemc_storage_to_litter_patch .= 0.02e-6
    cf_veg.hrv_livecrootc_storage_to_litter_patch .= 0.015e-6
    cf_veg.hrv_deadcrootc_storage_to_litter_patch .= 0.01e-6
    cf_veg.hrv_leafc_xfer_to_litter_patch       .= 0.025e-6
    cf_veg.hrv_frootc_xfer_to_litter_patch      .= 0.02e-6
    cf_veg.hrv_livestemc_xfer_to_litter_patch   .= 0.015e-6
    cf_veg.hrv_deadstemc_xfer_to_litter_patch   .= 0.01e-6
    cf_veg.hrv_livecrootc_xfer_to_litter_patch  .= 0.008e-6
    cf_veg.hrv_deadcrootc_xfer_to_litter_patch  .= 0.005e-6
    cf_veg.hrv_gresp_storage_to_litter_patch    .= 0.005e-6
    cf_veg.hrv_gresp_xfer_to_litter_patch       .= 0.003e-6
    cf_veg.hrv_xsmrpool_to_atm_patch            .= 0.01e-6
    cf_veg.harvest_c_to_litr_c_col              .= 1.5e-6
    cf_veg.harvest_c_to_cwdc_col                .= 0.8e-6

    # --- Gross unrepresented landcover change fluxes (CStateUpdate2g) ---
    cf_veg.gru_leafc_to_litter_patch            .= 0.3e-6
    cf_veg.gru_frootc_to_litter_patch           .= 0.25e-6
    cf_veg.gru_livestemc_to_atm_patch           .= 0.2e-6
    cf_veg.gru_deadstemc_to_atm_patch           .= 0.15e-6
    cf_veg.gru_wood_productc_gain_patch         .= 0.1e-6
    cf_veg.gru_livecrootc_to_litter_patch       .= 0.12e-6
    cf_veg.gru_deadcrootc_to_litter_patch       .= 0.08e-6
    cf_veg.gru_leafc_storage_to_atm_patch       .= 0.03e-6
    cf_veg.gru_frootc_storage_to_atm_patch      .= 0.025e-6
    cf_veg.gru_livestemc_storage_to_atm_patch   .= 0.02e-6
    cf_veg.gru_deadstemc_storage_to_atm_patch   .= 0.015e-6
    cf_veg.gru_livecrootc_storage_to_atm_patch  .= 0.012e-6
    cf_veg.gru_deadcrootc_storage_to_atm_patch  .= 0.008e-6
    cf_veg.gru_leafc_xfer_to_atm_patch          .= 0.015e-6
    cf_veg.gru_frootc_xfer_to_atm_patch         .= 0.012e-6
    cf_veg.gru_livestemc_xfer_to_atm_patch      .= 0.01e-6
    cf_veg.gru_deadstemc_xfer_to_atm_patch      .= 0.008e-6
    cf_veg.gru_livecrootc_xfer_to_atm_patch     .= 0.006e-6
    cf_veg.gru_deadcrootc_xfer_to_atm_patch     .= 0.004e-6
    cf_veg.gru_gresp_storage_to_atm_patch       .= 0.003e-6
    cf_veg.gru_gresp_xfer_to_atm_patch          .= 0.002e-6
    cf_veg.gru_xsmrpool_to_atm_patch            .= 0.005e-6
    cf_veg.gru_c_to_litr_c_col                  .= 1.0e-6
    cf_veg.gru_c_to_cwdc_col                    .= 0.5e-6

    # --- Soil biogeochem carbon state (column-kernel target; RMW accumulation) ---
    cs_soil.decomp_cpools_vr_col .= 100.0

    # --- Patch/column topology: woody + non-woody + crop in one build ---
    # (mortality funcs are not ivt/woody-branchy; kept for parity with the template)
    ivt        = [1, 2, 1, 2, NPCROPMIN, NPCROPMIN]
    mask_soilc = trues(NC)
    mask_soilp = trues(NP)

    return (; cs_veg, cf_veg, cs_soil, ivt, mask_soilc, mask_soilp)
end

run_csu2!(S, m_c, m_p, dt) =
    CLM.c_state_update2!(S.cs_veg, S.cf_veg, S.cs_soil;
        mask_soilc=m_c, mask_soilp=m_p, bounds_col=1:NC, bounds_patch=1:NP,
        nlevdecomp=NLEVDECOMP, i_litr_min=I_LITR_MIN, i_litr_max=I_LITR_MAX, i_cwd=I_CWD, dt=dt)

run_csu2h!(S, m_c, m_p, dt) =
    CLM.c_state_update2h!(S.cs_veg, S.cf_veg, S.cs_soil;
        mask_soilc=m_c, mask_soilp=m_p, bounds_col=1:NC, bounds_patch=1:NP,
        nlevdecomp=NLEVDECOMP, i_litr_min=I_LITR_MIN, i_litr_max=I_LITR_MAX, i_cwd=I_CWD, dt=dt)

run_csu2g!(S, m_c, m_p, dt) =
    CLM.c_state_update2g!(S.cs_veg, S.cf_veg, S.cs_soil;
        mask_soilc=m_c, mask_soilp=m_p, bounds_col=1:NC, bounds_patch=1:NP,
        nlevdecomp=NLEVDECOMP, i_litr_min=I_LITR_MIN, i_litr_max=I_LITR_MAX, i_cwd=I_CWD, dt=dt)

function main(backend)
    println("=" ^ 70)
    println("END-TO-END GPU parity for c_state_update2!/2h!/2g! (BGC mortality)")
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

    if !(cs_d.leafc_patch isa device_array_type())
        println("  BLOCKED: a state struct did not move to the device under adapt.")
        return 2
    end

    dt = FT(1800.0)
    dmask(m) = to(collect(Bool, m))
    mc = dmask(B.mask_soilc); mp = dmask(B.mask_soilp)

    # Run all three mortality state-updates in sequence (CPU + device).
    run_csu2!(H, H.mask_soilc, H.mask_soilp, 1800.0)
    run_csu2!(Sd, mc, mp, dt)
    run_csu2h!(H, H.mask_soilc, H.mask_soilp, 1800.0)
    run_csu2h!(Sd, mc, mp, dt)
    run_csu2g!(H, H.mask_soilc, H.mask_soilp, 1800.0)
    run_csu2g!(Sd, mc, mp, dt)

    checks = [
        ("xsmrpool",          H.cs_veg.xsmrpool_patch,          Sd.cs_veg.xsmrpool_patch),
        ("leafc",             H.cs_veg.leafc_patch,             Sd.cs_veg.leafc_patch),
        ("leafc_storage",     H.cs_veg.leafc_storage_patch,     Sd.cs_veg.leafc_storage_patch),
        ("leafc_xfer",        H.cs_veg.leafc_xfer_patch,        Sd.cs_veg.leafc_xfer_patch),
        ("frootc",            H.cs_veg.frootc_patch,            Sd.cs_veg.frootc_patch),
        ("frootc_storage",    H.cs_veg.frootc_storage_patch,    Sd.cs_veg.frootc_storage_patch),
        ("frootc_xfer",       H.cs_veg.frootc_xfer_patch,       Sd.cs_veg.frootc_xfer_patch),
        ("livestemc",         H.cs_veg.livestemc_patch,         Sd.cs_veg.livestemc_patch),
        ("livestemc_storage", H.cs_veg.livestemc_storage_patch, Sd.cs_veg.livestemc_storage_patch),
        ("livestemc_xfer",    H.cs_veg.livestemc_xfer_patch,    Sd.cs_veg.livestemc_xfer_patch),
        ("deadstemc",         H.cs_veg.deadstemc_patch,         Sd.cs_veg.deadstemc_patch),
        ("deadstemc_storage", H.cs_veg.deadstemc_storage_patch, Sd.cs_veg.deadstemc_storage_patch),
        ("deadstemc_xfer",    H.cs_veg.deadstemc_xfer_patch,    Sd.cs_veg.deadstemc_xfer_patch),
        ("livecrootc",        H.cs_veg.livecrootc_patch,        Sd.cs_veg.livecrootc_patch),
        ("livecrootc_storage",H.cs_veg.livecrootc_storage_patch,Sd.cs_veg.livecrootc_storage_patch),
        ("livecrootc_xfer",   H.cs_veg.livecrootc_xfer_patch,   Sd.cs_veg.livecrootc_xfer_patch),
        ("deadcrootc",        H.cs_veg.deadcrootc_patch,        Sd.cs_veg.deadcrootc_patch),
        ("deadcrootc_storage",H.cs_veg.deadcrootc_storage_patch,Sd.cs_veg.deadcrootc_storage_patch),
        ("deadcrootc_xfer",   H.cs_veg.deadcrootc_xfer_patch,   Sd.cs_veg.deadcrootc_xfer_patch),
        ("gresp_storage",     H.cs_veg.gresp_storage_patch,     Sd.cs_veg.gresp_storage_patch),
        ("gresp_xfer",        H.cs_veg.gresp_xfer_patch,        Sd.cs_veg.gresp_xfer_patch),
        ("decomp_cpools_vr",  H.cs_soil.decomp_cpools_vr_col,   Sd.cs_soil.decomp_cpools_vr_col),
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
    println(nfail == 0 ? "  WHOLE c_state_update2!/2h!/2g! MATCHES CPU ON $name ($FT) ✓" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
