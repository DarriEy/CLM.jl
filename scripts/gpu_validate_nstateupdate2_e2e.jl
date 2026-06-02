# ==========================================================================
# gpu_validate_nstateupdate2_e2e.jl — end-to-end Metal parity for the WHOLE
# n_state_update2! / n_state_update2h! / n_state_update2g! BGC N-state mortality
# updates (gap-phase / harvest / gross-unrepresented-LCC nitrogen).
#
# Builds a small multi-column / multi-patch instance mirroring
# test/test_n_state_update2.jl, runs each function on the CPU, converts every
# state struct to Float32 + adapts to Metal, runs the SAME calls on the device, and
# compares the mutated outputs field-by-field. Each function exercises its own
# column kernel (j/i decomp RMW into decomp_npools_vr_col) and patch kernel
# (per-patch displayed/storage/transfer N pools).
#
#   julia --project=scripts scripts/gpu_validate_nstateupdate2_e2e.jl
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
# reconstructs the struct (integer/bool arrays move as-is).
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

# Build the CPU reference state (Float64). Returns a NamedTuple of the structs +
# the masks. Mirrors make_nstate_update2_data() in test/test_n_state_update2.jl.
function build()
    ns_veg = CLM.CNVegNitrogenStateData()
    CLM.cnveg_nitrogen_state_init!(ns_veg, NP, NC, NG; nrepr=NREPR)
    ns_veg.leafn_patch             .= 5.0
    ns_veg.leafn_storage_patch     .= 1.0
    ns_veg.leafn_xfer_patch        .= 0.5
    ns_veg.frootn_patch            .= 4.0
    ns_veg.frootn_storage_patch    .= 0.8
    ns_veg.frootn_xfer_patch       .= 0.4
    ns_veg.livestemn_patch         .= 3.0
    ns_veg.livestemn_storage_patch .= 0.6
    ns_veg.livestemn_xfer_patch    .= 0.3
    ns_veg.deadstemn_patch         .= 20.0
    ns_veg.deadstemn_storage_patch .= 0.4
    ns_veg.deadstemn_xfer_patch    .= 0.2
    ns_veg.livecrootn_patch        .= 2.0
    ns_veg.livecrootn_storage_patch .= 0.4
    ns_veg.livecrootn_xfer_patch   .= 0.2
    ns_veg.deadcrootn_patch        .= 10.0
    ns_veg.deadcrootn_storage_patch .= 0.2
    ns_veg.deadcrootn_xfer_patch   .= 0.1
    ns_veg.retransn_patch          .= 2.0

    nf_veg = CLM.CNVegNitrogenFluxData()
    CLM.cnveg_nitrogen_flux_init!(nf_veg, NP, NC, NG;
                                   nrepr=NREPR,
                                   nlevdecomp_full=NLEVDECOMP,
                                   ndecomp_pools=NDECOMP_POOLS,
                                   i_litr_max=I_LITR_MAX)

    # Gap mortality fluxes (NStateUpdate2)
    nf_veg.m_leafn_to_litter_patch              .= 1.0e-6
    nf_veg.m_frootn_to_litter_patch             .= 0.8e-6
    nf_veg.m_livestemn_to_litter_patch          .= 0.5e-6
    nf_veg.m_deadstemn_to_litter_patch          .= 0.3e-6
    nf_veg.m_livecrootn_to_litter_patch         .= 0.4e-6
    nf_veg.m_deadcrootn_to_litter_patch         .= 0.2e-6
    nf_veg.m_retransn_to_litter_patch           .= 0.15e-6
    nf_veg.m_leafn_storage_to_litter_patch      .= 0.1e-6
    nf_veg.m_frootn_storage_to_litter_patch     .= 0.08e-6
    nf_veg.m_livestemn_storage_to_litter_patch  .= 0.05e-6
    nf_veg.m_deadstemn_storage_to_litter_patch  .= 0.03e-6
    nf_veg.m_livecrootn_storage_to_litter_patch .= 0.04e-6
    nf_veg.m_deadcrootn_storage_to_litter_patch .= 0.02e-6
    nf_veg.m_leafn_xfer_to_litter_patch         .= 0.05e-6
    nf_veg.m_frootn_xfer_to_litter_patch        .= 0.04e-6
    nf_veg.m_livestemn_xfer_to_litter_patch     .= 0.03e-6
    nf_veg.m_deadstemn_xfer_to_litter_patch     .= 0.02e-6
    nf_veg.m_livecrootn_xfer_to_litter_patch    .= 0.015e-6
    nf_veg.m_deadcrootn_xfer_to_litter_patch    .= 0.01e-6
    nf_veg.gap_mortality_n_to_litr_n_col        .= 2.0e-6
    nf_veg.gap_mortality_n_to_cwdn_col          .= 1.0e-6

    # Harvest mortality fluxes (NStateUpdate2h)
    nf_veg.hrv_leafn_to_litter_patch             .= 0.5e-6
    nf_veg.hrv_frootn_to_litter_patch            .= 0.4e-6
    nf_veg.hrv_livestemn_to_litter_patch         .= 0.3e-6
    nf_veg.hrv_livecrootn_to_litter_patch        .= 0.2e-6
    nf_veg.hrv_deadcrootn_to_litter_patch        .= 0.1e-6
    nf_veg.wood_harvestn_patch                   .= 0.6e-6
    nf_veg.hrv_retransn_to_litter_patch          .= 0.12e-6
    nf_veg.hrv_leafn_storage_to_litter_patch     .= 0.05e-6
    nf_veg.hrv_frootn_storage_to_litter_patch    .= 0.04e-6
    nf_veg.hrv_livestemn_storage_to_litter_patch .= 0.03e-6
    nf_veg.hrv_deadstemn_storage_to_litter_patch .= 0.02e-6
    nf_veg.hrv_livecrootn_storage_to_litter_patch .= 0.015e-6
    nf_veg.hrv_deadcrootn_storage_to_litter_patch .= 0.01e-6
    nf_veg.hrv_leafn_xfer_to_litter_patch        .= 0.025e-6
    nf_veg.hrv_frootn_xfer_to_litter_patch       .= 0.02e-6
    nf_veg.hrv_livestemn_xfer_to_litter_patch    .= 0.015e-6
    nf_veg.hrv_deadstemn_xfer_to_litter_patch    .= 0.01e-6
    nf_veg.hrv_livecrootn_xfer_to_litter_patch   .= 0.008e-6
    nf_veg.hrv_deadcrootn_xfer_to_litter_patch   .= 0.005e-6
    nf_veg.harvest_n_to_litr_n_col               .= 1.5e-6
    nf_veg.harvest_n_to_cwdn_col                 .= 0.8e-6

    # Gross unrepresented landcover change fluxes (NStateUpdate2g)
    nf_veg.gru_leafn_to_litter_patch            .= 0.3e-6
    nf_veg.gru_frootn_to_litter_patch           .= 0.25e-6
    nf_veg.gru_livestemn_to_atm_patch           .= 0.2e-6
    nf_veg.gru_deadstemn_to_atm_patch           .= 0.15e-6
    nf_veg.gru_wood_productn_gain_patch         .= 0.1e-6
    nf_veg.gru_livecrootn_to_litter_patch       .= 0.12e-6
    nf_veg.gru_deadcrootn_to_litter_patch       .= 0.08e-6
    nf_veg.gru_retransn_to_litter_patch         .= 0.06e-6
    nf_veg.gru_leafn_storage_to_atm_patch       .= 0.03e-6
    nf_veg.gru_frootn_storage_to_atm_patch      .= 0.025e-6
    nf_veg.gru_livestemn_storage_to_atm_patch   .= 0.02e-6
    nf_veg.gru_deadstemn_storage_to_atm_patch   .= 0.015e-6
    nf_veg.gru_livecrootn_storage_to_atm_patch  .= 0.012e-6
    nf_veg.gru_deadcrootn_storage_to_atm_patch  .= 0.008e-6
    nf_veg.gru_leafn_xfer_to_atm_patch          .= 0.015e-6
    nf_veg.gru_frootn_xfer_to_atm_patch         .= 0.012e-6
    nf_veg.gru_livestemn_xfer_to_atm_patch      .= 0.01e-6
    nf_veg.gru_deadstemn_xfer_to_atm_patch      .= 0.008e-6
    nf_veg.gru_livecrootn_xfer_to_atm_patch     .= 0.006e-6
    nf_veg.gru_deadcrootn_xfer_to_atm_patch     .= 0.004e-6
    nf_veg.gru_n_to_litr_n_col                  .= 1.0e-6
    nf_veg.gru_n_to_cwdn_col                    .= 0.5e-6

    ns_soil = CLM.SoilBiogeochemNitrogenStateData()
    CLM.soil_bgc_nitrogen_state_init!(ns_soil, NC, NG, NLEVDECOMP, NDECOMP_POOLS)
    ns_soil.decomp_npools_vr_col .= 10.0

    mask_soilc = trues(NC)
    mask_soilp = trues(NP)

    return (; ns_veg, nf_veg, ns_soil, mask_soilc, mask_soilp)
end

# --- Per-function runners (CPU + device share these; only the args differ) ---
run_nsu2!(S, m_c, m_p, dt) = CLM.n_state_update2!(S.ns_veg, S.nf_veg, S.ns_soil;
    mask_soilc=m_c, mask_soilp=m_p, bounds_col=1:NC, bounds_patch=1:NP,
    nlevdecomp=NLEVDECOMP, i_litr_min=I_LITR_MIN, i_litr_max=I_LITR_MAX, i_cwd=I_CWD, dt=dt)

run_nsu2h!(S, m_c, m_p, dt) = CLM.n_state_update2h!(S.ns_veg, S.nf_veg, S.ns_soil;
    mask_soilc=m_c, mask_soilp=m_p, bounds_col=1:NC, bounds_patch=1:NP,
    nlevdecomp=NLEVDECOMP, i_litr_min=I_LITR_MIN, i_litr_max=I_LITR_MAX, i_cwd=I_CWD, dt=dt)

run_nsu2g!(S, m_c, m_p, dt) = CLM.n_state_update2g!(S.ns_veg, S.nf_veg, S.ns_soil;
    mask_soilc=m_c, mask_soilp=m_p, bounds_col=1:NC, bounds_patch=1:NP,
    nlevdecomp=NLEVDECOMP, i_litr_min=I_LITR_MIN, i_litr_max=I_LITR_MAX, i_cwd=I_CWD, dt=dt)

# Common set of mutated outputs to compare for any of the three functions.
checks_for(H, Sd) = [
    ("leafn",            H.ns_veg.leafn_patch,            Sd.ns_veg.leafn_patch),
    ("frootn",           H.ns_veg.frootn_patch,           Sd.ns_veg.frootn_patch),
    ("livestemn",        H.ns_veg.livestemn_patch,        Sd.ns_veg.livestemn_patch),
    ("deadstemn",        H.ns_veg.deadstemn_patch,        Sd.ns_veg.deadstemn_patch),
    ("livecrootn",       H.ns_veg.livecrootn_patch,       Sd.ns_veg.livecrootn_patch),
    ("deadcrootn",       H.ns_veg.deadcrootn_patch,       Sd.ns_veg.deadcrootn_patch),
    ("retransn",         H.ns_veg.retransn_patch,         Sd.ns_veg.retransn_patch),
    ("leafn_storage",    H.ns_veg.leafn_storage_patch,    Sd.ns_veg.leafn_storage_patch),
    ("frootn_storage",   H.ns_veg.frootn_storage_patch,   Sd.ns_veg.frootn_storage_patch),
    ("livestemn_storage",H.ns_veg.livestemn_storage_patch,Sd.ns_veg.livestemn_storage_patch),
    ("deadstemn_storage",H.ns_veg.deadstemn_storage_patch,Sd.ns_veg.deadstemn_storage_patch),
    ("livecrootn_storage",H.ns_veg.livecrootn_storage_patch,Sd.ns_veg.livecrootn_storage_patch),
    ("deadcrootn_storage",H.ns_veg.deadcrootn_storage_patch,Sd.ns_veg.deadcrootn_storage_patch),
    ("leafn_xfer",       H.ns_veg.leafn_xfer_patch,       Sd.ns_veg.leafn_xfer_patch),
    ("frootn_xfer",      H.ns_veg.frootn_xfer_patch,      Sd.ns_veg.frootn_xfer_patch),
    ("livestemn_xfer",   H.ns_veg.livestemn_xfer_patch,   Sd.ns_veg.livestemn_xfer_patch),
    ("deadstemn_xfer",   H.ns_veg.deadstemn_xfer_patch,   Sd.ns_veg.deadstemn_xfer_patch),
    ("livecrootn_xfer",  H.ns_veg.livecrootn_xfer_patch,  Sd.ns_veg.livecrootn_xfer_patch),
    ("deadcrootn_xfer",  H.ns_veg.deadcrootn_xfer_patch,  Sd.ns_veg.deadcrootn_xfer_patch),
    ("decomp_npools_vr", H.ns_soil.decomp_npools_vr_col,  Sd.ns_soil.decomp_npools_vr_col),
]

# Run one function on CPU + device and report field-by-field parity.
function check_fn!(label, runfn!, to, FT)
    H  = build()    # CPU reference (Float64)
    B  = build()    # source for the device snapshot

    ad(x) = CLM.Adapt.adapt(MetalF32(), x)
    Sd = (; ns_veg = ad(B.ns_veg), nf_veg = ad(B.nf_veg), ns_soil = ad(B.ns_soil))

    if !(Sd.ns_veg.leafn_patch isa Metal.MtlArray)
        println("  BLOCKED: a state struct did not move to the device under adapt.")
        return 2
    end

    dt = FT(1800.0)
    dmask(m) = to(collect(Bool, m))

    runfn!(H, H.mask_soilc, H.mask_soilp, 1800.0)
    runfn!(Sd, dmask(B.mask_soilc), dmask(B.mask_soilp), dt)

    @printf("  --- %s ---\n", label)
    nfail = 0
    for (nm, a, b) in checks_for(H, Sd)
        if !cpu_has_finite(a)
            @printf("    [WARN] %-20s CPU reference is all-NaN/Inf — skipping\n", nm)
            continue
        end
        d = reldiff(a, b); ok = d < 1f-3
        @printf("    [%s] %-20s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end
    return nfail
end

function main(backend)
    println("=" ^ 70)
    println("END-TO-END Metal parity for n_state_update2!/2h!/2g! (BGC N-state mortality)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU path exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    nfail = 0
    nfail += check_fn!("n_state_update2!  (gap-phase mortality)",      run_nsu2!,  to, FT)
    nfail += check_fn!("n_state_update2h! (harvest mortality)",        run_nsu2h!, to, FT)
    nfail += check_fn!("n_state_update2g! (gross unrep LCC mortality)", run_nsu2g!, to, FT)

    println()
    println(nfail == 0 ? "  WHOLE n_state_update2 family MATCHES CPU ON $name ($FT) ✓" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
