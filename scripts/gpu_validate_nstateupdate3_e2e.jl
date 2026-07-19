# ==========================================================================
# gpu_validate_nstateupdate3_e2e.jl — end-to-end GPU parity for the WHOLE
# n_state_update3! (fire N) and n_state_update_leaching! (sminn leaching) BGC
# drivers (the Phase B C/N state-update cascade — nitrogen fire/leaching arm).
#
# Builds a small multi-column / multi-patch instance mirroring
# test/test_n_state_update3.jl, runs both routines on the CPU, converts every
# state struct to Float32 + adapts to the GPU, runs the SAME calls on the device,
# and compares the mutated outputs field-by-field. The scenario exercises the
# column-kernel CWD+litter inputs, the per-pool fire losses, the patch-kernel
# vegetation fire updates (displayed/storage/transfer/retransn pools, plus the
# live-to-dead fire transfers), and both leaching branches (the
# use_nitrif_denitrif=false sminn path and the =true no3+nh4 max-clamped path).
#
#   julia --project=scripts scripts/gpu_validate_nstateupdate3_e2e.jl
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

# Float32-down-converting Metal adaptor. The CN state structs are concretely
# Float64-typed (default ctor pins ::Vector{Float64}/::Matrix{Float64}) and the
# *_init! routines fill Float64, so we can't setfield! a Float32 array. Instead we
# adapt with a custom storage rule that down-converts float arrays to Float32 as it
# reconstructs the struct (integer/bool arrays move as-is). @adapt_structure rebuilds
# each struct positionally, inferring the new {Float32,…} params from the adapted fields.
struct MetalF32 end
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{<:AbstractFloat}) = device_array_type()(Float32.(x))
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{<:Integer})       = device_array_type()(x)
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{Bool})            = device_array_type()(x)

const NC = 4
const NP = 6
const NG = 2
const NLEVDECOMP = 2
const NDECOMP_POOLS = 7
const NREPR = 1
const I_LITR_MIN = 1
const I_LITR_MAX = 3
const I_CWD = 4

# Build the CPU reference state (Float64) for the fire (n_state_update3!) scenario.
function build_fire()
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
    CLM.cnveg_nitrogen_flux_init!(nf_veg, NP, NC, NG; nrepr=NREPR,
                                   nlevdecomp_full=NLEVDECOMP,
                                   ndecomp_pools=NDECOMP_POOLS,
                                   i_litr_max=I_LITR_MAX)

    # Fire mortality fluxes — patch level to combustion
    nf_veg.m_leafn_to_fire_patch              .= 1.0e-6
    nf_veg.m_frootn_to_fire_patch             .= 0.8e-6
    nf_veg.m_livestemn_to_fire_patch          .= 0.5e-6
    nf_veg.m_deadstemn_to_fire_patch          .= 0.3e-6
    nf_veg.m_livecrootn_to_fire_patch         .= 0.4e-6
    nf_veg.m_deadcrootn_to_fire_patch         .= 0.2e-6
    nf_veg.m_leafn_storage_to_fire_patch      .= 0.1e-6
    nf_veg.m_frootn_storage_to_fire_patch     .= 0.08e-6
    nf_veg.m_livestemn_storage_to_fire_patch  .= 0.05e-6
    nf_veg.m_deadstemn_storage_to_fire_patch  .= 0.03e-6
    nf_veg.m_livecrootn_storage_to_fire_patch .= 0.04e-6
    nf_veg.m_deadcrootn_storage_to_fire_patch .= 0.02e-6
    nf_veg.m_leafn_xfer_to_fire_patch         .= 0.05e-6
    nf_veg.m_frootn_xfer_to_fire_patch        .= 0.04e-6
    nf_veg.m_livestemn_xfer_to_fire_patch     .= 0.03e-6
    nf_veg.m_deadstemn_xfer_to_fire_patch     .= 0.02e-6
    nf_veg.m_livecrootn_xfer_to_fire_patch    .= 0.015e-6
    nf_veg.m_deadcrootn_xfer_to_fire_patch    .= 0.01e-6
    nf_veg.m_retransn_to_fire_patch           .= 0.12e-6

    # Fire mortality fluxes — patch level to litter
    nf_veg.m_leafn_to_litter_fire_patch              .= 0.5e-6
    nf_veg.m_frootn_to_litter_fire_patch             .= 0.4e-6
    nf_veg.m_livestemn_to_litter_fire_patch          .= 0.3e-6
    nf_veg.m_deadstemn_to_litter_fire_patch          .= 0.15e-6
    nf_veg.m_livecrootn_to_litter_fire_patch         .= 0.2e-6
    nf_veg.m_deadcrootn_to_litter_fire_patch         .= 0.1e-6
    nf_veg.m_leafn_storage_to_litter_fire_patch      .= 0.05e-6
    nf_veg.m_frootn_storage_to_litter_fire_patch     .= 0.04e-6
    nf_veg.m_livestemn_storage_to_litter_fire_patch  .= 0.03e-6
    nf_veg.m_deadstemn_storage_to_litter_fire_patch  .= 0.015e-6
    nf_veg.m_livecrootn_storage_to_litter_fire_patch .= 0.02e-6
    nf_veg.m_deadcrootn_storage_to_litter_fire_patch .= 0.01e-6
    nf_veg.m_leafn_xfer_to_litter_fire_patch         .= 0.025e-6
    nf_veg.m_frootn_xfer_to_litter_fire_patch        .= 0.02e-6
    nf_veg.m_livestemn_xfer_to_litter_fire_patch     .= 0.015e-6
    nf_veg.m_deadstemn_xfer_to_litter_fire_patch     .= 0.01e-6
    nf_veg.m_livecrootn_xfer_to_litter_fire_patch    .= 0.008e-6
    nf_veg.m_deadcrootn_xfer_to_litter_fire_patch    .= 0.005e-6
    nf_veg.m_retransn_to_litter_fire_patch           .= 0.06e-6

    # Live-to-dead fire transfers
    nf_veg.m_livestemn_to_deadstemn_fire_patch       .= 0.2e-6
    nf_veg.m_livecrootn_to_deadcrootn_fire_patch     .= 0.15e-6

    # Column-level fire fluxes
    nf_veg.fire_mortality_n_to_cwdn_col   .= 1.0e-6
    nf_veg.m_n_to_litr_fire_col           .= 2.0e-6
    nf_veg.m_decomp_npools_to_fire_vr_col .= 0.5e-6

    ns_soil = CLM.SoilBiogeochemNitrogenStateData()
    CLM.soil_bgc_nitrogen_state_init!(ns_soil, NC, NG, NLEVDECOMP, NDECOMP_POOLS)
    ns_soil.decomp_npools_vr_col .= 10.0

    mask_soilc = trues(NC)
    mask_soilp = trues(NP)

    return (; ns_veg, nf_veg, ns_soil, mask_soilc, mask_soilp)
end

# Build the CPU reference state (Float64) for the leaching scenario, exercising
# both the use_nitrif_denitrif=false (sminn) and =true (no3+nh4, max-clamped) paths.
function build_leach(use_nitrif_denitrif::Bool)
    ns_soil = CLM.SoilBiogeochemNitrogenStateData()
    CLM.soil_bgc_nitrogen_state_init!(ns_soil, NC, NG, NLEVDECOMP, NDECOMP_POOLS)
    nf_soil = CLM.SoilBiogeochemNitrogenFluxData()
    CLM.soil_bgc_nitrogen_flux_init!(nf_soil, NC, NLEVDECOMP, NDECOMP_POOLS, 1)

    if !use_nitrif_denitrif
        ns_soil.sminn_vr_col .= 10.0
        nf_soil.sminn_leached_vr_col .= 1.0e-6
    else
        ns_soil.smin_no3_vr_col .= 5.0
        ns_soil.smin_nh4_vr_col .= 3.0
        ns_soil.sminn_vr_col    .= 8.0
        nf_soil.smin_no3_leached_vr_col .= 0.5e-6
        nf_soil.smin_no3_runoff_vr_col  .= 0.3e-6
    end

    mask_soilc = trues(NC)
    return (; ns_soil, nf_soil, mask_soilc)
end

run_nsu3!(S, m_c, m_p, dt) =
    CLM.n_state_update3!(S.ns_veg, S.nf_veg, S.ns_soil;
        mask_soilc=m_c, mask_soilp=m_p, bounds_col=1:NC, bounds_patch=1:NP,
        nlevdecomp=NLEVDECOMP, ndecomp_pools=NDECOMP_POOLS,
        i_litr_min=I_LITR_MIN, i_litr_max=I_LITR_MAX, i_cwd=I_CWD, dt=dt)

run_leach!(S, m_c, use_nd, dt) =
    CLM.n_state_update_leaching!(S.ns_soil, S.nf_soil;
        mask_soilc=m_c, bounds_col=1:NC, nlevdecomp=NLEVDECOMP,
        use_nitrif_denitrif=use_nd, dt=dt)

function main(backend)
    println("=" ^ 70)
    println("END-TO-END GPU parity for n_state_update3! + n_state_update_leaching!")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU path exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    ad(x) = CLM.Adapt.adapt(MetalF32(), x)
    dmask(m) = to(collect(Bool, m))

    # ---- Fire (n_state_update3!) ----
    Hf = build_fire(); Bf = build_fire()
    nsv_d = ad(Bf.ns_veg); nfv_d = ad(Bf.nf_veg); nss_d = ad(Bf.ns_soil)
    Sf = (; ns_veg=nsv_d, nf_veg=nfv_d, ns_soil=nss_d)

    if !(nsv_d.leafn_patch isa device_array_type())
        println("  BLOCKED: a state struct did not move to the device under adapt.")
        return 2
    end

    dt = FT(1800.0)
    run_nsu3!(Hf, Hf.mask_soilc, Hf.mask_soilp, 1800.0)
    run_nsu3!(Sf, dmask(Bf.mask_soilc), dmask(Bf.mask_soilp), dt)

    # ---- Leaching (n_state_update_leaching!), both branches ----
    Hl0 = build_leach(false); Bl0 = build_leach(false)
    Sl0 = (; ns_soil=ad(Bl0.ns_soil), nf_soil=ad(Bl0.nf_soil))
    run_leach!(Hl0, Hl0.mask_soilc, false, 1800.0)
    run_leach!(Sl0, dmask(Bl0.mask_soilc), false, dt)

    Hl1 = build_leach(true); Bl1 = build_leach(true)
    Sl1 = (; ns_soil=ad(Bl1.ns_soil), nf_soil=ad(Bl1.nf_soil))
    run_leach!(Hl1, Hl1.mask_soilc, true, 1800.0)
    run_leach!(Sl1, dmask(Bl1.mask_soilc), true, dt)

    checks = [
        # fire: vegetation N pools
        ("leafn",            Hf.ns_veg.leafn_patch,            Sf.ns_veg.leafn_patch),
        ("frootn",           Hf.ns_veg.frootn_patch,           Sf.ns_veg.frootn_patch),
        ("livestemn",        Hf.ns_veg.livestemn_patch,        Sf.ns_veg.livestemn_patch),
        ("deadstemn",        Hf.ns_veg.deadstemn_patch,        Sf.ns_veg.deadstemn_patch),
        ("livecrootn",       Hf.ns_veg.livecrootn_patch,       Sf.ns_veg.livecrootn_patch),
        ("deadcrootn",       Hf.ns_veg.deadcrootn_patch,       Sf.ns_veg.deadcrootn_patch),
        ("leafn_storage",    Hf.ns_veg.leafn_storage_patch,    Sf.ns_veg.leafn_storage_patch),
        ("frootn_storage",   Hf.ns_veg.frootn_storage_patch,   Sf.ns_veg.frootn_storage_patch),
        ("livestemn_storage",Hf.ns_veg.livestemn_storage_patch,Sf.ns_veg.livestemn_storage_patch),
        ("deadstemn_storage",Hf.ns_veg.deadstemn_storage_patch,Sf.ns_veg.deadstemn_storage_patch),
        ("livecrootn_storage",Hf.ns_veg.livecrootn_storage_patch,Sf.ns_veg.livecrootn_storage_patch),
        ("deadcrootn_storage",Hf.ns_veg.deadcrootn_storage_patch,Sf.ns_veg.deadcrootn_storage_patch),
        ("leafn_xfer",       Hf.ns_veg.leafn_xfer_patch,       Sf.ns_veg.leafn_xfer_patch),
        ("frootn_xfer",      Hf.ns_veg.frootn_xfer_patch,      Sf.ns_veg.frootn_xfer_patch),
        ("livestemn_xfer",   Hf.ns_veg.livestemn_xfer_patch,   Sf.ns_veg.livestemn_xfer_patch),
        ("deadstemn_xfer",   Hf.ns_veg.deadstemn_xfer_patch,   Sf.ns_veg.deadstemn_xfer_patch),
        ("livecrootn_xfer",  Hf.ns_veg.livecrootn_xfer_patch,  Sf.ns_veg.livecrootn_xfer_patch),
        ("deadcrootn_xfer",  Hf.ns_veg.deadcrootn_xfer_patch,  Sf.ns_veg.deadcrootn_xfer_patch),
        ("retransn",         Hf.ns_veg.retransn_patch,         Sf.ns_veg.retransn_patch),
        # fire: soil decomposition N pools (column kernels)
        ("decomp_npools_vr", Hf.ns_soil.decomp_npools_vr_col,  Sf.ns_soil.decomp_npools_vr_col),
        # leaching (no nitrif_denitrif)
        ("leach_sminn",      Hl0.ns_soil.sminn_vr_col,         Sl0.ns_soil.sminn_vr_col),
        # leaching (nitrif_denitrif)
        ("leach_no3",        Hl1.ns_soil.smin_no3_vr_col,      Sl1.ns_soil.smin_no3_vr_col),
        ("leach_nd_sminn",   Hl1.ns_soil.sminn_vr_col,         Sl1.ns_soil.sminn_vr_col),
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
    println(nfail == 0 ?
        "  WHOLE n_state_update3! + n_state_update_leaching! MATCH CPU ON $name ($FT) ✓" :
        "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
