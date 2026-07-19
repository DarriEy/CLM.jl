# ==========================================================================
# gpu_validate_nstateupdate1_e2e.jl — end-to-end GPU parity for the WHOLE
# n_state_update1! + n_state_update_dyn_patch! BGC drivers (the nitrogen analog
# of the Phase B C-state cascade entry point).
#
# Builds a small multi-column / multi-patch instance mirroring
# test/test_n_state_update1.jl, runs n_state_update_dyn_patch! then
# n_state_update1! on the CPU, converts every state struct to Float32 + adapts to
# Metal, runs the SAME calls on the device, and compares the mutated outputs
# field-by-field. The scenario deliberately exercises the branchy paths in one
# shot: a non-woody patch, a woody patch, two crop patches, a FATES column
# (column-kernel skip), use_fun=true (the free_retransn flux RMW), and the
# dyn-patch column/gridcell updates.
#
#   julia --project=scripts scripts/gpu_validate_nstateupdate1_e2e.jl
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
const NLEVDECOMP = 1
const NDECOMP_POOLS = 7
const NCASCADE = 5
const NREPR = 1
const I_LITR_MIN = 1
const I_LITR_MAX = 3
const I_CWD = 4
const NPCROPMIN = 17

# Build the CPU reference state (Float64). Returns a NamedTuple of the structs +
# the scalar/index params the N-update routines need.
function build()
    ns_veg = CLM.CNVegNitrogenStateData()
    CLM.cnveg_nitrogen_state_init!(ns_veg, NP, NC, NG; nrepr=NREPR)
    nf_veg = CLM.CNVegNitrogenFluxData()
    CLM.cnveg_nitrogen_flux_init!(nf_veg, NP, NC, NG; nrepr=NREPR,
                                   nlevdecomp_full=NLEVDECOMP, ndecomp_pools=NDECOMP_POOLS,
                                   i_litr_max=I_LITR_MAX)
    nf_soil = CLM.SoilBiogeochemNitrogenFluxData()
    CLM.soil_bgc_nitrogen_flux_init!(nf_soil, NC, NLEVDECOMP, NDECOMP_POOLS, NCASCADE)
    ns_soil = CLM.SoilBiogeochemNitrogenStateData()
    CLM.soil_bgc_nitrogen_state_init!(ns_soil, NC, NG, NLEVDECOMP, NDECOMP_POOLS)

    # --- Nitrogen-state pools (known nonzero so every flux moves something) ---
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
    ns_veg.npool_patch             .= 10.0
    ns_veg.seedn_grc               .= 100.0
    ns_veg.cropseedn_deficit_patch .= -0.5
    for k in 1:NREPR
        ns_veg.reproductiven_patch[:, k]         .= 1.0
        ns_veg.reproductiven_storage_patch[:, k] .= 0.2
        ns_veg.reproductiven_xfer_patch[:, k]    .= 0.1
    end
    ns_soil.decomp_npools_vr_col .= 10.0

    # --- Nitrogen fluxes (mirror test_n_state_update1.jl + extras so all branches
    #     move; reproductive/crop fluxes nonzero to exercise the k loops). ---
    nf_veg.leafn_xfer_to_leafn_patch        .= 2.0e-7
    nf_veg.frootn_xfer_to_frootn_patch      .= 1.0e-7
    nf_veg.livestemn_xfer_to_livestemn_patch .= 1.5e-7
    nf_veg.deadstemn_xfer_to_deadstemn_patch .= 1.0e-7
    nf_veg.livecrootn_xfer_to_livecrootn_patch .= 0.8e-7
    nf_veg.deadcrootn_xfer_to_deadcrootn_patch .= 0.5e-7
    nf_veg.leafn_to_litter_patch            .= 3.0e-7
    nf_veg.frootn_to_litter_patch           .= 2.0e-7
    nf_veg.leafn_to_retransn_patch          .= 1.0e-7
    nf_veg.livestemn_to_deadstemn_patch     .= 1.0e-8
    nf_veg.livestemn_to_retransn_patch      .= 0.5e-8
    nf_veg.livecrootn_to_deadcrootn_patch   .= 1.0e-8
    nf_veg.livecrootn_to_retransn_patch     .= 0.5e-8
    nf_veg.free_retransn_to_npool_patch     .= 0.0
    nf_veg.frootn_to_retransn_patch         .= 1.0e-8
    nf_veg.livestemn_to_litter_patch        .= 1.0e-9
    nf_veg.livestemn_to_biofueln_patch      .= 1.0e-9
    nf_veg.livestemn_to_removedresiduen_patch .= 1.0e-9
    nf_veg.leafn_to_biofueln_patch          .= 1.0e-9
    nf_veg.leafn_to_removedresiduen_patch   .= 1.0e-9
    nf_veg.crop_seedn_to_leaf_patch         .= 2.0e-9
    for k in 1:NREPR
        nf_veg.repr_grainn_to_food_patch[:, k]         .= 1.0e-9
        nf_veg.repr_grainn_to_seed_patch[:, k]         .= 1.0e-9
        nf_veg.repr_structuren_to_cropprod_patch[:, k] .= 1.0e-9
        nf_veg.repr_structuren_to_litter_patch[:, k]   .= 1.0e-9
        nf_veg.reproductiven_xfer_to_reproductiven_patch[:, k] .= 1.0e-9
        nf_veg.reproductiven_storage_to_xfer_patch[:, k] .= 1.0e-9
        nf_veg.npool_to_reproductiven_patch[:, k]      .= 1.0e-9
        nf_veg.npool_to_reproductiven_storage_patch[:, k] .= 1.0e-9
    end
    nf_veg.sminn_to_npool_patch             .= 5.0e-7
    nf_veg.retransn_to_npool_patch          .= 2.0e-7
    nf_veg.npool_to_leafn_patch             .= 5.0e-7
    nf_veg.npool_to_leafn_storage_patch     .= 2.0e-7
    nf_veg.npool_to_frootn_patch            .= 3.0e-7
    nf_veg.npool_to_frootn_storage_patch    .= 1.0e-7
    nf_veg.npool_to_livestemn_patch         .= 2.0e-7
    nf_veg.npool_to_livestemn_storage_patch .= 1.0e-7
    nf_veg.npool_to_deadstemn_patch         .= 1.5e-7
    nf_veg.npool_to_deadstemn_storage_patch .= 0.5e-7
    nf_veg.npool_to_livecrootn_patch        .= 1.0e-7
    nf_veg.npool_to_livecrootn_storage_patch .= 0.5e-7
    nf_veg.npool_to_deadcrootn_patch        .= 0.5e-7
    nf_veg.npool_to_deadcrootn_storage_patch .= 0.3e-7
    nf_veg.leafn_storage_to_xfer_patch      .= 1.0e-9
    nf_veg.frootn_storage_to_xfer_patch     .= 1.0e-9
    nf_veg.livestemn_storage_to_xfer_patch  .= 1.0e-9
    nf_veg.deadstemn_storage_to_xfer_patch  .= 1.0e-9
    nf_veg.livecrootn_storage_to_xfer_patch .= 1.0e-9
    nf_veg.deadcrootn_storage_to_xfer_patch .= 1.0e-9

    # --- dyn-patch fluxes (column + gridcell) ---
    for i in I_LITR_MIN:I_LITR_MAX
        nf_veg.dwt_frootn_to_litr_n_col[:, :, i] .= 1.0e-6
    end
    nf_veg.dwt_livecrootn_to_cwdn_col .= 2.0e-6
    nf_veg.dwt_deadcrootn_to_cwdn_col .= 3.0e-6
    nf_veg.dwt_seedn_to_leaf_grc      .= 0.5e-6
    nf_veg.dwt_seedn_to_deadstem_grc  .= 0.3e-6

    # --- column input flux for n_state_update1! ---
    nf_veg.phenology_n_to_litr_n_col  .= 1.0e-7

    # --- Soil BGC nitrogen flux (column-kernel output) ---
    nf_soil.decomp_npools_sourcesink_col .= 0.0

    # --- Patch/column topology + branch selectors ---
    # patches: 1 non-woody, 2 woody, 3 non-woody, 4 woody, 5 crop, 6 crop
    ivt          = [1, 2, 1, 2, NPCROPMIN, NPCROPMIN]
    woody        = zeros(Float64, 80); woody[2] = 1.0
    col_is_fates = [false, false, false, true]    # column 4 is FATES (col-kernel skip)
    mask_soilc   = trues(NC)
    mask_soilp   = trues(NP)
    mask_soilc_with_inactive = trues(NC)

    return (; ns_veg, nf_veg, nf_soil, ns_soil,
              ivt, woody, col_is_fates, mask_soilc, mask_soilp, mask_soilc_with_inactive)
end

# Run dyn_patch then update1 on a state set S, given device-or-host arrays.
function run_nsu1!(S, m_c, m_p, m_cinact, fates, ivt, woody, dt)
    CLM.n_state_update_dyn_patch!(S.ns_veg, S.nf_veg, S.ns_soil;
        mask_soilc_with_inactive=m_cinact, bounds_col=1:NC, bounds_grc=1:NG,
        nlevdecomp=NLEVDECOMP, i_litr_min=I_LITR_MIN, i_litr_max=I_LITR_MAX,
        i_cwd=I_CWD, dt=dt)
    CLM.n_state_update1!(S.ns_veg, S.nf_veg, S.nf_soil;
        mask_soilc=m_c, mask_soilp=m_p, bounds_col=1:NC, bounds_patch=1:NP,
        ivt=ivt, woody=woody, col_is_fates=fates,
        nlevdecomp=NLEVDECOMP, i_litr_min=I_LITR_MIN, i_litr_max=I_LITR_MAX,
        i_cwd=I_CWD, npcropmin=NPCROPMIN, nrepr=NREPR,
        repr_grain_min=1, repr_grain_max=NREPR,
        repr_structure_min=1, repr_structure_max=0,
        use_matrixcn=false, use_soil_matrixcn=false, use_fun=true, dt=dt)
end

function main(backend)
    println("=" ^ 70)
    println("END-TO-END GPU parity for n_state_update1! (BGC N-state cascade)")
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
    ns_d  = ad(B.ns_veg); nf_d = ad(B.nf_veg)
    nfs_d = ad(B.nf_soil); nss_d = ad(B.ns_soil)
    Sd = (; ns_veg=ns_d, nf_veg=nf_d, nf_soil=nfs_d, ns_soil=nss_d)

    if !(ns_d.npool_patch isa device_array_type())
        println("  BLOCKED: a state struct did not move to the device under adapt.")
        return 2
    end

    dt = FT(1800.0)
    dmask(m) = to(collect(Bool, m))

    # CPU
    run_nsu1!(H, H.mask_soilc, H.mask_soilp, H.mask_soilc_with_inactive,
              H.col_is_fates, H.ivt, H.woody, 1800.0)
    # Device
    run_nsu1!(Sd, dmask(B.mask_soilc), dmask(B.mask_soilp), dmask(B.mask_soilc_with_inactive),
              dmask(B.col_is_fates), to(B.ivt), to(Float32.(B.woody)), dt)

    checks = [
        ("npool",            H.ns_veg.npool_patch,            Sd.ns_veg.npool_patch),
        ("retransn",         H.ns_veg.retransn_patch,         Sd.ns_veg.retransn_patch),
        ("leafn",            H.ns_veg.leafn_patch,            Sd.ns_veg.leafn_patch),
        ("leafn_storage",    H.ns_veg.leafn_storage_patch,    Sd.ns_veg.leafn_storage_patch),
        ("leafn_xfer",       H.ns_veg.leafn_xfer_patch,       Sd.ns_veg.leafn_xfer_patch),
        ("frootn",           H.ns_veg.frootn_patch,           Sd.ns_veg.frootn_patch),
        ("frootn_storage",   H.ns_veg.frootn_storage_patch,   Sd.ns_veg.frootn_storage_patch),
        ("frootn_xfer",      H.ns_veg.frootn_xfer_patch,      Sd.ns_veg.frootn_xfer_patch),
        ("livestemn",        H.ns_veg.livestemn_patch,        Sd.ns_veg.livestemn_patch),
        ("livestemn_storage",H.ns_veg.livestemn_storage_patch,Sd.ns_veg.livestemn_storage_patch),
        ("livestemn_xfer",   H.ns_veg.livestemn_xfer_patch,   Sd.ns_veg.livestemn_xfer_patch),
        ("deadstemn",        H.ns_veg.deadstemn_patch,        Sd.ns_veg.deadstemn_patch),
        ("deadstemn_storage",H.ns_veg.deadstemn_storage_patch,Sd.ns_veg.deadstemn_storage_patch),
        ("livecrootn",       H.ns_veg.livecrootn_patch,       Sd.ns_veg.livecrootn_patch),
        ("deadcrootn",       H.ns_veg.deadcrootn_patch,       Sd.ns_veg.deadcrootn_patch),
        ("cropseedn_deficit",H.ns_veg.cropseedn_deficit_patch,Sd.ns_veg.cropseedn_deficit_patch),
        ("seedn_grc",        H.ns_veg.seedn_grc,              Sd.ns_veg.seedn_grc),
        ("reproductiven",    H.ns_veg.reproductiven_patch,    Sd.ns_veg.reproductiven_patch),
        ("reproductiven_xfer",H.ns_veg.reproductiven_xfer_patch,Sd.ns_veg.reproductiven_xfer_patch),
        ("free_retransn",    H.nf_veg.free_retransn_to_npool_patch, Sd.nf_veg.free_retransn_to_npool_patch),
        ("decomp_vr_col",    H.ns_soil.decomp_npools_vr_col,  Sd.ns_soil.decomp_npools_vr_col),
        ("sourcesink_col",   H.nf_soil.decomp_npools_sourcesink_col, Sd.nf_soil.decomp_npools_sourcesink_col),
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
    println(nfail == 0 ? "  WHOLE n_state_update1! MATCHES CPU ON $name ($FT) ✓" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
