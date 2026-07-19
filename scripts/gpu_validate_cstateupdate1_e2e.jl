# ==========================================================================
# gpu_validate_cstateupdate1_e2e.jl — end-to-end GPU parity for the WHOLE
# c_state_update1! BGC driver (the Phase B C/N state-update cascade entry point).
#
# Builds a small multi-column / multi-patch instance mirroring
# test/test_c_state_update1.jl, runs c_state_update1! on the CPU, converts every
# state struct to Float32 + adapts to the GPU, runs the SAME call on the device, and
# compares the mutated outputs field-by-field. The scenario deliberately exercises
# the branchy paths in one shot: a non-woody patch, a woody patch, two crop patches
# (one harvested), a FATES column (column-kernel skip), and a non-trivial
# decomposition cascade (donor/receiver RMW into decomp_cpools_sourcesink_col).
#
#   julia --project=scripts scripts/gpu_validate_cstateupdate1_e2e.jl
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
# the scalar/index params c_state_update1! needs.
function build()
    cs_veg = CLM.CNVegCarbonStateData()
    CLM.cnveg_carbon_state_init!(cs_veg, NP, NC, NG; nrepr=NREPR)
    cf_veg = CLM.CNVegCarbonFluxData()
    CLM.cnveg_carbon_flux_init!(cf_veg, NP, NC, NG; nrepr=NREPR,
                                nlevdecomp_full=NLEVDECOMP, ndecomp_pools=NDECOMP_POOLS)
    cf_soil = CLM.SoilBiogeochemCarbonFluxData()
    CLM.soil_bgc_carbon_flux_init!(cf_soil, NC, NLEVDECOMP, NDECOMP_POOLS, NCASCADE)

    # --- Carbon-state pools (known nonzero so every flux moves something) ---
    cs_veg.cpool_patch              .= 100.0
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
    cs_veg.xsmrpool_loss_patch      .= 0.0
    cs_veg.gresp_storage_patch      .= 5.0
    cs_veg.gresp_xfer_patch         .= 2.0
    cs_veg.cropseedc_deficit_patch  .= -5.0
    for k in 1:NREPR
        cs_veg.reproductivec_patch[:, k]         .= 10.0
        cs_veg.reproductivec_storage_patch[:, k] .= 2.0
        cs_veg.reproductivec_xfer_patch[:, k]    .= 1.0
    end

    # --- Carbon fluxes (mirror test_c_state_update1.jl + extras for the resp/crop
    #     and cascade branches). Set a broad swath nonzero so all branches move. ---
    cf_veg.leafc_xfer_to_leafc_patch         .= 2.0e-7
    cf_veg.frootc_xfer_to_frootc_patch       .= 1.0e-7
    cf_veg.livestemc_xfer_to_livestemc_patch .= 1.5e-7
    cf_veg.deadstemc_xfer_to_deadstemc_patch .= 1.0e-7
    cf_veg.livecrootc_xfer_to_livecrootc_patch .= 0.8e-7
    cf_veg.deadcrootc_xfer_to_deadcrootc_patch .= 0.5e-7
    cf_veg.leafc_to_litter_patch             .= 3.0e-7
    cf_veg.frootc_to_litter_patch            .= 2.0e-7
    cf_veg.livestemc_to_deadstemc_patch      .= 1.0e-8
    cf_veg.livecrootc_to_deadcrootc_patch    .= 1.0e-8
    cf_veg.livestemc_to_litter_patch         .= 1.0e-9
    cf_veg.livestemc_to_biofuelc_patch       .= 1.0e-9
    cf_veg.livestemc_to_removedresiduec_patch .= 1.0e-9
    cf_veg.leafc_to_biofuelc_patch           .= 1.0e-9
    cf_veg.leafc_to_removedresiduec_patch    .= 1.0e-9
    cf_veg.crop_seedc_to_leaf_patch          .= 2.0e-9
    for k in 1:NREPR
        cf_veg.repr_grainc_to_food_patch[:, k]          .= 1.0e-9
        cf_veg.repr_grainc_to_seed_patch[:, k]          .= 1.0e-9
        cf_veg.repr_structurec_to_cropprod_patch[:, k]  .= 1.0e-9
        cf_veg.repr_structurec_to_litter_patch[:, k]    .= 1.0e-9
        cf_veg.reproductivec_xfer_to_reproductivec_patch[:, k] .= 1.0e-9
    end
    cf_veg.cpool_to_xsmrpool_patch          .= 1.0e-7
    cf_veg.leaf_curmr_patch                 .= 5.0e-8
    cf_veg.froot_curmr_patch                .= 3.0e-8
    cf_veg.livestem_curmr_patch             .= 2.0e-8
    cf_veg.livecroot_curmr_patch            .= 1.0e-8
    cf_veg.cpool_to_resp_patch              .= 1.0e-9
    cf_veg.soilc_change_patch               .= 1.0e-9
    cf_veg.leaf_xsmr_patch                  .= 1.0e-8
    cf_veg.froot_xsmr_patch                 .= 0.5e-8
    cf_veg.livestem_xsmr_patch              .= 0.3e-8
    cf_veg.livecroot_xsmr_patch             .= 0.2e-8
    for k in 1:NREPR
        cf_veg.reproductive_curmr_patch[:, k] .= 1.0e-9
        cf_veg.reproductive_xsmr_patch[:, k]  .= 1.0e-9
    end
    cf_veg.cpool_to_leafc_patch             .= 5.0e-7
    cf_veg.cpool_to_leafc_storage_patch     .= 2.0e-7
    cf_veg.cpool_to_frootc_patch            .= 3.0e-7
    cf_veg.cpool_to_frootc_storage_patch    .= 1.0e-7
    cf_veg.cpool_to_leafc_resp_patch        .= 1.0e-9
    cf_veg.cpool_to_leafc_storage_resp_patch .= 1.0e-9
    cf_veg.cpool_to_frootc_resp_patch       .= 1.0e-9
    cf_veg.cpool_to_frootc_storage_resp_patch .= 1.0e-9
    cf_veg.cpool_to_livestemc_patch         .= 2.0e-7
    cf_veg.cpool_to_livestemc_storage_patch .= 1.0e-7
    cf_veg.cpool_to_deadstemc_patch         .= 1.5e-7
    cf_veg.cpool_to_deadstemc_storage_patch .= 0.5e-7
    cf_veg.cpool_to_livecrootc_patch        .= 1.0e-7
    cf_veg.cpool_to_livecrootc_storage_patch .= 0.5e-7
    cf_veg.cpool_to_deadcrootc_patch        .= 0.5e-7
    cf_veg.cpool_to_deadcrootc_storage_patch .= 0.3e-7
    cf_veg.cpool_to_livecrootc_resp_patch   .= 1.0e-9
    cf_veg.cpool_to_livecrootc_storage_resp_patch .= 1.0e-9
    cf_veg.cpool_to_livestemc_resp_patch    .= 1.0e-9
    cf_veg.cpool_to_livestemc_storage_resp_patch .= 1.0e-9
    for k in 1:NREPR
        cf_veg.cpool_to_reproductivec_patch[:, k]         .= 1.0e-9
        cf_veg.cpool_to_reproductivec_storage_patch[:, k] .= 1.0e-9
    end
    cf_veg.cpool_leaf_gr_patch              .= 1.0e-8
    cf_veg.cpool_froot_gr_patch             .= 0.8e-8
    cf_veg.cpool_livestem_gr_patch          .= 0.5e-8
    cf_veg.cpool_deadstem_gr_patch          .= 0.3e-8
    cf_veg.cpool_livecroot_gr_patch         .= 0.2e-8
    cf_veg.cpool_deadcroot_gr_patch         .= 0.1e-8
    for k in 1:NREPR
        cf_veg.cpool_reproductive_gr_patch[:, k]          .= 1.0e-9
        cf_veg.cpool_reproductive_storage_gr_patch[:, k]  .= 1.0e-9
        cf_veg.transfer_reproductive_gr_patch[:, k]       .= 1.0e-9
    end
    cf_veg.transfer_leaf_gr_patch           .= 0.5e-8
    cf_veg.transfer_froot_gr_patch          .= 0.3e-8
    cf_veg.transfer_livestem_gr_patch       .= 0.2e-8
    cf_veg.transfer_deadstem_gr_patch       .= 0.1e-8
    cf_veg.transfer_livecroot_gr_patch      .= 0.1e-8
    cf_veg.transfer_deadcroot_gr_patch      .= 0.05e-8
    cf_veg.cpool_leaf_storage_gr_patch      .= 0.5e-8
    cf_veg.cpool_froot_storage_gr_patch     .= 0.3e-8
    cf_veg.cpool_livestem_storage_gr_patch  .= 0.2e-8
    cf_veg.cpool_deadstem_storage_gr_patch  .= 0.1e-8
    cf_veg.cpool_livecroot_storage_gr_patch .= 0.1e-8
    cf_veg.cpool_deadcroot_storage_gr_patch .= 0.05e-8
    cf_veg.cpool_to_gresp_storage_patch     .= 1.0e-8
    cf_veg.leafc_storage_to_xfer_patch      .= 1.0e-9
    cf_veg.frootc_storage_to_xfer_patch     .= 1.0e-9
    cf_veg.gresp_storage_to_xfer_patch      .= 1.0e-9
    cf_veg.livestemc_storage_to_xfer_patch  .= 1.0e-9
    cf_veg.deadstemc_storage_to_xfer_patch  .= 1.0e-9
    cf_veg.livecrootc_storage_to_xfer_patch .= 1.0e-9
    cf_veg.deadcrootc_storage_to_xfer_patch .= 1.0e-9
    for k in 1:NREPR
        cf_veg.reproductivec_storage_to_xfer_patch[:, k] .= 1.0e-9
    end
    cf_veg.xsmrpool_to_atm_patch            .= 0.0
    cf_veg.phenology_c_to_litr_c_col        .= 1.0e-7

    # --- Soil BGC carbon flux (column-kernel inputs; nonzero cascade) ---
    cf_soil.decomp_cpools_sourcesink_col    .= 0.0
    cf_soil.decomp_cascade_hr_vr_col        .= 2.0e-8
    cf_soil.decomp_cascade_ctransfer_vr_col .= 3.0e-8

    # --- Decomposition cascade (donor/receiver; 0 = terminal) ---
    cascade_donor_pool    = [1, 2, 3, 1, 2]
    cascade_receiver_pool = [2, 3, 0, 4, 4]

    # --- Patch/column topology + branch selectors ---
    # patches: 1 non-woody, 2 woody, 3 non-woody, 4 woody, 5 crop(harvested), 6 crop
    ivt          = [1, 2, 1, 2, NPCROPMIN, NPCROPMIN]
    woody        = zeros(Float64, 80); woody[2] = 1.0
    harvdate     = [999, 999, 999, 999, 100, 999]   # patch 5 harvested
    patch_column = [1, 1, 2, 2, 3, 4]
    col_is_fates = [false, false, false, true]      # column 4 is FATES (col-kernel skip)
    mask_soilc   = trues(NC)
    mask_soilp   = trues(NP)

    return (; cs_veg, cf_veg, cf_soil, cascade_donor_pool, cascade_receiver_pool,
              ivt, woody, harvdate, patch_column, col_is_fates, mask_soilc, mask_soilp)
end

run_csu1!(S, m_c, m_p, donor, recv, harv, fates, ivt, woody, patchcol, dt) =
    CLM.c_state_update1!(S.cs_veg, S.cf_veg, S.cf_soil;
        mask_soilc=m_c, mask_soilp=m_p, bounds_col=1:NC, bounds_patch=1:NP,
        patch_column=patchcol, ivt=ivt, woody=woody,
        cascade_donor_pool=donor, cascade_receiver_pool=recv,
        harvdate=harv, col_is_fates=fates,
        nlevdecomp=NLEVDECOMP, ndecomp_cascade_transitions=NCASCADE,
        i_litr_min=I_LITR_MIN, i_litr_max=I_LITR_MAX, i_cwd=I_CWD,
        npcropmin=NPCROPMIN, nrepr=NREPR,
        repr_grain_min=1, repr_grain_max=NREPR,
        repr_structure_min=1, repr_structure_max=0,
        carbon_resp_opt=1, dribble_crophrv_xsmrpool_2atm=false, dt=dt)

function main(backend)
    println("=" ^ 70)
    println("END-TO-END GPU parity for c_state_update1! (BGC C-state cascade)")
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
    cs_d  = ad(B.cs_veg); cf_d = ad(B.cf_veg); cfs_d = ad(B.cf_soil)
    Sd = (; cs_veg=cs_d, cf_veg=cf_d, cf_soil=cfs_d)

    if !(cs_d.cpool_patch isa device_array_type())
        println("  BLOCKED: a state struct did not move to the device under adapt.")
        return 2
    end

    dt = FT(1800.0)
    dmask(m) = to(collect(Bool, m))

    # CPU
    run_csu1!(H, H.mask_soilc, H.mask_soilp, H.cascade_donor_pool, H.cascade_receiver_pool,
              H.harvdate, H.col_is_fates, H.ivt, H.woody, H.patch_column, 1800.0)
    # Device
    run_csu1!(Sd, dmask(B.mask_soilc), dmask(B.mask_soilp),
              to(B.cascade_donor_pool), to(B.cascade_receiver_pool),
              to(B.harvdate), dmask(B.col_is_fates), to(B.ivt),
              to(Float32.(B.woody)), to(B.patch_column), dt)

    checks = [
        ("cpool",            H.cs_veg.cpool_patch,            Sd.cs_veg.cpool_patch),
        ("xsmrpool",         H.cs_veg.xsmrpool_patch,         Sd.cs_veg.xsmrpool_patch),
        ("xsmrpool_loss",    H.cs_veg.xsmrpool_loss_patch,    Sd.cs_veg.xsmrpool_loss_patch),
        ("leafc",            H.cs_veg.leafc_patch,            Sd.cs_veg.leafc_patch),
        ("leafc_storage",    H.cs_veg.leafc_storage_patch,    Sd.cs_veg.leafc_storage_patch),
        ("leafc_xfer",       H.cs_veg.leafc_xfer_patch,       Sd.cs_veg.leafc_xfer_patch),
        ("frootc",           H.cs_veg.frootc_patch,           Sd.cs_veg.frootc_patch),
        ("frootc_storage",   H.cs_veg.frootc_storage_patch,   Sd.cs_veg.frootc_storage_patch),
        ("frootc_xfer",      H.cs_veg.frootc_xfer_patch,      Sd.cs_veg.frootc_xfer_patch),
        ("livestemc",        H.cs_veg.livestemc_patch,        Sd.cs_veg.livestemc_patch),
        ("livestemc_storage",H.cs_veg.livestemc_storage_patch,Sd.cs_veg.livestemc_storage_patch),
        ("livestemc_xfer",   H.cs_veg.livestemc_xfer_patch,   Sd.cs_veg.livestemc_xfer_patch),
        ("deadstemc",        H.cs_veg.deadstemc_patch,        Sd.cs_veg.deadstemc_patch),
        ("livecrootc",       H.cs_veg.livecrootc_patch,       Sd.cs_veg.livecrootc_patch),
        ("deadcrootc",       H.cs_veg.deadcrootc_patch,       Sd.cs_veg.deadcrootc_patch),
        ("gresp_storage",    H.cs_veg.gresp_storage_patch,    Sd.cs_veg.gresp_storage_patch),
        ("gresp_xfer",       H.cs_veg.gresp_xfer_patch,       Sd.cs_veg.gresp_xfer_patch),
        ("cropseedc_deficit",H.cs_veg.cropseedc_deficit_patch,Sd.cs_veg.cropseedc_deficit_patch),
        ("reproductivec",    H.cs_veg.reproductivec_patch,    Sd.cs_veg.reproductivec_patch),
        ("reproductivec_xfer",H.cs_veg.reproductivec_xfer_patch,Sd.cs_veg.reproductivec_xfer_patch),
        ("xsmrpool_to_atm",  H.cf_veg.xsmrpool_to_atm_patch,  Sd.cf_veg.xsmrpool_to_atm_patch),
        ("sourcesink_col",   H.cf_soil.decomp_cpools_sourcesink_col, Sd.cf_soil.decomp_cpools_sourcesink_col),
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
    println(nfail == 0 ? "  WHOLE c_state_update1! MATCHES CPU ON $name ($FT) ✓" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
