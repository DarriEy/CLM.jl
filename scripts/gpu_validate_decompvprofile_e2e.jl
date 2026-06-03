# ==========================================================================
# gpu_validate_decompvprofile_e2e.jl — end-to-end Metal parity for the WHOLE
# soil_bgc_vertical_profile! decomposition vertical-profile driver.
#
# soil_bgc_vertical_profile! runs whole-function on Metal: per-patch and
# per-column profile kernels, plus the single TRUE patch->column scatter of
# the root fraction (`col_cinput_rootfr[col,j] += cinput_rootfr[p,j]*wtcol[p]`,
# col = patch_column[p]) via _scatter_add!. This harness builds the SAME inputs
# as test/test_decomp_vertical_profile.jl at Float64, runs the function on the
# CPU, adapts every state struct (Float64 -> Float32 device arrays, Int/Bool
# preserved, BitVector -> device Vector{Bool}) with a down-converting Metal
# adaptor, runs the SAME call on the device, and compares the mutated outputs.
#
# To actually exercise the patch->column scatter, the setup uses >= 2 patches
# mapping to ONE column (with weights summing over the column) and compares
# col_cinput_rootfr indirectly via nfixation_prof_col (the column-native
# profile derived from col_cinput_rootfr) AND the direct per-output profiles.
#
#   julia --project=scripts scripts/gpu_validate_decompvprofile_e2e.jl
#
# NOTE (banked lesson): reldiff PASSES when both sides are NaN. main() asserts
# the CPU reference fields are finite before trusting the parity numbers.
# ==========================================================================

using CLM
using Printf
import Metal
const Adapt = CLM.Adapt   # Adapt isn't in scripts/Project.toml; reach it via CLM
include(joinpath(@__DIR__, "gpu_backends.jl"))

function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(A[i] - B[i]) / (1.0 + max(abs(A[i]), abs(B[i]))))
    end
    return m
end

allfinite(a) = all(isfinite, Array(a))

# --------------------------------------------------------------------------
# MetalF32 down-convert adaptor: build state at Float64, move to the device as
# Float32 (Metal has no f64). Integer index/filter vectors keep their Int
# eltype (coercing them to Float32 would break device indexing); Bool/BitArray
# masks become device Vector{Bool}.
# --------------------------------------------------------------------------
struct MetalF32 end
Adapt.adapt_storage(::MetalF32, x::AbstractArray{Float64}) = Metal.MtlArray(Float32.(x))
Adapt.adapt_storage(::MetalF32, x::AbstractArray{Float32}) = Metal.MtlArray(x)
Adapt.adapt_storage(::MetalF32, x::AbstractArray{<:Integer}) = Metal.MtlArray(collect(x))
Adapt.adapt_storage(::MetalF32, x::AbstractArray{Bool}) = Metal.MtlArray(collect(Bool, x))
Adapt.adapt_storage(::MetalF32, x::BitArray) = Metal.MtlArray(collect(Bool, x))
adF32(x) = Adapt.adapt(MetalF32(), x)

# --------------------------------------------------------------------------
# Build a Float64 instance mirroring test/test_decomp_vertical_profile.jl's
# make_test_setup, but with TWO vegetated patches mapping to ONE column so the
# patch->column scatter is genuinely exercised, plus a partial active layer and
# distinct per-patch root distributions to make the branches live.
# --------------------------------------------------------------------------
function build()
    nc = 1; np = 2
    nlevdecomp = 4; nlevdecomp_full = 4

    dz = fill(0.1, nlevdecomp)
    zs = [(j - 0.5) * 0.1 for j in 1:nlevdecomp]

    bgc_state = CLM.SoilBiogeochemStateData()
    bgc_state.leaf_prof_patch    = zeros(np, nlevdecomp_full)
    bgc_state.froot_prof_patch   = zeros(np, nlevdecomp_full)
    bgc_state.croot_prof_patch   = zeros(np, nlevdecomp_full)
    bgc_state.stem_prof_patch    = zeros(np, nlevdecomp_full)
    bgc_state.nfixation_prof_col = zeros(nc, nlevdecomp_full)
    bgc_state.ndep_prof_col      = zeros(nc, nlevdecomp_full)

    active_layer = CLM.ActiveLayerData()
    active_layer.altmax_lastyear_indx_col = fill(nlevdecomp, nc)

    soilstate = CLM.SoilStateData()
    # Patch 1: roots in layer 1; patch 2: roots in layer 3 (distinct profiles).
    rootfr = zeros(np, nlevdecomp_full)
    rootfr[1, 1] = 1.0
    rootfr[2, 3] = 1.0
    soilstate.crootfr_patch = rootfr

    col_data = CLM.ColumnData()
    col_data.nbedrock = fill(nlevdecomp, nc)
    col_data.is_fates = fill(false, nc)

    patch_data = CLM.PatchData()
    patch_data.itype  = fill(1, np)        # all vegetated
    patch_data.column = fill(1, np)        # BOTH patches -> column 1 (scatter!)
    patch_data.wtcol  = [0.6, 0.4]         # weights over the shared column

    mask_soilc = trues(nc)
    mask_vegp  = trues(np)

    return (; bgc_state, active_layer, soilstate, col=col_data, patch=patch_data,
              mask_soilc, mask_vegp, nlevdecomp, nlevdecomp_full,
              dzsoi_decomp=dz, zsoi=zs)
end

run_dvp!(S) = CLM.soil_bgc_vertical_profile!(
    S.bgc_state, S.active_layer, S.soilstate, S.col, S.patch;
    mask_bgc_soilc = S.mask_soilc,
    mask_bgc_vegp  = S.mask_vegp,
    nlevdecomp      = S.nlevdecomp,
    nlevdecomp_full = S.nlevdecomp_full,
    dzsoi_decomp    = S.dzsoi_decomp,
    zsoi            = S.zsoi)

function main(backend)
    println("=" ^ 70)
    println("END-TO-END Metal parity for soil_bgc_vertical_profile! (whole driver)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU driver exercised by the suite).")
        return 0
    end
    name, _to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    H = build()          # CPU reference (Float64)
    B = build()          # to be adapted to device

    # Adapt state structs (Float64 -> Float32 device) + the loose vectors.
    Sd = (; bgc_state    = adF32(B.bgc_state),
            active_layer = adF32(B.active_layer),
            soilstate    = adF32(B.soilstate),
            col          = adF32(B.col),
            patch        = adF32(B.patch),
            mask_soilc   = adF32(B.mask_soilc),
            mask_vegp    = adF32(B.mask_vegp),
            nlevdecomp      = B.nlevdecomp,
            nlevdecomp_full = B.nlevdecomp_full,
            dzsoi_decomp    = adF32(B.dzsoi_decomp),
            zsoi            = adF32(B.zsoi))

    if !(Sd.bgc_state.nfixation_prof_col isa Metal.MtlArray)
        println("  BLOCKED: SoilBiogeochemStateData did not move to the device under adapt.")
        return 2
    end
    if !(Sd.patch.column isa Metal.MtlArray) || eltype(Sd.patch.column) != Int
        println("  BLOCKED: patch.column index vector lost its Int eltype on the device.")
        return 2
    end

    run_dvp!(H)    # CPU
    run_dvp!(Sd)   # device

    Hs = H.bgc_state; Ds = Sd.bgc_state

    cpu_fields = [
        ("leaf_prof",      Hs.leaf_prof_patch),
        ("froot_prof",     Hs.froot_prof_patch),
        ("croot_prof",     Hs.croot_prof_patch),
        ("stem_prof",      Hs.stem_prof_patch),
        ("nfixation_prof", Hs.nfixation_prof_col),
        ("ndep_prof",      Hs.ndep_prof_col),
    ]
    nbad = 0
    for (nm, a) in cpu_fields
        if !allfinite(a)
            @printf("  [NONFINITE] CPU reference %-16s = %s\n", nm, string(Array(a)))
            nbad += 1
        end
    end
    if nbad > 0
        println("\n  BLOCKED: CPU reference has non-finite fields — parity numbers untrustworthy.")
        return 2
    end

    # nfixation_prof_col is the column-native profile derived directly from the
    # scattered col_cinput_rootfr, so comparing it validates the patch->column
    # scatter accumulation across the two patches on the shared column.
    checks = [
        ("leaf_prof",      Hs.leaf_prof_patch,    Ds.leaf_prof_patch),
        ("froot_prof",     Hs.froot_prof_patch,   Ds.froot_prof_patch),
        ("croot_prof",     Hs.croot_prof_patch,   Ds.croot_prof_patch),
        ("stem_prof",      Hs.stem_prof_patch,    Ds.stem_prof_patch),
        ("nfixation_prof", Hs.nfixation_prof_col, Ds.nfixation_prof_col),
        ("ndep_prof",      Hs.ndep_prof_col,      Ds.ndep_prof_col),
    ]
    nfail = 0
    for (nm, a, b) in checks
        d = reldiff(a, b); ok = d < 1f-2
        @printf("  [%s] %-16s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end
    println()
    println(nfail == 0 ?
        "  WHOLE soil_bgc_vertical_profile! MATCHES CPU ON $name ($FT) ✓" :
        "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
