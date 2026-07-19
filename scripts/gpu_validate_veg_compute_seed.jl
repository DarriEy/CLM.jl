# ==========================================================================
# gpu_validate_veg_compute_seed.jl — GPU parity for compute_seed_amounts!
# (dyn_subgrid seed amounts for newly-grown patch area). Per-patch own-index
# kernel with leaf_proportions + species_type_multiplier inlined.
#
# pftcon_data stays HOST (the const global holds concrete host Vectors); the
# wrapper copies the needed pftcon arrays onto the output backend. The per-patch
# inputs/outputs are adapted to Metal/Float32 and compared to the CPU reference,
# for both a carbon species (C12) and nitrogen (N, uses leafcn/deadwdcn).
#
#   julia --project=scripts scripts/gpu_validate_veg_compute_seed.jl
# ==========================================================================

using CLM
using Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))

struct _F32 end
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:AbstractFloat}) = Float32.(x)
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:Integer}) = x
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{Bool}) = x
CLM.Adapt.adapt_storage(::_F32, x::Float64) = Float32(x)

reldiff(a, b) = begin
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + max(abs(Float64(A[i])), abs(Float64(B[i])))))
    end
    m
end

function make_pftcon(; npft=15)
    p = CLM.PftconType(); n = npft + 1
    p.c3psn = zeros(n); p.evergreen = zeros(n); p.woody = zeros(n)
    p.leafcn = fill(25.0, n); p.deadwdcn = fill(100.0, n)
    p.c3psn[2]=1.0; p.evergreen[2]=1.0; p.woody[2]=1.0; p.leafcn[2]=30.0; p.deadwdcn[2]=120.0
    p.c3psn[5]=1.0; p.evergreen[5]=1.0; p.woody[5]=1.0; p.leafcn[5]=40.0; p.deadwdcn[5]=200.0
    p.c3psn[8]=1.0; p.evergreen[8]=0.0; p.woody[8]=1.0; p.leafcn[8]=25.0; p.deadwdcn[8]=100.0
    p.c3psn[13]=1.0; p.evergreen[13]=0.0; p.woody[13]=0.0; p.leafcn[13]=20.0
    p.c3psn[15]=0.0; p.evergreen[15]=0.0; p.woody[15]=0.0; p.leafcn[15]=22.0
    return p
end

const NP = 6
# PFTs: bare(0), needleleaf-evergreen(1), broadleaf-evergreen-trop(4),
# broadleaf-decid(7), c3 grass(12), c4 grass(14).
mkpatch() = (pd = CLM.PatchData{Float64}(); pd.itype = [0, 1, 4, 7, 12, 14]; pd)
const LEAF   = [1.0, 5.0, 8.0, 0.0, 3.0, 2.0]
const LSTOR  = [0.5, 3.0, 2.0, 0.0, 1.0, 1.0]
const LXFER  = [0.2, 2.0, 1.0, 0.0, 0.5, 0.5]

function run!(species, mask, ch, ig, leaf, lstor, lxfer, sl, sls, slx, sd, patch, pc)
    CLM.compute_seed_amounts!(mask, 1:NP, patch, pc;
        species=species, leafc_seed=10.0, deadstemc_seed=4.0,
        leaf_patch=leaf, leaf_storage_patch=lstor, leaf_xfer_patch=lxfer,
        compute_here_patch=ch, ignore_current_state_patch=ig,
        seed_leaf_patch=sl, seed_leaf_storage_patch=sls, seed_leaf_xfer_patch=slx,
        seed_deadstem_patch=sd, noveg_val=0)
end

function main(backend)
    println("=" ^ 72)
    println("GPU parity for compute_seed_amounts! (dyn_subgrid seed amounts)")
    println("=" ^ 72)
    if backend === nothing
        println("  No GPU backend — nothing to validate."); return 0
    end
    name, _, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n", name, FT)

    pc = make_pftcon()
    mf(x) = CLM.Adapt.adapt(device_array_type(), CLM.Adapt.adapt(_F32(), x))
    mb(x) = device_array_type()(collect(x))
    nfail = 0

    for (sp, spname) in ((CLM.CN_SPECIES_C12, "C12"), (CLM.CN_SPECIES_N, "N"))
        # CPU reference
        maskH = trues(NP); chH = trues(NP); igH = falses(NP)
        slH=zeros(NP); slsH=zeros(NP); slxH=zeros(NP); sdH=zeros(NP)
        patchH = mkpatch()
        run!(sp, maskH, chH, igH, LEAF, LSTOR, LXFER, slH, slsH, slxH, sdH, patchH, pc)

        # Device (pftcon pc stays host; the wrapper copies its arrays to device)
        patchD = mf(mkpatch())
        slD=mf(zeros(NP)); slsD=mf(zeros(NP)); slxD=mf(zeros(NP)); sdD=mf(zeros(NP))
        run!(sp, mb(trues(NP)), mb(trues(NP)), mb(falses(NP)),
             mf(LEAF), mf(LSTOR), mf(LXFER), slD, slsD, slxD, sdD, patchD, pc)

        for (nm, a, b) in (("seed_leaf", slH, slD), ("seed_leaf_storage", slsH, slsD),
                           ("seed_leaf_xfer", slxH, slxD), ("seed_deadstem", sdH, sdD))
            dd = reldiff(a, b); ok = dd < 1f-3
            @printf("  [%s] %-8s %-18s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", spname, nm, dd)
            ok || (nfail += 1)
        end
    end

    println()
    println(nfail == 0 ? "  compute_seed_amounts! MATCHES CPU ON $name ($FT)" :
                         "  DIVERGENCE — investigate ($nfail failed).")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
