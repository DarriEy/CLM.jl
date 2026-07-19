# ==========================================================================
# gpu_validate_dyn_cnbal_patch_carbon_e2e.jl — whole-function GPU parity for
# dynamic_patch_adjustments_carbon! (the dyn_cnbal carbon patch-adjustment path):
# compute_seed_amounts! + the ~25 conservative update_patch_state! calls + the
# deadstem partition-flux-by-type, all composed on-device.
#
# Mirrors test_dyn_cons_biogeochem's 3-patch subgrid (shrink / grow / initiate).
# Runs on CPU, adapts the CN carbon state + updater + patch + flux vectors to
# Metal/Float32 (pftcon stays host; the orchestrator copies its arrays), runs the
# SAME call on the device, and compares the mutated patch C pools + dwt fluxes.
#
#   julia --project=scripts scripts/gpu_validate_dyn_cnbal_patch_carbon_e2e.jl
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

const NP = 3
const NC = 3
const NG = 1
const PWT_OLD = [0.40, 0.10, 0.00]
const PWT_NEW = [0.20, 0.30, 0.25]

const CPOOL_FIELDS = (:leafc_patch, :leafc_storage_patch, :leafc_xfer_patch,
    :frootc_patch, :frootc_storage_patch, :frootc_xfer_patch,
    :livestemc_patch, :livestemc_storage_patch, :livestemc_xfer_patch,
    :deadstemc_patch, :deadstemc_storage_patch, :deadstemc_xfer_patch,
    :livecrootc_patch, :livecrootc_storage_patch, :livecrootc_xfer_patch,
    :deadcrootc_patch, :deadcrootc_storage_patch, :deadcrootc_xfer_patch,
    :gresp_storage_patch, :gresp_xfer_patch, :cpool_patch,
    :xsmrpool_patch, :ctrunc_patch)

function make_pftcon(n, ndecomp_pools)
    p = CLM.PftconType()
    p.c3psn = ones(n); p.evergreen = zeros(n); p.woody = zeros(n)
    p.leafcn = fill(25.0, n); p.deadwdcn = fill(100.0, n)
    p.noveg = fill(typemax(Int), n); p.pconv = fill(0.5, n)
    p.fr_f = zeros(n, ndecomp_pools)
    for j in 1:n; p.fr_f[j,1]=0.6; p.fr_f[j,2]=0.4; end
    return p
end

# Build the host fixture: (updater, patch, col, cs, fluxes...). set_old/new are
# host per-step setup; cs pools filled with distinct nonzero values.
function build()
    CLM.varpar_init!(CLM.varpar, 1, 15, 2, 5); CLM.varcon_init!()
    bounds = CLM.BoundsType(begg=1, endg=NG, begl=1, endl=1, begc=1, endc=NC,
                            begp=1, endp=NP, level=CLM.BOUNDS_LEVEL_PROC)
    pch = CLM.PatchData{Float64}()
    pch.column = collect(1:NP); pch.landunit = fill(1, NP); pch.gridcell = fill(1, NP)
    pch.itype = [7, 13, 8]; pch.wtgcell = copy(PWT_OLD)
    col = CLM.ColumnData{Float64}(); col.gridcell = fill(1, NC); col.wtgcell = copy(PWT_OLD)
    updater = CLM.PatchStateUpdater(bounds)
    CLM.set_old_weights!(updater, bounds, pch, col)
    pch.wtgcell .= PWT_NEW; col.wtgcell .= PWT_NEW
    CLM.set_new_weights!(updater, bounds, pch)

    cs = CLM.CNVegCarbonStateData{Float64}(); CLM.cnveg_carbon_state_init!(cs, NP, NC, NG)
    v = 1.0
    for f in CPOOL_FIELDS
        arr = getproperty(cs, f)
        for p in 1:NP; arr[p] = v; v += 1.0; end
    end
    return bounds, pch, col, updater, cs
end

# the 8 per-patch flux/seed accumulators the carbon adjustment writes into
newfluxes() = (conv=zeros(NP), wood=zeros(NP), crop=zeros(NP),
               froot=zeros(NP), livecroot=zeros(NP), deadcroot=zeros(NP),
               leafseed=zeros(NP), deadstemseed=zeros(NP))

function run!(cs, updater, bounds, pch, pftcon, mask, filterp, fx)
    CLM.dynamic_patch_adjustments_carbon!(cs, updater, bounds, pch, pftcon, mask, filterp;
        leafc_seed=5.0, deadstemc_seed=2.0,
        conv_cflux=fx.conv, wood_product_cflux=fx.wood, crop_product_cflux=fx.crop,
        dwt_frootc_to_litter=fx.froot, dwt_livecrootc_to_litter=fx.livecroot,
        dwt_deadcrootc_to_litter=fx.deadcroot,
        dwt_leafc_seed=fx.leafseed, dwt_deadstemc_seed=fx.deadstemseed)
end

function main(backend)
    println("=" ^ 72)
    println("WHOLE-FN GPU parity for dynamic_patch_adjustments_carbon! (dyn_cnbal C path)")
    println("=" ^ 72)
    if backend === nothing
        println("  No GPU backend — nothing to validate."); return 0
    end
    name, _, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n", name, FT)

    pftcon = make_pftcon(15, 3)   # host (orchestrator copies pconv/c3psn/... to backend)

    # CPU reference
    bH, pchH, _, updH, csH = build()
    fxH = newfluxes(); maskH = trues(NP); filterH = collect(1:NP)
    run!(csH, updH, bH, pchH, pftcon, maskH, filterH, fxH)

    # Device
    mf(x) = CLM.Adapt.adapt(device_array_type(), CLM.Adapt.adapt(_F32(), x))
    mi(x) = device_array_type()(collect(Int32.(x)))
    bD, pchD, _, updD, csD = build()
    csD = mf(csD); updD = mf(updD); pchD = mf(pchD)
    fxD = map(mf, newfluxes())
    if !(csD.leafc_patch isa device_array_type())
        println("  BLOCKED: CN carbon state did not move to the device."); return 2
    end
    # mask stays a host BitVector (the orchestrator collects+copies it to backend)
    run!(csD, updD, bD, pchD, pftcon, trues(NP), mi(collect(1:NP)), fxD)

    nfail = 0
    for f in CPOOL_FIELDS
        a = getproperty(csH, f); b = getproperty(csD, f)
        dd = reldiff(a, b); ok = dd < 1f-3
        ok || (@printf("  [FAIL] cs.%-24s rel = %.3e\n", f, dd); nfail += 1)
    end
    @printf("  [%s] %d/%d patch C pools match\n", nfail == 0 ? "PASS" : "FAIL",
            length(CPOOL_FIELDS) - nfail, length(CPOOL_FIELDS))
    for (nm, hv, dv) in (("conv_cflux", fxH.conv, fxD.conv),
                         ("wood_product_cflux", fxH.wood, fxD.wood),
                         ("dwt_frootc_to_litter", fxH.froot, fxD.froot),
                         ("dwt_livecrootc_to_litter", fxH.livecroot, fxD.livecroot),
                         ("dwt_deadcrootc_to_litter", fxH.deadcroot, fxD.deadcroot),
                         ("dwt_leafc_seed", fxH.leafseed, fxD.leafseed),
                         ("dwt_deadstemc_seed", fxH.deadstemseed, fxD.deadstemseed))
        dd = reldiff(hv, dv); ok = dd < 1f-3
        @printf("  [%s] %-26s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, dd)
        ok || (nfail += 1)
    end
    println()
    println(nfail == 0 ? "  WHOLE dynamic_patch_adjustments_carbon! MATCHES CPU ON $name ($FT)" :
                         "  DIVERGENCE — investigate ($nfail failed).")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
