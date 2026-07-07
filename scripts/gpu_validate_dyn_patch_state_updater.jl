# ==========================================================================
# gpu_validate_dyn_patch_state_updater.jl — Metal parity for the shared
# conservative patch-state updater infra (dyn_subgrid): set_old_weights!,
# set_new_weights!, update_patch_state! (with seed + col/grc fluxes), and
# update_patch_state_partition_flux_by_type!.
#
# Mirrors test_dyn_patch_state_updater's 4-patch subgrid (shrinking / growing /
# initiating / unchanged). Runs the full updater sequence on the CPU, adapts the
# updater + patch/col + state/flux arrays to Metal/Float32, runs the SAME
# sequence on the device, and compares.
#
#   julia --project=scripts scripts/gpu_validate_dyn_patch_state_updater.jl
# ==========================================================================

using CLM
using Printf
import Metal
include(joinpath(@__DIR__, "gpu_backends.jl"))

struct _F32 end
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:AbstractFloat}) = Float32.(x)
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:Integer}) = x
CLM.Adapt.adapt_storage(::_F32, x::Float64) = Float32(x)

reldiff(a, b) = begin
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + max(abs(Float64(A[i])), abs(Float64(B[i])))))
    end
    m
end

const NP = 4
const NC = 4
const PWT_OLD = [0.40, 0.10, 0.00, 0.25]
const PWT_NEW = [0.20, 0.30, 0.25, 0.25]

bounds() = CLM.BoundsType(begg=1, endg=1, begl=1, endl=1, begc=1, endc=NC,
                          begp=1, endp=NP, level=CLM.BOUNDS_LEVEL_PROC)

function make()
    pch = CLM.PatchData{Float64}()
    pch.column  = collect(1:NP)
    pch.itype   = [7, 13, 1, 0]
    pch.wtgcell = copy(PWT_OLD)
    col = CLM.ColumnData{Float64}()
    col.wtgcell = copy(PWT_OLD)
    return pch, col
end

# Full updater sequence on whatever backend the arrays live on. `dev` adapts
# per-call scratch/inputs (filter, frac) to the same backend as the updater.
function run_seq!(updater, b, pch, col, var, flux_col, flux_grc, seed, seed_add,
                  flux1, flux2, filterp, frac)
    CLM.set_old_weights!(updater, b, pch, col)
    pch.wtgcell .= (pch.wtgcell isa Metal.MtlArray ? Metal.MtlArray(Float32.(PWT_NEW)) : PWT_NEW)
    CLM.set_new_weights!(updater, b, pch)
    CLM.update_patch_state!(updater, b, pch, filterp, var;
        flux_out_col_area=flux_col, flux_out_grc_area=flux_grc,
        seed=seed, seed_addition=seed_add)
    var2 = copy(var)  # partition uses a fresh var to keep the flux math independent
    CLM.update_patch_state_partition_flux_by_type!(updater, b, pch, filterp, frac,
        var2, flux1, flux2)
    return nothing
end

function main(backend)
    println("=" ^ 72)
    println("Metal parity for the dyn_subgrid patch-state updater (conservative infra)")
    println("=" ^ 72)
    if backend === nothing
        println("  No GPU backend — nothing to validate."); return 0
    end
    name, _, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n", name, FT)

    b = bounds()

    # --- CPU reference ---
    pchH, colH = make()
    uH = CLM.PatchStateUpdater(b)
    varH = [10.0, 20.0, 0.0, 5.0]
    fcolH = zeros(NP); fgrcH = zeros(NP); seedH = [0.0, 100.0, 50.0, 0.0]
    saddH = zeros(NP); f1H = zeros(NP); f2H = zeros(NP)
    fracH = fill(0.25, CLM.MXPFT + 1); fracH[pchH.itype[1] + 1] = 0.75
    filterH = collect(1:NP)
    run_seq!(uH, b, pchH, colH, varH, fcolH, fgrcH, seedH, saddH, f1H, f2H, filterH, fracH)

    # --- Device run ---
    mf(x) = CLM.Adapt.adapt(Metal.MtlArray, CLM.Adapt.adapt(_F32(), x))
    mi(x) = Metal.MtlArray(collect(Int32.(x)))
    pchD, colD = make()
    pchD = mf(pchD); colD = mf(colD)
    uD = mf(CLM.PatchStateUpdater(b))
    varD  = mf([10.0, 20.0, 0.0, 5.0])
    fcolD = mf(zeros(NP)); fgrcD = mf(zeros(NP)); seedD = mf([0.0, 100.0, 50.0, 0.0])
    saddD = mf(zeros(NP)); f1D = mf(zeros(NP)); f2D = mf(zeros(NP))
    fracD = mf(fill(0.25, CLM.MXPFT + 1)); # set the 0.75 after adapt (device scalar set is awkward)
    fracH2 = fill(0.25, CLM.MXPFT + 1); fracH2[pchH.itype[1] + 1] = 0.75
    fracD = mf(fracH2)
    filterD = mi(collect(1:NP))

    if !(uD.dwt isa Metal.MtlArray)
        println("  BLOCKED: PatchStateUpdater did not move to the device."); return 2
    end

    run_seq!(uD, b, pchD, colD, varD, fcolD, fgrcD, seedD, saddD, f1D, f2D, filterD, fracD)

    checks = [
        ("updater.dwt",                  uH.dwt,                  uD.dwt),
        ("updater.growing_old_fraction", uH.growing_old_fraction, uD.growing_old_fraction),
        ("updater.growing_new_fraction", uH.growing_new_fraction, uD.growing_new_fraction),
        ("updater.pwtgcell_old",         uH.pwtgcell_old,         uD.pwtgcell_old),
        ("updater.cwtgcell_old",         uH.cwtgcell_old,         uD.cwtgcell_old),
        ("var (grow/seed/dilute)",       varH, varD),
        ("flux_out_col_area",            fcolH, fcolD),
        ("flux_out_grc_area",            fgrcH, fgrcD),
        ("seed_addition",                saddH, saddD),
        ("flux1_out (partition)",        f1H, f1D),
        ("flux2_out (partition)",        f2H, f2D),
    ]
    nfail = 0
    for (nm, a, c) in checks
        dd = reldiff(a, c); ok = dd < 1f-3
        @printf("  [%s] %-30s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, dd)
        ok || (nfail += 1)
    end
    println()
    println(nfail == 0 ? "  patch-state updater infra MATCHES CPU ON $name ($FT)" :
                         "  DIVERGENCE — investigate ($nfail failed).")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
