# ==========================================================================
# gpu_validate_dyn_cnbal_col_e2e.jl — whole-function Metal parity for
# dyn_cnbal_col! (the dyn_cnbal SOIL column-state conservation path): for every
# decomp pool / level it conservatively redistributes the vertically-resolved
# soil C/N pools across shrinking/growing columns (via the column-state updater)
# and depth-integrates the apparent per-column change — all composed on-device.
#
# Mirrors test_dyn_cons_biogeochem's 2-natveg-column area swap (c1 shrinks 0.2,
# c2 grows 0.2). Runs on CPU, adapts the soil BGC C/N state + updater + column to
# Metal/Float32, runs the SAME call on the device, and compares the redistributed
# vr pools + the depth-integrated adjustment diagnostics.
#
#   julia --project=scripts scripts/gpu_validate_dyn_cnbal_col_e2e.jl
# ==========================================================================

using CLM
using Printf
import Metal
include(joinpath(@__DIR__, "gpu_backends.jl"))

struct _F32 end
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:AbstractFloat}) = Float32.(x)
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:Integer}) = x
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{Bool}) = x
# NOTE: leave scalar Float64 fields (e.g. totvegcthresh) as-is — coercing them to
# Float32 breaks the Adapt positional reconstruction of the soil BGC state (its
# concrete ::Float64 field can't accept a Float32).

reldiff(a, b) = begin
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + max(abs(Float64(A[i])), abs(Float64(B[i])))))
    end
    m
end

const NLEV = (CLM.varpar_init!(CLM.varpar, 1, 15, 2, 5); CLM.varcon_init!(); CLM.varpar.nlevdecomp)
const NDP = 3

function build()
    CLM.varpar_init!(CLM.varpar, 1, 15, 2, 5); CLM.varcon_init!()
    ng, nl, nc, np = 1, 1, 2, 2
    grc = CLM.GridcellData(); lun = CLM.LandunitData()
    col = CLM.ColumnData();   pch = CLM.PatchData()
    CLM.gridcell_init!(grc, ng); CLM.landunit_init!(lun, nl)
    CLM.column_init!(col, nc);   CLM.patch_init!(pch, np)
    li = Ref(0); ci = Ref(0); pi = Ref(0)
    CLM.add_landunit!(lun, li, 1, CLM.ISTSOIL, 1.0)
    CLM.add_column!(col, lun, ci, li[], CLM.ISTSOIL, 0.5); CLM.add_patch!(pch, col, lun, pi, ci[], 1, 1.0)
    CLM.add_column!(col, lun, ci, li[], CLM.ISTSOIL, 0.5); CLM.add_patch!(pch, col, lun, pi, ci[], 1, 1.0)
    bounds = CLM.BoundsType(begg=1, endg=ng, begl=1, endl=nl, begc=1, endc=nc,
                            begp=1, endp=np, level=CLM.BOUNDS_LEVEL_PROC)
    CLM.clm_ptrs_compdown!(bounds, grc, lun, col, pch)

    csu = CLM.ColumnStateUpdater(bounds, 1)
    col.active .= [true, true]; col.wtgcell .= [0.5, 0.5]
    CLM.set_old_weights!(csu, bounds, grc, lun, col)
    col.wtgcell .= [0.3, 0.7]
    CLM.set_new_weights!(csu, bounds, 1, col)

    sc = CLM.SoilBiogeochemCarbonStateData{Float64}()
    CLM.soil_bgc_carbon_state_init!(sc, nc, ng, NLEV, NDP)
    sn = CLM.SoilBiogeochemNitrogenStateData{Float64}()
    CLM.soil_bgc_nitrogen_state_init!(sn, nc, ng, NLEV, NDP)
    for l in 1:NDP, j in 1:NLEV
        sc.decomp_cpools_vr_col[1, j, l] = 10.0 + j + l
        sc.decomp_cpools_vr_col[2, j, l] = 3.0  + j + l
        sn.decomp_npools_vr_col[1, j, l] = (10.0 + j + l) * 0.05
        sn.decomp_npools_vr_col[2, j, l] = (3.0  + j + l) * 0.05
    end
    for j in 1:NLEV
        sc.ctrunc_vr_col[1, j] = 0.5;  sc.ctrunc_vr_col[2, j] = 0.2
        sn.ntrunc_vr_col[1, j] = 0.01; sn.ntrunc_vr_col[2, j] = 0.005
        sn.sminn_vr_col[1, j]  = 1.0;  sn.sminn_vr_col[2, j]  = 0.4
    end
    return bounds, col, csu, sc, sn
end

function main(backend)
    println("=" ^ 72)
    println("WHOLE-FN Metal parity for dyn_cnbal_col! (dyn_cnbal soil-column C/N path)")
    println("=" ^ 72)
    if backend === nothing
        println("  No GPU backend — nothing to validate."); return 0
    end
    name, _, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n", name, FT)

    bH, colH, csuH, scH, snH = build()
    CLM.dyn_cnbal_col!(bH, 1, csuH, colH, scH, snH)

    mf(x) = CLM.Adapt.adapt(Metal.MtlArray, CLM.Adapt.adapt(_F32(), x))
    bD, colD, csuD, scD, snD = build()
    csuD = mf(csuD); colD = mf(colD); scD = mf(scD); snD = mf(snD)
    if !(scD.decomp_cpools_vr_col isa Metal.MtlArray)
        println("  BLOCKED: soil BGC state did not move to the device."); return 2
    end
    CLM.dyn_cnbal_col!(bD, 1, csuD, colD, scD, snD)

    nfail = 0
    checks = [
        ("sc.decomp_cpools_vr_col",  scH.decomp_cpools_vr_col,  scD.decomp_cpools_vr_col),
        ("sc.ctrunc_vr_col",         scH.ctrunc_vr_col,         scD.ctrunc_vr_col),
        ("sc.dyn_cbal_adjustments",  scH.dyn_cbal_adjustments_col, scD.dyn_cbal_adjustments_col),
        ("sn.decomp_npools_vr_col",  snH.decomp_npools_vr_col,  snD.decomp_npools_vr_col),
        ("sn.ntrunc_vr_col",         snH.ntrunc_vr_col,         snD.ntrunc_vr_col),
        ("sn.sminn_vr_col",          snH.sminn_vr_col,          snD.sminn_vr_col),
        ("sn.dyn_nbal_adjustments",  snH.dyn_nbal_adjustments_col, snD.dyn_nbal_adjustments_col),
    ]
    for (nm, a, b) in checks
        dd = reldiff(a, b); ok = dd < 1f-3
        @printf("  [%s] %-28s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, dd)
        ok || (nfail += 1)
    end
    println()
    println(nfail == 0 ? "  WHOLE dyn_cnbal_col! MATCHES CPU ON $name ($FT)" :
                         "  DIVERGENCE — investigate ($nfail failed).")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
