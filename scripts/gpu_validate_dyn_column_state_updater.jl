# ==========================================================================
# gpu_validate_dyn_column_state_updater.jl — GPU parity for the dyn_subgrid
# COLUMN-state conservative updater (the update path: gridcell reduction-scatter
# of shrinking-column mass, per-gridcell finalize, distribute-to-growing-columns).
#
# Mirrors test_dyn_column_state_updater's 3-column single-gridcell subgrid (2
# natveg + 1 special lake). set_old/new_weights! are host per-step setup; the
# update variants (no_special_handling / fill_special_using_natveg /
# fill_using_fixed_values) run on the device after adapting the updater + subgrid
# + state to Metal/Float32, and are compared to the CPU reference.
#
#   julia --project=scripts scripts/gpu_validate_dyn_column_state_updater.jl
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

function build_subgrid()
    CLM.varpar_init!(CLM.varpar, 1, 15, 2, 5)
    CLM.varcon_init!()
    ng, nl, nc, np = 1, 2, 3, 3
    grc = CLM.GridcellData(); lun = CLM.LandunitData()
    col = CLM.ColumnData();   pch = CLM.PatchData()
    CLM.gridcell_init!(grc, ng); CLM.landunit_init!(lun, nl)
    CLM.column_init!(col, nc);   CLM.patch_init!(pch, np)
    li = Ref(0); ci = Ref(0); pi = Ref(0)
    CLM.add_landunit!(lun, li, 1, CLM.ISTSOIL, 0.6)
    CLM.add_column!(col, lun, ci, li[], CLM.ISTSOIL, 0.5); CLM.add_patch!(pch, col, lun, pi, ci[], 1, 1.0)
    CLM.add_column!(col, lun, ci, li[], CLM.ISTSOIL, 0.5); CLM.add_patch!(pch, col, lun, pi, ci[], 1, 1.0)
    CLM.add_landunit!(lun, li, 1, CLM.ISTDLAK, 0.4)
    CLM.add_column!(col, lun, ci, li[], CLM.ISTDLAK, 1.0); CLM.add_patch!(pch, col, lun, pi, ci[], 1, 1.0)
    bounds = CLM.BoundsType(begg=1, endg=ng, begl=1, endl=nl, begc=1, endc=nc,
                            begp=1, endp=np, level=CLM.BOUNDS_LEVEL_PROC)
    CLM.clm_ptrs_compdown!(bounds, grc, lun, col, pch)
    return bounds, grc, lun, col, pch
end

# Prepare an updater with old=[0.3,0.3,0.4], new=[c1,c2,c3] (host set_old/new).
function prep_updater(bounds, grc, lun, col, new_w)
    csu = CLM.ColumnStateUpdater(bounds, 1)
    col.active  .= [true, true, true]
    col.wtgcell .= [0.3, 0.3, 0.4]
    CLM.set_old_weights!(csu, bounds, grc, lun, col)
    col.wtgcell .= new_w
    CLM.set_new_weights!(csu, bounds, 1, col)
    return csu
end

function main(backend)
    println("=" ^ 72)
    println("GPU parity for the dyn_subgrid COLUMN-state updater (conservation core)")
    println("=" ^ 72)
    if backend === nothing
        println("  No GPU backend — nothing to validate."); return 0
    end
    name, _, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n", name, FT)

    mf(x) = CLM.Adapt.adapt(device_array_type(), CLM.Adapt.adapt(_F32(), x))
    nfail = 0
    report(nm, a, b) = begin
        dd = reldiff(a, b); ok = dd < 1f-3
        @printf("  [%s] %-40s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, dd)
        ok || (global nfail += 1)
    end

    # --- Variant 1: no_special_handling (c1 shrinks, c2 grows) ---
    b, grc, lun, col, _ = build_subgrid()
    csuH = prep_updater(b, grc, lun, col, [0.1, 0.5, 0.4])
    varH = [10.0, 4.0, 7.0]; adjH = zeros(3)
    CLM.update_column_state_no_special_handling!(csuH, b, 1, col, varH; adjustment=adjH)

    bD, grcD, lunD, colD, _ = build_subgrid()
    csuD = mf(prep_updater(bD, grcD, lunD, colD, [0.1, 0.5, 0.4]))
    colDm = mf(colD)
    varD = mf([10.0, 4.0, 7.0]); adjD = mf(zeros(3))
    if !(csuD.area_gained_col isa device_array_type())
        println("  BLOCKED: ColumnStateUpdater did not move to the device."); return 2
    end
    CLM.update_column_state_no_special_handling!(csuD, b, 1, colDm, varD; adjustment=adjD)
    report("variant1 var (conservative redistribute)", varH, varD)
    report("variant1 adjustment", adjH, adjD)

    # --- Variant 2: fill_special_using_natveg (lake c3 shrinks, natveg c2 grows) ---
    b2, grc2, lun2, col2, _ = build_subgrid()
    csuH2 = prep_updater(b2, grc2, lun2, col2, [0.3, 0.5, 0.2])
    varH2 = [10.0, 4.0, 999.0]; ncH2 = zeros(1); adjH2 = zeros(3)
    CLM.update_column_state_fill_special_using_natveg!(csuH2, b2, 1, grc2, lun2, col2,
        varH2, ncH2; adjustment=adjH2)

    b2D, grc2D, lun2D, col2D, _ = build_subgrid()
    csuD2 = mf(prep_updater(b2D, grc2D, lun2D, col2D, [0.3, 0.5, 0.2]))
    varD2 = mf([10.0, 4.0, 999.0]); ncD2 = mf(zeros(1)); adjD2 = mf(zeros(3))
    CLM.update_column_state_fill_special_using_natveg!(csuD2, b2D, 1, mf(grc2D), mf(lun2D),
        mf(col2D), varD2, ncD2; adjustment=adjD2)
    report("variant2 var (natveg-seeded special)", varH2, varD2)
    report("variant2 non_conserved_mass", ncH2, ncD2)

    # --- Variant 3: fill_using_fixed_values (lake c3 shrinks -> fixed 5.0; c1 grows) ---
    b3, grc3, lun3, col3, _ = build_subgrid()
    csuH3 = prep_updater(b3, grc3, lun3, col3, [0.5, 0.3, 0.2])
    lv = fill(CLM.FILLVAL_USE_EXISTING_VALUE, CLM.MAX_LUNIT); lv[CLM.ISTDLAK] = 5.0
    varH3 = [10.0, 4.0, 999.0]; ncH3 = zeros(1)
    CLM.update_column_state_fill_using_fixed_values!(csuH3, b3, 1, lun3, col3, varH3, lv, ncH3)

    b3D, grc3D, lun3D, col3D, _ = build_subgrid()
    csuD3 = mf(prep_updater(b3D, grc3D, lun3D, col3D, [0.5, 0.3, 0.2]))
    varD3 = mf([10.0, 4.0, 999.0]); ncD3 = mf(zeros(1))
    CLM.update_column_state_fill_using_fixed_values!(csuD3, b3D, 1, mf(lun3D), mf(col3D),
        varD3, mf(lv), ncD3)
    report("variant3 var (fixed-value special)", varH3, varD3)
    report("variant3 non_conserved_mass", ncH3, ncD3)

    println()
    println(nfail == 0 ? "  COLUMN-state updater conservation core MATCHES CPU ON $name ($FT)" :
                         "  DIVERGENCE — investigate ($nfail failed).")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
