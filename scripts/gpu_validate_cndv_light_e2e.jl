# ==========================================================================
# gpu_validate_cndv_light_e2e.jl — end-to-end Metal parity for cndv_light!
# (CNDV light competition: per-patch crown/LAI/FPC + patch->gridcell scatter
# of tree/shrub/grass FPC totals + per-gridcell max + per-patch competition).
# eco/PftconType stay host (their fields are extracted into a device bundle);
# the DGVS state + patch go to the device.
#
#   julia --project=scripts scripts/gpu_validate_cndv_light_e2e.jl
# ==========================================================================
using CLM
using Printf
import Metal
include(joinpath(@__DIR__, "gpu_backends.jl"))

struct _F32 end
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:AbstractFloat}) = Float32.(x)
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:Integer}) = x
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{Bool}) = x

reldiff(a, b) = begin
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + max(abs(Float64(A[i])), abs(Float64(B[i])))))
    end
    m
end

function make_pft()
    pft = CLM.PftconType(); CLM.pftcon_allocate!(pft)
    for i in [2, 3, 4]
        pft.woody[i] = 1.0; pft.is_tree[i] = true; pft.is_shrub[i] = false
        pft.dwood[i] = 2.5e5; pft.slatop[i] = 0.012; pft.dsladlai[i] = 0.0
    end
    for i in [13, 14]
        pft.woody[i] = 0.0; pft.is_tree[i] = false; pft.is_shrub[i] = false
        pft.dwood[i] = 0.0; pft.slatop[i] = 0.030; pft.dsladlai[i] = 0.0
    end
    pft.pftpar20 .= 15.0; pft.pftpar28 .= -35.0; pft.pftpar29 .= 28.0
    pft.pftpar30 .= 350.0; pft.pftpar31 .= 1000.0
    return pft
end

function make_data()
    np = 4
    pft = make_pft()
    eco = CLM.DGVEcophysCon(); CLM.dgv_ecophyscon_init!(eco, pft)
    dgvs = CLM.DGVSData(); CLM.dgvs_init!(dgvs, np)
    dgvs.present_patch .= [false, true, true, true]
    dgvs.nind_patch .= [0.0, 0.05, 0.03, 1.0]
    dgvs.fpcgrid_patch .= [0.0, 0.6, 0.6, 0.3]   # tree total > 0.95 -> competition fires
    dgvs.fpcgridold_patch .= [0.0, 0.5, 0.5, 0.2]
    dgvs.crownarea_patch .= [0.0, 8.0, 8.0, 1.0]
    dgvs.agdd_patch .= [10.0, 20.0, 30.0, 40.0]
    dgvs.agddtw_patch .= [1.0, 2.0, 3.0, 4.0]
    dgvs.agdd20_patch .= [400.0, 500.0, 600.0, 700.0]
    dgvs.tmomin20_patch .= [CLM.TFRZ + 5.0, CLM.TFRZ + 8.0, CLM.TFRZ + 6.0, CLM.TFRZ + 3.0]
    dgvs.pftmayexist_patch .= [true, true, true, true]
    patch = CLM.PatchData(); CLM.patch_init!(patch, np)
    patch.itype .= [CLM.noveg, 1, 2, 12]; patch.gridcell .= 1; patch.column .= 1
    patch.landunit .= 1; patch.wtcol .= [0.1, 0.3, 0.3, 0.3]
    lun = CLM.LandunitData(); CLM.landunit_init!(lun, 1)
    lun.itype[1] = CLM.ISTSOIL; lun.wtgcell[1] = 1.0
    return (; dgvs, eco, pft, patch, lun, np, ng=1,
            deadstemc=[0.0, 5000.0, 3000.0, 0.0], leafcmax=[0.0, 200.0, 150.0, 50.0],
            t_a10=[280.0, 290.0, 300.0, 310.0], t_ref2m=[281.0, 291.0, 301.0, 311.0],
            prec365=[1.0e-5], annsum_npp=[0.0, 300.0, 250.0, 100.0],
            annsum_litfall=[0.0, 100.0, 80.0, 40.0], mask=trues(np))
end

run_light!(d) = CLM.cndv_light!(d.dgvs, d.eco, d.deadstemc, d.leafcmax,
    d.pft, d.patch, d.mask, 1:d.np; bounds_gridcell=1:d.ng)
run_init!(d) = CLM.dyn_cndv_init!(d.dgvs, d.patch, 1:d.np)
run_interp!(d) = CLM.dyn_cndv_interp!(d.dgvs, d.patch, d.lun, 1:d.np; wt1=0.3, is_beg_curr_year=true)
run_acc!(d) = CLM.cndv_update_acc_vars!(d.dgvs, d.eco, d.t_a10, d.t_ref2m, 1:d.np, 1800.0; month=6, day=15, secs=1800)
run_estab!(d) = CLM.cndv_establishment!(d.dgvs, d.eco, d.prec365, d.annsum_npp, d.annsum_litfall,
    d.deadstemc, d.leafcmax, d.pft, d.patch, d.lun, 1:d.np; bounds_gridcell=1:d.ng)

function main(backend)
    println("="^66); println("END-TO-END Metal parity for cndv_light!"); println("="^66)
    if backend === nothing; println("  No GPU backend."); return 0; end
    name, _, FT = backend
    @printf("  Backend: %s   (precision: %s)\n", name, FT)
    mf(x) = CLM.Adapt.adapt(Metal.MtlArray, CLM.Adapt.adapt(_F32(), x))
    movedev(B) = (; dgvs=mf(B.dgvs), eco=B.eco, pft=B.pft, patch=mf(B.patch), lun=mf(B.lun),
        np=B.np, ng=B.ng, deadstemc=mf(B.deadstemc), leafcmax=mf(B.leafcmax),
        t_a10=mf(B.t_a10), t_ref2m=mf(B.t_ref2m), prec365=mf(B.prec365),
        annsum_npp=mf(B.annsum_npp), annsum_litfall=mf(B.annsum_litfall),
        mask=Metal.MtlArray(collect(B.mask)))
    checks = Tuple{String,Any,Any}[]

    # cndv_light! (the scatter/competition kernel) — run on its own fresh data
    Hl = make_data(); Dl = movedev(make_data())
    if !(Dl.dgvs.fpcgrid_patch isa Metal.MtlArray); println("  BLOCKED."); return 2; end
    run_light!(Hl); run_light!(Dl)
    push!(checks, ("light: crownarea", Hl.dgvs.crownarea_patch, Dl.dgvs.crownarea_patch))
    push!(checks, ("light: fpcgrid", Hl.dgvs.fpcgrid_patch, Dl.dgvs.fpcgrid_patch))
    push!(checks, ("light: nind", Hl.dgvs.nind_patch, Dl.dgvs.nind_patch))

    # dyn_cndv_init!
    Hi = make_data(); Di = movedev(make_data()); run_init!(Hi); run_init!(Di)
    push!(checks, ("init: fpcgrid", Hi.dgvs.fpcgrid_patch, Di.dgvs.fpcgrid_patch))
    push!(checks, ("init: fpcgridold", Hi.dgvs.fpcgridold_patch, Di.dgvs.fpcgridold_patch))

    # dyn_cndv_interp!
    Hp = make_data(); Dp = movedev(make_data()); run_interp!(Hp); run_interp!(Dp)
    push!(checks, ("interp: wtcol", Hp.patch.wtcol, Dp.patch.wtcol))
    push!(checks, ("interp: fpcgridold", Hp.dgvs.fpcgridold_patch, Dp.dgvs.fpcgridold_patch))

    # cndv_update_acc_vars!
    Ha = make_data(); Da = movedev(make_data()); run_acc!(Ha); run_acc!(Da)
    push!(checks, ("acc: agddtw", Ha.dgvs.agddtw_patch, Da.dgvs.agddtw_patch))
    push!(checks, ("acc: agdd", Ha.dgvs.agdd_patch, Da.dgvs.agdd_patch))

    # cndv_establishment! (per-gridcell kernel)
    He = make_data(); De = movedev(make_data()); run_estab!(He); run_estab!(De)
    push!(checks, ("estab: present", He.dgvs.present_patch, De.dgvs.present_patch))
    push!(checks, ("estab: nind", He.dgvs.nind_patch, De.dgvs.nind_patch))
    push!(checks, ("estab: fpcgrid", He.dgvs.fpcgrid_patch, De.dgvs.fpcgrid_patch))
    push!(checks, ("estab: crownarea", He.dgvs.crownarea_patch, De.dgvs.crownarea_patch))
    push!(checks, ("estab: greffic", He.dgvs.greffic_patch, De.dgvs.greffic_patch))
    push!(checks, ("estab: heatstress", He.dgvs.heatstress_patch, De.dgvs.heatstress_patch))
    push!(checks, ("estab: leafcmax", He.leafcmax, De.leafcmax))
    nfail = 0
    for (nm, a, b) in checks
        dd = reldiff(a, b); ok = dd < 1f-3
        @printf("  [%s] %-18s rel = %.3e\n", ok ? "PASS" : "FAIL", nm, dd); ok || (nfail += 1)
    end
    println(); println(nfail == 0 ? "  cndv per-patch fns MATCH CPU ON $name ($FT)" : "  DIVERGENCE ($nfail).")
    return nfail == 0 ? 0 : 1
end
const BACKEND = detect_backend()
exit(main(BACKEND))
