# ==========================================================================
# gpu_validate_n_leaching_e2e.jl — end-to-end GPU parity for the WHOLE
# n_leaching! driver (mineral-N leaching + runoff): per-column level reductions
# (tot_water/surface_water via in-thread j-loops) + per-(c,j) leaching/runoff
# physics, depth-conditional, for both use_nitrif_denitrif modes.
#
#   julia --project=scripts scripts/gpu_validate_n_leaching_e2e.jl
# ==========================================================================

using CLM
using Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))

struct _F32 end
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:AbstractFloat}) = Float32.(x)
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:Integer}) = x
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{Bool}) = x

function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + max(abs(Float64(A[i])), abs(Float64(B[i])))))
    end
    return m
end
allfinite(a) = all(isfinite, Array(a))

function build(; nc=4, nlevdecomp=3, nlevsoi=3)
    ndecomp_pools=7; ndecomp_cascade_transitions=10
    params = CLM.NLeachingParams(sf=0.1, sf_no3=1.0)
    nf = CLM.SoilBiogeochemNitrogenFluxData()
    CLM.soil_bgc_nitrogen_flux_init!(nf, nc, nlevsoi, ndecomp_pools, ndecomp_cascade_transitions)
    ns = CLM.SoilBiogeochemNitrogenStateData()
    CLM.soil_bgc_nitrogen_state_init!(ns, nc, 1, nlevsoi, ndecomp_pools)
    for j in 1:nlevsoi, c in 1:nc
        ns.sminn_vr_col[c,j]=5.0; ns.smin_no3_vr_col[c,j]=2.0
    end
    col_dz = fill(0.1, nc, nlevsoi)
    zisoi = zeros(nlevsoi+1); for j in 1:nlevsoi; zisoi[j+1]=zisoi[j]+col_dz[1,j]; end
    return (; params, nf, ns, h2osoi_liq=fill(30.0,nc,nlevsoi),
        qflx_drain=fill(1.0e-5,nc), qflx_surf=fill(5.0e-6,nc), col_dz, zisoi, dt=1800.0,
        mask_bgc_soilc=trues(nc), nlevdecomp, nlevsoi, nc)
end

run_leach!(d, zisoi; use_nitrif) = CLM.n_leaching!(d.nf, d.ns, d.params;
    mask_bgc_soilc=d.mask_bgc_soilc, bounds=1:d.nc, nlevdecomp=d.nlevdecomp, nlevsoi=d.nlevsoi,
    dt=d.dt, h2osoi_liq=d.h2osoi_liq, qflx_drain=d.qflx_drain, qflx_surf=d.qflx_surf,
    col_dz=d.col_dz, zisoi=zisoi, use_nitrif_denitrif=use_nitrif)

function check_case(name, FT, use_nitrif)
    @printf("\n  --- use_nitrif_denitrif=%s ---\n", use_nitrif)
    H = build(); B = build()
    run_leach!(H, H.zisoi; use_nitrif=use_nitrif)
    mf(x) = CLM.Adapt.adapt(device_array_type(), CLM.Adapt.adapt(_F32(), x))
    D = (; params=B.params, nf=mf(B.nf), ns=mf(B.ns), h2osoi_liq=mf(B.h2osoi_liq),
        qflx_drain=mf(B.qflx_drain), qflx_surf=mf(B.qflx_surf), col_dz=mf(B.col_dz),
        dt=B.dt, mask_bgc_soilc=device_array_type()(collect(B.mask_bgc_soilc)),
        nlevdecomp=B.nlevdecomp, nlevsoi=B.nlevsoi, nc=B.nc)
    zisoi_d = device_array_type()(Float32.(B.zisoi))
    if !(D.nf.sminn_leached_vr_col isa device_array_type()); println("  BLOCKED."); return 2; end
    run_leach!(D, zisoi_d; use_nitrif=use_nitrif)
    outs = use_nitrif ? (:smin_no3_leached_vr_col, :smin_no3_runoff_vr_col) : (:sminn_leached_vr_col,)
    nfail = 0
    for f in outs
        a = getfield(H.nf, f); b = getfield(D.nf, f)
        allfinite(a) || (@printf("  [WARN] %-26s non-finite\n", f); continue)
        dd = reldiff(a, b); ok = dd < 1f-3
        @printf("  [%s] %-26s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", f, dd)
        ok || (nfail += 1)
    end
    return nfail
end

function main(backend)
    println("=" ^ 72); println("END-TO-END GPU parity for n_leaching!"); println("=" ^ 72)
    if backend === nothing; println("  No GPU backend."); return 0; end
    name, _, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n", name, FT)
    nfail = check_case(name, FT, false) + check_case(name, FT, true)
    println()
    println(nfail == 0 ? "  WHOLE n_leaching! MATCHES CPU ON $name ($FT)" :
                         "  DIVERGENCE — investigate ($nfail failed).")
    return nfail == 0 ? 0 : 1
end
const BACKEND = detect_backend()
exit(main(BACKEND))
