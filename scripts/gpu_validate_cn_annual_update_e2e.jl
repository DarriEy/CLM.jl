# ==========================================================================
# gpu_validate_cn_annual_update_e2e.jl — end-to-end GPU parity for
# cn_annual_update! (end-of-year column flag kernel + patch annual-rollover
# kernel + the kernelized p2c_1d_filter! averaging). Triggers end-of-year so
# the p2c averaging path actually runs on the device.
#
#   julia --project=scripts scripts/gpu_validate_cn_annual_update_e2e.jl
# ==========================================================================
using CLM
using Printf
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

function make_data()
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    nc = 2; np = 4
    days_per_year = 365.0; secspday = CLM.SECSPDAY
    secspyear = days_per_year * secspday
    patch_column = [1, 1, 2, 2]
    col = CLM.ColumnData(); CLM.column_init!(col, nc)
    col.is_fates[1:nc] .= false
    patch = CLM.PatchData(); CLM.patch_init!(patch, np)
    patch.column[1:np] .= patch_column
    patch.active[1:np] .= true
    for p in 1:np; patch.wtcol[p] = 1.0 / count(==(patch_column[p]), patch_column); end
    for c in 1:nc
        plist = findall(==(c), patch_column)
        col.patchi[c] = first(plist); col.patchf[c] = last(plist)
    end
    vs = CLM.CNVegStateData(); CLM.cnveg_state_init!(vs, np, nc)
    vs.annsum_counter_col[1:nc] .= secspyear          # force end-of-year
    vs.tempsum_potential_gpp_patch[1:np] .= [10.0, 20.0, 30.0, 40.0]
    vs.annsum_potential_gpp_patch[1:np] .= 0.0
    vs.tempmax_retransn_patch[1:np] .= [1.0, 2.0, 3.0, 4.0]
    vs.annmax_retransn_patch[1:np] .= 0.0
    vs.tempavg_t2m_patch[1:np] .= [280.0, 281.0, 282.0, 283.0]
    vs.annavg_t2m_patch[1:np] .= 0.0; vs.annavg_t2m_col[1:nc] .= 0.0
    cf = CLM.CNVegCarbonFluxData(); CLM.cnveg_carbon_flux_init!(cf, np, nc, 1)
    cf.tempsum_npp_patch[1:np] .= [100.0, 200.0, 300.0, 400.0]
    cf.annsum_npp_patch[1:np] .= 0.0
    cf.tempsum_litfall_patch[1:np] .= [5.0, 6.0, 7.0, 8.0]
    cf.annsum_litfall_patch[1:np] .= 0.0; cf.annsum_npp_col[1:nc] .= 0.0
    return (; col, patch, vs, cf, mask_c=trues(nc), mask_p=trues(np),
            bounds_c=1:nc, bounds_p=1:np, dt=1800.0, days_per_year, secspday)
end

run_au!(d) = CLM.cn_annual_update!(d.mask_c, d.mask_p, d.bounds_c, d.bounds_p,
    d.col, d.patch, d.vs, d.cf; dt=d.dt, days_per_year=d.days_per_year, secspday=d.secspday)

function main(backend)
    println("="^66); println("END-TO-END GPU parity for cn_annual_update!"); println("="^66)
    if backend === nothing; println("  No GPU backend."); return 0; end
    name, _, FT = backend
    @printf("  Backend: %s   (precision: %s)\n", name, FT)
    H = make_data(); B = make_data()
    run_au!(H)
    mf(x) = CLM.Adapt.adapt(device_array_type(), CLM.Adapt.adapt(_F32(), x))
    mm(b) = device_array_type()(collect(b))
    D = (; col=mf(B.col), patch=mf(B.patch), vs=mf(B.vs), cf=mf(B.cf),
         mask_c=mm(B.mask_c), mask_p=mm(B.mask_p), bounds_c=B.bounds_c, bounds_p=B.bounds_p,
         dt=B.dt, days_per_year=B.days_per_year, secspday=B.secspday)
    if !(D.vs.annsum_potential_gpp_patch isa device_array_type()); println("  BLOCKED."); return 2; end
    run_au!(D)
    checks = [("annsum_potential_gpp", H.vs.annsum_potential_gpp_patch, D.vs.annsum_potential_gpp_patch),
              ("tempsum_potential_gpp", H.vs.tempsum_potential_gpp_patch, D.vs.tempsum_potential_gpp_patch),
              ("annmax_retransn", H.vs.annmax_retransn_patch, D.vs.annmax_retransn_patch),
              ("tempmax_retransn", H.vs.tempmax_retransn_patch, D.vs.tempmax_retransn_patch),
              ("annavg_t2m_patch", H.vs.annavg_t2m_patch, D.vs.annavg_t2m_patch),
              ("annsum_npp_patch", H.cf.annsum_npp_patch, D.cf.annsum_npp_patch),
              ("annsum_litfall_patch", H.cf.annsum_litfall_patch, D.cf.annsum_litfall_patch),
              ("annsum_npp_col (p2c)", H.cf.annsum_npp_col, D.cf.annsum_npp_col),
              ("annavg_t2m_col (p2c)", H.vs.annavg_t2m_col, D.vs.annavg_t2m_col),
              ("annsum_counter_col", H.vs.annsum_counter_col, D.vs.annsum_counter_col)]
    nfail = 0
    for (nm, a, b) in checks
        dd = reldiff(a, b); ok = dd < 1f-3
        @printf("  [%s] %-22s rel = %.3e\n", ok ? "PASS" : "FAIL", nm, dd); ok || (nfail += 1)
    end
    println(); println(nfail == 0 ? "  cn_annual_update! MATCHES CPU ON $name ($FT)" : "  DIVERGENCE ($nfail).")
    return nfail == 0 ? 0 : 1
end
const BACKEND = detect_backend()
exit(main(BACKEND))
