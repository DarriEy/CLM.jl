# ==========================================================================
# gpu_validate_c_iso_flux_e2e.jl — end-to-end Metal parity for the carbon-
# isotope flux module (CNCIsoFluxMod).  Exercises the two grouped patch->
# column scatter wrappers that needed device-view struct grouping to clear
# Metal's ~31-arg limit:
#   (1) cn_c_iso_gap_pft_to_column!     — gap-mortality scatter (9 bundles)
#   (2) cn_c_iso_harvest_pft_to_column! — harvest scatter (6 bundles) + the
#       loose wood-harvest 1D scatter launch
#
#   julia --project=scripts scripts/gpu_validate_c_iso_flux_e2e.jl
# ==========================================================================
using CLM
using Printf
import Metal
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

# ---- gap patch->column scatter setup (mirrors test_c_iso_flux) ---------
function build_gap()
    np = 4; nc = 2; nlevdecomp = 3; i_litr_min = 1; i_litr_max = 2; i_met_lit = 1
    cf = CLM.CNVegCarbonFluxData()
    flds = (:m_leafc_to_litter_patch, :m_frootc_to_litter_patch, :m_livestemc_to_litter_patch,
            :m_deadstemc_to_litter_patch, :m_livecrootc_to_litter_patch, :m_deadcrootc_to_litter_patch,
            :m_leafc_storage_to_litter_patch, :m_frootc_storage_to_litter_patch,
            :m_livestemc_storage_to_litter_patch, :m_deadstemc_storage_to_litter_patch,
            :m_livecrootc_storage_to_litter_patch, :m_deadcrootc_storage_to_litter_patch,
            :m_gresp_storage_to_litter_patch, :m_leafc_xfer_to_litter_patch,
            :m_frootc_xfer_to_litter_patch, :m_livestemc_xfer_to_litter_patch,
            :m_deadstemc_xfer_to_litter_patch, :m_livecrootc_xfer_to_litter_patch,
            :m_deadcrootc_xfer_to_litter_patch, :m_gresp_xfer_to_litter_patch)
    for (k, f) in enumerate(flds)
        setfield!(cf, f, collect(range(0.01 * k, 0.01 * k + 0.03, length=np)))
    end
    cf.gap_mortality_c_to_litr_c_col = zeros(nc, nlevdecomp, i_litr_max)
    cf.gap_mortality_c_to_cwdc_col = zeros(nc, nlevdecomp)
    sb = CLM.SoilBiogeochemStateData()
    sb.leaf_prof_patch = fill(0.4, np, nlevdecomp)
    sb.froot_prof_patch = fill(0.3, np, nlevdecomp)
    sb.croot_prof_patch = fill(0.2, np, nlevdecomp)
    sb.stem_prof_patch = fill(0.1, np, nlevdecomp)
    return (; cf, sb, mask=fill(true, np), bounds=1:np, nlevdecomp, i_litr_min, i_litr_max, i_met_lit,
            patch_column=[1, 1, 2, 2], patch_itype=[1, 1, 1, 1], patch_wtcol=fill(0.5, np),
            lf_f=ones(20, i_litr_max), fr_f=ones(20, i_litr_max), np)
end

# ---- harvest patch->column scatter setup -------------------------------
function build_harvest()
    np = 4; nc = 2; nlevdecomp = 3; i_litr_min = 1; i_litr_max = 2; i_met_lit = 1
    cf = CLM.CNVegCarbonFluxData()
    hflds = (:hrv_leafc_to_litter_patch, :hrv_frootc_to_litter_patch, :hrv_livestemc_to_litter_patch,
             :hrv_livecrootc_to_litter_patch, :hrv_deadcrootc_to_litter_patch,
             :hrv_leafc_storage_to_litter_patch, :hrv_frootc_storage_to_litter_patch,
             :hrv_livestemc_storage_to_litter_patch, :hrv_deadstemc_storage_to_litter_patch,
             :hrv_livecrootc_storage_to_litter_patch, :hrv_deadcrootc_storage_to_litter_patch,
             :hrv_gresp_storage_to_litter_patch, :hrv_leafc_xfer_to_litter_patch,
             :hrv_frootc_xfer_to_litter_patch, :hrv_livestemc_xfer_to_litter_patch,
             :hrv_deadstemc_xfer_to_litter_patch, :hrv_livecrootc_xfer_to_litter_patch,
             :hrv_deadcrootc_xfer_to_litter_patch, :hrv_gresp_xfer_to_litter_patch,
             :wood_harvestc_patch)
    for (k, f) in enumerate(hflds)
        setfield!(cf, f, collect(range(0.01 * k, 0.01 * k + 0.03, length=np)))
    end
    cf.harvest_c_to_litr_c_col = zeros(nc, nlevdecomp, i_litr_max)
    cf.harvest_c_to_cwdc_col = zeros(nc, nlevdecomp)
    cf.wood_harvestc_col = zeros(nc)
    sb = CLM.SoilBiogeochemStateData()
    sb.leaf_prof_patch = fill(0.4, np, nlevdecomp)
    sb.froot_prof_patch = fill(0.3, np, nlevdecomp)
    sb.croot_prof_patch = fill(0.2, np, nlevdecomp)
    sb.stem_prof_patch = fill(0.1, np, nlevdecomp)
    return (; cf, sb, mask=fill(true, np), bounds=1:np, nlevdecomp, i_litr_min, i_litr_max, i_met_lit,
            patch_column=[1, 1, 2, 2], patch_itype=[1, 1, 1, 1], patch_wtcol=fill(0.5, np),
            lf_f=ones(20, i_litr_max), fr_f=ones(20, i_litr_max), np)
end

run_gap!(d) = CLM.cn_c_iso_gap_pft_to_column!(d.cf, d.sb, d.mask, d.bounds, d.nlevdecomp,
                                              patch_column=d.patch_column, patch_itype=d.patch_itype,
                                              patch_wtcol=d.patch_wtcol, lf_f=d.lf_f, fr_f=d.fr_f,
                                              i_litr_min=d.i_litr_min, i_litr_max=d.i_litr_max,
                                              i_met_lit=d.i_met_lit)
run_harvest!(d) = CLM.cn_c_iso_harvest_pft_to_column!(d.cf, d.sb, d.mask, d.bounds, d.nlevdecomp,
                                              patch_column=d.patch_column, patch_itype=d.patch_itype,
                                              patch_wtcol=d.patch_wtcol, lf_f=d.lf_f, fr_f=d.fr_f,
                                              i_litr_min=d.i_litr_min, i_litr_max=d.i_litr_max,
                                              i_met_lit=d.i_met_lit)

function main(backend)
    println("=" ^ 66); println("END-TO-END Metal parity for CNCIsoFluxMod"); println("=" ^ 66)
    if backend === nothing; println("  No GPU backend."); return 0; end
    name, _, FT = backend
    @printf("  Backend: %s   (precision: %s)\n", name, FT)
    mf(x) = CLM.Adapt.adapt(Metal.MtlArray, CLM.Adapt.adapt(_F32(), x))
    mmask(b) = Metal.MtlArray(collect(b))
    checks = Tuple{String,Any,Any}[]

    movedev(d) = (; cf=mf(d.cf), sb=mf(d.sb), mask=mmask(d.mask), bounds=d.bounds, nlevdecomp=d.nlevdecomp,
                  i_litr_min=d.i_litr_min, i_litr_max=d.i_litr_max, i_met_lit=d.i_met_lit,
                  patch_column=mf(d.patch_column), patch_itype=mf(d.patch_itype), patch_wtcol=mf(d.patch_wtcol),
                  lf_f=mf(d.lf_f), fr_f=mf(d.fr_f), np=d.np)

    # (1) gap-mortality patch->column scatter (grouped device-view kernel)
    Hg = build_gap(); Bg = build_gap(); run_gap!(Hg); Dg = movedev(Bg)
    if !(Dg.cf.gap_mortality_c_to_litr_c_col isa Metal.MtlArray); println("  BLOCKED (gap)."); return 2; end
    run_gap!(Dg)
    push!(checks, ("gap_mortality_c_to_litr_c_col", Hg.cf.gap_mortality_c_to_litr_c_col, Dg.cf.gap_mortality_c_to_litr_c_col))
    push!(checks, ("gap_mortality_c_to_cwdc_col", Hg.cf.gap_mortality_c_to_cwdc_col, Dg.cf.gap_mortality_c_to_cwdc_col))

    # (2) harvest patch->column scatter (grouped kernel) + wood-harvest 1D launch
    Hh = build_harvest(); Bh = build_harvest(); run_harvest!(Hh); Dh = movedev(Bh)
    run_harvest!(Dh)
    push!(checks, ("harvest_c_to_litr_c_col", Hh.cf.harvest_c_to_litr_c_col, Dh.cf.harvest_c_to_litr_c_col))
    push!(checks, ("harvest_c_to_cwdc_col", Hh.cf.harvest_c_to_cwdc_col, Dh.cf.harvest_c_to_cwdc_col))
    push!(checks, ("wood_harvestc_col", Hh.cf.wood_harvestc_col, Dh.cf.wood_harvestc_col))

    nfail = 0
    for (nm, a, b) in checks
        dd = reldiff(a, b); ok = dd < 1f-3
        @printf("  [%s] %-34s rel = %.3e\n", ok ? "PASS" : "FAIL", nm, dd); ok || (nfail += 1)
    end
    println(); println(nfail == 0 ? "  CNCIsoFluxMod kernels MATCH CPU ON $name ($FT)" : "  DIVERGENCE ($nfail).")
    return nfail == 0 ? 0 : 1
end
const BACKEND = detect_backend()
exit(main(BACKEND))
