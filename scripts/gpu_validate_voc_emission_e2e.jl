# ==========================================================================
# gpu_validate_voc_emission_e2e.jl — end-to-end Metal parity for voc_emission!
# (MEGAN VOC emissions). The per-patch kernel calls the 6 hardened get_gamma_*
# helpers + get_map_EF and reads the flattened compound metadata + the
# device-movable MEGANFactors. The Vector{MEGANCompound/MechComp} metadata stays
# on the host (the function flattens it to numeric arrays internally).
#
#   julia --project=scripts scripts/gpu_validate_voc_emission_e2e.jl
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

function make_data()
    np = 4; nc = 3; ng = 2
    iso = CLM.MEGANCompound(name="isoprene", index=1, class_number=1,
        molec_weight=68.12, coeff=1.0, emis_factors=fill(600.0, 20))
    mono = CLM.MEGANCompound(name="myrcene", index=2, class_number=2,
        molec_weight=136.23, coeff=1.0, emis_factors=fill(100.0, 20))
    meg_compounds = [iso, mono]
    mech_comps = [CLM.MEGANMechComp(name="ISOP", n_megan_comps=1, megan_indices=[1]),
                  CLM.MEGANMechComp(name="TERP", n_megan_comps=1, megan_indices=[2])]
    mfac = CLM.MEGANFactors(); CLM.megan_factors_init!(mfac, 20)
    voc = CLM.VOCEmisData(); CLM.vocemis_init!(voc, np, ng, 2, 2)
    voc.efisop_grc[:, 1] .= 600.0; voc.efisop_grc[:, 2] .= 500.0
    patch = CLM.PatchData(); CLM.patch_init!(patch, np)
    patch.itype .= [0, 1, 6, 12]; patch.gridcell .= [1, 1, 1, 2]; patch.column .= [1, 1, 2, 3]
    ci_val = 0.7 * 400.0e-6 * 101325.0
    return (; voc, meg_compounds, mech_comps, mfac, patch, bounds_p=1:np,
            mask=trues(np),
            forc_solad_col=fill(200.0, nc, 2), forc_solai_grc=fill(100.0, ng, 2),
            forc_pbot_col=fill(101325.0, nc), forc_pco2_grc=fill(40.53, ng),
            fsd24=fill(200.0, np), fsd240=fill(200.0, np), fsi24=fill(100.0, np), fsi240=fill(100.0, np),
            fsun=fill(0.5, np), fsun24=fill(0.5, np), fsun240=fill(0.5, np),
            elai=fill(2.0, np), elai240=fill(2.0, np),
            cisun_z=fill(ci_val, np, 1), cisha_z=fill(ci_val, np, 1),
            t_veg=fill(300.0, np), t_veg24=fill(300.0, np), t_veg240=fill(297.0, np),
            btran=fill(0.8, np))
end

run_voc!(d) = CLM.voc_emission!(d.voc, d.meg_compounds, d.mech_comps, d.mfac,
    d.patch, d.bounds_p, d.mask,
    d.forc_solad_col, d.forc_solai_grc, d.forc_pbot_col, d.forc_pco2_grc,
    d.fsd24, d.fsd240, d.fsi24, d.fsi240,
    d.fsun, d.fsun24, d.fsun240, d.elai, d.elai240,
    d.cisun_z, d.cisha_z, d.t_veg, d.t_veg24, d.t_veg240, d.btran)

function main(backend)
    println("="^66); println("END-TO-END Metal parity for voc_emission! (MEGAN)"); println("="^66)
    if backend === nothing; println("  No GPU backend."); return 0; end
    name, _, FT = backend
    @printf("  Backend: %s   (precision: %s)\n", name, FT)
    H = make_data(); B = make_data()
    run_voc!(H)
    mf(x) = CLM.Adapt.adapt(Metal.MtlArray, CLM.Adapt.adapt(_F32(), x))
    D = (; voc=mf(B.voc), meg_compounds=B.meg_compounds, mech_comps=B.mech_comps,
         mfac=mf(B.mfac), patch=mf(B.patch), bounds_p=B.bounds_p,
         mask=Metal.MtlArray(collect(B.mask)),
         forc_solad_col=mf(B.forc_solad_col), forc_solai_grc=mf(B.forc_solai_grc),
         forc_pbot_col=mf(B.forc_pbot_col), forc_pco2_grc=mf(B.forc_pco2_grc),
         fsd24=mf(B.fsd24), fsd240=mf(B.fsd240), fsi24=mf(B.fsi24), fsi240=mf(B.fsi240),
         fsun=mf(B.fsun), fsun24=mf(B.fsun24), fsun240=mf(B.fsun240),
         elai=mf(B.elai), elai240=mf(B.elai240),
         cisun_z=mf(B.cisun_z), cisha_z=mf(B.cisha_z),
         t_veg=mf(B.t_veg), t_veg24=mf(B.t_veg24), t_veg240=mf(B.t_veg240), btran=mf(B.btran))
    if !(D.voc.vocflx_tot_patch isa Metal.MtlArray); println("  BLOCKED."); return 2; end
    run_voc!(D)
    checks = [("vocflx_tot_patch", H.voc.vocflx_tot_patch, D.voc.vocflx_tot_patch),
              ("vocflx_patch", H.voc.vocflx_patch, D.voc.vocflx_patch),
              ("meg_flux_out", H.voc.meg_flux_out, D.voc.meg_flux_out),
              ("gamma_out", H.voc.gamma_out_patch, D.voc.gamma_out_patch),
              ("gammaP_out", H.voc.gammaP_out_patch, D.voc.gammaP_out_patch),
              ("gammaT_out", H.voc.gammaT_out_patch, D.voc.gammaT_out_patch),
              ("gammaA_out", H.voc.gammaA_out_patch, D.voc.gammaA_out_patch),
              ("gammaC_out", H.voc.gammaC_out_patch, D.voc.gammaC_out_patch),
              ("paru_out", H.voc.paru_out_patch, D.voc.paru_out_patch),
              ("Eopt_out", H.voc.Eopt_out_patch, D.voc.Eopt_out_patch),
              ("topt_out", H.voc.topt_out_patch, D.voc.topt_out_patch)]
    nfail = 0
    for (nm, a, b) in checks
        dd = reldiff(a, b); ok = dd < 1f-3
        @printf("  [%s] %-18s rel = %.3e\n", ok ? "PASS" : "FAIL", nm, dd); ok || (nfail += 1)
    end
    println(); println(nfail == 0 ? "  voc_emission! MATCHES CPU ON $name ($FT)" : "  DIVERGENCE ($nfail).")
    return nfail == 0 ? 0 : 1
end
const BACKEND = detect_backend()
exit(main(BACKEND))
