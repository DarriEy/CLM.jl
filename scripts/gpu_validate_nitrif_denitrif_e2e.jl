# ==========================================================================
# gpu_validate_nitrif_denitrif_e2e.jl — end-to-end Metal parity for the WHOLE
# nitrif_denitrif! driver (nitrification + denitrification: soil diffusivity,
# anaerobic fraction, CENTURY nitrification, del Grosso denitrification, N2:N2O).
# One per-column kernel with an internal level loop; exp/atan/sqrt/pow-heavy.
#
# Builds a nitrif/denitrif dataset (clone of test_nitrif_denitrif's
# make_nitrif_denitrif_data) at Float64, runs nitrif_denitrif! on the CPU, adapts
# nf/ns/cf + the soil/water/temp col arrays to Metal/Float32, runs the SAME call on
# the device, and compares the mutated per-(col,level) flux fields.
#
#   julia --project=scripts scripts/gpu_validate_nitrif_denitrif_e2e.jl
# ==========================================================================

using CLM
using Printf
import Metal
include(joinpath(@__DIR__, "gpu_backends.jl"))

struct _F32 end
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:AbstractFloat}) = Float32.(x)
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:Integer}) = x
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{Bool}) = x
# NB: do NOT convert scalar Float64 fields — the Soil/CNVeg N structs have concrete
# ::Float64 scalar fields the parametric ctor pins as Float64.

function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + max(abs(Float64(A[i])), abs(Float64(B[i])))))
    end
    return m
end
allfinite(a) = all(isfinite, Array(a))

function build(; nc=4, nlevdecomp=1)
    ndecomp_pools = 7; ndecomp_cascade_transitions = 10
    params = CLM.NitrifDenitrifParams(
        k_nitr_max_perday=0.1, surface_tension_water=0.073, rij_kro_a=1.5e-10,
        rij_kro_alpha=1.26, rij_kro_beta=0.6, rij_kro_gamma=0.6, rij_kro_delta=0.85,
        denitrif_respiration_coefficient=0.1, denitrif_respiration_exponent=1.3,
        denitrif_nitrateconc_coefficient=0.1, denitrif_nitrateconc_exponent=1.3, om_frac_sf=1.0)
    cn_params = CLM.CNSharedParamsData(organic_max=130.0)
    nf = CLM.SoilBiogeochemNitrogenFluxData()
    CLM.soil_bgc_nitrogen_flux_init!(nf, nc, nlevdecomp, ndecomp_pools, ndecomp_cascade_transitions)
    ns = CLM.SoilBiogeochemNitrogenStateData()
    CLM.soil_bgc_nitrogen_state_init!(ns, nc, 1, nlevdecomp, ndecomp_pools)
    for j in 1:nlevdecomp, c in 1:nc
        ns.smin_nh4_vr_col[c,j] = 1.0; ns.smin_no3_vr_col[c,j] = 2.0
    end
    cf = CLM.SoilBiogeochemCarbonFluxData()
    CLM.soil_bgc_carbon_flux_init!(cf, nc, nlevdecomp, ndecomp_pools, ndecomp_cascade_transitions)
    for j in 1:nlevdecomp, c in 1:nc
        cf.t_scalar_col[c,j]=0.8; cf.w_scalar_col[c,j]=0.6; cf.phr_vr_col[c,j]=1.0e-6
    end
    M(v) = fill(Float64(v), nc, nlevdecomp)
    return (; params, cn_params, nf, ns, cf,
        watsat=M(0.45), watfc=M(0.30), bd=M(1300.0), bsw=M(5.0), cellorg=M(20.0),
        sucsat=M(200.0), soilpsi=M(-1.0), h2osoi_vol=M(0.30), h2osoi_liq=M(30.0),
        t_soisno=fill(CLM.TFRZ+15.0, nc, nlevdecomp), col_dz=M(0.1),
        o2_decomp_depth_unsat=M(1.0e-6), conc_o2_unsat=M(0.2),
        mask_bgc_soilc=trues(nc), nlevdecomp, nc)
end

run_nd!(d; use_lch4) = CLM.nitrif_denitrif!(d.nf, d.ns, d.cf, d.params, d.cn_params;
    mask_bgc_soilc=d.mask_bgc_soilc, bounds=1:d.nc, nlevdecomp=d.nlevdecomp,
    watsat=d.watsat, watfc=d.watfc, bd=d.bd, bsw=d.bsw, cellorg=d.cellorg,
    sucsat=d.sucsat, soilpsi=d.soilpsi, h2osoi_vol=d.h2osoi_vol, h2osoi_liq=d.h2osoi_liq,
    t_soisno=d.t_soisno, col_dz=d.col_dz,
    o2_decomp_depth_unsat=d.o2_decomp_depth_unsat, conc_o2_unsat=d.conc_o2_unsat,
    use_lch4=use_lch4, no_frozen_nitrif_denitrif=false)

const OUT = (:diffus_col, :anaerobic_frac_col, :k_nitr_t_vr_col, :k_nitr_ph_vr_col,
             :k_nitr_h2o_vr_col, :k_nitr_vr_col, :pot_f_nit_vr_col, :soil_bulkdensity_col,
             :smin_no3_massdens_vr_col, :pot_f_denit_vr_col, :n2_n2o_ratio_denit_vr_col,
             :r_psi_col, :fmax_denit_nitrate_vr_col)

function check_case(name, FT, use_lch4)
    @printf("\n  --- use_lch4=%s ---\n", use_lch4)
    H = build(); B = build()
    run_nd!(H; use_lch4=use_lch4)
    mf(x) = CLM.Adapt.adapt(Metal.MtlArray, CLM.Adapt.adapt(_F32(), x))
    D = (; params=B.params, cn_params=B.cn_params, nf=mf(B.nf), ns=mf(B.ns), cf=mf(B.cf),
        watsat=mf(B.watsat), watfc=mf(B.watfc), bd=mf(B.bd), bsw=mf(B.bsw), cellorg=mf(B.cellorg),
        sucsat=mf(B.sucsat), soilpsi=mf(B.soilpsi), h2osoi_vol=mf(B.h2osoi_vol), h2osoi_liq=mf(B.h2osoi_liq),
        t_soisno=mf(B.t_soisno), col_dz=mf(B.col_dz),
        o2_decomp_depth_unsat=mf(B.o2_decomp_depth_unsat), conc_o2_unsat=mf(B.conc_o2_unsat),
        mask_bgc_soilc=Metal.MtlArray(collect(B.mask_bgc_soilc)), nlevdecomp=B.nlevdecomp, nc=B.nc)
    if !(D.nf.diffus_col isa Metal.MtlArray); println("  BLOCKED: nf not on device."); return 2; end
    run_nd!(D; use_lch4=use_lch4)
    nfail = 0
    for f in OUT
        a = getfield(H.nf, f); b = getfield(D.nf, f)
        allfinite(a) || (@printf("  [WARN] %-28s CPU non-finite\n", f); continue)
        dd = reldiff(a, b); ok = dd < 1f-3
        @printf("  [%s] %-28s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", f, dd)
        ok || (nfail += 1)
    end
    return nfail
end

function main(backend)
    println("=" ^ 72); println("END-TO-END Metal parity for nitrif_denitrif!"); println("=" ^ 72)
    if backend === nothing; println("  No GPU backend."); return 0; end
    name, _, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n", name, FT)
    nfail = check_case(name, FT, false) + check_case(name, FT, true)
    println()
    println(nfail == 0 ? "  WHOLE nitrif_denitrif! MATCHES CPU ON $name ($FT)" :
                         "  DIVERGENCE — investigate ($nfail failed).")
    return nfail == 0 ? 0 : 1
end
const BACKEND = detect_backend()
exit(main(BACKEND))
