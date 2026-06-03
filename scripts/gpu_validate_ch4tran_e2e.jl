# ==========================================================================
# gpu_validate_ch4tran_e2e.jl — end-to-end Metal parity for the WHOLE
# ch4_tran! routine (the 2-species CH4/O2 reaction-diffusion vertical transport,
# the hardest methane function). One per-column kernel runs competition, the
# ebullition/aerenchyma reductions, the 2-species source/epsilon setup, and for
# s=1,2 the Patankar tridiagonal assembly + an in-thread Thomas solve.
#
# Builds a CH4Data + ch4_tran! arg set at Float64, runs ch4_tran! on the CPU,
# adapts the whole CH4Data (now @adapt_structure'd) + forcing to Float32/Metal,
# runs the SAME call on the device, and compares the mutated conc/flux outputs.
#
#   julia --project=scripts scripts/gpu_validate_ch4tran_e2e.jl
# ==========================================================================

using CLM
using Printf
import Metal
include(joinpath(@__DIR__, "gpu_backends.jl"))

function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + max(abs(Float64(A[i])), abs(Float64(B[i])))))
    end
    return m
end
allfinite(a) = all(isfinite, Array(a))

# Build a CH4Data + ch4_tran! inputs at precision FT. nc columns, nlevsoi levels.
function build(::Type{FT}; nc=3, nlevsoi=6, nlevsno=5) where {FT}
    ng = nc
    M(v) = fill(FT(v), nc, nlevsoi)
    V(v) = fill(FT(v), nc)
    ch4 = CLM.CH4Data{FT}(
        # selected (sat=0 / unsat) depth fluxes + concentrations + stresses
        ch4_prod_depth_unsat_col = M(2.0e-6), ch4_oxid_depth_unsat_col = M(1.0e-6),
        ch4_aere_depth_unsat_col = M(5.0e-7), ch4_ebul_depth_unsat_col = M(3.0e-7),
        o2_decomp_depth_unsat_col = M(2.0e-6), o2_oxid_depth_unsat_col = M(1.5e-6),
        o2_aere_depth_unsat_col = M(4.0e-6),
        conc_ch4_unsat_col = M(0.1), conc_o2_unsat_col = M(8.0),
        ch4stress_unsat_col = M(1.0), o2stress_unsat_col = M(1.0),
        co2_decomp_depth_unsat_col = M(1.0e-6),
        ch4_ebul_total_unsat_col = V(0.0), ch4_surf_aere_unsat_col = V(0.0),
        ch4_surf_ebul_unsat_col = V(0.0), ch4_surf_diff_unsat_col = V(0.0),
        # also the sat variants (so the sat=1 case has data)
        ch4_prod_depth_sat_col = M(2.0e-6), ch4_oxid_depth_sat_col = M(1.0e-6),
        ch4_aere_depth_sat_col = M(5.0e-7), ch4_ebul_depth_sat_col = M(3.0e-7),
        o2_decomp_depth_sat_col = M(2.0e-6), o2_oxid_depth_sat_col = M(1.5e-6),
        o2_aere_depth_sat_col = M(4.0e-6),
        conc_ch4_sat_col = M(0.1), conc_o2_sat_col = M(8.0),
        ch4stress_sat_col = M(1.0), o2stress_sat_col = M(1.0),
        co2_decomp_depth_sat_col = M(1.0e-6),
        ch4_ebul_total_sat_col = V(0.0), ch4_surf_aere_sat_col = V(0.0),
        ch4_surf_ebul_sat_col = V(0.0), ch4_surf_diff_sat_col = V(0.0),
        grnd_ch4_cond_col = V(1.0e-3),
        c_atm_grc = FT[1.7e-3 8.0 0.4; 1.7e-3 8.0 0.4; 1.7e-3 8.0 0.4][1:ng, :],
    )

    params = CLM.CH4Params()
    ch4vc  = CLM.CH4VarCon()

    return (; ch4, params, ch4vc,
        mask_soil = trues(nc),
        col_gridcell = collect(1:nc),
        watsat = M(0.5), h2osoi_vol = M(0.4), h2osoi_liq = M(300.0), h2osoi_ice = M(10.0),
        h2osfc = V(5.0), bsw = M(6.0), cellorg = M(20.0), t_soisno = fill(FT(CLM.TFRZ + 8.0), nc, nlevsoi),
        t_grnd = V(CLM.TFRZ + 9.0), t_h2osfc = V(CLM.TFRZ + 7.0), frac_h2osfc = V(0.2),
        snow_depth = V(0.0), snl = fill(0, nc),
        z = M(0.2), dz = M(0.1), zi = M(0.25),
        jwt = fill(3, nc),    # water table at level 3 (exercises jwt / jwt+1 branches)
        nlevsoi, nlevsno, organic_max = 130.0)
end

run_tran!(S, ch4; sat, lake) = CLM.ch4_tran!(
    ch4, S.params, S.ch4vc, S.mask_soil, S.col_gridcell,
    S.watsat, S.h2osoi_vol, S.h2osoi_liq, S.h2osoi_ice, S.h2osfc,
    S.bsw, S.cellorg, S.t_soisno, S.t_grnd, S.t_h2osfc, S.frac_h2osfc,
    S.snow_depth, S.snl, S.z, S.dz, S.zi, S.jwt, sat, lake,
    S.nlevsoi, S.nlevsno, 1800.0, S.organic_max)

function check_case(to, FT, label; sat, lake)
    @printf("\n  --- %s (sat=%d, lake=%s) ---\n", label, sat, lake)
    H = build(Float64)   # CPU reference (uniform Float64)
    B = build(Float64)   # device source

    ad(x) = to(x isa AbstractArray ? (eltype(x) <: AbstractFloat ? FT.(x) : x) : x)
    ch4_d = CLM.Adapt.adapt(Metal.MtlArray, CLM.Adapt.adapt(_F32(), B.ch4))
    Sd = (; B.params, B.ch4vc,
        mask_soil = B.mask_soil,   # ch4_tran! moves it to backend internally
        col_gridcell = B.col_gridcell,
        watsat=ad(B.watsat), h2osoi_vol=ad(B.h2osoi_vol), h2osoi_liq=ad(B.h2osoi_liq),
        h2osoi_ice=ad(B.h2osoi_ice), h2osfc=ad(B.h2osfc), bsw=ad(B.bsw), cellorg=ad(B.cellorg),
        t_soisno=ad(B.t_soisno), t_grnd=ad(B.t_grnd), t_h2osfc=ad(B.t_h2osfc),
        frac_h2osfc=ad(B.frac_h2osfc), snow_depth=ad(B.snow_depth), snl=B.snl,
        z=ad(B.z), dz=ad(B.dz), zi=ad(B.zi), jwt=B.jwt, B.nlevsoi, B.nlevsno, B.organic_max)

    if !(ch4_d.conc_ch4_sat_col isa Metal.MtlArray)
        println("  BLOCKED: CH4Data did not move to the device under adapt.")
        return 2
    end

    run_tran!(H, H.ch4; sat=sat, lake=lake)
    run_tran!(Sd, ch4_d; sat=sat, lake=lake)

    cc = sat == 0 ? (:conc_ch4_unsat_col, :conc_o2_unsat_col, :ch4_surf_diff_unsat_col) :
                    (:conc_ch4_sat_col, :conc_o2_sat_col, :ch4_surf_diff_sat_col)
    checks = [
        ("conc_ch4",      getfield(H.ch4, cc[1]),            getfield(ch4_d, cc[1])),
        ("conc_o2",       getfield(H.ch4, cc[2]),            getfield(ch4_d, cc[2])),
        ("ch4_surf_diff", getfield(H.ch4, cc[3]),            getfield(ch4_d, cc[3])),
        ("grnd_ch4_cond", H.ch4.grnd_ch4_cond_col,           ch4_d.grnd_ch4_cond_col),
    ]
    nfail = 0
    for (nm, a, b) in checks
        if !allfinite(a)
            @printf("  [WARN] %-14s CPU reference non-finite — skipping\n", nm); continue
        end
        d = reldiff(a, b); ok = d < 1f-3
        @printf("  [%s] %-14s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end
    return nfail
end

# Float32 down-convert adaptor (CH4Data fields are Float64 on the CPU build).
struct _F32 end
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:AbstractFloat}) = Float32.(x)
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:Integer}) = x
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{Bool}) = x

function main(backend)
    println("=" ^ 70)
    println("END-TO-END Metal parity for ch4_tran! (2-species reaction-diffusion)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n", name, FT)
    nfail = 0
    nfail += check_case(to, FT, "unsat"; sat=0, lake=false)
    nfail += check_case(to, FT, "sat";   sat=1, lake=false)
    println()
    println(nfail == 0 ? "  WHOLE ch4_tran! MATCHES CPU ON $name ($FT)" :
                         "  DIVERGENCE — investigate ($nfail failed).")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
