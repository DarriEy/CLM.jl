# ==========================================================================
# gpu_validate_soilwater.jl — real-hardware GPU parity for the soil_water!
# (Zeng-Decker 2009) kernels. Same spirit as gpu_validate_soiltemp.jl: run each
# @kernel on a CPU Array and on a device array of the SAME inputs at the SAME
# precision and compare. Parity ("does the device execute the kernel the way the
# CPU does?"), not accuracy-vs-Fortran.
#
#   julia --project=scripts scripts/gpu_validate_soilwater.jl
#
# Kernels covered (swm-zd1, increment 1 of soilwater_zengdecker2009!):
#   _soilwm_voleq_zq_kernel!  : equilibrium water content vol_eq + matric pot. zq
#   _soilwm_hk_smp_kernel!    : hydraulic conductivity hk, matric pot. smp, derivs
# ==========================================================================

using CLM
using Printf
using Random

include(joinpath(@__DIR__, "gpu_backends.jl"))

maxabsdiff(a, b) = maximum(abs.(Array(a) .- Array(b)); init = 0.0)
# Mixed abs/rel error: |a-b| / (1 + max|a|,|b|). Relative for large-magnitude
# quantities (matric potential smp/zq clamp near smpmin ~ -1e8; dsmpdw ~ 1e10),
# absolute for small ones. The right metric here: CPU-Float32 vs Metal-Float32
# differ only in pow()'s last bits, so a 1e8 quantity shows ~1e2 absolute but
# ~1e-6 relative round-off. Same spirit as the relaxed phase_change_beta bound.
function reldiff(a, b)
    A = Array(a); B = Array(b)
    maximum(abs.(A .- B) ./ (one(eltype(A)) .+ max.(abs.(A), abs.(B))); init = 0.0)
end
# Float32 relative bound generous enough for power/division-heavy math (s^(-bsw)).
TOL(::Type{Float32}) = 1f-2
TOL(::Type{T}) where {T} = 1e-9

const NC = 8
const NLEVSOI = 10

# Same _parity contract as gpu_validate_soiltemp.jl: out0 is the primary output;
# mk(f) returns (positional kernel args after `out`, tuple of extra written arrays
# to also diff). f maps each array to its cpu (identity) or device form.
function _parity(name, to, FT, K, ndr, out0, mk; tol = TOL(FT))
    out_cpu = copy(out0)
    args_cpu, extras_cpu = mk(identity)
    CLM._launch!(K, out_cpu, args_cpu...; ndrange = ndr)
    out_dev = to(copy(out0))
    args_dev, extras_dev = mk(x -> _dev(to, x))
    CLM._launch!(K, out_dev, args_dev...; ndrange = ndr)
    d = max(reldiff(out_cpu, out_dev),
            maximum((reldiff(a, b) for (a, b) in zip(extras_cpu, extras_dev)); init = 0.0))
    return (name, d, d < tol)
end

# Branch-exercising inputs at precision FT. zwtmm is spread across the per-column
# interface depths zimm_arr so the three vol_eq branches (above / within / below
# the water table) all fire across the 8 columns.
function build_inputs(::Type{FT}) where {FT}
    rng = MersenneTwister(19)
    rnd(d...) = FT.(rand(rng, d...))
    mask = trues(NC)

    # interface depths in mm, strictly increasing in j (index 1 = Fortran zimm(c,0)=0)
    zimm_arr = FT[100.0 * (j - 1) for c in 1:NC, j in 1:(NLEVSOI + 1)]
    # one water-table depth per column, spanning below-top .. below-bottom
    zwtmm = FT[50.0 + 110.0 * (c - 1) for c in 1:NC]

    watsat = FT(0.35) .+ rnd(NC, NLEVSOI) .* FT(0.1)
    sucsat = FT(100.0) .+ rnd(NC, NLEVSOI) .* FT(100.0)   # mm
    bsw    = FT(4.0)   .+ rnd(NC, NLEVSOI) .* FT(3.0)     # > 1 so (1-1/bsw) in (0,1)
    smpmin = FT(-1.0e8) .* (FT(1.0) .+ rnd(NC))
    hksat  = FT(1.0e-5) .+ rnd(NC, NLEVSOI) .* FT(1.0e-5)
    icefrac = rnd(NC, NLEVSOI) .* FT(0.5)
    # vwc_liq has nlevsoi+1 columns (aquifer slot); kernels read j and min(nlevsoi,j+1)
    vwc_liq = FT(0.05) .+ rnd(NC, NLEVSOI + 1) .* FT(0.25)

    return (; mask, zimm_arr, zwtmm, watsat, sucsat, bsw, smpmin,
            hksat, icefrac, vwc_liq, e_ice = FT(6.0))
end

function tests(to, ::Type{FT}) where {FT}
    I = build_inputs(FT)
    results = Tuple{String,Float64,Bool}[]

    # ---- vol_eq / zq (out = vol_eq, extra = zq) ----
    voleq0 = zeros(FT, NC, NLEVSOI + 1)
    push!(results, _parity("voleq_zq", to, FT, CLM._soilwm_voleq_zq_kernel!,
        (NC, NLEVSOI), voleq0,
        f -> begin
            zq = f(zeros(FT, NC, NLEVSOI + 1))
            ((zq, f(I.mask), f(I.zwtmm), f(I.zimm_arr), f(I.watsat), f(I.sucsat),
              f(I.bsw), f(I.smpmin)), (zq,))
        end))

    # ---- hk / smp (out = hk, extras = the other six written arrays) ----
    hk0 = zeros(FT, NC, NLEVSOI)
    push!(results, _parity("hk_smp", to, FT, CLM._soilwm_hk_smp_kernel!,
        (NC, NLEVSOI), hk0,
        f -> begin
            dhkdw   = f(zeros(FT, NC, NLEVSOI))
            imped   = f(zeros(FT, NC, NLEVSOI))
            smp     = f(zeros(FT, NC, NLEVSOI))
            dsmpdw  = f(zeros(FT, NC, NLEVSOI))
            smp_l   = f(zeros(FT, NC, NLEVSOI))
            hk_l    = f(zeros(FT, NC, NLEVSOI))
            ((dhkdw, imped, smp, dsmpdw, smp_l, hk_l, f(I.mask), f(I.vwc_liq),
              f(I.watsat), f(I.hksat), f(I.bsw), f(I.icefrac), f(I.sucsat),
              f(I.smpmin), NLEVSOI, FT(I.e_ice)),
             (dhkdw, imped, smp, dsmpdw, smp_l, hk_l))
        end))

    return results
end

function main(backend)
    println("=" ^ 70)
    println("METAL PARITY for soilwater_zengdecker2009! kernels (swm-zd1)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend detected — nothing to validate (kernels run on KA CPU in the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Detected backend: %s   (working precision: %s)\n\n", name, FT)
    res = tests(to, FT)
    npass = nfail = 0
    for (knm, d, ok) in res
        @printf("  [%s] %-12s  rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", knm, d)
        ok ? (npass += 1) : (nfail += 1)
    end
    @printf("\n  %d passed, %d failed\n", npass, nfail)
    println(nfail == 0 ? "\n  ALL ZD09 swm-zd1 kernels MATCH CPU ON $name ($FT) ✓" :
                          "\n  SOME KERNELS DIVERGED — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
