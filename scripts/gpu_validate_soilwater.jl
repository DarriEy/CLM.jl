# ==========================================================================
# gpu_validate_soilwater.jl — real-hardware GPU parity for the soil_water!
# (Zeng-Decker 2009) kernels. Same spirit as gpu_validate_soiltemp.jl: run each
# @kernel on a CPU Array and on a device array of the SAME inputs at the SAME
# precision and compare. Parity ("does the device execute the kernel the way the
# CPU does?"), not accuracy-vs-Fortran.
#
#   julia --project=scripts scripts/gpu_validate_soilwater.jl
#
# Kernels covered (swm-zd1..zd4 increments of soilwater_zengdecker2009!):
#   _soilwm_voleq_zq_kernel!      : equilibrium water content vol_eq + matric pot. zq
#   _soilwm_hk_smp_kernel!        : hydraulic conductivity hk, matric pot. smp, derivs
#   _soilwm_jwt_vwczwt_kernel!    : water-table layer index jwt + vwc at WT depth
#   _soilwm_voleq_zq_aqu_kernel!  : aquifer-layer vol_eq/zq (water table below column)
#   _soilwm_aqu_zmm_kernel!       : aquifer-layer geometry zmm/dzmm
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
const NLEVSNO = 5
const NLEVGRND = 15          # > NLEVSOI so t_soisno[joff+nlevsoi+1] is in bounds
const JOFF = NLEVSNO         # offset for z/dz/t_soisno/h2osoi (snow+soil)
const JOFF_ZI = NLEVSNO + 1  # offset for zi (extra top element)

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

    # interface depths in mm, strictly increasing in j (index 1 = Fortran zimm(c,0)=0):
    # zimm_arr[c, j+1] = 100*j, so zimm_arr[c, NLEVSOI+1] = 1000.
    zimm_arr = FT[100.0 * (j - 1) for c in 1:NC, j in 1:(NLEVSOI + 1)]
    # Water-table depth per column (mm), CONSISTENT with jwt_in below (jwt is
    # derived from zwt in the real driver). Aquifer columns (jwt==NLEVSOI) must
    # have zwtmm > 1000 so the aquifer kernel's (suc+zwt-zimm_j) base stays > 0;
    # in-column columns sit between their zimm interfaces.
    #             jwt:  10     3    10      7   10    0     5    10
    zwtmm = FT[1050.0, 350.0, 1075.0, 750.0, 1100.0, 50.0, 550.0, 1125.0]

    watsat = FT(0.35) .+ rnd(NC, NLEVSOI) .* FT(0.1)
    sucsat = FT(100.0) .+ rnd(NC, NLEVSOI) .* FT(100.0)   # mm
    bsw    = FT(4.0)   .+ rnd(NC, NLEVSOI) .* FT(3.0)     # > 1 so (1-1/bsw) in (0,1)
    smpmin = FT(-1.0e8) .* (FT(1.0) .+ rnd(NC))
    hksat  = FT(1.0e-5) .+ rnd(NC, NLEVSOI) .* FT(1.0e-5)
    icefrac = rnd(NC, NLEVSOI) .* FT(0.5)
    # vwc_liq has nlevsoi+1 columns (aquifer slot); kernels read j and min(nlevsoi,j+1)
    vwc_liq = FT(0.05) .+ rnd(NC, NLEVSOI + 1) .* FT(0.25)

    # --- extra inputs for the per-column search/aquifer kernels (swm-zd2) ---
    # interface depths zi in m, strictly increasing in k; zi[c, JOFF_ZI+j] = 0.1*(NLEVSNO+j)
    zi = FT[0.1 * (k - 1) for c in 1:NC, k in 1:(NLEVSNO + 1 + NLEVGRND)]
    # water-table depth (m) spread so some columns sit within the column, some below
    zwt = FT[0.4 + 0.18 * (c - 1) for c in 1:NC]
    # t_soisno (snow+ground) spanning TFRZ so the frozen-water-table branch fires
    t_soisno = FT[273.15 + 1.2 * sinpi((c + k) / 6) for c in 1:NC, k in 1:(NLEVSNO + NLEVGRND)]
    h2osoi_vol = FT(0.2) .+ rnd(NC, NLEVSOI) .* FT(0.15)
    # jwt mix: several == NLEVSOI (water table below column → aquifer branch active)
    jwt_in = Int[NLEVSOI, 3, NLEVSOI, 7, NLEVSOI, 0, 5, NLEVSOI]
    # zmm/dzmm pre-filled (aquifer kernel reads [c,nlevsoi], writes [c,nlevsoi+1]);
    # zmm strictly increasing in j so den = zmm[j+1]-zmm[j] = 100 > 0 (no /0).
    zmm_in  = FT[100.0 * j for c in 1:NC, j in 1:(NLEVSOI + 1)]
    dzmm_in = fill(FT(100.0), NC, NLEVSOI + 1)

    # --- extra inputs for the tridiagonal-assembly kernel (swm-zd3) ---
    zq_in     = FT[-(50.0 + 5.0 * j) for c in 1:NC, j in 1:(NLEVSOI + 1)]
    smp_in    = FT[-(200.0 + 7.0 * j) - 3.0 * c for c in 1:NC, j in 1:(NLEVSOI + 1)]
    hk_in     = FT(1.0e-6) .+ rnd(NC, NLEVSOI) .* FT(1.0e-6)
    dhkdw_in  = FT(1.0e-7) .+ rnd(NC, NLEVSOI) .* FT(1.0e-7)
    dsmpdw_in = FT(-5.0e3) .- rnd(NC, NLEVSOI) .* FT(1.0e3)
    qflx_rootsoi_in = rnd(NC, NLEVSOI) .* FT(1.0e-4)
    qflx_infl_in    = rnd(NC) .* FT(1.0e-3)
    vwc_zwt_in      = FT(0.2) .+ rnd(NC) .* FT(0.15)

    # --- extra inputs for the post-solve kernels (swm-zd4) ---
    # h2osoi_liq over snow+soil (JOFF+NLEVSOI cols); soil layers carry a few
    # slightly negative values so the deficit smooth_ifelse branch fires.
    h2osoi_liq_in = FT[ (k > JOFF && (c + k) % 7 == 0) ? -0.5 : 2.0 + 0.3 * k
                        for c in 1:NC, k in 1:(JOFF + NLEVSOI) ]
    z_in     = FT[0.1 * k for c in 1:NC, k in 1:(JOFF + NLEVSOI)]   # increasing depth (m)
    dwat2_in = (rnd(NC, NLEVSOI + 1) .- FT(0.5)) .* FT(1.0e-3)
    imped_in = FT(0.3) .+ rnd(NC, NLEVSOI) .* FT(0.7)              # ice impedance in (0,1]

    return (; mask, zimm_arr, zwtmm, watsat, sucsat, bsw, smpmin,
            hksat, icefrac, vwc_liq, e_ice = FT(6.0),
            zi, zwt, t_soisno, h2osoi_vol, jwt_in, zmm_in, dzmm_in,
            zq_in, smp_in, hk_in, dhkdw_in, dsmpdw_in, qflx_rootsoi_in,
            qflx_infl_in, vwc_zwt_in, sdamp = FT(0.0), dtime = FT(1800.0),
            h2osoi_liq_in, z_in, dwat2_in, imped_in)
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

    # ---- jwt / vwc_zwt (out = vwc_zwt, extra = jwt[Int]) ----
    vwc_zwt0 = zeros(FT, NC)
    push!(results, _parity("jwt_vwczwt", to, FT, CLM._soilwm_jwt_vwczwt_kernel!,
        NC, vwc_zwt0,
        f -> begin
            jwt = f(zeros(Int, NC))
            ((jwt, f(I.mask), f(I.zwt), f(I.zi), f(I.watsat), f(I.vwc_liq),
              f(I.t_soisno), f(I.sucsat), f(I.bsw), f(I.h2osoi_vol),
              JOFF, JOFF_ZI, NLEVSOI, NLEVGRND), (jwt,))
        end))

    # ---- aquifer vol_eq/zq (out = vol_eq, extra = zq); needs jwt input ----
    voleqA0 = zeros(FT, NC, NLEVSOI + 1)
    push!(results, _parity("voleq_zq_aqu", to, FT, CLM._soilwm_voleq_zq_aqu_kernel!,
        NC, voleqA0,
        f -> begin
            zq = f(zeros(FT, NC, NLEVSOI + 1))
            ((zq, f(I.mask), f(I.jwt_in), f(I.zwtmm), f(I.zimm_arr), f(I.watsat),
              f(I.sucsat), f(I.bsw), f(I.smpmin), NLEVSOI), (zq,))
        end))

    # ---- aquifer zmm/dzmm (out = zmm pre-filled, extra = dzmm pre-filled) ----
    push!(results, _parity("aqu_zmm", to, FT, CLM._soilwm_aqu_zmm_kernel!,
        NC, I.zmm_in,
        f -> begin
            dzmm = f(copy(I.dzmm_in))
            ((dzmm, f(I.mask), f(I.jwt_in), f(I.zwt), NLEVSOI), (dzmm,))
        end))

    # ---- tridiagonal assembly (device-view structs SwmTriOut/SwmTriIn) ----
    push!(results, test_tridiag(to, FT, I))

    # ---- post-solve: renew liquid + qcharge (out = h2osoi_liq, extra = qcharge) ----
    push!(results, _parity("renew_qcharge", to, FT, CLM._soilwm_renew_qcharge_kernel!,
        NC, I.h2osoi_liq_in,
        f -> begin
            qcharge = f(zeros(FT, NC))
            ((qcharge, f(I.mask), f(I.jwt_in), f(I.dwat2_in), f(I.dzmm_in),
              f(I.h2osoi_vol), f(I.watsat), f(I.imped_in), f(I.hksat), f(I.bsw),
              f(I.smpmin), f(I.smp_in), f(I.zq_in), f(I.zwt), f(I.z_in),
              JOFF, NLEVSOI, FT(I.dtime)), (qcharge,))
        end))

    # ---- post-solve: water deficit reduction (out = qflx_deficit) ----
    push!(results, _parity("deficit", to, FT, CLM._soilwm_deficit_kernel!,
        NC, zeros(FT, NC),
        f -> ((f(I.mask), f(I.h2osoi_liq_in), JOFF, NLEVSOI), ())))

    return results
end

# tridiagonal-assembly kernel: build the SwmTriOut/SwmTriIn structs on each
# backend, run, and compare every output array. Mirrors test_pc_beta in
# gpu_validate_soiltemp.jl (struct-arg kernel → manual backend + synchronize).
function test_tridiag(to, ::Type{FT}, I) where {FT}
    z2() = zeros(FT, NC, NLEVSOI + 1)
    mkout(cv) = CLM.SwmTriOut(; qin = cv(z2()), qout = cv(z2()), dqidw0 = cv(z2()),
        dqidw1 = cv(z2()), dqodw1 = cv(z2()), dqodw2 = cv(z2()), rmx = cv(z2()),
        amx = cv(z2()), bmx = cv(z2()), cmx = cv(z2()))
    mkin(cv) = CLM.SwmTriIn(; zmm = cv(copy(I.zmm_in)), dzmm = cv(copy(I.dzmm_in)),
        zq = cv(copy(I.zq_in)), smp = cv(copy(I.smp_in)), hk = cv(copy(I.hk_in)),
        dhkdw = cv(copy(I.dhkdw_in)), dsmpdw = cv(copy(I.dsmpdw_in)),
        qflx_rootsoi = cv(copy(I.qflx_rootsoi_in)), watsat = cv(copy(I.watsat)),
        sucsat = cv(copy(I.sucsat)), bsw = cv(copy(I.bsw)), vwc_liq = cv(copy(I.vwc_liq)),
        qflx_infl = cv(copy(I.qflx_infl_in)), vwc_zwt = cv(copy(I.vwc_zwt_in)),
        smpmin = cv(copy(I.smpmin)))

    run!(out, tin, to_) = CLM.soilwm_tridiag_assemble!(out, tin,
        _dev(to_, I.mask), _dev(to_, I.jwt_in), I.sdamp, I.dtime, NLEVSOI)

    oc = mkout(identity); ic = mkin(identity); run!(oc, ic, identity)
    cv = x -> _dev(to, x)
    od = mkout(cv); id = mkin(cv); run!(od, id, to)

    d = maximum(reldiff(getfield(oc, f), getfield(od, f)) for f in fieldnames(CLM.SwmTriOut))
    return ("tridiag_assemble", d, d < TOL(FT))
end

function main(backend)
    println("=" ^ 70)
    println("METAL PARITY for soilwater_zengdecker2009! kernels (swm-zd1..zd4)")
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
    println(nfail == 0 ? "\n  ALL ZD09 swm-zd1..zd4 kernels MATCH CPU ON $name ($FT) ✓" :
                          "\n  SOME KERNELS DIVERGED — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
