# ==========================================================================
# gpu_ad_validate.jl — real-hardware AUTODIFF validation for CLM.jl's kernels.
#
# Phase 5: verify that automatic differentiation works THROUGH the Phase-4 kernels
# when they run on the GPU. Forward-mode (ForwardDiff `Dual` numbers) flows through
# a kernel as a custom isbits element type: an array of `Dual`s on the device runs
# the same kernel and the derivative rides along in the partial component.
#
# For each kernel we differentiate the output w.r.t. a scalar θ added to one input
# (so every output carries a partial), then check:
#   1. PARITY  — device-AD partials vs CPU-AD partials (the core question: does AD
#                give the SAME answer on the GPU as on the CPU?). PASS/FAIL on this.
#   2. CORRECT — AD partials vs a finite-difference estimate (reported; confirms the
#                derivative itself is right, not just that CPU and GPU agree).
#
# RUN
#   julia --project=scripts scripts/gpu_ad_validate.jl
#
# REVERSE MODE (Enzyme): forward-mode is what is validated here. Reverse-mode
# through these kernels works on the CPU backend, but Enzyme does not differentiate
# the Metal device kernel launch (no Enzyme GPU support for the Metal/Apple backend);
# see the note printed at the end. On CUDA, Enzyme-over-KA is the intended reverse path.
# ==========================================================================

using CLM
using Printf
using Random
using ForwardDiff
import KernelAbstractions  # synchronize() after the phase_change_beta struct-kernel launch

include(joinpath(@__DIR__, "gpu_backends.jl"))

# --- forward-mode helpers ---------------------------------------------------
_dualT(::Type{FT}) where {FT} = ForwardDiff.Dual{Nothing,FT,1}
_seed(v, p) = ForwardDiff.Dual(v, p)                 # value v, single partial p
_todual(::Type{FT}, A) where {FT} = _dualT(FT).(A)   # lift to Dual, partial 0
_part(A) = ForwardDiff.partials.(Array(A), 1)        # extract the 1st partial
_val(A)  = ForwardDiff.value.(Array(A))

maxabs(x) = maximum(abs.(x); init = 0.0)
relerr(a, b) = maxabs(a .- b) / max(maxabs(b), eps(Float64))

# PASS thresholds. Parity (CPU-AD vs device-AD) must be tight; the finite-difference
# correctness check is generous (Float32 FD is noisy, esp. through the solver sweeps).
_parity_tol(::Type{Float32}) = 1f-3
_parity_tol(::Type{Float64}) = 1e-9
const FD_REL_TOL = 0.05    # 5% relative — finite-difference sanity, not a tight bound

# Each test returns (name, parity, fd_rel, ok).
function adtest_forc_q(to, ::Type{FT}) where {FT}
    rng = MersenneTwister(1); nc = 32
    gc = collect(1:nc)
    vp  = FT.(500 .+ 600 .* rand(rng, nc))
    pb  = FT.(80000 .+ 10000 .* rand(rng, nc))
    # θ added to every vp → each out carries d(out)/d(vp)
    vpd = _todual(FT, vp); for i in eachindex(vpd); vpd[i] = _seed(vp[i], one(FT)); end
    pbd = _todual(FT, pb)
    oc = zeros(_dualT(FT), nc);    CLM.compute_forc_q!(oc, gc, vpd, pbd)
    od = to(zeros(_dualT(FT), nc)); CLM.compute_forc_q!(od, _dev(to,gc), _dev(to,vpd), _dev(to,pbd))
    parity = maxabs(_part(oc) .- _part(od))
    # finite difference wrt θ: perturb all vp
    h = FT(1f-1)
    o0 = zeros(FT, nc); CLM.compute_forc_q!(o0, gc, vp, pb)
    o1 = zeros(FT, nc); CLM.compute_forc_q!(o1, gc, vp .+ h, pb)
    fd = (o1 .- o0) ./ h
    fdrel = relerr(_part(od), fd)
    return ("compute_forc_q!", parity, fdrel, parity < _parity_tol(FT) && fdrel < FD_REL_TOL)
end

function adtest_tridiagonal(to, ::Type{FT}) where {FT}
    rng = MersenneTwister(2); nc, nl = 24, 16
    a = FT.(-0.5 .* ones(nc, nl)); c = FT.(-0.5 .* ones(nc, nl))
    b = FT.(2 .+ 0.5 .* rand(rng, nc, nl)); r = FT.(rand(rng, nc, nl))
    a[:,1] .= 0; c[:,nl] .= 0; jt = ones(Int, nc); mk = trues(nc)
    rd = _todual(FT, r); for i in eachindex(rd); rd[i] = _seed(r[i], one(FT)); end
    ad = _todual(FT, a); bd = _todual(FT, b); cd = _todual(FT, c)
    uc = zeros(_dualT(FT), nc, nl); CLM.tridiagonal_multi!(uc, ad, bd, cd, rd, jt, mk, nc, nl)
    ud = to(zeros(_dualT(FT), nc, nl))
    CLM.tridiagonal_multi!(ud, _dev(to,ad), _dev(to,bd), _dev(to,cd), _dev(to,rd), _dev(to,jt), _dev(to,mk), nc, nl)
    parity = maxabs(_part(uc) .- _part(ud))
    h = FT(1f-3)
    u0 = zeros(FT, nc, nl); CLM.tridiagonal_multi!(u0, a, b, c, r,        jt, mk, nc, nl)
    u1 = zeros(FT, nc, nl); CLM.tridiagonal_multi!(u1, a, b, c, r .+ h,   jt, mk, nc, nl)
    fd = (u1 .- u0) ./ h
    fdrel = relerr(_part(ud), fd)
    return ("tridiagonal_multi! (batched Thomas)", parity, fdrel, parity < _parity_tol(FT) && fdrel < FD_REL_TOL)
end

function adtest_band(to, ::Type{FT}) where {FT}
    rng = MersenneTwister(3); nc = 24; kl = ku = 2; nband = 5; nlevsno = 5; nlev = 25
    jt = ones(Int, nc); jb = fill(15, nc); mk = trues(nc)
    bm = zeros(FT, nc, nband, nlev)
    for c in 1:nc, lev in 1:nlev
        bm[c, kl+1, lev] = FT(4) + rand(rng, FT)
        bm[c, kl, lev] = FT(-0.7); bm[c, kl+2, lev] = FT(-0.7)
        bm[c, 1, lev] = FT(-0.3);  bm[c, nband, lev] = FT(-0.3)
    end
    rv = FT.(rand(rng, nc, nlev))
    rvd = _todual(FT, rv); for i in eachindex(rvd); rvd[i] = _seed(rv[i], one(FT)); end
    bmd = _todual(FT, bm)
    tc = zeros(_dualT(FT), nc, nlev); CLM.batched_band_solve!(tc, bmd, rvd, jt, jb, mk, kl, ku, nlevsno)
    td = to(zeros(_dualT(FT), nc, nlev))
    CLM.batched_band_solve!(td, _dev(to,bmd), _dev(to,rvd), _dev(to,jt), _dev(to,jb), _dev(to,mk), kl, ku, nlevsno)
    parity = maxabs(_part(tc) .- _part(td))
    h = FT(1f-3)
    t0 = zeros(FT, nc, nlev); CLM.batched_band_solve!(t0, bm, rv,      jt, jb, mk, kl, ku, nlevsno)
    t1 = zeros(FT, nc, nlev); CLM.batched_band_solve!(t1, bm, rv .+ h, jt, jb, mk, kl, ku, nlevsno)
    fd = (t1 .- t0) ./ h
    fdrel = relerr(_part(td), fd)
    return ("batched_band_solve! (pentadiagonal GE)", parity, fdrel, parity < _parity_tol(FT) && fdrel < FD_REL_TOL)
end

# ==========================================================================
# soil_temperature! kernels (Metal pipeline fn 4-8) — forward-mode AD parity.
# Same idea as above: seed ONE float input with partial 1, lift the other floats to
# Dual (partial 0), leave Int/Bool index/mask arrays plain. Compare device-AD vs
# CPU-AD partials (PASS) and the device-AD partial vs a finite difference (sanity).
# ==========================================================================
const STNLEVSNO = 5; const STNLEVGRND = 10; const STNLEVURB = 5; const STNLEVMAX = 10
const STNLEV = STNLEVSNO + STNLEVMAX            # 15
const STNLEVTOT = STNLEVSNO + 1 + STNLEVMAX     # 16
const STNBAND = 5

# Branch-exercising Float32 inputs (mirror gpu_validate_soiltemp.jl).
function _stbuild(::Type{FT}) where {FT}
    rng = MersenneTwister(7); nc = 8
    rnd(d...) = FT.(rand(rng, d...))
    (; nc,
        itype     = Int[1,1,1,1, CLM.ICOL_SUNWALL, CLM.ICOL_ROOF, CLM.ICOL_ROAD_PERV, CLM.ICOL_ROAD_IMPERV],
        urbpoi    = Bool[false,false,false,false,true,true,true,true],
        lun_itype = Int[CLM.ISTSOIL, CLM.ISTCROP, CLM.ISTSOIL, CLM.ISTSOIL, 70,70,70,70],
        landunit  = collect(1:nc),
        snl       = Int[0,-1,-2,-3,-1,0,-2,0],
        mask      = trues(nc),
        z  = FT[0.1*jj for c in 1:nc, jj in 1:STNLEV],
        zi = FT[0.1*jj-0.05 for c in 1:nc, jj in 1:(STNLEV+1)],
        dz = fill(FT(0.1), nc, STNLEV),
        t_soisno = FT[273.15 + 1.5*sinpi((c+jj)/5) for c in 1:nc, jj in 1:STNLEV],
        t_h2osfc = FT[272.0 + 0.5*c for c in 1:nc],
        tk = FT(1.5) .+ rnd(nc, STNLEV), fn = rnd(nc, STNLEV),
        fact = FT(1.0e-3) .+ rnd(nc, STNLEV).*FT(1e-4),
        dhsdT = FT(-2) .- rnd(nc),
        hs_top = rnd(nc).*50, hs_top_snow = rnd(nc).*50, hs_soil = rnd(nc).*50, hs_h2osfc = rnd(nc).*50,
        sabg_lyr_col = rnd(nc, STNLEVSNO+1),
        tk_h2osfc = FT(0.5) .+ rnd(nc), c_h2osfc = FT(4.0e4) .+ rnd(nc).*1000,
        dz_h2osfc = FT(0.01) .+ rnd(nc).*FT(0.01),
        frac_h2osfc  = FT[0.0,0.3,0.0,0.4,0.0,0.0,0.2,0.0],
        frac_sno_eff = FT[0.0,0.5,0.8,1.0,0.5,0.0,0.6,0.0],
        snow_depth = FT(0.1) .+ rnd(nc), h2osfc = FT(0.5) .+ rnd(nc),
        h2osno_no_layers = FT[0.0,2.0,0.0,5.0,1.0,0.0,3.0,0.0],
        int_snow = rnd(nc).*10,
        h2osoi_ice = FT(2.0) .+ rnd(nc, STNLEV).*5, h2osoi_liq = FT(2.0) .+ rnd(nc, STNLEV).*5,
        excess_ice = rnd(nc, STNLEVMAX).*3,
        bsw = FT(4.0) .+ rnd(nc, STNLEVMAX), sucsat = FT(100.0) .+ rnd(nc, STNLEVMAX).*100,
        watsat = FT(0.4) .+ rnd(nc, STNLEVMAX).*FT(0.1))
end

_sd(::Type{FT}, A) where {FT} = (D = _todual(FT, A); @inbounds for i in eachindex(D); D[i] = _seed(A[i], one(FT)); end; D)  # seeded partial 1
_devmap(to, args) = map(a -> a isa AbstractArray ? _dev(to, a) : a, args)

# Generic `_launch!`-kernel AD test. `mk(sfn, ofn)` returns the positional args after `out`,
# applying sfn to the SEEDED float array and ofn to the other float arrays (Int/Bool passed
# plain). dims/ndr describe `out`. `out0val` (Float32) seeds a nonzero output (for += / -= kernels).
function _adkernel(name, to, ::Type{FT}, K, dims, ndr, mk, h; out0val = nothing) where {FT}
    args_d = mk(x -> _sd(FT, x), x -> _todual(FT, x))
    oc = out0val === nothing ? zeros(_dualT(FT), dims...) : _todual(FT, out0val)
    CLM._launch!(K, oc, args_d...; ndrange = ndr)
    od = out0val === nothing ? to(zeros(_dualT(FT), dims...)) : to(_todual(FT, out0val))
    CLM._launch!(K, od, _devmap(to, args_d)...; ndrange = ndr)
    parity = maxabs(_part(oc) .- _part(od))
    o0 = out0val === nothing ? zeros(FT, dims...) : copy(out0val)
    CLM._launch!(K, o0, mk(identity, identity)...; ndrange = ndr)
    o1 = out0val === nothing ? zeros(FT, dims...) : copy(out0val)
    CLM._launch!(K, o1, mk(x -> x .+ h, identity)...; ndrange = ndr)
    fd = (o1 .- o0) ./ h
    fdrel = relerr(_part(od), fd)
    return (name, parity, fdrel, parity < _parity_tol(FT) && fdrel < FD_REL_TOL)
end

adtest_rhs_snow(to, ::Type{FT}) where {FT} = (B = _stbuild(FT);
    _adkernel("rhs_snow", to, FT, CLM._rhs_snow_kernel!, (B.nc, STNLEVTOT), (B.nc, STNLEVSNO),
        (s, o) -> (B.mask, B.itype, B.snl, s(B.t_soisno), o(B.fact), o(B.fn), o(B.dhsdT),
                   o(B.hs_top_snow), o(B.hs_top), o(B.sabg_lyr_col), STNLEVSNO), FT(1f-2)))

adtest_rhs_ssw(to, ::Type{FT}) where {FT} = (B = _stbuild(FT); dt = FT(1800);
    _adkernel("rhs_ssw", to, FT, CLM._rhs_ssw_kernel!, (B.nc, STNLEVTOT), B.nc,
        (s, o) -> (B.mask, B.itype, o(B.t_soisno), s(B.t_h2osfc), o(B.z), o(B.tk_h2osfc), o(B.dz_h2osfc),
                   o(B.c_h2osfc), o(B.hs_h2osfc), o(B.dhsdT), o(zeros(FT, B.nc)), dt, STNLEVSNO), FT(1f-2)))

adtest_rhs_soil_urban(to, ::Type{FT}) where {FT} = (B = _stbuild(FT);
    _adkernel("rhs_soil_urban", to, FT, CLM._rhs_soil_urban_kernel!, (B.nc, STNLEVTOT), (B.nc, STNLEVURB),
        (s, o) -> (B.mask, B.itype, B.snl, s(B.t_soisno), o(B.fact), o(B.fn), o(B.dhsdT),
                   o(B.hs_top), o(B.sabg_lyr_col), STNLEVSNO, STNLEVURB), FT(1f-2)))

adtest_rhs_soil(to, ::Type{FT}) where {FT} = (B = _stbuild(FT);
    _adkernel("rhs_soil", to, FT, CLM._rhs_soil_kernel!, (B.nc, STNLEVTOT), (B.nc, STNLEVGRND),
        (s, o) -> (B.mask, B.itype, B.snl, B.landunit, B.urbpoi, s(B.t_soisno), o(B.fact), o(B.fn),
                   o(B.dhsdT), o(B.frac_sno_eff), o(B.hs_top_snow), o(B.hs_soil), o(B.sabg_lyr_col),
                   STNLEVSNO, STNLEVGRND), FT(1f-2)))

adtest_rhs_h2osfc_corr(to, ::Type{FT}) where {FT} = (B = _stbuild(FT);
    rv0 = FT.(1.0 .+ rand(MersenneTwister(11), B.nc, STNLEVTOT));
    _adkernel("rhs_h2osfc_corr", to, FT, CLM._rhs_h2osfc_corr_kernel!, (B.nc, STNLEVTOT), B.nc,
        (s, o) -> (B.mask, o(B.frac_h2osfc), o(B.fact), o(B.dhsdT), s(B.t_soisno), o(B.hs_soil),
                   o(B.tk_h2osfc), STNLEVSNO), FT(1f-2); out0val = rv0))

adtest_mat_snow(to, ::Type{FT}) where {FT} = (B = _stbuild(FT);
    _adkernel("mat_snow", to, FT, CLM._mat_snow_kernel!, (B.nc, STNBAND, STNLEVTOT), (B.nc, STNLEVSNO),
        (s, o) -> (B.mask, B.snl, o(B.z), s(B.tk), o(B.fact), o(B.dhsdT), STNLEVSNO), FT(1f-2)))

adtest_mat_ssw(to, ::Type{FT}) where {FT} = (B = _stbuild(FT); dt = FT(1800);
    _adkernel("mat_ssw", to, FT, CLM._mat_ssw_kernel!, (B.nc, STNBAND, STNLEVTOT), B.nc,
        (s, o) -> (B.mask, B.itype, o(B.z), s(B.tk_h2osfc), o(B.dz_h2osfc), o(B.c_h2osfc), o(B.dhsdT),
                   dt, STNLEVSNO), FT(1f-2)))

adtest_mat_soil_urban(to, ::Type{FT}) where {FT} = (B = _stbuild(FT);
    _adkernel("mat_soil_urban", to, FT, CLM._mat_soil_urban_kernel!, (B.nc, STNBAND, STNLEVTOT), (B.nc, STNLEVURB),
        (s, o) -> (B.mask, B.itype, B.snl, o(B.z), o(B.zi), s(B.tk), o(B.fact), o(B.dhsdT),
                   STNLEVSNO, STNLEVURB), FT(1f-2)))

adtest_mat_soil(to, ::Type{FT}) where {FT} = (B = _stbuild(FT);
    _adkernel("mat_soil", to, FT, CLM._mat_soil_kernel!, (B.nc, STNBAND, STNLEVTOT), (B.nc, STNLEVGRND),
        (s, o) -> (B.mask, B.itype, B.snl, B.landunit, B.urbpoi, o(B.z), s(B.tk), o(B.fact),
                   o(B.dhsdT), o(B.frac_sno_eff), STNLEVSNO, STNLEVGRND), FT(1f-2)))

adtest_mat_h2osfc_corr(to, ::Type{FT}) where {FT} = (B = _stbuild(FT);
    bm0 = FT.(1.0 .+ rand(MersenneTwister(13), B.nc, STNBAND, STNLEVTOT));
    _adkernel("mat_h2osfc_corr", to, FT, CLM._mat_h2osfc_corr_kernel!, (B.nc, STNBAND, STNLEVTOT), B.nc,
        (s, o) -> (B.mask, o(B.frac_h2osfc), o(B.z), s(B.tk_h2osfc), o(B.dz_h2osfc), o(B.fact),
                   o(B.dhsdT), STNLEVSNO), FT(1f-2); out0val = bm0))

# phase_change_h2osfc!: seed t_h2osfc; check the t_soisno output partial (CPU vs device) + FD.
# Every WRITABLE float array receives Dual-derived values, so all floats are lifted to Dual
# (L = _todual / identity); only t_h2osfc is seeded (SF). dm moves to device (or identity on CPU).
function adtest_pc_h2osfc(to, ::Type{FT}) where {FT}
    B = _stbuild(FT); nc = B.nc; dt = FT(1800); pck = FT(CLM.PHASE_CHANGE_MASS_K[])
    function go(L, SF, dm)
        th = dm(SF(B.t_h2osfc)); ts = dm(L(B.t_soisno))
        CLM._launch!(CLM._phase_change_h2osfc_kernel!, th, ts, dm(L(B.h2osfc)),
            dm(L(B.h2osno_no_layers)), dm(L(B.h2osoi_ice)), dm(L(B.snow_depth)), dm(L(B.int_snow)),
            dm(L(zeros(FT, nc))), dm(L(zeros(FT, nc))), dm(L(zeros(FT, nc))), dm(B.mask), dm(B.snl),
            dm(L(B.fact)), dm(L(B.c_h2osfc)), dm(L(B.frac_sno_eff)), dm(L(B.frac_h2osfc)),
            dm(L(B.h2osoi_liq)), dm(L(B.dhsdT)), dt, STNLEVSNO, pck)
        return ts
    end
    td(x) = _todual(FT, x); sdl(x) = _sd(FT, x); h = FT(1f-2)
    tsc = go(td, sdl, identity)
    tsd = go(td, sdl, x -> _dev(to, x))
    parity = maxabs(_part(tsc) .- _part(tsd))
    ts0 = go(identity, identity, identity)
    ts1 = go(identity, x -> x .+ h, identity)
    fd = (Array(ts1) .- Array(ts0)) ./ h
    fdrel = relerr(_part(tsd), fd)
    # Phase change has hard melt/freeze branch decisions: a finite-h perturbation can flip a
    # cell across TFRZ and produce a spurious FD jump, so FD-vs-AD is not a reliable gate here.
    # Gate on PARITY (device-AD == CPU-AD) plus a finite, live (nonzero) partial.
    live = isfinite(parity) && maxabs(_part(tsd)) > 0
    return ("phase_change_h2osfc*", parity, fdrel, parity < _parity_tol(FT) && live)
end

# phase_change_beta!: seed t_soisno; check the h2osoi_ice output partial (smooth via smooth_max/min).
# All struct float fields share one element type, so every float is lifted to Dual (L); only
# t_soisno is seeded (SF). dm moves each array to device (or identity on CPU).
function adtest_pc_beta(to, ::Type{FT}) where {FT}
    B = _stbuild(FT); nc = B.nc; dt = FT(1800); pck = FT(CLM.PHASE_CHANGE_MASS_K[])
    Z(d2) = zeros(FT, nc, d2)
    function go(L, SF, dm)
        Lyr = CLM.PcbLyr(; t_soisno = dm(SF(B.t_soisno)), h2osoi_ice = dm(L(B.h2osoi_ice)),
            h2osoi_liq = dm(L(B.h2osoi_liq)), excess_ice = dm(L(B.excess_ice)),
            exice_subs = dm(L(Z(STNLEVMAX))), qflx_snomelt_lyr = dm(L(Z(STNLEV))),
            qflx_snofrz_lyr = dm(L(Z(STNLEV))))
        Col = CLM.PcbCol(; xmf = dm(L(zeros(FT, nc))), qflx_snomelt = dm(L(zeros(FT, nc))),
            qflx_snofrz = dm(L(zeros(FT, nc))), qflx_snow_drain = dm(L(zeros(FT, nc))),
            h2osno_no_layers = dm(L(B.h2osno_no_layers)), snow_depth = dm(L(B.snow_depth)),
            snomelt_accum = dm(L(zeros(FT, nc))), eflx_snomelt = dm(L(zeros(FT, nc))),
            eflx_snomelt_r = dm(L(zeros(FT, nc))), eflx_snomelt_u = dm(L(zeros(FT, nc))))
        Pin = CLM.PcbIn(; dz = dm(L(B.dz)), fact = dm(L(B.fact)), bsw = dm(L(B.bsw)),
            sucsat = dm(L(B.sucsat)), watsat = dm(L(B.watsat)), dhsdT = dm(L(B.dhsdT)),
            frac_sno_eff = dm(L(B.frac_sno_eff)), frac_h2osfc = dm(L(B.frac_h2osfc)))
        Tmp = CLM.PcbTmp(; hm = dm(L(Z(STNLEV))), xm = dm(L(Z(STNLEV))), xm2 = dm(L(Z(STNLEV))),
            wice0 = dm(L(Z(STNLEV))), wliq0 = dm(L(Z(STNLEV))), wexice0 = dm(L(Z(STNLEV))),
            wmass0 = dm(L(Z(STNLEV))), supercool = dm(L(Z(STNLEVMAX))), tinc = dm(L(Z(STNLEV))))
        be = CLM._kernel_backend(Lyr.t_soisno)
        CLM._phase_change_beta_kernel!(be)(Lyr, Col, Pin, Tmp, dm(zeros(Int, nc, STNLEV)),
            dm(B.mask), dm(B.urbpoi), dm(B.snl), dm(B.landunit), dm(B.itype), dm(B.lun_itype),
            dt, STNLEVSNO, STNLEVGRND, STNLEVURB, STNLEVMAX, pck; ndrange = nc)
        KernelAbstractions.synchronize(be)
        return Lyr
    end
    td(x) = _todual(FT, x); sdl(x) = _sd(FT, x); h = FT(1f-3)
    Lc = go(td, sdl, identity)
    Ld = go(td, sdl, x -> _dev(to, x))
    parity = maxabs(_part(Lc.h2osoi_ice) .- _part(Ld.h2osoi_ice))
    L0 = go(identity, identity, identity)
    L1 = go(identity, x -> x .+ h, identity)
    fd = (Array(L1.h2osoi_ice) .- Array(L0.h2osoi_ice)) ./ h
    fdrel = relerr(_part(Ld.h2osoi_ice), fd)
    # See phase_change_h2osfc: FD across hard melt/freeze boundaries is unreliable; gate on
    # PARITY + a finite, live partial.
    live = isfinite(parity) && maxabs(_part(Ld.h2osoi_ice)) > 0
    return ("phase_change_beta*", parity, fdrel, parity < _parity_tol(FT) && live)
end

const AD_TESTS = [adtest_forc_q, adtest_tridiagonal, adtest_band,
    adtest_rhs_snow, adtest_rhs_ssw, adtest_rhs_soil_urban, adtest_rhs_soil, adtest_rhs_h2osfc_corr,
    adtest_mat_snow, adtest_mat_ssw, adtest_mat_soil_urban, adtest_mat_soil, adtest_mat_h2osfc_corr,
    adtest_pc_h2osfc, adtest_pc_beta]

function main(backend)
    println("=" ^ 74)
    println("GPU AUTODIFF VALIDATION for CLM.jl kernels (forward-mode, CPU vs device)")
    println("=" ^ 74)
    if backend === nothing
        println("""
  No GPU backend detected (CUDA / AMDGPU / oneAPI / Metal).
  Forward-mode AD through the kernels is exercised on the CPU by the test suite.
  Install a backend package and re-run on real GPU hardware.
""")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n", name, FT)
    println("  PARITY = max|device-AD partial − CPU-AD partial|   (must be ~0)")
    println("  FD-rel = relative error of AD partial vs finite difference (sanity)\n")

    npass = 0; nfail = 0
    for t in AD_TESTS
        nm, parity, fdrel, ok = try
            t(to, FT)
        catch e
            (string(t) * "  (ERROR: " * first(split(sprint(showerror, e), Char(10))) * ")", NaN, NaN, false)
        end
        @printf("  [%s] %-42s  PARITY=%.2e  FD-rel=%.2e\n", ok ? "PASS" : "FAIL", nm, parity, fdrel)
        ok ? (npass += 1) : (nfail += 1)
    end
    @printf("\n  %d passed, %d failed\n", npass, nfail)
    println("  (*) phase_change_* gated on PARITY only; their FD-rel is large by design —")
    println("      a finite-h probe flips melt/freeze branches at TFRZ, so FD≠AD there.")
    if nfail == 0
        println("\n  FORWARD-MODE AD MATCHES CPU ON $(name) ($(FT)) ✓")
        println("  Reverse-mode (Enzyme) through these kernels works on CPU; Enzyme does not")
        println("  differentiate the Metal device launch (no Enzyme GPU support for Apple).")
    else
        println("\n  SOME AD CHECKS FAILED — investigate before trusting GPU gradients.")
    end
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()   # top-level: advance world age before main()
exit(main(BACKEND))
