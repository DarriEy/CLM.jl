# ==========================================================================
# gpu_validate_soiltemp.jl — real-hardware GPU parity for the soil_temperature!
# kernels (Metal pipeline, functions 4-8). Same spirit as gpu_validate.jl:
# run each @kernel on a CPU Array and on a device array of the SAME inputs at the
# SAME precision and compare. Parity ("does the device execute the kernel the way
# the CPU does?"), not accuracy-vs-Fortran.
#
#   julia --project=scripts scripts/gpu_validate_soiltemp.jl
#
# Kernels covered:
#   fn4 set_rhs_vec!         : _rhs_snow/_rhs_ssw/_rhs_soil_urban/_rhs_soil/_rhs_h2osfc_corr
#   fn5 set_matrix!          : _mat_snow/_mat_ssw/_mat_soil_urban/_mat_soil/_mat_h2osfc_corr
#   fn6 building_hac!        : CPU-pinned by design (BitVector outputs) — N/A on device
#   fn7 phase_change_h2osfc! : _phase_change_h2osfc
#   fn8 phase_change_beta!   : _phase_change_beta (device-view structs PcbIn/Lyr/Col/Tmp)
# ==========================================================================

using CLM
using Printf
using Random
import KernelAbstractions  # for synchronize() after the struct-kernel launch

include(joinpath(@__DIR__, "gpu_backends.jl"))

maxabsdiff(a, b) = maximum(abs.(Array(a) .- Array(b)); init = 0.0)
combinediff(pairs...) = maximum(maxabsdiff(a, b) for (a, b) in pairs)
TOL(::Type{Float32}) = 5f-3
TOL(::Type{T}) where {T} = 1e-9

# Dimensions (mirror test_soil_temperature.jl).
const NLEVSNO = 5
const NLEVGRND = 10
const NLEVURB = 5
const NLEVMAX = 10          # nlevmaxurbgrnd
const JOFF = NLEVSNO
const NLEV = NLEVSNO + NLEVMAX            # layer-array 2nd dim (15)
const NLEVTOT = NLEVSNO + 1 + NLEVMAX     # rvector/bmatrix 2nd/3rd dim (16)
const NBAND = 5

# Build a realistic, branch-exercising input set on the CPU at precision FT.
# nc = 8 columns: 1-4 non-urban soil (snl 0,-1,-2,-3), 5 sunwall, 6 roof, 7 road-perv, 8 road-imperv.
function build_inputs(::Type{FT}) where {FT}
    rng = MersenneTwister(7)
    nc = 8
    itype     = Int[1, 1, 1, 1, CLM.ICOL_SUNWALL, CLM.ICOL_ROOF, CLM.ICOL_ROAD_PERV, CLM.ICOL_ROAD_IMPERV]
    urbpoi    = Bool[false, false, false, false, true, true, true, true]
    lun_itype = Int[CLM.ISTSOIL, CLM.ISTCROP, CLM.ISTSOIL, CLM.ISTSOIL, 70, 70, 70, 70]
    landunit  = collect(1:nc)
    snl       = Int[0, -1, -2, -3, -1, 0, -2, 0]
    mask      = trues(nc)

    rnd(d...) = FT.(rand(rng, d...))
    # strictly increasing depths so all (z[jj+1]-z[jj]) and (z[jj]-z[jj-1]) are > 0
    z  = FT[0.1 * jj for c in 1:nc, jj in 1:NLEV]
    zi = FT[0.1 * jj - 0.05 for c in 1:nc, jj in 1:(NLEV + 1)]
    dz = fill(FT(0.1), nc, NLEV)

    # temperatures spanning TFRZ to exercise melt (>TFRZ) and freeze (<TFRZ) branches
    t_soisno = FT[273.15 + 1.5 * sinpi((c + jj) / 5) for c in 1:nc, jj in 1:NLEV]
    t_h2osfc = FT[272.0 + 0.5 * c for c in 1:nc]

    return (; nc, itype, urbpoi, lun_itype, landunit, snl, mask,
        z, zi, dz, t_soisno, t_h2osfc,
        tk   = FT(1.5) .+ rnd(nc, NLEV),
        cv   = FT(1.0e6) .+ rnd(nc, NLEV),
        fn   = rnd(nc, NLEV),
        fact = FT(1.0e-3) .+ rnd(nc, NLEV) .* FT(1e-4),
        dhsdT = FT(-2) .- rnd(nc),
        hs_top = rnd(nc) .* 50, hs_top_snow = rnd(nc) .* 50,
        hs_soil = rnd(nc) .* 50, hs_h2osfc = rnd(nc) .* 50,
        sabg_lyr_col = rnd(nc, NLEVSNO + 1),
        tk_h2osfc = FT(0.5) .+ rnd(nc),
        c_h2osfc  = FT(4.0e4) .+ rnd(nc) .* 1000,
        dz_h2osfc = FT(0.01) .+ rnd(nc) .* FT(0.01),
        frac_h2osfc  = FT[0.0, 0.3, 0.0, 0.4, 0.0, 0.0, 0.2, 0.0],
        frac_sno_eff = FT[0.0, 0.5, 0.8, 1.0, 0.5, 0.0, 0.6, 0.0],
        snow_depth = FT(0.1) .+ rnd(nc),
        h2osfc = FT(0.5) .+ rnd(nc),
        h2osno_no_layers = FT[0.0, 2.0, 0.0, 5.0, 1.0, 0.0, 3.0, 0.0],
        int_snow = rnd(nc) .* 10,
        h2osoi_ice = FT(2.0) .+ rnd(nc, NLEV) .* 5,
        h2osoi_liq = FT(2.0) .+ rnd(nc, NLEV) .* 5,
        excess_ice = rnd(nc, NLEVMAX) .* 3,
        bsw = FT(4.0) .+ rnd(nc, NLEVMAX),
        sucsat = FT(100.0) .+ rnd(nc, NLEVMAX) .* 100,
        watsat = FT(0.4) .+ rnd(nc, NLEVMAX) .* FT(0.1),
        FT = FT)
end

# Run a `_launch!`-style kernel on CPU then device (inputs moved with _dev) and diff `out`
# plus any extra written arrays. `mkargs(arrs)` returns the positional kernel args given a
# function mapping each named array to its (cpu or device) version.
function _parity(name, to, FT, K, ndr, out0, mk; extra_out = (), tol = TOL(FT))
    # CPU run
    out_cpu = copy(out0)
    cpu(x) = x
    args_cpu, extras_cpu = mk(cpu)
    CLM._launch!(K, out_cpu, args_cpu...; ndrange = ndr)
    # device run
    out_dev = to(copy(out0))
    args_dev, extras_dev = mk(x -> _dev(to, x))
    CLM._launch!(K, out_dev, args_dev...; ndrange = ndr)
    d = max(maxabsdiff(out_cpu, out_dev),
            maximum((maxabsdiff(a, b) for (a, b) in zip(extras_cpu, extras_dev)); init = 0.0))
    return (name, d, d < tol)
end

function tests(to, ::Type{FT}) where {FT}
    I = build_inputs(FT)
    nc = I.nc
    dt = FT(1800)
    results = Tuple{String,Float64,Bool}[]

    # ---- fn4: set_rhs_vec! kernels (rvector starts at 0 here so unwritten cells match) ----
    rv0 = zeros(FT, nc, NLEVTOT)
    push!(results, _parity("rhs_snow", to, FT, CLM._rhs_snow_kernel!, (nc, NLEVSNO), rv0,
        f -> ((f(I.mask), f(I.itype), f(I.snl), f(I.t_soisno), f(I.fact), f(I.fn), f(I.dhsdT),
               f(I.hs_top_snow), f(I.hs_top), f(I.sabg_lyr_col), NLEVSNO), ())))
    # ssw also writes fn_h2osfc
    push!(results, _parity("rhs_ssw", to, FT, CLM._rhs_ssw_kernel!, nc, rv0,
        f -> begin fh = f(zeros(FT, nc));
            # `itype` was added to _rhs_ssw_kernel! (soil_temperature.jl:1324) to
            # decouple the urban roof/wall standing-water row, and this call was never
            # updated — 13 args against a 14-arg kernel, so it died with a MethodError
            # before running. Its two sibling calls just below already pass f(I.itype).
            ((f(I.mask), f(I.itype), f(I.t_soisno), f(I.t_h2osfc), f(I.z), f(I.tk_h2osfc), f(I.dz_h2osfc),
              f(I.c_h2osfc), f(I.hs_h2osfc), f(I.dhsdT), fh, dt, NLEVSNO), (fh,)) end))
    push!(results, _parity("rhs_soil_urban", to, FT, CLM._rhs_soil_urban_kernel!, (nc, NLEVURB), rv0,
        f -> ((f(I.mask), f(I.itype), f(I.snl), f(I.t_soisno), f(I.fact), f(I.fn), f(I.dhsdT),
               f(I.hs_top), f(I.sabg_lyr_col), NLEVSNO, NLEVURB), ())))
    push!(results, _parity("rhs_soil", to, FT, CLM._rhs_soil_kernel!, (nc, NLEVGRND), rv0,
        f -> ((f(I.mask), f(I.itype), f(I.snl), f(I.landunit), f(I.urbpoi), f(I.t_soisno),
               f(I.fact), f(I.fn), f(I.dhsdT), f(I.frac_sno_eff), f(I.hs_top_snow), f(I.hs_soil),
               f(I.sabg_lyr_col), NLEVSNO, NLEVGRND), ())))
    rv_corr = FT.(1.0 .+ rand(MersenneTwister(11), nc, NLEVTOT))  # nonzero so -= is exercised
    fnh_in = FT.(rand(MersenneTwister(12), nc))
    push!(results, _parity("rhs_h2osfc_corr", to, FT, CLM._rhs_h2osfc_corr_kernel!, nc, rv_corr,
        f -> ((f(I.mask), f(I.frac_h2osfc), f(I.fact), f(I.dhsdT), f(I.t_soisno), f(I.hs_soil),
               f(fnh_in), NLEVSNO), ())))

    # ---- fn5: set_matrix! kernels (bmatrix starts at 0) ----
    bm0 = zeros(FT, nc, NBAND, NLEVTOT)
    push!(results, _parity("mat_snow", to, FT, CLM._mat_snow_kernel!, (nc, NLEVSNO), bm0,
        f -> ((f(I.mask), f(I.snl), f(I.z), f(I.tk), f(I.fact), f(I.dhsdT), NLEVSNO), ())))
    push!(results, _parity("mat_ssw", to, FT, CLM._mat_ssw_kernel!, nc, bm0,
        # Same missing `itype` as the rhs_ssw call above (soil_temperature.jl:1513).
        f -> ((f(I.mask), f(I.itype), f(I.z), f(I.tk_h2osfc), f(I.dz_h2osfc), f(I.c_h2osfc), f(I.dhsdT),
               dt, NLEVSNO), ())))
    push!(results, _parity("mat_soil_urban", to, FT, CLM._mat_soil_urban_kernel!, (nc, NLEVURB), bm0,
        f -> ((f(I.mask), f(I.itype), f(I.snl), f(I.z), f(I.zi), f(I.tk), f(I.fact), f(I.dhsdT),
               NLEVSNO, NLEVURB), ())))
    push!(results, _parity("mat_soil", to, FT, CLM._mat_soil_kernel!, (nc, NLEVGRND), bm0,
        f -> ((f(I.mask), f(I.itype), f(I.snl), f(I.landunit), f(I.urbpoi), f(I.z), f(I.tk),
               f(I.fact), f(I.dhsdT), f(I.frac_sno_eff), NLEVSNO, NLEVGRND), ())))
    bm_corr = FT.(1.0 .+ rand(MersenneTwister(13), nc, NBAND, NLEVTOT))
    push!(results, _parity("mat_h2osfc_corr", to, FT, CLM._mat_h2osfc_corr_kernel!, nc, bm_corr,
        f -> ((f(I.mask), f(I.frac_h2osfc), f(I.z), f(I.tk_h2osfc), f(I.dz_h2osfc), f(I.fact),
               f(I.dhsdT), NLEVSNO), ())))

    # ---- fn7: phase_change_h2osfc! (in-place on many column/layer arrays) ----
    push!(results, test_pc_h2osfc(to, FT, I, dt))
    # ---- fn8: phase_change_beta! (device-view structs) ----
    push!(results, test_pc_beta(to, FT, I, dt))

    return results
end

# phase_change_h2osfc!: build cpu + device copies of every written array, run, compare.
function test_pc_h2osfc(to, ::Type{FT}, I, dt) where {FT}
    nc = I.nc
    mkset(cv) = (t_h2osfc = cv(copy(I.t_h2osfc)), t_soisno = cv(copy(I.t_soisno)),
        h2osfc = cv(copy(I.h2osfc)), h2osno = cv(copy(I.h2osno_no_layers)),
        ice = cv(copy(I.h2osoi_ice)), sd = cv(copy(I.snow_depth)), isnow = cv(copy(I.int_snow)),
        xmf = cv(zeros(FT, nc)), q2i = cv(zeros(FT, nc)), e2s = cv(zeros(FT, nc)))
    run!(s, to_) = (CLM._launch!(CLM._phase_change_h2osfc_kernel!, s.t_h2osfc, s.t_soisno, s.h2osfc,
        s.h2osno, s.ice, s.sd, s.isnow, s.xmf, s.q2i, s.e2s, _dev(to_, I.mask), _dev(to_, I.snl),
        _dev(to_, I.fact), _dev(to_, I.c_h2osfc), _dev(to_, I.frac_sno_eff), _dev(to_, I.frac_h2osfc),
        # Trailing `pck`: the PHASE_CHANGE_MASS_K smoothing constant, read on the HOST
        # and passed in as a scalar (soil_temperature.jl:1822) precisely so the kernel
        # never dereferences a host Ref on a device. Added to the kernel after this
        # harness was written — 20 args against a 21-arg kernel, MethodError before
        # running. Convert at the working eltype, exactly as the real launch does.
        _dev(to_, I.h2osoi_liq), _dev(to_, I.dhsdT), dt, NLEVSNO,
        convert(FT, CLM.PHASE_CHANGE_MASS_K[])); s)
    cpu = mkset(identity); run!(cpu, identity)
    dev = mkset(x -> _dev(to, x)); run!(dev, to)
    d = combinediff((cpu.t_h2osfc, dev.t_h2osfc), (cpu.t_soisno, dev.t_soisno),
        (cpu.h2osfc, dev.h2osfc), (cpu.h2osno, dev.h2osno), (cpu.ice, dev.ice),
        (cpu.sd, dev.sd), (cpu.isnow, dev.isnow), (cpu.xmf, dev.xmf),
        (cpu.q2i, dev.q2i), (cpu.e2s, dev.e2s))
    return ("phase_change_h2osfc", d, d < TOL(FT))
end

# phase_change_beta!: build the four device-view structs on each backend, run, compare.
function test_pc_beta(to, ::Type{FT}, I, dt) where {FT}
    nc = I.nc
    mklyr(cv) = CLM.PcbLyr(; t_soisno = cv(copy(I.t_soisno)), h2osoi_ice = cv(copy(I.h2osoi_ice)),
        h2osoi_liq = cv(copy(I.h2osoi_liq)), excess_ice = cv(copy(I.excess_ice)),
        exice_subs = cv(zeros(FT, nc, NLEVMAX)),
        qflx_snomelt_lyr = cv(zeros(FT, nc, NLEV)), qflx_snofrz_lyr = cv(zeros(FT, nc, NLEV)))
    mkcol(cv) = CLM.PcbCol(; xmf = cv(zeros(FT, nc)), qflx_snomelt = cv(zeros(FT, nc)),
        qflx_snofrz = cv(zeros(FT, nc)), qflx_snow_drain = cv(zeros(FT, nc)),
        h2osno_no_layers = cv(copy(I.h2osno_no_layers)), snow_depth = cv(copy(I.snow_depth)),
        snomelt_accum = cv(zeros(FT, nc)), eflx_snomelt = cv(zeros(FT, nc)),
        eflx_snomelt_r = cv(zeros(FT, nc)), eflx_snomelt_u = cv(zeros(FT, nc)))
    mkin(cv, cm) = CLM.PcbIn(; dz = cm(I.dz), fact = cm(I.fact), bsw = cm(I.bsw),
        sucsat = cm(I.sucsat), watsat = cm(I.watsat), dhsdT = cv(I.dhsdT),
        frac_sno_eff = cv(I.frac_sno_eff), frac_h2osfc = cv(I.frac_h2osfc))
    mktmp(cv) = CLM.PcbTmp(; hm = cv(zeros(FT, nc, NLEV)), xm = cv(zeros(FT, nc, NLEV)),
        xm2 = cv(zeros(FT, nc, NLEV)), wice0 = cv(zeros(FT, nc, NLEV)), wliq0 = cv(zeros(FT, nc, NLEV)),
        wexice0 = cv(zeros(FT, nc, NLEV)), wmass0 = cv(zeros(FT, nc, NLEV)),
        supercool = cv(zeros(FT, nc, NLEVMAX)), tinc = cv(zeros(FT, nc, NLEV)))

    function run!(lyr, colv, pin, tmp, imelt, to_)
        backend = CLM._kernel_backend(lyr.t_soisno)
        CLM._phase_change_beta_kernel!(backend)(lyr, colv, pin, tmp, imelt,
            _dev(to_, I.mask), _dev(to_, I.urbpoi), _dev(to_, I.snl), _dev(to_, I.landunit),
            _dev(to_, I.itype), _dev(to_, I.lun_itype), dt, NLEVSNO, NLEVGRND, NLEVURB, NLEVMAX,
            # Same host-hoisted PHASE_CHANGE_MASS_K trailing scalar as the fn7 call
            # above (real launch: soil_temperature.jl:2192).
            convert(FT, CLM.PHASE_CHANGE_MASS_K[]);
            ndrange = nc)
        KernelAbstractions.synchronize(backend)
    end

    lyr_c, col_c, in_c, tmp_c = mklyr(identity), mkcol(identity), mkin(identity, identity), mktmp(identity)
    im_c = zeros(Int, nc, NLEV)
    run!(lyr_c, col_c, in_c, tmp_c, im_c, identity)
    cv = x -> _dev(to, x)
    lyr_d, col_d, in_d, tmp_d, im_d = mklyr(cv), mkcol(cv), mkin(cv, cv), mktmp(cv), to(zeros(Int, nc, NLEV))
    run!(lyr_d, col_d, in_d, tmp_d, im_d, to)

    d = combinediff((lyr_c.t_soisno, lyr_d.t_soisno), (lyr_c.h2osoi_ice, lyr_d.h2osoi_ice),
        (lyr_c.h2osoi_liq, lyr_d.h2osoi_liq), (lyr_c.excess_ice, lyr_d.excess_ice),
        (lyr_c.exice_subs, lyr_d.exice_subs), (lyr_c.qflx_snomelt_lyr, lyr_d.qflx_snomelt_lyr),
        (lyr_c.qflx_snofrz_lyr, lyr_d.qflx_snofrz_lyr), (col_c.xmf, col_d.xmf),
        (col_c.qflx_snomelt, col_d.qflx_snomelt), (col_c.qflx_snofrz, col_d.qflx_snofrz),
        (col_c.h2osno_no_layers, col_d.h2osno_no_layers), (col_c.snow_depth, col_d.snow_depth),
        (col_c.snomelt_accum, col_d.snomelt_accum), (col_c.eflx_snomelt, col_d.eflx_snomelt),
        (Array(im_c), Array(im_d)))
    return ("phase_change_beta", d, d < TOL(FT))
end

function main(backend)
    println("=" ^ 70)
    println("METAL PARITY for soil_temperature! kernels (fn 4-8)")
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
        @printf("  [%s] %-26s  max|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", knm, d)
        ok ? (npass += 1) : (nfail += 1)
    end
    println()
    @printf("  building_hac          : CPU-pinned by design (BitVector outputs) — N/A on device\n")
    @printf("\n  %d passed, %d failed\n", npass, nfail)
    println(nfail == 0 ? "\n  ALL soil_temperature! kernels MATCH CPU ON $name ($FT) ✓" :
                          "\n  SOME KERNELS DIVERGED — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
