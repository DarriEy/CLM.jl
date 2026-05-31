# ==========================================================================
# gpu_validate_lake.jl — real-hardware GPU parity for the lake_temperature!
# kernels. Same spirit as gpu_validate_soiltemp.jl / gpu_validate_soilwater.jl:
# run each @kernel on a CPU Array and on a device array of the SAME inputs at the
# SAME precision and compare. Parity (does the device execute the kernel the way
# the CPU does?), not accuracy-vs-Fortran.
#
#   julia --project=scripts scripts/gpu_validate_lake.jl
#
# Kernels covered (grows per increment):
#   lt1: _lake_soil_tk_kernel! / _lake_tk_interface_kernel! / _lake_soil_cv_kernel! /
#        _lake_snow_cv_kernel! / _lake_total_h2osno_kernel!
# ==========================================================================

using CLM
using Printf
using Random
include(joinpath(@__DIR__, "gpu_backends.jl"))

function reldiff(a, b)
    A = Array(a); B = Array(b)
    maximum(abs.(A .- B) ./ (one(eltype(A)) .+ max.(abs.(A), abs.(B))); init = 0.0)
end
TOL(::Type{Float32}) = 1f-2
TOL(::Type{T}) where {T} = 1e-9

const NC = 8
const NS = 5          # nlevsno
const NG = 10         # nlevgrnd
const NSOI = 10       # nlevsoi
const NLAK = 10       # nlevlak
const NLEV = NS + NG  # layer-array 2nd dim (15)
const HFUS = CLM.HFUS; const TFRZ = CLM.TFRZ

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

# Branch-exercising inputs at precision FT. snl spans 0..-3 so snow layers fire;
# t_soisno spans TFRZ so frozen/unfrozen soil-tk branches both run.
function build_inputs(::Type{FT}) where {FT}
    rng = MersenneTwister(23)
    rnd(d...) = FT.(rand(rng, d...))
    mask = trues(NC)
    snl  = Int[0, -1, -2, -3, -1, -2, 0, -3]
    # strictly increasing depths so interface denominators are well-posed
    z  = FT[0.1 * jj for c in 1:NC, jj in 1:NLEV]
    zi = FT[0.1 * jj - 0.05 for c in 1:NC, jj in 1:(NLEV + 1)]
    dz = fill(FT(0.1), NC, NLEV)
    t_soisno = FT[273.15 + 1.5 * sinpi((c + jj) / 5) for c in 1:NC, jj in 1:NLEV]
    return (; mask, snl, z, zi, dz, t_soisno,
        watsat = FT(0.4) .+ rnd(NC, NG) .* FT(0.1),
        tksatu = FT(1.5) .+ rnd(NC, NG),
        tkmg   = FT(2.0) .+ rnd(NC, NG),
        tkdry  = FT(0.25) .+ rnd(NC, NG) .* FT(0.1),
        csol   = FT(2.0e6) .+ rnd(NC, NG) .* FT(1.0e5),
        h2osoi_ice = FT(2.0) .+ rnd(NC, NLEV) .* 5,
        h2osoi_liq = FT(2.0) .+ rnd(NC, NLEV) .* 5,
        h2osno_no_layers = FT[0.0, 2.0, 0.0, 5.0, 1.0, 0.0, 3.0, 0.0])
end

function tests(to, ::Type{FT}) where {FT}
    I = build_inputs(FT)
    results = Tuple{String,Float64,Bool}[]

    thk0 = zeros(FT, NC, NLEV)
    push!(results, _parity("soil_tk", to, FT, CLM._lake_soil_tk_kernel!, (NC, NLEV), thk0,
        f -> ((f(I.mask), f(I.snl), f(I.dz), f(I.t_soisno), f(I.h2osoi_liq), f(I.h2osoi_ice),
               f(I.watsat), f(I.tksatu), f(I.tkmg), f(I.tkdry), NS, NSOI), ())))

    # interface tk needs a populated thk: run the CPU soil-tk once to get it.
    thk_cpu = zeros(FT, NC, NLEV)
    CLM._launch!(CLM._lake_soil_tk_kernel!, thk_cpu, I.mask, I.snl, I.dz, I.t_soisno,
        I.h2osoi_liq, I.h2osoi_ice, I.watsat, I.tksatu, I.tkmg, I.tkdry, NS, NSOI;
        ndrange = (NC, NLEV))
    tk0 = zeros(FT, NC, NLEV)
    push!(results, _parity("tk_interface", to, FT, CLM._lake_tk_interface_kernel!, (NC, NLEV), tk0,
        f -> begin
            ttop = f(zeros(FT, NC))
            ((ttop, f(thk_cpu), f(I.mask), f(I.snl), f(I.z), f(I.zi), NS, NG), (ttop,))
        end))

    cv0 = zeros(FT, NC, NLEV)
    push!(results, _parity("soil_cv", to, FT, CLM._lake_soil_cv_kernel!, (NC, NG), cv0,
        f -> ((f(I.mask), f(I.csol), f(I.watsat), f(I.dz), f(I.h2osoi_ice), f(I.h2osoi_liq),
               NS, NSOI), ())))
    push!(results, _parity("snow_cv", to, FT, CLM._lake_snow_cv_kernel!, (NC, NS), cv0,
        f -> ((f(I.mask), f(I.snl), f(I.h2osoi_liq), f(I.h2osoi_ice), NS), ())))

    h0 = zeros(FT, NC)
    push!(results, _parity("total_h2osno", to, FT, CLM._lake_total_h2osno_kernel!, NC, h0,
        f -> ((f(I.mask), f(I.snl), f(I.h2osno_no_layers), f(I.h2osoi_ice), f(I.h2osoi_liq),
               NS), ())))

    push!(results, test_phase_change(to, FT))

    # ---- lt3: lake density + diffusivity (the Metal-risk ones: smooth_abs, exp, pow) ----
    rng2 = MersenneTwister(37)
    rnd2(d...) = FT.(rand(rng2, d...))
    icefrac = rnd2(NC, NLAK) .* FT(0.5)
    t_lake  = FT[272.0 + 4.0 * sinpi((c + j) / 5) for c in 1:NC, j in 1:NLAK]
    z_lake  = FT[0.5 * j for c in 1:NC, j in 1:NLAK]
    rhow_in = FT(999.0) .+ rnd2(NC, NLAK) .* 2
    ws_in = FT(0.05) .+ rnd2(NC); ks_in = FT(0.1) .+ rnd2(NC) .* FT(0.1)
    t_grnd = FT[274.0, 271.0, 275.0, 272.0, 274.0, 270.0, 276.0, 273.0]
    snl = Int[0, 0, 0, -1, 0, -2, 0, 0]
    lakedepth = FT[2.0, 60.0, 5.0, 80.0, 3.0, 50.0, 10.0, 100.0]
    cwat = FT(CLM.CPLIQ * CLM.DENH2O); tkice_eff = FT(CLM.TKICE * CLM.DENICE / CLM.DENH2O)
    km = FT(CLM.TKWAT) / cwat; p0 = FT(1.0)

    rho0 = zeros(FT, NC, NLAK)
    push!(results, _parity("lake_density", to, FT, CLM._lake_density_kernel!, (NC, NLAK), rho0,
        f -> ((f(trues(NC)), f(icefrac), f(t_lake), NLAK), ())))

    kme0 = zeros(FT, NC, NLAK)
    push!(results, _parity("diffusivity", to, FT, CLM._lake_diffusivity_kernel!,
        (NC, NLAK - 1), kme0,
        f -> begin
            tkl = f(zeros(FT, NC, NLAK)); frzn = f(fill(false, NC))
            ((tkl, frzn, f(trues(NC)), f(rhow_in), f(z_lake), f(ws_in), f(ks_in), f(t_grnd),
              f(t_lake), f(snl), f(lakedepth), f(icefrac), NLAK, km, p0, cwat, tkice_eff),
             (tkl,))
        end))
    return results
end

# phase_change_lake!: run BOTH per-column kernels (water then soil) on CPU and on
# device with identical inputs, compare every mutated array. Inputs trigger melt
# (t>TFRZ + ice) and freeze (t<TFRZ + liquid) on both lake and snow/soil sides.
function test_phase_change(to, ::Type{FT}) where {FT}
    rng = MersenneTwister(31)
    rnd(d...) = FT.(rand(rng, d...))
    mask = trues(NC)
    snl  = Int[0, -1, -2, -3, -1, -2, 0, -3]
    dz_lake = fill(FT(1.0), NC, NLAK)
    # alternate warm/cold lake layers; ice present so both branches fire
    t_lake0 = FT[(c + j) % 2 == 0 ? 274.0 : 271.0 for c in 1:NC, j in 1:NLAK]
    icefrac0 = fill(FT(0.5), NC, NLAK)
    cvlake0 = FT(4.0e6) .+ rnd(NC, NLAK) .* FT(1.0e5)
    snowd0 = FT(0.1) .+ rnd(NC)
    h2osno0 = FT[0.0, 2.0, 0.0, 5.0, 1.0, 0.0, 3.0, 0.0]
    t_soisno0 = FT[273.15 + 1.5 * sinpi((c + jj) / 4) for c in 1:NC, jj in 1:NLEV]
    ice0 = FT(2.0) .+ rnd(NC, NLEV) .* 3
    liq0 = FT(2.0) .+ rnd(NC, NLEV) .* 3
    cv0  = FT(2.0e6) .+ rnd(NC, NLEV) .* FT(1.0e5)
    dt = FT(1800)

    # mutable set: copies of every read/written array, on backend `cv`(=identity/device)
    mkset(cv) = (
        t_lake = cv(copy(t_lake0)), snow_depth = cv(copy(snowd0)),
        h2osno = cv(copy(h2osno0)), cv_lake = cv(copy(cvlake0)),
        icefrac = cv(copy(icefrac0)), qsm = cv(zeros(FT, NC)), qsf = cv(zeros(FT, NC)),
        qsd = cv(zeros(FT, NC)), lhabs = cv(zeros(FT, NC)),
        qsm_lyr = cv(zeros(FT, NC, NLEV)), qsf_lyr = cv(zeros(FT, NC, NLEV)),
        imelt = cv(zeros(Int, NC, NLEV)), t_soisno = cv(copy(t_soisno0)),
        ice = cv(copy(ice0)), liq = cv(copy(liq0)), cv = cv(copy(cv0)),
        eflx = cv(zeros(FT, NC)))
    run!(s, to_) = begin
        m = _dev(to_, mask); sn = _dev(to_, snl); dzl = _dev(to_, dz_lake)
        CLM._launch!(CLM._phase_change_lake_water_kernel!, s.t_lake, s.snow_depth, s.h2osno,
            s.cv_lake, s.icefrac, s.qsm, s.qsf, s.qsd, s.lhabs, s.qsm_lyr, s.qsf_lyr,
            s.imelt, m, dzl, NS, NLAK, dt; ndrange = NC)
        CLM._launch!(CLM._phase_change_lake_soil_kernel!, s.t_soisno, s.ice, s.liq, s.cv,
            s.lhabs, s.qsm, s.qsf, s.eflx, s.imelt, s.qsm_lyr, s.qsf_lyr, m, sn, NS, NG, dt;
            ndrange = NC)
        s
    end
    cpu = run!(mkset(identity), identity)
    dev = run!(mkset(x -> _dev(to, x)), to)
    flds = (:t_lake, :snow_depth, :h2osno, :cv_lake, :icefrac, :qsm, :qsf, :qsd, :lhabs,
            :qsm_lyr, :qsf_lyr, :imelt, :t_soisno, :ice, :liq, :cv, :eflx)
    d = maximum(reldiff(getfield(cpu, f), getfield(dev, f)) for f in flds)
    return ("phase_change", d, d < TOL(FT))
end

function main(backend)
    println("=" ^ 70)
    println("METAL PARITY for lake_temperature! kernels")
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
        @printf("  [%s] %-14s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", knm, d)
        ok ? (npass += 1) : (nfail += 1)
    end
    @printf("\n  %d passed, %d failed\n", npass, nfail)
    println(nfail == 0 ? "\n  ALL lake kernels MATCH CPU ON $name ($FT) ✓" :
                          "\n  SOME KERNELS DIVERGED — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
