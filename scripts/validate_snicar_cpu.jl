# ==========================================================================
# validate_snicar_cpu.jl — CPU-backend parity check of snicar_rt_device!
# against the host snicar_rt!. Both run on identical Float64 CPU arrays (the
# device wrapper uses the KernelAbstractions CPU backend), so any divergence is
# a port bug, not a precision artifact. This de-risks the kernel BEFORE Metal.
#
#   julia --project=. scripts/validate_snicar_cpu.jl
# ==========================================================================
using CLM, Printf

const FSNOWOPTICS = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_c013122.nc"

CLM.snow_optics_init!(CLM.snicar_optics; fsnowoptics=FSNOWOPTICS)

const nlevsno = 12
const NR = CLM.NUMRAD
const AER = CLM.SNO_NBR_AER

# Synthetic columns covering every branch:
#   1: snl=-3, sun, 3 snow layers, aerosols           (multi-layer RT)
#   2: snl=-1, low sun (mu_not<0.2588 → SZA adjust)    (1 layer + SZA fix)
#   3: snl= 0, sun, snow mass > MIN_SNW (flg_nosnl=1)   (lumped layer)
#   4: snl= 0, sun, 0 < snow < MIN_SNW                  (surface-albedo branch)
#   5: snl=-2, no sun (coszen=0)                        (zero branch)
const nc = 5
snl   = Int[-3, -1, 0, 0, -2]
coszen = Float64[0.6, 0.15, 0.5, 0.5, 0.0]

h2osno_liq = zeros(nc, nlevsno)
h2osno_ice = zeros(nc, nlevsno)
snw_rds    = fill(round(Int, CLM.snicar_params.snw_rds_min), nc, nlevsno)
mss_cnc    = zeros(nc, nlevsno, AER)
albsfc     = zeros(nc, NR)
h2osno_tot = zeros(nc)
frac_sno   = fill(0.9, nc)

for c in 1:nc
    albsfc[c, CLM.IVIS] = 0.18 + 0.02c
    albsfc[c, CLM.INIR] = 0.34 + 0.02c
end

# snow-bearing layered columns (1,2,5)
function fill_layers!(c, s)
    layers = (nlevsno + s + 1):nlevsno
    tot = 0.0
    for (k, L) in enumerate(layers)
        h2osno_ice[c, L] = 18.0 + 6.0k
        h2osno_liq[c, L] = 0.8 + 0.2k
        snw_rds[c, L]    = 100 + 25k
        # aerosols (BC/OC/dust) — small, exercise the aerosol-optics path
        mss_cnc[c, L, 1] = 5.0e-9
        mss_cnc[c, L, 2] = 3.0e-9
        mss_cnc[c, L, 3] = 2.0e-9
        mss_cnc[c, L, 5] = 1.0e-8
        mss_cnc[c, L, 6] = 5.0e-9
        tot += h2osno_ice[c, L] + h2osno_liq[c, L]
    end
    h2osno_tot[c] = tot
end
fill_layers!(1, -3)
fill_layers!(2, -1)
fill_layers!(5, -2)         # no sun, but give it mass anyway
h2osno_tot[3] = 30.0        # flg_nosnl=1 lumped path
h2osno_tot[4] = 1.0e-40     # tiny snow, surface-albedo branch

reldiff(a, b) = begin
    A = Array(a); B = Array(b); m = 0.0; n = 0
    for i in eachindex(A, B)
        (isfinite(A[i]) && isfinite(B[i])) || continue
        n += 1
        m = max(m, abs(A[i] - B[i]) / (1.0 + max(abs(A[i]), abs(B[i]))))
    end
    (m, n)
end

function run_case(flg)
    aH = zeros(nc, NR); fH = zeros(nc, nlevsno + 1, NR)
    aD = zeros(nc, NR); fD = zeros(nc, nlevsno + 1, NR)
    CLM.snicar_rt!(copy(coszen), flg, copy(h2osno_liq), copy(h2osno_ice), copy(h2osno_tot),
                   copy(snw_rds), copy(mss_cnc), copy(albsfc), copy(snl), copy(frac_sno),
                   aH, fH, nlevsno; mask_nourbanc=trues(nc))
    CLM.snicar_rt_device!(copy(coszen), flg, copy(h2osno_liq), copy(h2osno_ice), copy(h2osno_tot),
                          copy(snw_rds), copy(mss_cnc), copy(albsfc), copy(snl), copy(frac_sno),
                          aD, fD, nlevsno; mask=nothing)
    da, na = reldiff(aH, aD)
    df, nf = reldiff(fH, fD)
    @printf("  flg_slr=%d  albout rel=%.3e (%d)   flx_abs rel=%.3e (%d)\n", flg, da, na, df, nf)
    return aH, aD, fH, fD, max(da, df)
end

println("SNICAR CPU-backend parity: host snicar_rt! vs snicar_rt_device!")
aH1, aD1, fH1, fD1, m1 = run_case(1)
aH2, aD2, fH2, fD2, m2 = run_case(2)
gmax = max(m1, m2)
@printf("\n  max rel over both beams = %.3e — %s\n", gmax,
        gmax < 1e-10 ? "MATCH" : "DIVERGENCE")

if gmax >= 1e-10
    println("\n  per-column albout (host | device), flg_slr=1:")
    for c in 1:nc
        @printf("   col %d snl=%2d cosz=%.2f  H=[%.6f %.6f]  D=[%.6f %.6f]\n",
                c, snl[c], coszen[c], aH1[c,1], aH1[c,2], aD1[c,1], aD1[c,2])
    end
end
exit(gmax < 1e-10 ? 0 : 1)
