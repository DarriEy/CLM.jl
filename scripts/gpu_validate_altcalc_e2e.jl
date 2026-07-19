# ==========================================================================
# gpu_validate_altcalc_e2e.jl — end-to-end GPU parity for the WHOLE alt_calc!
# active-layer-thickness driver (annual-maxima reset kernel on the set annual
# timestep + the per-column ALT search/interpolation kernel).
#
# Builds a multi-column ActiveLayerData/TemperatureData/ColumnData/GridcellData
# at Float32 exercising all branches:
#   col 1 — fully thawed (deepest layer > TFRZ)            -> alt = zsoi[nlevgrnd]
#   col 2 — fully frozen (all layers <= TFRZ)              -> alt = 0
#   col 3 — partial thaw (warm shallow, frozen deep)       -> linear interp branch
# and both hemispheres (NH lat>0 / SH lat<=0) for the Jan-1 reset.
#
# Runs alt_calc! on the CPU, adapts the structs to the GPU, runs the SAME call on
# the device, and compares alt / alt_indx / altmax / altmax_indx /
# altmax_lastyear / altmax_lastyear_indx with reldiff.
#
# CRITICAL: t_soisno is set to real temperatures so the CPU reference is FINITE;
# we assert finiteness before trusting parity.
#
#   julia --project=scripts scripts/gpu_validate_altcalc_e2e.jl
# ==========================================================================

using CLM
using Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))

function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(A[i] - B[i]) / (1.0 + max(abs(A[i]), abs(B[i]))))
    end
    return m
end

function setup_params!()
    vp = CLM.varpar
    vp.nlevsno = 5
    vp.nlevsoi = 10
    vp.nlevgrnd = 15
    vp.nlevlak = 10
    vp.nlevurb = 5
    vp.nlevmaxurbgrnd = max(vp.nlevurb, vp.nlevgrnd)
    CLM.varcon_init!()
    return vp
end

function build(::Type{FT}) where {FT}
    setup_params!()
    nlevgrnd = CLM.varpar.nlevgrnd
    nlevsno  = CLM.varpar.nlevsno
    joff = nlevsno
    nc = 3; ng = 2

    al = CLM.ActiveLayerData{FT}()
    CLM.active_layer_init!(al, nc)
    for c in 1:nc
        al.alt_col[c]                  = FT(0.0)
        al.altmax_col[c]               = FT(0.3)   # nonzero prior max to test reset/update
        al.altmax_lastyear_col[c]      = FT(0.1)
        al.alt_indx_col[c]             = 0
        al.altmax_indx_col[c]          = 5
        al.altmax_lastyear_indx_col[c] = 2
    end

    temperature = CLM.TemperatureData{FT}()
    CLM.temperature_init!(temperature, nc, nc, 1, ng)
    # initialize entire profile cold, then set per-column branch temperatures
    for c in 1:nc, j in 1:size(temperature.t_soisno_col, 2)
        temperature.t_soisno_col[c, j] = FT(260.0)
    end
    # col 1: fully thawed (all soil > TFRZ)
    for j in 1:nlevgrnd
        temperature.t_soisno_col[1, j + joff] = FT(280.0)
    end
    # col 2: fully frozen (all soil <= TFRZ) — already 260
    # col 3: warm shallow (j=1..6 > TFRZ), frozen deeper -> interp at the crossing
    for j in 1:6
        temperature.t_soisno_col[3, j + joff] = FT(278.0)
    end
    for j in 7:nlevgrnd
        temperature.t_soisno_col[3, j + joff] = FT(268.0)
    end

    col_data = CLM.ColumnData{FT}()
    CLM.column_init!(col_data, nc)
    col_data.gridcell[1] = 1   # NH
    col_data.gridcell[2] = 2   # SH
    col_data.gridcell[3] = 1   # NH

    grc = CLM.GridcellData{FT}()
    CLM.gridcell_init!(grc, ng)
    grc.lat[1] = FT(0.8)    # NH
    grc.lat[2] = FT(-0.8)   # SH

    mask = trues(nc)
    return (; nc, al, temperature, col_data, grc, mask)
end

# Use the Jan-1 annual-reset timestep (mon=1,day=1, div(sec,dtime)==1) so the
# reset kernel (NH branch) is exercised in addition to the ALT search kernel.
run_alt!(H) = CLM.alt_calc!(H.al, H.mask, H.temperature, H.col_data, H.grc;
                            mon = 1, day = 1, sec = 1800, dtime = 1800)

function main(backend)
    println("=" ^ 70)
    println("END-TO-END GPU parity for alt_calc! (whole driver)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU driver exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    H = build(FT)
    B = build(FT)

    ad(x) = CLM.Adapt.adapt(device_array_type(), x)
    al_d   = ad(B.al)
    temp_d = ad(B.temperature)
    col_d  = ad(B.col_data)
    grc_d  = ad(B.grc)
    # The whole alt_calc! runs on the device, so its mask must be a device Bool
    # vector (BitVector can't reach a GPU kernel); convert BitArray -> Bool.
    mask_d = to(collect(Bool, B.mask))

    if !(al_d.alt_col isa device_array_type())
        println("  BLOCKED: ActiveLayerData did not move to the device under adapt.")
        return 2
    end

    run_alt!(H)
    CLM.alt_calc!(al_d, mask_d, temp_d, col_d, grc_d;
                  mon = 1, day = 1, sec = 1800, dtime = 1800)

    # Guard: CPU reference alt/altmax must be FINITE.
    if !(all(isfinite, Array(H.al.alt_col)) && all(isfinite, Array(H.al.altmax_col)))
        println("  BLOCKED: CPU alt/altmax not finite — parity would be a false PASS.")
        return 2
    end

    checks = [
        ("alt",                  H.al.alt_col,                  al_d.alt_col),
        ("alt_indx",             H.al.alt_indx_col,             al_d.alt_indx_col),
        ("altmax",               H.al.altmax_col,               al_d.altmax_col),
        ("altmax_indx",          H.al.altmax_indx_col,          al_d.altmax_indx_col),
        ("altmax_lastyear",      H.al.altmax_lastyear_col,      al_d.altmax_lastyear_col),
        ("altmax_lastyear_indx", H.al.altmax_lastyear_indx_col, al_d.altmax_lastyear_indx_col),
    ]
    nfail = 0; worst = 0.0
    for (nm, a, b) in checks
        d = reldiff(a, b); ok = d < 1f-2
        worst = max(worst, d)
        @printf("  [%s] %-22s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end
    @printf("\n  worst reldiff = %.3e\n", worst)
    println(nfail == 0 ? "  WHOLE alt_calc! MATCHES CPU ON $name ($FT)" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
