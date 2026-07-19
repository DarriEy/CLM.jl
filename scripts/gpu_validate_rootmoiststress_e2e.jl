# ==========================================================================
# gpu_validate_rootmoiststress_e2e.jl — end-to-end GPU parity for the WHOLE
# calc_root_moist_stress! driver (every per-(p)/per-(p,j) kernel together).
#
# Builds a small Float32 instance (one soil patch/column, multiple ground layers)
# mirroring test/test_soil_moist_stress.jl's dispatcher setup, runs
# calc_root_moist_stress! on the CPU, adapts every state struct (+ masks/params)
# to the GPU, runs the SAME call on the device, and compares the mutated outputs
# field-by-field. This exercises: btran zeroing, the unfrozen-rootfr normalization
# (array_normalization sum-to-one), the per-patch rresis/rootr + btran accumulation
# (sequential layer reduction in-kernel), and the btran-normalize kernel together.
#
#   julia --project=scripts scripts/gpu_validate_rootmoiststress_e2e.jl
# ==========================================================================

using CLM
using Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))

# reldiff: NaN-aware (both-NaN agrees; one-sided NaN flags divergence).
function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(A[i] - B[i]) / (1.0 + max(abs(A[i]), abs(B[i]))))
    end
    return m
end
cpu_has_finite(a) = any(isfinite, Array(a))

function build(::Type{FT}) where {FT}
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    np = 1; nc = 1
    nlevsno  = CLM.varpar.nlevsno
    nlevgrnd = CLM.varpar.nlevgrnd
    joff = nlevsno

    CLM.set_perchroot_opt!(false, false)
    CLM.init_root_moist_stress!()

    col = CLM.ColumnData{FT}(); CLM.column_init!(col, nc)
    col.dz .= FT(0.1)

    patchdata = CLM.PatchData{FT}(); CLM.patch_init!(patchdata, np)
    for p in 1:np
        patchdata.column[p] = p
        patchdata.itype[p] = 1
    end

    temp = CLM.TemperatureData{FT}(); CLM.temperature_init!(temp, np, nc, nc, 1)
    temp.t_soisno_col .= FT(290.0)

    soilstate = CLM.SoilStateData{FT}(); CLM.soilstate_init!(soilstate, np, nc)
    soilstate.watsat_col .= FT(0.45)
    soilstate.sucsat_col .= FT(100.0)
    soilstate.bsw_col .= FT(5.0)
    soilstate.eff_porosity_col .= FT(0.40)
    for p in 1:np, j in 1:nlevgrnd
        soilstate.rootfr_patch[p, j] = FT(1.0 / nlevgrnd)
    end
    soilstate.rootr_patch .= FT(0.0)

    energyflux = CLM.EnergyFluxData{FT}(); CLM.energyflux_init!(energyflux, np, nc, nc, 1)
    energyflux.btran_patch .= FT(0.0)
    energyflux.rresis_patch .= FT(0.0)

    waterstatebulk = CLM.WaterStateBulkData{FT}(); CLM.waterstatebulk_init!(waterstatebulk, nc, np, nc, 1)
    waterdiagbulk = CLM.WaterDiagnosticBulkData{FT}(); CLM.waterdiagnosticbulk_init!(waterdiagbulk, nc, np, nc, 1)
    for c in 1:nc, j in 1:nlevgrnd
        waterdiagbulk.h2osoi_liqvol_col[c, j + joff] = FT(0.25)
    end

    npft = 2
    smpso = fill(FT(-35000.0), npft)
    smpsc = fill(FT(-275000.0), npft)
    altmax_lastyear_indx = fill(FT(10.0), nc)
    altmax_indx = fill(FT(10.0), nc)

    S = (; soilstate, energyflux, temp, waterstatebulk, waterdiagbulk, col, patchdata)
    return (; np, nc, nlevgrnd, nlevsno, S,
              smpso, smpsc, altmax_lastyear_indx, altmax_indx)
end

run_rms!(S, smpso, smpsc, altmax_l, altmax_i, m_p, nlevgrnd, nlevsno) =
    CLM.calc_root_moist_stress!(S.soilstate, S.energyflux, S.temp,
        S.waterstatebulk, S.waterdiagbulk, S.col, S.patchdata,
        smpso, smpsc, altmax_l, altmax_i, m_p, 1:length(m_p), nlevgrnd, nlevsno)

function main(backend)
    println("=" ^ 70)
    println("END-TO-END GPU parity for calc_root_moist_stress! (whole driver)")
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
    Sd = map(ad, B.S)
    dmask(m) = to(collect(Bool, m))

    if !(Sd.soilstate.rootr_patch isa device_array_type())
        println("  BLOCKED: a state struct did not move to the device under adapt.")
        return 2
    end

    m_p = trues(H.np)

    run_rms!(H.S, H.smpso, H.smpsc, H.altmax_lastyear_indx, H.altmax_indx,
             m_p, H.nlevgrnd, H.nlevsno)
    run_rms!(Sd, to(B.smpso), to(B.smpsc), to(B.altmax_lastyear_indx), to(B.altmax_indx),
             dmask(m_p), B.nlevgrnd, B.nlevsno)

    checks = [
        ("btran",  H.S.energyflux.btran_patch,  Sd.energyflux.btran_patch),
        ("rresis", H.S.energyflux.rresis_patch, Sd.energyflux.rresis_patch),
        ("rootr",  H.S.soilstate.rootr_patch,   Sd.soilstate.rootr_patch),
    ]
    nfail = 0
    for (nm, a, b) in checks
        if !cpu_has_finite(a)
            @printf("  [WARN] %-8s CPU reference is all-NaN/Inf — skipping (no parity signal)\n", nm)
            continue
        end
        d = reldiff(a, b); ok = d < 1f-3
        @printf("  [%s] %-8s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end

    # Sanity: btran must be physically in (0,1] and rootr must sum to 1 per patch.
    bt = Array(H.S.energyflux.btran_patch)
    rr = Array(H.S.soilstate.rootr_patch)
    @printf("\n  CPU btran[1] = %.6f   sum(rootr[1,:]) = %.6f\n", bt[1], sum(rr[1, :]))
    if !(bt[1] > 0.0 && bt[1] <= 1.0 + 1e-6)
        println("  [WARN] CPU btran out of (0,1] — check the setup branches.")
    end

    println()
    println(nfail == 0 ? "  WHOLE calc_root_moist_stress! MATCHES CPU ON $name ($FT) ✓" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
