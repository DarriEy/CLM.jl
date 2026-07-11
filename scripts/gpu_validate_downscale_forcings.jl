# ==========================================================================
# gpu_validate_downscale_forcings.jl — parity for the kernelized atmospheric
# forcing downscaling (atm2lndMod.F90 → downscale_forcings.jl).
#
# Proves TWO things:
#   1. BYTE-IDENTITY of the default host path: the KA.CPU kernels (Float64)
#      reproduce the ORIGINAL explicit-loop reference bit-for-bit.
#   2. FLOAT32 PARITY on Metal: the same kernels on device match the host
#      Float64 result within Float32 tolerance.
#
#   julia --project=scripts scripts/gpu_validate_downscale_forcings.jl
# ==========================================================================
using CLM
using Printf
import Metal
include(joinpath(@__DIR__, "gpu_backends.jl"))
include(joinpath(@__DIR__, "gpu_adapt.jl"))
mf(x) = mf(Metal.MtlArray, x)

const NR = CLM.NUMRAD

function reldiff(H, D)
    A = Array(H); B = Array(D); m = 0.0; n = 0
    for i in eachindex(A, B)
        (isfinite(A[i]) && isfinite(B[i])) || continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + abs(Float64(A[i])))); n += 1
    end
    return m, n
end

maxabs(H, R) = (A = Array(H); B = Array(R); maximum(abs.(Float64.(A) .- Float64.(B)); init = 0.0))

# --------------------------------------------------------------------------
# Fixture: 3 gridcells / 3 columns / 3 landunits (1:1:1), with elevation
# offsets and one glacier (ISTICE) landunit. Repartition + LW downscaling on.
# Built fresh each call so reference / host / device runs are independent.
# --------------------------------------------------------------------------
function build()
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    ng = 3; nc = 3; nl = 3; np = 3
    ISTSOIL = CLM.ISTSOIL; ISTICE = CLM.ISTICE

    a2l = CLM.Atm2LndData{Float64}()
    CLM.atm2lnd_init!(a2l, ng, nc, np)
    CLM.atm2lnd_params_init!(a2l.params;
        repartition_rain_snow = true,
        glcmec_downscale_longwave = true,
        lapse_rate = 0.006,
        lapse_rate_longwave = 0.032,
        longwave_downscaling_limit = 0.5,
        precip_repartition_glc_all_snow_t = -5.0,
        precip_repartition_glc_all_rain_t = 5.0,
        precip_repartition_nonglc_all_snow_t = -5.0,
        precip_repartition_nonglc_all_rain_t = 5.0)

    # gridcell not-downscaled forcing
    a2l.forc_topo_grc              .= [0.0, 100.0, 500.0]
    a2l.forc_hgt_grc               .= [30.0, 30.0, 30.0]
    a2l.forc_t_not_downscaled_grc  .= [274.0, 272.0, 278.0]
    a2l.forc_th_not_downscaled_grc .= [274.5, 272.5, 278.5]
    a2l.forc_pbot_not_downscaled_grc .= [1.0e5, 0.98e5, 0.95e5]
    a2l.forc_rho_not_downscaled_grc  .= [1.20, 1.22, 1.15]
    a2l.forc_lwrad_not_downscaled_grc .= [300.0, 290.0, 310.0]
    a2l.forc_solar_not_downscaled_grc .= [400.0, 350.0, 500.0]
    a2l.forc_vp_grc                 .= [900.0, 800.0, 1100.0]
    a2l.forc_pco2_grc               .= [37.0, 36.3, 35.2]
    a2l.forc_po2_grc                .= [2.1e4, 2.06e4, 2.0e4]
    a2l.forc_rain_not_downscaled_grc .= [1.0e-4, 0.5e-4, 2.0e-4]
    a2l.forc_snow_not_downscaled_grc .= [0.5e-4, 2.0e-4, 0.2e-4]
    for b in 1:NR
        a2l.forc_solad_not_downscaled_grc[:, b] .= [200.0 + 10b, 150.0 + 10b, 250.0 + 10b]
    end

    col = CLM.ColumnData{Float64}()
    CLM.column_init!(col, nc)
    col.gridcell .= [1, 2, 3]
    col.landunit .= [1, 2, 3]

    lun = CLM.LandunitData{Float64}()
    CLM.landunit_init!(lun, nl)
    lun.itype .= [ISTSOIL, ISTICE, ISTSOIL]

    topo = CLM.TopoData{Float64}()
    CLM.topo_init_allocate!(topo, nc)
    topo.topo_col .= [200.0, 100.0, 800.0]   # elevation offsets from gridcell

    bounds = CLM.BoundsType(begg = 1, endg = ng, begl = 1, endl = nl,
                            begc = 1, endc = nc, begp = 1, endp = np)
    return (; a2l, col, lun, topo, bounds)
end

# --------------------------------------------------------------------------
# Reference: the ORIGINAL explicit-loop implementation (verbatim from the
# pre-kernelization downscale_forcings.jl), Float64. Used only to prove the
# host KA path is byte-identical.
# --------------------------------------------------------------------------
function reference!(f)
    a2l = f.a2l; col = f.col; lun = f.lun; topo = f.topo; bounds = f.bounds
    RAIR = CLM.RAIR; GRAV = CLM.GRAV; CPAIR = CLM.CPAIR; ISTICE = CLM.ISTICE
    bc_col = bounds.begc:bounds.endc
    lapse_rate = a2l.params.lapse_rate
    isnan(lapse_rate) && (lapse_rate = 0.006)

    for c in bc_col
        g = col.gridcell[c]
        a2l.forc_t_downscaled_col[c]     = a2l.forc_t_not_downscaled_grc[g]
        a2l.forc_th_downscaled_col[c]    = a2l.forc_th_not_downscaled_grc[g]
        a2l.forc_pbot_downscaled_col[c]  = a2l.forc_pbot_not_downscaled_grc[g]
        a2l.forc_rho_downscaled_col[c]   = a2l.forc_rho_not_downscaled_grc[g]
        a2l.forc_lwrad_downscaled_col[c] = a2l.forc_lwrad_not_downscaled_grc[g]
        for b in 1:CLM.NUMRAD
            a2l.forc_solad_downscaled_col[c, b] = a2l.forc_solad_not_downscaled_grc[g, b]
        end
        a2l.forc_solar_downscaled_col[c] = a2l.forc_solar_not_downscaled_grc[g]
    end
    for c in bc_col
        g = col.gridcell[c]
        hsurf_g = a2l.forc_topo_grc[g]; hsurf_c = topo.topo_col[c]
        (isnan(hsurf_c) || isnan(hsurf_g)) && continue
        tbot_g = a2l.forc_t_not_downscaled_grc[g]
        tbot_c = tbot_g - lapse_rate * (hsurf_c - hsurf_g)
        a2l.forc_t_downscaled_col[c] = tbot_c
        zbot = a2l.forc_hgt_grc[g]
        Hbot = RAIR * 0.5 * (tbot_g + tbot_c) / GRAV
        if Hbot > 0.0
            a2l.forc_th_downscaled_col[c] = a2l.forc_th_not_downscaled_grc[g] +
                (tbot_c - tbot_g) * exp((zbot / Hbot) * (RAIR / CPAIR))
        end
        if Hbot > 0.0
            a2l.forc_pbot_downscaled_col[c] = a2l.forc_pbot_not_downscaled_grc[g] *
                exp(-(hsurf_c - hsurf_g) / Hbot)
        end
        pbot_g = a2l.forc_pbot_not_downscaled_grc[g]; pbot_c = a2l.forc_pbot_downscaled_col[c]
        rho_est_g = pbot_g / (RAIR * tbot_g); rho_est_c = pbot_c / (RAIR * tbot_c)
        if rho_est_g > 0.0
            a2l.forc_rho_downscaled_col[c] = a2l.forc_rho_not_downscaled_grc[g] *
                (rho_est_c / rho_est_g)
        end
    end
    for c in bc_col
        g = col.gridcell[c]
        pbot_nd = a2l.forc_pbot_not_downscaled_grc[g]; pbot_ds = a2l.forc_pbot_downscaled_col[c]
        if pbot_nd > 0.0 && pbot_ds > 0.0
            a2l.forc_pco2_grc[g] = (a2l.forc_pco2_grc[g] / pbot_nd) * pbot_ds
            a2l.forc_po2_grc[g]  = (a2l.forc_po2_grc[g]  / pbot_nd) * pbot_ds
        end
    end
    # partition_precip
    for c in bc_col
        g = col.gridcell[c]
        total_rain = a2l.forc_rain_not_downscaled_grc[g]
        total_snow = a2l.forc_snow_not_downscaled_grc[g]
        total_precip = total_rain + total_snow
        t_col = a2l.forc_t_downscaled_col[c]
        if a2l.params.repartition_rain_snow && total_precip > 0.0
            l = col.landunit[c]; is_glc = lun.itype[l] == ISTICE
            all_snow_t      = is_glc ? a2l.params.precip_repartition_glc_all_snow_t :
                                       a2l.params.precip_repartition_nonglc_all_snow_t
            frac_rain_slope = is_glc ? a2l.params.precip_repartition_glc_frac_rain_slope :
                                       a2l.params.precip_repartition_nonglc_frac_rain_slope
            frac_rain = clamp((t_col - all_snow_t) * frac_rain_slope, 0.0, 1.0)
            a2l.forc_rain_downscaled_col[c] = total_precip * frac_rain
            a2l.forc_snow_downscaled_col[c] = total_precip * (1.0 - frac_rain)
        else
            a2l.forc_rain_downscaled_col[c] = total_rain
            a2l.forc_snow_downscaled_col[c] = total_snow
        end
    end
    for c in bc_col
        g = col.gridcell[c]
        vp = a2l.forc_vp_grc[g]; pbot = a2l.forc_pbot_downscaled_col[c]
        a2l.forc_q_downscaled_col[c] = 0.622 * vp / max(pbot - 0.378 * vp, 1.0)
    end
    # downscale_longwave
    lr_lw = a2l.params.lapse_rate_longwave; limit = a2l.params.longwave_downscaling_limit
    for c in bc_col
        g = col.gridcell[c]
        hsurf_g = a2l.forc_topo_grc[g]; hsurf_c = topo.topo_col[c]
        lwrad_g = a2l.forc_lwrad_not_downscaled_grc[g]
        lwrad_c = lwrad_g - lr_lw * (hsurf_c - hsurf_g)
        lwrad_c = clamp(lwrad_c, lwrad_g * (1.0 - limit), lwrad_g * (1.0 + limit))
        a2l.forc_lwrad_downscaled_col[c] = max(lwrad_c, 0.0)
    end
    return a2l
end

# Move a whole fixture to device (Float32) and run the kernelized routine.
function run_device!(f)
    a2l = mf(f.a2l); col = mf(f.col); lun = mf(f.lun); topo = mf(f.topo)
    CLM.downscale_forcings!(f.bounds, a2l, col, lun, topo)
    Metal.synchronize()
    return a2l
end

const OUT = (:forc_t_downscaled_col, :forc_th_downscaled_col, :forc_pbot_downscaled_col,
             :forc_rho_downscaled_col, :forc_lwrad_downscaled_col, :forc_solad_downscaled_col,
             :forc_solar_downscaled_col, :forc_rain_downscaled_col, :forc_snow_downscaled_col,
             :forc_q_downscaled_col, :forc_pco2_grc, :forc_po2_grc)

function main(backend)
    println("="^66); println("  downscale_forcings! — host byte-identity + device parity"); println("="^66)

    # (1) byte-identity: KA.CPU host path == original explicit-loop reference
    R = reference!(build())
    H = build(); CLM.downscale_forcings!(H.bounds, H.a2l, H.col, H.lun, H.topo)
    nbad = 0
    println("\n  [byte-identity] host KA.CPU (Float64) vs original loop reference:")
    for nm in OUT
        d = maxabs(getfield(H.a2l, nm), getfield(R, nm))
        ok = d == 0.0
        @printf("    [%s] %-28s maxabs=%.3e\n", ok ? "OK" : "XX", nm, d)
        ok || (nbad += 1)
    end
    if nbad != 0
        println("\n  HOST PATH NOT BYTE-IDENTICAL ($nbad fields) — default path changed!"); return 1
    end
    println("    → host default path is BYTE-IDENTICAL to the pre-kernelization loop.")

    # (2) device Float32 parity vs host Float64
    if backend === nothing; println("\n  No GPU backend — skipping device parity."); return 0; end
    name, _, FT = backend
    @printf("\n  Backend: %s (%s)\n", name, FT)
    D = run_device!(build())
    nfail = 0; ncmp = 0
    for nm in OUT
        r, k = reldiff(getfield(H.a2l, nm), getfield(D, nm)); ncmp += k
        ok = r < 1f-4
        @printf("    [%s] %-28s rel=%.2e over %d finite\n", ok ? "PASS" : "FAIL", nm, r, k)
        ok || (nfail += 1)
    end
    println()
    println(nfail == 0 ?
        "  downscale_forcings! kernels MATCH host on $name over $ncmp finite outputs" :
        "  DIVERGENCE ($nfail fields).")
    return nfail == 0 ? 0 : 1
end

exit(main(detect_backend()))
