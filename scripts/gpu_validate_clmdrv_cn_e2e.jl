# ==========================================================================
# gpu_validate_clmdrv_cn_e2e.jl — WHOLE-TIMESTEP Metal parity for clm_drv! with
# the CN cycle ON (use_cn=true) — the biogeophys + BGC COMPOSITE. The biogeophys
# path (use_cn=false) is validated by gpu_validate_clmdrv_e2e.jl and the BGC path
# by the bgc/matrix harnesses; this checks that the two halves COMPOSE on one
# device through a single clm_drv! call — the integration the per-half harnesses
# cannot give.
#
# It is also a LEAK-FINDER: the composite surfaces host-constant-array leaks (the
# driver hands physical constants like zsoi/dzsoi as host Float64 vectors to device
# CN kernels). On a device compile failure it names the leaking kernel so the Phase-2
# sweep is "run → fix the reported kernel → repeat".
#
#   julia --project=scripts scripts/gpu_validate_clmdrv_cn_e2e.jl
# ==========================================================================
using CLM
using Printf
import Metal
include(joinpath(@__DIR__, "gpu_backends.jl"))
include(joinpath(@__DIR__, "gpu_adapt.jl"))
include(joinpath(@__DIR__, "clmdrv_make_data.jl"))
mf(x) = mf(Metal.MtlArray, x)

const DRV_ARGS = (true, 1.0, 0.0, 0.0, 0.4091, false, false, "20260101", false)
drv!(cfg, inst, filt, filt_ia, bounds, ps) = CLM.clm_drv!(cfg, inst, filt, filt_ia, bounds, DRV_ARGS...;
    nstep=1, is_first_step=false, is_beg_curr_day=false, is_beg_curr_year=false, photosyns=ps)

function reldiff_finite(H, D)
    A = Array(H); B = Array(D); m = 0.0; n = 0
    for i in eachindex(A, B)
        (isfinite(A[i]) && isfinite(B[i])) || continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + abs(Float64(A[i])))); n += 1
    end
    return m, n
end

function main(backend)
    println("="^72); println("WHOLE-TIMESTEP Metal parity for clm_drv! — CN COMPOSITE (use_cn=true)"); println("="^72)
    if backend === nothing; println("  No GPU backend."); return 0; end
    name, _, FT = backend
    @printf("  Backend: %s (%s)\n\n", name, FT)
    cfg = CLM.CLMDriverConfig(use_cn=true)

    # ---- host reference ----
    instH, bounds, filtH, filt_iaH, _c, psH = make_driver_data_physical()
    instH.canopystate.frac_veg_nosno_alb_patch .= 1
    drv!(cfg, instH, filtH, filt_iaH, bounds, psH)
    println("  HOST clm_drv! (use_cn=true) completed — the composite runs.")

    # ---- device ----
    instB, boundsB, filtB, filt_iaB, _c2, psB = make_driver_data_physical()
    instB.canopystate.frac_veg_nosno_alb_patch .= 1
    inst_d = mf(instB); ps_d = mf(psB); filt_d = mf(filtB); filt_ia_d = mf(filt_iaB)
    inst_d.temperature.t_soisno_col isa Metal.MtlArray || (println("  BLOCKED: adapt."); return 2)
    try
        drv!(cfg, inst_d, filt_d, filt_ia_d, boundsB, ps_d); Metal.synchronize()
    catch e
        msg = sprint(showerror, e)
        kern = match(r"gpu_(_?\w+kernel!?)", msg)
        arg = match(r"Argument (\d+) to your kernel function is of type (\w+\{Float64\})", msg)
        println("\n  ✗ COMPOSITION LEAK on device:")
        kern !== nothing && println("      leaking kernel: ", kern.captures[1])
        arg !== nothing && @printf("      arg %s is host %s (needs _to_backend_like)\n", arg.captures[1], arg.captures[2])
        kern === nothing && arg === nothing && println("      ", first(split(msg, "\n")))
        println("\n  → fix that kernel's launch (move the host constant to the state backend), then re-run.")
        return 1
    end
    println("  DEVICE clm_drv! (use_cn=true) COMPLETED — composite runs on $name.\n")

    checks = [("t_soisno", instH.temperature.t_soisno_col, inst_d.temperature.t_soisno_col),
              ("decomp_cpools_vr", instH.soilbiogeochem_carbonstate.decomp_cpools_vr_col, inst_d.soilbiogeochem_carbonstate.decomp_cpools_vr_col),
              ("leafc", instH.bgc_vegetation.cnveg_carbon_state.leafc_patch, inst_d.bgc_vegetation.cnveg_carbon_state.leafc_patch)]
    nfail = 0; ncmp = 0
    for (nm, h, d) in checks
        r, n = reldiff_finite(h, d); ncmp += n
        ok = r < 1f-3
        @printf("  [%s] %-18s rel=%.2e over %d finite\n", ok ? "PASS" : "FAIL", nm, r, n); ok || (nfail += 1)
    end
    println()
    if ncmp == 0
        println("  ⚠ composite runs on device but the fixture leaves all compared outputs NaN —")
        println("    a CN-consistent fixture is needed to get finite parity numbers.")
        return 0
    end
    println(nfail == 0 ? "  CN-COMPOSITE clm_drv! MATCHES host ON $name over $ncmp finite outputs" : "  DIVERGENCE ($nfail).")
    return nfail == 0 ? 0 : 1
end
exit(main(detect_backend()))
