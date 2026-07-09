# ==========================================================================
# gpu_validate_depvel.jl — host (Float64) vs device (Float32) parity for the
# kernelized dry-deposition-velocity routine `depvel_compute!` (feature-gated on
# n_drydep>0). Exercises both Wesely hemispheres, several PFT->landtype branches,
# snow/no-snow, day/night, warm/cold, and reactive/unreactive species.
#
#   julia --project=scripts scripts/gpu_validate_depvel.jl
# ==========================================================================
using CLM
using Printf
import Metal
include(joinpath(@__DIR__, "gpu_backends.jl"))
include(joinpath(@__DIR__, "gpu_adapt.jl"))
mf(x) = mf(Metal.MtlArray, x)

function reldiff(H, D)
    A = Array(H); B = Array(D); m = 0.0; n = 0
    for i in eachindex(A, B)
        (isfinite(A[i]) && isfinite(B[i])) || continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + abs(Float64(A[i])))); n += 1
    end
    return m, n
end

# Fresh DryDepVelocityData (host, Float64) — velocity_patch is mutated by the run.
function build_dd(np)
    n_drydep = 3
    dd = CLM.DryDepVelocityData{Float64}()
    CLM.drydep_init!(dd, np, n_drydep;
        foxd = [1.0, 0.5, 0.0],       # reactive / mid / unreactive
        dv   = [0.25, 0.15, 0.10],    # diffusivity [cm^2/s]
        mapping = [1, 1, 1])
    return dd
end

# Loose forcing arrays (host, Float64), shared by host + device runs.
function build_forcings()
    np = 6; nc = 4; ng = 2
    mask_patch     = Bool[true, true, true, true, false, true]  # p5 masked out
    patch_gridcell = [1, 1, 2, 2, 1, 2]        # g1 NH (+45), g2 SH (-30)
    patch_column   = [1, 2, 3, 4, 1, 2]
    patch_landunit = [1, 1, 1, 1, 1, 1]        # unused by depvel
    patch_itype    = [0, 2, 6, 10, 15, 3]      # barren/conifer/decid/range/ag/conifer
    ram1_patch     = [50.0, 80.0, 30.0, 120.0, 60.0, 90.0]
    rb1_patch      = fill(10.0, np)            # unused by depvel
    fv_patch       = [0.30, 0.10, 0.50, 0.05, 0.20, 0.40]
    elai_patch     = [3.0, 1.5, 5.0, 0.5, 2.0, 4.0]
    forc_t_col     = [298.0, 260.0, 285.0, 275.0]   # warm/cold/mild/near-freeze
    forc_solar_col = [400.0, 0.0, 200.0, 100.0]     # day/night/…
    frac_sno       = [0.0, 0.8, 0.0, 0.3]
    lat_grc        = [45.0, -30.0]
    month          = 7                               # NH summer / SH winter
    return (; np, nc, ng, mask_patch, patch_gridcell, patch_column, patch_landunit,
            patch_itype, ram1_patch, rb1_patch, fv_patch, elai_patch,
            forc_t_col, forc_solar_col, frac_sno, lat_grc, month)
end

function run_depvel!(dd, f, dev)
    id(x) = dev ? mf(x) : x
    CLM.depvel_compute!(dd, id(f.mask_patch), 1:f.np,
        id(f.patch_gridcell), id(f.patch_column), id(f.patch_landunit), id(f.patch_itype),
        id(f.ram1_patch), id(f.rb1_patch), id(f.fv_patch), id(f.elai_patch),
        id(f.forc_t_col), id(f.forc_solar_col), id(f.frac_sno), id(f.lat_grc), f.month)
    dev && Metal.synchronize()
    return dd.velocity_patch
end

function main(backend)
    println("="^64); println("depvel_compute! — host vs device parity"); println("="^64)
    if backend === nothing; println("  No GPU backend."); return 0; end
    name, _, FT = backend
    @printf("  Backend: %s (%s)\n\n", name, FT)

    f = build_forcings()
    Hvel = run_depvel!(build_dd(f.np), f, false)
    Dvel = run_depvel!(mf(build_dd(f.np)), f, true)

    r, n = reldiff(Hvel, Dvel)
    ok = r < 1f-4
    @printf("  [%s] velocity_patch   rel=%.2e over %d finite\n", ok ? "PASS" : "FAIL", r, n)
    println()
    println(ok ? "  depvel_compute! kernel MATCHES host on $name over $n finite outputs" :
                 "  DIVERGENCE.")
    return ok ? 0 : 1
end
exit(main(detect_backend()))
