# ==========================================================================
# gpu_validate_phs_e2e.jl — WHOLE-FUNCTION GPU parity for photosynthesis_hydrstress!
# (the PHS / plant-hydraulic-stress photosynthesis path), all 4 passes on device.
#
# Builds the same small Float32 PHS case as the CPU oracle
# (gpu_validate_phs_cpu_oracle.jl), runs photosynthesis_hydrstress! on the CPU,
# adapts ps + every input array to the GPU, runs the SAME call on device, and
# compares the mutated outputs. This is the parity gate for Passes 1-4 together.
# Passes 1/2/4 run on the GPU; Pass 3's fused Newton solve runs via the HYBRID
# fallback (on the CPU over host copies, results copied back to the device) because
# that kernel exceeds the Apple Metal shader compiler's capacity. So vegwp/qflx_tran
# here exercise the hybrid path, while k_soil_root/vcmax_z/cisun_z exercise on-GPU.
#
#   julia --project=scripts scripts/gpu_validate_phs_e2e.jl
# ==========================================================================

using CLM
using Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))

reldiff(a, b) = begin
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(A[i] - B[i]) / (1.0 + max(abs(A[i]), abs(B[i]))))
    end
    m
end

# Reuse the oracle harness's setparams!/build/run_phs! verbatim (same inputs).
include(joinpath(@__DIR__, "gpu_validate_phs_cpu_oracle.jl"))  # defines setparams!, build, run_phs!

function main_e2e()
    println("=" ^ 72)
    println("WHOLE-FN GPU parity for photosynthesis_hydrstress! (PHS, all 4 passes)")
    println("=" ^ 72)
    if !gpu_functional()
        println("  No GPU backend detected — nothing to validate."); return 0
    end
    CLM.varpar_init!(CLM.varpar, 5, 16, 0, 5)
    nlevsoi = 5
    FT = Float32

    setparams!()
    H = build(FT, nlevsoi)        # CPU reference
    setparams!()
    B = build(FT, nlevsoi)        # to be adapted to device

    ad(x) = CLM.Adapt.adapt(device_array_type(), x)
    # Adapt ps + every input array + vegwp/vegwp_ln; mask -> device Bool.
    Bd = (; np = B.np, nc = B.nc, nlevsoi = B.nlevsoi, nlevcan = B.nlevcan,
            ps = ad(B.ps),
            vegwp = ad(B.vegwp), vegwp_ln = ad(B.vegwp_ln),
            inp = (; (k => ad(getfield(B.inp, k)) for k in keys(B.inp))...),
            mask = device_array_type()(collect(Bool, B.mask)))

    if !(Bd.ps.ac_phs_patch isa device_array_type())
        println("  BLOCKED: PhotosynthesisData did not move to the device under adapt."); return 2
    end

    run_phs!(H)          # CPU
    run_phs!(Bd)         # device

    checks = [
        ("psnsun",     H.ps.psnsun_patch,     Bd.ps.psnsun_patch),
        ("psnsha",     H.ps.psnsha_patch,     Bd.ps.psnsha_patch),
        ("an_sun",     H.ps.an_sun_patch,     Bd.ps.an_sun_patch),
        ("an_sha",     H.ps.an_sha_patch,     Bd.ps.an_sha_patch),
        ("gs_mol_sun", H.ps.gs_mol_sun_patch, Bd.ps.gs_mol_sun_patch),
        ("gs_mol_sha", H.ps.gs_mol_sha_patch, Bd.ps.gs_mol_sha_patch),
        ("rssun",      H.ps.rssun_patch,      Bd.ps.rssun_patch),
        ("rssha",      H.ps.rssha_patch,      Bd.ps.rssha_patch),
        ("lmrsun",     H.ps.lmrsun_patch,     Bd.ps.lmrsun_patch),
        ("vcmax_z_phs",H.ps.vcmax_z_phs_patch,Bd.ps.vcmax_z_phs_patch),  # Pass 2
        ("cisun_z",    H.ps.cisun_z_patch,    Bd.ps.cisun_z_patch),
        ("vegwp",      H.vegwp,               Bd.vegwp),                 # Pass 3 Newton
        ("btran",      H.inp.btran,           Bd.inp.btran),             # Pass 4
        ("qflx_tran",  H.inp.qflx_tran_veg,   Bd.inp.qflx_tran_veg),
        ("k_soil_root",H.inp.k_soil_root,     Bd.inp.k_soil_root),       # Pass 1
    ]
    nfail = 0
    for (nm, a, b) in checks
        d = reldiff(a, b); ok = d < 1f-2
        @printf("  [%s] %-12s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end
    println()
    println(nfail == 0 ? "  WHOLE photosynthesis_hydrstress! MATCHES CPU ON GPU ✓" :
                         "  DIVERGENCE ($nfail) — investigate.")
    return nfail == 0 ? 0 : 1
end

exit(main_e2e())
