# ==========================================================================
# gpu_validate_veg_struct_update.jl — GPU parity for cn_veg_struct_update!
# (the per-patch vegetation-structure diagnosis kernelized in veg_struct_update.jl).
#
# Builds the mixed-PFT fixture from test/test_veg_struct_update.jl (noveg + woody
# tree + grass + prognostic crop, so every branch runs), runs cn_veg_struct_update!
# on the CPU (Float64 reference), adapts every state struct + the mask to Metal
# (Float32) via the shared mf adaptor, runs the SAME call on the device (pftcon
# stays a host-global — the function moves its params onto the backend internally),
# and compares the mutated canopy/cnveg outputs with reldiff.
#
#   julia --project=scripts scripts/gpu_validate_veg_struct_update.jl
# ==========================================================================

using CLM
using Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))
include(joinpath(@__DIR__, "gpu_adapt.jl"))
mf(x) = mf(device_array_type(), x)

function reldiff_finite(H, D)
    A = Array(H); B = Array(D); m = 0.0; n = 0
    for i in eachindex(A, B)
        (isfinite(A[i]) && isfinite(B[i])) || continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + abs(Float64(A[i])))); n += 1
    end
    return m, n
end

# Mixed-PFT fixture (mirrors make_veg_struct_data + test 24): 4 patches spanning
# noveg / woody tree / grass / prognostic crop, 2 columns. Float64 host structs.
function build_fixture()
    np = 4; nc = 2; npft = 80

    patch = CLM.PatchData()
    patch.itype    = [0, 2, 13, 17]        # noveg, needleleaf tree, c3 grass, temperate corn
    patch.column   = [1, 1, 2, 2]
    patch.gridcell = ones(Int, np)
    patch.active   = trues(np)

    cs = CLM.CanopyStateData()
    cs.tlai_patch               = [0.0, 0.5, 0.0, 0.0]
    cs.tsai_patch               = [0.0, 1.0, 1.0, 0.5]
    cs.elai_patch               = zeros(np)
    cs.esai_patch               = zeros(np)
    cs.htop_patch               = zeros(np)
    cs.hbot_patch               = zeros(np)
    cs.stem_biomass_patch       = zeros(np)
    cs.leaf_biomass_patch       = zeros(np)
    cs.frac_veg_nosno_alb_patch = zeros(Int, np)

    cnveg_cs = CLM.CNVegCarbonStateData()
    cnveg_cs.leafc_patch     = [0.0, 80.0, 50.0, 100.0]
    cnveg_cs.deadstemc_patch = [0.0, 3000.0, 0.0, 0.0]
    cnveg_cs.livestemc_patch = [0.0, 500.0, 0.0, 0.0]

    wd = CLM.WaterDiagnosticBulkData()
    wd.frac_sno_col   = [0.3, 0.5]
    wd.snow_depth_col = [0.2, 0.05]

    fv = CLM.FrictionVelocityData()
    fv.forc_hgt_u_patch = fill(30.0, np)

    vs = CLM.CNVegStateData()
    vs.farea_burned_col = [0.1, 0.0]
    vs.htmx_patch       = zeros(np)
    vs.peaklai_patch    = zeros(Int, np)

    cr = CLM.CropData()
    cr.harvdate_patch = fill(999, np)

    pft = CLM.PftconType()
    pft.woody    = zeros(npft);      pft.woody[3] = 1.0
    pft.slatop   = fill(0.02, npft); pft.slatop[3] = 0.012
    pft.dsladlai = zeros(npft)
    pft.z0mr     = fill(0.055, npft)
    pft.displar  = fill(0.67, npft)
    pft.dwood    = fill(2.5e5, npft)
    pft.ztopmx   = fill(1.5, npft);  pft.ztopmx[18] = 2.5
    pft.laimx    = fill(6.0, npft)
    pft.nstem    = fill(0.005, npft)
    pft.taper    = fill(200.0, npft)
    pft.fbw      = fill(0.5, npft)

    return (; patch, cs, cnveg_cs, wd, fv, vs, cr, pft, mask = trues(np), bounds = 1:np)
end

const KW = (dt=1800.0, noveg=0, npcropmin=17, ntmp_corn=17,
            use_biomass_heat_storage=true, spinup_factor_deadwood=1.0)

run_it(d, mask) = CLM.cn_veg_struct_update!(mask, d.bounds, d.patch, d.cs,
    d.cnveg_cs, d.wd, d.fv, d.vs, d.cr, d.pft; KW...)

function main(backend)
    println("="^70); println("GPU parity — cn_veg_struct_update! (per-patch kernel)"); println("="^70)
    if backend === nothing; println("  No GPU backend."); return 0; end
    name, _, FT = backend
    @printf("  Backend: %s (%s)\n\n", name, FT)

    # ---- host reference (Float64) ----
    dH = build_fixture(); run_it(dH, dH.mask)

    # ---- device (Float32) ----
    dD = build_fixture()
    dev = (patch = mf(dD.patch), cs = mf(dD.cs), cnveg_cs = mf(dD.cnveg_cs),
           wd = mf(dD.wd), fv = mf(dD.fv), vs = mf(dD.vs), cr = mf(dD.cr),
           pft = dD.pft, bounds = dD.bounds)          # pft stays host (moved internally)
    mask_d = mf(dD.mask)
    dev.cs.tlai_patch isa device_array_type() || (println("  BLOCKED: adapt did not reach device."); return 2)
    CLM.cn_veg_struct_update!(mask_d, dev.bounds, dev.patch, dev.cs, dev.cnveg_cs,
        dev.wd, dev.fv, dev.vs, dev.cr, dev.pft; KW...)
    device_synchronize()
    println("  cn_veg_struct_update! completed on $name.\n")

    checks = [("tlai", dH.cs.tlai_patch, dev.cs.tlai_patch),
              ("tsai", dH.cs.tsai_patch, dev.cs.tsai_patch),
              ("htop", dH.cs.htop_patch, dev.cs.htop_patch),
              ("hbot", dH.cs.hbot_patch, dev.cs.hbot_patch),
              ("elai", dH.cs.elai_patch, dev.cs.elai_patch),
              ("esai", dH.cs.esai_patch, dev.cs.esai_patch),
              ("stem_biomass", dH.cs.stem_biomass_patch, dev.cs.stem_biomass_patch),
              ("leaf_biomass", dH.cs.leaf_biomass_patch, dev.cs.leaf_biomass_patch),
              ("frac_veg_nosno_alb", dH.cs.frac_veg_nosno_alb_patch, dev.cs.frac_veg_nosno_alb_patch),
              ("htmx", dH.vs.htmx_patch, dev.vs.htmx_patch),
              ("peaklai", dH.vs.peaklai_patch, dev.vs.peaklai_patch)]
    nfail = 0; ncmp = 0
    for (nm, h, d) in checks
        r, n = reldiff_finite(h, d); ncmp += n
        ok = r < 1f-3
        @printf("  [%s] %-20s rel=%.2e over %d finite\n", ok ? "PASS" : "FAIL", nm, r, n)
        ok || (nfail += 1)
    end
    println()
    if ncmp == 0
        println("  ⚠ all comparisons non-finite — fixture problem."); return 1
    end
    println(nfail == 0 ? "  cn_veg_struct_update! MATCHES host ON $name over $ncmp finite outputs" :
                         "  DIVERGENCE ($nfail fields).")
    return nfail == 0 ? 0 : 1
end
exit(main(detect_backend()))
