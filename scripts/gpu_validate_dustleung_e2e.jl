# gpu_validate_dustleung_e2e.jl — Leung (2023) dust emission on Metal parity.
#
# Dry erodible high-wind scene (1 landunit / column / patch). Runs
# dust_emission_leung2023! on the CPU, then adapts the DustEmisBaseData struct +
# every input array to Metal (Float32) and reruns on-device, comparing the
# per-bin flux and the total.
#
# Run: julia +1.12 --project=scripts scripts/gpu_validate_dustleung_e2e.jl

using CLM, Printf
import Metal
include(joinpath(@__DIR__, "gpu_adapt.jl"))
ad(x)    = mf(Metal.MtlArray, x)
to(x)    = Metal.MtlArray(x)
dmask(m) = to(collect(Bool, m))
tof(x)   = to(Float32.(x))

function scene()
    dust = CLM.DustEmisBaseData(); CLM.dust_emis_init!(dust, 1; nc=1)
    dust.dpfct_rock_patch[1] = 1.0
    return (; dust, nolakep_mask = trues(1), patch_active = trues(1), patch_column = [1],
        patch_landunit = [1], patch_itype = [CLM.noveg], patch_wtlunit = [1.0],
        lun_itype = [CLM.ISTSOIL], forc_rho = [1.225], gwc_thr = [0.17],
        mss_frc_cly_vld = [0.15], watsat = reshape([0.45], 1, 1), tlai = [0.0], tsai = [0.0],
        frac_sno = [0.0], h2osoi_vol = reshape([0.02], 1, 1), h2osoi_liq = reshape([30.0], 1, 1),
        h2osoi_ice = reshape([0.0], 1, 1), fv = [0.9], obu = [-50.0])
end

runcpu(s) = CLM.dust_emission_leung2023!(s.dust, s.nolakep_mask, s.patch_active,
    s.patch_column, s.patch_landunit, s.patch_itype, s.patch_wtlunit, s.lun_itype, 1:1, 1,
    s.forc_rho, s.gwc_thr, s.mss_frc_cly_vld, s.watsat, s.tlai, s.tsai, s.frac_sno,
    s.h2osoi_vol, s.h2osoi_liq, s.h2osoi_ice, s.fv, s.obu)

# ---- CPU reference --------------------------------------------------------
cpu = scene(); runcpu(cpu)
cpu_bin = copy(cpu.dust.flx_mss_vrt_dst_patch); cpu_tot = copy(cpu.dust.flx_mss_vrt_dst_tot_patch)
@printf("  CPU total = %.6e  (bins nonzero: %d)\n", cpu_tot[1], count(>(0.0), cpu_bin[1, :]))

# ---- Metal device run -----------------------------------------------------
if !Metal.functional()
    println("Metal not functional — skipping device leg."); exit(0)
end
s = scene()
dust_d = ad(s.dust)
CLM.dust_emission_leung2023!(dust_d, dmask(s.nolakep_mask), dmask(s.patch_active),
    to(s.patch_column), to(s.patch_landunit), to(s.patch_itype), tof(s.patch_wtlunit),
    to(s.lun_itype), 1:1, 1, tof(s.forc_rho), tof(s.gwc_thr), tof(s.mss_frc_cly_vld),
    tof(s.watsat), tof(s.tlai), tof(s.tsai), tof(s.frac_sno), tof(s.h2osoi_vol),
    tof(s.h2osoi_liq), tof(s.h2osoi_ice), tof(s.fv), tof(s.obu))
if !(dust_d.flx_mss_vrt_dst_patch isa Metal.MtlArray)
    println("  BLOCKED: DustEmisBaseData did not move to the device."); exit(1)
end
dev_bin = Array(dust_d.flx_mss_vrt_dst_patch); dev_tot = Array(dust_d.flx_mss_vrt_dst_tot_patch)

# ---- Compare (relative, since fluxes are ~1e-9) ---------------------------
scale = max(abs(cpu_tot[1]), 1e-30)
dbin = maximum(abs.(Float64.(vec(cpu_bin)) .- Float64.(vec(dev_bin))); init = 0.0)
dtot = abs(cpu_tot[1] - Float64(dev_tot[1]))
@printf("  Metal total = %.6e\n", dev_tot[1])
@printf("  max|Δbin| = %.2e   |Δtot| = %.2e   rel = %.2e\n", dbin, dtot, dtot / scale)
println(dtot / scale < 1e-4 ? "★★ dust_emission_leung2023! RUNS ON METAL — MATCHES CPU" :
                              "RAN on Metal but rel Δ=$(dtot/scale)")
