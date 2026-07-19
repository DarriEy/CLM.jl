# ==========================================================================
# gpu_validate_clusterC_misc.jl — GPU parity for the Cluster-C default-hot
# kernelizations (Phase 3 Batch 2):
#   1. dust_emission_zender2003!         (DEFAULT dust scheme, per-patch)
#   2. soilbiogeochem_n_state_update1!   (mineral NH4/NO3 state update, per-column)
#
# Each builds a small Float64 fixture, runs the routine on the CPU, adapts the
# state + inputs to GPU, reruns on-device, and compares the mutated
# outputs. (depvel_compute! is intentionally excluded — it returns early when
# n_drydep==0, the default, and depends on non-device-callable Wesely-table
# helpers; see the fork report.)
#
#   julia --project=scripts scripts/gpu_validate_clusterC_misc.jl
# ==========================================================================

using CLM, Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))
include(joinpath(@__DIR__, "gpu_adapt.jl"))   # shared Float32 struct adaptor `mf`

to(x)    = device_array_type()(x)
dmask(m) = to(collect(Bool, m))
tof(x)   = to(Float32.(x))
mfs(x)   = mf(device_array_type(), x)   # handles DustEmisBaseData + Soil-BGC N structs

reldiff(a, b) = begin
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + max(abs(Float64(A[i])), abs(Float64(B[i])))))
    end
    m
end

# --------------------------------------------------------------------------
# 1. dust_emission_zender2003!
# --------------------------------------------------------------------------
function dust_scene()
    dust = CLM.DustEmisBaseData(); CLM.dust_emis_init!(dust, 1; nc=1)
    return (; dust, nolakep_mask=trues(1), patch_active=trues(1), patch_column=[1],
        patch_landunit=[1], patch_wtlunit=[1.0], lun_itype=[CLM.ISTSOIL],
        forc_rho=[1.225], gwc_thr=[0.17], mss_frc_cly_vld=[0.15],
        watsat=reshape([0.45], 1, 1), tlai=[0.0], tsai=[0.0], frac_sno=[0.0],
        h2osoi_vol=reshape([0.02], 1, 1), h2osoi_liq=reshape([30.0], 1, 1),
        h2osoi_ice=reshape([0.0], 1, 1), fv=[1.2], u10=[15.0])
end
runcpu_dust(s) = CLM.dust_emission_zender2003!(s.dust, s.nolakep_mask, s.patch_active,
    s.patch_column, s.patch_landunit, s.patch_wtlunit, s.lun_itype, 1:1, 1,
    s.forc_rho, s.gwc_thr, s.mss_frc_cly_vld, s.watsat, s.tlai, s.tsai, s.frac_sno,
    s.h2osoi_vol, s.h2osoi_liq, s.h2osoi_ice, s.fv, s.u10)

function check_dust()
    println("\n  --- dust_emission_zender2003! ---")
    cpu = dust_scene(); runcpu_dust(cpu)
    cbin = copy(cpu.dust.flx_mss_vrt_dst_patch); ctot = copy(cpu.dust.flx_mss_vrt_dst_tot_patch)
    @printf("    CPU total = %.6e  (nonzero bins: %d)\n", ctot[1], count(>(0.0), cbin[1, :]))
    s = dust_scene(); dust_d = mfs(s.dust)
    CLM.dust_emission_zender2003!(dust_d, dmask(s.nolakep_mask), dmask(s.patch_active),
        to(s.patch_column), to(s.patch_landunit), tof(s.patch_wtlunit), to(s.lun_itype),
        1:1, 1, tof(s.forc_rho), tof(s.gwc_thr), tof(s.mss_frc_cly_vld), tof(s.watsat),
        tof(s.tlai), tof(s.tsai), tof(s.frac_sno), tof(s.h2osoi_vol), tof(s.h2osoi_liq),
        tof(s.h2osoi_ice), tof(s.fv), tof(s.u10))
    dust_d.flx_mss_vrt_dst_patch isa device_array_type() || (println("    BLOCKED: dust not on device."); return 1)
    db = reldiff(cbin, dust_d.flx_mss_vrt_dst_patch)
    dt = reldiff(ctot, dust_d.flx_mss_vrt_dst_tot_patch)
    ok = db < 1e-4 && dt < 1e-4
    @printf("    [%s] rel|Δbin|=%.2e  rel|Δtot|=%.2e\n", ok ? "PASS" : "FAIL", db, dt)
    return ok ? 0 : 1
end

# --------------------------------------------------------------------------
# 2. soilbiogeochem_n_state_update1!
# --------------------------------------------------------------------------
function nsu_scene(; nc=4, nlevdecomp=4)
    ndecomp_pools = 7; ndct = 10
    nf = CLM.SoilBiogeochemNitrogenFluxData()
    CLM.soil_bgc_nitrogen_flux_init!(nf, nc, nlevdecomp, ndecomp_pools, ndct)
    ns = CLM.SoilBiogeochemNitrogenStateData()
    CLM.soil_bgc_nitrogen_state_init!(ns, nc, 1, nlevdecomp, ndecomp_pools)
    st = CLM.SoilBiogeochemStateData()
    CLM.soil_bgc_state_init!(st, nc, 1, nlevdecomp, ndct)
    for c in 1:nc
        nf.ndep_to_sminn_col[c] = 1.0e-6; nf.nfix_to_sminn_col[c] = 5.0e-7
        for j in 1:nlevdecomp
            ns.smin_nh4_vr_col[c, j] = 1.0; ns.smin_no3_vr_col[c, j] = 2.0
            ns.sminn_vr_col[c, j]    = 3.0
            nf.gross_nmin_vr_col[c, j]           = 1.0e-7
            nf.actual_immob_nh4_vr_col[c, j]     = 5.0e-8
            nf.actual_immob_no3_vr_col[c, j]     = 3.0e-8
            nf.smin_nh4_to_plant_vr_col[c, j]    = 2.0e-8
            nf.smin_no3_to_plant_vr_col[c, j]    = 1.0e-8
            nf.sminn_to_plant_fun_nh4_vr_col[c, j] = 1.5e-8
            nf.sminn_to_plant_fun_no3_vr_col[c, j] = 0.8e-8
            nf.f_nit_vr_col[c, j]                = 4.0e-8
            nf.f_denit_vr_col[c, j]              = 2.0e-8
            nf.supplement_to_sminn_vr_col[c, j]  = 1.0e-8
            st.ndep_prof_col[c, j]       = 0.5
            st.nfixation_prof_col[c, j]  = 0.3
        end
    end
    return (; ns, nf, st, mask=trues(nc), nc, nlevdecomp)
end
run_nsu!(d) = CLM.soilbiogeochem_n_state_update1!(d.ns, d.nf, d.st;
    mask_bgc_soilc=d.mask, bounds_col=1:d.nc, nlevdecomp=d.nlevdecomp, dt=1800.0)

function check_nsu()
    println("\n  --- soilbiogeochem_n_state_update1! ---")
    H = nsu_scene(); run_nsu!(H)
    B = nsu_scene()
    D = (; ns=mfs(B.ns), nf=mfs(B.nf), st=mfs(B.st),
         mask=device_array_type()(collect(B.mask)), nc=B.nc, nlevdecomp=B.nlevdecomp)
    D.ns.smin_nh4_vr_col isa device_array_type() || (println("    BLOCKED: ns not on device."); return 1)
    run_nsu!(D)
    nfail = 0
    for f in (:smin_nh4_vr_col, :smin_no3_vr_col, :sminn_vr_col)
        dd = reldiff(getfield(H.ns, f), getfield(D.ns, f)); ok = dd < 1f-3
        @printf("    [%s] %-18s rel=%.3e\n", ok ? "PASS" : "FAIL", f, dd); ok || (nfail += 1)
    end
    return nfail
end

function main(backend)
    println("=" ^ 72); println("Cluster-C misc kernelizations — GPU parity"); println("=" ^ 72)
    if backend === nothing; println("  No GPU backend."); return 0; end
    name, _, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n", name, FT)
    n = check_dust() + check_nsu()
    println()
    println(n == 0 ? "  ★★ Cluster-C misc kernels MATCH CPU on $name" : "  DIVERGENCE ($n failed).")
    return n == 0 ? 0 : 1
end
exit(main(detect_backend()))
