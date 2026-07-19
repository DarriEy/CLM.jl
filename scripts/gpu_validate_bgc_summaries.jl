# ==========================================================================
# gpu_validate_bgc_summaries.jl — host-vs-GPU parity for the BGC state-summary
# reductions (soil_bgc_carbon_state_summary! + soil_bgc_nitrogen_state_summary!),
# kernelized as one-thread-per-column fused reductions.
#
#   julia --project=scripts scripts/gpu_validate_bgc_summaries.jl
# ==========================================================================
using CLM
using Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))
include(joinpath(@__DIR__, "gpu_adapt.jl"))
mf(x) = mf(device_array_type(), x)

function reldiff(H, D)
    A = Array(H); B = Array(D); m = 0.0; n = 0
    for i in eachindex(A, B)
        (isfinite(A[i]) && isfinite(B[i])) || continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + abs(Float64(A[i])))); n += 1
    end
    return m, n
end

# Build a carbon-state fixture with varied (non-uniform) vr pools so a dropped
# level/pool would show up as a mismatch.
function make_cs(nc, ng, nlev, npools)
    cs = CLM.SoilBiogeochemCarbonStateData()
    CLM.soil_bgc_carbon_state_init!(cs, nc, ng, nlev, npools)
    for c in 1:nc, j in 1:nlev, l in 1:npools
        cs.decomp_cpools_vr_col[c, j, l] = 1.0 + 0.1c + 0.5j + 2.0l
    end
    cs.ctrunc_vr_col .= 0.3
    return cs
end
function make_ns(nc, ng, nlev, npools)
    ns = CLM.SoilBiogeochemNitrogenStateData()
    CLM.soil_bgc_nitrogen_state_init!(ns, nc, ng, nlev, npools)
    for c in 1:nc, j in 1:nlev, l in 1:npools
        ns.decomp_npools_vr_col[c, j, l] = 0.2 + 0.05c + 0.1j + 0.3l
    end
    ns.sminn_vr_col   .= 0.4
    ns.ntrunc_vr_col  .= 0.02
    ns.smin_no3_vr_col .= 0.15
    ns.smin_nh4_vr_col .= 0.25
    return ns
end

function main(backend)
    println("="^64); println("BGC summary reductions — host vs device parity"); println("="^64)
    if backend === nothing; println("  No GPU backend."); return 0; end
    name, _, FT = backend
    @printf("  Backend: %s (%s)\n", name, FT)

    nc = 4; ng = 1; nlev = 3; npools = 7
    dzsoi = [0.1, 0.3, 0.6]
    zisoi = [0.0, 0.1, 0.4, 1.0]
    is_litter  = BitVector([true, true, true, false, false, false, false])
    is_soil    = BitVector([false, false, false, true, true, true, false])
    is_microbe = BitVector([false, false, false, false, false, true, false])
    is_cwd     = BitVector([false, false, false, false, false, false, true])
    mask = trues(nc)

    nfail = 0; ncmp = 0

    # ---- carbon ----
    csH = make_cs(nc, ng, nlev, npools); csD = mf(make_cs(nc, ng, nlev, npools))
    csD.decomp_cpools_vr_col isa device_array_type() || (println("  BLOCKED: cs adapt"); return 2)
    kw = (; nlevdecomp=nlev, nlevdecomp_full=nlev, ndecomp_pools=npools,
          dzsoi_decomp_vals=dzsoi, zisoi_vals=zisoi,
          is_litter=is_litter, is_soil=is_soil, is_microbe=is_microbe, is_cwd=is_cwd)
    CLM.soil_bgc_carbon_state_summary!(csH, mask, 1:nc; kw...)
    CLM.soil_bgc_carbon_state_summary!(csD, mf(mask), 1:nc; kw...); device_synchronize()
    for f in (:decomp_cpools_col, :decomp_cpools_1m_col, :decomp_soilc_vr_col, :ctrunc_col,
              :totlitc_col, :totsomc_col, :totmicc_col, :totlitc_1m_col, :totsomc_1m_col,
              :cwdc_col, :totecosysc_col, :totc_col)
        r, n = reldiff(getfield(csH, f), getfield(csD, f)); ncmp += n
        ok = r < 1f-3; ok || (nfail += 1)
        @printf("  [%s] C.%-22s rel=%.2e over %d\n", ok ? "PASS" : "FAIL", f, r, n)
    end

    # ---- nitrogen ----
    nsH = make_ns(nc, ng, nlev, npools); nsD = mf(make_ns(nc, ng, nlev, npools))
    kwn = (; nlevdecomp=nlev, nlevdecomp_full=nlev, ndecomp_pools=npools,
          dzsoi_decomp_vals=dzsoi, zisoi_vals=zisoi,
          is_litter=is_litter, is_soil=is_soil, is_microbe=is_microbe, is_cwd=is_cwd,
          use_nitrif_denitrif=true)
    CLM.soil_bgc_nitrogen_state_summary!(nsH, mask, 1:nc; kwn...)
    CLM.soil_bgc_nitrogen_state_summary!(nsD, mf(mask), 1:nc; kwn...); device_synchronize()
    for f in (:decomp_npools_col, :decomp_npools_1m_col, :decomp_soiln_vr_col, :sminn_col,
              :ntrunc_col, :smin_no3_col, :smin_nh4_col, :totlitn_col, :totsomn_col,
              :totmicn_col, :totlitn_1m_col, :totsomn_1m_col, :cwdn_col, :totecosysn_col, :totn_col)
        r, n = reldiff(getfield(nsH, f), getfield(nsD, f)); ncmp += n
        ok = r < 1f-3; ok || (nfail += 1)
        @printf("  [%s] N.%-22s rel=%.2e over %d\n", ok ? "PASS" : "FAIL", f, r, n)
    end

    println()
    println(nfail == 0 ? "  BGC SUMMARIES MATCH host on $name over $ncmp finite outputs" : "  DIVERGENCE ($nfail).")
    return nfail == 0 ? 0 : 1
end
exit(main(detect_backend()))
