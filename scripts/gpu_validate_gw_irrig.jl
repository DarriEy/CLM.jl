# ==========================================================================
# gpu_validate_gw_irrig.jl — Metal parity for the use_groundwater_irrigation
# path in soil_hydrology.jl. These were already kernelized (calc_irrig_withdrawals!
# via _soilhyd_irrig_withdrawals_kernel! + a device-guarded host-only negative-demand
# pre-check; withdraw_groundwater_irrigation! -> soilhyd_withdraw_gw_lyr!/_con! kernels).
# This confirms all three kernels run on Metal at parity vs the CPU (KA) path.
#
#   julia +1.12 --project=scripts scripts/gpu_validate_gw_irrig.jl
# ==========================================================================
using CLM, Printf
import Metal
include(joinpath(@__DIR__, "gpu_backends.jl"))

df(x) = Metal.MtlArray(Float32.(x))          # float array -> device Float32
di(x) = Metal.MtlArray(x)                    # int/bool array -> device (keep type)
function reldiff(H, D)
    A = Array(H); B = Array(D); m = 0.0
    for i in eachindex(A, B)
        (isfinite(A[i]) && isfinite(B[i])) || continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + abs(Float64(A[i]))))
    end
    return m
end

function main()
    Metal.functional() || (println("Metal not functional."); return 0)
    println("="^60); println("  Groundwater-irrigation kernels — Metal parity"); println("="^60)
    nc = 8; nlev = 10; dtime = 1800.0
    mask = fill(true, nc)
    # --- fixtures (realistic-ish, all columns irrigating) ---
    h2osoi_liq = fill(30.0, nc, nlev + 5)          # padded like the real col array
    wa         = fill(4800.0, nc)
    qflx_uncon = fill(1.0e-4, nc, nlev + 5)
    qflx_con   = fill(2.0e-4, nc)
    demand     = fill(5.0e-4, nc)
    nbedrock   = fill(nlev, nc)
    zi   = repeat(collect(0.0:0.2:0.2*(nlev+4))', nc, 1)
    zwt  = fill(1.5, nc)
    watsat = fill(0.4, nc, nlev + 5); sucsat = fill(100.0, nc, nlev + 5); bsw = fill(5.0, nc, nlev + 5)
    aqy = 0.02

    checks = []
    # ---- 1. soilhyd_withdraw_gw_lyr! (per col,layer) ----
    hl_c = copy(h2osoi_liq); CLM.soilhyd_withdraw_gw_lyr!(hl_c, mask, qflx_uncon, nlev, dtime)
    hl_d = df(h2osoi_liq);   CLM.soilhyd_withdraw_gw_lyr!(hl_d, di(mask), df(qflx_uncon), nlev, dtime); Metal.synchronize()
    push!(checks, ("withdraw_gw_lyr h2osoi_liq", reldiff(hl_c, hl_d)))
    # ---- 2. soilhyd_withdraw_gw_con! (per col) ----
    wa_c = copy(wa); CLM.soilhyd_withdraw_gw_con!(wa_c, mask, qflx_con, dtime)
    wa_d = df(wa);   CLM.soilhyd_withdraw_gw_con!(wa_d, di(mask), df(qflx_con), dtime); Metal.synchronize()
    push!(checks, ("withdraw_gw_con wa", reldiff(wa_c, wa_d)))
    # ---- 3. calc_irrig_withdrawals kernel (per col; water-table split) ----
    lyr_c = zeros(nc, nlev + 5); con_c = zeros(nc)
    CLM._launch!(CLM._soilhyd_irrig_withdrawals_kernel!, lyr_c, con_c, mask, demand,
        nbedrock, zi, zwt, watsat, sucsat, bsw, nlev, dtime, aqy; ndrange = nc)
    lyr_d = df(zeros(nc, nlev + 5)); con_d = df(zeros(nc))
    CLM._launch!(CLM._soilhyd_irrig_withdrawals_kernel!, lyr_d, con_d, di(mask), df(demand),
        di(nbedrock), df(zi), df(zwt), df(watsat), df(sucsat), df(bsw), nlev, Float32(dtime), Float32(aqy); ndrange = nc)
    Metal.synchronize()
    push!(checks, ("calc_irrig_withdrawals lyr", reldiff(lyr_c, lyr_d)))
    push!(checks, ("calc_irrig_withdrawals con", reldiff(con_c, con_d)))

    nfail = 0
    for (nm, r) in checks
        ok = r < 1f-3; ok || (nfail += 1)
        @printf("  [%s] %-30s rel=%.2e\n", ok ? "PASS" : "FAIL", nm, r)
    end
    println(nfail == 0 ? "  GROUNDWATER-IRRIGATION kernels MATCH host on Metal" : "  DIVERGENCE ($nfail)")
    return nfail
end
exit(main())
