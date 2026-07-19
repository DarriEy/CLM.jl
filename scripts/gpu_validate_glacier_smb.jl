# ==========================================================================
# gpu_validate_glacier_smb.jl — host (Float64) vs device (Float32) parity for the
# three kernelized glacier surface-mass-balance routines:
#   handle_ice_melt! · compute_surface_mass_balance! · adjust_runoff_terms!
#
#   julia --project=scripts scripts/gpu_validate_glacier_smb.jl
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

# Build the fixture fresh (host, Float64) so host and device runs are independent.
function build()
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    nlevsno = CLM.varpar.nlevsno; nlevgrnd = CLM.varpar.nlevgrnd
    ncol_lyr = nlevsno + nlevgrnd
    nc = 3
    ISTICE = CLM.ISTICE; ISTSOIL = CLM.ISTSOIL
    # col1 istice w/ meltwater in 2 layers, col2 istsoil (persistent snow), col3 istice no melt
    col_landunit = [1, 2, 3]
    col_gridcell = [1, 1, 1]
    lun_itype    = [ISTICE, ISTSOIL, ISTICE]
    mask_do_smb  = Bool[true, true, true]
    mask_allc    = Bool[true, true, true]
    h2osoi_liq = zeros(nc, ncol_lyr); h2osoi_ice = zeros(nc, ncol_lyr)
    j1 = 1 + nlevsno; j2 = 2 + nlevsno
    h2osoi_liq[1, j1] = 5.0;  h2osoi_ice[1, j1] = 100.0
    h2osoi_liq[1, j2] = 3.0;  h2osoi_ice[1, j2] = 50.0
    h2osoi_liq[2, j1] = 7.0                     # non-istice: untouched
    max_days = 7300; thr = max_days * CLM.SECSPDAY
    snow_persistence = [0.0, thr + 1.0, 0.0]
    qflx_snwcp_ice   = [0.4, 0.5, 0.6]
    qflx_glcice_melt = fill(NaN, nc)
    qflx_glcice      = fill(NaN, nc)
    qflx_glcice_frz  = fill(NaN, nc)
    qflx_dyn         = fill(NaN, nc)
    glc_dyn_routing  = [0.7]
    qflx_qrgwl       = [0.1, 0.2, 0.3]
    qflx_ice_runoff  = [2.0, 3.0, 4.0]
    return (; nc, nlevsno, nlevgrnd, max_days, col_landunit, col_gridcell, lun_itype,
            mask_do_smb, mask_allc, h2osoi_liq, h2osoi_ice, snow_persistence,
            qflx_snwcp_ice, qflx_glcice_melt, qflx_glcice, qflx_glcice_frz, qflx_dyn,
            glc_dyn_routing, qflx_qrgwl, qflx_ice_runoff)
end

# Run all three routines in driver order over a fixture (host arrays or device dict).
function run_smb!(f, dev)
    id(x) = dev ? mf(x) : x
    nc = f.nc; bounds = 1:nc; dtime = 1800.0
    liq = id(f.h2osoi_liq); ice = id(f.h2osoi_ice); melt = id(f.qflx_glcice_melt)
    cl = id(f.col_landunit); cg = id(f.col_gridcell); li = id(f.lun_itype)
    smb = id(f.mask_do_smb); allc = id(f.mask_allc)
    snwcp = id(f.qflx_snwcp_ice); persist = id(f.snow_persistence)
    glc = id(f.qflx_glcice); frz = id(f.qflx_glcice_frz); dyn = id(f.qflx_dyn)
    routing = id(f.glc_dyn_routing); qrgwl = id(f.qflx_qrgwl); icerun = id(f.qflx_ice_runoff)

    CLM.handle_ice_melt!(liq, ice, melt, cl, li, smb, bounds, dtime, f.nlevsno, f.nlevgrnd)
    CLM.compute_surface_mass_balance!(glc, frz, dyn, snwcp, melt, persist, routing,
        cl, cg, li, allc, smb, bounds; glc_snow_persistence_max_days=f.max_days)
    CLM.adjust_runoff_terms!(qrgwl, icerun, frz, melt, routing, cg, smb, bounds)
    dev && device_synchronize()
    return (; melt, ice, liq, glc, frz, dyn, qrgwl, icerun)
end

function main(backend)
    println("="^64); println("Glacier SMB — host vs device parity"); println("="^64)
    if backend === nothing; println("  No GPU backend."); return 0; end
    name, _, FT = backend
    @printf("  Backend: %s (%s)\n\n", name, FT)

    H = run_smb!(build(), false)
    D = run_smb!(build(), true)

    fields = (:melt, :ice, :liq, :glc, :frz, :dyn, :qrgwl, :icerun)
    nfail = 0; ncmp = 0
    for nm in fields
        r, n = reldiff(getfield(H, nm), getfield(D, nm)); ncmp += n
        ok = r < 1f-4
        @printf("  [%s] %-8s rel=%.2e over %d finite\n", ok ? "PASS" : "FAIL", nm, r, n)
        ok || (nfail += 1)
    end
    println()
    println(nfail == 0 ? "  GLACIER SMB kernels MATCH host on $name over $ncmp finite outputs" :
                         "  DIVERGENCE ($nfail fields).")
    return nfail == 0 ? 0 : 1
end
exit(main(detect_backend()))
