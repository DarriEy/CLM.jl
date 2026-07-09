# ==========================================================================
# gpu_validate_dyn_cnbal_patch_e2e.jl — WHOLE-FUNCTION Metal parity for the
# dyn_cnbal_patch! ORCHESTRATOR (the transient-land-use patch C/N balance).
#
# The sub-functions dynamic_patch_adjustments_carbon!/_nitrogen! (see
# gpu_validate_dyn_cnbal_patch_carbon_e2e.jl) and dyn_cnbal_col! already have
# Metal harnesses; this covers the orchestrator's own loops that were the last
# host-only piece: the dwt/re-init kernel, the sign-flips, and the patch->gridcell
# / patch->column atomic scatters (seeding, litter/CWD, product, conversion, slash).
#
# Mirrors test_dyn_cons_biogeochem's 3-patch subgrid (shrink / grow / initiate),
# runs on CPU, adapts every device-resident struct to Metal/Float32, runs the SAME
# call on-device, and compares the mutated flux + re-init + veg-pool fields.
#
#   julia --project=scripts scripts/gpu_validate_dyn_cnbal_patch_e2e.jl
# ==========================================================================

using CLM
using Printf
import Metal
include(joinpath(@__DIR__, "gpu_backends.jl"))
include(joinpath(@__DIR__, "gpu_adapt.jl"))   # shared mf: keeps pinned ::Float64 scalars (leaf_mr_vcm)

reldiff(a, b) = begin
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + max(abs(Float64(A[i])), abs(Float64(B[i])))))
    end
    m
end

const NP = 3
const NC = 3
const NG = 1
const PWT_OLD = [0.40, 0.10, 0.00]
const PWT_NEW = [0.20, 0.30, 0.25]
const NDECOMP_POOLS = 3
const I_LITR_MIN = 1
const I_LITR_MAX = 2

function make_pftcon(n, ndecomp_pools)
    p = CLM.PftconType()
    p.c3psn = ones(n); p.evergreen = zeros(n); p.woody = zeros(n)
    p.leafcn = fill(25.0, n); p.deadwdcn = fill(100.0, n)
    p.noveg = fill(typemax(Int), n); p.pconv = fill(0.5, n)
    p.fr_f = zeros(n, ndecomp_pools)
    for j in 1:n; p.fr_f[j,1]=0.6; p.fr_f[j,2]=0.4; end
    return p
end

# Build the full host fixture; returns the named tuple of everything dyn_cnbal_patch!
# needs (+ the mask / prior weights). Mirrors test_dyn_cons_biogeochem.
function build()
    CLM.varpar_init!(CLM.varpar, 1, 15, 2, 5); CLM.varcon_init!()
    nlevdecomp = CLM.varpar.nlevdecomp
    bounds = CLM.BoundsType(begg=1, endg=NG, begl=1, endl=1, begc=1, endc=NC,
                            begp=1, endp=NP, level=CLM.BOUNDS_LEVEL_PROC)
    pch = CLM.PatchData{Float64}()
    pch.column = collect(1:NP); pch.landunit = fill(1, NP); pch.gridcell = fill(1, NP)
    pch.itype = [7, 13, 8]; pch.wtgcell = copy(PWT_OLD)
    col = CLM.ColumnData{Float64}(); col.gridcell = fill(1, NC); col.wtgcell = copy(PWT_OLD)
    lun = CLM.LandunitData(); lun.itype = [CLM.ISTSOIL]

    updater = CLM.PatchStateUpdater(bounds)
    CLM.set_old_weights!(updater, bounds, pch, col)
    pch.wtgcell .= PWT_NEW; col.wtgcell .= PWT_NEW
    CLM.set_new_weights!(updater, bounds, pch)

    cs = CLM.CNVegCarbonStateData{Float64}(); CLM.cnveg_carbon_state_init!(cs, NP, NC, NG)
    ns = CLM.CNVegNitrogenStateData{Float64}(); CLM.cnveg_nitrogen_state_init!(ns, NP, NC, NG)
    cpool_fields = (:leafc_patch, :leafc_storage_patch, :leafc_xfer_patch,
        :frootc_patch, :frootc_storage_patch, :frootc_xfer_patch,
        :livestemc_patch, :livestemc_storage_patch, :livestemc_xfer_patch,
        :deadstemc_patch, :deadstemc_storage_patch, :deadstemc_xfer_patch,
        :livecrootc_patch, :livecrootc_storage_patch, :livecrootc_xfer_patch,
        :deadcrootc_patch, :deadcrootc_storage_patch, :deadcrootc_xfer_patch,
        :gresp_storage_patch, :gresp_xfer_patch, :cpool_patch, :xsmrpool_patch, :ctrunc_patch)
    npool_fields = (:leafn_patch, :leafn_storage_patch, :leafn_xfer_patch,
        :frootn_patch, :frootn_storage_patch, :frootn_xfer_patch,
        :livestemn_patch, :livestemn_storage_patch, :livestemn_xfer_patch,
        :deadstemn_patch, :deadstemn_storage_patch, :deadstemn_xfer_patch,
        :livecrootn_patch, :livecrootn_storage_patch, :livecrootn_xfer_patch,
        :deadcrootn_patch, :deadcrootn_storage_patch, :deadcrootn_xfer_patch,
        :retransn_patch, :npool_patch, :ntrunc_patch)
    v = 1.0
    for f in cpool_fields; arr = getproperty(cs, f); for p in 1:NP; arr[p] = v; v += 1.0; end; end
    for f in npool_fields; arr = getproperty(ns, f); for p in 1:NP; arr[p] = v*0.01; v += 1.0; end; end

    cf = CLM.CNVegCarbonFluxData{Float64}()
    CLM.cnveg_carbon_flux_init!(cf, NP, NC, NG; nlevdecomp_full=nlevdecomp, ndecomp_pools=NDECOMP_POOLS)
    nf = CLM.CNVegNitrogenFluxData{Float64}()
    CLM.cnveg_nitrogen_flux_init!(nf, NP, NC, NG; nlevdecomp_full=nlevdecomp, i_litr_max=I_LITR_MAX)
    for fld in (:dwt_seedc_to_leaf_patch, :dwt_seedc_to_leaf_grc,
                :dwt_seedc_to_deadstem_patch, :dwt_seedc_to_deadstem_grc,
                :dwt_conv_cflux_patch, :dwt_conv_cflux_grc,
                :dwt_wood_productc_gain_patch, :dwt_crop_productc_gain_patch,
                :dwt_slash_cflux_patch, :dwt_slash_cflux_grc)
        getproperty(cf, fld) .= 0.0
    end
    cf.dwt_frootc_to_litr_c_col .= 0.0; cf.dwt_livecrootc_to_cwdc_col .= 0.0; cf.dwt_deadcrootc_to_cwdc_col .= 0.0
    for fld in (:dwt_seedn_to_leaf_patch, :dwt_seedn_to_leaf_grc,
                :dwt_seedn_to_deadstem_patch, :dwt_seedn_to_deadstem_grc,
                :dwt_conv_nflux_patch, :dwt_conv_nflux_grc,
                :dwt_wood_productn_gain_patch, :dwt_crop_productn_gain_patch)
        getproperty(nf, fld) .= 0.0
    end
    nf.dwt_frootn_to_litr_n_col .= 0.0; nf.dwt_livecrootn_to_cwdn_col .= 0.0; nf.dwt_deadcrootn_to_cwdn_col .= 0.0

    cnst = CLM.CNVegStateData{Float64}(); CLM.cnveg_state_init!(cnst, NP, NC)
    cnst.dwt_smoothed_patch .= 0.0; cnst.annavg_t2m_col .= 285.0
    canopy = CLM.CanopyStateData{Float64}(); CLM.canopystate_init!(canopy, NP)
    sbgc = CLM.SoilBiogeochemStateData{Float64}(); CLM.soil_bgc_state_init!(sbgc, NC, NP, nlevdecomp, 0)
    sbgc.froot_prof_patch .= 1.0 / nlevdecomp; sbgc.croot_prof_patch .= 1.0 / nlevdecomp

    return (dynbal=CLM.DynConsBiogeochemState(), bounds=bounds, updater=updater,
            pch=pch, lun=lun, col=col, cs=cs, cf=cf, ns=ns, nf=nf, cnst=cnst,
            canopy=canopy, sbgc=sbgc)
end

run!(f, pftcon, mask) = CLM.dyn_cnbal_patch!(f.dynbal, f.bounds, mask, PWT_OLD, f.updater,
    f.pch, f.lun, f.col, pftcon, f.canopy, f.cnst, f.cs, f.cf, f.ns, f.nf, f.sbgc;
    dt=1800.0, i_litr_min=I_LITR_MIN, i_litr_max=I_LITR_MAX)

function main(backend)
    println("=" ^ 72)
    println("WHOLE-FN Metal parity for dyn_cnbal_patch! (transient land-use orchestrator)")
    println("=" ^ 72)
    if backend === nothing; println("  No GPU backend — nothing to validate."); return 0; end
    name, _, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n", name, FT)
    pftcon = make_pftcon(15, NDECOMP_POOLS)

    # CPU reference
    H = build(); run!(H, pftcon, trues(NP))

    # Device
    dev(x) = mf(Metal.MtlArray, x)
    D = build()
    D = (dynbal=D.dynbal, bounds=D.bounds, updater=dev(D.updater),
         pch=dev(D.pch), lun=dev(D.lun), col=dev(D.col), cs=dev(D.cs), cf=dev(D.cf),
         ns=dev(D.ns), nf=dev(D.nf), cnst=dev(D.cnst), canopy=dev(D.canopy), sbgc=dev(D.sbgc))
    if !(D.cs.leafc_patch isa Metal.MtlArray)
        println("  BLOCKED: CN carbon state did not move to the device."); return 2
    end
    run!(D, pftcon, trues(NP)); Metal.synchronize()

    # compare the mutated fields
    checks = [
        ("cf.dwt_seedc_to_leaf_patch", H.cf.dwt_seedc_to_leaf_patch, D.cf.dwt_seedc_to_leaf_patch),
        ("cf.dwt_seedc_to_leaf_grc", H.cf.dwt_seedc_to_leaf_grc, D.cf.dwt_seedc_to_leaf_grc),
        ("cf.dwt_seedc_to_deadstem_grc", H.cf.dwt_seedc_to_deadstem_grc, D.cf.dwt_seedc_to_deadstem_grc),
        ("cf.dwt_conv_cflux_patch", H.cf.dwt_conv_cflux_patch, D.cf.dwt_conv_cflux_patch),
        ("cf.dwt_conv_cflux_grc", H.cf.dwt_conv_cflux_grc, D.cf.dwt_conv_cflux_grc),
        ("cf.dwt_wood_productc_gain_patch", H.cf.dwt_wood_productc_gain_patch, D.cf.dwt_wood_productc_gain_patch),
        ("cf.dwt_slash_cflux_patch", H.cf.dwt_slash_cflux_patch, D.cf.dwt_slash_cflux_patch),
        ("cf.dwt_slash_cflux_grc", H.cf.dwt_slash_cflux_grc, D.cf.dwt_slash_cflux_grc),
        ("cf.dwt_frootc_to_litr_c_col", H.cf.dwt_frootc_to_litr_c_col, D.cf.dwt_frootc_to_litr_c_col),
        ("cf.dwt_livecrootc_to_cwdc_col", H.cf.dwt_livecrootc_to_cwdc_col, D.cf.dwt_livecrootc_to_cwdc_col),
        ("cf.dwt_deadcrootc_to_cwdc_col", H.cf.dwt_deadcrootc_to_cwdc_col, D.cf.dwt_deadcrootc_to_cwdc_col),
        ("nf.dwt_seedn_to_leaf_grc", H.nf.dwt_seedn_to_leaf_grc, D.nf.dwt_seedn_to_leaf_grc),
        ("nf.dwt_conv_nflux_grc", H.nf.dwt_conv_nflux_grc, D.nf.dwt_conv_nflux_grc),
        ("nf.dwt_frootn_to_litr_n_col", H.nf.dwt_frootn_to_litr_n_col, D.nf.dwt_frootn_to_litr_n_col),
        ("nf.dwt_deadcrootn_to_cwdn_col", H.nf.dwt_deadcrootn_to_cwdn_col, D.nf.dwt_deadcrootn_to_cwdn_col),
        # re-init fields (initiating patch p=3)
        ("cnst.dwt_smoothed_patch", H.cnst.dwt_smoothed_patch, D.cnst.dwt_smoothed_patch),
        ("cnst.dormant_flag_patch", H.cnst.dormant_flag_patch, D.cnst.dormant_flag_patch),
        ("cnst.annavg_t2m_patch", H.cnst.annavg_t2m_patch, D.cnst.annavg_t2m_patch),
        ("cf.availc_patch", H.cf.availc_patch, D.cf.availc_patch),
        ("canopy.laisun_patch", H.canopy.laisun_patch, D.canopy.laisun_patch),
        # veg pools (mutated by the adjustment subfns)
        ("cs.leafc_patch", H.cs.leafc_patch, D.cs.leafc_patch),
        ("cs.deadstemc_patch", H.cs.deadstemc_patch, D.cs.deadstemc_patch),
        ("ns.leafn_patch", H.ns.leafn_patch, D.ns.leafn_patch),
    ]
    nfail = 0; ncmp = 0
    for (nm, h, d) in checks
        dd = reldiff(h, d); ok = dd < 1f-3; ncmp += 1
        @printf("  [%s] %-32s rel = %.3e\n", ok ? "PASS" : "FAIL", nm, dd)
        ok || (nfail += 1)
    end
    println()
    println(nfail == 0 ? "  WHOLE dyn_cnbal_patch! MATCHES CPU ON $name over $ncmp fields" :
                         "  DIVERGENCE — investigate ($nfail failed).")
    return nfail == 0 ? 0 : 1
end

exit(main(detect_backend()))
