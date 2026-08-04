# Setups shared by the default-bounds conservation checks. Kept beside
# run_default_bounds.jl rather than in test/ proper so it is obvious these are
# not part of the `Pkg.test` suite (which forces --check-bounds=yes and is
# therefore blind to what they test — see the header of run_default_bounds.jl).

using CLM

"""
    _lvt_conservation_case(nc, nlevdecomp; ndecomp_pools = 7)

Drive `litter_vert_transp!` on a synthetic column set and return the worst
violation of its exact mass identity.

The Patankar discretisation with the zero-flux top (a=0, b=1, c=-1, r=0) and
bottom (a=-1, b=1, c=0, r=0) boundary rows this routine builds conserves mass
EXACTLY, so for every non-CWD pool and every column

    Σⱼ (conc_after[j] − conc_before[j])·dz[j]  ==  Σⱼ source[j]·dz[j]

Advective leakage is floored at `epsilon = 1e-30`, so there is no legitimate
residual: anything above roundoff means the in-thread Thomas sweep was
miscompiled.
"""
function _lvt_conservation_case(nc::Int, nlevdecomp::Int; ndecomp_pools::Int = 7)
    scalez = 0.025
    nlevgrnd = nlevdecomp + 5
    zsoi = [scalez * (exp(0.5 * (j - 0.5)) - 1.0) for j in 1:nlevgrnd]
    zisoi = zeros(nlevgrnd + 1)
    for j in 1:nlevgrnd
        zisoi[j+1] = j < nlevgrnd ? 0.5 * (zsoi[j] + zsoi[j+1]) :
                                    zsoi[j] + 0.5 * (zsoi[j] - zisoi[j])
    end
    dzs = [zisoi[j+1] - zisoi[j] for j in 1:nlevdecomp]

    cs = CLM.SoilBiogeochemCarbonStateData()
    cs.decomp_cpools_vr_col =
        [100.0 * exp(-0.5 * j) * (1 + 0.1 * s) for c in 1:nc, j in 1:nlevdecomp, s in 1:ndecomp_pools]
    cf = CLM.SoilBiogeochemCarbonFluxData()
    cf.decomp_cpools_sourcesink_col = [0.01 * j * s for c in 1:nc, j in 1:nlevdecomp, s in 1:ndecomp_pools]
    cf.decomp_cpools_transport_tendency_col = zeros(nc, nlevdecomp, ndecomp_pools)
    cf.tri_ma_vr = zeros(nc, 1)

    ns = CLM.SoilBiogeochemNitrogenStateData()
    ns.decomp_npools_vr_col = cs.decomp_cpools_vr_col ./ 15.0
    nf = CLM.SoilBiogeochemNitrogenFluxData()
    nf.decomp_npools_sourcesink_col = zeros(nc, nlevdecomp, ndecomp_pools)
    nf.decomp_npools_transport_tendency_col = zeros(nc, nlevdecomp, ndecomp_pools)

    st = CLM.SoilBiogeochemStateData()
    st.som_adv_coef_col = zeros(nc, nlevdecomp + 1)
    st.som_diffus_coef_col = zeros(nc, nlevdecomp + 1)

    col = CLM.ColumnData()
    col.nbedrock = fill(nlevdecomp, nc); col.gridcell = ones(Int, nc)
    grc = CLM.GridcellData(); grc.latdeg = fill(45.0, 1)
    al = CLM.ActiveLayerData()
    al.altmax_col = fill(3.0, nc); al.altmax_lastyear_col = fill(3.0, nc)

    cc = CLM.DecompCascadeConData()
    is_cwd = falses(ndecomp_pools); is_cwd[5] = true
    cc.is_cwd = is_cwd; cc.spinup_factor = ones(ndecomp_pools)

    pr = CLM.LitterVertTranspParams()
    pr.som_diffus = 1.0e-4
    pr.cryoturb_diffusion_k = 5.0e-4
    pr.max_altdepth_cryoturbation = 2.0

    c0  = copy(cs.decomp_cpools_vr_col)
    src = copy(cf.decomp_cpools_sourcesink_col)

    CLM.litter_vert_transp!(cs, cf, ns, nf, st, col, grc, cc, al, pr;
        mask_bgc_soilc = trues(nc), bounds = 1:nc, dtime = 1800.0,
        nlevdecomp = nlevdecomp, ndecomp_pools = ndecomp_pools,
        zsoi_vals = zsoi, dzsoi_decomp_vals = dzs, zisoi_vals = zisoi,
        spinup_state = 0, use_soil_matrixcn = false,
        som_adv_flux_val = 0.0, max_depth_cryoturb_val = 3.0)

    max_abs = 0.0; max_rel = 0.0
    for s in 1:ndecomp_pools
        is_cwd[s] && continue
        for c in 1:nc
            dmass = sum((cs.decomp_cpools_vr_col[c, j, s] - c0[c, j, s]) * dzs[j] for j in 1:nlevdecomp)
            smass = sum(src[c, j, s] * dzs[j] for j in 1:nlevdecomp)
            err = abs(dmass - smass)
            max_abs = max(max_abs, err)
            max_rel = max(max_rel, err / max(abs(smass), 1e-30))
        end
    end
    return (max_abs_err = max_abs, max_rel_err = max_rel)
end

"""
    _snicar_host_vs_device(flg_slr_in; ncols = 12, nlevsno = 5)

Run the SNICAR radiative transfer through BOTH implementations and return their
worst disagreement.

`snicar_rt!` is a plain per-column host loop, not a `@kernel`, so it is an
independent reference for `snicar_rt_device!`. The two agreed bit-exactly under
`--check-bounds=yes` even while the device kernel's adding-doubling sweeps were
being miscompiled — which is precisely what proved the transliteration faithful
and the optimizer responsible. Five resolved snow layers (snl = −5) so both
sweeps run full depth.
"""
function _snicar_host_vs_device(flg_slr_in::Int; ncols::Int = 12, nlevsno::Int = 5)
    numrad_snw = 5
    od = CLM.SnicarOpticsData()
    CLM.snow_optics_init!(od; numrad_snw = numrad_snw)
    fill!(od.ss_alb_snw_drc, 0.9999); fill!(od.asm_prm_snw_drc, 0.85)
    fill!(od.ext_cff_mss_snw_drc, 10.0)
    fill!(od.ss_alb_snw_dfs, 0.9999); fill!(od.asm_prm_snw_dfs, 0.85)
    fill!(od.ext_cff_mss_snw_dfs, 10.0)
    od.flx_wgt_dir .= [1.0, 0.5, 0.3, 0.15, 0.05]
    od.flx_wgt_dif .= [1.0, 0.5, 0.3, 0.15, 0.05]

    old_numrad = CLM.varctl.snicar_numrad_snw
    old_shape  = CLM.varctl.snicar_snw_shape
    CLM.varctl.snicar_numrad_snw = numrad_snw
    CLM.varctl.snicar_snw_shape  = "sphere"
    try
        coszen       = fill(0.5, ncols)
        h2osno_liq   = [1.0 + 0.1 * j for c in 1:ncols, j in 1:nlevsno]
        h2osno_ice   = [20.0 + 2.0 * j for c in 1:ncols, j in 1:nlevsno]
        h2osno_total = [sum(h2osno_liq[c, :]) + sum(h2osno_ice[c, :]) for c in 1:ncols]
        snw_rds      = fill(100, ncols, nlevsno)
        mss_cnc_aer  = fill(1e-8, ncols, nlevsno, CLM.SNO_NBR_AER)
        albsfc       = fill(0.2, ncols, CLM.NUMRAD)
        snl_vec      = fill(-5, ncols)
        frac_sno     = fill(1.0, ncols)

        ah = zeros(ncols, CLM.NUMRAD); fh = zeros(ncols, nlevsno + 1, CLM.NUMRAD)
        ad = zeros(ncols, CLM.NUMRAD); fd = zeros(ncols, nlevsno + 1, CLM.NUMRAD)

        CLM.snicar_rt!(coszen, flg_slr_in, h2osno_liq, h2osno_ice, h2osno_total,
            snw_rds, mss_cnc_aer, albsfc, snl_vec, frac_sno, ah, fh, nlevsno; optics = od)
        CLM.snicar_rt_device!(coszen, flg_slr_in, h2osno_liq, h2osno_ice, h2osno_total,
            snw_rds, mss_cnc_aer, albsfc, snl_vec, frac_sno, ad, fd, nlevsno; optics = od)

        return (max_albedo_diff = maximum(abs, ah .- ad),
                max_flx_abs_diff = maximum(abs, fh .- fd))
    finally
        CLM.varctl.snicar_numrad_snw = old_numrad
        CLM.varctl.snicar_snw_shape  = old_shape
    end
end
