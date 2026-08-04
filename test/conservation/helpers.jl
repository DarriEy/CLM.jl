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

"""
    _ch4_tran_case()

Drive `ch4_tran!` on a synthetic saturated column set (water table at level 7)
and return its own CH4 balance residual plus the resulting concentrations.

`jwt = 7` matters: ebullition is added to the CH4 source at the water-table
layer, so that is where a corrupted `ch4_ebul_total` reduction enters the
tracer solve.
"""
function _ch4_tran_case()
    nc = 12; np = 12; ng = 2; nlevsoi = 15
    params = CLM.CH4Params(); ch4vc = CLM.CH4VarCon()
    z2(a...) = zeros(a...)
    ch4 = CLM.CH4Data(
        ch4_prod_depth_sat_col = fill(1e-8, nc, nlevsoi),
        ch4_prod_depth_unsat_col = fill(1e-9, nc, nlevsoi),
        ch4_prod_depth_lake_col = z2(nc, nlevsoi),
        ch4_oxid_depth_sat_col = fill(1e-9, nc, nlevsoi),
        ch4_oxid_depth_unsat_col = fill(1e-10, nc, nlevsoi),
        ch4_oxid_depth_lake_col = z2(nc, nlevsoi),
        ch4_aere_depth_sat_col = fill(1e-10, nc, nlevsoi),
        ch4_aere_depth_unsat_col = fill(1e-10, nc, nlevsoi),
        ch4_tran_depth_sat_col = z2(nc, nlevsoi), ch4_tran_depth_unsat_col = z2(nc, nlevsoi),
        ch4_ebul_depth_sat_col = fill(1e-10, nc, nlevsoi), ch4_ebul_depth_unsat_col = z2(nc, nlevsoi),
        o2_oxid_depth_sat_col = fill(2e-9, nc, nlevsoi), o2_oxid_depth_unsat_col = fill(2e-10, nc, nlevsoi),
        o2_aere_depth_sat_col = fill(1e-9, nc, nlevsoi), o2_aere_depth_unsat_col = fill(1e-9, nc, nlevsoi),
        co2_decomp_depth_sat_col = z2(nc, nlevsoi), co2_decomp_depth_unsat_col = z2(nc, nlevsoi),
        co2_oxid_depth_sat_col = z2(nc, nlevsoi), co2_oxid_depth_unsat_col = z2(nc, nlevsoi),
        co2_aere_depth_sat_col = z2(nc, nlevsoi), co2_aere_depth_unsat_col = z2(nc, nlevsoi),
        ch4_ebul_total_sat_col = z2(nc), ch4_ebul_total_unsat_col = z2(nc),
        ch4_surf_aere_sat_col = z2(nc), ch4_surf_aere_unsat_col = z2(nc),
        ch4_surf_ebul_sat_col = z2(nc), ch4_surf_ebul_unsat_col = z2(nc), ch4_surf_ebul_lake_col = z2(nc),
        ch4_surf_diff_sat_col = z2(nc), ch4_surf_diff_unsat_col = z2(nc), ch4_surf_diff_lake_col = z2(nc),
        ch4_dfsat_flux_col = z2(nc), ch4_surf_flux_tot_col = z2(nc),
        conc_ch4_sat_col = fill(1.0e-4, nc, nlevsoi), conc_ch4_unsat_col = fill(1.0e-5, nc, nlevsoi),
        conc_ch4_lake_col = z2(nc, nlevsoi),
        conc_o2_sat_col = fill(0.01, nc, nlevsoi), conc_o2_unsat_col = fill(0.05, nc, nlevsoi),
        conc_o2_lake_col = z2(nc, nlevsoi),
        o2_decomp_depth_sat_col = fill(1e-9, nc, nlevsoi), o2_decomp_depth_unsat_col = fill(1e-9, nc, nlevsoi),
        o2stress_sat_col = ones(nc, nlevsoi), o2stress_unsat_col = ones(nc, nlevsoi),
        ch4stress_sat_col = ones(nc, nlevsoi), ch4stress_unsat_col = ones(nc, nlevsoi),
        zwt_ch4_unsat_col = z2(nc), lake_soilc_col = fill(100.0, nc, nlevsoi),
        totcolch4_col = z2(nc), totcolch4_bef_col = z2(nc), annsum_counter_col = z2(nc),
        tempavg_somhr_col = z2(nc), annavg_somhr_col = fill(1.0e-6, nc),
        tempavg_finrw_col = z2(nc), annavg_finrw_col = fill(0.1, nc), sif_col = ones(nc),
        qflx_surf_lag_col = z2(nc), finundated_col = fill(0.1, nc),
        finundated_pre_snow_col = fill(0.1, nc), finundated_lag_col = fill(0.1, nc),
        layer_sat_lag_col = fill(0.5, nc, nlevsoi), pH_col = fill(6.5, nc),
        c_atm_grc = fill(0.03, ng, 3), ch4co2f_grc = z2(ng), ch4prodg_grc = z2(ng),
        totcolch4_grc = z2(ng), totcolch4_bef_grc = z2(ng),
        annavg_agnpp_patch = fill(1.0e-5, np), annavg_bgnpp_patch = fill(1.0e-5, np),
        tempavg_agnpp_patch = fill(1.0e-6, np), tempavg_bgnpp_patch = fill(1.0e-6, np),
        grnd_ch4_cond_patch = fill(0.01, np), grnd_ch4_cond_col = fill(0.01, nc),
        ch4_first_time_grc = fill(true, ng))

    watsat = fill(0.45, nc, nlevsoi); h2osoi_vol = fill(0.30, nc, nlevsoi)
    h2osoi_liq = fill(30.0, nc, nlevsoi); h2osoi_ice = zeros(nc, nlevsoi)
    bsw = fill(5.0, nc, nlevsoi); cellorg = fill(10.0, nc, nlevsoi)
    t_soisno = fill(CLM.TFRZ + 15.0, nc, nlevsoi)
    t_grnd = fill(CLM.TFRZ + 15.0, nc); t_h2osfc = fill(CLM.TFRZ + 15.0, nc)
    frac_h2osfc = fill(0.1, nc); h2osfc = fill(1.0, nc); snow_depth = zeros(nc)
    snl = zeros(Int, nc)
    dz = fill(0.1, nc, nlevsoi); z = [(j-0.5)*0.1 for c in 1:nc, j in 1:nlevsoi]
    zi = [j*0.1 for c in 1:nc, j in 1:nlevsoi]
    jwt = fill(7, nc)
    col_gridcell = [1 + (c-1) % ng for c in 1:nc]

    conc_ch4_0 = copy(ch4.conc_ch4_sat_col)

    CLM.ch4_tran!(ch4, params, ch4vc, trues(nc), col_gridcell, watsat, h2osoi_vol,
        h2osoi_liq, h2osoi_ice, h2osfc, bsw, cellorg, t_soisno, t_grnd, t_h2osfc,
        frac_h2osfc, snow_depth, snl, z, dz, zi, jwt, 1, false, nlevsoi, 5, 1800.0, 130.0)

    dtime = 1800.0
    worst = 0.0
    for c in 1:nc
        e = 0.0
        for j in 1:nlevsoi
            e += (ch4.conc_ch4_sat_col[c,j] - conc_ch4_0[c,j]) * dz[c,j]
            e -= ch4.ch4_prod_depth_sat_col[c,j] * dz[c,j] * dtime
            e += ch4.ch4_oxid_depth_sat_col[c,j] * dz[c,j] * dtime
        end
        e += (ch4.ch4_surf_aere_sat_col[c] + ch4.ch4_surf_ebul_sat_col[c] +
              ch4.ch4_surf_diff_sat_col[c]) * dtime
        worst = max(worst, abs(e))
    end
    return (max_errch4 = worst,
            conc_ch4 = copy(ch4.conc_ch4_sat_col),
            conc_o2 = copy(ch4.conc_o2_sat_col))
end
