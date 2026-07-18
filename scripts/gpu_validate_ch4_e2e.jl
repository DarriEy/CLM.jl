# ==========================================================================
# gpu_validate_ch4_e2e.jl — end-to-end Metal parity for the WHOLE ch4! methane
# orchestrator: its per-column/per-(c,j) pre/post loops + the grnd_ch4_cond
# patch->column and the column->gridcell flux scatters are kernels; the sat/lake
# dispatch calls the 6 device process kernels; conservation @warn checks run host.
#
# Clones test/test_methane.jl's "ch4! main driver" setup at Float64, runs ch4! on
# the CPU, adapts the whole CH4Data + every forcing array to Float32/Metal, runs the
# SAME call on the device, and compares the mutated column/gridcell outputs.
#
# Covers TWO finundation configurations:
#   (1) h2osfc      (CTSM default; _ch4o_finund_kernel!)
#   (2) TWS_inversion (#238; _ch4o_finund_inv_kernel! — the satellite-regression
#       CalcFinundated device kernel, previously NOT Metal-validated: the default
#       ch4! harness only ever exercised the h2osfc route where finundated≡frac_h2osfc).
#
#   julia --project=scripts scripts/gpu_validate_ch4_e2e.jl
# ==========================================================================

using CLM
using Printf
import Metal
include(joinpath(@__DIR__, "gpu_backends.jl"))

function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + max(abs(Float64(A[i])), abs(Float64(B[i])))))
    end
    return m
end
allfinite(a) = all(isfinite, Array(a))

# Float32 down-convert adaptor for the whole CH4Data + arrays.
struct MF32 end
CLM.Adapt.adapt_storage(::MF32, x::AbstractArray{<:AbstractFloat}) = Metal.MtlArray(Float32.(x))
CLM.Adapt.adapt_storage(::MF32, x::AbstractArray{<:Integer})       = Metal.MtlArray(collect(Int, x))
CLM.Adapt.adapt_storage(::MF32, x::AbstractArray{Bool})            = Metal.MtlArray(collect(Bool, x))

function build(::Type{FT}; finmtd=CLM.FINUNDATION_MTD_H2OSFC) where {FT}
    nc=3; np=4; ng=2; nlevsoi=5
    M(v) = fill(FT(v), nc, nlevsoi); Mp(v) = fill(FT(v), np, nlevsoi)
    V(v) = fill(FT(v), nc); Vp(v) = fill(FT(v), np); Vg(v) = fill(FT(v), ng)
    ch4 = CLM.CH4Data{FT}(
        ch4_prod_depth_sat_col=M(0), ch4_prod_depth_unsat_col=M(0), ch4_prod_depth_lake_col=M(0),
        ch4_oxid_depth_sat_col=M(0), ch4_oxid_depth_unsat_col=M(0), ch4_oxid_depth_lake_col=M(0),
        ch4_aere_depth_sat_col=M(0), ch4_aere_depth_unsat_col=M(0),
        ch4_tran_depth_sat_col=M(0), ch4_tran_depth_unsat_col=M(0),
        ch4_ebul_depth_sat_col=M(0), ch4_ebul_depth_unsat_col=M(0),
        o2_oxid_depth_sat_col=M(0), o2_oxid_depth_unsat_col=M(0),
        o2_aere_depth_sat_col=M(0), o2_aere_depth_unsat_col=M(0),
        co2_decomp_depth_sat_col=M(0), co2_decomp_depth_unsat_col=M(0),
        co2_oxid_depth_sat_col=M(0), co2_oxid_depth_unsat_col=M(0),
        co2_aere_depth_sat_col=M(0), co2_aere_depth_unsat_col=M(0),
        ch4_ebul_total_sat_col=V(0), ch4_ebul_total_unsat_col=V(0),
        ch4_surf_aere_sat_col=V(0), ch4_surf_aere_unsat_col=V(0),
        ch4_surf_ebul_sat_col=V(0), ch4_surf_ebul_unsat_col=V(0), ch4_surf_ebul_lake_col=V(0),
        ch4_surf_diff_sat_col=V(0), ch4_surf_diff_unsat_col=V(0), ch4_surf_diff_lake_col=V(0),
        ch4_dfsat_flux_col=V(0), ch4_surf_flux_tot_col=V(0),
        conc_ch4_sat_col=M(1.0e-4), conc_ch4_unsat_col=M(1.0e-5), conc_ch4_lake_col=M(0),
        conc_o2_sat_col=M(0.01), conc_o2_unsat_col=M(0.05), conc_o2_lake_col=M(0),
        o2_decomp_depth_sat_col=M(0), o2_decomp_depth_unsat_col=M(0),
        o2stress_sat_col=M(1), o2stress_unsat_col=M(1), ch4stress_sat_col=M(1), ch4stress_unsat_col=M(1),
        zwt_ch4_unsat_col=V(0), lake_soilc_col=M(100), totcolch4_col=V(0), totcolch4_bef_col=V(0),
        annsum_counter_col=V(0), tempavg_somhr_col=V(0), annavg_somhr_col=V(1.0e-6),
        tempavg_finrw_col=V(0), annavg_finrw_col=V(0.1), sif_col=V(1), qflx_surf_lag_col=V(0),
        finundated_col=V(0.1), finundated_pre_snow_col=V(0.1), finundated_lag_col=V(0.1),
        layer_sat_lag_col=M(0.5), pH_col=V(6.5),
        c_atm_grc=fill(FT(0.03), ng, 3), ch4co2f_grc=Vg(0), ch4prodg_grc=Vg(0),
        totcolch4_grc=Vg(0), totcolch4_bef_grc=Vg(0),
        annavg_agnpp_patch=Vp(1.0e-5), annavg_bgnpp_patch=Vp(1.0e-5),
        tempavg_agnpp_patch=Vp(1.0e-6), tempavg_bgnpp_patch=Vp(1.0e-6),
        grnd_ch4_cond_patch=Vp(0.01), grnd_ch4_cond_col=V(0.01),
        ch4_first_time_grc=fill(true, ng))
    # TWS_inversion needs a satellite-regression stream (fws_slope*TWS + fws_intercept,
    # per gridcell) + a tws gridcell vector; the coeffs vary by gridcell so finundated
    # ends up gridcell-distinct and far from the cold-start 0.1 (a meaningful device test).
    # Stays a HOST CH4FInundatedStream (Float64); ch4! copies its coeff vectors to the
    # backend internally (the _dz/_dvec path) — exactly the production device flow.
    stream = finmtd == CLM.FINUNDATION_MTD_TWS_INVERSION ?
        CLM.CH4FInundatedStream(active=true,
            fws_slope_gdc     = FT[0.05, 0.03],
            fws_intercept_gdc = FT[0.10, 0.20]) : nothing
    tws = finmtd == CLM.FINUNDATION_MTD_TWS_INVERSION ? FT[2.0, 8.0] : nothing
    return (; ch4, params=CLM.CH4Params(), ch4vc=CLM.CH4VarCon(finundation_mtd=finmtd),
        stream, tws,
        nc, np, ng, nlevsoi,
        mask_soil=trues(nc), mask_soilp=trues(np), mask_lake=falses(nc), mask_nolake=trues(nc),
        col_gridcell=[1,1,2], col_wtgcell=V(0.5), patch_column=[1,1,2,3], patch_itype=[1,2,1,1],
        patch_wtcol=Vp(0.5), is_fates=falses(nc), latdeg=FT[45.0, 30.0],
        forc_pbot=fill(FT(101325.0), ng), forc_t=fill(FT(CLM.TFRZ + 15.0), ng),
        forc_po2=fill(FT(101325.0 * 0.209), ng), forc_pco2=fill(FT(40.0), ng), forc_pch4=fill(FT(101325.0 * 1.7e-6), ng),
        watsat=M(0.5), h2osoi_vol=M(0.4), h2osoi_liq=M(0.3 * CLM.DENH2O * 0.1), h2osoi_ice=M(0.0), h2osfc=V(0.0),
        bsw=M(5.0), cellorg=M(10.0), smp_l=M(-50.0), t_soisno=fill(FT(CLM.TFRZ + 8.0), nc, nlevsoi),
        t_grnd=V(CLM.TFRZ + 15.0), t_h2osfc=V(CLM.TFRZ + 15.0), frac_h2osfc=V(0.0), snow_depth=V(0.0), snl=zeros(Int, nc),
        qflx_surf=V(0.0), rootfr_p=Mp(0.2), rootfr_col=M(0.2), crootfr=Mp(0.2), rootr_p=Mp(0.2), elai=Vp(2.0),
        qflx_tran_veg=Vp(1.0e-3), annsum_npp=Vp(500.0), rr=Vp(1.0e-6),
        somhr=V(1.0e-6), lithr=V(1.0e-7), hr_vr=M(1.0e-8), o_scalar=M(1.0), fphr=M(1.0), pot_f_nit_vr=M(0.0),
        lake_icefrac=M(0.0), lakedepth=V(5.0), z=M(0.2), dz=M(0.1), zi=M(0.25), agnpp=Vp(1.0e-5), bgnpp=Vp(1.0e-5))
end

# The TWS_inversion coeff stream + tws stay HOST (Float64); ch4! copies them to the
# backend internally, so the same host stream feeds both the CPU and device runs.
run_ch4!(S, ch4) = CLM.ch4!(ch4, S.params, S.ch4vc,
    S.mask_soil, S.mask_soilp, S.mask_lake, S.mask_nolake,
    S.col_gridcell, S.col_wtgcell, S.patch_column, S.patch_itype, S.patch_wtcol, S.is_fates, S.latdeg,
    S.forc_pbot, S.forc_t, S.forc_po2, S.forc_pco2, S.forc_pch4,
    S.watsat, S.h2osoi_vol, S.h2osoi_liq, S.h2osoi_ice, S.h2osfc, S.bsw, S.cellorg, S.smp_l, S.t_soisno,
    S.t_grnd, S.t_h2osfc, S.frac_h2osfc, S.snow_depth, S.snl, S.qflx_surf,
    S.rootfr_p, S.rootfr_col, S.crootfr, S.rootr_p, S.elai, S.qflx_tran_veg, S.annsum_npp, S.rr,
    S.somhr, S.lithr, S.hr_vr, S.o_scalar, S.fphr, S.pot_f_nit_vr, S.lake_icefrac, S.lakedepth,
    S.z, S.dz, S.zi, S.nlevsoi, 5, 1, S.nlevsoi, 5, 0.01, 130.0, S.ng, 0, 1800.0,
    true, false, false, S.agnpp, S.bgnpp, 365.0 * 86400.0;
    finundated_stream = S.stream, tws = S.tws)

function check_config(name, FT, label; finmtd)
    @printf("\n  --- config: %s ---\n", label)
    H = build(Float64; finmtd=finmtd); B = build(Float64; finmtd=finmtd)
    ad(x) = CLM.Adapt.adapt(MF32(), x)
    # Device snapshot: whole CH4Data + every array arg. The stream/tws stay HOST
    # (ch4! copies their coeff vectors to the backend internally).
    Sd = (; B.params, B.ch4vc, B.stream, B.tws, B.nc, B.np, B.ng, B.nlevsoi,
        mask_soil=ad(B.mask_soil), mask_soilp=ad(B.mask_soilp), mask_lake=ad(B.mask_lake), mask_nolake=ad(B.mask_nolake),
        col_gridcell=ad(B.col_gridcell), col_wtgcell=ad(B.col_wtgcell), patch_column=ad(B.patch_column),
        patch_itype=ad(B.patch_itype), patch_wtcol=ad(B.patch_wtcol), is_fates=ad(B.is_fates), latdeg=ad(B.latdeg),
        forc_pbot=ad(B.forc_pbot), forc_t=ad(B.forc_t), forc_po2=ad(B.forc_po2), forc_pco2=ad(B.forc_pco2), forc_pch4=ad(B.forc_pch4),
        watsat=ad(B.watsat), h2osoi_vol=ad(B.h2osoi_vol), h2osoi_liq=ad(B.h2osoi_liq), h2osoi_ice=ad(B.h2osoi_ice), h2osfc=ad(B.h2osfc),
        bsw=ad(B.bsw), cellorg=ad(B.cellorg), smp_l=ad(B.smp_l), t_soisno=ad(B.t_soisno), t_grnd=ad(B.t_grnd), t_h2osfc=ad(B.t_h2osfc),
        frac_h2osfc=ad(B.frac_h2osfc), snow_depth=ad(B.snow_depth), snl=ad(B.snl), qflx_surf=ad(B.qflx_surf),
        rootfr_p=ad(B.rootfr_p), rootfr_col=ad(B.rootfr_col), crootfr=ad(B.crootfr), rootr_p=ad(B.rootr_p), elai=ad(B.elai),
        qflx_tran_veg=ad(B.qflx_tran_veg), annsum_npp=ad(B.annsum_npp), rr=ad(B.rr),
        somhr=ad(B.somhr), lithr=ad(B.lithr), hr_vr=ad(B.hr_vr), o_scalar=ad(B.o_scalar), fphr=ad(B.fphr), pot_f_nit_vr=ad(B.pot_f_nit_vr),
        lake_icefrac=ad(B.lake_icefrac), lakedepth=ad(B.lakedepth), z=ad(B.z), dz=ad(B.dz), zi=ad(B.zi),
        agnpp=ad(B.agnpp), bgnpp=ad(B.bgnpp))
    ch4_d = CLM.Adapt.adapt(MF32(), B.ch4)

    run_ch4!(H, H.ch4); run_ch4!(Sd, ch4_d)

    checks = [
        ("ch4_surf_flux_tot", H.ch4.ch4_surf_flux_tot_col, ch4_d.ch4_surf_flux_tot_col),
        ("totcolch4_col",     H.ch4.totcolch4_col,         ch4_d.totcolch4_col),
        ("conc_ch4_sat",      H.ch4.conc_ch4_sat_col,      ch4_d.conc_ch4_sat_col),
        ("conc_ch4_unsat",    H.ch4.conc_ch4_unsat_col,    ch4_d.conc_ch4_unsat_col),
        ("c_atm_grc",         H.ch4.c_atm_grc,             ch4_d.c_atm_grc),
        ("finundated",        H.ch4.finundated_col,        ch4_d.finundated_col),
        ("grnd_ch4_cond",     H.ch4.grnd_ch4_cond_col,     ch4_d.grnd_ch4_cond_col),
        ("ch4co2f_grc",       H.ch4.ch4co2f_grc,           ch4_d.ch4co2f_grc),
        ("ch4prodg_grc",      H.ch4.ch4prodg_grc,          ch4_d.ch4prodg_grc),
        ("totcolch4_grc",     H.ch4.totcolch4_grc,         ch4_d.totcolch4_grc),
    ]
    nfail = 0
    for (nm, a, b) in checks
        if !allfinite(a)
            @printf("  [WARN] %-18s CPU non-finite — skipping\n", nm); continue
        end
        d = reldiff(a, b); ok = d < 1f-3
        @printf("  [%s] %-18s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end
    # For the inversion config, prove the kernel actually did regression work rather
    # than sitting at the cold-start finundated=0.1 fixed point (a no-op would be a
    # trivially-identical, meaningless "pass"). finundated = fws_slope*tws + fws_intercept
    # then lagged toward it; with the coeffs above the host value must be well off 0.1.
    if finmtd == CLM.FINUNDATION_MTD_TWS_INVERSION
        fh = Array(H.ch4.finundated_col)
        @printf("  [info] host finundated (TWS inversion) = %s\n", fh)
        if all(f -> abs(f - 0.1) < 1e-6, fh)
            println("  [FAIL] inversion kernel was a NO-OP (finundated pinned at cold-start 0.1)")
            nfail += 1
        end
    end
    return nfail
end

function main(backend)
    println("=" ^ 70); println("END-TO-END Metal parity for ch4! (methane orchestrator)"); println("=" ^ 70)
    backend === nothing && (println("  No GPU backend."); return 0)
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n", name, FT)

    nfail = 0
    nfail += check_config(name, FT, "h2osfc (default)";       finmtd=CLM.FINUNDATION_MTD_H2OSFC)
    nfail += check_config(name, FT, "TWS_inversion (#238)";   finmtd=CLM.FINUNDATION_MTD_TWS_INVERSION)
    println()
    println(nfail == 0 ? "  WHOLE ch4! MATCHES CPU ON $name ($FT) — both finundation methods" :
                         "  DIVERGENCE — investigate ($nfail failed).")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
