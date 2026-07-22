# Regression guard for the MIMICS live-path WIRING (decomp_method==2 activation).
# Unlike test_decomp_mimics.jl (qualitative kernel unit tests), this asserts:
#   (a) decomp_rates_mimics! is BIT-EXACT vs an independent inline scalar oracle of
#       the Fortran decomp_rates_mimics kinetics, with the real clm50_params values
#       (densdep=1.2, vint=6.6, ...) via mimics_default_read_params! — non-vacuous.
#   (b) the activation config derivation: CLMDriverConfig(decomp_method=2) sizes the
#       MIMICS 8-pool/15-transition layout, CENTURY (default) is untouched, and the
#       CN-vegetation facade propagates decomp_method through _sync_driver_config!.
# Fast + portable: no NetCDF / Fortran / Bow-fixture deps.
@testset "MIMICS wiring (decomp_method==2 activation)" begin
    _SPHR = 3600.0; _SPDAY = 86400.0; _TFRZ = 273.15
    _G2MG = 1.0e3; _CM32M3 = 1.0e-6; _P2F = 1.0e-2

    # -- inline independent scalar oracle (multi-level, non-FATES, non-anoxia) --
    function _oracle(p; t_soisno, soilpsi, cpools, dz, cellclay, lignin, annpp,
                     minpsi, maxpsi, tau_cwd, dpy, icop, ioli)
        nc, nlev = size(t_soisno)
        dk = fill(NaN, nc, nlev, 8); pf = fill(NaN, nc, nlev, 15)
        rf = fill(NaN, nc, nlev, 15); ws = fill(NaN, nc, nlev); os = fill(NaN, nc, nlev)
        cnc = fill(NaN, nc); cno = fill(NaN, nc)
        vmod=p.mimics_vmod; vint=p.mimics_vint; vsl=p.mimics_vslope
        kmod=p.mimics_kmod; kint=p.mimics_kint; ksl=p.mimics_kslope; mge=p.mimics_mge
        f1,f2,f3,f4 = p.mimics_fmet[1:4]
        fcr1,fcr2=p.mimics_fchem_r[1:2]; fck1,fck2=p.mimics_fchem_k[1:2]
        fpr1,fpr2=p.mimics_fphys_r[1:2]; fpk1,fpk2=p.mimics_fphys_k[1:2]
        ds1,ds2=p.mimics_desorp[1:2]; ps1,ps2=p.mimics_p_scalar[1:2]
        tr1,tr2=p.mimics_tau_r[1:2]; tk1,tk2=p.mimics_tau_k[1:2]
        tmn,tmx,tmf=p.mimics_tau_mod_min,p.mimics_tau_mod_max,p.mimics_tau_mod_factor
        kor,kok=p.mimics_ko_r,p.mimics_ko_k; dd=p.mimics_densdep
        dq=p.mimics_desorpQ10; tref=p.mimics_t_soi_ref
        cnm=p.mimics_cn_mod_num; cnr=p.mimics_cn_r; cnk=p.mimics_cn_k
        kfrag = 1.0/(_SPDAY*dpy*tau_cwd)
        r_l1m1=1-mge[1]; r_l2m1=1-mge[2]; r_s1m1=1-mge[3]
        r_l1m2=1-mge[4]; r_l2m2=1-mge[5]; r_s1m2=1-mge[6]
        for c in 1:nc
            fmet=f1*(f2-f3*min(f4,lignin[c]))
            tmod=min(tmx,max(tmn,sqrt(tmf*max(0.0,annpp[c]))))
            tm1=tr1*exp(tr2*fmet)*tmod/_SPHR; tm2=tk1*exp(tk2*fmet)*tmod/_SPHR
            cnc[c]=cnr*sqrt(cnm/fmet); cno[c]=cnk*sqrt(cnm/fmet)
            fcm1=min(1.0,max(0.0,fcr1*exp(fcr2*fmet)))
            fcm2=min(1.0,max(0.0,fck1*exp(fck2*fmet)))
            for j in 1:nlev
                psi=min(soilpsi[c,j],maxpsi)
                w = psi>minpsi ? log(minpsi/psi)/log(minpsi/maxpsi) : 0.0
                ws[c,j]=w; os[c,j]=1.0; wdo=w
                clay=_P2F*min(100.0,cellclay[c,j])
                dsp=ds1*exp(ds2*clay); fpm1=min(1.0,fpr1*exp(fpr2*clay))
                fpm2=min(1.0,fpk1*exp(fpk2*clay)); psc=1.0/(ps1*exp(ps2*sqrt(clay)))
                td=t_soisno[c,j]-_TFRZ
                vl1m1=exp(vsl[1]*td+vint[1])*vmod[1]/_SPHR
                vl2m1=exp(vsl[2]*td+vint[2])*vmod[2]/_SPHR
                vs1m1=exp(vsl[3]*td+vint[3])*vmod[3]/_SPHR
                vl1m2=exp(vsl[4]*td+vint[4])*vmod[4]/_SPHR
                vl2m2=exp(vsl[5]*td+vint[5])*vmod[5]/_SPHR
                vs1m2=exp(vsl[6]*td+vint[6])*vmod[6]/_SPHR
                kl1m1=exp(ksl[1]*td+kint[1])*kmod[1]
                kl2m1=exp(ksl[2]*td+kint[2])*kmod[2]
                ks1m1=exp(ksl[3]*td+kint[3])*kmod[3]*psc
                kl1m2=exp(ksl[4]*td+kint[4])*kmod[4]
                kl2m2=exp(ksl[5]*td+kint[5])*kmod[5]
                ks1m2=exp(ksl[6]*td+kint[6])*kmod[6]*psc
                dsn=(dsp/_SPHR)*dq*exp((td-tref)/10.0)
                m1=(cpools[c,j,icop]/dz[c,j])*_G2MG*_CM32M3
                m2=(cpools[c,j,ioli]/dz[c,j])*_G2MG*_CM32M3
                t1=vl1m1*m1/(kl1m1+m1); t2=vl1m2*m2/(kl1m2+m2)
                dk[c,j,1]=(t1+t2)*wdo
                if (t1+t2)!=0; pf[c,j,1]=t1/(t1+t2); pf[c,j,2]=t2/(t1+t2) else; pf[c,j,1]=0.0; pf[c,j,2]=0.0 end
                t1=vl2m1*m1/(kl2m1+m1); t2=vl2m2*m2/(kl2m2+m2)
                dk[c,j,2]=(t1+t2)*wdo
                if (t1+t2)!=0; pf[c,j,3]=t1/(t1+t2); pf[c,j,4]=t2/(t1+t2) else; pf[c,j,3]=0.0; pf[c,j,4]=0.0 end
                t1=vs1m1*m1/(ks1m1+m1); t2=vs1m2*m2/(ks1m2+m2)
                dk[c,j,3]=(t1+t2)*wdo
                if (t1+t2)!=0; pf[c,j,5]=t1/(t1+t2); pf[c,j,6]=t2/(t1+t2) else; pf[c,j,5]=0.0; pf[c,j,6]=0.0 end
                dk[c,j,5]=dsn
                t1=vl2m1*m1/(kor*kl2m1+m1); t2=vl2m2*m2/(kok*kl2m2+m2)
                dk[c,j,4]=(t1+t2)*wdo
                dk[c,j,6]=tm1*m1^(dd-1.0)*wdo
                pf[c,j,9]=min(1.0,max(0.0,1.0-fpm1-fcm1)); pf[c,j,10]=fcm1
                dk[c,j,7]=tm2*m2^(dd-1.0)*wdo
                pf[c,j,12]=min(1.0,max(0.0,1.0-fpm2-fcm2)); pf[c,j,13]=fcm2
                dk[c,j,8]=kfrag*wdo
                pf[c,j,7]=1.0; pf[c,j,8]=1.0; pf[c,j,11]=fpm1; pf[c,j,14]=fpm2; pf[c,j,15]=1.0
                rf[c,j,1]=r_l1m1; rf[c,j,2]=r_l1m2; rf[c,j,3]=r_l2m1
                rf[c,j,4]=r_l2m2; rf[c,j,5]=r_s1m1; rf[c,j,6]=r_s1m2
                for k in (7,8,9,10,11,12,13,14,15); rf[c,j,k]=0.0 end
            end
        end
        return (; dk, pf, rf, ws, os, cnc, cno)
    end
    _mrd(a,b) = begin m=0.0; for i in eachindex(a,b)
        isnan(b[i]) && continue; @test isfinite(b[i])
        m=max(m, abs(a[i]-b[i])/(1.0+max(abs(a[i]),abs(b[i])))) end; m end

    @testset "(a) decomp_rates_mimics! bit-exact vs scalar oracle" begin
        nc=3; nlev=4
        p = CLM.DecompMIMICSParams(); CLM.mimics_default_read_params!(p)
        @test p.mimics_densdep == 1.2   # confirms clm50_params values are in effect
        st = CLM.DecompMIMICSState(); casc = CLM.DecompCascadeConData()
        casc.initial_stock_soildepth = 0.0
        cnp = CLM.CNSharedParamsData()
        CLM.cn_shared_params_read!(cnp; q10_mr=1.5, minpsi_hr=-2.0, maxpsi_hr=-0.002,
            rf_cwdl2=0.0, tau_cwd=3.3333333, cwd_flig=0.24, decomp_depth_efolding=10.0,
            froz_q10=1.5, mino2lim=0.2, organic_max=130.0)
        cellclay = [12.0+11.0*(c-1)+3.0*(j-1) for c in 1:nc, j in 1:nlev]
        CLM.init_decompcascade_mimics!(st, casc, p, cnp; cellclay=cellclay,
            bounds=1:nc, nlevdecomp=nlev, ndecomp_pools_max=8,
            ndecomp_cascade_transitions_max=15, spinup_state=0, use_fates=false)
        cf = CLM.SoilBiogeochemCarbonFluxData()
        CLM.soil_bgc_carbon_flux_init!(cf, nc, nlev, 8, 15; nlevdecomp=nlev)
        t_soisno = [276.0+6.0*sin(0.6c+0.4j) for c in 1:nc, j in 1:nlev]
        soilpsi  = [-0.03-0.4*((c+j)%5)/5 for c in 1:nc, j in 1:nlev]
        dz       = [0.03+0.04*(j-1) for c in 1:nc, j in 1:nlev]
        cpools   = [8.0+3.0*c+2.0*j+1.5*k for c in 1:nc, j in 1:nlev, k in 1:8]
        lignin   = [18.0+6.0*c for c in 1:nc]
        annpp    = [80.0+50.0*c for c in 1:nc]
        CLM.decomp_rates_mimics!(cf, st, p, cnp, casc;
            mask_bgc_soilc=trues(nc), bounds=1:nc, nlevdecomp=nlev,
            t_soisno=t_soisno, soilpsi=soilpsi, decomp_cpools_vr=cpools, col_dz=dz,
            ligninNratioAvg=lignin, annsum_npp_col=annpp, days_per_year=365.0, dt=1800.0,
            spinup_state=0, use_lch4=false, anoxia=false, use_fates=false)
        ref = _oracle(p; t_soisno, soilpsi, cpools, dz, cellclay, lignin, annpp,
            minpsi=cnp.minpsi, maxpsi=cnp.maxpsi, tau_cwd=cnp.tau_cwd, dpy=365.0,
            icop=st.i_cop_mic, ioli=st.i_oli_mic)
        @test _mrd(Array(cf.decomp_k_col), ref.dk) == 0.0
        @test _mrd(Array(cf.pathfrac_decomp_cascade_col), ref.pf) == 0.0
        @test _mrd(Array(cf.rf_decomp_cascade_col), ref.rf) == 0.0
        @test _mrd(Array(cf.w_scalar_col), ref.ws) == 0.0
        @test _mrd(Array(cf.o_scalar_col), ref.os) == 0.0
        @test _mrd(Array(cf.cn_col)[:, st.i_cop_mic], ref.cnc) == 0.0
        @test _mrd(Array(cf.cn_col)[:, st.i_oli_mic], ref.cno) == 0.0
        # non-vacuous: decomp_k actually populated and positive
        dkv = Array(cf.decomp_k_col)
        @test all(isfinite, dkv)
        @test count(>(0.0), dkv) == length(dkv)
    end

    @testset "(b) activation config derivation + facade sync" begin
        m = CLM.CLMDriverConfig(use_cn=true, decomp_method=2)
        @test m.decomp_method == 2
        @test m.ndecomp_pools == 8
        @test m.ndecomp_cascade_transitions == 15
        @test m.i_litr_min == 1
        @test m.i_litr_max == 2
        @test m.i_cwd == 8
        # CENTURY default (decomp_method=1) untouched
        c = CLM.CLMDriverConfig(use_cn=true)
        @test c.decomp_method == 1
        @test c.ndecomp_pools == 7
        @test c.ndecomp_cascade_transitions == 10
        @test c.i_litr_max == 3
        @test c.i_cwd == 7
        # facade propagates decomp_method into the CN driver config
        veg = CLM.CNVegetationData(); veg.config.use_cn = true; veg.config.decomp_method = 2
        CLM._sync_driver_config!(veg)
        @test veg.driver_config.decomp_method == 2
    end

    @testset "(c) litr_lig_c_to_n_col (fmet input) vs scalar oracle" begin
        # Ports SoilBiogeochemCarbonFluxType.F90:961-1016 (mimics_decomp block).
        # Fixture: 4 columns / 5 patches covering all branches —
        #   col1: full case (litterfall + cwdn>0 + cwd_ctransfer>0)
        #   col2: cwdn==0 branch (cwd term forced to 0)
        #   col3: FATES column (skipped; litr_lig_c_to_n_col left untouched)
        #   col4: floor case (zero litterfall + zero cwd_ctransfer ⇒ max(1e-3,·) denom)
        ncol = 4; ntrans = 15; i_cwdl2 = 15
        cwd_flig = 0.24

        # pftcon vectors, indexed by ivt = itype + 1 (1-based)
        lf_flig = [0.10, 0.20, 0.25]
        lflitcn = [50.0, 60.0, 70.0]
        fr_flig = [0.15, 0.30, 0.35]
        frootcn = [40.0, 45.0, 42.0]

        # patches: 1,2 → col1 ; 3 → col2 ; 4 → col3(FATES) ; 5 → col4
        itype   = [1, 2, 1, 0, 2]          # ivt = itype+1
        wtcol   = [0.6, 0.4, 1.0, 1.0, 1.0]
        leafc_to_litter_patch  = [0.03, 0.05, 0.02, 0.99, 0.0]
        frootc_to_litter_patch = [0.02, 0.01, 0.04, 0.88, 0.0]
        mask_vegp = Bool[true, true, true, false, true]   # patch4 inactive

        cwdc_col = [3.0, 2.0, 9.9, 1.0]
        cwdn_col = [0.1, 0.0, 9.9, 0.05]                  # col2 → 0 branch
        cwd_ctransfer_by_col = [0.008, 0.004, 9.9, 0.0]   # col4 → floor branch

        col = CLM.ColumnData()
        col.patchi   = [1, 3, 4, 5]
        col.patchf   = [2, 3, 4, 5]
        col.is_fates = [false, false, true, false]

        patch = CLM.PatchData()
        patch.itype = itype
        patch.wtcol = wtcol

        cf = CLM.SoilBiogeochemCarbonFluxData()
        cf.litr_lig_c_to_n_col = fill(-999.0, ncol)       # sentinel (col3 must survive)
        cf.decomp_cascade_ctransfer_col = zeros(ncol, ntrans)
        for c in 1:ncol
            cf.decomp_cascade_ctransfer_col[c, i_cwdl2] = cwd_ctransfer_by_col[c]
        end

        CLM.soil_bgc_carbon_flux_lignin_n_ratio!(cf, cwdc_col, cwdn_col,
            mask_vegp, 1:ncol; col=col, patch=patch,
            lf_flig=lf_flig, lflitcn=lflitcn, fr_flig=fr_flig, frootcn=frootcn,
            leafc_to_litter_patch=leafc_to_litter_patch,
            frootc_to_litter_patch=frootc_to_litter_patch,
            cwd_flig=cwd_flig, i_cwdl2=i_cwdl2)

        # --- independent plain-scalar oracle ---
        function oracle(patches)
            lnl = 0.0; lnf = 0.0; lc = 0.0; fc = 0.0
            for p in patches
                mask_vegp[p] || continue
                ivt = itype[p] + 1; wt = wtcol[p]
                lnl += lf_flig[ivt] * lflitcn[ivt] * leafc_to_litter_patch[p] * wt
                lnf += fr_flig[ivt] * frootcn[ivt] * frootc_to_litter_patch[p] * wt
                lc  += leafc_to_litter_patch[p] * wt
                fc  += frootc_to_litter_patch[p] * wt
            end
            (lnl, lnf, lc, fc)
        end
        function ratio(patches, c)
            lnl, lnf, lc, fc = oracle(patches)
            ct = cwd_ctransfer_by_col[c]
            lncwd = cwdn_col[c] > 0.0 ? cwd_flig * (cwdc_col[c] / cwdn_col[c]) * ct : 0.0
            (lnl + lnf + lncwd) / max(1.0e-3, lc + fc + ct)
        end

        @test cf.litr_lig_c_to_n_col[1] ≈ ratio(1:2, 1) rtol=1e-12
        @test cf.litr_lig_c_to_n_col[2] ≈ ratio(3:3, 2) rtol=1e-12
        @test cf.litr_lig_c_to_n_col[3] == -999.0            # FATES column untouched
        @test cf.litr_lig_c_to_n_col[4] ≈ ratio(5:5, 4) rtol=1e-12

        # Non-vacuous: col1 is a real nonzero ratio; col4 exercises the 1e-3 floor
        # (numerator 0, denom would be 0 without the floor).
        @test cf.litr_lig_c_to_n_col[1] > 0.0
        @test cf.litr_lig_c_to_n_col[2] > 0.0
        @test cf.litr_lig_c_to_n_col[4] == 0.0
        # cwdn==0 branch: col2 ratio must NOT include a cwd term
        lnl2, lnf2, lc2, fc2 = oracle(3:3)
        @test cf.litr_lig_c_to_n_col[2] ≈ (lnl2 + lnf2) /
            max(1.0e-3, lc2 + fc2 + cwd_ctransfer_by_col[2]) rtol=1e-12
    end
end
