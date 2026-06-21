# Fire burned-area reverse phase. cnfire_area_li2014! is a THRESHOLD-discrete phase but is fully
# written with the smooth primitives (smooth_max/smooth_min throughout → continuous farea_burned,
# no FSM flag), so under SMOOTH_MODE=:always it reverse-differentiates: perturb a continuous fire
# driver (forc_rh, through the fire-weather thresholds), seed farea_burned, FD-validate. Set
# SMOOTH=always|auto via ENV (one mode per process — Enzyme caches the adjoint by function type).
#   SMOOTH=always julia +1.10 --project=/tmp/clm_jl10_<id> scripts/enzyme_fire.jl
using CLM, Enzyme, Printf
const C = CLM
C.SMOOTH_MODE[] = (get(ENV,"SMOOTH","always")=="always") ? :always : :auto
const NP=4; const NC=2; const NG=1; const NLG=3; const NLD=1; const NDP=4; const NLITR=3; const NPFT=20

# --- read-only Const aux (params, PFT cons, geometry, forcing) ---
const PFTCON = C.PftConFireBase(woody=vcat(fill(1.0,8),fill(0.0,12)), cc_leaf=fill(0.4,NPFT),
    cc_lstem=fill(0.2,NPFT), cc_dstem=fill(0.1,NPFT), cc_other=fill(0.3,NPFT), fm_leaf=fill(0.6,NPFT),
    fm_lstem=fill(0.5,NPFT), fm_other=fill(0.4,NPFT), fm_root=fill(0.3,NPFT), fm_lroot=fill(0.5,NPFT),
    fm_droot=fill(0.2,NPFT), lf_f=fill(1.0/NLITR,NPFT,NLITR), fr_f=fill(1.0/NLITR,NPFT,NLITR),
    smpso=fill(-66000.0,NPFT), smpsc=fill(-275000.0,NPFT))
const PFTLI = C.PftConFireLi2014(fsr_pft=fill(0.2,NPFT), fd_pft=fill(1.0,NPFT))
const FCONST = C.CNFireConstData(); const FPARAMS = C.CNFireParams(prh30=0.05, ignition_efficiency=0.02)
const PATCH = (p=C.PatchData(); p.itype=[2,10,5,15]; p.column=[1,1,2,2]; p.wtcol=[0.5,0.5,0.6,0.4]; p)
const COL = (c=C.ColumnData(); c.gridcell=[1,1]; c)
const GRC = (g=C.GridcellData(); g.latdeg=[45.0]; g.lat=[45.0*pi/180]; g)
const SS = (s=C.SoilStateData(); s.watsat_col=fill(0.45,NC,NLG); s.rootfr_patch=fill(1.0/NLG,NP,NLG);
    s.sucsat_col=fill(200.0,NC,NLG); s.bsw_col=fill(5.0,NC,NLG); s)
const H2OVOL = fill(0.1, NC, NLG)
const CASC = (c=C.DecompCascadeConData(); c.is_litter=BitVector([true,true,true,false]); c.is_cwd=BitVector([false,false,false,true]); c)
const TOTLITC = [200.0,150.0]; const DCPOOLS = fill(100.0,NC,NLD,NDP); const TSOI17 = fill(280.0,NC)
const WIND = [3.0]; const FORCT = fill(290.0,NC); const RAIN = fill(0.0,NC); const SNOW = fill(0.0,NC)
const PREC60 = fill(2.0e-5,NP); const PREC10 = fill(3.0e-5,NP); const FSAT = fill(0.1,NC); const WF=fill(0.3,NC); const WF2=fill(0.25,NC)
const MSC = trues(NC); const MSP = trues(NP); const MEV = trues(NP); const MNEV = falses(NP)
const AUX = (;)

# --- per-call bundle: structs fire mutates + the perturbed input forc_rh ---
function fresh()
    fli = C.CNFireLi2014Data(forc_hdm=[50.0], forc_lnfm=[0.05], gdp_lf_col=[10.0,10.0], peatf_lf_col=[0.0,0.0], abm_lf_col=[6,6])
    fdata = C.CNFireBaseData(btran2_patch=zeros(NP))
    cst = C.CNVegStateData()
    for f in (:dwt_smoothed_patch,); setfield!(cst,f,zeros(NP)); end
    for f in (:cropf_col,:baf_crop_col,:baf_peatf_col,:fbac_col,:fbac1_col,:farea_burned_col,:nfire_col,
              :fsr_col,:fd_col,:lgdp_col,:lgdp1_col,:lpop_col,:lfwt_col,:trotr1_col,:trotr2_col,:dtrotr_col,:lfc_col,:wtlf_col)
        setfield!(cst,f,zeros(NC))
    end
    cst.burndate_patch=fill(10000,NP)
    cs = C.CNVegCarbonStateData()
    cs.totvegc_col=[500.0,400.0]; cs.rootc_col=zeros(NC); cs.leafc_col=zeros(NC); cs.fuelc_col=zeros(NC); cs.fuelc_crop_col=zeros(NC)
    cs.leafc_patch=[10.0,5.0,8.0,3.0]; cs.leafc_storage_patch=[1.0,0.5,0.8,0.3]; cs.leafc_xfer_patch=[0.5,0.25,0.4,0.15]
    cs.frootc_patch=[4.0,2.0,3.0,1.0]; cs.frootc_storage_patch=[0.5,0.2,0.3,0.1]; cs.frootc_xfer_patch=[0.2,0.1,0.15,0.05]
    cs.deadcrootc_patch=[25.0,0.0,20.0,0.0]; cs.deadcrootc_storage_patch=[1.5,0.0,1.2,0.0]; cs.deadcrootc_xfer_patch=[0.8,0.0,0.6,0.0]
    cs.livecrootc_patch=[10.0,0.0,8.0,0.0]; cs.livecrootc_storage_patch=[1.0,0.0,0.8,0.0]; cs.livecrootc_xfer_patch=[0.5,0.0,0.4,0.0]
    return (; fli=fli, fdata=fdata, cst=cst, cs=cs, rh=[50.0])
end
phase!(b, a) = (C.cnfire_area_li2014!(b.fli, PFTLI, b.fdata, FCONST, FPARAMS, PFTCON,
    MSC, MSP, MEV, MNEV, 1:NC, 1:NP, PATCH, COL, GRC, SS, H2OVOL, b.cst, b.cs, CASC,
    TOTLITC, DCPOOLS, TSOI17; forc_rh_grc=b.rh, forc_wind_grc=WIND, forc_t_col=FORCT,
    forc_rain_col=RAIN, forc_snow_col=SNOW, prec60_patch=PREC60, prec10_patch=PREC10,
    fsat_col=FSAT, wf_col=WF, wf2_col=WF2, dt=1800.0, dayspyr=365.0, kmo=6, kda=15, mcsec=3600,
    nstep=10, nlevgrnd=NLG, nlevdecomp=NLD, ndecomp_pools=NDP); nothing)

farea(rh)=(b=fresh(); b.rh[1]=rh; phase!(b,AUX); copy(b.cst.farea_burned_col))
let b=fresh(); phase!(b,AUX)
    @printf("primal: farea_burned=%s finite=%s  (SMOOTH_MODE=%s)\n",
        string(round.(b.cst.farea_burned_col;sigdigits=4)), all(isfinite,b.cst.farea_burned_col), C.SMOOTH_MODE[])
end
L(b)=sum(abs2, b.cst.farea_burned_col)
richardson(f)=(cfd(h)=(f(h)-f(-h))/(2h); (4*cfd(5e-3)-cfd(1e-2))/3)
g_fd = richardson(δ -> (b=fresh(); b.rh[1]+=δ; phase!(b,AUX); L(b)))
db = C.compositional_reverse!(Any[(phase!,(AUX,))], fresh(),
    (db,b)->(db.cst.farea_burned_col .= 2 .* b.cst.farea_burned_col))
g_rev = db.rh[1]; rl = abs(g_rev-g_fd)/max(abs(g_fd),1e-30)
println("="^72)
@printf("FIRE burned-area REVERSE [%s]  dL/d(forc_rh)  rev=% .5e  FD=% .5e  rel=%.2e  %s\n",
    C.SMOOTH_MODE[], g_rev, g_fd, rl, rl<1e-3 ? "PASS ✓ (fire reverses under Enzyme)" : "FAIL ✗")
println("="^72)
