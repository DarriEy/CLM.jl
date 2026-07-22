#!/usr/bin/env julia
# =============================================================================
# validate_voc_megan_activation.jl
#
# ACTIVATED-PATH end-to-end harness for the VOC/MEGAN activation plumbing.
#
# The kernel physics is already oracle-validated bit-exact in
# scripts/validate_voc_megan_fortran_parity.jl. THIS harness proves the thing that
# harness could NOT reach: that the LIVE INIT PLUMBING actually feeds the kernel.
# Specifically it drives voc_emission! through descriptors + factors + efisop built
# ENTIRELY by the new live initializers — NOT hand-built:
#   * megan_config_init         (shr_megan_readnl + megan_factors_init analogue):
#                                reads megan_factors_file -> table + per-class mf,
#                                parses megan_specifier -> compound/mech descriptors
#   * read_efisop_from_surfdata! (VOCEmission::iniTimeConst analogue): EF1_* -> efisop
#   * CLMDriverConfig(use_voc=true, megan_specifier=…, megan_factors_file=…)
#
# and asserts:
#   1. NONZERO isoprene/monoterpene flux matching an INDEPENDENT scalar oracle
#      (transcribed from VOCEmissionMod.F90) to ~1e-10, using the SAME per-compound
#      EFs + per-class coefficients the init path loaded from the file.
#   2. efisop_grc is finite/populated (not the NaN it stays without iniTimeConst).
#   3. the DEFAULT (VOC-off) CLMDriverConfig is inert / byte-identical (=== empty).
#
# Uses the REAL megan_factors_file + a REAL surfdata carrying EF1_* when present on
# the local inputdata tree; otherwise writes CTSM-layout synthetic fixtures and
# SAYS SO (never fabricates a value silently, never disables TLS to fetch data).
#
# Run:  julia +1.12 --project=. scripts/validate_voc_megan_activation.jl
# Exit: 0 = all activation-path checks pass, 1 = a divergence / wiring regression.
# =============================================================================
using CLM
using NCDatasets
using Printf

const TOL  = 1e-10
const TFRZ = 273.15

# ---- independent scalar oracle (transcribed from VOCEmissionMod.F90) ---------
ref_gamma_L(f240, elai) = (f240 > 0.0 && f240 < 1.0e30) ? 0.30*elai : 0.24*elai
ref_gamma_SM(btran) = btran >= 1.0 ? 1.0 : 1.0/(1.0 + 3.2552*exp(-7.4463*(btran-0.2)))
function ref_gamma_P(par_sun,par24_sun,par240_sun,par_sha,par24_sha,par240_sha,
                     fsun,fsun240,solad240,solai240,LDF)
    ca1,ca2,ca3 = 0.004,0.0005,0.0468
    if fsun240 > 0.0 && fsun240 < 1.0 && solad240 > 0.0 && solai240 > 0.0
        a=ca1-ca2*log(par240_sun); cp=ca3*exp(ca2*(par24_sun-200.0))*par240_sun^0.6
        gp=fsun*(cp*a*par_sun*(1.0+a*a*par_sun*par_sun)^(-0.5))
        a=ca1-ca2*log(par240_sha); cp=ca3*exp(ca2*(par_sha-50.0))*par240_sha^0.6
        gp+= (1.0-fsun)*(cp*a*par_sha*(1.0+a*a*par_sha*par_sha)^(-0.5))
    else
        a,cp=0.001,1.21
        gp=fsun*(cp*a*par_sun*(1.0+a*a*par_sun*par_sun)^(-0.5))
        gp+=(1.0-fsun)*(cp*a*par_sha*(1.0+a*a*par_sha*par_sha)^(-0.5))
    end
    (1.0-LDF)+LDF*gp
end
function ref_gamma_T(t240,t24,tv,ct1,ct2,betaT,LDF,Ceo,ivt)
    if t240 > 0.0 && t240 < 1.0e30
        topt=313.0+0.6*(t240-297.0)
        if ivt==11;      Eopt=7.9*exp(0.217*(t24-TFRZ-24.0))
        elseif ivt==12;  Eopt=exp(0.12*(t240-TFRZ-15.0))
        else;            Eopt=Ceo*exp(0.05*(t24-297.0))*exp(0.05*(t240-297.0)); end
    else
        topt,Eopt=317.0,2.26
    end
    x=((1.0/topt)-(1.0/tv))/0.00831
    if ivt==12
        bet=min(95.0+9.49*exp(0.53*(TFRZ+15.0-t240)),300.0)
        gLDF=Eopt*exp(bet*((1.0/(TFRZ+30.0)-1.0/tv)/0.00831))
    else
        gLDF=Eopt*(ct2*exp(ct1*x)/(ct2-ct1*(1.0-exp(ct2*x))))
    end
    gLIF=exp(betaT*(tv-303.15))
    (1.0-LDF)*gLIF+LDF*gLDF
end
function ref_gamma_A(ivt,e240,e,nc,Anew,Agro,Amat,Aold)
    (ivt==3 || ivt>=6) || return 1.0
    (e240>0.0 && e240<1.0e30) || return 1.0
    ep=2.0*e240-e
    if ep==e; fn,fg,fm,fo=0.0,0.0,1.0,0.0
    elseif ep>e; fn,fg=0.0,0.0; fm=1.0-(ep-e)/ep; fo=(ep-e)/ep
    else; fn=1.0-(ep/e); fg=0.0; fm=ep/e; fo=0.0; end
    fn*Anew[nc]+fg*Agro[nc]+fm*Amat[nc]+fo*Aold[nc]
end
function ref_gamma_C(cisun,cisha,pbot,fsun,co2)
    gca=1.344-(1.344*(0.7*co2)^1.4614)/(585.0^1.4614+(0.7*co2)^1.4614)
    if co2<400.0; Is,h,Cs=1.072,1.70,1218.0
    elseif co2<600.0; t=(co2-400.0)/200.0; Is=t*1.036+(1-t)*1.072; h=t*2.0125+(1-t)*1.70; Cs=t*1150.0+(1-t)*1218.0
    elseif co2<800.0; t=(co2-600.0)/200.0; Is=t*1.046+(1-t)*1.036; h=t*1.5380+(1-t)*2.0125; Cs=t*2025.0+(1-t)*1150.0
    else; Is,h,Cs=1.014,2.861,1525.0; end
    if !isnan(cisun) && !isnan(cisha) && pbot>0.0 && fsun>0.0
        ci=(fsun*cisun+(1.0-fsun)*cisha)/pbot*1.0e6; gci=Is-(Is*ci^h)/(Cs^h+ci^h)
    else; gci=1.0; end
    gci*gca
end
ref_map_EF(ivt,g,ef)= ivt in (1,2) ? ef[2,g] : ivt==3 ? ef[3,g] : (4<=ivt<=8) ? ef[1,g] :
    (9<=ivt<=11) ? ef[4,g] : (12<=ivt<=14) ? ef[5,g] : ivt>=15 ? ef[6,g] : 0.0

const NFAIL = Ref(0)
function check(name, got, exp; tol=TOL)
    d = abs(got-exp); ok = d <= tol*max(1.0,abs(exp)); ok || (NFAIL[]+=1)
    @printf("  %-46s got=% .8e ref=% .8e |Δ|=% .1e  %s\n", name, got, exp, d, ok ? "ok" : "**FAIL**")
    ok
end

# ---- locate REAL fixtures, else synthesize (and say so) ----------------------
const _REAL_MFF = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/installs/cesm-inputdata/atm/cam/chem/trop_mozart/emis/megan21_emis_factors_78pft_c20161108.nc"
const _REAL_SD_CANDIDATES = [
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Peatland_MerBleue_Canada/settings/CLM/parameters/surfdata_clm.nc",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_lake_run/surfdata_lake100.nc",
]

function _synth_mff(path)
    names=["isoprene","myrcene","limonene"]; npft=20; nclass=3
    NCDataset(path,"c") do ds
        defDim(ds,"Comp_Num",3); defDim(ds,"Class_Num",nclass); defDim(ds,"PFT_Num",npft); defDim(ds,"char",40)
        cn=fill('\0',40,3); for (j,nm) in enumerate(names), (i,c) in enumerate(collect(nm)); cn[i,j]=c; end
        defVar(ds,"Comp_Name",Char,("char","Comp_Num"))[:,:]=cn
        defVar(ds,"Comp_MW",Float64,("Comp_Num",))[:]=[68.12,136.23,136.23]
        defVar(ds,"Class_Num",Int16,("Comp_Num",))[:]=Int16[1,2,3]
        defVar(ds,"Comp_EF",Float64,("Comp_Num","PFT_Num"))[:,:]=[Float64(2+i) for i in 1:3, k in 1:npft]
        defVar(ds,"Class_EF",Float64,("Class_Num","PFT_Num"))[:,:]=[Float64(10c) for c in 1:nclass, k in 1:npft]
        defVar(ds,"LDF",Float64,("Class_Num",))[:]=[0.999,0.4,0.2]
        for (nm,val) in (("Agro",0.6),("Amat",1.125),("Anew",0.05),("Aold",1.0),("betaT",0.10),("ct1",95.0),("ct2",230.0),("Ceo",2.0))
            defVar(ds,nm,Float64,("Class_Num",))[:]=fill(val,nclass)
        end
    end
    path
end
function _synth_sd(path)
    NCDataset(path,"c") do ds
        defDim(ds,"lsmlon",1); defDim(ds,"lsmlat",1)
        for (k,vn) in enumerate(("EF1_BTR","EF1_FET","EF1_FDT","EF1_SHR","EF1_GRS","EF1_CRP"))
            defVar(ds,vn,Float64,("lsmlon","lsmlat"))[:,:]=fill((10000.0,2000.0,2000.0,4000.0,800.0,1.0)[k],1,1)
        end
    end
    path
end

println("="^78)
println("VOC/MEGAN ACTIVATED-PATH harness (init plumbing -> kernel, indep oracle)")
println("="^78)

tmp = mktempdir()
if isfile(_REAL_MFF)
    mff = _REAL_MFF; println("  megan_factors_file : REAL  $(mff)")
else
    mff = _synth_mff(joinpath(tmp,"megan_factors.nc"))
    println("  megan_factors_file : SYNTHETIC (real file not on local inputdata) $(mff)")
end
_real_sd = findfirst(isfile, _REAL_SD_CANDIDATES)
if _real_sd !== nothing && (NCDataset(ds->haskey(ds,"EF1_BTR"), _REAL_SD_CANDIDATES[_real_sd], "r"))
    sd = _REAL_SD_CANDIDATES[_real_sd]; println("  surfdata (EF1_*)   : REAL  $(sd)")
else
    sd = _synth_sd(joinpath(tmp,"surfdata_ef.nc"))
    println("  surfdata (EF1_*)   : SYNTHETIC (no local surfdata with EF1_*) $(sd)")
end

# =============================================================================
# [1] Build EVERYTHING through the live init path (no hand-built descriptors).
# =============================================================================
println("\n[1] Live init path: CLMDriverConfig(use_voc=true, specifier, factors_file)")
spec = ["ISOP = isoprene", "BIGENE = isoprene + 0.5*myrcene"]
cfg  = CLM.CLMDriverConfig(use_voc=true, megan_specifier=spec, megan_factors_file=mff)
mc   = cfg.megan
@printf("  use_voc=%s meg_compounds=%d mech_comps=%d mapped_emisfctrs=%s (CTSM default .false.)\n",
        cfg.use_voc, length(mc.meg_compounds), length(mc.mech_comps), mc.mapped_emisfctrs)
(cfg.use_voc && length(mc.meg_compounds)==2 && length(mc.mech_comps)==2) ||
    (NFAIL[]+=1; println("  **FAIL** init path did not build the descriptors"))

# efisop via the iniTimeConst analogue
np,nc,ng = 5,4,2
voc = CLM.VOCEmisData(); CLM.vocemis_init!(voc, np, ng, length(mc.meg_compounds), length(mc.mech_comps))
@assert all(isnan, voc.efisop_grc)   # NaN until the live read
efok = CLM.read_efisop_from_surfdata!(voc, sd, ones(Int, ng))   # map both gridcells to cell 1
@printf("  efisop read ok=%s finite=%s efisop[:,1]=%s\n", efok, all(isfinite,voc.efisop_grc), voc.efisop_grc[:,1])
(all(isfinite, voc.efisop_grc)) || (NFAIL[]+=1; println("  **FAIL** efisop not populated"))

# =============================================================================
# [2] Drive the kernel through the init-built config; match independent oracle.
# =============================================================================
println("\n[2] voc_emission! through init-built descriptors vs independent oracle")
patch = CLM.PatchData(); CLM.patch_init!(patch, np)
patch.itype    .= [0, 1, 6, 12, 11]     # noveg, needleleaf-evr, broadleaf-decid-trop, arctic-C3, boreal-shrub
patch.gridcell .= [1, 1, 1, 2, 2]
patch.column   .= [1, 1, 2, 3, 4]
forc_solad = fill(220.0, nc, 2); forc_solad[:,1] .= [210.0,215.0,225.0,205.0]
forc_solai = fill(95.0, ng, 2);  forc_solai[:,1] .= [90.0,100.0]
forc_pbot  = fill(101325.0, nc); forc_pco2 = fill(45.0, ng)
fsd24=[180.,190.,200.,175.,185.]; fsd240=[170.,178.,188.,168.,176.]
fsi24=[85.,88.,92.,83.,87.];      fsi240=[80.,84.,89.,79.,83.]
fsun=[0.5,0.55,0.6,0.45,0.5]; fsun24=[0.5,0.52,0.58,0.47,0.5]; fsun240=[0.5,0.53,0.57,0.48,0.5]
elai=[2.0,2.2,1.8,1.2,1.5];   elai240=[2.0,2.1,1.6,1.3,1.4]
civ=0.7*45.0; cisun=fill(civ,np,1); cisha=fill(0.9*civ,np,1)
t_veg=[300.,301.,299.,293.,297.]; t_veg24=[300.,300.5,298.5,292.,296.]; t_veg240=[297.,298.,296.,290.,294.]
btran=[0.8,0.7,0.9,0.6,0.85]

CLM.voc_emission!(voc, mc.meg_compounds, mc.mech_comps, mc.megan_factors,
    patch, 1:np, trues(np),
    forc_solad, forc_solai, forc_pbot, forc_pco2,
    fsd24, fsd240, fsi24, fsi240, fsun, fsun24, fsun240, elai, elai240,
    cisun, cisha, t_veg, t_veg24, t_veg240, btran;
    use_mapped_emisfctrs = mc.mapped_emisfctrs)

# independent scalar recomputation of the whole pipeline using the SAME descriptors
mf = mc.megan_factors
megfac = 1.0/3600.0/1.0e6
exp_tot = zeros(np)
for p in 1:np
    ivt = patch.itype[p]; ivt==0 && continue
    g=patch.gridcell[p]; c=patch.column[p]
    par_sun    = (forc_solad[c,1]+fsun[p]*forc_solai[g,1])*4.6
    par24_sun  = (fsd24[p]+fsun24[p]*fsi24[p])*4.6
    par240_sun = (fsd240[p]+fsun240[p]*fsi240[p])*4.6
    par_sha    = ((1.0-fsun[p])*forc_solai[g,1])*4.6
    par24_sha  = ((1.0-fsun24[p])*fsi24[p])*4.6
    par240_sha = ((1.0-fsun240[p])*fsi240[p])*4.6
    gl=ref_gamma_L(fsun240[p],elai[p]); gsm=ref_gamma_SM(btran[p])
    vmeg=zeros(length(mc.meg_compounds))
    for (im,cmp) in enumerate(mc.meg_compounds)
        cn=cmp.class_number; is_iso = cmp.name=="isoprene"
        eps = (is_iso && mc.mapped_emisfctrs) ? ref_map_EF(ivt,g,voc.efisop_grc) : cmp.emis_factors[ivt]
        gp=ref_gamma_P(par_sun,par24_sun,par240_sun,par_sha,par24_sha,par240_sha,fsun[p],fsun240[p],fsd240[p],fsi240[p],mf.LDF[cn])
        gt=ref_gamma_T(t_veg240[p],t_veg24[p],t_veg[p],mf.ct1[cn],mf.ct2[cn],mf.betaT[cn],mf.LDF[cn],mf.Ceo[cn],ivt)
        ga=ref_gamma_A(ivt,elai240[p],elai[p],cn,mf.Anew,mf.Agro,mf.Amat,mf.Aold)
        gc= is_iso ? ref_gamma_C(cisun[p,1],cisha[p,1],forc_pbot[c],fsun[p],1.0e6*forc_pco2[g]/forc_pbot[c]) : 1.0
        gamma=gl*gsm*ga*gp*gt*gc
        (gamma>=0.0 && gamma<100.0) && (vmeg[im]=cmp.coeff*eps*gamma*megfac/cmp.molec_weight)
    end
    for imech in 1:length(mc.mech_comps)
        acc=0.0
        for ii in 1:mc.mech_comps[imech].n_megan_comps; acc+=vmeg[mc.mech_comps[imech].megan_indices[ii]]; end
        exp_tot[p]+=acc
    end
end
for p in 1:np
    check("vocflx_tot[$p]", voc.vocflx_tot_patch[p], exp_tot[p])
end
if any(>(0.0), exp_tot) && any(p->voc.vocflx_tot_patch[p]>0.0, 1:np)
    println("  NON-VACUITY ok: at least one vegetated patch carries real nonzero flux.")
else
    NFAIL[]+=1; println("  **FAIL** vacuous — no nonzero flux produced")
end
(voc.vocflx_tot_patch[1]==0.0) || (NFAIL[]+=1; println("  **FAIL** noveg patch must stay zero"))

# =============================================================================
# [3] DEFAULT config remains VOC-off / inert (byte-identical).
# =============================================================================
println("\n[3] Default CLMDriverConfig is inert (byte-identical)")
d = CLM.CLMDriverConfig()
gate = d.use_voc && !isempty(d.megan.mech_comps) && !isempty(d.megan.meg_compounds)
@printf("  use_voc=%s meg=%d mech=%d mapped=%s gate=%s (expect false / 0 / 0 / false / false)\n",
        d.use_voc, length(d.megan.meg_compounds), length(d.megan.mech_comps), d.megan.mapped_emisfctrs, gate)
(d.use_voc===false && isempty(d.megan.meg_compounds) && isempty(d.megan.mech_comps) &&
 d.megan.mapped_emisfctrs===false && gate===false) ||
    (NFAIL[]+=1; println("  **FAIL** default config not inert"))

println("\n" * "="^78)
if NFAIL[]==0
    println("RESULT: VOC/MEGAN ACTIVATION path WIRED + validated end-to-end.")
    println("  Init plumbing (megan_config_init + read_efisop_from_surfdata!) feeds the")
    println("  kernel: nonzero isoprene/monoterpene flux matches the independent oracle")
    println("  to <1e-10, efisop populated finite, default run stays VOC-off / inert.")
    exit(0)
else
    println("RESULT: $(NFAIL[]) divergence(s) — see **FAIL** lines above."); exit(1)
end
