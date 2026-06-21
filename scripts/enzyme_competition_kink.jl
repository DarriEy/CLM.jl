# Validate the soil_bgc_competition! fpi smoothing AT THE KINK (supply == demand). The hard fpi
# = min(1, sminn/(demand·dt)) has a derivative KINK at supply==demand: slope 1/(demand·dt) from
# the N-limited side, 0 from the supply side. Enzyme reverse on the HARD form (SMOOTH_MODE=:auto)
# returns a one-sided SUBGRADIENT that disagrees with the central FD; the SMOOTH form (:always,
# smooth_min) returns the FD-consistent averaged derivative. This is the cross-domain-residual
# kink, now removed for reverse-AD while the forward stays byte-identical under :auto.
#   julia +1.10 --project=/tmp/clm_jl10_<id> scripts/enzyme_competition_kink.jl
using CLM, Enzyme, Printf
const C = CLM
# ONE mode per process (SMOOTH=auto|always) — set BEFORE any Enzyme compile. Enzyme caches the
# adjoint by function type, not by the runtime SMOOTH_MODE global, so compiling phase! once under
# :auto would poison it for :always within the same process. Run this script twice.
C.SMOOTH_MODE[] = (get(ENV, "SMOOTH", "auto") == "always") ? :always : :auto

# Minimal competition state placed EXACTLY at the kink, then perturb the mineral-N supply sminn.
const NC=1; const NLD=1; const NCT=2
const STATE = C.SoilBGCCompetitionState(); const PARAMS = C.SoilBGCCompetitionParams()
const DZ = fill(1.0, NLD); const PMNF = zeros(NC,NLD,NCT); const PCNG = zeros(NC,NLD,NCT); const CRP = ones(Int,NCT)
const AUX = (; mask=BitVector([true]), bounds=1:NC, nld=NLD, nct=NCT)

# choose inputs so sum_ndemand·dt == sminn (the kink): one column/level, plant demand + immob.
const DT = 1800.0
function fresh(sminn)
    st = C.SoilBiogeochemStateData()
    st.fpg_col=zeros(NC); st.fpi_col=zeros(NC); st.fpi_vr_col=zeros(NC,NLD)
    st.nfixation_prof_col=fill(0.5,NC,NLD); st.plant_ndemand_col=fill(1.0e-6,NC)
    ns = C.SoilBiogeochemNitrogenStateData()
    ns.sminn_vr_col=fill(sminn,NC,NLD); ns.smin_nh4_vr_col=fill(sminn/2,NC,NLD); ns.smin_no3_vr_col=fill(sminn/2,NC,NLD)
    nf = C.SoilBiogeochemNitrogenFluxData()
    for f in (:potential_immob_vr_col,:actual_immob_vr_col,:sminn_to_plant_vr_col,:supplement_to_sminn_vr_col,
              :sminn_to_denit_excess_vr_col,:actual_immob_no3_vr_col,:actual_immob_nh4_vr_col,:smin_no3_to_plant_vr_col,
              :smin_nh4_to_plant_vr_col,:pot_f_nit_vr_col,:pot_f_denit_vr_col,:f_nit_vr_col,:f_denit_vr_col,
              :n2_n2o_ratio_denit_vr_col,:f_n2o_denit_vr_col,:f_n2o_nit_vr_col,:sminn_to_plant_fun_vr_col,
              :sminn_to_plant_fun_no3_vr_col,:sminn_to_plant_fun_nh4_vr_col)
        setfield!(nf, f, zeros(NC,NLD))
    end
    nf.potential_immob_vr_col=fill(1.0e-6,NC,NLD)         # immob demand
    for f in (:sminn_to_plant_col,:actual_immob_col,:potential_immob_col); setfield!(nf,f,zeros(NC)); end
    cf = C.SoilBiogeochemCarbonFluxData(); cf.c_overflow_vr=zeros(NC,NLD,NCT)
    return (; st=st, nf=nf, cf=cf, ns=ns)
end
phase!(b,a)=(C.soil_bgc_competition!(b.st,b.nf,b.cf,b.ns,STATE,PARAMS; mask_bgc_soilc=a.mask,bounds=a.bounds,
    nlevdecomp=a.nld,ndecomp_cascade_transitions=a.nct,dzsoi_decomp=DZ,pmnf_decomp_cascade=PMNF,
    p_decomp_cn_gain=PCNG,cascade_receiver_pool=CRP,use_nitrif_denitrif=false,use_fun=false); nothing)
fpi(sminn)=(b=fresh(sminn); phase!(b,AUX); b.st.fpi_vr_col[1,1])

# find the kink sminn (where fpi crosses ~1): plant_demand·prof + immob ≈ sum_ndemand; pick sminn = that·dt
let lo=1e-9, hi=1e-2
    for _ in 1:50; mid=(lo+hi)/2; if fpi(mid) < 0.999; lo = mid; else; hi = mid; end; end
    global SMINN_KINK = (lo+hi)/2
end
@printf("kink at sminn≈%.4e   fpi(0.98·kink)=%.4f  fpi(kink)=%.4f  fpi(1.02·kink)=%.4f\n",
    SMINN_KINK, fpi(0.98*SMINN_KINK), fpi(SMINN_KINK), fpi(1.02*SMINN_KINK))

rev(sminn)=(b=fresh(sminn); db=C.compositional_reverse!(Any[(phase!,(AUX,))], fresh(sminn),
    (db,b)->(db.st.fpi_vr_col[1,1]=1.0)); db.ns.sminn_vr_col[1,1])
fd(sminn;h)= (fpi(sminn+h)-fpi(sminn-h))/(2h)

println("="^72)
let mode = C.SMOOTH_MODE[], h = 5e-4*SMINN_KINK   # << transition width (~kink/k) so FD resolves it
    r = rev(SMINN_KINK); f = fd(SMINN_KINK; h=h); rl = abs(r-f)/max(abs(f),1e-30)
    @printf("[%-7s] reverse d(fpi)/d(sminn)=% .5e   central-FD=% .5e   rel=%.2e  %s\n",
        mode, r, f, rl,
        mode===:auto  ? "(HARD: reverse = one-sided SUBGRADIENT, ~2× the averaged FD)" :
                        (rl < 0.05 ? "PASS ✓ (SMOOTH: reverse matches central FD at the kink)" : "FAIL ✗"))
end
println("="^72)
