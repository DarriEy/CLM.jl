# Gap mortality reverse phase: cn_gap_mortality! is a CONTINUOUS/linear BGC phase
# (m_leafc_to_litter = leafc · r_mort_rate; its only branches are on the STATIC `woody` PFT flag
# and `>0` true-zeros), so it reverse-differentiates at machine precision with no smoothing needed
# — extending the reverse chain to the mortality flux. Perturb leafc, seed m_leafc_to_litter.
#   julia +1.10 --project=/tmp/clm_jl10_<id> scripts/enzyme_gapmort.jl
using CLM, Enzyme, Printf
const C = CLM
const NP=3; const NC=2; const NLD=1; const NLITR=3

const PARAMS = C.GapMortalityParams(k_mort=0.3, r_mort=fill(0.02, 21))
const PFTCON = C.PftConGapMort(woody=vcat(fill(1.0,9),fill(0.0,12)), leafcn=fill(25.0,21),
    livewdcn=fill(50.0,21), lf_f=fill(1.0/NLITR,21,NLITR), fr_f=fill(1.0/NLITR,21,NLITR))
const DGVS = C.DgvsGapMortData(greffic_patch=fill(0.5,NP), heatstress_patch=fill(0.0,NP), nind_patch=fill(100.0,NP))
const PATCH = (p=C.PatchData(); p.itype=[2,10,5]; p.column=[1,1,2]; p.wtcol=[0.5,0.3,0.2]; p)
const CST = (c=C.CanopyStateData(); C.canopystate_init!(c, NP); c.laisun_patch.=2.0; c.laisha_patch.=1.0; c)
const AUX = (; mask=BitVector([true,true,true]), bounds=1:NP)

function fresh()
    cs = C.CNVegCarbonStateData(); C.cnveg_carbon_state_init!(cs, NP, NC, 1)
    cf = C.CNVegCarbonFluxData(); C.cnveg_carbon_flux_init!(cf, NP, NC, 1)
    ns = C.CNVegNitrogenStateData(); C.cnveg_nitrogen_state_init!(ns, NP, NC, 1)
    nf = C.CNVegNitrogenFluxData(); C.cnveg_nitrogen_flux_init!(nf, NP, NC, 1)
    for s in (cs, cf, ns, nf), f in fieldnames(typeof(s))   # init NaN-fills → zero (else 0·NaN poisons reverse)
        v = getfield(s, f); v isa AbstractArray{<:AbstractFloat} && fill!(v, 0.0)
    end
    cs.leafc_patch .= [10.0, 5.0, 8.0]
    cs.frootc_patch .= [4.0, 2.0, 3.0]; cs.livestemc_patch .= [20.0,0.0,15.0]; cs.deadstemc_patch .= [50.0,0.0,40.0]
    cs.livecrootc_patch .= [10.0,0.0,8.0]; cs.deadcrootc_patch .= [25.0,0.0,20.0]
    ns.leafn_patch .= [0.4,0.2,0.32]
    return (; cs=cs, cf=cf, ns=ns, nf=nf)
end

phase!(b, a) = (C.cn_gap_mortality!(a.mask, a.bounds, PARAMS, PFTCON, DGVS, PATCH, CST,
    b.cs, b.cf, b.ns, b.nf; dt=1800.0, days_per_year=365.0, use_cndv=false, use_matrixcn=false,
    spinup_state=0, npcropmin=17, spinup_factor_deadwood=1.0); nothing)

m = 0.02 / (365.0 * C.SECSPDAY)
let b=fresh(); phase!(b,AUX)
    @printf("primal: m_leafc_to_litter=%s (expect leafc·m = %s)\n",
        string(round.(b.cf.m_leafc_to_litter_patch; sigdigits=4)), string(round.([10,5,8].*m; sigdigits=4)))
end
L(b) = sum(abs2, b.cf.m_leafc_to_litter_patch)
richardson(f) = (cfd(h)=(f(h)-f(-h))/(2h); (4*cfd(5e-3)-cfd(1e-2))/3)
g_fd = richardson(δ -> (b=fresh(); b.cs.leafc_patch[1]+=δ; phase!(b,AUX); L(b)))
db = C.compositional_reverse!(Any[(phase!,(AUX,))], fresh(),
    (db,b)->(db.cf.m_leafc_to_litter_patch .= 2 .* b.cf.m_leafc_to_litter_patch))
g_rev = db.cs.leafc_patch[1]
rl = abs(g_rev-g_fd)/max(abs(g_fd),1e-30)
println("="^66)
@printf("gap mortality REVERSE  dL/d(leafc)  rev=% .6e  FD=% .6e  rel=%.2e  %s\n",
    g_rev, g_fd, rl, rl<1e-5 ? "PASS ✓ (machine precision)" : "FAIL ✗")
println("="^66)
