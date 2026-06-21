# Phenology onset-growth flux reverse. REVISES the earlier "phenology blocked" conclusion: the
# FSM flag (onset_flag) is a FIXED per-step INPUT state, so the onset-growth flux
# (leafc_xfer_to_leafc = t1·leafc_xfer, t1=2/onset_counter) is PIECEWISE-LINEAR in the pools given
# the flag — it reverse-differentiates fine. Only the THRESHOLD-CROSSING gradient (d/d(crit_onset_gdd),
# which flips the flag — a one-time triggered event) is discontinuous and would need event-smoothing.
# Here: onset active (flag=1), perturb leafc_xfer, seed leafc_xfer_to_leafc → reverse = FD.
#   julia +1.10 --project=/tmp/clm_jl10_<id> scripts/enzyme_phenology.jl
using CLM, Enzyme, Printf
const C = CLM
const NP=1; const DT=1800.0
# Const aux: the FSM state (flag fixed = onset active) + indices.
const MASK=BitVector([true]); const ITYPE=[1]; const WOODY=zeros(80)
const ONSET_FLAG=[1.0]; const ONSET_COUNTER=[10.0*C.SECSPDAY]; const BGTR=[0.0]
z() = zeros(NP)
# bundle = all 12 output flux arrays + 12 input xfer arrays (only leaf/froot used for non-woody)
function fresh()
    o = (; lc=z(),fc=z(),ln=z(),fn=z(), lsc=z(),dsc=z(),lcc=z(),dcc=z(), lsn=z(),dsn=z(),lcn=z(),dcn=z())
    x = (; lc=[5.0],fc=[3.0],ln=[0.2],fn=[0.12], lsc=z(),dsc=z(),lcc=z(),dcc=z(), lsn=z(),dsn=z(),lcn=z(),dcn=z())
    return (; o=o, x=x)
end
phase!(b, a) = (C.phen_onset_growth!(
    b.o.lc,b.o.fc,b.o.ln,b.o.fn, b.o.lsc,b.o.dsc,b.o.lcc,b.o.dcc, b.o.lsn,b.o.dsn,b.o.lcn,b.o.dcn,
    MASK, ITYPE, WOODY, ONSET_FLAG, ONSET_COUNTER, BGTR,
    b.x.lc,b.x.fc,b.x.ln,b.x.fn, b.x.lsc,b.x.dsc,b.x.lcc,b.x.dcc, b.x.lsn,b.x.dsn,b.x.lcn,b.x.dcn, DT); nothing)

let b=fresh(); phase!(b,nothing); t1=2.0/ONSET_COUNTER[1]
    @printf("primal: leafc_xfer_to_leafc=%.6e (expect t1·leafc_xfer = %.6e, t1=%.4e)\n",
        b.o.lc[1], t1*5.0, t1)
end
L(b)=sum(abs2, b.o.lc)
richardson(f)=(cfd(h)=(f(h)-f(-h))/(2h); (4*cfd(5e-3)-cfd(1e-2))/3)
g_fd = richardson(δ -> (b=fresh(); b.x.lc[1]+=δ; phase!(b,nothing); L(b)))
db = C.compositional_reverse!(Any[(phase!,(nothing,))], fresh(), (db,b)->(db.o.lc .= 2 .* b.o.lc))
g_rev = db.x.lc[1]; rl = abs(g_rev-g_fd)/max(abs(g_fd),1e-30)
println("="^70)
@printf("PHENOLOGY onset-growth REVERSE  dL/d(leafc_xfer)  rev=% .6e  FD=% .6e  rel=%.2e  %s\n",
    g_rev, g_fd, rl, rl<1e-5 ? "PASS ✓ (phenology reverses for pool gradients; flag = fixed input)" : "FAIL ✗")
println("="^70)
