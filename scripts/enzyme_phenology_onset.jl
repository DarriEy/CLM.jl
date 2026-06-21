# Phenology ONSET-THRESHOLD one-time-event smoothing — validates the reverse gradient
# d(leaf-out burst)/d(crit_onset_gdd) that the hard event LACKED. The onset burst (phenology.jl
# season_decid kernel) is now leafc_storage_to_xfer = smooth_heaviside(og − crit)·fstor2tran·
# storage/dt: under :auto smooth_heaviside is the exact step (ow=1 ⇒ byte-identical, suite 190/190),
# under :always it ramps as og crosses crit → a FINITE d/d(crit) where the instant event had a
# delta-function gradient. This exercises the EXACT expression added to the kernel, one mode per
# process (Enzyme caches the adjoint by function type).
#   SMOOTH=always|auto julia +1.10 --project=/tmp/clm_jl10_<id> scripts/enzyme_phenology_onset.jl
using CLM, Enzyme, Printf
const C = CLM
C.SMOOTH_MODE[] = (get(ENV,"SMOOTH","always")=="always") ? :always : :auto
const FSTOR2TRAN = 0.5; const DT = 1800.0
const OG = 100.1          # onset_gdd just past the threshold (og−crit=0.1): inside the k=50 smoothing
                          # band (saturates at |og−crit|>36/50≈0.72), but off the exact step so the HARD
                          # :auto branch is cleanly flat (FD=0) while :always keeps a finite gradient.
const STORAGE = 5.0       # leafc_storage

# the exact onset-burst expression from the kernel: leafc_storage_to_xfer = ow·fstor2tran·storage/dt
burst(crit) = (b=[0.0]; b[1] = C.smooth_heaviside(OG - crit) * FSTOR2TRAN * STORAGE / DT; b)
L(b) = b[1]^2
phase!(b, a) = (b.o[1] = C.smooth_heaviside(OG - b.crit[1]) * FSTOR2TRAN * STORAGE / DT; nothing)
fresh(crit) = (; o=[0.0], crit=[crit])

const CRIT = 100.0        # crit_onset_gdd; og−crit = 2 (near threshold under :always)
let b=fresh(CRIT); phase!(b,nothing)
    @printf("primal: leafc_storage_to_xfer=%.6e  (SMOOTH_MODE=%s; og−crit=%.1f)\n", b.o[1], C.SMOOTH_MODE[], OG-CRIT)
end
g_fd = (b1=fresh(CRIT); b1.crit[1]+=1e-4; phase!(b1,nothing); a=b1.o[1]^2;
        b2=fresh(CRIT); b2.crit[1]-=1e-4; phase!(b2,nothing); (a-b2.o[1]^2)/2e-4)
db = C.compositional_reverse!(Any[(phase!,(nothing,))], fresh(CRIT), (db,b)->(db.o[1]=2*b.o[1]))
g_rev = db.crit[1]; rl = abs(g_rev-g_fd)/max(abs(g_fd),1e-30)
println("="^74)
@printf("ONSET-THRESHOLD REVERSE [%s]  dL/d(crit_onset_gdd)  rev=% .5e  FD=% .5e  rel=%.2e  %s\n",
    C.SMOOTH_MODE[], g_rev, g_fd, rl,
    C.SMOOTH_MODE[]===:auto ? "(HARD: gradient = 0, the delta-function event is invisible)" :
                              (rl<1e-3 ? "PASS ✓ (SMOOTH: finite gradient, reverse matches FD)" : "FAIL ✗"))
println("="^74)
