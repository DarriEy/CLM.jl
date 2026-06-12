# Compositional reverse-AD over REAL coupled canopy sub-phases — NUMERICAL proof.
#
# enzyme_localize_canopy.jl proved: the canopy_fluxes_core! monolith hits the Enzyme
# sret_ty assertion, but EVERY constituent sub-phase differentiates cleanly in
# isolation. The unblock is COMPOSITIONAL reverse-AD — differentiate each sub-phase
# with a SEPARATE Enzyme.autodiff call (each compiles its own smaller thunk, dodging
# the monolith's sret bug) and chain the adjoints, exactly as
# compositional_enzyme_proto.jl validates on toy phases.
#
# Here we prove it NUMERICALLY on two genuinely-coupled real CLM sub-phases, on the
# finite single-patch synthetic state from gpu_validate_canopy_e2e.jl (the driver's
# real warmup leaves NaN radiation on this tiny domain — irrelevant to the sret
# COMPILE bug, but a numerical check needs finite inputs):
#
#     phase1 = friction_velocity!   writes frictionvel.ustar_patch, um_patch
#     phase2 = cf_resist_update!     reads ustar/um, writes frictionvel.ram1_patch
#
# The coupling crosses the phase boundary through the inst state fields ustar/um.
# We perturb the upstream Monin-Obukhov length obu, and validate the chained
# input-adjoint  dL/d(obu),  L = sum(abs2, ram1),  against central finite
# differences of the SAME primal chain.
#
#   julia +1.10 --project=/tmp/clm_jl10_env scripts/enzyme_compositional_canopy.jl
#
using CLM, Enzyme, Printf

# --- config (mirrors gpu_validate_canopy_e2e.jl setconfig!) ---
CLM.soil_resistance_read_nl!(soil_resis_method = CLM.SOIL_RESIS_LEEPIELKE_1992)
CLM.canopy_fluxes_read_nml!(use_undercanopy_stability = false,
    use_biomass_heat_storage = false, itmax_canopy_fluxes = 40)
CLM.canopy_fluxes_read_params!()

# --- finite single-patch state (trimmed from gpu_validate_canopy_e2e.jl build()) ---
const FT = Float64
CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
const NP = 1; const NC = 1; const NG = 1; const NL = 1

function build_state()
    canopystate = CLM.CanopyStateData{FT}(); CLM.canopystate_init!(canopystate, NP)
    canopystate.elai_patch[1] = 2.0;  canopystate.esai_patch[1] = 0.5
    canopystate.laisun_patch[1] = 1.2; canopystate.laisha_patch[1] = 0.8
    canopystate.displa_patch[1] = 5.0; canopystate.htop_patch[1] = 10.0
    canopystate.frac_veg_nosno_patch[1] = 1; canopystate.dleaf_patch[1] = 0.04
    canopystate.stem_biomass_patch[1] = 0.0; canopystate.leaf_biomass_patch[1] = 0.0

    energyflux = CLM.EnergyFluxData{FT}(); CLM.energyflux_init!(energyflux, NP, NC, NL, NG)
    energyflux.btran_patch[1] = 0.5; energyflux.bsun_patch[1] = 0.5; energyflux.bsha_patch[1] = 0.5
    energyflux.htvp_col[1] = CLM.HVAP

    frictionvel = CLM.FrictionVelocityData{FT}(); CLM.frictionvel_init!(frictionvel, NP, NC)
    frictionvel.zetamaxstable = 0.5; frictionvel.zsno = 0.00085; frictionvel.zlnd = 0.000775
    frictionvel.z0mv_patch[1] = 0.5; frictionvel.z0hv_patch[1] = 0.5; frictionvel.z0qv_patch[1] = 0.5
    frictionvel.z0mg_col[1] = 0.01; frictionvel.z0hg_col[1] = 0.01; frictionvel.z0qg_col[1] = 0.01
    frictionvel.forc_hgt_u_patch[1] = 30.0; frictionvel.forc_hgt_t_patch[1] = 30.0
    frictionvel.forc_hgt_q_patch[1] = 30.0
    frictionvel.ustar_patch[1] = 0.5; frictionvel.um_patch[1] = 5.0; frictionvel.uaf_patch[1] = 3.0
    frictionvel.taf_patch[1] = 290.0; frictionvel.qaf_patch[1] = 0.008
    frictionvel.obu_patch[1] = -100.0; frictionvel.zeta_patch[1] = -0.1
    frictionvel.vpd_patch[1] = 1.0; frictionvel.rb1_patch[1] = 50.0
    frictionvel.ram1_patch[1] = 50.0; frictionvel.num_iter_patch[1] = 0.0

    temperature = CLM.TemperatureData{FT}(); CLM.temperature_init!(temperature, NP, NC, NL, NG)
    temperature.t_veg_patch[1] = 290.0; temperature.t_stem_patch[1] = 289.0
    temperature.t_skin_patch[1] = 290.0; temperature.thm_patch[1] = 290.0
    temperature.t_grnd_col[1] = 288.0; temperature.t_h2osfc_col[1] = 288.0
    temperature.thv_col[1] = 291.0; temperature.emv_patch[1] = 0.97; temperature.emg_col[1] = 0.96
    for j in 1:size(temperature.t_soisno_col, 2); temperature.t_soisno_col[1, j] = 288.0; end

    waterdiagbulk = CLM.WaterDiagnosticBulkData{FT}(); CLM.waterdiagnosticbulk_init!(waterdiagbulk, NC, NP, NL, NG)
    waterdiagbulk.fwet_patch[1] = 0.1; waterdiagbulk.fdry_patch[1] = 0.8
    waterdiagbulk.frac_sno_eff_col[1] = 0.0; waterdiagbulk.frac_h2osfc_col[1] = 0.0
    waterdiagbulk.snow_depth_col[1] = 0.0; waterdiagbulk.qg_col[1] = 0.005
    waterdiagbulk.qg_snow_col[1] = 0.005; waterdiagbulk.qg_soil_col[1] = 0.005
    waterdiagbulk.qg_h2osfc_col[1] = 0.005; waterdiagbulk.dqgdT_col[1] = 0.0003
    waterdiagbulk.rh_af_patch[1] = 0.6

    return (; frictionvel, canopystate, temperature, waterdiagbulk)
end

# Const (non-differentiated) auxiliaries captured by the phase closures.
patch_data = CLM.PatchData{FT}(); CLM.patch_init!(patch_data, NP)
patch_data.column[1] = 1; patch_data.gridcell[1] = 1; patch_data.landunit[1] = 1
patch_data.itype[1] = 1; patch_data.active[1] = true
const PATCH = patch_data
const FILTERP = Int[1]
const FN = 1
const ACTIVE = Bool[true]
const FORC_PBOT_COL = FT[101325.0]
const DLEAF_PFT = fill(FT(0.04), CLM.MXPFT + 1)
const GRND_CH4 = fill(FT(0.0), NP)
const PARAMS_CF = CLM.canopy_fluxes_params
const CSOILC = PARAMS_CF.csoilc
mkc(v) = fill(FT(v), NP)

# --- The two real coupled sub-phases (mutate bits, return nothing) ---
function phase1!(b)   # friction_velocity!: obu -> ustar, um
    fv = b.frictionvel
    ur = mkc(3.2); temp1 = mkc(0.0); temp2 = mkc(0.0)
    temp12m = mkc(0.0); temp22m = mkc(0.0); fm = mkc(0.0)
    CLM.friction_velocity!(fv, FN, FILTERP,
        b.canopystate.displa_patch, fv.z0mv_patch, fv.z0hv_patch, fv.z0qv_patch,
        fv.obu_patch, 1, ur, fv.um_patch, fv.ustar_patch,
        temp1, temp2, temp12m, temp22m, fm; active = ACTIVE)
    return nothing
end

function phase2!(b)   # cf_resist_update!: ustar,um -> ram1
    z(a...) = zeros(FT, a...)
    # temp1/temp2 are the φ-profile factors friction_velocity! produces; they appear
    # in rah=1/(temp1·ustar) (a dead output here — L depends on ram1 only). Held at a
    # fixed nonzero value: zero would make rah=1/0=Inf with an Inf derivative that
    # reverse-mode multiplies by the dead 0 adjoint → 0·Inf = NaN. The objective is
    # independent of temp1/temp2, so any nonzero value yields the correct dL/d(ustar).
    temp1 = mkc(0.1); temp2 = mkc(0.1); tlbef = mkc(289.0); del2 = z(NP); del_arr = z(NP)
    rah = z(NP, 2); raw = z(NP, 2); uuc = z(NP); rb = mkc(40.0)
    svpts = mkc(1900.0); eah = mkc(1400.0); el = mkc(1900.0)
    CLM.cf_resist_update!(b.frictionvel, b.canopystate, b.temperature, b.waterdiagbulk,
        PATCH, FILTERP, FN, ACTIVE,
        temp1, temp2, tlbef, del2, del_arr, rah, raw, uuc, rb, svpts, eah, el,
        DLEAF_PFT, GRND_CH4, FORC_PBOT_COL,
        CSOILC, false, false, false, PARAMS_CF)
    return nothing
end

chain!(b) = (phase1!(b); phase2!(b); nothing)
Lval(b) = sum(abs2, b.frictionvel.ram1_patch)

bits0 = build_state()

# Sanity: the primal chain must be finite on this state.
let s = deepcopy(bits0); chain!(s)
    @printf("primal: ustar=%.6g  um=%.6g  ram1=%.6g  L=%.6g\n",
        s.frictionvel.ustar_patch[1], s.frictionvel.um_patch[1],
        s.frictionvel.ram1_patch[1], Lval(s))
end

# ---------------------------------------------------------------------------
# Reference: central finite differences of L w.r.t. initial obu[1].
# ---------------------------------------------------------------------------
function L_perturbed(δ)
    s = deepcopy(bits0)
    s.frictionvel.obu_patch[1] += δ
    chain!(s)
    return Lval(s)
end
h = 1e-3
g_fd = (L_perturbed(h) - L_perturbed(-h)) / (2h)
@printf("\nFD   dL/d(obu)        = % .8e\n", g_fd)

# ---------------------------------------------------------------------------
# Compositional reverse-AD: per-phase separate Enzyme calls, chained adjoint.
# ---------------------------------------------------------------------------
Enzyme.API.strictAliasing!(false)
revmode = Enzyme.set_runtime_activity(Enzyme.Reverse)

# --- Per-phase isolation: localize any NaN to a single phase's adjoint ---
println("\n--- per-phase reverse vs FD ---")
# phase1 alone: L1 = sum(abs2, ustar);  d L1 / d(obu)
let
    f1(δ) = (s = deepcopy(bits0); s.frictionvel.obu_patch[1] += δ; phase1!(s);
             sum(abs2, s.frictionvel.ustar_patch))
    fd1 = (f1(h) - f1(-h)) / (2h)
    sf = deepcopy(bits0); phase1!(sf); ustar1 = copy(sf.frictionvel.ustar_patch)
    d = Enzyme.make_zero(bits0); d.frictionvel.ustar_patch .= 2.0 .* ustar1
    b = deepcopy(bits0)
    Enzyme.autodiff(revmode, (x) -> (phase1!(x); nothing), Enzyme.Const, Enzyme.Duplicated(b, d))
    @printf("  phase1  dL1/d(obu):   FD=% .6e  rev=% .6e\n", fd1, d.frictionvel.obu_patch[1])
end
# phase2 alone: L2 = sum(abs2, ram1);  d L2 / d(ustar_in), input = post-phase1 state
let
    base = deepcopy(bits0); phase1!(base)        # the input state to phase2
    f2(δ) = (s = deepcopy(base); s.frictionvel.ustar_patch[1] += δ; phase2!(s);
             sum(abs2, s.frictionvel.ram1_patch))
    fd2 = (f2(h) - f2(-h)) / (2h)
    sf = deepcopy(base); phase2!(sf); ram1b = copy(sf.frictionvel.ram1_patch)
    d = Enzyme.make_zero(bits0); d.frictionvel.ram1_patch .= 2.0 .* ram1b
    b = deepcopy(base)
    Enzyme.autodiff(revmode, (x) -> (phase2!(x); nothing), Enzyme.Const, Enzyme.Duplicated(b, d))
    @printf("  phase2  dL2/d(ustar): FD=% .6e  rev=% .6e\n", fd2, d.frictionvel.ustar_patch[1])
end

# Forward sweep with input checkpoints.
s = deepcopy(bits0)
snap1 = deepcopy(s);  phase1!(s)
snap2 = deepcopy(s);  phase2!(s)
ram1_final = copy(s.frictionvel.ram1_patch)

# Seed the running adjoint at the chain output: dL/d(ram1) = 2 ram1.
dbits = Enzyme.make_zero(bits0)
dbits.frictionvel.ram1_patch .= 2.0 .* ram1_final

# Reverse phase 2 (dbits: output-adjoint -> input-adjoint of phase2).
let b2 = deepcopy(snap2)
    Enzyme.autodiff(revmode, (x) -> (phase2!(x); nothing), Enzyme.Const,
                    Enzyme.Duplicated(b2, dbits))
end
# Reverse phase 1 (dbits -> input-adjoint of phase1 == chain-start adjoint).
let b1 = deepcopy(snap1)
    Enzyme.autodiff(revmode, (x) -> (phase1!(x); nothing), Enzyme.Const,
                    Enzyme.Duplicated(b1, dbits))
end
g_comp = dbits.frictionvel.obu_patch[1]
@printf("comp dL/d(obu)        = % .8e\n", g_comp)

# ---------------------------------------------------------------------------
# Compare.
# ---------------------------------------------------------------------------
println("\n", "="^62)
relerr = abs(g_comp - g_fd) / max(abs(g_fd), 1e-10)
@printf("abs error  = %.3e\n", abs(g_comp - g_fd))
@printf("rel error  = %.3e\n", relerr)
if isfinite(g_comp) && isfinite(g_fd) && relerr < 1e-3
    println("\nCOMPOSITIONAL CANOPY REVERSE-AD VALIDATED ✓")
    println("two real coupled sub-phases, chained via separate Enzyme calls,")
    println("gradient matches finite differences.")
else
    println("\nMISMATCH ✗ — investigate")
end
println("="^62)
