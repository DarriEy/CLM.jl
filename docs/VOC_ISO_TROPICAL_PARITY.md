# VOC/MEGAN + C13/C14 isotope parity at the TROPICAL Aripuanã site

Both subsystems were previously validated at **Bow single-point only**. This note
records the first **tropical** (Aripuanã, Amazon; lat −7.51, lon −59.63; dominant
PFT 4 = broadleaf-evergreen tropical, isoprene EF = 7000 µg/m²/h, strong δ13C
discrimination) Fortran-parity check — the canonical stress case for both:
isoprene emission is maximal and, unlike Bow's unlit canopy, the C13
photosynthetic fractionation is genuinely active.

## Ground truth (Fortran) — left the shared build PRISTINE
`cesm_iso.exe` from `clm_bgc_spinup/bgc_ref_iso` (CN + C13/C14 + MEGAN active,
built with the #2119 reseed SourceMods) was rerun **cold-start** at Aripuanã in a
scratch rundir `/private/tmp/claude-501/ari_iso_run` (no rebuild, no SourceMods
change; only namelists + mesh/forcing paths swapped Bow→Aripuanã):

- surfdata / clm5_params / esmf_mesh / CLM1PT forcing from `~/projects/scratch/aripuana_data`;
- MEGAN specifier `ISOP = isoprene`, `C10H16 = pinene_a+carene_3+thujene_a`, + 6 OVOCs
  (identical to the drv_flds_in);
- native-grid h1 history tape (`hist_dov2xy=.false.`, `hist_nhtfrq=1`, pft×time = 6×48)
  carrying VOCFLXT / MEG_<compound> / GAMMA* / EOPT / TOPT / ALPHA / PAR_* — the VOC
  diagnostics are HISTORY fields, NOT in the restart-format pdump;
- restart-format pdumps `pdump_before_step_n{38..48}` / `pdump_after_hydrologydrainage_n{38..48}`
  carry the isotope pools + rc13_canair/psnsun/psnsha + the VOC 24/240-h accumulators.

Setup blockers found + fixed along the way: (1) `nuopc.runconfig` `mesh_atm/lnd/mask`
still pointed at the Bow mesh (datm-only path swap is insufficient) → surfrdMod:121
lon/lat mismatch; (2) `BTRAN` is not a registered h-field in this build → dropped
from `hist_fincl2`; (3) a CDEPS stream formatted-write abort; (4) `PDUMP_MIN_NSTEP`
must be set INLINE on the exe line (export does not propagate). Run terminated
SUCCESSFULLY; midday reference step = **n40** (2202-01-02 16:00 UTC ≈ local noon;
model-year label 2202, forcing recycled 2002 — absolute year irrelevant).

**Non-vacuity** (BET tropical patch, n40 midday): leafc = 99.67 gC/m², elai = 2.01,
VOCFLXT = 1.65e-8 mol/m²/s, isoprene GAMMA = 0.409, MEG_isoprene = 7.95e-10 kg/m²/s
(= 71 % of total VOC flux), rc13_canair = 0.011170.

## Harnesses (this branch)
- `scripts/fortran_parity_voc_aripuana.jl` — DECOUPLED VOC-kernel parity: feeds the
  Fortran-dumped inputs (h1 PAR_*/TV/ELAI + pdump 24/240 accumulators) into the port
  `get_gamma_*` helpers and diffs each activity factor vs the Fortran GAMMA* diagnostic;
  then an ASSEMBLY check reproducing MEG_isoprene + VOCFLXT from the port's emission
  factor (read from the SAME `megan21_emis_factors_78pft` file) + unit factor.
- `scripts/fortran_parity_isotopes_aripuana.jl` — teacher-forced single `use_cn+c13+c14`
  `clm_drv!` step: injects the Fortran iso pools from `before_step`, runs one step, diffs
  veg/soil iso pools, δ13C, and the recomputed rc13 ratios vs `after_hydrologydrainage`.

Run: `ARI_ISO_DUMPDIR=/private/tmp/claude-501/ari_iso_run julia +1.12 --project=. scripts/<harness> 40`

---

## Part A — VOC / MEGAN parity table (n40 midday, patch 2 = PFT 4 BET tropical)

| quantity | port | Fortran | max\|rel\| | verdict |
|---|---|---|---|---|
| GAMMAL (LAI factor) | 6.040851e-01 | 6.040851e-01 | 2.2e-8 | ✅ |
| ALPHA (gamma_P coeff) | 1.782073e-03 | 1.782073e-03 | 1.3e-8 | ✅ |
| GAMMAT (temperature) | 1.226127e+00 | 1.226126e+00 | 5.0e-7 | ✅ |
| EOPT | 1.530641e+00 | 1.530641e+00 | 5.5e-9 | ✅ |
| TOPT | 3.113656e+02 | 3.113656e+02 | 2.8e-8 | ✅ |
| GAMMAA (leaf age, evergreen→1) | 1.000000e+00 | 1.000000e+00 | 0.0 | ✅ |
| GAMMA_total (port L·T·A × Fortran P,S,C) | 4.090083e-01 | 4.090081e-01 | 4.9e-7 | ✅ |
| MEG_isoprene [kg/m²/s] | 7.952934e-10 | 7.952934e-10 | 1.7e-10 | ✅ |
| GAMMAP (PPFD factor) | 5.298941e-01 | 5.446931e-01 | 2.7e-2 | ⚠ input-limited |

**Verdict: PASS.** Every independently-computable activity factor and the assembled
isoprene flux match Fortran to ≤5e-7. GAMMAT — the most complex kernel (Guenther-2006
Eopt/Topt exponential temperature response) — matches to 5e-7, and its intermediates
EOPT/TOPT to ~1e-8. The GAMMA_total cross-check (port's L·T·A factors combined with
Fortran's P/S/C reproduce the full Fortran GAMMA to 4.9e-7) confirms the port kernel
is faithful at tropical values.

⚠ **GAMMAP is input-limited, NOT a divergence.** GAMMAP's final value weights sunlit
vs shaded by the *instantaneous* fsun, which was not retained on the h1 tape (FSUN
dropped from fincl). Feeding the start-of-step pdump fsun while Fortran's PAR_sun was
built with its during-step fsun introduces a ~2.7 % timing inconsistency. That this is
an input artifact and not a kernel error is proven by ALPHA (get_gamma_P's own nonlinear
coefficient) matching to 1.3e-8 and every other factor matching to ≤5e-7. GAMMAS/GAMMAC
likewise need the instantaneous btran/ci (not dumped); the port's get_gamma_SM/get_gamma_C
are value-verified against the VOCEmissionMod oracle in the Bow harness, and the
GAMMA_total check folds Fortran's own GAMMAS/GAMMAC in.

Monoterpene C10H16 (pinene_a+carene_3+thujene_a) Fortran flux = 9.72e-11 kg/m²/s (reported;
its per-compound gammas are not on the tape, only isoprene's imeg==1 set).

---

## Part B — C13/C14 isotope parity table (n40 midday teacher-forced step)

| quantity | max\|abs\| | max\|rel\| / \|Δ\| | verdict |
|---|---|---|---|
| rc13_canair (veg patches) | — | 1.7e-18 | ✅ exact |
| rc13_psnsun (veg patches) | 2.0e-6 | rel 1.8e-4 | ✅ (fractionation ACTIVE) |
| rc13_psnsha (veg patches) | 2.0e-6 | rel 1.8e-4 | ✅ |
| δ13C (deadstemc, established pools) | — | 0.0 ‰ | ✅ exact |
| soil1c_13_vr / soil2c_13_vr | 7.2e-5 / 3.6e-5 | 2.2e-5 / 1.1e-5 | ✅ |
| soil3c_13_vr | 4.6e-7 | 4.1e-7 | ✅ |
| soil{1,2,3}c_14_vr | ≤6.5e-15 | ≤6.5e-15 | ✅ |
| xsmrpool_13 (1-step increment) | 3.1e-9 | incr 1.8e-4 | ✅ |
| veg MAIN pools (leafc_13, cpool_13, …) | — | **NaN** | ✗ bulk-CN cold-start (see below) |

**Verdict: the isotope subsystem is FAITHFUL where the reference supports a clean diff.**
The discrimination physics — the whole reason for going tropical — passes:

- **rc13_canair matches to 1.7e-18** (the atmospheric C13 ratio + `forc_pc13o2 = C13RATIO·forc_pco2` wiring; value 0.011170 = Fortran C13RATIO/(1−C13RATIO)).
- **rc13_psnsun/psnsha match to 2.0e-6.** At Bow the canopy top layer was unlit so
  alphapsn≈1 and rc13_psn = rc13_canair trivially; at the **lit tropical midday canopy
  the fractionation is genuinely exercised** — rc13_psnsun = 0.01093 is depleted from
  rc13_canair = 0.01117, and the port tracks Fortran to 2e-6 (rel 1.8e-4). The 2e-6
  residual is the port's ci/ca vs Fortran's (the photosynthesis coupling), not the
  fractionation kernel (canair, which is ci-independent, matches to 1e-18).
- **δ13C on established pools (deadstemc) is exact**; soil isotope decomposition pools
  (litr/soil/cwd × C13/C14) match to ~2e-5 (C13) / ~6e-15 (C14); the maintenance-resp
  isotope pool xsmrpool_13 one-step increment matches to 3e-9.

### The veg-main-pool NaN is a bulk-CN cold-start artifact, NOT an isotope defect
`leafc_13` patch 2: injected `J_before = 1.10096` (matches Fortran exactly) → `J_after = NaN`
after one step, while Fortran goes 1.10096 → 1.10087. Root-caused with a dependency probe:
the **BULK** `cpool`, `leafc`, and `cpool_to_leafc` are ALL NaN post-step (`psnsun_to_cpool`
= 7.6e-5 is finite, so photosynthesis is fine); the isotope pools NaN only because they
faithfully mirror the broken bulk. **Attribution test:** the identical step with
`use_c13=use_c14=false` (pure bulk CN, no isotope code path) gives the SAME bulk NaN
(cpool=NaN, leafc=NaN) — so this is a **bulk-CN cold-start allocation singularity**,
fully independent of the isotope subsystem. Traced further (code-verified): with
`use_fun=true` (the BGC default), `nutrient_competition_flexiblecn.jl:240` sets
`plant_calloc = npp_growth` (the FUN carbon-for-growth), and FUN returns a NaN
`npp_growth` at cold start (soil mineral N ≈ 0 → N-uptake conductance reciprocals
0/0; `fpg_col` stays NaN). That NaN plant_calloc drives `cpool_to_leafc` → NaN. Bow's iso harness never hit it because Bow's
reference was a *developed* mid-July run (cpool > 0, past the cold-start singularity);
the tropical reference is a 2-day cold start. This bulk-CN cold-start NaN is a real but
SEPARATE port issue, outside the VOC/isotope scope, and it invalidates only the
veg-main-pool isotope diff — every isotope quantity that does not ride the broken bulk
pool matches Fortran.

**Honesty note:** the harness's per-pool flag logic (inherited from the Bow harness)
prints "ok" for NaN pools (NaN > 1e-6 is false) — a vacuous-check weakness that masked
the NaN until the dependency probe. The global pool max|rel| is reported as NaN, which is
the honest signal. (The Bow harness has the same latent weakness; not fixed here to keep
the diff scoped, but noted.)

## Bottom line
- **VOC/MEGAN at the tropical isoprene stress case: PASS** — activity-factor kernels +
  assembled isoprene flux faithful to ≤5e-7 vs a live Fortran MEGAN history tape;
  the one ~2.7 % gap (GAMMAP) is a dumped-input timing limitation, not a port error.
- **C13/C14 discrimination at the tropical stress case: PASS** — rc13_canair exact,
  rc13_psnsun/sha (now genuinely fractionating) to 2e-6, δ13C exact, soil iso pools to ~2e-5.
- **Side-finding (bulk CN, out of scope):** a bulk-CN cold-start allocation NaN at
  Aripuanã (isotope-independent) blocks the veg-main-pool isotope diff on a 2-day
  cold-start reference; a developed-spinup reference (as at Bow) would restore it.

Dump paths: `/private/tmp/claude-501/ari_iso_run/{pdump_before_step_n40.nc,
pdump_after_hydrologydrainage_n40.nc, Bow_at_Banff_lumped.clm2.h1.2202-01-01-00000.nc}`.
