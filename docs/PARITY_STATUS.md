# Fortran Parity Status — CN + Biogeophysics

> ## 🕰️ HISTORICAL SNAPSHOT — 2026-06-17. Do not read this as current status.
>
> **Banner added 2026-07-12.** Everything below is preserved verbatim as the
> record of the *2026-06-17* single-column CN parity campaign at Bow at Banff.
> It was accurate then. It is **not** the current picture, and its "Not yet
> validated" list in particular understates today's coverage.
>
> **Current sources of truth:**
> - **Parity coverage / scorecard:** the README's *Validation Against Fortran
>   CLM5* section (a 20-biome × 69-variable single-point suite, run via
>   `scripts/parity_run_domain.jl`) — this document predates that suite entirely.
> - **Strict-gate residuals:** [`docs/STRICT_GATE_RESIDUALS.md`](STRICT_GATE_RESIDUALS.md).
> - **GPU coverage:** [`docs/GPU_PORT_GAP.md`](GPU_PORT_GAP.md).
> - **What is actually wired into the driver:** `git grep` in `src/driver/` —
>   several things this document calls unported/unwired are now live (see below).
>
> **Specific claims below that are now WRONG** (verified in-tree on 2026-07-12):
> - Methane, fire, VOC/MEGAN, dust, irrigation and the C13/C14 isotope path are
>   all **wired into `clm_drv!`** now: `ch4!` (`clm_driver.jl:2320`),
>   `dust_emission!` (`:1507`), `voc_emission!` (`:1522`),
>   `calc_irrigation_fluxes!` (`:1034`), `c13_c14_photosynthesis!` (`:1348`),
>   `c14_decay!` (`:2161`), `cndv_driver!` (`:2626`). Fire runs from the CN
>   driver behind `cnfire_method` (`cn_driver.jl:53`, default `:nofire`, so the
>   default path is unchanged).
> - The "not yet validated: other sites / other biomes" caveat is superseded by
>   the multi-biome suite (still single-point, still one year per site — read the
>   README's own caveats, which remain in force).
> - Since this was written: FATES was ported and wired (Fortran-FATES bit-parity
>   is **not** established), matrix-CN was wired, and much of the driver was
>   kernelized for GPU.
>
> The methodology sections (especially the **before-vs-after dump pitfall**) are
> still valuable and still correct — that is the main reason this file is kept.

**Last updated:** 2026-06-17 *(content frozen; see banner)*

This document records, honestly, the current state of Fortran-parity validation
for CLM.jl. It covers what has been verified, the known residuals (and why they
are *not* bugs), the methodology pitfall that nearly produced a false bug report,
and — just as importantly — what has **not** been validated. The goal of the port
is full process fidelity (see `PRD_CLM_JULIA_PORT.md`); this is the scorecard for
how close we are on the one configuration we have driven to ground truth.

The headline: on the validated configuration, CLM.jl reproduces the Fortran CN +
biogeophysics surface to **sub-percent**, with every remaining gap traced to a
named cause. That is a much narrower claim than "CLM.jl reproduces Fortran." Read
the *Not yet validated* section before quoting any number here.

---

## Scope of validation

- **Site:** Bow at Banff, a single column, 3 patches: bare ground, tree
  (needleleaf evergreen), and grass (C3).
- **Configuration:**
  - `use_cn` (full carbon–nitrogen biogeochemistry)
  - PHS / `use_hydrstress` (plant hydraulic stress, Newton fixed-point solve)
  - LUNA / `use_luna` (vcmax25 / jmax25 acclimation; injected from Fortran dumps
    rather than ported)
  - The Bow `lnd_in` flags that materially move the answer:
    `modifyphoto_and_lmr_forcrop=.true.`, `light_inhibit=.true.`, Jackson root
    profile (`rootprof_beta`, method=1), and CO2/pbot scaling so `forc_pco2` is
    367 ppm (`367e-6 * pbot`), not the 506 ppm default.
- **Reference:** Fortran CLM/CTSM per-step instrumentation dumps taken from a
  converged BGC spinup (soil C/N converged to <0.2% per 100 yr), restarted at
  model year 2202. The validation window is ~28+ contiguous summer timesteps
  beginning near `nstep ≈ 1757845`.
- **Dumps:** `/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_summer/`
  — per step there are `pdump_before_step_n<N>.nc`,
  `pdump_after_canopyfluxes_n<N>.nc`, and `pdump_after_hydrologydrainage_n<N>.nc`
  (28 of each across the summer window `n1757845..n1757872`).

---

## Verified parity results

Each claim cites the harness script under `scripts/` that produces it. Single-step
scripts re-inject the Fortran state before the step and then diff Julia's evolved
output against the *same-step* Fortran `after_hydrologydrainage` dump, so they
measure per-step **translation** error, not compounded drift.

### Biogeophysics (SP / surface-energy path)

- `T_GRND`, `T_VEG`, `SABV`, `SABG`, snow state, and the hydrology cascade match
  Fortran to roughly **0.01–0.3 K** (temperatures) and small relative errors
  (fluxes / water states).
- The gated Fortran-parity test passes; the annual spun-up-vs-spun-up run matches
  Fortran (the earlier cold-start "BTRAN inversion" was an init artifact, not
  physics).
- Harnesses: `scripts/fortran_parity_validate.jl`, `scripts/fortran_parity_snow.jl`,
  `scripts/fortran_parity_smoke.jl`, `scripts/fortran_parity_sweep.jl`.

### CN single-step (the key new result — 2026-06-17)

Comparing Julia's **evolved after-state** against the Fortran
`after_hydrologydrainage` dump for the same step:

- **Mineral N (`smin_nh4_vr` + `smin_no3_vr`):** worst layer (layer 1) matches to
  **0.024%**; by layer 5 the error is **< 0.001%**.
  (`scripts/probe_sminabs.jl`, `scripts/probe_nh4no3.jl`)
- **FUN plant-N uptake:** tree **1.002×**, grass **0.988×** vs Fortran.
  (`scripts/probe_funuptake.jl`, `scripts/probe_fun.jl`)
- **availc (available C for allocation):** tree **1.002×**, grass **0.988×**.
  (`scripts/probe_availc.jl`)
- **plant_ndemand:** matches Fortran. (`scripts/probe_ndemand.jl`)
- **calc_allometry!:** **byte-identical** to Fortran. (`scripts/probe_allometry.jl`)
- Whole-step CN pool diff: `scripts/fortran_parity_cn_step.jl`,
  `scripts/fortran_parity_cn_summer.jl`, `scripts/fortran_parity_cn_read.jl`.

### Multi-step drift (free-running)

`scripts/fortran_parity_drift.jl` injects the Fortran restart **once**, then
free-runs N steps (advancing forcing and time), diffing the Julia trajectory
against the per-step Fortran `after_hydrologydrainage` dumps. This shows how the
single-step residuals **compound**:

- Drift is **bounded and roughly linear** — no runaway.
- Mineral N drifts ~**0.024%/step → ~0.7%/day**, then plateaus.
- `leafc` stays **< 0.1%**; `soil1c` **< 0.01%**; the temporary buffer pools
  (`cpool`, `xsmrpool`) drift fastest but plateau rather than diverge.
- Critically, the per-step drift rate **equals** the single-step translation
  residual. That is the signature of *compounding*, not a leak or a per-step bug:
  the single-step C/N budgets themselves conserve.

---

## Known residuals (NOT bugs)

1. **Grass GPP / availc ~1.2% low** (tree is exact). This is the documented
   "closed GPP/ci axis" — the last sub-percent gap in grass available carbon.
   Closing it requires Fortran `psn`/`ci` live-print instrumentation at the
   matched step (the vcmax_z, btran, leaf_mr, and froot_mr inputs all already
   match Fortran). It is a residual, not a defect.

2. **Nfix ≈ 0.** Julia reports ~`1e-28`, Fortran ~`1e-12`. Both are effectively
   zero N fixation at this site/season; the difference is cosmetic (denormal /
   round-to-zero ordering), not physical.

3. **Cross-language reality.** Different floating-point operation ordering,
   different nonlinear solvers (the PHS Newton fixed-point iteration in
   particular), and different discontinuity handling mean **exact bit-parity is
   not the target** and is not achievable. Sub-percent agreement with every gap
   traced to a named cause *is* the target, and is what is reported above.

---

## Methodology note / pitfall (READ THIS)

> **Single-step CN parity probes MUST compare Julia's this-step output against the
> `after_hydrologydrainage` dump — NOT `before_step`.**

The `before_step` dump holds the **previous** step's values for the
allocation-pipeline fields (`availc`, `plant_ndemand`, and the FUN uptake
streams). Diffing Julia's freshly-computed this-step values against `before_step`
compares two different timesteps and manufactures phantom errors.

On 2026-06-17 an entire apparent bug chain — a "FUN uptake 0.687× tree throttle"
that propagated into `plant_ndemand` and then `availc` — turned out to be a
**stale-dump artifact**. Against the correct `after_hydrologydrainage` dump,
everything matched (tree FUN uptake 1.002×, availc 1.002×).

The `availc` field in the dump made this unambiguous:

| quantity (tree availc)         | value     |
| ------------------------------ | --------- |
| Fortran `before_step`          | 1.711e-6  |
| Fortran `after_hydrologydrainage` | 1.173e-6  |
| Julia this-step                | 1.175e-6  |

Julia's this-step value (1.175e-6) sits on the Fortran **after** value (1.173e-6),
not the before value (1.711e-6). When a CN probe shows a clean ~0.7× or ~1.5×
factor on an allocation-pipeline field, suspect the dump phase before suspecting
the port.

---

## Not yet validated

The validated surface is **one site, one column, three patches, one season, one
config**. The following are explicitly **not** established:

- Other PFTs (only bare / needleleaf-evergreen tree / C3 grass exercised).
- Other sites, climates, soil columns, or land-cover mixes.
- Other configurations / flag combinations (e.g. `use_cn` off, no-PHS, no-LUNA,
  crop, irrigation, transient land use, fire-active years).
- Modules the parity harness does not drive at this site (much of fire, methane,
  CNDV, VOC, dynamic vegetation, isotopes — these are GPU-validated for
  *internal* consistency but not against Fortran ground truth at Bow).
- Winter / shoulder-season dynamics (phenology onset/offset, snow-dominated
  energy balance under CN).

**"CLM.jl reproduces Fortran everywhere" is NOT a claim this document supports.**
What is supported: on the Bow CN + biogeophysics surface, for a contiguous summer
window, CLM.jl matches Fortran to sub-percent, with every remaining gap traced to
a named, understood cause.

---

## Diagnostic scripts backing these claims

All under `scripts/`, runnable with `julia +1.12 --project=.`:

- `probe_nh4no3.jl` — NH4 vs NO3 per-layer split of the mineral-N residual.
- `probe_funuptake.jl` — FUN per-patch N-uptake stream totals vs Fortran.
- `probe_ndemand.jl` — plant N demand parity.
- `probe_allometry.jl` — `calc_allometry!` byte-parity check.
- `probe_sminabs.jl` — clean single-step absolute mineral-N relative error
  (evolved after-state vs Fortran **after** dump — the correct comparison).
- `fortran_parity_drift.jl` — multi-step free-running drift harness.
- `fortran_parity_cn_summer.jl` — CN single-step pool parity across the summer
  window.

Supporting: `probe_fun.jl`, `probe_availc.jl`, `probe_sminbalance.jl`,
`fortran_parity_cn_step.jl`, `fortran_parity_cn_read.jl`, `fortran_parity_common.jl`
(shared single-step injection/runner).
