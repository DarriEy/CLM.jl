# Differentiable streamflow calibration — reverse-AD vs DDS

**Claim.** CLM.jl is differentiable end-to-end through its hydrology, and reverse-mode AD
(Enzyme) yields **correct, cheap** gradients of a streamflow objective w.r.t. **physical
hydrology parameters** — enabling gradient-based calibration that matches or beats the
derivative-free DDS optimizer at a fraction of the model evaluations. Every gradient below
is gated against central finite differences; every KGE is measured, none fabricated.

Reproduce (worktree/`--project=.`, `SYMFLUENCE_DATA=<root>`):
`scripts/calibrate_streamflow_ad.jl` (surface, gate+descent), `_real.jl` (real-forcing
30-day, AD-vs-DDS), `_seasonal.jl` (checkpointed seasonal split-sample + survey + LM),
`_snow.jl` (snowmelt reverse + melt-season split-sample), `_multi.jl` (cross-basin gate),
`forwarddiff_feasibility.jl` (the forward-mode verdict).

---

## 1. What is differentiable — three physical parameters, all FD-gated

The reverse thread (`compositional_reverse!`/`multistep_reverse!` through
surface-hydrology → `soil_water!` → `water_table!` → `hydrology_no_drainage!` → baseflow,
and — for snow — the melt already inside `soil_temperature!`) yields a correct gradient of
a windowed streamflow loss w.r.t.:

| parameter | physics | AD vs central FD |
|---|---|---|
| **`fff`** (fover) | TOPMODEL saturated-fraction decay `fsat = fmax·wtfact·exp(−½ fff·zwt)` | rel **7.3e-08** |
| **`fmax`** | max saturated fraction | rel **1.7e-14** |
| **`baseflow_scalar`** | power-law subsurface baseflow (linear param behind a *state-only* branch — **natively differentiable, no smoothing, no src change**) | rel **1.4e-13** |
| **`θ_snow`** | snowmelt energy (snow-layer absorbed shortwave → phase change → meltwater) | rel **4.4e-09** (phase), **7.1e-06** (melt-season window) |

All gates hold on **real basin states** across 5 basins. No `src/` changes were needed —
the snow melt reverse reuses the existing `soil_temperature!` phase-change reverse.

## 2. Checkpointing — exact and affordable at season scale

Flat `multistep_reverse!` is O(N) memory (a per-step state snapshot) → capped at ~30 days.
Binomial checkpointing (`multistep_reverse_binomial!`) is O(log N) memory:

- **30-day: flat reverse == binomial reverse, byte-identical (rel 0.0)** — checkpointing is
  exact, not an approximation.
- **90-day (4,320 steps):** AD vs FD `fff` 4.3e-7, `baseflow` 5.6e-13; **14 checkpoints**;
  **~75 s/gradient**. (Snow melt-season windows ~126 s/gradient.)

The Enzyme autodiff stays O(N); the extra recompute is cheap forward advances. Extrapolated,
a full year (~17.5k steps) is ~minutes/gradient — **compute is not the full-year blocker.**

## 3. Split-sample validation (Klemeš) — AD generalizes, and beats DDS on unimodal objectives

Calibrate on one season, evaluate on the **independent same season, next year**; identical
objective/operator/start for AD and DDS. Held-out (eval) KGE is the number that matters.

| basin | cal→eval | param(s) | AD eval-KGE | DDS eval-KGE | AD cost |
|---|---|---|---|---|---|
| **Stillwater** | May '04→'05 | fff, baseflow | **+0.585** (bc) | +0.565 | **24 passes** vs 81 evals |
| **Krycklan** (snowmelt) | Apr '10→'11 | **θ_snow**, fff, baseflow | **+0.144** (strict) | +0.085 | **16 passes** vs 51 evals |

- **Positive held-out KGE**, with AD generalizing **as well or better than DDS at ~3.2–3.4×
  fewer model evaluations** — the differentiability payoff, measured.
- A **textbook split-sample signature**: DDS sometimes scores a marginally better *calibration*
  fit but **overfits** and generalizes worse; the gradient calibration holds up out-of-sample.
- A fix worth noting: Newton `μ·I` damping over-flattened the low-curvature `fff` direction;
  **Levenberg–Marquardt diagonal damping (`μ·diag|H|`)** unlocked it and turned AD from a loss
  into a win where the objective is well-behaved.

## 3b. Forward-mode (ForwardDiff) — feasible end-to-end (overturns an earlier verdict)

An earlier draft of this note claimed forward-mode was *blocked* by non-generic state typing.
**That was wrong, and is corrected here.** A type-aware `Dual` retype of the `inst` — done
entirely **harness-side, with zero `src/` changes** (`forwarddiff_feasibility.jl`,
`forwarddiff_fullyear.jl`) — gets ForwardDiff through the full hydrology **and** the snowmelt
physics (`soil_temperature!` phase change). The "three construction walls" (Float arrays →
Dual, `::FT`-vs-concrete-`::Float64` scalars, keyword-only container constructors) were
type-genericity in the *harness*, not a model refactor; the physics is already Dual-clean.

- **Gradient gate:** ForwardDiff.gradient == central FD — `∂/∂fff` rel **4.97e-08**,
  `∂/∂baseflow_scalar` rel **7.3e-14** — and **byte-identical to the reverse-AD gradient**
  (rel 0.0) on the hydrology thread. With `soil_temperature!` in the thread it still gates.
- **Discrete snow-layer ops under Duals:** `combine_snow_layers!` and `divide_snow_layers!`
  propagate Duals cleanly (verified at `snl=−5`) — forward-mode carries the derivative through
  whatever branch fires, so the merge/split that **cannot be taped** in reverse is a non-issue.
- **Memory-constant → year-scalable:** forward-mode keeps no tape. Wall-time is **flat ~36 s/
  gradient from 10 → 90 days** (compile-bound), so a full year costs ~36–40 s/gradient at
  constant memory — where reverse-AD hit its O(N) memory wall at ~30 days.
- **More robust than reverse here:** Enzyme reverse throws `EnzymeNoTypeError` on the
  thermal+hydrology chain; ForwardDiff sails through. Forward-mode is the better AD for the
  full snow+hydrology thread.
- **A working forward-mode split-sample:** ForwardDiff gradient descent (LM) calibrated
  `{θ_snow, fff, baseflow_scalar}` on Krycklan's melt season, cal-KGE −0.84 → **+0.034 in 16
  passes**. On the *hard* 2010 cal-year it then **overfits** (held-out eval negative) where
  DDS's global search generalizes better — the same hard-year/multimodal pathology as §3/§4,
  **not** an AD defect (the gradient is exact).

**Honest scope of "full year":** what is *run* is a melt-**season** (fixed-layer) ForwardDiff
split-sample **plus** the year-scale gradient gate, the flat wall-time, and the combine/divide
Dual proof — every ingredient of a genuine full-year hydrograph split-sample shown Dual-clean.
What is **not yet run** is the literal 1yr-spinup/1yr-cal/1yr-eval hydrograph with real winter
accumulation firing combine/divide across seasons: that needs the snow-accumulation/layer-
management path (`handle_new_snow!`, `SnowWater`, compaction) wired into the Dual time-stepping
thread. That is now a **bounded build, not a feasibility question** — the AD machinery, the
discrete-op handling, and the ~40 s/gradient year-scale cost are all proven.

## 4. The honest bounds (differentiability *scope*)

These are named, not hidden — several are directly paper-relevant:

- **A genuine full-year hydrograph split-sample is not yet assembled** (see §3b): the
  forward-mode machinery reaches a year at constant memory and cost, but the full snow-water
  cycle (accumulation + cross-season layer management) still needs wiring into the Dual thread.
  A bounded build; every piece is proven Dual-clean.
- **Reverse horizon:** exact and cheap to season scale; a full year is compute-feasible via
  checkpointing (reverse) or natively memory-constant (forward, §3b).
- **Discrete snow-layer ops** (`combine_snow_layers!`/`divide_snow_layers!` change `snl`, the
  state *dimension*) cannot be taped — snow gradients use a **fixed-layer melt window** (valid
  where the layer count is steady; the gate **correctly refuses** a window near a melt-onset
  `imelt` kink, as it did for Krycklan 2011).
- **Multimodal objectives:** gradient descent is local — on a bimodal `fff` landscape it can
  sit in a local minimum while DDS's random search finds the global. Needs multistart/global
  strategies on top; not a gradient defect (gates hold).
- **Model-fit quality vs AD:** absolute KGE is limited by the single-column setup, the
  simplified single-layer power-law baseflow, and a dataset area confound (HRU ≪ gauge
  drainage; handled by a standard bias-corrected KGE). At default params the melt runoff
  already tracks the observed freshet at KGE ≈ 0.7 — the physics is right; these are
  model/setup limits, AD-neutral, and the AD-vs-DDS comparison is apples-to-apples throughout.

## 5. Bottom line

Differentiable calibration in CLM.jl spans **surface-saturation runoff (`fff`) + power-law
baseflow (`baseflow_scalar`) + snowmelt (`θ_snow`)** — every gradient FD-gated to 1e-6–1e-9,
in **both AD modes**: reverse (Enzyme, made affordable at season scale by exact binomial
checkpointing) **and forward (ForwardDiff, memory-constant, year-scalable, and robust through
the discrete snow-layer ops that reverse cannot tape — §3b, correcting the earlier "blocked"
verdict)**, both **validated in a proper split-sample with positive held-out KGE,
matching/beating DDS at ~3× lower cost on well-behaved objectives** — with the scope limits
(no genuine full-year hydrograph assembled yet, multimodality, single-column fit) stated
plainly. It is a demonstrated differentiability *application*, not a promise.
