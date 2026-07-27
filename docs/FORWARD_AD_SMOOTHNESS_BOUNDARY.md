# CLM.jl's forward-AD boundary is *smoothness*, not type-genericity

**Thesis.** ForwardDiff propagates cleanly through the *entire* per-step CLM.jl physics —
hydrology, the full snow-water cycle (accumulation, melt, compaction, layer combine/divide),
and the surface energy balance. There is **no type wall**. But a clean forward pass proves only
that the code is *element-type-generic* (`Dual`-clean); it does **not** prove the gradient is
correct. The gradient is correct only where the model is a *smooth* function of its parameters,
and CLM's snow physics is full of discrete events (`imelt` freeze/thaw flags, `combine`/`divide`
layer changes) that make the streamflow loss **piecewise-smooth**. So CLM.jl's real forward-AD
boundary is **differentiability smoothness at those discrete events** — a mathematical property
of the model, not a software property of the port. This document states the distinction, gives
the measured evidence (independently reproduced), and names the remedy and the lesson.

This corrects two overstated claims made earlier in `AD_STREAMFLOW_CALIBRATION.md` (§3b): that
the full-year hydrograph is "a large engineering assembly, not a feasibility barrier," and that
"Farquhar photosynthesis is Dual-clean." Both are revised below.

---

## 1. Two independent meanings of "differentiable"

A correct AD gradient requires **both** of the following. They are routinely conflated; this
whole investigation is the story of conflating them.

**(A) Type-genericity — "Dual-clean" — a property of the *code*.**
The code must accept and propagate `ForwardDiff.Dual` numbers without erroring: no `Float64(x)`
casts, no `::Float64` field pins, no type-narrowing containers. This is a *software* property,
removable by making structs element-type-generic. **CLM.jl satisfies it across the whole per-step
physics** (proven with zero `src` changes via a harness-side `Dual` retype — see §3b of the
calibration doc). This is the wall we kept looking for and kept correctly knocking down.

**(B) Mathematical smoothness — "Dual-smooth" — a property of the *function*.**
AD computes the exact derivative *of whichever branch executes*. If the underlying mathematical
function has a kink (slope discontinuity) or a jump (value discontinuity) at the operating point,
AD returns a one-sided (sub)derivative that does **not** equal the true derivative — because there
isn't one. **No amount of type-plumbing fixes this.** You cannot type-refactor past a discontinuity.

**The trap.** A forward pass that runs clean is evidence of (A) and says *nothing* about (B).
"It ran without erroring" is type-genericity masquerading as a valid gradient. Every over-claim
in this arc was this substitution.

---

## 2. Where CLM's snow physics violates smoothness

Two discrete constructs, both fundamental to the model, both non-smooth in the parameters:

- **`imelt` — the per-layer freeze/thaw flag.** Each step, each layer, the phase-change routine
  decides *discretely* whether a layer melts/freezes by comparing temperature to `tfrz`. A
  parameter nudge that pushes a layer across `tfrz` **flips the branch** — energy routed to phase
  change vs. sensible heating changes discontinuously. The loss gains a small jump at every
  threshold crossing.

- **`combine_snow_layers!` / `divide_snow_layers!` — changing `snl`, the state *dimension*.**
  A thinning layer is merged with a neighbour (combine); a thick layer is split (divide). At the
  threshold the **number of snow layers changes** and mass/energy is repartitioned across a
  different-length state vector — a genuine discontinuity in the state trajectory as a function
  of the parameters.

---

## 3. The measured evidence (independently reproduced)

`scripts/forwarddiff_fullyear.jl` (`CLM_FY_STAGE=sweep CLM_LAYERS=1`) sweeps `θ_snow` over a
range that drives the Krycklan pack through combine/divide and `imelt` flips, and at each `θ`
compares the **ForwardDiff gradient** of the streamflow loss against **central finite difference**
(FD re-runs the whole perturbed loop, so FD is the ground truth for the *actual*, kinked function).
It flags a **straddle** when the `snl` trajectory differs across the FD bracket `θ ± h`.

| regime | FwdDiff-vs-FD rel error | interpretation |
|---|---|---|
| transition-free control (no-melt hydrology) | **~1e-8** | AD **is** the exact derivative where the function is smooth |
| smooth melt steps (`imelt` active, `snl` fixed) | **~1e-3** (6e-4 – 1.4e-3, pervasive) | discrete `imelt` determination injects a per-step kink-error |
| **straddled `snl` transition** (combine/divide fires within `θ ± h`) | **11 – 85 %** (my sweep worst **2.06e-1** at `θ=0.907`, `snl −5→−1→0`) | ForwardDiff returns the one-sided subgradient of the active branch; FD straddles the kink |

Reproduced on two independent runs (agent + coordinator). Worst over a 28-point sweep:
**rel 2.06e-1 at the `snl` straddle**; a second spike **1.08e-1** at `θ=4.574` (a within-bracket
combine). The magnitude at a straddle depends on how squarely `θ ± h` brackets the discontinuity
(the agent's finer-resolved straddle reached 85%).

**Crucially, the pervasive ~1e-3 is the same residual** earlier dismissed as "FD kink-error" in the
90-day gate. It was never FD's error — it is the `imelt` discreteness, and it is present on *every*
melt step. The transition-free control at ~1e-8 proves the ~1e-3 is genuinely the discrete events,
not the AD.

---

## 4. Why this is a graceful, density-dependent boundary — not a cliff

The forward-AD gradient degrades *gracefully* with transition density:

| transition density | gradient quality | evidence |
|---|---|---|
| none (smooth region) | **exact** (~1e-8) | the AD gradient equals the derivative |
| sparse (short melt window, one freshet) | **usable** — kinks are measure-zero, subgradient descent tolerates them | the earlier split-samples **calibrated and generalized**: Stillwater eval-KGE **+0.585**, Krycklan **+0.144** |
| dense (a full year: winter divides + spring combines + thousands of diurnal `imelt` flips) | **piecewise-valid but kink-dominated near transitions** | the sweep above; a year fires hundreds–thousands of transitions |

The middle row is the reconciliation the earlier over-claims skipped: **the kinked gradient
empirically works on sparse-transition problems** — gradient descent calibrated real basins to
positive held-out KGE across melt. It is *usable*, just not globally C¹. Robustness on a
*full-year* (dense-transition) problem is where smoothing becomes necessary.

---

## 5. The remedy — and it is already the plan

This is the standard problem in differentiable physics (differentiable simulators, DiffTaichi,
JAX/PyTorch physics engines). Standard remedies, in CLM.jl terms:

1. **Smoothing / relaxation (the Phase-3 plan).** Replace hard `min`/`max`/threshold/step with
   smooth surrogates. CLM.jl already ships `smooth_max`/`smooth_min` wrappers — they currently
   reduce to *exact* `min`/`max` on the Float64 path (so defaults are byte-identical), but they
   exist precisely to be switched on at the `imelt` and combine/divide thresholds, making the loss
   C¹ and the gradient globally valid.
   **The catch (documented):** the smoothing constant *k* is **dimensional** — `smooth_max(0,x)`
   overshoots by `log(2)/k` *in the axis's units*. Smoothing trades a kink for a small, controllable
   physics bias; choosing *k* is a modeling decision, not a porting fix.
2. **Subgradient descent (no code change).** Use the one-sided gradient AD already gives. Works for
   piecewise-smooth functions (kinks are measure-zero); this is *why* the sparse-transition
   calibrations succeeded. Cheapest option; degrades on dense-transition problems.
3. **Straight-through / fixed-branch estimators.** Hold the discrete choice fixed for the derivative
   (literally what ForwardDiff does at a kink) — biased but usable.

---

## 5b. Phase-3 executed — the `imelt` label smoothed, quantified

Two things sharpen §3/§5 after actually doing the work.

**(i) The continuous smoothing already existed, and a `:auto` artifact inflated the residual.**
`soil_temperature.jl`/`snow_hydrology.jl` already smooth the *continuous* heat/mass partitions
(gated by `SMOOTH_MODE`); the authors deliberately left the integer `imelt` **label** and the
combine/divide **layer count** discrete. My original sweep ran under `SMOOTH_MODE=:auto`, where
ForwardDiff (Dual) smooths but central-FD (Float64) does not — so part of the "~1e-3 pervasive"
was a **smooth-vs-hard mismatch**, not non-smoothness. Re-run under `:always` (both paths on the
same smoothed function), the settled-melt residual tightens to **~2–6e-4** (FD-truncation-limited),
and the genuine non-smoothness concentrates at the discrete-label flips.

**(ii) The `imelt` melt label is now smoothed (the last mile).** The snow melt identification
`t_soisno > tfrz` (with ice) — a hard threshold — is replaced, under a smooth mode, by a C¹ melt
driver `tinc = −smooth_max(0, t − tfrz; k = PHASE_CHANGE_TEMP_K)`. The key is the **axis**: this is
a **temperature in Kelvin**, so width `log(2)/k` at `k=1e3` is `6.9e-4 K`, and the `smooth_max`
saturation guard returns **exactly** 0 for `t < tfrz − 0.036 K` — a genuinely cold layer is
byte-identical to the hard path; only a ~0.036 K band around freezing is rounded, a negligible
melt-mass bias (unlike the kg/m² mass axis, where a comparable width would misassign 0.0139 mm of
ice). Default (`:auto` + Float64) runs the original hard branch verbatim.

**Measured effect** (`CLM_SMOOTH=always` sweep, before → after the label smoothing):

| θ_snow | before | after | what it is |
|---|---|---|---|
| 0.500 | **2.44** (sign-flipped) | **0.015** | pure `imelt` label flip → **fixed** |
| 1.926 | **0.427** | **0.004** | pure `imelt` label flip → **fixed** |
| 0.704 / 0.907 | 0.9 / 0.05 | 0.9 / 0.6 | `snl` combine/divide straddle → **remains** |
| settled melt | 2–6e-4 | 2–6e-4 | unchanged (FD-truncation-limited) |

So smoothing the temperature-threshold label removes the label-flip spikes; **the residual
non-smoothness is now purely the `snl` combine/divide layer-count change** — the genuine
dimension-change discontinuity that cannot be smoothed without a structural fixed-layer
reformulation. Validated: **507 default-path tests byte-identical** (`_pc_smooth(T)` is false for
`Float64`+`:auto`); the melt driver is C¹ across `tfrz` under Dual and `→ hard as k → ∞`
(`test_phase_snow_smoothing.jl`). Freeze and soil-layer labels are left hard (same technique
applies; they do not drive the spring-melt streamflow gradient).

## 6. The photosynthesis claim, corrected

The earlier "Farquhar photosynthesis is Dual-clean" claim was **unproven — potentially vacuous**.
The 6.57e-7 gate was `d(Σ t_veg)/d(air-temp)`, which flows through the surface **energy balance**
(genuinely non-degenerate `Dual`-clean — that part stands) but **not demonstrably through the
photosynthesis solve**. A direct `d(GPP)/d(vcmax)` and `d(GPP)/d(absorbed-PAR)` gate returned
**0.0 (vacuous)**: `psnsun` was ~0 in the probe state, so the derivative never entered the Farquhar
Newton solve. **Honest status: the canopy energy balance is non-degenerately `Dual`-clean; whether
Farquhar photosynthesis is non-degenerately `Dual`-clean is UNCONFIRMED** (it runs without a
type-break, but a non-zero psn gradient — requiring a genuinely photosynthesizing state — was not
produced). "It ran clean" was, again, not sufficient evidence.

---

## 7. The lesson

**"It ran clean" ≠ "the gradient is correct."** A forward pass propagating `Dual`s without erroring
proves *type-genericity*, not *mathematical smoothness* — and the gradient is only as valid as the
function is smooth at the operating point. To trust an AD gradient you must gate it against finite
differences **at the operating point**, especially **across the discrete events** a physical model
is full of, **with a smooth-region control** to prove any error is the transitions and not the AD.
This is the same class as the port's other hard-won lessons — *conservation is not accuracy*,
*vacuous checks*: assert the **value** (the gradient, across a transition), never that it merely ran.

**Bottom line.** CLM.jl has **no forward-AD type wall** — the whole per-step physics is `Dual`-clean.
Its forward-AD boundary is **smoothness**: the discrete snow events make the streamflow loss
piecewise-smooth. The continuous heat/mass partitions were already smoothed; **Phase-3 now also
smooths the `imelt` melt label** on the temperature axis (negligible bias, default byte-identical),
removing the label-flip gradient spikes (§5b). The gradient is **exact away from transitions and
usable across sparse ones (and empirically calibrates real basins)**; the **one remaining
discontinuity is the `snl` combine/divide layer-count** — a genuine dimension change that needs a
structural fixed-layer reformulation, not a threshold smoothing. `Dual`-clean is necessary but not
sufficient; smoothness is the deeper requirement, and it is now closed except for the integer
layer count.
