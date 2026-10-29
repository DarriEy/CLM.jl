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

## 5c. Phase-3 continued — the integer layer count `snl` (fixed-dimension, Stage 1)

`snl` — the number of snow layers — is the last discrete residual, and it is *categorically*
harder than the `imelt` label: it is a **dimension**, not a threshold on a continuous quantity,
so no sigmoid smooths it. Diagnosis first (this mattered): the melt-time `snl` reduction is **not**
the thickness/`dzmin` combine — it is the **ice-mass removal** in `combine_snow_layers!`, which
drops a layer once its ice falls below `0.01 kg/m²` (a *finite* mass). A three-order-of-magnitude
sweep of the collapse threshold (`edz` = 1e-5 → 1e-12) was byte-identical until *that* trigger was
found and gated — a reminder to locate the actual discontinuity, not the assumed one.

**Approach (fixed-dimension without the 621-site rewrite):** hold `snl` *fixed* through the melt
window by making all three removal triggers (ice-mass, thickness/`dzmin`, and the `snl→0`
all-snow-gone switch) wait, under a smooth mode, until the layer is **essentially empty** (ice or
thickness `< SNOW_COMBINE_EMPTY_DZ`, 1e-8). Removing a ≈0 layer is continuous — the same
"discrete event at zero-contribution" principle as the `imelt` fix — and because `snl` simply
stays fixed, none of the 621 `snl`-bound loops need to change.

**Result (`CLM_SMOOTH=always` sweep):**

| θ_snow | before | after (fixed-`snl`) | what it is |
|---|---|---|---|
| 0.500 | 1.5e-2 | **7.1e-5** | thin-layer removal → continuous |
| **0.907** | **5.9e-1** | **2.9e-3** | the **layered→bulk mode switch** → continuous |
| settled melt | 3–5e-4 | 3–5e-4 | unchanged |
| **0.704** | 9.3e-1 | **8.7e-1** | **whole-pack melt-out** — remains |

Two structural risks cleared: the **thermal solve survives** the resulting persistent thin layers
(no NaN — a naive fixed-dimension change would divide by dz→0), and the **default path stays
byte-identical** (217 snow/robustness tests; gates are behind `_pc_smooth(T)`, false for
`Float64`+`:auto`). The sharpest spike — the `snl→0` layered→bulk representation switch — is
essentially removed.

**Stage 2 — the whole-pack melt-out residual closed.** The last spike (θ=0.704) was probed and
turned out **not** to be a layer-count *difference* (the `snl` trajectory was identical across
θ±h) but the whole-pack *collapse itself*: melted layers empty in **mass** but retain **thickness**,
so reindexing a finite-`dz` layer when the pack drops 6→1 is a restructuring AD cannot
differentiate (FwdDiff 6.8 vs FD 52 at near-identical states). The fix: require **both** ice *and*
thickness `< SNOW_COMBINE_EMPTY_DZ` before collapsing a layer, so a melted-but-thick layer
**persists in place** and `snl` stays genuinely fixed — no reindex, ever. Result: across the whole
sweep `snl` is now **constant (`−6`)**, θ=0.704 drops **8.7e-1 → 1.5e-4**, the worst point over the
entire sweep is **4.5e-2** (isolated, not `snl`-related), the thermal solve **survives the persistent
empty-thick layers (no NaN)**, and the sweep's own verdict flips to **"gradient-BENIGN
(FwdDiff==FD)."** Default byte-identical (gates behind `_pc_smooth(T)`). The smooth-mode physics
differs slightly from the hard path (empty slots linger — the analog of the dimensional-`k` trade),
which is AD-neutral for calibration (AD and DDS use the same smooth path). Net: **the integer
layer-count — the last discrete residual of the forward-AD smoothness boundary — is resolved for the
melt window; the gradient is now benign across the full θ sweep, default-exact.**

**Stage 3: the θ=4.981 residual was not a kink. It was a ForwardDiff tie-break.** The isolated
4.5e-2 point was reproduced (FwdDiff 15.68 vs FD 16.42) and bracketed. It is **not** FD truncation
or a straddle: FD is **h-independent** (16.4250437 for every h from 1e-4 to 1e-6), and the bad
region is an *interval* (θ ∈ ~[4.978, 4.989]) where the AD slope is wrong while the loss itself is
smooth. Locating it took three probes in turn. A Dual-vs-FD field walk found nothing until step 14.
A **primal** walk (Dual value vs Float64 at the *same* θ) put the first divergence at step 14
inside `soil_temperature!`: top-soil `thk` was 0.943 on the Dual path and 1.398 on Float64, with
`t_h2osfc` 0.26 K off. A fork test then showed that every Dual-typed field was bit-identical to
Float64 going in, and that a Dual built from the Float64 state reproduced Float64 exactly. So the
only difference left was the **partials**.

The cause is ForwardDiff ≥ 1.0 semantics (1.4.3 here). An exact-value tie in `<`/`<=`/`>=`/`==` is
broken by comparing **partials**: `Dual(tfrz, ∂) >= tfrz` is `∂ >= 0`. CLM pins state **exactly**
at thresholds by assignment. A melting layer gets `t_soisno = tfrz`. After melt,
`h2osoi_liq = max(0, wmass0 − ice)` is exactly 0. The phase-change residual
`heatr = hm − hfus·(wice0 − ice)/dtime` is identically 0 whenever the melt is not mass-limited.
At those ties the value is exact but the partials are **roundoff noise**, so the noise sign chose
the branch. Two chains were found:
(i) `abs(heatr) > 0` passed on the noise, adding `fact·heatr` noise to `∂t_soisno` of the pinned top
soil layer. The next step's Kersten test `t_soisno >= tfrz` then took the **frozen** branch
(thk 0.94 instead of 1.40 W/m/K).
(ii) the snow test `h2osoi_liq·h2osoi_ice > 0` passed on a zero-valued, noise-partialled liquid
mass and reset the bottom snow layer to `tfrz` (Dual 273.150 K vs Float64 273.116 K).
Central FD runs in Float64, which never takes these branches, so the AD differentiated a function
the primal never evaluates.

**Fix: decide branches on the primal value.** `_primal(x)` (`smooth_ad.jl`) strips partials and is
the **identity on every `AbstractFloat`**. It now drives every branch test in the phase-change kernel
(identification, `hm`/`xm`/`heatr` guards, the `liq·ice` reset) and the three `t_soisno`-vs-`TFRZ`
tests in the soil/snow thermal-conductivity kernel. It is not a smoothing, so there is no `k` and no
physics bias. Under a Dual the branch is now exactly the Float64 branch. Float64 and Float32 are
byte-identical by construction, which is stronger than a `_pc_smooth(T)` gate.

| `CLM_SMOOTH=always` sweep, 28 points | before | after |
|---|---|---|
| θ = 4.981 | **4.54e-2** | **2.8e-11** |
| worst point over the sweep | 4.54e-2 (θ=4.981) | **6.2e-3** (θ=2.537) |
| points with rel ≤ 1e-8 | 0 / 28 | **21 / 28** |
| typical settled point | ~1e-3 | **~1e-10** |

The second-to-last row corrects §5b. The "~2–6e-4 / ~1e-3 FD-truncation-limited settled-melt
residual" was **mostly this same tie artifact, not truncation**. With branches on the primal,
settled points agree with central FD to ~1e-10. Seven points stay above 1e-8 (θ = 0.50, 0.91,
1.32, 2.33, 2.54, 2.74, 3.35; worst 6.2e-3). The worst one, θ=2.537, was checked directly and is
**not** a tie-break: a primal walk finds no Dual/Float64 divergence over all 480 steps. It is an FD
straddle of a genuine kink (§3/§4) sitting 1e-5 to 3e-5 from θ. The FD error is **non-monotone in h**
(6.0e-4, 2.0e-3, 6.2e-3, 2.1e-2 at h = 1e-3, 3e-4, 1e-4, 3e-5). It then collapses to **1.5e-10 at
h = 1e-5**, once the bracket no longer reaches the kink. Neighbours at θ ± 5e-4 agree to ~1e-11.
The AD derivative at θ is exact; the sweep's h = 1e-4 FD is the inexact side.

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

**A gate can also fail for a reason that is not a kink (§5c, Stage 3).** Before attributing an
AD-vs-FD gap to non-smoothness, check two things. First, FD must be **h-stable**. Second, the Dual
**primal** must equal the Float64 primal at the same θ. If FD is h-independent and the primals
diverge, the AD is differentiating a different branch from the one the model takes. In ForwardDiff
≥ 1.0 this happens at any **exact tie**: comparisons break value ties on partials, and a model that
pins state at thresholds by assignment hands them roundoff-noise partials. The θ=4.981 "residual"
was exactly this case. So was most of the "~1e-3 FD-truncation floor" that had been accepted as
irreducible. Branching on `_primal(x)` removes it without smoothing anything.

**Bottom line.** CLM.jl has **no forward-AD type wall**: the whole per-step physics is `Dual`-clean.
Its forward-AD boundary is **smoothness**, because the discrete snow events make the streamflow loss
piecewise-smooth. The continuous heat/mass partitions were already smoothed. Phase-3 then smoothed
the `imelt` melt label (§5b) and held the `snl` layer count fixed through melt (§5c, Stages 1–2).
Stage 3 made the phase-change and thermal-conductivity branches decide on the primal value, so a
Dual follows the Float64 branch at exact ties. The `CLM_SMOOTH=always` θ_snow sweep now has
**FwdDiff == central FD to ≤ 1e-8 at 21 of 28 points**, with a worst case of **6.2e-3**. That worst
case is at points where the primal is identical and the remaining gap is a genuine piecewise-smooth
kink inside the FD bracket. All of this is default byte-identical. `Dual`-clean is necessary but not
sufficient. Smoothness is the deeper requirement, **and so is branch-consistency**: the Dual must
take the same branch as the primal it is differentiating.
