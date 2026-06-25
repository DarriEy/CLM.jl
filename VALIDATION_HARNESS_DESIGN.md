# CLM.jl — Full-Port Validation Oracle & Harness Design

Goal: a single, scalable validation system that drives **multi-year stability + Fortran
parity across every configuration path of the port**, and runs the matrix in the optimal
parallel shape. This document is the design; the build order is at the end.

---

## 0. The central problem — the *oracle*, not the runner

Running the model in every configuration is easy. Knowing whether each run is **correct**
is the hard part, and the answer is different per config:

- We have **Fortran reference dumps for exactly one domain (Bow) and a handful of configs**
  (biogeophys single-step, CN/BGC step + drift, LUNA). Generating a Fortran reference for
  an *arbitrary* config (crop on, lake landunit, Leung dust, …) means building+running
  Fortran CTSM with that namelist and dumping state — **minutes-to-hours of serial work per
  config**, and impossible for paths the Fortran build at hand doesn't enable.
- The configuration space is **combinatorial** (3 BGC modes × ~20 independent toggles ×
  9 domains × 7 landunits × 4 backends × parallelism modes). A Fortran oracle for the full
  cartesian product is infeasible and largely meaningless.

**Therefore the design is a layered oracle.** Fortran parity anchors the *core physics on
reference domains*; everything else is validated by oracles that need **no external truth**:
conservation laws and metamorphic/differential relations. This is what makes full-coverage
validation tractable.

---

## 1. The four oracle tiers

| Tier | Oracle | Needs Fortran? | Strength | Coverage role |
|------|--------|----------------|----------|---------------|
| **T1 — Parity** | Julia state == Fortran dump (per-field tol) | **Yes** (dumps) | Strongest — absolute correctness | Core physics on the configs we have/generate dumps for |
| **T2 — Conservation** | water/energy/C/N closure each step; multi-year finite + bounded drift | No | Strong — catches leaks, blow-ups, instability | **The workhorse** — runs for *any* config path |
| **T3 — Metamorphic / differential** | invariant relations between two runs | No | Strong for regressions | Scales to the whole flag lattice |
| **T4 — Observational** | streamflow KGE/NSE vs gauges | No (real obs) | Realism, not bit-correctness | Hydrologic plausibility across domains |

### T1 — Fortran parity (anchor)
Reuse `scripts/fortran_parity_common.jl` primitives (`build_bow_inst`, dump injection via
`read_fortran_restart!`, single-step run, per-field diff against `after_hydrologydrainage`).
Already covers: biogeophys (T_SOISNO/H2OSOI/T_GRND/ZWT/… ≤1e-10 rel where clean, banded tol
in `test_fortran_parity.jl`), CN/BGC step + 28-step drift, LUNA vcmax/jmax.
**Extend to:** an annual parity loop (inject Jan-1 spun-up IC → run a year → compare daily
aggregates + end-state to Fortran), and a *small, deliberately chosen* set of new Fortran
dumps for the highest-value config paths (see §4).

### T2 — Conservation (scalable workhorse)
Already in-model: `balance_check.jl` (water, errh2o ≤1e-5 mm/step) + `cn_balance_check.jl`
(C ≤1e-7, N ≤1e-3 gN/m²/step). `longhorizon_conservation.jl` does the multi-year free-run.
**Promote these to the universal pass/fail for every config**: for any config the harness can
construct, run N steps/years and assert (a) no NaN/Inf anywhere, (b) per-step balances within
threshold, (c) annual budget closes, (d) seasonal-cycle drift bounded under repeat-year
forcing. No Fortran needed → this is what gives *every* config a verdict.

### T3 — Metamorphic / differential (regression net over the flag lattice)
Pairs of runs whose relationship is known a-priori — no external truth:
- **Gated-off byte-identity:** every feature flag OFF ⇒ bit-identical to the baseline binary.
  (We already assert this per-PR; generalize to *all* flags systematically.)
- **Equivalent-method identity:** `use_matrixcn` == sequential cascade (≤1e-9, already proven);
  CENTURY vs MIMICS only where they should agree.
- **Invariance:** MPI 2-rank gather == serial (bit-identical, already a CI lane); threaded-clumps
  == serial; GPU(CPU-backend proxy / real device) == CPU.
- **AD vs FD:** ForwardDiff/Enzyme gradient == finite-difference (already in ad_robustness).
- **Restart round-trip:** write → read → continue == uninterrupted run (already in test_run_clm).
- **Symmetry/scaling:** ncopies=N replicate ⇒ N identical columns; spatial permutation invariance.

### T4 — Observational (realism)
`run_clm_streamflow.jl` across the 7 wired gauges (Bow/Krycklan/Abisko/Tagus/Massa/Baltimore/
Iceland): KGE/NSE bands per domain. Not parity — guards against physically-wrong-but-stable runs.

---

## 2. The configuration matrix — coverage, not cartesian

The flag lattice is too large for cartesian and most flags are **independent** (don't interact),
so full cartesian wastes effort. Use a **coverage strategy**:

**Axes** (from the config-surface map):
- **BGC mode:** SP, CN, FATES(sp), FATES(bgc) — *mutually exclusive* (pick 1).
- **CN sub-flags** (CN only, gated on use_cn): lch4, c13, c14, crop, crop_agsys, matrixcn, cndv,
  fun, flexiblecn, nitrif_denitrif; `cnfire_method∈{nofire,li2014,li2016,li2021,li2024}`;
  `decomp_method∈{CENTURY,MIMICS}`.
- **Driver flags** (largely independent): irrigate, hydrstress, luna, voc, ozone, aquifer_layer,
  soil_moisture_streams, lai_streams, threaded_clumps; dust∈{off,Zender2003,Leung2023};
  h2osfcflag∈{0,1}.
- **Landunit:** soil, crop, glacier(istice), deep-lake, wetland, urban{tbd,hd,md} — set by surfdata.
- **Domain:** Bow, Aripuanã, Stillwater, Krycklan, Abisko, Tagus, Massa, Baltimore, Iceland.
- **Backend:** cpu (always), metal, cuda, amdgpu (device-gated).
- **Parallelism:** serial, threaded-clumps, MPI-2rank, multi-GPU.
- **Init:** cold-start, Fortran-warm-restart.

**Coverage tactics (instead of 2^N):**
1. **Single-axis sweep** — baseline (SP/Bow/cpu/cold) with **exactly one flag flipped on**, for
   every flag. Catches "this flag alone breaks something." ~30 configs.
2. **t-wise (pairwise) combinations** — for the independent driver+CN flags, a pairwise-covering
   set (every pair of flag values co-occurs in ≥1 config). Covers all 2-way interactions in
   ~20–40 configs instead of thousands. (Generate the covering array offline; store as a table.)
3. **Stacked "realistic" profiles** — the flag *bundles* real science uses, e.g.
   `CN+crop+nitrif_denitrif+fun+flexiblecn+li2021-fire+lch4` (a full BGC production config),
   `SP+hydrstress+luna+ozone+voc` (full biophysics), FATES-bgc+spitfire. ~6–10 configs.
4. **Landunit coverage** — each landunit type exercised on ≥1 domain (glacier/snow/lake/urban
   smokes already exist; add crop + wetland).
5. **Domain coverage** — each of the 9 domains run at least at smoke depth; Bow + 2 contrasting
   (Aripuanã, Stillwater) at multi-year depth.
6. **Cross-cutting** — each backend and each parallelism mode run on ≥1 nontrivial config via
   the metamorphic invariants (== cpu/serial).

**Result:** ~120–160 distinct config runs give *single-flag + all-pairwise + realistic-bundle +
landunit + domain + backend* coverage — a tractable, defensible matrix instead of millions.

---

## 3. Harness architecture — one runner, declarative matrix

```
                    configs.toml  (the matrix: id → {mode, flags, domain, backend, depth, oracles})
                          │
                    ┌─────▼──────────────┐
                    │  validate.jl        │   one entry point
                    │  (config → verdict) │
                    └─────┬──────────────┘
        ┌─────────────────┼───────────────────┬──────────────────┐
        ▼                 ▼                   ▼                  ▼
   build_inst(cfg)   run(depth)          oracles[cfg]        record
   (domain+flags →   (N steps / N yrs,   T1 parity?          {id, pass/fail,
    CLMInstances)     repeat-year)        T2 conservation     per-field resid,
                                          T3 metamorphic      drift quantiles,
                                          T4 streamflow       wallclock}
                          │
                    results.jsonl  →  report (per-config verdict + tolerance histograms + trend)
```

- **Declarative matrix** (`test/validation/configs.toml`): each row is one config with its
  applicable oracle tiers and run depth (smoke=6 steps / day=48 / season / year / multi-year).
  New config = new row, not new code.
- **`validate.jl`**: maps a config row → builds the instance (reusing `build_bow_inst` /
  multisite / streamflow domain builders), runs to depth, applies the oracle set, emits a
  structured verdict to `results.jsonl`. Pure data out — no human-facing prints in the hot loop.
- **Reuse, don't rebuild:** every primitive exists (parity common, longhorizon, streamflow,
  balance checks, multisite builders, restart). The harness is the **orchestration layer** that
  unifies them under one matrix + one report — it should add ~no new physics.
- **Tolerances as data:** per-field/per-oracle tolerances live in the matrix, fitted from
  observed residual quantiles (e.g. 95th-percentile rel-diff over a clean baseline year) so the
  harness flags *regressions* (residual moves) rather than re-litigating known drift.

---

## 4. The Fortran-reference bottleneck (minimize it deliberately)

T1 is the only tier that costs Fortran runs. Strategy: **spend Fortran-dump generation only
where T2/T3 can't substitute** — i.e. where absolute physics correctness on a new path matters
and isn't implied by an invariant.

- **Already have:** Bow biogeophys (per-step), Bow CN/BGC (step+drift), LUNA. Keep as anchors.
- **Generate (one-time, batched, highest value):** a Bow **annual** reference (Jan-1 IC → 1 yr
  daily aggregates) for biogeophys + CN; and step references for **crop**, **lch4/methane**, and
  one **non-soil landunit** (lake or glacier) if the Fortran build enables them. ~5–8 new dump
  sets. Script it once as `scripts/gen_fortran_refs.sh` (serial, run on the box with the Fortran
  build) so it's reproducible and not a manual ritual.
- **Everything else** rides T2 (conservation) + T3 (metamorphic). A config with no Fortran dump
  is still *validated* — just by invariants + stability, not bit-parity. State this explicitly in
  the report so coverage is honest (parity-verified vs invariant-verified).

---

## 5. Optimal parallel execution plan

Config runs are **independent** ⇒ embarrassingly parallel. Three concurrency layers:

1. **Within a run:** the existing `use_threaded_clumps` / MPI / GPU paths (not needed for small
   validation domains; relevant only for the scaling-invariant configs).
2. **Across configs (the big win):** a process pool — `Distributed.jl`/`pmap` over the matrix
   rows, or a shell launcher (GNU-parallel-style) spawning `julia validate.jl <id>` per config.
   Process isolation is *required* anyway because parity/robustness tests are global-state
   sensitive (runtests.jl already runs them in subprocesses). Cap at `ncores-2`. ~150 configs ×
   (smoke ~30 s … year ~3–5 min) → **the whole smoke+single-axis matrix in minutes; the
   multi-year subset in ~an hour on a workstation.**
3. **Across machines / CI tiers:**
   - **Per-PR (fast, CI):** baseline byte-identity + single-axis smokes for touched subsystems +
     the MPI 2-rank lane. Minutes. (Mostly already what runtests does.)
   - **Nightly / on-demand (full):** the entire matrix at smoke+day depth + the metamorphic suite,
     parallel over cores. A few GPU/AMD configs only on a GPU runner (skip-guarded elsewhere).
   - **Weekly / release (deep):** multi-year stability + streamflow KGE across all domains +
     annual Fortran parity. Hours — the only genuinely long tier.

**Scheduling shape:** sort the matrix by cost (cheap smokes first → fail fast), run the pool,
stream verdicts to `results.jsonl`, render one report. The Fortran-ref *generation* is the only
inherently serial step and runs out-of-band (it's an input to T1, not part of the run loop).

**Driving it at scale:** for a one-shot full sweep, a `Workflow` (or a `pmap` script) fans the
matrix across workers, each returning a structured verdict — same pattern as the parallel
porting agents, but the "agents" are deterministic validation processes.

---

## 6. Coverage ledger (what "validated everything" means, honestly)

The report ends with a ledger per config: **parity-verified** (T1), **invariant-verified**
(T2+T3), **realism-checked** (T4), or **smoke-only** (ran finite, no oracle yet). The goal is
zero "smoke-only" rows: every config path reaches at least T2+T3, and the core/reference paths
reach T1. That's the concrete definition of "checked and validated every configuration path."

---

## 7. Build order

1. **Matrix + runner skeleton** — `configs.toml` + `validate.jl` that builds/runs/records for the
   SP/Bow/cpu/cold baseline and emits `results.jsonl` + a report. (Wires T2 conservation first —
   the universal verdict.)
2. **Single-axis sweep** — enumerate every flag flipped-on from baseline; T2+T3 (gated-off
   byte-identity) for each. Catches solo-flag breakage + proves default byte-identity systematically.
3. **Metamorphic suite (T3)** — fold in matrix==sequential, MPI==serial, restart round-trip,
   AD==FD, ncopies-replication as matrix rows.
4. **Domain + landunit coverage** — the 9 domains + 7 landunits at smoke/day depth (reuse
   multisite/streamflow/landunit builders).
5. **t-wise + realistic bundles** — generate the pairwise covering array; add the science bundles.
6. **Fortran T1 extension** — annual Bow parity loop + `gen_fortran_refs.sh` for the chosen new
   dumps (crop/lch4/landunit).
7. **Multi-year + streamflow deep tier** — wire T2 multi-year + T4 KGE across domains; schedule
   nightly/weekly. Add the parallel pool + CI tiering.

Each step is independently mergeable and adds coverage rows, never physics.
