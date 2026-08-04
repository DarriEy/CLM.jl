# CLM.jl — Port Status & Validation Scorecard

An honest, current account of how complete and how *validated* the Julia port of CLM5
(CTSM) is — across the four goals: **process fidelity, validation, differentiability (AD),
and GPU**. The aim is to separate "done and verified" from "done but unproven" from
"genuinely remaining," with no overclaiming.

_Last updated: 2026-08-02._

---

## TL;DR

A **functionally complete CLM5 port** with a working **layered-oracle validation harness**,
**Fortran-parity-anchored** and **invariant-validated across 8 contrasting domains**,
**forward- and reverse-differentiable** (whole-driver reverse incl. photosynthesis now
FD-validated on **both Julia 1.10 and 1.12**), and **Metal-GPU-validated at 0.0 parity**. The
honest asterisks: CUDA is now validated on real NVIDIA silicon, while AMDGPU remains
CPU-proxy-only;
the PHS-coupled (`use_hydrstress`) reverse path and whole-*function* canopy autodiff remain
1.12-blocked (the latter avoided by decomposition); and Fortran parity is anchored on Bow
reference configs by design (full-cartesian parity is out of scope).

---

## 1. Process fidelity (the port itself) — ✅ complete

- **Tiers A–F complete**: single-point CLM5 biogeophysics + hydrology, BGC (CN, FUN,
  FlexibleCN, matrix-CN, CNDV, methane, N-cycling, fire), crop, transient/dynamic subgrid,
  MPI/distributed, and FATES (SP + BGC + SPITFIRE).
- **~23,400-test suite green** on Julia 1.10 (min) and latest stable, every PR, plus an
  MPI 2-rank bit-identity lane.
- **Fortran single-step parity** on the Bow reference configs (biogeophysics, CN/BGC step +
  drift, LUNA) to ≤1e-10 relative where clean, banded tolerances elsewhere.

## 2. Validation harness — ✅ complete (this is the new backbone)

A **layered oracle** (no single config needs Fortran ground truth): T1 parity anchors core
physics; T2 conservation + T3 metamorphic give *every* config a verdict; T4 streamflow checks
realism. Built in 7 steps (PRs #123–#132), ~54 configs, run **subprocess-isolated** with
`--jobs` parallelism and `--tier` CI cadence (pr/nightly/weekly).

| Tier | Oracle | Status |
|------|--------|--------|
| **T1 parity** | Julia state == Fortran dump (16 fields, fitted band) | ✅ Bow anchor (n13461); new-config refs = instrumentation written + compile-verified (§5) |
| **T2 conservation** | finiteness + water/energy/C-N closure, smoke→multi-year | ✅ the universal verdict |
| **T3 metamorphic** | determinism · matrix==sequential · restart **bit-zero** · MPI==serial · AD==FD | ✅ all wired |
| **T4 streamflow** | in-harness KGE/NSE vs gauge | ✅ real, **multi-year across all 7 gauges** (best Abisko KGE=+0.31; modest skill = calibration, not wiring — §2.5) |

- **8 domains exercised** (not just wired): Bow + Aripuanã, Stillwater, Krycklan, Abisko,
  Tagus, Massa, Iceland, Baltimore — the last three **built from cloud ERA5 via Symfluence**
  this session, then run finite with zero harness changes.
- **Restart**: round-trip bit-exact; and a **complete checkpoint** (`CLM.write_checkpoint`)
  gives **bit-zero continue==uninterrupted from a file** (the curated NetCDF restart omits
  forward-feeding state; the checkpoint serializes all of it).
- Harness building **found and fixed a real CNDV cold-start bug** and a real harness-design
  flaw (in-process global-state leakage → false verdicts) — it earns its keep.

### 2.5 Multi-year T4 streamflow — real uncalibrated baselines

All 7 wired gauges run to completion at their **full multi-year registry windows**
(multi-year spin-up + multi-year scoring), in ~3 min wall under 7-way parallelism; every
domain scored ≥1,000 overlapping gauge-days with all simulated runoff finite (nbad=0).
These are **uncalibrated cold-start** baselines — the honest floor before any calibration.

| Domain (climate) | Scored window | KGE | NSE | r | α | β | Note |
|---|---|---:|---:|---:|---:|---:|---|
| Abisko (arctic) | 2010–18 | **+0.31** | −0.65 | 0.48 | 1.41 | 1.20 | best skill; only positive KGE |
| Baltimore (urban) | 2010–18 | −0.18 | **+0.12** | 0.58 | 0.18 | 0.26 | positive NSE; under-predicts volume 4× |
| Bow (alpine) | 2003–04 | −0.32 | −3.4 | 0.18 | 2.03 | 0.97 | volume right, snowmelt timing/variance off |
| Krycklan (boreal) | 2010–18 | −0.59 | −5.0 | 0.05 | 2.27 | 1.09 | same: right mean, 2× too variable |
| Massa (glaciated) | 2010–18 | −0.64 | −0.5 | 0.02 | 0.09 | 0.06 | glacial melt discharge not represented |
| Iceland Jökulsá (glacial) | 2016–18 | −0.69 | −3.6 | 0.05 | 0.02 | 0.01 | glacial + gauge-area mismatch |
| Tagus (Mediterranean) | 2010–18 | −36 | −1417 | 0.08 | 35.3 | 15.8 | basin-area/gauge mismatch (16× volume) — config, not physics |

**Honest read:** volume is often right (β≈1 for Bow/Krycklan/Abisko); timing/variance is the
uncalibrated error (α≈2 — classic snowmelt mistiming). Tagus/Iceland/Massa are **domain-config
issues** (auto-delineated gauge-area mismatch; no glacial-melt routing), not physics faults —
they need explicit per-gauge `area_km2` overrides in the `DOMAINS` registry before their KGE is
meaningful. Real positive KGE requires a **calibration pass** (DDS over baseflow/snow/soil-
hydrology params + routing lag + area corrections) — that belongs to the separate Symfluence
calibration pipeline and was deliberately **not** run; these are the baseline-skill numbers.

## 3. Differentiability (AD)

- **Forward-mode (ForwardDiff)** — ✅ broad; validated vs finite-difference across kernels
  (`test_ad_robustness.jl`, `test_ad_e2e.jl`), including on the Metal device path.
- **Reverse-mode (Enzyme)** — ✅ **whole-driver energy-balance reverse now validates on BOTH
  Julia 1.10 and 1.12.** The top-level `clm_drv_reverse!` (convergence-aware, `use_psn=0`:
  canopy(decomposed) → soil_temperature! → surface-hydrology → soil_water! → water_table! →
  hydrology_no_drainage) gives `dL/d(t_grnd)` matching finite-difference — **1.12 rel 4.4e-7,
  1.10 rel 2.0e-7, gradient value bit-identical across versions**. The 1.12 blocker was found
  and fixed: scratch-allocation closures capturing a `Type`-valued `FT` (`similar(arr, FT, …)`)
  emitted a `has_free_typevars` runtime check Enzyme 1.12 can't differentiate; dropping the
  captured `FT` (`similar(arr, …)` + `zero(eltype(arr))`) is byte-identical and GPU-safe
  (`soil_water_movement.jl`, `surface_water.jl`). Decomposed canopy, multistep/checkpointed
  composition, and the hydrology phases all pass on 1.12.
  - **Photosynthesis-coupled (`use_psn=1`) whole-driver reverse now ALSO validates on 1.12**
    (rel 4.443e-7 vs FD; 1.10 unchanged). The 4 photosynthesis KA kernels segfaulted under
    Enzyme reverse on 1.12; fixed with an **AD-mode host-loop fallback** — each kernel body
    factored into a shared `@inline` function the KA `@kernel` still inlines (GPU path
    unchanged), with a host `for`-loop branch taken only under AD (a `Ref{Bool}` flipped by
    the reverse engine in `try/finally`). **Primal is byte-identical** (`max|KA−hostloop| = 0.0`).
  - **Remaining 1.12 gaps:** (1) **PHS-coupled** photosynthesis (`use_hydrstress=true`, its own
    `_psn_phs_*` KA kernels) would still segfault under reverse on 1.12 — same fix extends it;
    (2) the **whole-function** (non-decomposed) canopy autodiff hits the historical "sret wall"
    segfault — already avoided by the decomposed production path. 1.10 retains everything.
- **AD-smoothing** of discontinuities (btran, phase-change, snow merge/split) — ✅ done, and
  byte-identical when smoothing is off.

## 4. GPU

- **Apple Metal** — ✅ **validated at 0.0 parity** across the entire biogeophys/hydrology
  driver + BGC (decomposition, MIMICS, methane, fire, N-cycling, phenology, allocation).
  This is real hardware, not a proxy.
- **CUDA (NVIDIA)** — ✅ **validated on real silicon** (2026-08-02, NVIDIA GeForce RTX
  5070 Ti Laptop GPU, compute capability 12.0, Julia 1.12.6, CUDA.jl 5.11.3). Backend
  registration/device transfer passed; core kernels passed at ≤1.11e-16 absolute CPU
  difference; forward AD passed 15/15; the whole biogeophysics/hydrology `clm_drv!`
  timestep matched over 410 finite outputs at ≤4.701e-08 relative; and the `use_cn=true`
  composite matched over 754 finite outputs at ≤2.23e-08 relative. The synthetic driver
  fixture still contains non-finite uninitialized fields, so those entries are explicitly
  skipped rather than counted as parity. A portable real-NetCDF initialization gate now also
  passes four timesteps over mixed soil/lake columns: 36 fields and 427 finite values compared,
  with maximum relative difference 1.862e-03 (integration tolerance 1e-02).
  `scripts/gpu_validate_cuda.jl` includes all three whole-driver gates.
- **AMDGPU (ROCm)** — ⚠️ CPU-proxy-validated only; not yet run on real AMD silicon.
- **Reverse-AD on GPU** — ✅ the CUDA kernel-shape suite passes **4/4**. Atomic scatter is
  differentiated at the host launch boundary: `scatter_add_1d!` keeps the portable Atomix
  primal, while its reverse rule launches a race-free gather kernel instead of asking Enzyme
  to differentiate Atomix's CAS loop. Forward-vs-reverse and CUDA-vs-CPU errors are both
  exactly 0 for scatter; the other probes remain at ≤2.22e-16 cross-AD and ≤1.11e-16 device
  parity. Whole-driver CUDA reverse is still pending migration of embedded per-thread scatter
  sites to launch-level operations and an end-to-end gradient gate. Enzyme has no Metal device
  support; Metal's primal Atomix path is unchanged.

## 5. Fortran parity breadth

T1 is the only tier needing external truth, and it's **anchored, not exhaustive** — by design
(a Fortran build+run per config is infeasible and largely meaningless). The Bow anchor is
validated; broadening to new configs (e.g. `use_cn`) needs the **pdump instrumentation** that
had been stripped from the build. That instrumentation is **reconstructed and compile-verified**
(`scripts/validation/fortran_pdump/`: `pdumpMod.F90` + driver-hooks diff + recipe) — it builds
into CTSM cleanly. The final **byte-verify run did not complete**: the CTSM install tree
vanished from the filesystem mid-run (an environment/drive issue, not an instrumentation fault).
Re-run the README recipe once the CTSM case is restored.

---

## The honest scorecard

| Goal | State | The asterisk |
|------|-------|--------------|
| Process fidelity | ✅ complete | parity anchored on reference configs, invariant-validated elsewhere (by design) |
| Validation harness | ✅ complete | T4 now multi-year across all 7 gauges (§2.5); skill is modest by design (calibration, not wiring) |
| Forward AD | ✅ validated | — |
| Reverse AD | ✅ whole-driver incl. **photosynthesis** (`use_psn=1`) on 1.10 **and 1.12** | PHS-coupled (`use_hydrstress`) reverse + whole-*function* canopy still 1.12-blocked (latter avoided by decomposition) |
| GPU Metal | ✅ validated 0.0 | — |
| GPU CUDA | ✅ real-hardware validated | synthetic + CN + portable real-NetCDF mixed soil/lake full-driver gates pass |
| GPU AMD | ⚠️ proxy only | unrun on real AMD hardware |

## Concrete remaining work to claim "fully validated AD + GPU"

1. **Reverse-AD on 1.12** — energy-balance AND photosynthesis-coupled (`use_psn=1`) whole-driver
   reverse are both FD-validated on 1.12 (done). Only the PHS-coupled path (`use_hydrstress`, its
   own `_psn_phs_*` KA kernels) remains — the same AD-mode host-loop transform extends to it.
2. ~~Run `scripts/gpu_validate_cuda.jl` on real NVIDIA hardware~~ — **done** (§4). Remaining
   GPU breadth is real AMD hardware, multi-GPU/MPI execution, and CUDA reverse-mode AD.
3. **Restore the CTSM case + run the pdump verify** → broaden T1 parity beyond the Bow anchor.
4. ~~Multi-year T4 streamflow at depth~~ — **done** (§2.5: all 7 gauges, full multi-year
   windows, real uncalibrated KGE/NSE). Remaining there is *calibration* for real KGE skill
   (DDS/optimization over hydrology params + per-gauge `area_km2` overrides) — a separate
   Symfluence-pipeline effort, not a harness gap.

Everything except the remaining hardware-dependent GPU breadth and (3)'s external CTSM case
dependency is code that exists and runs today.
