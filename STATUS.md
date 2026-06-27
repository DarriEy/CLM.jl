# CLM.jl — Port Status & Validation Scorecard

An honest, current account of how complete and how *validated* the Julia port of CLM5
(CTSM) is — across the four goals: **process fidelity, validation, differentiability (AD),
and GPU**. The aim is to separate "done and verified" from "done but unproven" from
"genuinely remaining," with no overclaiming.

_Last updated: 2026-06-27._

---

## TL;DR

A **functionally complete CLM5 port** with a working **layered-oracle validation harness**,
**Fortran-parity-anchored** and **invariant-validated across 8 contrasting domains**,
**forward- and reverse-differentiable** (whole-driver energy-balance reverse now FD-validated
on **both Julia 1.10 and 1.12**), and **Metal-GPU-validated at 0.0 parity**. The honest
asterisks: the *photosynthesis-coupled* reverse path is still 1.10-only on 1.12, CUDA/AMD GPU
is validated by CPU proxy but **not yet run on real silicon**, and Fortran parity is anchored
on Bow reference configs by design (full-cartesian parity is out of scope).

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
| **T4 streamflow** | in-harness KGE/NSE vs gauge | ✅ real (Bow KGE=−0.282; modest skill = calibration, not wiring) |

- **8 domains exercised** (not just wired): Bow + Aripuanã, Stillwater, Krycklan, Abisko,
  Tagus, Massa, Iceland, Baltimore — the last three **built from cloud ERA5 via Symfluence**
  this session, then run finite with zero harness changes.
- **Restart**: round-trip bit-exact; and a **complete checkpoint** (`CLM.write_checkpoint`)
  gives **bit-zero continue==uninterrupted from a file** (the curated NetCDF restart omits
  forward-feeding state; the checkpoint serializes all of it).
- Harness building **found and fixed a real CNDV cold-start bug** and a real harness-design
  flaw (in-process global-state leakage → false verdicts) — it earns its keep.

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
  - **Remaining 1.12 gaps:** (1) the **photosynthesis-coupled** whole-driver (`use_psn=1`)
    segfaults in the KA-kernel `photosynth_ci_solve!` under Enzyme reverse on 1.12 (works on
    1.10) — closeable with an AD-mode host-loop fallback for that one ci-solve; (2) the
    **whole-function** (non-decomposed) canopy autodiff still hits the historical "sret wall"
    segfault — already avoided by the decomposed production path. 1.10 retains everything
    including the photosynthesis-coupled path.
- **AD-smoothing** of discontinuities (btran, phase-change, snow merge/split) — ✅ done, and
  byte-identical when smoothing is off.

## 4. GPU

- **Apple Metal** — ✅ **validated at 0.0 parity** across the entire biogeophys/hydrology
  driver + BGC (decomposition, MIMICS, methane, fire, N-cycling, phenology, allocation).
  This is real hardware, not a proxy.
- **CUDA (NVIDIA) / AMDGPU (ROCm)** — ⚠️ **CPU-proxy-validated only; never run on real
  silicon.** The extension backends + `clm_set_backend(:cuda|:amdgpu)` registry exist and
  pass the CPU-proxy path. A one-command validator `scripts/gpu_validate_cuda.jl` is ready;
  the remaining step is to **run it on an actual NVIDIA GPU** (rented cloud box via SSH, or
  the JuliaGPU Buildkite `juliagpu` queue). Until then, "CUDA-optimized" is unproven.
- **Reverse-AD on GPU** — deferred (Enzyme has no Metal device support; CUDA path untried).

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
| Validation harness | ✅ complete | T4 skill is modest (calibration, not wiring); deep multi-year tier opt-in |
| Forward AD | ✅ validated | — |
| Reverse AD | ✅ whole-driver (energy-balance) on 1.10 **and 1.12** | photosynthesis-coupled (`use_psn=1`) reverse still 1.10-only on 1.12 (KA-kernel ci-solve) |
| GPU Metal | ✅ validated 0.0 | — |
| GPU CUDA/AMD | ⚠️ proxy only | **unrun on real hardware** — `gpu_validate_cuda.jl` ready for a GPU box |

## Concrete remaining work to claim "fully validated AD + GPU"

1. **Photosynthesis-coupled reverse-AD on 1.12** — the energy-balance whole-driver reverse is
   now FD-validated on 1.12 (done); the `use_psn=1` path needs an AD-mode host-loop fallback for
   the KA-kernel `photosynth_ci_solve!` to clear its Enzyme-on-1.12 segfault.
2. **Run `scripts/gpu_validate_cuda.jl` on real NVIDIA hardware** (rented GPU over SSH, or
   Buildkite `juliagpu`) — the one external dependency for "CUDA-validated".
3. **Restore the CTSM case + run the pdump verify** → broaden T1 parity beyond the Bow anchor.
4. (optional) Multi-year T4 streamflow at depth + calibration for real KGE skill.

Everything except (2) and (3)'s external dependencies is code that exists and runs today.
