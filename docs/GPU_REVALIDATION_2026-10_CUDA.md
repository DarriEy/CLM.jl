# CUDA full-fleet re-validation — 2026-10

First full-fleet CUDA run since #263, on the same box: NVIDIA RTX 5070 Ti Laptop
(sm_120, 12 GB), driver 610.78 (CUDA 13.3), Windows 11, Julia 1.12.6, CUDA.jl 5.11.3,
KernelAbstractions 0.9.45. All device runs are **Float64** (`CuArray`, no `cu()`).

Every `scripts/gpu_validate_*.jl` and `scripts/gpu_ad_*.jl` harness (129) ran serially
at `main` `e5c3c18`, one fresh Julia process each, 25-min per-harness timeout (6.8 h in
all). The non-passes (except the known-pending reverse-AD one) were then rerun on this
branch with an 80-min limit, together with every harness that uses `scripts/gpu_adapt.jl` (see defect 3) — those had never actually run in Float64 on CUDA.

Site data: the small Bow / Aripuana `clm5_params.nc` + `surfdata_clm.nc` and the two
SNICAR optics files (6.7 MB total) were copied to the server and pointed at with
`SYMFLUENCE_DATA` / `CLM_PARAMFILE` / `CLM_FSNOWOPTICS` / `CLM_FSNOWAGING`, so no
harness had to be skipped for missing data.

## Result

| | `main` (e5c3c18) | this branch |
|---|---|---|
| harnesses | 129 | 129 (24 rerun, 105 unchanged from main) |
| pass | 120 | **128** |
| fail | 6 | 1 (known-pending) |
| timeout (25 min) | 3 | 0 |
| skipped | 0 | 0 |

Non-pass harnesses on `main` and their outcome here:

| harness | main | cause | here |
|---|---|---|---|
| `gpu_validate_clmdrv_e2e` | FAIL | `update_lnd2glc!` / `lnd2atm_rof!` scalar-index device arrays (defect 1) | PASS, all finite fields 0.0 |
| `gpu_validate_clmdrv_realinit_e2e` | FAIL | same | PASS, max rel 3.1e-14 over 427 entries |
| `gpu_validate_computewf2_e2e` | FAIL | hard-coded `import Metal` (missed by the #263 migration) | PASS, 0.0 |
| `gpu_validate_downscale_forcings` | FAIL | stale pco2 host reference (pre-5d27c43 per-column rescale) | PASS, 0.0 |
| `gpu_validate_soiltemp` | FAIL | arity: `PHASE_CHANGE_TEMP_K` arg added by #321 | PASS, 12/12, max 4.5e-13 |
| `gpu_ad_reverse_driver_validate` | FAIL | whole-driver CUDA reverse AD — known pending (STATUS.md §4); `canopy_rev_aux` scalar-iterates a device mask | not addressed |
| `gpu_validate_clmdrv_cn_e2e` | TIMEOUT | slow whole-driver CN run (and, on main, defect 1) | PASS, 0.0 over 754 entries (26 min); the #263 composition leak is gone |
| `gpu_validate_clmdrv_fates_e2e` | TIMEOUT | slow whole-driver FATES run | PASS, 7/7 fields, max rel 3.1e-16 (25 min) |
| `gpu_validate_cuda` | TIMEOUT | wrapper that runs `gpu_validate` + `gpu_ad_validate` + three clmdrv harnesses | PASS, 6/6 suites (72 min; run with an 80-min limit) |

`gpu_validate_ch4tran_e2e`, the open CH4-transport divergence in #263, now passes.

## Defects

### 1. `clm_drv!` could not run on any GPU (src)

The CESM coupling wiring in #328 added two host scalar loops to the default driver
path, after the Metal full-fleet run:

- `update_lnd2glc!` (`lnd2glc.jl`) — a per-column fill with a data-dependent
  duplicate-slot error check. Now runs on host copies of its (small) inputs when the
  state is on a device and copies the three `(nec+1, ng)` outputs back; the host path
  is the same loop on the same arrays.
- `lnd2atm_rof!` (`lnd2atm.jl`) — the stream-zero and dynbal-correction gridcell
  loops. Now broadcasts; bit-identical (a non-finite dynbal contributes exactly 0).

On CUDA this is `Scalar indexing is disallowed`; Metal fails the same way.

### 2. The CUDA "real GPU" test skipped itself on a CUDA box (test)

`test/test_backend.jl`'s fake-backend seam test registered a fake `:cuda` and then
`delete!`d it, which also removed the real `CLMCUDAExt` registration, so the following
"Real GPU execution" testset always took the `@test_skip` branch. Under `main` with
CUDA loaded: 28 pass / 3 broken. The seam test now restores the prior registration,
and the real-GPU testset runs `compute_forc_q!` on the device and checks it
bit-identical to the CPU result with a Float64 eltype. Under CUDA: **33 pass / 2 broken**
(AMDGPU, Metal); on a CPU-only box unchanged (32 pass / 3 broken).

### 3. `mf` demoted Float64 to Float32 on every backend (scripts)

`scripts/gpu_adapt.jl`'s `mf`/`mfs` always ran the `_F32` down-adaptor, so 23
harnesses printed "working precision: Float64" on CUDA while comparing a Float32
device run to the Float64 CPU reference — the same silent demotion #263 removed from
`cu()`. `mf` now down-converts only when the detected backend is Float32 (Metal);
`GPU_ADAPT_F32=1` restores the old path. All 23 were rerun in Float64.

### 4. Harness clone drift (scripts)

`computewf2` (Metal-only), `downscale_forcings` (stale reference), `soiltemp`
(arity) — see table.

## Still open

- **`gpu_ad_reverse_driver_validate`**: whole-driver CUDA reverse AD is pending
  (STATUS.md §4). On this server there is an uncommitted continuation of that work,
  now saved on the local branch `backup/server-wip-2026-08-06`.
- The CUDA.jl artifact runtime is 13.2 against a 13.3 driver; CUDA.jl logs
  `You are using CUDA 13.3.0, but CUDA.jl was precompiled for CUDA 13.2.0` at
  precompile. It does not affect results.
- The fleet runner glob `gpu_(validate|ad)_*` misses the core `scripts/gpu_validate.jl`;
  it is covered via `gpu_validate_cuda.jl` (PASS).

## Reproduce

```
julia --project=scripts -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=scripts scripts/gpu_validate_<name>.jl        # one harness
julia -e 'using CUDA, Test, CLM; include("test/test_backend.jl")'   # in an env with CLM+CUDA+Adapt+KA
```
