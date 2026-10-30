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
| `gpu_ad_reverse_driver_validate` | FAIL | whole-driver CUDA reverse AD — known pending (STATUS.md §4); `canopy_rev_aux` scalar-iterates a device mask | not addressed here; see the follow-up section below |
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

## Follow-up: whole-driver CUDA reverse AD (`ad/cuda-reverse-driver`)

`gpu_ad_reverse_driver_validate` was the one failure left above. It reverse-
differentiates one timestep of the compositional driver chain (canopy energy balance
with 2 Newton iterations → soil temperature → surface hydrology → soil water → water
table → hydrology_no_drainage; 33 phases) on CUDA, then compares `dL/d(t_grnd)` with
the CPU reverse and with central finite differences of the CUDA forward chain,
L = Σ (soil-layer t_soisno)². Same server, Julia 1.12.6, Enzyme 0.13.216,
EnzymeCore 0.8.26, KernelAbstractions 0.9.45, CUDA.jl 5.11.3, GPUCompiler 1.15.2.

**Status: runs end to end on CUDA; gradient still wrong.** Full gate at 5102a2e
(`CLM_REVERSE_CPU=1`, FD on; 6.0 h wall, nearly all of it Enzyme/GPU compilation):

| | dL/d(t_grnd[1]) |
|---|---|
| CPU reverse (host phase list, unchanged) | -3.50231505e+00 |
| Central FD of the CUDA forward chain, h = 1e-2 / 5e-3 / 2.5e-3 | -3.50231505 / -3.50231504 / -3.50231505 (Richardson -3.50231505) |
| CUDA reverse (device phase list) | -1.34647340e-14 — **FAIL** |

CPU reverse and FD agree to < 1e-8 relative, so the CPU path and the device forward
chain are right; the device reverse loses the gradient. Per-phase CUDA-vs-CPU reverse
checks (`phase_diff`, same checkpoint, every shadow seeded with ones, CPU running the
same device phase list) agree to ≤ 1.3e-15 relative for tk_interface, lwrad, heatdiff,
gnet rural / rural_contrib / reduce, SNICAR, set_qflx, infiltration. The soil-temperature
**solve** and **post** phases disagree by O(1), in the surface-water (h2osfc) row: z,
tk_h2osfc, dz_h2osfc, c_h2osfc, frac_h2osfc, fact, dhsdT shadows. The batched band
solve was one cause (item 7 below); after its rule the same h2osfc-row mismatch remains,
so a second cause is still open. Not yet checked per phase: the canopy phases, cv,
tk_layers, setup, and the hydrology `*_dev` phases.

The soil-temperature block alone carries this gradient, and its device decomposition is
right on CPU: reversing only soil temperature with the same seed gives
dL/d(t_grnd[1]) = -3.5023150483 for the production `soil_temperature!` and
-3.5023150483 for the 13-phase decomposition on the CPU (shadow maxima of t_soisno,
t_h2osfc, frac_h2osfc, h2osfc identical to 7 digits). The open defect is CUDA-specific.

### Blockers found, in the order they surfaced

1. **`canopy_rev_aux` scalar-indexed device masks** (the failure in the table above).
   The August server WIP (`backup/server-wip-2026-08-06`, 775df61) already built that
   aux on the host and adapted it; it is applied on this branch (minus its CUDA 6
   compat bump) together with its allocation-free device decomposition of the soil-
   temperature, surface-hydrology, soil-water and water-table phases.
2. **`jtop`/`jbot` never set in the decomposed soil temperature.** The decomposition
   skipped `soil_temperature!`'s band-extent kernel, so `batched_band_solve!` indexed
   the -9999 scratch sentinel: `ERROR_ILLEGAL_ADDRESS` in the forward sweep. Fixed by a
   setup phase; a CPU test checks the decomposition against `soil_temperature!` exactly.
3. **Scalar kernel arguments under Enzyme on CUDA (upstream).** KernelAbstractions'
   Enzyme extension rejects any GPU kernel argument annotated Active ("Active kernel
   arguments not supported on GPU"), and Enzyme annotates every Float64 scalar kernel
   argument Active: one read from a global `Ref`, one read from a `Const` argument, and
   one passed through `EnzymeCore.ignore_derivatives` all fail; only a literal that LLVM
   constant-folds gets through. The same code differentiates exactly on the CPU
   backend. Reduced to `scripts/gpu_ad_active_scalar_repro.jl` (8-element axpy kernel).
   Every CLM kernel taking `dtime` or a physical constant hit it.
   **Workaround:** GPU launches go through `_kcall!` → `_gpu_launch!`, whose Enzyme rule
   calls KernelAbstractions' own kernel rule (`augmented_primal`/`reverse` on the
   `Kernel` object) with those scalars re-annotated `Const`. The kernel adjoint, tape
   and caching stay KernelAbstractions'. (Re-annotating through a nested `autodiff`
   call does not work: the nested call marks the same scalars Active again.)
   **Cost:** scalar kernel arguments get a zero derivative on GPU, so a quantity to be
   differentiated on device has to reach kernels as an array. Every scalar kernel
   argument on the driver reverse path is `dtime`, a physical constant or a
   configuration value. CPU launches never reach the rule.
4. **GPUArrays' scalar fill/broadcast kernels** (`bmatrix .= 0`, `fill!(bmatrix, 0)` in
   `set_matrix!`) hit the same error. Those two device scratch fills now use a CLM fill
   kernel (`_fill_dev!`; `fill!` unchanged on host arrays).
5. **Atomic scatter inside a kernel.** `compute_gnet_snicar!` scatters patch → column
   with Atomix atomics inside its kernel; Enzyme cannot lower the reverse of the CAS loop
   on CUDA (`LLVM ERROR: Cannot select ... AtomicLoad`, a hard process abort). The device
   SNICAR phase writes per-patch values to scratch and reduces them with
   `scatter_add_1d!`, whose launch-level rule is a race-free gather. The scatter rules
   now also match a `Duplicated` index array.
6. **Reverse-rule return type.** Enzyme checks the inferred return type of a reverse
   rule; the `_gpu_launch!` rule now builds it by dispatch so it is concrete.
7. **Wrong adjoint from the batched band solve (no error raised).** With 1–6 fixed the
   whole chain ran on CUDA but returned dL/d(t_grnd) = -1.3e-14 against -3.5023 on CPU
   and by FD. A per-phase CUDA-vs-CPU reverse comparison (same checkpoint, all shadows
   seeded with ones) matched every soil-temperature phase checked to ≤1.3e-15 relative
   except the solve phase, which was off by O(1). Enzyme's reverse of
   `_batched_band_solve_kernel!` (per-thread pivoted Gaussian elimination with a
   data-dependent system size, overwriting its scratch in place) is wrong on CUDA. A
   reverse rule for `batched_band_solve!` on GPU arrays now supplies the exact
   linear-solve adjoint (λ = A⁻ᵀ x̄ from the same batched solver on the transposed
   band; r̄ += λ; Ā -= λ xᵀ), FD-checked in `test/test_band_solve_adjoint.jl`. Not yet
   reduced to a standalone upstream reproduction.

The CPU reverse is unchanged: `driver_rev_phases` gives host filters the same phase
list as before (production `soil_temperature!`, `water_table!`, the canopy block, …)
and only device filters get the decomposition (`device` keyword; `*_dev_rev_phase!`).
`test/test_driver_reverse.jl` runs the device list on the host and checks its forward
state against the host list (t_soisno, t_grnd, t_veg, h2osoi_liq/vol, h2osfc, zwt).

Compile time is the practical cost: a fresh process spends hours compiling the 33
reverse thunks for CUDA (per-phase sweep: soil water 19 min, canopy energy 37 min,
urban ground-heat-flux kernel 100 min). The urban ground-heat-flux phase is skipped
when the domain has no urban column (an exact no-op there).

## Still open

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
# whole-driver CUDA reverse (follow-up section): gate, upstream repro, per-phase diagnostic
CLM_REVERSE_CPU=1 julia --project=scripts scripts/gpu_ad_reverse_driver_validate.jl
julia --project=scripts scripts/gpu_ad_active_scalar_repro.jl
PHASES=solve,post julia --project=scripts scripts/cuda_reverse_phase_diff.jl
```
