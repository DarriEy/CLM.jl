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
table → hydrology_no_drainage; 36 device phases on this domain) on CUDA, then compares
`dL/d(t_grnd)` with the CPU reverse and with central finite differences of the CUDA
forward chain, L = Σ (soil-layer t_soisno)². Same server, Julia 1.12.6, Enzyme
0.13.216, EnzymeCore 0.8.26, KernelAbstractions 0.9.45, CUDA.jl 5.11.3, GPUCompiler
1.15.2.

**Result: PASS** (`gate2`, commit 665a90b, `CLM_REVERSE_CPU=1 CLM_REVERSE_FD=1`,
16466 s wall including compilation):

| | dL/d(t_grnd[1]) |
|---|---|
| CUDA reverse (36-phase device list) | -3.50231505e+00 |
| CPU reverse (15-phase host list) | -3.50231505e+00 (absdiff 2.8e-14, tol 7.0e-9) |
| CUDA forward FD, h = 1e-2 / 5e-3 / 2.5e-3 | -3.50231505 / -3.50231504 / -3.50231505 |
| Richardson FD | -3.50231505e+00 (rel 4.8e-10 vs CUDA reverse, bracketed) |

The CPU reverse runs the host phase list (canopy block, `soil_temperature!`, 10
surface-hydrology phases, soil water, water table, hydrology_no_drainage = 15 phases);
the CUDA reverse runs the device decomposition of the same chain (36 phases).

### Blockers, in the order they surfaced

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
   rule; the launch rule builds it by dispatch so it is concrete.
7. **Stale primal inputs in the GPU kernel reverse (silent wrong gradient).** With 1–6
   fixed the whole chain ran on CUDA but returned dL/d(t_grnd) = -1.3e-14 against
   -3.5023 on CPU and by FD. Localised by
   `scripts/cuda_reverse_phase_diff.jl` (same checkpoint, every shadow seeded with ones,
   CUDA vs the same device phase on CPU): every phase matched to ≤1.3e-15 except the
   soil-temperature solve and post phases; splitting the solve into RHS / matrix / init
   / band solve / extract made every piece match. KernelAbstractions' GPU reverse
   kernel re-reads its primal inputs from the arrays it is handed at reverse time; when
   a later launch in the same differentiated call has overwritten one of them (the
   extraction overwrites the `t_soisno` the RHS kernels read; the phase-change kernels
   overwrite inputs of `_fn1_kernel!`), the coefficients' adjoints come out wrong while
   the linear-path adjoints stay right. The CPU path is not affected (Enzyme's own
   overwritten analysis caches what it needs).

   **Workaround (3 + 7):** GPU launches go through `_kcall!` → `_gpu_launch!`, whose
   Enzyme rule (a) snapshots the launch's arguments (one `deepcopy`, aliasing kept) and
   (b) calls KernelAbstractions' own kernel rule (`augmented_primal`/`reverse` on the
   `Kernel` object) with scalar arguments re-annotated `Const`, reversing against the
   snapshot. The kernel adjoint and tape stay KernelAbstractions'. Re-annotating through
   a nested `autodiff` call does not work: the nested call marks the scalars Active
   again.
   **Cost:** one device copy of each launch's arguments, held until the reverse sweep
   reaches it (per phase — each compositional phase is a separate Enzyme call); and
   scalar kernel arguments get a zero derivative on GPU, so a quantity to be
   differentiated on device has to reach kernels as an array. Every scalar kernel
   argument on the driver reverse path is `dtime`, a physical constant or a
   configuration value. CPU launches never reach the rule.
   *Correction (follow-up 2 below):* as merged, the snapshot broke runtime activity's
   `shadow === primal` marker for Const members of kernel-argument structs. Canopy
   gradients on CUDA came out wrong and Const inputs were corrupted; the shadows are now
   remapped onto the snapshot.

An intermediate suspicion that Enzyme's reverse of the batched band-solve kernel was
itself wrong on CUDA did not hold: `scripts/cuda_band_solve_reverse_check.jl` (raw
kernel launch vs the exact linear-solve adjoint, which is FD-checked in
`test/test_band_solve_adjoint.jl`) gives relative errors 2.1e-16 / 3.2e-16 on CUDA,
the same as CPU. No rule is used for it.

The CPU reverse is unchanged: `driver_rev_phases` gives host filters the same phase
list as before (production `soil_temperature!`, `water_table!`, the canopy block, …)
and only device filters get the decomposition (`device` keyword; `*_dev_rev_phase!`).
`test/test_driver_reverse.jl` runs the device list on the host and checks its forward
state against the host list, and a soil-temperature-only reverse gives
dL/d(t_grnd[1]) = -3.5023150483 for production `soil_temperature!`, for the device
decomposition on CPU, and for the device decomposition on CUDA.

Compile time is the practical cost: a fresh process spends hours compiling the reverse
thunks for CUDA (canopy energy 37 min, soil water 19 min, ground-heat-flux reduction
20 min). The urban ground-heat-flux phase (a 93 min compile) is left out of the phase
list when the domain has no urban column, where it is an exact no-op.

## Follow-up 2: per-phase checks and compile time (`ad/cuda-reverse-phase-checks`)

Same server and stack. All runs from one checkpoint on this domain (3 patches, 2 of
them in the exposed-vegetation filter).

### Per-phase checks found a CUDA-only bug in the #373 launch rule

`scripts/cuda_reverse_phase_diff.jl PHASES=all` reverses every one of the 36 device
phases from the same checkpoint on CUDA and on CPU (same device phase list forced onto
the host). Each Float array leaf of the phase's view gets an independent uniform(0.5, 1.5)
seed per entry, so every output carries a cotangent and every input gets a gradient. The
loss is not one scalar along one path. The diff then compares every gradient leaf.
`DETAIL=1` also reruns the CUDA reverse on the same input and arbitrates the worst
entries with a central FD of the seeded loss on the CPU forward.

Before the fix, the canopy resist and energy phases were wrong on CUDA, with scaled
differences of 5.5 and 7.6e3 on the first iteration and up to 3e16 on the second. They
were also not reproducible: reruns on identical input differed by up to 1.6e11. FD sided
with the CPU to 8–10 digits, for example `esai_patch[2]`: CPU −16.808422889, CUDA
0.853358, FD −16.808422855. Patch 1, which is outside the filter, was always right. The
soil-temperature and hydrology phases matched when run on their own (≤4e-15,
deterministic). In the full sequence they showed 1e-6 to 1e-2 differences, inherited
from inputs that the canopy reverse had corrupted. The forward states of the two
backends agree to ≤4e-14 after every phase (`scripts/cuda_reverse_forward_diff.jl`),
and no array aliasing differs between the host and the adapted device bundle.

**Cause.** Under runtime activity, a kernel-argument struct built inside the
differentiated code (`CfResIn`, `CfEnC`, the gnet `pin`/`cin` NamedTuples, …) can have
members that come from a Const argument: the canopy aux's forcing arrays, `dleaf_pft`,
`forc_pbot`. For those members Enzyme's shadow *is* the primal array, and
`shadow === primal` is how runtime activity marks them inactive. The #373 launch rule
reverses KernelAbstractions' kernel against `Duplicated(snapshot, live_shadow)`. The
primal is now the snapshot copy while the shadow is still the original array, so the
kernel reverse treats those members as active and accumulates their adjoints into the
original primal. That corrupts the Const aux and the live forcing arrays: the canopy
gradients come out wrong, and inputs read by later phases are corrupted. The CPU path
never reaches the rule.

**Fix** (`src/infrastructure/enzyme_rules.jl`). The augmented primal keeps the
deepcopy's `IdDict`. The reverse maps every shadow leaf that is one of the snapshotted
primal objects to its copy (`_remap_shadow`), which restores `shadow === primal`. A host
unit test covers the helper (`test/test_driver_reverse.jl`).

After the fix (`pdfix1`), all 36 phases agree and every CUDA rerun is bit-identical:

| # | phase | max abs Δ | max rel Δ (per leaf) | before fix: max scaled Δ | CUDA 1st call, lattice 6/512 (s) | CUDA 1st call, lattice 3/64 (s) |
|---|---|---|---|---|---|---|
| 1 | `cf_rev_init!` | 5.684e-14 | 2.179e-16 | 2.179e-16 | 1438 | 417 |
| 2 | `cf_rev_friction!` | 8.882e-16 | 1.431e-16 | 1.431e-16 | 437 | 108 |
| 3 | `cf_rev_resist!` | 1.455e-11 | 4.614e-16 | 5.497e+00 | 1103 | 83 |
| 4 | `cf_rev_energy!` | 5.722e-06 | 2.623e-14 | 7.633e+03 | 2318 | 110 |
| 5 | `cf_rev_friction!` | 1.776e-15 | 6.640e-16 | 6.403e-01 | 0 | 0 |
| 6 | `cf_rev_resist!` | 5.821e-11 | 3.908e-16 | 9.734e+06 | 0 | 0 |
| 7 | `cf_rev_energy!` | 5.722e-06 | 6.964e-14 | 3.271e+16 | 0 | 0 |
| 8 | `soiltemp_setup_rev_phase!` | 0 | 0 | 0 | 139 | 89 |
| 9 | `soiltemp_tk_layers_array_phase!` | 7.105e-15 | 5.137e-16 | 5.137e-16 | 313 | 98 |
| 10 | `soiltemp_tk_interface_array_phase!` | 1.776e-14 | 4.561e-16 | 4.561e-16 | 99 | 92 |
| 11 | `soiltemp_cv_array_phase!` | 4.657e-10 | 1.593e-16 | 1.593e-16 | 299 | 134 |
| 12 | `soiltemp_lwrad_array_phase!` | 0 | 0 | 0 | 122 | 68 |
| 13 | `soiltemp_gnet_patch_rural_array_phase!` | 0 | 0 | 1.158e-02 | 296 | 70 |
| 14 | `soiltemp_gnet_patch_rural_contrib_array_phase!` | 2.910e-11 | 2.041e-16 | 1.122e-02 | 395 | 67 |
| 15 | `soiltemp_gnet_patch_reduce_array_phase!` | 0 | 0 | 0 | 73 | 58 |
| 16 | `soiltemp_gnet_snicar_rev_phase!` | 4.657e-10 | 2.088e-16 | 1.140e-02 | 383 | 101 |
| 17 | `soiltemp_heatdiff_rev_phase!` | 1.137e-13 | 3.623e-16 | 3.623e-16 | 171 | 94 |
| 18 | `soiltemp_rhs_rev_phase!` | 7.812e-03 | 1.975e-16 | 4.264e-05 | 241 | 193 |
| 19 | `soiltemp_matrix_rev_phase!` | 1.600e+01 | 2.969e-15 | 2.969e-15 | 221 | 192 |
| 20 | `soiltemp_tvinit_rev_phase!` | 0 | 0 | 0 | 81 | 58 |
| 21 | `soiltemp_bandsolve_rev_phase!` | 2.274e-13 | 5.347e-16 | 2.290e-05 | 92 | 68 |
| 22 | `soiltemp_tvextract_rev_phase!` | 0 | 0 | 0 | 86 | 66 |
| 23 | `soiltemp_post_rev_phase!` | 7.913e-11 | 1.811e-14 | 6.536e-04 | 994 | 394 |
| 24 | `setsoilfrac_rev_phase!` | 8.438e-15 | 2.339e-15 | 5.198e-05 | 67 | 59 |
| 25 | `setfloodc_rev_phase!` | 0 | 0 | 0 | 66 | 68 |
| 26 | `satexcess_dev_rev_phase!` | 0 | 0 | 0 | 110 | 74 |
| 27 | `setqflx_rev_phase!` | 0 | 0 | 0 | 68 | 63 |
| 28 | `inflexcess_dev_rev_phase!` | 2.220e-16 | 1.492e-16 | 1.069e-06 | 89 | 82 |
| 29 | `routeinfl_rev_phase!` | 0 | 0 | 0 | 68 | 58 |
| 30 | `updateh2osfc_dev_rev_phase!` | 1.998e-13 | 7.125e-14 | 7.125e-14 | 213 | 203 |
| 31 | `infil_rev_phase!` | 0 | 0 | 0 | 62 | 57 |
| 32 | `totalrunoff_dev_rev_phase!` | 0 | 0 | 0 | 129 | 124 |
| 33 | `rootsink_dev_rev_phase!` | 0 | 0 | 0 | 163 | 138 |
| 34 | `soilwater_dev_rev_phase!` | 2.048e+03 | 3.849e-15 | 2.100e-03 | 1120 | 102 |
| 35 | `watertable_dev_rev_phase!` | 6.661e-16 | 3.975e-16 | 5.850e-06 | 250 | 192 |
| 36 | `hydnodrain_rev_phase!` | 7.451e-09 | 3.037e-16 | 2.420e-05 | 439 | 398 |
| | **total** | | | | **12145** | **4178** |

"max rel Δ (per leaf)" is max|Δ| / max|g_cpu| for the worst leaf. The large absolute
values (2e3 in soil water, 16 in the matrix phase) are on leaves whose gradients are
themselves very large (up to ~1e18 in soil water). "before fix" is the same diff at the pre-fix commit (scaled by
max(max|g_cpu|, 1)). Phases 5–7 reuse the thunks compiled for 2–4.

### Gradients with respect to other inputs (whole chain)

`CLM_REVERSE_EXTRA=1` in the gate adds three checks.

- **(a)** For the default loss, every one of the 3484 `db.inst` input-gradient leaves is
  compared, CUDA (device list) vs CPU (host list).
- **(b)** A second loss puts random weights on soil/snow temperature, t_grnd, t_veg, soil
  liquid and ice, h2osfc, sensible, latent and ground heat flux, surface runoff and
  infiltration. For that loss every leaf is compared again.
- **(c)** Central FD of that loss on the CUDA forward with respect to t_grnd, a soil-layer
  t_soisno, t_veg, h2osoi_liq, `watsat` and `bsw` (soil parameters) and
  `forc_lwrad_downscaled`. The canopy aux is held fixed: the reverse treats it as Const,
  so a forcing read through the aux gets no adjoint, and the FD has to differentiate that
  same function.

**Gate result (`gate4`, 9bbf2b2, `CLM_REVERSE_CPU=1 CLM_REVERSE_FD=1 CLM_REVERSE_EXTRA=1`,
7190 s wall): PASS.** The reference value moved from −3.5023 to −3.4433 because current
main reads `e_ice` from the params file (#376) and the litter-mixing params (#378); it is
the value the soil-temperature-only check gives on cpu-host, cpu-device and cuda-device.

| check | result |
|---|---|
| dL/d(t_grnd[1]), CUDA vs CPU | −3.44327925 both, absdiff 2.1e-14 |
| dL/d(t_grnd[1]), Richardson FD of the CUDA forward | rel 1.4e-9, bracketed |
| (a) default loss, all 3484 db.inst leaves, CUDA vs CPU | worst leaf rel 2.9e-14 (219 nonzero entries) |
| (b) wide loss, all 3484 db.inst leaves, CUDA vs CPU | worst leaf rel 1.1e-14 (478 nonzero entries) |
| (c) FD wide loss / t_grnd | rel 4.5e-10 |
| (c) FD / t_soisno[c0, 13] (276.3 K) | rel 3.5e-10 |
| (c) FD / t_veg[2] | rel 4.1e-10 |
| (c) FD / h2osoi_liq[c0, 14] | rel 7.9e-9 |
| (c) FD / watsat[c0, 2] | rel 9.9e-11 |
| (c) FD / bsw[c0, 3] | rel 1.7e-8 |
| (c) FD / forc_lwrad_downscaled[c0] | rel 3.8e-10 |

In gate3, the first run of these checks, the t_soisno probe sat at exactly 273.15 K. The
central difference there straddles the phase-change switch (FD ∝ 1/h), so the probe now
picks a layer away from TFRZ. The all-leaf comparison treats a leaf whose largest
difference is ≤1e-15 as agreeing; gate3 had 1e-20 gradients there that were exactly 0
on the other backend.

**Finding (CPU host list), fixed in `ad/host-watertable-production`:** the host reverse
phase `watertable_rev_phase!` called the aquifer `water_table!` unconditionally and
skipped the bedrock clip. Production `clm_drv!` uses `perched_water_table!` +
`theta_based_water_table!` when `use_aquifer_layer = false` (the default, and what this
domain uses), and clips zwt at bedrock either way. So host-list gradients through `zwt` /
`qflx_drain` differed from production. In gate3 the wide loss with `zwt_col` and
`qflx_drain_col` seeded showed this as rel 1 in `qflx_drain_col`, `zwt_col` and
`column.zi`.

Both phases now follow production. The device aquifer branch also gains the bedrock
clip, which it lacked. `test/test_watertable_reverse.jl` checks, for both
`use_aquifer_layer` settings, that the host phase, the device phase (on the host) and
the production block give bit-identical forward state, that the host and device Enzyme
gradients through zwt / zwt_perched / qflx_drain agree to 1e-12, and that the host
gradient matches FD. `zwt_col` and `qflx_drain_col` are back in the gate's wide loss.
The default loss (soil temperature) does not depend on the water table, so
dL/d(t_grnd) is unchanged.

Gate on this branch (`gatewt1`, 37f1c6a, `CLM_REVERSE_CPU=1 CLM_REVERSE_FD=1
CLM_REVERSE_EXTRA=1`, 6759 s wall): **PASS**.
- dL/d(t_grnd[1]) is −3.44327925 on CUDA and CPU (absdiff 2.1e-14), unchanged; FD rel
  1.4e-9.
- Default-loss all-leaf CUDA vs CPU: 2.9e-14.
- Wide loss, now with `zwt_col` and `qflx_drain_col` seeded: all 3484 leaves CUDA vs
  CPU (host list) worst 1.1e-14, 481 nonzero entries. The seven FD probes are unchanged
  from gate4 (≤1.7e-8). With `use_aquifer_layer = false` the perched/theta-based water
  table does not write `qflx_drain`, so its seed passes straight through on both lists.
- CUDA phase diff of `watertable_dev_rev_phase!`: rel 4.0e-16, deterministic.

### Compile time

**Where it goes.** The launch rule now has an opt-in profile (`CLM._GPU_LAUNCH_PROFILE`,
reported per kernel by the phase diff) and the diff reports Julia's compile-time
counter. Device-kernel compilation is small: KernelAbstractions' Enzyme tape-type, aug
and rev kernels take 5–30 s per kernel (≈60 s of `cf_rev_init!`'s 1438 s; 16+9 s of
`cf_rev_friction!`'s 437 s). 95–99% of each phase's first call is compilation on the
host. A sampling profile of `cf_rev_friction!` (`PROFILE_PHASES=…`) puts 81% of the
first call in Enzyme's thunk generator (`get_staged`). About half of that is Enzyme
TypeAnalysis (`TypeTree::insert` / std::map iteration) over the argument structs. A
CuArray's nested `DataRef`/`Managed` fields give a much deeper type tree than an
Array's, which is why the same phase compiles 3–15× faster on CPU. Identical kernels
and repeated phases are already cached in-process (phases 5–7 cost 0 s; every kernel
compiles one signature).

**Changes** (numerics unchanged, see the table):

1. *Device-scoped Enzyme type lattice.* `_phase_autodiff!` sets
   `Enzyme.API.maxtypedepth!/maxtypeoffset!` to `CLM.GPU_REVERSE_TYPE_LATTICE` = (3, 64)
   (Enzyme defaults 6 / 512) around each per-phase Enzyme call whose bundle lives on a
   device, then restores the defaults. Host reverse calls are untouched. Set it to
   `nothing`, or set `CLM_GPU_REVERSE_TYPE_LATTICE=default`, to opt out. The full 36-phase
   diff at (3, 64) and at (6, 512) gave the same per-phase differences to every printed
   digit, both before the fix (pd364 vs pdall4) and after it (pdfix1 ran at 3/64). On
   `cf_rev_friction!` alone: 602 s (6/512), 343 s (4/512), 354 s (6/128), 292 s (3/64),
   with identical gradients (`scripts/cuda_reverse_typeknob_sweep.ps1`).
2. *Narrow view for the scatter-only gnet reduce phase* (`soiltemp_gnet_patch_reduce_view`).
   It contains only the arrays its four `scatter_add_1d!` calls touch, instead of every
   Gnet*In/Out struct: 1210 s → 73 s in the gate.

Not done, and why:
- **Caching across runs.** GPUCompiler 1.15's disk cache (`disk_cache` preference)
  only stores kernels whose CodeInstance comes from a package image; KA's Enzyme
  aug/rev kernels and Enzyme's host thunks are runtime compilations, so they never hit
  it. Enzyme has no persistent thunk cache. A precompile workload in a package
  extension would need a GPU at precompile time and would redo the hours-long compile
  on every source change.
- **Splitting the canopy energy kernel.** Its device compile is 18+10 s; the cost was
  host-side.

**Before/after (whole gate, fresh process, same server).**

| | gate2 (#373, 665a90b) | gate3 (this branch, 5ea6749) | gate4 (9bbf2b2) |
|---|---|---|---|
| CUDA reverse sweep (36 phases) | 13425 s | 4564 s | 4254 s |
| CPU reverse sweep (15 phases, unchanged) | 2014 s | 1956 s | 1963 s |
| whole run (incl. precompile, CPU leg, FD) | 16466 s | 7466 s, including the EXTRA legs | 7190 s, including the EXTRA legs |

Per CUDA phase in the gate, before → after: canopy energy 2324 → 118 s, canopy init
1193 → 206 s, canopy resist 1064 → 82 s, soil water 1096 → 108 s, gnet reduce
1210 → 74 s, soil-temperature post 1006 → 400 s, hydrology_no_drainage 836 → 797 s.

## Follow-up 3: CPU reverse-chain scripts after #382

All fifteen `scripts/calibrate_streamflow_ad*.jl`, `enzyme_driver_reverse*.jl`,
`enzyme_multistep_*.jl` and `enzyme_clm_drv_reverse.jl` were re-run on a laptop
(Apple Silicon, Julia 1.12.6, CPU only) at 8bcccbf. Every script whose chain contains
`watertable_rev_phase!` was also run at 822fda2, the commit before #382. All of them use
`CLMDriverConfig()`, so `use_aquifer_layer = false`.

**Unchanged by #382 (bit-identical printed gradients and FD):** `enzyme_driver_reverse_whole`
(dL/d(t_grnd) = −2.86850579e−8), `_fullstep` and `enzyme_clm_drv_reverse` (same value),
`enzyme_multistep_reverse` (2.36114321e+2), `_forced` (2.77574493e+2), `_gpp` Part A
(1.65082275e+2), `calibrate_streamflow_ad` (gates 1 and 2) and `_multi` (5/5 basins).
These scripts start from a warmed state in which the water table sits at bedrock. That
holds because production `clm_drv!` re-diagnoses zwt from an unsaturated profile and
clips it at bedrock, and leaves `qcharge = 0`. At that state the old aquifer
`water_table!` changes nothing. The new perched/θ-based diagnosis gives bedrock again,
and the clip has no effect. None of these losses reads `qflx_drain`. The perched water
table and the frost table only enter `satexcess_cal!` through a branch that is off in
these windows.

**Changed by #382: scripts that set zwt by hand.** `calibrate_streamflow_ad_real`,
`_seasonal` and `_snow` write `zwt = zi_bedrock − 1 m` to switch on the power-law
baseflow. The old host phase kept that zwt: with `qcharge = 0` the aquifer routine left
it alone. The production diagnosis overwrites it on the first step, so zwt returns to
bedrock, the head is 0, and dL/d(baseflow_scalar) is exactly 0 for both AD and FD. All
three gates then failed as VACUOUS. The old nonzero baseflow gradients were therefore
gradients of a water table frozen at the hand-set depth, not of the production model.
The scripts now set the wet state on the moisture profile (`wet_profile!`: saturate the
layers below zi_bedrock − head, then re-diagnose zwt). The baseflow sink drains the
first saturated layer (jwt+1, as the script comments always said) instead of the bedrock
layer, so the diagnosed water table falls as the aquifer drains. In the 7-day window
zwt goes from 1.325 to 1.554 m, and baseflow stays on in 336/336 steps.

| script | before #382 | after #382, unfixed | after fix |
|---|---|---|---|
| `_real` gate, ∂L/∂fff | −4.344527 (FD rel 3.9e−8) | 22.18363 (FD rel 2.2e−6) | −4.884137 (FD rel 6.5e−8) |
| `_real` gate, ∂L/∂baseflow_scalar | 9.598658e4 (FD rel 1.2e−13) | 0, VACUOUS → gate FAIL | 1.029781e5 (FD rel 9.5e−8) |
| `_real` 7-day calibration, KGE: start → AD / DDS | −1.668 → 0.723 / 0.490 | not run | −9.911 → −0.320 / 0.622 |
| `_snow` gate, ∂L/∂(θ_snow, fff, bf) | 4.678e−7, −2.909e−7, 2.646e−4 | 2.174e−7, −1.981e−7, 0 → FAIL | 4.356e−7, −2.800e−7, 2.141e−4 (FD rel ≤ 9.4e−9) |
| `_seasonal` Stage A, 30 d | −1.669739e2, 8.460133e5 | — | 1.183595e2, 1.439947e6 (FD rel 2.1e−10) |
| `_seasonal` Stage A, 90 d | 2.410578e2, 6.781069e6 | — | −2.523242e3, 2.365639e7 (FD rel ≤ 2.9e−10) |

The seasonal gate needed a smaller FD step. Over 30–90 days the production objective is
piecewise smooth in baseflow_scalar. The θ-diagnosed zwt is a piecewise-linear
interpolation between layer nodes, and the sink moves between layers. Within ±1 % of θ
the slope varies from 1.4e6 to 4.7e6. With the old 1e−3 relative step, FD differed from
AD by 0.41 (30 d) and 1.4e−2 (90 d). Steps of 1e−4 to 1e−6 agree with AD to 9 digits.
`CLM_FD_REL` now defaults to 1e−5. The AD calibration in `_real` now stalls at
KGE −0.32, behind DDS: it fits SSE, and the kinked objective favours derivative-free
search. That is a finding about this objective, not about the gradient.

**Other fixes (independent of #382).**
- `enzyme_driver_reverse.jl` [A]: `canopy_fluxes!` has seeded `z0mv_patch` from
  `canopystate.z0m_patch` since 0e23b3e. The fixture left `z0m_patch` NaN, so the
  real-canopy reference returned t_veg = NaN. With `z0m_patch = 0.5` (equal to its z0mv),
  [A] rel 9.8e−10, [B] 3.3e−9 and [C] 3.8e−7 all pass. `enzyme_canopy_mirror.jl`,
  `enzyme_canopy_psn.jl`, `enzyme_compositional_canopy.jl` and `enzyme_phs_canopy.jl`
  have the same fixture and were not re-run.
- `enzyme_driver_reverse_hydro.jl` [G]: the warmed state has `qcharge = 0` with zwt at
  bedrock. That is the kink between the aquifer routine's rising and deepening branches,
  so the reverse returned the one-sided value 0 while central FD gave −55.67. At
  `qcharge = 1e−4` (`CLM_QCHARGE0`) the check gives −55.6449479 for both (rel 1.1e−11).
- The default data paths in seven of the scripts now honour `SYMFLUENCE_DATA`.

**Still looks wrong (not fixed here).** The decomposed canopy reverse
(`cf_rev_init!` / `cf_rev_energy!`) differs from production `canopy_fluxes!` in two
ways.
- It never copies `z0m_patch` into `z0mv_patch` before the roughness blend.
- It never updates `efeb`, which production updates in its convergence kernel.

As a result `canopy_rev_converged_n` never sees `dele` shrink. On the Bow state in
`enzyme_clm_drv_reverse.jl` it returns itmax = 40, where production converges in
11 / 8 iterations. The sign-change limiter (`efe·efeb < 0`) also never fires in the
decomposed forward. The script FD checks still pass because they difference the same
decomposed forward. Fixing this changes the CUDA-validated canopy phases, so it needs
its own change and a GPU re-run.

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
PHASES=soiltemp_rhs,post julia --project=scripts scripts/cuda_reverse_phase_diff.jl
julia --project=scripts scripts/cuda_band_solve_reverse_check.jl
# follow-up 2: all-phase diff with determinism + FD arbitration, forward/alias diff,
# wider gate checks, compile profile of one phase
PHASES=all DETAIL=1 julia --project=scripts scripts/cuda_reverse_phase_diff.jl
julia --project=scripts scripts/cuda_reverse_forward_diff.jl
CLM_REVERSE_CPU=1 CLM_REVERSE_EXTRA=1 julia --project=scripts scripts/gpu_ad_reverse_driver_validate.jl
PHASES=cf_rev_friction PROFILE=1 julia --project=scripts scripts/cuda_reverse_phase_diff.jl
```
