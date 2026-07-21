# Metal full-fleet re-validation — 2026-07-20

**Supersedes the scope (not the results) of `GPU_REVALIDATION_2026-07.md` (#258).**

## Why this run

#258 validated the device paths whose `src/` file had changed since #242 and
concluded "no regressions found; no `src/` change was required." That verdict was
correct *for what it ran*. But changed-files scoping cannot see a harness that was
**already broken** — a harness that is never invoked is indistinguishable from a
passing one in a summary that never lists it.

So this run executed **every** `gpu_validate_*.jl` harness (123, excluding the
NVIDIA-only `gpu_validate_cuda.jl`) on the M4 Max at HEAD. Each takes ~15-25 s;
the whole fleet is ~50 minutes. It is cheap, and it should be the default.

Motivation was also concrete: 19 `src/` files changed after #258 across 20 commits,
including five PRs that flipped **default flags** (`use_aquifer_layer` #259,
`create_crop_landunit` #265, `use_luna` #267, `use_hydrstress` #273, `use_bedrock`).
A flipped default silently re-points a harness at a different code path.

## Result

| | count |
|---|---|
| harnesses run | 123 |
| passed first time | 114 |
| **failed** | **9** |

Of the 9: **8 fixed here**, 1 blocked on missing input data (not a defect). The last
one (`watertablevariants`) was initially logged as an open device divergence but
turned out to be a harness fixture bug and is now also fixed — see the resolved item
below. A separate, softer CPU-side issue (`combine_snow_layers!` optimization
dependence) remains open and is not a device-parity failure.

Three of the failures were already failing at the #258 commit (`a99a4ce`) — verified
by re-running them there in a detached worktree. They are not regressions from the
recent commits; they were simply never run.

## The three real `src/` defects

### 1. `combine_snow_layers!` could not compile on Metal at all

`_snowhyd_combine_kernel!` received **raw Float64 `dtime`** and did `dt = T(dtime)`
*on device*. Metal has no double at all, so this was
`InvalidIRError: unsupported use of double value` — the kernel never compiled.

The kernel header comment claimed literals were converted "so no Float64 reaches a
Float32-only backend". It covered literals and missed a scalar **argument**.

Fix: `eltype(dz)(dtime)` at the launch, matching the sibling launches already doing
exactly this in the same file (lines 396, 663, 825). Identity on Float64, so the CPU
path is byte-identical.

### 2 & 3. #268 broke the non-PHS photosynthesis kernel on device (two ways)

#268 ("port the LUNA branch into the NON-PHS photosynthesis!") added three arguments
to `_psn_pass2_kernel!` — `crop_pft`, `vcmaxcint_sun`, `use_luna`. Both new arrays
default to a host `Float64[]` (photosynthesis.jl:2350, :2366), and a kernel argument
is rejected for being non-bitstype **whether or not the body reads it**. Then, once
that was fixed, the added arguments pushed the kernel to 32 indirect argument
buffers against Apple's hard limit of 31:

```
NSError: Total number of indirect argument buffer resources exceeded
         for buffers (32/31) (AGXMetalG16X, code 3)
```

Net effect: **the default (non-PHS) photosynthesis path could not run on Metal.**

Fixes: substitute a backend-matched zero array when either is empty (the guard #273
added for `crop_pft` one file over in `canopy_fluxes.jl`), and fold the six loose
per-PFT vectors into one `_PsnP2Pft` bundle — the same device-view-bundle trick the
PHS pass (`_Psn2Pft`) and `fire_base.jl` already use. `_psn_pass2_body!`, which the
AD host loop calls directly, is unchanged.

## The dominant harness bug class: CLONE DRIFT `test/` → `scripts/`

Four failures were a `src/` change that updated the `test/` fixture and not the
`scripts/` GPU clone:

- **`canopy_e2e`** — NaN on every field. Bisected (852 commits, ~10 steps) to
  `0e23b3e`, which made `canopy_fluxes.jl:1709` re-seed
  `z0mv_patch .= canopystate.z0m_patch` each step. It updated
  `test/test_canopy_fluxes.jl:102` and `test_canopy_reverse.jl` but not this clone,
  so the harness's `z0mv = 0.5` was overwritten by `z0m_patch`'s **NaN init**
  (`canopy_state.jl:124`); the `max(z0mv, z0mg)` floor cannot rescue a NaN. After the
  fix the harness returns byte-identical pre-regression numbers (t_veg 0.0).
- **`clmdrv_e2e` + `clmdrv_cn_e2e`** — both died before any device comparison because
  #273 flipped `use_hydrstress` on and `scripts/clmdrv_make_data.jl` lacked the PHS
  pftcon params. `test/test_clm_driver.jl:177-179` had been updated; the clone —
  whose own header says "Kept in sync with the test" — had not.
- **`gpu_validate_lake`** (18 args vs 24) and **`gpu_validate_soiltemp`** (four
  separate stale calls) were plain arity mismatches against evolved kernel
  signatures. Both had been dead since those signatures changed, pre-#258.

**Takeaway:** when a `src/` change adds a fixture obligation or a kernel argument,
grep `scripts/` for the same struct/kernel, not just `test/`.

## Coverage added

- **`scripts/gpu_validate_tridiag_jtop.jl`** — `tridiagonal_multi!` on the ragged
  `jtop > 1` path. This is a live production configuration
  (`lake_temperature.jl:1219` passes `jtop .+ nlevsno`) with **no** direct device
  coverage: `gpu_validate.jl:88` pins `jtop = ones()`. #263's message states the
  pre-fix code "corrupted the jtop > 1 path outright".
  It compares against an **independent dense `LinearAlgebra.Tridiagonal` oracle, not
  the KA-CPU result** — the whole lesson of #263 is that the CPU side can itself be
  wrong, so a device-vs-host check measures both against a possibly-broken yardstick.
  Host 1.3e-16, Metal 8.8e-08 at nlevs 21 and 25.
- **`scripts/gpu_validate_fire_transient_e2e.jl`** — closes the single gap #258 left
  open. `transient_landcover` appeared in **zero** `gpu_validate_*` scripts, and the
  existing fire harness both omits the kwarg and builds `fbac_col` as zeros, so
  flipping it there would prove nothing. The new harness asserts the two regimes
  disagree on the host, pins device == host-transient (~1e-10), and separately
  asserts device != host-non-transient (2.1e-03) so a wrong-branch device cannot pass.

Both harnesses also pass under `--check-bounds=yes`.

## Still open

1. ~~`gpu_validate_watertablevariants_e2e` — genuine device divergence.~~
   **RESOLVED — it was a harness fixture bug, not a device/port bug.** The initial
   read here ("device NaN, host finite; not a fixture-shape problem") was wrong on
   both counts. Reproducing it showed **host = NaN, device = 0.0**, and the cause was
   fixture layout: both water-table kernels index the snow+soil column arrays with the
   standard snow offset (soil layer k at padded index `k + nlevsno`, interfaces at
   `k + nlevsno + 1`) — what the real caller passes (`col.dz/col.z/col.zi` are
   nc×(nlevsno+nlevgrnd), init_vertical.jl:152-160) and how every state field is laid
   out. The harness built `col_dz/col_z/col_zi` as soil-only nc×nlevsoi and wrote
   h2osoi_liq/ice and t_soisno at `[c,k]` instead of `[c,k+nlevsno]`, so the geometry
   reads ran off the end of the 20-wide arrays: CPU `@inbounds` returned adjacent-heap
   garbage (NaN), Metal's OOB read returned 0.0. h2osoi_liq/ice were mis-written but
   stayed in bounds of their 37-wide arrays, so `h2osoi_vol` matched device==host on
   shifted-but-consistent values — a green 0.0 over physically wrong inputs.
   Fix (harness only, no `src/` change): build geometry snow-padded, write soil
   quantities at `k+nlevsno`, keep `watsat_col` soil-indexed. All four columns now pass
   with finite, per-column-distinct depths (zwt 0.05 / 1.39 / 27.1 / 104.5 m across the
   four branches), and it passes under `--check-bounds=yes` — the mode that would have
   thrown a BoundsError on the original OOB instead of silently reading garbage.

2. ~~`combine_snow_layers!` gives an optimization-dependent answer on the KA CPU
   backend.~~ **RESOLVED — a real #263-class store→load miscompile, now fixed.** Under
   `--check-bounds=yes` the host returned `snl = [0, 0, -3, -3, 0]` vs `[0, -2, -3, -3, 0]`
   under normal bounds; the device was stable and correct in both. The physically-right
   answer is unambiguous (col 2's two snow layers have ice 10.0 and 0.5, both far above
   the 0.01 removal threshold → snl = -2), so no oracle was needed. Root cause: the
   snow-depth recompute (`snow_hydrology.jl`, `_snowhyd_combine_kernel!`, ~line 1652)
   accumulated into the array cell — `cv.snow_depth[c] = 0; for jl: cv.snow_depth[c] += m.dz[c,jj]`
   — the textbook reduction-into-a-memory-location the KA CPU backend miscompiles.
   Under `--check-bounds=yes` it kept only the last layer's dz, so a healthy 0.054 m
   two-layer pack read as 0.004 m, tripped the `frac_sno_eff*snow_depth < dzmin`
   all-snow-gone branch, and collapsed the column to snl = 0. **This is a genuine
   production CPU correctness hazard, not a CI artifact** — any compile that
   reassociates the reduction can wrongly collapse a multi-layer snowpack; the device
   was always correct. Fix: accumulate into a local scalar, write the cell once at the
   end (the #263 register-carry pattern). All three CPU modes + device now agree on -2,
   and the harness now asserts the expected `snl` vector, not just device==host parity.
   A sweep of the file's other `[c,jj] +=` sites found only downward mass-transfers
   (distinct source/dest cells, RHS a local) — not reductions, so LLVM won't reassociate
   them; left as-is.

3. **`gpu_validate_clmdrv_fates_e2e` cannot run here** — needs
   `domain_Aripuana_Amazon/settings/CLM/parameters/surfdata_clm.nc`, absent on this
   box. Environment, not a port defect.

4. **`gpu_validate_clmdrv_e2e` passes but is partly blind.** It reports 4 fields as
   `[skip] (no finite entries on both)`, including **`t_veg_patch`** — the
   whole-driver harness does not currently exercise the canopy path, which is exactly
   where defect (1) above and the `canopy_e2e` NaN lived.

5. **The "~15 latent host `Ref` globals" figure from #263 does not hold up.** A full
   call-graph census found exactly **2** kernel-reachable (`SNOW_THRESH_K`
   snow_hydrology.jl:1975 and `SNOW_COMBO_TEMP_K` :2069, both from
   `_snowhyd_divide_kernel!`); the rest are host-only, statically eliminated behind
   `_use_smooth`/`_is_ad_type`, or FATES host code. The original count was probably a
   plain grep for `Ref`. Moreover, once `combinedivide` compiled, **neither predicted
   deref faulted on Metal** — the divide kernel ran to exact 0.0 parity. Still worth
   checking on CUDA/AMDGPU; not a confirmed Metal defect.

## CPU regression check

Full suite, **all three `src/` changes** (snow_hydrology.jl + both photosynthesis.jl):
**26588 pass / 0 fail / 3 broken** (pre-existing `@test_broken`), exit 0, 27m40s —
byte-identical pass count to the pre-change baseline (26588). The two new harnesses
and every fixed harness also pass under `julia --check-bounds=yes`.
