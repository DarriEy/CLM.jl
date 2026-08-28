# Accelerator qualification (E8/E9) — 2026-08-28

Status: **E8 device-correctness fleet PASS on Apple Metal — 127/129 exit-zero, the two
non-zero rows CUDA-hardware-gated; E9 scaling run recorded separately.**

## Hardware and environment

- Apple M5, Metal 4, device `AGXG17GDevice` — a **different Metal device generation** than
  the July 2026 fleet's `AGXG16CDevice`, so the accelerator claim now has two distinct
  Metal devices behind it.
- macOS 26.5 (Darwin 25.5.0), Julia 1.12.6, Metal.jl 1.9.3; harness project
  `scripts/Manifest.toml` at the qualification commit.
- Working precision on device: Float32 (the production Metal path); host Float64.

## E8 — device==CPU correctness fleet

Runner: `scripts/gmd/run_metal_fleet.sh` (structural exit-status classification only —
no log-text parsing; every harness always runs and appears as a row). Summary:
`repro/results/metal-fleet/metal_fleet_1dcd5c9.json`, produced by
`scripts/gmd/summarize_metal_fleet.py`.

Result at commit `1dcd5c9` (clean tree):

- **129 harnesses ran: 127 exit-zero, 0 unexplained failures.**
- The 2 non-zero rows are hardware-gated by design on a Metal-only box:
  `gpu_validate_cuda` and `gpu_ad_reverse_driver_validate` (both CUDA-targeted).
- The whole-driver composites all pass, including `clmdrv_cn_e2e`
  (device `clm_drv!` with `use_cn=true` matching host over **754 finite outputs**, max
  rel 2.23e-8 on `t_soisno`, `decomp_cpools_vr` and `leafc` bit-identical) and
  `clmdrv_fates_e2e`, which had been input-blocked on every previous box.

## What the qualification (shakedown) pass caught first

The first full-fleet pass on this box was a deliberate shakedown; it flushed four real
defects, all introduced after the July fleet and invisible without real Metal hardware:

1. `update_lnd2glc!` scalar-indexed device arrays (defaults loops + host column loop) —
   faulted all whole-driver Metal composites. Kernelized behind the `isa Array` host gate.
2. `lnd2atm_rof!` stream-zeroing and dynbal-correction loops scalar-indexed device
   arrays. Replaced with view broadcasts (identical semantics, unset dynbal ≡ 0).
3. The `downscale_forcings` harness's byte-identity reference still used the pre-#328
   per-column pco2 rescale (1 ulp off; compounding for multi-column gridcells) — the
   reference was stale, the src behaviour intended. Reference updated; harness green
   over 39 finite outputs.
4. The `soiltemp` harness called `_phase_change_beta_kernel!` without the
   `PHASE_CHANGE_TEMP_K` scalar added by #321 (MethodError before anything ran). With
   the real launch's argument list it passes 12/0.

(1) and (2) are src fixes with byte-identical CPU behaviour; (3) and (4) are harness
repairs of clone/signature drift — the documented `test/`→`scripts/` drift class.

## E9 — performance scaling

Runner: `scripts/gmd/run_metal_scaling.jl` (median + IQR over 11 repetitions after one
discarded warm-up; compilation excluded; host-device transfers excluded — state resident
each side, device timing includes a full synchronize; device, threads, and precisions
recorded in the JSON). Result: `repro/results/metal-fleet/metal_scaling_1dcd5c9.json`.
The crossover column count, not only the peak speedup, is part of the record.

## Claim consequence

- C06 (accelerator reproduces CPU pathway): qualification-pass on Apple Metal at commit
  `1dcd5c9`. The claim scope is the tested devices (M-series Metal, two generations);
  CUDA remains implemented-but-untested on this hardware and must be claimed only as
  backend-capable design unless separately qualified.
- C07 (crossover/scaling): evidence recorded by the E9 runner; threshold-free (the
  record reports the measured crossover and dispersion, no pass/fail tuning).
