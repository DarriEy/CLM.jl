# Calibration requalification — 2026-08-27 (evening)

Status: **joint Medlyn/Vcmax twin experiment PASS — C05 gate repaired at the cause**.

## What the morning FAIL actually was

The morning record (`calibration_qualification_2026-08-27.md`) diagnosed the two-parameter
twin experiment as lacking an identifiable observation design: objective and gradient were
exactly zero at the defaults against observations generated with the declared non-default
truth. The observation design was not the defect. The truth run and the default run were
**bit-identical** because the two parameters never reached the physics in the default
(`use_hydrstress`, `use_luna`) configuration:

1. **PHS call-site drop.** `canopy_fluxes_core!`'s `use_hydrstress` branch called
   `photosynthesis_hydrstress!` without forwarding the `CalibrationOverrides`, so the PHS
   solve rebuilt its default all-NaN overrides and `medlyn_slope`/`vcmax25_scale`/
   `jmax25top_sf` were silently inert. The non-PHS branch forwarded them correctly —
   which is why `csoilc` (consumed inside `canopy_fluxes_core!` itself, which held the
   real overrides) calibrated fine while the photosynthesis parameters did not.
2. **LUNA bypass.** With `use_luna` on, the operative vcmax25 is the LUNA-acclimated
   `vcmx25_z`; `vcmax25_scale` only multiplied the static `vcmax25top`, whose scaled value
   is stored diagnostically and then overridden. A ×3 scale left every flux bit-identical
   (its only live consumers were the C4 `kp` path and the `kn` profile, both inert for a
   C3 tree under LUNA's lmr override). Both pass-2 bodies now scale the acclimated
   vcmax25 under the same `!isnan` sentinel guard.

Verification of the mechanism (Bow fixture, patch 2, PFT 1, 3 warmup steps, T=285 K):
medlyn 6→20 moves rssun 281→78 s/m and LH 152.6→194.7 W/m²; vcmax25_scale ×3 moves
`vcmax_z_phs` 10.60→31.24 and psnsun 4.28→5.60 (sub-linear through the `aj`
light-limitation, as expected physically). Before the fix all of these were bit-identical.

Default-path safety: both fixes sit behind the NaN-sentinel guard (`NaN` = use default),
so the default configuration is unchanged — the same guard convention every other
override uses.

## Gate result after the fix

`test/test_parameter_recovery.jl` (the strict gate introduced in `96609ab`, unchanged):

- **19 pass / 0 fail** (macOS, Julia 1.12.6, commit `6f54aee`);
- override gradient wiring 4/4, AD/FD ratio 1.0;
- single-parameter csoilc recovery 4/4 (recovers 0.012 from 0.004);
- **two-parameter Medlyn/Vcmax recovery: all recovery assertions pass**, including the
  previously failing finite/nonzero starting objective and gradient, and both parameters
  recovered within 10%;
- multi-timestep evaluation 4/4.

The Linux locked-environment rerun of this file remains part of the definitive campaign.

## Claim consequence

C05 ("Gradient calibration recovers synthetic parameters") moves from FAIL-QUALIFICATION
to qualification-pass on macOS, pending the locked-Linux confirmation. The morning record
is retained; this record supersedes its diagnosis, not its data.
