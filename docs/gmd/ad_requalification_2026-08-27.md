# AD winter-scenario requalification — 2026-08-27 (evening)

Status: **winter FAIL root-caused to test-harness state aliasing, fixed in
`_calib_dual_copy`; awaiting full `test_ad_robustness.jl` confirmation run.**

## What the morning FAIL actually was

`ad_qualification_2026-08-27.md` recorded the winter (250 K) scenario failing all nine
finite-reference assertions: the AD run at the nominal temperature was finite, but both
finite-difference endpoint trajectories emitted non-finite water/snow state. The record
interpreted this as the perturbed smoothed-model trajectories being non-finite.

Localization on the Bow fixture shows the failure has **nothing to do with the
perturbation**:

- Three IDENTICAL dual runs (dT = 0 each) from one warmed-up Float64 state: run 1 is
  finite, runs 2 and 3 are NaN. Perturbations of ±1e-6 … ±1 K behave identically —
  whichever dual run happens first is finite, every later one NaN.
- Phase probes show run 2's state is finite through the canopy solve and goes all-NaN
  inside `soil_temperature!` (28 `t_soisno` entries at once), then propagates to the
  water state.

Mechanism: `_calib_dual_copy` copied only `Array{Float64}` fields into a dual instance;
**integer, Bool, and Bit arrays were shared by reference** across every dual instance
built from the same Float64 base. In the winter scenario, dual run 1's snow-layer
initiation flips the shared `col.snl` from 0 to −1 (writing the new layer's dz/z/zi and
temperature only into run 1's own dual arrays). Dual run 2 then starts from a fresh copy
of the Float64 base — whose snow-layer slots are NaN-padded because that base never
initiated a layer — but with the mutated shared `snl = −1`, so `soil_temperature!` reads
a snow layer that was never initialized. Only the winter scenario snows, which is why
only winter failed; and the AD (first) run always passed while both FD (later) runs
always failed, which reproduced exactly the morning's 10-pass/9-fail pattern.

## The fix

`_calib_dual_copy` now copies `Array{<:Integer}`, `Array{Bool}`, and `BitArray` fields
(`copy(sv)`) instead of aliasing them — the same isolation the Float64 state always had.
Default-path physics is untouched: the function only builds dual instances for AD.

With isolated state, repeated identical dual runs and perturbed runs (±1e-3, ±0.01 K)
all remain finite on the previously failing fixture.

## Claim consequence (pending confirmation)

The winter failure was an artifact of the qualification harness's shared-state pattern,
not a model-physics fragility of the smoothed winter trajectory. C03's exclusion of the
winter regime is lifted on macOS and awaits the locked-Linux confirmation:

- [x] `test_ad_robustness.jl` full-file result on macOS (this box, Apple M5,
      Julia 1.12.6): **150 pass / 0 fail, all six scenarios** (3m28.8s). Winter
      (250 K) AD vs smoothed-FD agreement: d(LH)/dT = 0.036 / 0.036,
      d(SH)/dT = −11.8525 / −11.8525, d(T_grnd)/dT = 0.4691 / 0.4691 —
      0.0% relative error on all three, versus the morning's 9 non-finite
      references. (The total rises from 144 to 150 because winter's guarded
      sign/agreement comparisons now run instead of being skipped on NaN.)
- [ ] locked-Linux rerun as part of the definitive campaign
