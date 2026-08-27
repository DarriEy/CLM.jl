# Automatic-differentiation qualification — 2026-08-27

Status: **winter scenario FAIL; differentiable seasonal-cycle claim blocked**.

The Linux strict-suite audit showed that `test_ad_robustness.jl` printed the winter scenario
as passed even though all three finite-difference references were `NaN`. The comparison
guards only ran when `abs(fd_deriv)` exceeded a threshold; comparisons with NaN are false,
so the assertions were skipped accidentally.

The test now requires each positive/negative finite-difference endpoint and the resulting
derivative to be finite before testing sign or relative agreement. The focused locked-Linux
run produced:

- overall: 135 pass, 9 fail, 144 total;
- winter cold (250 K): 10 pass, 9 fail;
- near-freezing, spring, summer hot, summer wet, and dry stress: 125/125 pass.

The nine winter failures are the positive endpoint, negative endpoint, and central
difference for each of latent heat, sensible heat, and ground temperature. The AD values
and derivatives themselves are finite, but the perturbed smoothed-model trajectories emit
non-finite water/snow balance state. Without a finite independent numerical reference, this
scenario does not validate the AD derivative.

The failure is deliberately retained. Submission text may claim demonstrated AD only for
the five passing forcing regimes unless the winter trajectory is repaired and independently
requalified. It may not claim seasonal-cycle robustness from the current experiment.
