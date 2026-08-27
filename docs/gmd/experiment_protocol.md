# CLM.jl GMD experiment protocol

Status: **draft — thresholds and experiment membership are not frozen**.

## 1. Scientific questions

1. Does the scoped Julia implementation reproduce the configured Fortran CLM5 reference
   across process boundaries and complete annual trajectories at contrasting sites?
2. What error is introduced by numerical reformulation, floating-point precision, and the
   smoothing used to enable differentiation?
3. Are derivatives of selected model responses correct away from physical discontinuities,
   and are they useful in a controlled calibration problem?
4. Does the same scoped implementation execute correctly on the demonstrated accelerator,
   and at what workload size does acceleration become beneficial?

## 2. Scope

### In scope

- single-site/single-column CLM5 biogeophysics configuration;
- energy, water, snow, soil-state, and photosynthesis diagnostics exercised by that path;
- full declared simulation period with normal conservation checks active;
- Julia CPU Float64 as the primary translated implementation;
- differentiable variant and its quantified smoothing bias;
- accelerator backend(s) actually run in the definitive campaign.

### Out of scope unless promoted before protocol freeze

- global or coupled-CESM production readiness;
- universal CTSM configuration equivalence;
- biogeochemistry, methane, fire, crops, dynamic vegetation, FATES, and other optional
  subsystems not exercised by the primary experiment matrix;
- observational validation of CLM5 physics itself (the core comparison is implementation
  fidelity to the declared source model);
- claims about all NVIDIA/AMD/Apple devices based on one tested device.

## 3. Implementations

- **Reference:** CTSM tag `ctsm5.3.012`, commit
  `ab466d6f9789ca3df2c72bda46cf7afed2d04102`; the current macOS build and oracle
  configuration are recorded in `repro/configs/bow_btran_oracle.toml`. The prescribed
  macOS GNU Fortran 15.3/MPICH 4.3.2 `-O`, Linux GNU Fortran 14.3 `-O`, and Linux GNU
  Fortran 14.3 trailing-`-O0` sensitivity cells completed on 2026-08-27. The locked Linux
  image digest and exact cell results are recorded there; the acceptance threshold remains
  unfrozen pending independent scientific and numerical review.
- **Candidate:** exact CLM.jl submission commit TBD; Julia CPU Float64.
- **Differentiable candidate:** same commit/configuration with explicitly listed smoothing
  and AD mode.
- **Accelerator candidate:** same commit/configuration, backend/device and precision TBD.

## 4. Site matrix

Qualification must include at least:

- one cold/snow-dominated site;
- one wet tropical site;
- one dry or water-limited site.

The final breadth matrix may retain the existing biome collection only after all inputs,
initial states, licences, site-selection logic, and completion status are audited. Sites
are selected by prespecified climate/land-cover axes, not by whether they pass.

For each site record coordinates at policy-compliant precision, land cover/PFT composition,
forcing source, period, timestep, calendar, surface/parameter data, spin-up, restart, and
all checksums.

## 5. Equivalence controls

Before comparing outputs, prove that reference and candidate share:

- forcing values after any unit conversion and temporal interpolation;
- surface, parameter, and initial-state values;
- calendar, timestep, start/end timestamps, and leap-day handling;
- enabled process options and numerical tolerances;
- output definitions, units, sampling frequency, and aggregation semantics.

Include input-equivalence checks in the machine-readable results rather than relying on
manual inspection.

## 6. Metrics

Metrics are computed per site and variable at native comparison frequency and for declared
aggregates:

- mean bias and normalized mean bias;
- RMSE and a robustly normalized RMSE;
- mean absolute error;
- Pearson correlation where variance is nonzero;
- errors in annual total/mean, extrema, and timing where physically relevant;
- water and energy conservation residuals;
- maximum process-boundary absolute and relative error for oracle experiments.

Near-zero denominators use variable- and unit-specific absolute floors. Temperature errors
are reported in kelvin, not relative percent. Flux totals use physically meaningful units
and time integration. Missing, NaN, Inf, and constant-series cases receive explicit
verdicts and cannot disappear from aggregates.

Exact thresholds remain TBD until reviewed by a land-model scientist and frozen before the
definitive run. Both passing and failing cells will be published.

For process-oracle tolerances, exact agreement is required for discrete control flow,
configuration, injected state, and declared inputs. Numerical output tolerances will be
set only after running the same oracle with the locked Linux compiler and an independent
compiler/optimization configuration. That required matrix is complete: the grass absolute
error spans only `1.0272300575098203e-5` to `1.0272300575597804e-5` across the three cells.
The tolerance study reports the full range and must not choose a threshold from whichever
configuration most closely matches Julia. Independent review is now the remaining decision.

## 7. Experiments

### E1 — input equivalence

Compare all inputs presented to the first model timestep and checksum configuration files.

### E2 — process-boundary oracle

At representative seasonal states, compare inputs/outputs around major radiation,
turbulence, hydrology, snow, soil-temperature, and photosynthesis boundaries. Reference
dumps must be generated independently of Julia and instrumentation must be state-neutral.

### E3 — annual qualification

Run the three contrasting qualification sites for the complete declared year. Conservation
checks remain active. Any early termination fails qualification.

### E4 — multi-biome breadth

Run the frozen site matrix and report every site-variable verdict plus distributions of
daily and aggregate error. No site is removed after results are seen.

### E5 — smoothing bias

Run exact and smoothed Julia variants from identical states. Report state/flux differences,
conservation impact, and sensitivity to smoothing parameters at transitions and annually.

### E6 — derivative verification

For prespecified parameters and representative states, compare AD directional derivatives
with a centred finite-difference step-size sweep. Separate smooth interior states from
known physical discontinuities. Report convergence curves, absolute/relative derivative
error, non-finite results, and zero-gradient pathologies.

### E7 — calibration demonstration

Required minimum: a synthetic-twin recovery experiment with held-out responses. A real-data
case is included only if its scientific setup, uncertainty treatment, and data licence are
sound. Compare gradient-based calibration with a declared derivative-free baseline using
matched model-evaluation budgets.

### E8 — accelerator correctness

Compare CPU and accelerator results at matched precision first; compare each with the
Fortran reference separately. Report conservation and every discrepancy outside the frozen
precision-appropriate thresholds.

### E9 — performance scaling

Use logarithmically increasing column counts. Record hardware, power mode, software,
precision, threads, warm-up, compilation exclusion, host-device transfers, repetitions,
statistic, dispersion, and peak memory. Report the crossover, not only maximum speedup.

## 8. Reproducibility and results schema

Every run emits a record containing:

- experiment/run ID and UTC timestamp;
- Git commit and dirty-tree assertion;
- implementation/backend, precision, threads/ranks, device;
- complete configuration and input checksums;
- Julia/compiler/package environment identifiers;
- command and random seed, if applicable;
- start/end status, elapsed time, and output checksums;
- metrics and verdicts with threshold-set version.

One generated results index is the only source for manuscript numbers and plots. Analysis
and plotting must fail if expected runs are missing or provenance checks fail.

## 9. Freeze procedure

The protocol freezes only after scientific and numerical review. Freeze by committing:

1. exact experiment membership and configurations;
2. metric definitions and unit conversions;
3. threshold table and scientific justification;
4. analysis version and expected result schema;
5. release-candidate commit and locked environments.

Any later scientific change increments the protocol version and requires a new definitive
campaign. Superseded results remain archived.
