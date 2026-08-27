# Calibration qualification — 2026-08-27

Status: **single-parameter twin experiment PASS; two-parameter recovery FAIL**.

The original two-parameter test accepted optimizer termination and objective reduction but
did not assert recovery of its generating parameters. It therefore printed passed after
returning Medlyn slope 6 instead of 8 and Vcmax scale 1.0 instead of 1.3.

The gate now requires:

- a finite, nonzero objective at the starting parameters;
- a finite, nonzero starting gradient;
- each recovered physical parameter within 10% of its generating value;
- the existing objective-reduction and iteration bounds.

The locked Linux focused run produced 15 pass and 4 fail across 19 assertions:

- override gradient wiring: 4/4 pass, with AD/FD gradient ratio 1.0;
- single-parameter `csoilc`: 4/4 pass, recovering 0.012 from a 0.004 start to
  effectively zero relative error in four iterations;
- two-parameter Medlyn/Vcmax: 3 pass, 4 fail;
- multi-timestep evaluation: 4/4 pass.

For the failed twin experiment, the defaults already produce objective 0 and gradient 0
against observations generated using the declared non-default truth. Optimization therefore
terminates after one iteration at Medlyn 6.0 and Vcmax scale 1.0, with relative parameter
errors of 25.0% and 23.1%. The current LH/GPP target extraction does not identify these
parameters in this initialized four-step configuration.

The failure is deliberately retained. The submission may demonstrate calibration plumbing
and the successful scalar recovery, but it may not claim joint Medlyn/Vcmax identifiability
or recovery until a prespecified, signal-bearing observation design passes this gate.
