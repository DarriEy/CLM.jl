# Three-site annual qualification — 2026-08-25

Status: **annual completion passes; strict scientific-parity gate fails**.

## Source and environment

- CLM.jl release base: `e5c3c183ebb318bdaaa89cf16fa1f458a39fbd20`
- Output-routing-only commit under test: `e60bda7`
- Julia: 1.12.6
- Platform: macOS 26.5, Apple M1, CPU Float64
- Inputs: staged `SYMFLUENCE_DATA` archive
- Configuration: the existing PHS+LUNA domain runner and Fortran restart initialization
- Normal balance checks remained active

These are qualification results, not yet definitive paper artifacts: the input archive and
full environment still need published checksums and licensing review.

## Annual completion

| Site | Role | Steps | Result | Runtime | Output SHA-256 |
|---|---|---:|---|---:|---|
| Bow at Banff | cold/snow | 8,760 | PASS | 191.7 s | `810ef3a8aaffa412abcd9fd42ef53b90c2da6ebdd87aa007fe6659db6324fbcb` |
| Aripuanã | wet tropical | 8,784 | PASS | 189.9 s | `2ef865c631021ab5609eb46791ce50cd4aae4d319655539ee55ebe76669aebde` |
| Stillwater | dry/semi-arid | 8,760 | PASS | 184.8 s | `88bd8c007c42129982a3f904d5e80de5e4bf612099c78a4737f9b5ef6e63af34` |

All three runs reached their declared final timestep and exited successfully. This
independently confirms that the earlier Aripuanã and Stillwater mid-year balance aborts are
closed on the current release base.

## Strict parity result

The existing strict gate compares aligned daily series and requires both annual-mean and
daily trajectory criteria. It was executed in enforceable `--check` mode for exactly the
three qualification sites.

- Passing cells: **204/207**
- Fully strict sites: **1/3** (Aripuanã 69/69)
- Undocumented failing cells: **3**
- Gate process exit status: **1**, as required for an unmet gate

| Site | Variable | Annual metric | Daily nRMSE | Failed component |
|---|---|---:|---:|---|
| Bow | BTRAN | +4.2 % | 0.21 | annual and daily |
| Bow | SNOW_DEPTH | +0.9 % | 0.12 | daily |
| Stillwater | SNOW_DEPTH | +15.8 % | 0.50 | daily |

Therefore claim C01 (strict agreement with the declared Fortran reference) remains TBD and
cannot be asserted for the qualification matrix. Annual completion claim C02 has supporting
qualification evidence but is not final until inputs and run metadata are fully archived.

## Immediate investigation

1. Confirm date alignment, units, aggregation, fill handling, and the exact mapped Fortran
   diagnostics for `BTRAN`/`BTRANMN` and `SNOW_DEPTH` before changing physics.
2. Plot daily residuals and locate season/state transitions responsible for each nRMSE.
3. Determine whether BTRAN compares unlike diagnostics (instantaneous/mean/minimum,
   patch/column aggregation) or reflects a model-path difference.
4. For snow depth, compare SWE, snow ice/liquid, snow fraction, layer count, density, melt,
   and the distinction between snow-covered-area and gridcell-average depth.
5. Use process-boundary dumps around the first divergence; fix only demonstrated defects.
6. Keep thresholds frozen during this investigation. Do not add exceptions until a physical
   mechanism and acceptability argument have independent scientific review.
