# Three-site annual qualification — 2026-08-25

Status: **annual completion passes; snow residuals are diagnostically classified; Bow
BTRAN has a reproducible process oracle but the annual parity cell remains open**.

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

Therefore claim C01 (strict agreement with the declared Fortran reference) remains failed
until the Bow annual BTRAN trajectory is regenerated with the now-audited controls and
meets a prespecified criterion. The target shared-state oracle is now reproducible, but it
does not retroactively qualify the existing annual output. Annual completion claim C02 has supporting
qualification evidence but is not final until inputs and run metadata are fully archived.

## Fresh residual audit

The submission branch adds `scripts/gmd/audit_qualification_residuals.py`, which reads only the
new qualification outputs and staged Fortran references. It independently reproduced the
residual mechanisms rather than importing historical paper results.

`SNOW_DEPTH` is CTSM's height over the **snow-covered area**. CTSM's gridcell-mean snow
height is the separate `SNOWDP` diagnostic. The two failed covered-area cells have passing
gridcell-mean, snow-mass, and physical-pack comparisons:

| Site | Covered-area depth | Covered nRMSE | Gridcell `SNOWDP` | Grid nRMSE | SWE | Physical-pack result |
|---|---:|---:|---:|---:|---:|---|
| Bow | +0.851% | 0.123 | -0.385% | 0.0188 | -0.501% | covered depth +0.852%; density -0.531% |
| Stillwater | +15.757% | 0.503 | -0.243% | 0.00386 | -0.003% | 12 days: depth -0.170%, nRMSE 0.0070; density +0.136% |

Stillwater has 60 days where SWE is positive but at or below the gate's 0.02 mm physical
floor. Bow's largest covered-area depth differences occur from 26 September through 4
October when snow cover is only about 2–3.5%; multiplying covered-area height by snow cover
reproduces `SNOWDP` and removes the daily failure. Both `SNOW_DEPTH` cells are therefore
retained visibly as documented secondary-diagnostic exceptions, while `SNOWDP` is the
appropriate gridcell-scale paper variable.

The fresh Bow audit localizes BTRAN to 29 April–6 May (largest difference on 30 April:
Julia 0.33648, Fortran 0.21410). Yet the annual consumers pass: QVEGT and FCTR -0.674%,
latent heat +0.574%; total soil liquid differs +0.087% and SWE -0.501%. This supports the
daily-minimum amplification hypothesis but does not prove instantaneous implementation
parity. The historical n11616 oracle remains excluded. A replacement oracle has now been
generated from CTSM tag `ctsm5.3.012` (commit
`ab466d6f9789ca3df2c72bda46cf7afed2d04102`) at step 2901 in the worst target
window. Its configuration and checksums are in `repro/configs/bow_btran_oracle.toml`.
After complete state/control replay, tree BTRAN agrees to `9.78e-13` absolute error and
grass differs by `1.03e-5`, amplified from approximately `1e-8` canopy-energy arithmetic.

## Immediate investigation

1. Completed 2026-08-27: repeat the target oracle with locked Linux GNU Fortran 14.3
   release and independent trailing-`-O0` builds. Across macOS release, Linux release, and
   Linux `-O0`, grass absolute error spans only `1.0272300575098203e-5` to
   `1.0272300575597804e-5`; optimization sensitivity is not the cause.
2. Obtain independent land-model and numerical-method review and prespecify the disposition
   of the stable grass residual before deciding whether it is within tolerance. Do not tune
   a threshold to any one cell.
3. Add `SNOWDP` explicitly to the frozen paper variable table and label `SNOW_DEPTH` as
   snow-covered-area height everywhere.
4. Regenerate the Bow annual trajectory with the audited restart fields and namelist controls,
   then rerun the unchanged 207-cell gate. Do not add a BTRAN exception without independent
   scientific and numerical review.
