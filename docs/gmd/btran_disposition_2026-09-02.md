# BTRAN residual disposition — 2026-09-02 (author decision)

Status: **C01 rescoped by author decision; no gate exception added; independence of
review deferred to GMD interactive open peer review and disclosed as such.**

## The constraint this records

The qualification records required "independent scientific and numerical review" before
dispositioning the Bow BTRAN residual. **No independent pre-submission reviewers are
available to this project.** This memo records the author-level decision taken in that
situation and the safeguards that replace pretended independence:

1. **Nothing is exempted.** `DOCUMENTED_EXCEPTIONS` is untouched; `parity_gate.py`
   still reports the Bow BTRAN cell as an undocumented strict failure and still exits
   nonzero on `--check`. No threshold was changed at any point.
2. **The claim is rescoped to what the evidence shows**, rather than the evidence being
   promoted to the original claim (see below).
3. **The residual is presented, not buried.** The manuscript must show the failing cell
   in the main results (scorecard figure and text), with its full characterization.
4. **Independence comes from the journal.** GMD uses interactive public peer review;
   the referees and open discussion are the independent scientific and numerical review
   this residual receives. The cover letter and manuscript must state explicitly that
   pre-submission review was internal only and direct referee attention to this cell.

## The rescoped C01 claim

Old (never achieved): strict agreement in all 207 cells of the frozen three-site matrix.

New: *the scoped Julia CPU implementation agrees with the declared Fortran CLM5
reference in **204 of 207** strict cells of the frozen three-site matrix; the three
non-conforming cells are fully characterized — two are snow-covered-area
diagnostic-definition effects whose gridcell-mean counterparts (`SNOWDP`) pass, and one
is a stable 1.03e-5 internal-diagnostic residual (`BTRAN`, a mid-solve stress factor)
whose consuming fluxes all pass strict.*

Supporting characterization (all previously recorded):

- per-step oracle: tree patch at 9.8e-13; grass patch at a **stable** 1.0272e-5
  absolute offset, invariant across macOS release, Linux release, and Linux `-O0`
  (span 5e-16), amplified from ~1e-8 canopy-energy arithmetic
  (`repro/results/btran-sensitivity/matrix.json`, `repro/configs/bow_btran_oracle.toml`);
- the audited annual rerun with corrected restart/namelist controls reproduces the cell
  unchanged (`bow_audited_rerun_2026-08-31.md`) — restart, namelist, and compiler
  sensitivity are all excluded as causes;
- consumers pass strict: QVEGT and FCTR −0.674%, latent heat +0.574%
  (`qualification_2026-08-25.md`);
- the five-axis jitter screen and the daily-minimum amplification analysis
  (`docs/PARITY_FAITHFULNESS_AUDIT.md`) find no translation-bug fingerprint.

What the characterization does **not** establish — and the manuscript must say so —
is instantaneous implementation parity of the grass-patch BTRAN itself: a stable 1e-5
formulation-level difference (most plausibly order-of-operations in the coupled
canopy-energy/PHS solve) remains possible and unresolved.

## Ledger and verifier consequence

- `claims.csv` C01 → `PASS-QUALIFICATION-SCOPED`, evidence this memo; the strict
  204/207 gate output remains the recorded measurement.
- `verify_results_index.py`'s C01 rule changes from "must be FAIL" to "must be the
  SCOPED status and this memo must exist" — updated together with this record, per the
  script's own contract. The strict-parity gate's `claim_supported` stays **false**:
  204/207 cannot and does not support the unscoped parity claim.
- The manuscript may not use the word "parity" for C01 without the 204/207
  qualification attached.
