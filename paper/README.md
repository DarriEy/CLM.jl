# Fresh GMD manuscript workspace

No historical manuscript text or figure is an input to this directory. This file defines
the manuscript architecture; prose and headline numbers must wait for the frozen definitive
campaign.

## Working identity

- Article type: model description paper
- Working title: **CLM.jl v0.1.0: verification-gated translation of CLM5
  biogeophysics into differentiable, accelerator-capable Julia**
- Primary boundary: configured single-site/single-column biogeophysics
- Evidence authority: `repro/results/index.json` after it is marked definitive
- Claim authority: `docs/gmd/claims.csv`

“Differentiable” and “accelerator-capable” remain working-title terms, not validated
performance claims. Remove or narrow them before freeze if C03-C07 do not pass their
prespecified experiments.

## Planned manuscript

1. **Abstract** — problem, scoped implementation, frozen methods, only passed headline
   results, limitations, and availability.
2. **Introduction** — motivation for a Julia translation, distinction between translation
   fidelity and Earth-system validation, and explicit contribution list.
3. **Source model and scope** — exact CTSM version, retained processes, omitted processes,
   configuration boundary, and terminology.
4. **Translation architecture** — state representation, process ordering, numerical
   reformulations, optional differentiation and accelerator design, and traceability to
   upstream code.
5. **Verification protocol** — input equivalence, process oracles, annual experiments,
   metrics, prespecified thresholds, conservation checks, hardware, and environments.
6. **Results** — CPU process/annual fidelity first; smoothing/AD/calibration/accelerator
   subsections only for claims that pass the frozen campaign. Failed cells remain visible.
7. **Discussion** — interpretation of residuals, portability and maintainability,
   limitations of the site/configuration matrix, and distinction from observational model
   validation.
8. **Conclusions** — claims no broader than the passed ledger.
9. **Code and data availability** — exact software and experiment archive citations,
   development URL, licences, restrictions, and reviewer access where applicable.
10. **Author contributions, competing interests, acknowledgements, financial support, and
    AI-use disclosure** — completed and approved by the authors.

Appendices should carry technical material that belongs with the article, including the
process crosswalk, complete threshold table, supplementary pass/fail matrices, and
environment details. Large code, inputs, outputs, and ordinary analysis products belong in
the persistent archives rather than a journal supplement.

## Initial figure and table plan

Every item is provisional until its source runs exist in the definitive results index.

| Item | Purpose | Required claim/evidence |
|---|---|---|
| Fig. 1 | Scope and process-order architecture | implementation/source crosswalk |
| Fig. 2 | Verification workflow and evidence chain | E1-E4 provenance |
| Fig. 3 | Process-boundary Fortran–Julia discrepancies | C01 / E2 |
| Fig. 4 | Annual state/flux agreement across frozen sites | C01-C02 / E3-E4 |
| Table 1 | Exact implementations, environments, sites, inputs | E1 provenance |
| Table 2 | Prespecified metrics, thresholds, and verdicts | frozen protocol |
| Table 3 | Complete retained-claim results and limitations | definitive index |

Optional figures for smoothing, derivatives, calibration, correctness, or performance may
be added only if C03-C07 pass. Qualification plots are not manuscript figures.

## Drafting guardrails

- Do not type a numerical result directly from a log or qualification note.
- Do not convert a qualification pass into a paper result.
- Do not hide failed, missing, non-finite, or constant-series cells in aggregation.
- Use “bit-identical”, “within numerical tolerance”, and “scientifically acceptable” only
  with their distinct frozen definitions.
- Do not describe excluded inputs as reproducible until the acquisition workflow and exact
  checksums are independently exercised.
- Do not create the abstract, short summary, conclusions, or graphical abstract before the
  retained claim set is frozen.

