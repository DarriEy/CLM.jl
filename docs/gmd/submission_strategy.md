# GMD submission strategy

Status: **working strategy — submission is not yet authorized**  
Policy checked: 2026-08-27

This strategy starts from the current evidence workspace. Historical drafts and figures
are not manuscript inputs.

## Recommended article and claim boundary

Target a GMD **model description paper** centred on the scoped CLM.jl v0.1.0 translation
and its verification methodology. The defensible core is narrower than the repository:

- single-site/single-column CLM5 biogeophysics on Julia CPU Float64;
- an exact, named CTSM source version and explicitly matched configurations;
- process-boundary and annual verification with conservation checks active;
- differentiability and accelerator results only if their dedicated gates pass before
  protocol freeze.

Do not claim general CLM5 equivalence, production/global readiness, or observational
skill. Seasonal-cycle differentiability (C03), joint parameter recovery (C05), Metal
accelerator correctness (C06), and the measured scaling crossover (C07) reached
qualification passes on 2026-08-28/31 (see `claims.csv` and the requalification
records), so they are candidates for retention — but like every claim they enter the
manuscript only if they pass the frozen definitive campaign, and the accelerator claim
is scoped to the tested M-series Metal devices (CUDA stays backend-capable design
only).

## Current journal requirements that shape the package

The current GMD instructions require a review manuscript PDF in one-column portrait form
with page and line numbers, embedded fonts, and figures and tables placed near first
mention. The upload also needs separately entered abstract text, a non-technical short
summary of at most 500 characters including spaces, funding information, and any justified
supplement. Figure files must use Arabic numbering, normally be 300 dpi, combine panels,
and use colour-vision-accessible schemes.

GMD's code and data policy requires a final section entitled **Code and data availability**
before acknowledgements. It must cite persistent public archives of the precise code and
data versions, identify the licence and generic development location, and provide the
configuration, boundary conditions, inputs, run control, preprocessing, postprocessing,
and plotting material used for reported results. Access must exist when the preprint is
submitted; an acceptance embargo is not permitted. GitHub is a development location, not
a frozen archive. Restrictions outside the authors' control must be stated, while
confidential editor/reviewer access must still be arranged without compromising reviewer
anonymity.

The submission instructions also require disclosure in Methods or Acknowledgements when AI
tools generated part of the manuscript. The corresponding author must review the final
wording and make that disclosure accurate before submission.

Official sources:

- https://www.geoscientific-model-development.net/submission.html
- https://www.geoscientific-model-development.net/about/manuscript_types.html
- https://www.geoscientific-model-development.net/policies/code_and_data_policy.html

## Release sequence

### Stage A — narrow and freeze

1. Obtain domain-science and numerical review of the experiment matrix, metrics, and the
   stable Bow grass BTRAN residual.
2. Decide whether E5-E9 are repaired and retained or explicitly removed from v0.1.0 paper
   claims. Exclusion is preferable to an indefinitely expanding campaign.
3. Freeze the protocol, threshold table, site membership, release-candidate commit, and
   Linux environment. Replace every `TBD` in the frozen configuration.
4. Tag the exact software release only after the freeze commit passes the source, licence,
   and clean-room gates.

### Stage B — definitive campaign

1. Run the complete Linux campaign from a clean checkout using one fail-closed command.
2. Preserve every pass, failure, environment record, input hash, raw-output hash, and
   analysis version in the results index.
3. Generate all tables and figures exclusively from that index and archived raw outputs.
4. Independently reproduce the documented workflow on a second clean environment.

### Stage C — archive before writing headline results

Create two linked, immutable deposits:

- **software deposit:** exact tagged source, licences, documentation, environment lock,
  run/analysis/plot scripts, and the source-archive manifest;
- **experiment deposit:** policy-cleared inputs, configurations, raw outputs, result index,
  checksums, and generated figure/table data.

Mint version-specific persistent identifiers. If restricted inputs remain essential,
document the restriction and create a reviewer-access route before registration. The
current source-only archive is a candidate build, not the required final deposit.

### Stage D — manuscript and internal review

Write into the architecture in `paper/README.md`. Every quantitative statement must cite a
claim ID and resolve through the definitive results index. Run four sign-offs:

1. scientific scope and threshold review;
2. independent numerical/reproducibility review;
3. code/data/licence and privacy review;
4. author, funding, contribution, conflict, acknowledgement, and AI-disclosure review.

Only after these sign-offs should the Copernicus LaTeX template be vendored from the
current official distribution and the review PDF compiled.

### Stage E — submission

Upload the review PDF, abstract, short summary, and any justified supplement; register the
model-description manuscript in the most appropriate GMD subject area. Cite both archive
deposits in the text and bibliography. Retain the exact submitted PDF and archive versions
as a submission record. Any scientific change during review requires a new archive version
and re-generated results, figures, and manuscript values.

## Go/no-go decision

Submission is **NO-GO** today. The minimum transition to GO is:

- protocol and numerical thresholds frozen after recorded human review;
- primary CPU equivalence claim passes the definitive campaign;
- all retained claims pass and all failed/TBD claims are removed from title, abstract,
  conclusions, figures, and tables;
- inputs and outputs are publicly archived or have policy-compliant reviewer access;
- clean-room reproduction passes;
- exact software and experiment deposits have persistent identifiers;
- manuscript assets pass journal-format, licence, accessibility, and authorship checks.

