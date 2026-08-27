# Upstream data-terms review — 2026-08-27

Status: **four provider terms support reuse with conditions; two exact-source terms remain
unresolved; no held fixture is promoted**.

This is a provenance assessment, not legal advice. Only primary provider pages were used.
The machine-readable record is `repro/manifests/upstream_terms.json`.

## Findings

| Source | Finding | Release consequence |
|---|---|---|
| Copernicus DEM | GLO-30/GLO-90 are described as freely licensed for distribution and reuse, with product-specific notices for unmodified and adapted products. | Potentially clear after confirming which DEM instance SYMFLUENCE downloaded and packaging the exact required notice. |
| NASA MODIS | NASA-led Earthdata are unrestricted or CC0 unless specifically marked otherwise; citation is strongly requested. | The underlying MCD12Q1 product is reusable, subject to exact-product confirmation and citation. |
| OpenGeoHub Zenodo mosaic | Record 8367523 identifies a mosaic derived from MCD12Q1 v061, but its displayed Rights/License field is blank. | Do not infer that NASA's terms cover the mosaic producer's contribution; request clarification or bypass it with a direct NASA acquisition recipe. |
| ECCC RDRS | ECCC's data-server licence permits copying, modification, publication, adaptation, and distribution, including commercial use, with attribution and third-party-right exceptions. | RDRS is conditionally clear when the exact accessed product is confirmed and the ECCC notice is packaged. |
| Water Survey of Canada | The official hydrometric dataset record states Open Government Licence—Canada. | WSC data are conditionally clear with required attribution. |
| CTSM input data | CTSM documents public automatic download of standard runtime inputs, but no file-level licence was found for exact S03. The source-code licence is not evidence for every input dataset. | Keep S03 and derived products on hold; seek NCAR confirmation or publish acquisition/checksum instructions without repackaging it. |

Primary records reviewed:

- <https://dataspace.copernicus.eu/explore-data/data-collections/copernicus-contributing-missions/collections-description/COP-DEM>
- <https://www.earthdata.nasa.gov/engage/open-data-services-software/data-use-policy>
- <https://zenodo.org/records/8367523>
- <https://eccc-msc.github.io/open-data/licence/readme_en/>
- <https://open.canada.ca/data/en/dataset/46763060-e859-4812-8da5-2361d99b4c34>
- <https://escomp.github.io/CTSM/lilac/obtaining-building-and-running/setting-ctsm-runtime-options.html>

## Recommended release route

1. Replace the OpenGeoHub mosaic dependency with direct NASA MCD12Q1 acquisition if the
   transformation can be reproduced without changing scientific inputs; otherwise obtain
   written licence clarification for Zenodo record 8367523.
2. Ask NCAR/CTSM maintainers for the redistribution status of exact S03. If it is not
   affirmatively redistributable, exclude D03 and provide the standard CTSM downloader
   path plus checksum.
3. D04's exact extraction recipe is recovered. Package ECCC/WSC attributions, the
   applicable Copernicus notice, NASA citations, and the complete transformation record
   only after T03 and T06 are resolved or bypassed.

Until these actions are complete, all six NetCDF artifacts remain `HOLD` and the GMD
redistribution checklist remains open.
