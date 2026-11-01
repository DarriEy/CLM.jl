# Single-point LUH2 land-use data (Bow at Banff)

These two files drive the FATES land-use tests (`test/test_fates_landuse.jl`) and the
LUH2 Fortran-parity harness (`scripts/validation/fates_fortran_parity/`) offline,
without the 339 MB source file.

| File | Derived from |
|---|---|
| `LUH2_1pt_Bow_hist_simyr1850-2015.nc` | CESM inputdata `lnd/clm2/surfdata_map/fates-sci.1.68.3_api.31.0.0_tools.1.0.1/LUH2_states_transitions_management.timeseries_4x5_hist_simyr1850-2015_c231101.nc` |
| `fates_landuse_pft_map_1pt_Bow.nc` | CESM inputdata `lnd/clm2/surfdata_map/fates-sci.1.77.0_api.36.0.0/fates_landuse_pft_map_4x5_240206.nc` |

Both source files come from the CESM input data repository,
<https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/>.

## Changes made

Each file is a single-grid-cell subset of the 4°×5° CESM file. It is the cell nearest
the Bow at Banff site (50°N, 245°E), written on a 1×1 grid with the same variables. No
values were modified. The subsets were produced by
[`scripts/validation/fates_fortran_parity/make_luh_1pt.jl`](../../../scripts/validation/fates_fortran_parity/make_luh_1pt.jl).
To regenerate them, run:

```bash
julia --project=. scripts/validation/fates_fortran_parity/make_luh_1pt.jl \
    <LUH2_states_transitions_management...4x5...c231101.nc> \
    <fates_landuse_pft_map_4x5_240206.nc> <Bow surfdata_clm.nc> data/fates/luh_1pt
```

## Source data and attribution

The LUH2 states, transitions and harvest fields are CESM's 4°×5° regridding of the
**Land-Use Harmonization 2 (LUH2)** dataset for CMIP6. That dataset was produced by the
University of Maryland and distributed through input4MIPs as
`input4MIPs.CMIP6.CMIP.UofMD`
([doi:10.22033/ESGF/input4MIPs.10454](https://doi.org/10.22033/ESGF/input4MIPs.10454)).

Please cite:

> Hurtt, G. C., Chini, L., Sahajpal, R., Frolking, S., Bodirsky, B. L., Calvin, K.,
> et al. (2020). Harmonization of global land use change and management for the period
> 850–2100 (LUH2) for CMIP6. *Geoscientific Model Development*, 13, 5425–5464.
> <https://doi.org/10.5194/gmd-13-5425-2020>

Use of CMIP6 forcing data is also subject to the
[CMIP6 terms of use](https://pcmdi.llnl.gov/CMIP6/TermsOfUse).

The landuse × PFT map is a FATES/CTSM input file prepared by the NCAR/FATES
development teams for the CESM input data repository.

## License

The input4MIPs LUH2 dataset is licensed under the
[Creative Commons Attribution-ShareAlike 4.0 International License (CC BY-SA 4.0)](https://creativecommons.org/licenses/by-sa/4.0/).
**These derived files are distributed under the same CC BY-SA 4.0 license.** They are
not covered by the MIT License that applies to the CLM.jl source code.
