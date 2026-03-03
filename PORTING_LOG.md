# CLM → Julia Porting Log

## Status

| Module | Tier | Status | Attempts | Notes |
|--------|------|--------|----------|-------|
| Gridcell type definition | 2 | ✓ PASSED | 1 | — |
| Landunit type definition | 2 | ✓ PASSED | 2 | — |
| Column type definition | 2 | ✓ PASSED | 2 | — |
| Patch type definition | 2 | ✓ PASSED | 2 | — |
| Domain decomposition | 2 | ✓ PASSED | 2 | — |
| Filter to mask conversion | 2 | ✓ PASSED | 2 | — |
| Temperature state | 3 | ✓ PASSED | 2 | — |
| Energy flux state | 3 | ✓ PASSED | 2 | — |
| Soil state | 3 | ✓ PASSED | 2 | — |
| Soil hydrology state | 3 | ✓ PASSED | 2 | — |
| Canopy state | 3 | ✓ PASSED | 2 | — |
| Lake state | 3 | ✓ PASSED | 2 | — |
| Albedo state | 3 | ✓ PASSED | 2 | — |
| Solar absorption | 3 | ✓ PASSED | 2 | — |
| Urban parameters | 3 | ✓ PASSED | 2 | — |
| Water info base | 3 | ✓ PASSED | 2 | — |
| Water state | 3 | ✓ PASSED | 2 | — |
| Water flux | 3 | ✓ PASSED | 2 | — |
| Water diagnostic bulk | 3 | ✓ PASSED | 2 | — |
| Water state bulk | 3 | ✓ PASSED | 2 | — |
| Water flux bulk | 3 | ✓ PASSED | 2 | — |
| Water balance | 3 | ✓ PASSED | 2 | — |
| Master water container | 3 | ✓ PASSED | 2 | — |
| Friction velocity | 3 | ✓ PASSED | 2 | — |
| Veg CN state | 4 | ✓ PASSED | 1 | — |
| Veg carbon state | 4 | ✓ PASSED | 1 | — |
| Veg nitrogen state | 4 | ✓ PASSED | 2 | — |
| Veg carbon flux | 4 | ✓ PASSED | 2 | — |
| Veg nitrogen flux | 4 | ✓ PASSED | 1 | — |
| Soil C state | 4 | ✓ PASSED | 1 | — |
| Soil N state | 4 | ✓ PASSED | 2 | — |
| Soil C flux | 4 | ✓ PASSED | 2 | — |
| Soil N flux | 4 | ✓ PASSED | 2 | — |
| Soil BGC state | 4 | ✓ PASSED | 2 | — |
| Crop type | 4 | ✓ PASSED | 2 | — |
| CN shared params | 4 | ✓ PASSED | 2 | — |
| Day length | 5 | ✓ PASSED | 2 | — |
| Surface albedo | 5 | ✓ PASSED | 2 | — |
| Urban albedo | 5 | ✓ PASSED | 2 | — |
| Surface radiation | 5 | ✓ PASSED | 2 | — |
| Urban radiation | 5 | ✓ PASSED | 2 | — |
| Snow SNICAR | 5 | ✓ PASSED | 2 | — |
| Aerosol | 5 | ✓ PASSED | 2 | — |
| Surface humidity | 5 | ✓ PASSED | 2 | — |
| Surface resistance | 5 | ✓ PASSED | 2 | — |
| Soil moisture stress | 5 | ✓ PASSED | 2 | — |
| Soil temperature | 6 | ✓ PASSED | 2 | — |
| Lake temperature | 6 | ✓ PASSED | 2 | — |
| SWRC base | 6 | ✓ PASSED | 2 | — |
| SWRC Clapp-Hornberg | 6 | ✓ PASSED | 2 | — |
| SWRC Van Genuchten | 6 | ✓ PASSED | 2 | — |
| Soil water movement | 6 | ✓ PASSED | 2 | — |
| Soil hydrology | 6 | ✓ PASSED | 2 | — |
| Snow hydrology | 6 | ✓ PASSED | 2 | — |
| Hydrology no drain | 6 | ✓ PASSED | 2 | — |
| Hydrology drainage | 6 | ✓ PASSED | 2 | — |
| Lake hydrology | 6 | ✓ PASSED | 2 | — |
| Hillslope hydrology | 6 | ✓ PASSED | 2 | — |
| Photosynthesis | 7 | ✓ PASSED | 2 | — |
| Canopy fluxes | 7 | ✓ PASSED | 2 | — |
| Canopy hydrology | 7 | ✓ PASSED | 2 | — |
| Bare ground fluxes | 7 | ✓ PASSED | 2 | — |
| Soil fluxes | 7 | ✓ PASSED | 2 | — |
| Urban fluxes | 7 | ✓ PASSED | 2 | — |
| LUNA | 7 | ✓ PASSED | 2 | — |
| Ozone damage | 7 | ✓ PASSED | 2 | — |
| Phenology | 8 | ✓ PASSED | 2 | — |
| CN allocation | 8 | ✓ PASSED | 2 | — |
| Growth respiration | 8 | ✓ PASSED | 2 | — |
| Maintenance respiration | 8 | ✓ PASSED | 2 | — |
| FUN model | 8 | ✓ PASSED | 2 | — |
| Gap mortality | 8 | ✓ PASSED | 2 | — |
| Decomp BGC | 8 | ✓ PASSED | 2 | — |
| Decomp MIMICS | 8 | ✓ PASSED | 2 | — |
| Decomposition | 8 | ✓ PASSED | 2 | — |
| Nitrif-Denitrif | 8 | ✓ PASSED | 2 | — |
| N leaching | 8 | ✓ PASSED | 2 | — |
| N dynamics | 8 | ✓ PASSED | 2 | — |
| C state update 1 | 8 | ✓ PASSED | 2 | — |
| C state update 2 | 8 | ✓ PASSED | 2 | — |
| C state update 3 | 8 | ✓ PASSED | 2 | — |
| N state update 1 | 8 | ✓ PASSED | 2 | — |
| N state update 2 | 8 | ✓ PASSED | 2 | — |
| N state update 3 | 8 | ✓ PASSED | 2 | — |
| CN driver | 8 | ✓ PASSED | 2 | — |
| Vegetation facade | 8 | ✓ PASSED | 2 | — |
| Fire base | 9 | ✓ PASSED | 2 | — |
| Fire Li2014 | 9 | ✓ PASSED | 2 | — |
| Methane | 9 | ✓ PASSED | 2 | — |
| VOC emissions | 9 | ✓ PASSED | 2 | — |
| Dust emission | 9 | ✓ PASSED | 2 | — |
| Irrigation | 9 | ✓ PASSED | 2 | — |
