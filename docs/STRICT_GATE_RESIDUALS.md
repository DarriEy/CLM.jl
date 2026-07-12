# Strict-gate residual audit

This file records the evidence required before a strict-gate miss is classified
as a documented exception. An exception does not make the underlying series
numerically identical; it states why the selected diagnostic crosses a hard
threshold under two otherwise-close trajectories.

## Shortwave penetration at snow melt-out

`SABG_PEN` is discontinuous when the final snow layer disappears. The remaining
misses are concentrated at that boundary rather than distributed through the
snow season:

| Site | Dominant date | Julia `SABG_PEN` | Fortran `SABG_PEN` | Julia SWE | Fortran SWE |
|---|---:|---:|---:|---:|---:|
| Kherlen | 2017-05-07 | 62.325 W m⁻² | 0 | 0.039 mm | 0 |
| Hubbard Brook | 2017-04-16 | 0 | 28.955 W m⁻² | 4.270 mm | 7.200 mm |
| Massa Aletsch | 2013-07-09 | 85.126 W m⁻² | 0 | 12.585 mm | 0 |

Kherlen's single May 7 value is larger than its full annual-mean bias after
division by 365. Hubbard Brook similarly gets nearly its full mean bias from
April 16. Massa has several delayed melt-out events, including July 9–10, so its
glacier residual is broader but has the same mechanism. These cells are
documented exceptions; persistent absorbed/reflected shortwave fields remain
subject to the strict gate.

## Bow BTRAN

The annual miss is concentrated from 2003-04-29 through early May. On April 30,
Julia/Fortran daily values are:

| Quantity | Julia | Fortran |
|---|---:|---:|
| `BTRAN` | 0.3364 | 0.2141 |
| transpiration (`QVEGT`) | 6.358e-6 | 6.834e-6 |
| photosynthesis (`FPSN`) | 1.1479 | 1.3060 |
| total soil liquid | 260.080 mm | 259.717 mm |
| SWE | 199.885 mm | 201.069 mm |
| exposed LAI | 2.02062 | 2.02097 |

Thus the BTRAN residual is not explained by materially wetter Julia soil:
transpiration is lower despite the higher diagnostic.

The spring oracle was extended with `vegwp`, complete LUNA 10-day state,
`vcmx25_z` / `jmx25_z`, `SMP` / `HK`, and soil/root conductances. A full CTSM
restart at 2003-04-29 23:00 supplied the remaining prior-step state. Replaying
n11616 from that shared state is finite and gives:

| Patch | Julia BTRAN | Fortran BTRAN | Difference |
|---|---:|---:|---:|
| tree | 0.206950 | 0.206950 | +9.009e-7 |
| grass | 0.681003 | 0.681029 | −2.613e-5 |

This directly clears the instantaneous PHS/BTRAN implementation at the annual
worst episode. The annual +4.3% and daily nRMSE 0.21 are accumulated upstream
trajectory/history differences, so Bow `BTRAN` is a documented exception.
Snow, temperature, and soil-water state remain independently gated: this
classification does not exempt their observed one-step differences.

## Yakutia RH2M

The annual RH difference is confined primarily to the very cold months
(November–February: about +1.6 to +3.0 percentage points). The underlying
diagnostics remain close:

| Quantity | Annual RMSE |
|---|---:|
| `TSA` | 0.0083 K |
| `Q2M` | 8.2e-6 kg kg⁻¹ |
| `RH2M` | 1.47 percentage points |

At roughly −40 °C, saturation specific humidity is so small that the accepted
absolute Q2M residual is magnified in `Q2M / QSAT(TSA)`, and Julia reaches the
100% cap while Fortran commonly reports 95–96%. Both implementations use the
same Flatau ice polynomial and hard 100% cap. This derived cold-saturation
diagnostic is documented; its underlying temperature and humidity remain
strictly gated.

## Hubbard Brook SNOWLIQ

The annual +3.74% residual is localized to the same final spring melt episode
as Hubbard Brook `SABG_PEN`. April 12–14 accounts for 84% of total absolute
SNOWLIQ error; April 14 alone contributes 57% (Julia 11.258 mm, Fortran
5.179 mm). Values outside that transition are close. This is documented as a
melt/refreeze timing exception, while SWE, snow ice, snowmelt, and runoff remain
independently gated.

## Massa Aletsch soil ice

Soil ice is exactly zero in both trajectories through 2013-11-25. The divergence
begins at the November 26 freeze/snow boundary:

| Date | Julia soil ice | Fortran soil ice | Julia 10-cm T | Fortran 10-cm T |
|---|---:|---:|---:|---:|
| Nov 25 | 0 | 0 | 276.228 K | 276.190 K |
| Nov 26 | 0 | 1.199 mm | 276.142 K | 274.479 K |
| Nov 27 | 1.367 mm | 12.239 mm | 274.422 K | 272.718 K |
| Dec 20 | 57.532 mm | 61.554 mm | 272.323 K | 272.326 K |

After the event, temperature reconverges while the liquid/ice partition retains
the latent-heat memory. Total soil water stays close: on November 27 it is
751.586 mm in Julia and 752.230 mm in Fortran. The relative annual ice error is
large because ice exists only during this short late-year interval. This
glacier freeze-onset phase diagnostic is documented; total water, soil
temperature, and energy fluxes remain independently gated.
