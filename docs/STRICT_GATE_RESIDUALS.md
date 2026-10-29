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
| Kherlen (resolved, see below) | 2017-05-07 | 62.325 W m⁻² | 0 | 0.039 mm | 0 |
| Hubbard Brook | 2017-04-16 | 0 | 28.955 W m⁻² | 4.270 mm | 7.200 mm |
| Massa Aletsch | 2013-07-09 | 85.126 W m⁻² | 0 | 12.585 mm | 0 |

### 2026-10-28 re-audit: purely melt-out timing — and Kherlen's timing was a bug

`SABG_PEN = sabg − sabg_lyr(top)` (SurfaceRadiationMod.F90:851) is non-zero only
while explicit snow layers exist, so a daily mean flips on whichever day the last
layer goes. All three sites were rerun on this branch and on its base
(`infra/env-data-roots`) against the Fortran h0:

| Site | cell | base: ann / nRMSE | this PR: ann / nRMSE |
|---|---|---:|---:|
| Kherlen | `SABG_PEN` | +1.84 % / 0.214 **FAIL** | **+0.001 % / 0.0004 PASS** |
| Kherlen | `SNOW_DEPTH` | +7.48 % / 0.252 FAIL | +2.01 % / 0.169 FAIL |
| Hubbard Brook | `SABG_PEN` | −1.59 % / 0.167 FAIL | −1.58 % / 0.167 FAIL |
| Massa | `SABG_PEN` | +1.61 % / 0.167 FAIL | +1.30 % / 0.165 FAIL |

**Kherlen — fixed.** The base run kept 0.036–0.039 mm of snow through May 6–7
(SABG_PEN 1.26 and 62.3 W m⁻²) while Fortran melted out on May 6. With the CLMNCEP
datm humidity cap (`QBOT ≤ qsat(datm_eSat(TBOT))`, added for the lake work and
verified at Bow to reproduce the h0 `QBOT` to 2.5e-8) Julia's May 6 SWE is
0.035 mm = Fortran's, and SABG_PEN is 0 on May 7 exactly as Fortran. The cold,
dry steppe forcing has supersaturated hours; uncapped, they deposited frost the
reference never sees. Kherlen `SABG_PEN` is removed from the exception list.

**Hubbard Brook — timing, inherent so far.** Both runs melt out on Apr 16, but
Julia's last layer goes a few hours earlier (daily-mean SWE 4.26 vs 7.20 mm,
SABG_PEN 0 vs 29.0). The trajectory difference is the melt/refreeze partition of
the preceding warm spell: `SNOWLIQ` Apr 12–14 is 7.41/17.39/11.26 mm (Julia) vs
5.96/15.94/5.18 (Fortran) — the same episode as the documented SNOWLIQ exception —
with melt 33.1 vs 31.9 mm d⁻¹ on Apr 12. Unchanged by this PR (−1.59 → −1.58 %).
A fixable cause, if any, is in the liquid-retention / refreeze step of that
episode, not in the radiation diagnostic.

**Massa — timing, partly improved.** Julia carries 8–10 mm more SWE through the
first week of July (Jul 7: 43.3 vs 35.4 mm), so the glacier snow melts out on
Jul 10 instead of Jul 9. The datm cap removed the Jul 10 tail (base 14.2 W m⁻²,
now 0.008) and moved the cell +1.61 → +1.30 %; the remaining SWE excess builds
over the melt season and is not a melt-out-day artefact.

Hubbard Brook and Massa `SABG_PEN` remain documented exceptions, as melt-out
timing differences of hours (HB) and one day (Massa).

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
worst episode.

### 2026-10-28 re-audit: what actually drives the gap (verdict: inherent)

The hypothesis that a BTRAN gap this size with near-identical soil water must be
a root-stress defect (frozen-layer exclusion, `smpso`/`smpsc`, `rootr`, the
`t_soisno > tfrz` test) **does not apply to this configuration**. Bow runs PHS
(`use_hydrstress = .true.`): `calc_root_moist_stress` still runs, but its `btran`
is overwritten in `PhotosynthesisHydraulicStress` (PhotosynthesisMod.F90:3778,
`btran = bsun·f_sun + bsha·f_sha` from the leaf water potentials). Under PHS the
frozen-layer / `smpso` / `smpsc` path never reaches `BTRANMN`. What does is the
vegetation water potential, which tracks the root-weighted soil matric potential.

To see the hours the daily file hides, the Fortran reference was re-run from the
Jan 1 restart with an hourly per-patch tape (`hist_dov2xy=.false.`: `VEGWP`,
`BTRANMN`, `TSOI`, `SOILLIQ`, `SOILICE`, `SMP`, …). The rerun reproduces the
reference h0 **bit-for-bit** (BTRANMN, H2OSNO, TOTSOILLIQ, FSH, TSOI_10CM:
max|Δ| = 0 over Jan 1 – Apr 30). Julia was probed hourly over the same window.

Weekly means, grass patch (35 % of the gridcell; `liq` = relative difference
J/F−1 per layer, `SMP` and root `vegwp` as J/F ratios):

| ISO week | liq L1 | liq L2 | liq L3 | SMP L1 | SMP L2 | SMP L3 | vegwp root |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1–10 (Jan–early Mar) | ≤ +0.003 | ≤ +0.001 | ≤ 0.0001 | 0.99–1.00 | 1.00 | 1.00 | 0.990 |
| 14 (Apr 1–6) | +0.006 | +0.002 | −0.000 | 0.930 | 0.984 | 1.002 | 1.001 |
| 16 (Apr 14–20) | +0.017 | +0.001 | −0.000 | 0.913 | 0.993 | 1.000 | 0.994 |
| 17 (Apr 21–27) | +0.008 | +0.004 | +0.001 | 0.955 | 0.981 | 0.993 | 0.956 |
| 18 (Apr 28–May 4) | +0.004 | **+0.012** | **+0.007** | 0.986 | **0.925** | **0.962** | **0.873** |

Soil temperature agrees to a few hundredths of a kelvin and ice to < 2 % of a
mm throughout. The mechanism is the soil-water retention curve: at Bow's spring
potentials (−0.4 to −5 MPa) `SMP ∝ s^(−B)` with B ≈ 10, so the **+1.2 %**
layer-2 liquid (from snowmelt infiltration in the third week of April) becomes a
**7.5 %** less negative SMP, a **13 %** less negative root water potential, and —
because PHS stress is `2^(−(ψ/ψ50)^ck)` around the curve's knee — a daily-minimum
grass BTRAN of 0.684 versus Fortran's 0.342 for the Apr 29 UTC day (the h0
record stamped Apr 30), read per patch from the hourly tape; the tree patch
agrees (0.162 vs 0.157). The grass split opens day by day as the shallow layers
diverge: 1.000/0.9995 (Apr 26), 0.997/0.985, 0.907/0.648, 0.684/0.342 (Apr 29). Hourly on Apr 29 at 18 UTC the grass leaf
potential is −265 282 mm (Julia) vs −307 070 mm (Fortran) with layer-1/2 liquid
3.79/5.73 vs 3.76/5.67 mm. The gridcell gap starts on Apr 26 (0.030 vs 0.018) and
persists into mid-May (May 8: 0.113 vs 0.080), i.e. it is not confined to a
single flip of the exposed-vegetation filter — `frac_veg_nosno` is 1 for both
patches at every probed hour.

Two smaller observations, recorded rather than resolved:

* Through the winter the grass root `vegwp` sits **1 % less negative** in Julia
  (ratio 0.990) while SMP agrees to 0.1 % in layers 1–8. The shared-state replay
  above matches to 3e-5, so this is the PHS solve's carried `vegwp` (the Newton
  starts from the previous step's potentials), not its formula; it does not
  produce the spring gap, which appears only once the shallow layers diverge.
* The `QBOT` datm cap added in this PR (it reproduces the Bow h0 daily `QBOT` to
  2.5e-8 vs 1.3e-3 rms without it) leaves the spring soil state unchanged to the
  printed precision and moves annual BTRAN from +4.249 % to +4.261 %.

**Verdict: inherent, not a port defect.** The calculation is correct given its
state (replay) and the state agrees to ~1 %; a ~1 % shallow-soil-liquid
trajectory difference is amplified ~10× by the retention curve and again by the
PHS stress function at the curve's knee. The daily-*minimum* diagnostic then
selects the most amplified hour. A fixable cause would have to be found in the
spring snowmelt-infiltration partition (week 16–18), not in the BTRAN path; Bow
`BTRAN` stays a documented exception. Snow, temperature, and soil-water state
remain independently gated and pass (TOTSOILLIQ +0.09 %, SOILWATER_10CM +0.22 %).

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
