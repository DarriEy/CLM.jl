# AD-smoothing bias audit — the "conservation is not accuracy" sweep

**Status:** complete. 460 `smooth_*` call sites audited across `src/` (33 files).
**Material bias found & fixed:** 1 site-group (the methane aerenchyma/transpiration
flux floors — a 1.39e7× overshoot on a mol/m³/s axis). Everything else is either
already axis-scaled, saturated (far from its kink), or on an axis where `log(2)/k`
is a physically negligible fraction of the signal. Numbers below.

This document is the durable ranked record. The executable regression guard is
`test/test_smooth_axis_guard.jl` (error law + axis-scaled-constant registry +
"hardened guards stay hard" pins + per-file census pin).

---

## 1. The thesis and the error law

CLM.jl replaces ~460 discontinuous `max/min/clamp/heaviside/abs/ifelse` guards with
smooth surrogates (`src/infrastructure/smooth_ad.jl`) so ForwardDiff/Enzyme get
continuous gradients. Each takes a sharpness `k` (default 50). The load-bearing fact:

```
smooth_max(0, x)  overshoots  max(0, x)  by  log(2)/k   — IN THE UNITS OF THE AXIS.
smooth_min(a, b)  at  a ≈ b   biases     min(a, b)  by  log(2)/k   likewise.
```

`k = 50` is calibrated for an **O(1) axis** (`log(2)/50 = 0.013863`). On any axis
where 0.0139-in-those-units is NOT negligible, the surrogate **fabricates physics**,
and a `smooth_max(0, ·)` / `smooth_min(·, cap)` ReLU sits **AT its kink in the common
case**, so the bias fires every step — it is not a rare-branch rounding.

Water and energy have fatal conservation checks that eventually caught the water-axis
instances (PRs #211/#214/#219). **Temperature (K), resistances/conductances (s/m,
mol/m²/s), carbon/nitrogen pools & fluxes (gC/m², mol/m³/s), fractions [0,1], and
albedo have NO conservation check — a bias there is INVISIBLE.** This audit targets
exactly those axes.

### The critical gating fact (scope of every finding below)

The surrogates are **type-/mode-gated** (`_use_smooth`, `SMOOTH_MODE`):

- `SMOOTH_MODE = :auto` (the default) + Float64 host **→ exact `max/min`**.
- Float32 GPU kernels **→ exact `max/min`** (`_is_ad_type(Float32) == false`).
- Bias of `log(2)/k` materializes **only** for `ForwardDiff.Dual` (AD gradients) and
  for Float64 under `SMOOTH_MODE = :always` (calibration FD/AD-consistency; Enzyme
  reverse Float64).

**So the shipped forward model and every published parity result are byte-identical
to hard `max/min` and are unaffected by everything in this document.** "Material bias"
here means material to the **AD gradient / `:always`-mode smoothed physics** — which
is precisely the model a gradient-based calibration or an adjoint would optimize
against, and which has no balance check to flag it.

---

## 2. Method

1. **Static inventory.** Every `smooth_*` call site classified by (axis units, `k`,
   `log(2)/k` in those units, at-kink-vs-saturated). Saturated = the two arguments are
   routinely `|k·(a−b)| > 36` apart, where the surrogate returns the exact value.
2. **Empirical divergence probe** (`scripts/smooth_divergence_probe.jl`, extended):
   warm up N steps in exact mode, snapshot, then take ONE exact step and ONE `:always`
   step **from the identical state**, and rank all 1265 state arrays by
   `max|Δ|/max|exact|`. Run at a winter state (warmup 120) and a summer state
   (warmup 7900, active photosynthesis) on the Bow single-point domain.
3. **Rank** by (at-kink frequency) × (0.0139-in-axis-units as a fraction of the signal)
   × (downstream reach), focused on the no-conservation-check axes.

---

## 3. Empirical result — biogeophysics single-step bias is SMALL (post prior fixes)

One `:always` step vs one exact step from an identical warmed-up state, top divergences:

| state field | axis | summer rel `max|Δ|/max|exact|` | note |
|---|---|---|---|
| `energyflux.cgrndl_patch` | ∂latent/∂T | 6.7e-2 | near-freezing patch only; tiny absolute (1.1e-7) |
| `energyflux.cgrnd_patch` / `dgnetdT` | ∂net/∂T | 2.5e-2 / 1.8e-2 | surface-Newton derivative terms |
| `energyflux.eflx_snomelt` / `temperature.xmf` | latent heat of fusion | 2.1e-3 | melt-season phase change |
| `surfalb.flx_abs*` | absorbed SW | 6e-4 | two-stream absorbed-flux |
| `photosyns.gs_mol*` | stomatal conductance | **1.7e-4** | see §5: gs selection uses HARD max by design |
| `photosyns.vcmax_z` / `lmr*` | µmol/m²/s | 3.5e-5 | negligible |

**Verdict:** the biogeophysics smoothing is well-controlled. The largest divergence
(`cgrndl`, 6.7%) is the **qsat freezing blend** — `es/desdT = smooth_ifelse(T−TFRZ,
water_poly, ice_poly)` (`qsat.jl:91–96`), a `k=50` sigmoid on a °C axis → a ±0.72 °C
transition width blending the water (`desdT` slope 0.444) and ice (0.503) branches that
genuinely meet at 0 °C. This is a **legitimate symmetric phase-transition smoothing**,
not a one-sided dimensional-`k` floor bias, and it is confined to patches within ~0.7 °C
of freezing. **Not material.** (The winter probe: all fields rel < 1e-8 except a lake
diagnostic artifact, §6.)

The probe runs `use_cn=false`, so it does not exercise BGC carbon/N accumulation — that
axis was covered by static analysis (§4–§5).

---

## 4. THE MATERIAL FIND — methane aerenchyma/transpiration flux floors (FIXED)

`src/biogeochem/methane.jl`, the **live** `ch4_aere!` kernel (`_ch4aere_patch_kernel!`).
CH4 and O2 are transported through plant aerenchyma and by transpiration; those fluxes
live on a **mol/m³/s axis of magnitude ~1e-9**. The kernel floored/capped them with
`smooth_max`/`smooth_min`:

| line | expression | axis | bias under `:always`/AD |
|---|---|---|---|
| 1382 | `tranloss = smooth_max(tranloss, 0)` | CH4 transpiration flux, mol/m³/s ~1e-9 | → 0.0139, **1.39e7×** the flux |
| 1415 | `aere = smooth_max(aere, 0)` | CH4 aerenchyma flux, mol/m³/s | → 0.0139, ~7 orders too big |
| 1424 | `oxaere = smooth_max(oxaere, 0)` | O2 aerenchyma flux, mol/m³/s | → 0.0139, ~7 orders too big |
| 1435 | `aeretran = smooth_min(aere+tranloss, conc/dt+prod)` | scattered CH4 flux, mol/m³/s | biases transported CH4 by 0.0139 |
| 1437 | `smooth_min(tranloss, aeretran)` | scattered CH4 transpiration flux | same |

Aerenchyma is the **dominant CH4 conduit** (per the module's own comments), and these
are the quantities scattered into `ch4_aere_depth` / `ch4_tran_depth` / `o2_aere_depth`
— the actual transported gas. Under AD / `:always`, `smooth_max(1e-9, 0)` returns
**0.013862944** instead of `1e-9` — a measured **1.39e7×** overshoot (verified at the
primitive). There is no CH4 conservation check to catch it.

**This was a porting inconsistency, not a design choice.** The module's own site-level
reference `site_ox_aere!` already uses **hard** `max` on the identical axis, with the
explicit comments `# HARD: mol/m3/s axis (~1e-9)` / `# HARD: ... 0.0139 is 7 orders too
big` (methane.jl:1238/1271/1280). The kernel — the live runtime path — had diverged to
`smooth_max`.

**Fix:** the five kernel guards reverted to hard `max`/`min`, matching `site_ox_aere!`
and Fortran `ch4Mod`. This is:
- **byte-identical** on the default forward path (Float64:auto and Float32 GPU already
  evaluate `smooth_*` as exact `max/min`);
- **correct** under AD / `:always` (no 0.0139 injection);
- **AD-differentiable** (`max`/`min` are ForwardDiff/Enzyme-differentiable; and for
  `tranloss`, floored at 0, the guard is a provable no-op — `tranloss ≥ 0` always).

Per the established remedy ([[smooth-k-is-dimensional]]): a floor against a constant
carries no derivative on the clamped branch, so a smooth surrogate buys nothing for AD
and only injects the bias. A hard floor is legitimate exactly here.

Guarded by new "hardened guards stay hard" pins + the tightened census pin
(methane.jl 23 → 18) in `test/test_smooth_axis_guard.jl`. `test_methane.jl` 146/146 and
`test_smooth_axis_guard.jl` 101/101 pass under `--check-bounds=yes`.

---

## 5. Other at-kink sites on invisible axes — audited, NOT material

Ranked by residual concern. All are AD/`:always`-only; none warrant a code change.

### Photosynthesis / LUNA (resistance, conductance, Ci, N axes)
- **The dangerous conductance selection is HARD by design.** The `quadratic_solve`
  root selections that produce `gs_mol` in **mol/m²/s** use plain `max` per an in-code
  CONTRACT (`photosynthesis.jl:581–594`, pinned in the guard test) — a 0.0139 mol/m²/s
  floor there would be a 13 900 µmol/m²/s floor on gs. The `smooth_max` guards that
  touch gs do so only **after** the `×1e6` rescale, i.e. on a µmol/m²/s axis where
  0.0139 is negligible. Empirically `gs_mol` diverges only 1.7e-4 (§3).
- **LUNA RH floors** (`luna.jl:332` `smooth_max(0.65, relh)`; `829`/`908`
  `smooth_min(1.0, rh)`; `516` RH-deficit `/0.35`): ~1.4–2.1% bias on an RH fraction
  that feeds the Ball-Berry conductance inside the LUNA acclimation optimum. Small, and
  LUNA is a slow (days-integrated) accumulator. Highest residual concern in this group.
- **`smooth_max(ci − Γ*, 0)` in Pa** (photosynthesis ~24 sites + luna ~10): during
  active photosynthesis `ci − Γ* ≈ 20 Pa` (0.07%, negligible); near the CO₂ compensation
  point (dawn/dusk/deep-shade) it sits at the kink and fabricates ~0.0139 Pa → a small
  positive assimilation on patches whose true rate is 0. Confined to the compensation
  regime; the CONTRACT comment quantifies the sibling at 0.008–0.024 µmol/m²/s.
- **`smooth_clamp(t10−TFRZ, 11, 35)` (°C, ~15 sites)**, VPD/vapor-pressure guards (Pa),
  PAR floors (µmol photon), resistance caps `smooth_min(·, 2e4)` (s/m): 0.0139 in those
  units is < 0.1% of the signal, or the guard is saturated. Negligible.

### Fire (Li2014/2016/2021/2024 — fractions, burned area)
- **Fuel/RH/moisture/btran rails** (`fb`, `afuel`, `arh`, `arh30`, btran-stress): each a
  dimensionless [0,1] factor sitting at a 0/1 rail in the common case, leaking ~1.4% into
  the burned-area product. Pervasive (~90 sites) but each is 1.4% of a bounded factor.
  Moderate, no conservation check — the largest *fraction-axis* residual.
- **`smooth_min(0.0016, dtrotr_rate)`** (`fire_li2021.jl:157`, `fire_li2024.jl:246`):
  the cap value (0.0016, an annual deforestation fraction) is ~9× **smaller** than the
  0.0139 bias, so at the kink the min undershoots by more than the cap. **Material in
  magnitude**, but narrow reach: transient-landcover + tropical-closed-forest only, and
  it feeds the deforestation-fire term, not a conserved pool. Flagged; left for a
  fire-specific pass (no hard reference sibling exists, unlike methane).
- **Burned-fraction caps** `smooth_min(1, baf_crop+baf_peatf)`, `smooth_min(1, spread²+…)`,
  `smooth_min(60, |lat|)`: argument is a per-step rate « 1 → saturated → **exact**. Clean.

### Phenology / decomp / veg-struct (carbon & nitrogen pools)
- **The classic `smooth_max(0, carbon_pool)`-written-back-to-a-pool pattern does NOT
  occur.** The genuinely dangerous per-step N-flux guards in `decomp_competition.jl`
  (`supplement_to_sminn_vr`, `sminn_to_plant_vr`) and the C/LAI floors in
  `veg_struct_update.jl` are **hard `max`/`min`**, not smoothed.
- **`phenology.jl:2021/2022`** `smooth_max(0, leafc_storage − leafc)·bgtr` — the only
  smooth calls on a real gC/m² carbon axis that can rest at their kink. Throttled by
  `bgtr ≈ 3.2e-8/s` → a phantom source of ~**0.014 gC/m²/yr**. Negligible.
- **`decomp_competition.jl:190`** `fpi = smooth_min(1, supply/demand)` — sits at the
  supply==demand N-limitation transition (at kink by construction); 0.0139 ≈ **1.4%** of
  the N-partition fraction `fpi`. Moderate (a fraction, not a pool); the highest-leverage
  N-side smooth site since it directly gates an immobilization flux.
- **GDD/FDD accumulators** (soil-T − TFRZ, °C·day): sit at the kink every shoulder-season
  step and steer leaf-out/leaf-off *timing* (an indirect carbon-burst trigger, not a
  per-step leak). ≤ 0.0139 °C·day/day. Watch-item.
- All height/LAI/SAI/daylength/porosity/[0,1]-fraction sites: not C/N axes, or saturated.

### Temperature / snow / albedo
- **`surface_albedo.jl:1059/1061/1119/1121`** `smooth_max(d_fab*, 0)` — absorbed-PAR
  **flux fraction [0,1]**; ReLU at kink for every optically-thin canopy layer, ~1.4% into
  per-layer sunlit/shaded PAR absorption feeding photosynthesis. **The largest
  fraction-axis residual in biogeophysics.** Confined to thin layers; AD/`:always`-only.
- **`snow_hydrology.jl:265`** `smooth_max(0, snomelt_accum − newsnow·1e-3)` — a
  **metres-of-water axis**: 0.0139 m = **13.9 mm** of phantom accumulated snowmelt at the
  kink. Same class as the fixed `SNOW_PERC_WATER_K` leak, but limited reach — it drives
  only z0m/snow-roughness aging, **not** the water balance. High magnitude, low blast
  radius. Candidate for an axis-scaled `k` if snow roughness is ever calibrated via AD.
- **Snow grain radius (µm ~54–1500), snow density (kg/m³), temperature (K)**: 0.0139 in
  those units is negligible (a persistent but tiny K bias). Snow-mass `kg/m²` phase-change
  guards are **already axis-scaled** (`PHASE_CHANGE_MASS_K = 1e9`, `SNOW_PERC_WATER_K =
  1e9`). Most albedo/optical-depth/coszen guards are saturated.

---

## 6. Illustrative artifact (not a live bug) — lake `betaprime`

The winter probe flagged `lakestate.betaprime_col` jumping 0 → 1386 under `:always`.
Cause: `lake_temperature.jl:428–429` `sabg_nir = smooth_min(sabg_nir, sabg)` at
low/zero solar (both args ≈ 0) returns `−log(2)/50 = −0.0139`, then divided by the
`max(1e-5, sabg)` guard → −1386. **Energy-invisible**: `betaprime` is only ever
multiplied by `sabg`, and where the garbage appears `sabg ≈ 0` annihilates it, so
`eflx_soil_grnd` is unaffected (errsoi rel 8e-3). A corrupted diagnostic, not a material
flux. A defensible `smooth_clamp(betaprime, 0, 1)` (physical fraction, a no-op in exact
physics) would remove it; left as documented given zero energy impact.

---

## 7. Bottom line

- **460 sites audited.** The shipped forward model (Float64:auto, Float32 GPU) is
  byte-identical to hard `max/min` throughout — none of this touches published results.
- **1 material bias fixed:** the methane aerenchyma/transpiration mol/m³/s flux floors
  (1.39e7× overshoot under AD/`:always`), by reverting the live kernel to the hard
  `max/min` its own `site_ox_aere!` reference already used.
- **Everything else is negligible or already handled:** the worst residuals are
  fraction-axis rails at ~1.4% (fire fuel/RH factors; absorbed-PAR `d_fab` floors; LUNA
  RH; `fpi`), all AD/`:always`-only, all on axes with no conservation check but a bias
  small relative to the signal. Two magnitude-large-but-narrow watch-items are noted:
  the fire `smooth_min(0.0016, dtrotr)` deforestation-rate cap and the snow-roughness
  `snomelt_accum` metres-of-water floor.
- Guarded by `test/test_smooth_axis_guard.jl` (error law + axis registry + hardened-guard
  pins + per-file census pin), which now pins the methane guards hard.

See also: [[conservation-is-not-accuracy]], [[smooth-k-is-dimensional]],
[[vacuous-checks-bug-class]].
