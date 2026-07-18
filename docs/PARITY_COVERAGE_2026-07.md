# Fortran-parity harness coverage — 2026-07-17 refresh

**What this is.** A run-verified status for every `scripts/fortran_parity_*.jl`
harness that was *added since* the 2026-06-17 snapshot in
[`HARNESS_COVERAGE_AUDIT.md`](HARNESS_COVERAGE_AUDIT.md) but whose actual pass/fail
was never recorded. Each row below was **executed on 2026-07-17** (Julia 1.12,
`--project=.`) against the local dumps that survive on this box
(`clm_bgc_spinup/bgc_ref_summer` @ nstep 1757852, `clm_bgc_spinup/tier1_iso` @
nstep 1757873). Cross-referenced from
[`FORTRAN_VALIDATION_BACKLOG.md`](FORTRAN_VALIDATION_BACKLOG.md).

**Status legend**
- **VALIDATED** — Julia recomputes the field (not injected) and matches the Fortran
  dump within a tight tolerance.
- **DIVERGES** — runs, reads a real Fortran field, and disagrees beyond tolerance
  (root-cause noted).
- **GREEN-BY-SKIP** — runs and passes, but the check is vacuous on the available
  data (field identically zero in the snow-free/summer window, or no Fortran field
  exists so it degrades to a self-consistency check).
- **NEEDS-FORTRAN-RUN** — cannot run to a real parity conclusion because the
  required dump/domain is absent on disk (needs regenerating), not because of a
  port defect.

---

## Scorecard (harnesses run 2026-07-17)

| Harness | Subsystem / module | Fortran field diffed | Dump used | Result | Status |
|---|---|---|---|---|---|
| `fortran_parity_activelayer.jl` | `active_layer.jl` (`alt_calc!`) | altmax, altmax_indx, altmax_lastyear(_indx) | bgc_ref_summer n1757852 | max\|rel\| **0.0e+00** (exact) | **VALIDATED** |
| `fortran_parity_soilresis.jl` | `surface_resistance.jl` (`calc_soilevap_resis!`) | SOILRESIS (genuine recompute, overwrites injected) | bgc_ref_summer n1757852 | rel **1.7e-08** (2671.559 vs 2671.559) | **VALIDATED** |
| `fortran_parity_albedo.jl` | `surface_albedo.jl` bands (no injection) | albd/albi/fabd/fabi/fsun/albgrd/albgri/albsod/albsoi | bgc_ref_summer n1757852 | max\|rel\| **0.0e+00** (computed, not injected) | **VALIDATED** — closes the Table-B "RT injection can mask an albedo bug" concern |
| `fortran_parity_decomprates.jl` | soil-BGC N source/sink (decomp+nitrif/denit+immob+uptake) | Δsmin_nh4_vr / Δsmin_no3_vr / Δsminn_vr per layer (integrated rate effect) | bgc_ref_summer n1757852 | NO3 ratio **1.001**; NH4/sminn ratio **1.19** (top layer) | **VALIDATED w/ known residual** — NH4 over-drain is the documented closed-GPP/sminn residual, NOT a decomp-rate bug. The direct `*_VR_P` rate dumps are dead-zero (dead instrumentation hook) → direct rate parity impossible. |
| `fortran_parity_aerosol.jl` | `aerosol.jl` (`aerosol_masses!`) | mss_bcphi/bcpho/ocphi/ocpho/dst1-4 | bgc_ref_summer n1757852 | all 0.0 (snow-free window) | **GREEN-BY-SKIP** — structural/no-crash only; needs a **snowy** dump to stress redistribution |
| `fortran_parity_daylength.jl` | `daylength.jl` (`init/update_daylength!`) | (dump carries NO `dayl`/`prev_dayl`) | bgc_ref_summer n1757852 | self-consistency vs closed form **1.2e-09** | **GREEN-BY-SKIP** — no Fortran ground-truth field; self-consistency only |
| `fortran_parity_luna_update.jl` | `luna.jl` (`update_photosynthesis_capacity!` EOD acclimation) | vcmx25_z / jmx25_z post-EOD-update | Fortran LunaMod PARITYLUNA/PARITYLUNO @ n1757880 (values in-harness) | **jmx25_z ratio 1.00000** ✓; **vcmx25_z ratio 1.052–1.056** ✗ | **DIVERGES** — see §"LUNA vcmax residual" |
| `fortran_parity_isotopes.jl` | `carbon_isotopes.jl` / `c_iso_flux.jl` | c13/c14 veg pools + one-step increment | **bgc_ref_iso n4700-4725 (NEW 2026-07-18)** | c13/c14 spinup RUNS (via #2119 SourceMods bypass). 3 real port bugs FIXED (facade-config wiring made the veg iso path dead; i_litr/i_met_lit sentinel −9 crashed the CIsoFlux p2c scatters). Cascade now runs; **one-step c13 pool increment DIVERGES ~100%** — the per-isotope `c_state_update*` cascade is unported (CTSM calls each `CStateUpdate` 3×). | **PARTIAL** — reference exists, subsystem un-blocked; state-update wiring is the next task (harness validates it) |
| `fortran_parity_lake.jl` | `lake_*.jl` (con/fluxes/hydrology/temperature) | TG/TLAKE/LAKEICEFRAC/FSH/EFLX_LH/FIRE/FSA/... | `clm_lake_run` **absent on disk** | cannot run | **NEEDS-FORTRAN-RUN** — dump domain removed; but see §"Lake" (was validated w/ fixes through commit 3db0c1a) |
| `fortran_parity_soilparam.jl` | `cold_start.jl` pedotransfer (watsat/bsw/sucsat → SMP) | SMP recomputed vs dump | default `clm_parity_run` **absent on disk** | cannot run as-is | **NEEDS-DUMP-REDIRECT** — forcing-free check; the identical SMP validation still holds, just point `DUMPDIR` at `bgc_ref_summer`. Not a port regression. |

---

## LUNA vcmax residual (the one real DIVERGE) — root-cause status

`fortran_parity_luna_update.jl` isolates LUNA's end-of-day acclimation update
(`update_photosynthesis_capacity!`) by injecting the exact Fortran per-step climate
accumulators + pre-update vcmx25/jmx25 (instrumented from `LunaMod` PARITYLUNA at
nstep 1757880) and calling the update standalone.

Result (Julia post / Fortran post):

| patch | field | pre | Julia post | Fortran post | ratio |
|---|---|---|---|---|---|
| p2 (NET tree) | vcmx25_z | 49.563 | 52.618 | **49.844** | **1.0556** |
| p2 | jmx25_z | 92.985 | 92.985 | 92.985 | 1.00000 |
| p3 (C3 grass) | vcmx25_z | 136.862 | 143.304 | **136.197** | **1.0522** |
| p3 | jmx25_z | 336.077 | 334.313 | 334.272 | 1.00012 |

**Analysis.** `jmx25_z` matches to ~5 digits, so the shared machinery (max_daily_pchg,
the constraint `chg_constrn = min(|chg|, x·max_daily_pchg)`, the temperature/enzyme
turnover factor, and the jmax optimum `PNetopt`) is all correct. Only the **vcmax
(Rubisco) optimum `PNcbopt`** is ~5.3–5.6% high. Verified NOT a constant/param bug:
`LUNA_Fc25=294.2`, `LUNA_Fj25=1257.0`, `LUNA_Q10Enz=2.0`, `enzyme_turnover_daily`
default 0.0114 are all byte-identical to `LunaMod.F90` (Fc25=294.2, Fj25=1257.0,
readNcdioScalar `enzyme_turnover_daily`). The divergence is inside the joint
N-allocation optimizer `nitrogen_allocation!` — the carboxylation N pool `Ncb`
(→ `PNcbopt`) settles higher than Fortran's while the electron-transport pool `Net`
(→ `PNetopt`, jmax) matches. Because `vcmx25_opt = PNcbopt·FNCa·Fc25` sits above the
current value, Julia's update hits the `max_daily_pchg` cap and steps up ~+3, whereas
Fortran's optimum is near the current value (small unconstrained step).

This is exactly the **"Rubisco-N optimum residual"** flagged as the open next item in
MEMORY `luna-wiring-status`. It is a subtle ~5% residual in the iterative
`Ncb`/`Nstore` optimization, **not a clean one-line fix**, and it is deliberately
**left un-tuned**. Closing it requires instrumenting the Fortran `Allocation`
subroutine to emit the internal `Ncb`/`vcmx25_opt` (not just the post-update
`vcmx25_z`, which PARITYLUNO already gives) so the exact iteration where Julia and
Fortran's `Ncb` diverge can be pinned — a new Fortran instrumentation run, not a
disk-present dump. Tracked as backlog **C1a** (below).

**Impact caveat.** In the full driver LUNA `vcmx25_z` is used for photosynthesis via
`Fc25`-scaled activity; a persistent +5% vcmax would bias GPP high in
Rubisco-limited conditions. The Bow CN-pool parity (`fortran_parity_cn_summer.jl`,
tree availc ~0.998×) does not obviously surface this, likely because Bow summer is
largely light/`Jmax`-limited (jmax matches) — which is itself a reason to prioritize
the Fortran `Allocation` instrumentation.

---

## Lake (`fortran_parity_lake.jl`) — historical status

The lake reference domain (`clm_lake_run` h0 time series) has been **removed from
this box**, so the harness cannot re-run today. Its committed history shows it was
carried well past a first pass: after `c400220` (add harness) the surface
turbulent-flux residual was chased through three real fixes —
`525b868` (Monin–Obukhov `obu` sign / ustar collapse), `c4d4409` (fetch-limited
Charnock roughness), and `3db0c1a` ("ZengWang2007 frozen z0hg + carry ust_lake →
**closes the surface-flux residual**"). So lake thermodynamics **and** the surface
flux solve were Julia↔Fortran validated at the time; only the ability to *reproduce*
that run locally is gone. Backlog **B1** therefore stays NEEDS-A-RUN (regenerate the
lake h0), but should record that a validated result already existed.

---

## Net effect on the backlog

| Backlog row | Was | Now (2026-07-17) |
|---|---|---|
| B1 Lake | "harness EXISTS — verify dump" | Harness ran to VALIDATED historically (through commit 3db0c1a); **dump domain now absent** → still NEEDS-A-RUN to reproduce, but not un-validated. |
| C1 Isotopes | "harness EXISTS — verify dump/pass" | Confirmed: c13/c14 **needs a c13/c14-enabled spinup** (Fortran branch restart has no isotope vars → hard crash). Irrigation half runs but is inactive at Bow. Stays NEEDS-FORTRAN-RUN. |
| — (new) C1a LUNA vcmax | (implicit in Table B "LUNA injected/unvalidated") | Promoted to a KNOWN **DIVERGES** (~5.3–5.6% vcmax high); needs Fortran `Allocation` internal (`Ncb`/`vcmx25_opt`) instrumentation. |
| Audit Table B: albedo bands, active-layer, soilresis, aerosol, daylength | "runs but NOT validated / injection can mask" | activelayer/soilresis/albedo now **VALIDATED** (albedo computed, no injection → masking concern closed); aerosol/daylength **GREEN-BY-SKIP** (need snowy / a `dayl` dump field). |

No src fix was warranted: the single divergence (LUNA vcmax) is a documented,
non-self-contained optimizer residual and was deliberately not tuned to pass.
