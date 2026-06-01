# CLM.jl GPU Kernelization — Remaining Work Inventory & Plan

Status as of 2026-06-01. Goal: full `clm_drv!` timestep on GPU + AD-on-GPU.

**Done (whole-function GPU, Metal-validated at parity):** `soil_temperature!`,
`soil_water!` (ZD09 + moisture_form), `lake_temperature!`, `canopy_fluxes!`.
Plus ~156 element-wise ops, both batched linear solvers, and the
friction-velocity / photosynthesis-Ci-core kernels.

Established patterns (reuse these): per-(c,j) element-wise kernels with
eltype-converted literals + String→Int config; grouped device-view
`@adapt_structure` structs for field-heavy loops (alias struct fields to locals
so the body stays verbatim; **Metal caps total kernel args ~31** — bundle scalar
constants into one isbits struct, give Float64-kwarg fields their own type param
for AD); `_scatter_add!` (atomic on GPU, `+=` on Dual/CPU) for patch→column;
per-column kernels with internal sequential level/substep loops for loop-carried
deps; mask/`active`-based iteration instead of host filter compaction.

Per-function gate: full suite 15681/15681 + `--check-bounds=yes` + a
`scripts/gpu_validate_<fn>_e2e.jl` Metal harness at 0.0/round-off parity.

---

## Phase A — Biogeophys driver completion ("the other half" of clm_drv!)

Ordered by value × reuse-of-existing-patterns (do top-down):

| # | Target | Status now | Pattern / difficulty | Notes |
|---|--------|-----------|----------------------|-------|
| A1 | **photosynthesis! Ci solver** | PARTIAL (Ci core done; ~27 host loops outer) | per-patch + per-canopy-level loops; Ci Brent/hybrid core already GPU-callable | Completes the canopy chain (item #2). Self-contained. **Start here.** |
| A2 | **bareground_fluxes!** | HOST | stability iteration — SAME shape as canopy_fluxes Newton (active-mask rework + eltype literals + device-view structs) | Highest pattern reuse from the canopy work just landed. |
| A3 | **lake_fluxes!** | HOST | per-column Newton-Raphson | mirrors lake_temperature approach. |
| A4 | **soil_fluxes!** | HOST | per-col/per-patch energy partitioning | feeds soil_temperature outputs. |
| A5 | **surface_radiation! + canopy_sun_shade_fracs!** | HOST | per-patch/per-band/per-layer, branchy | radiation pair. |
| A6 | **canopy_interception_and_throughfall!** | PARTIAL | finish patch→column scatter (`_scatter_add!`) | mostly done. |
| A7 | **Snow suite**: handle_new_snow!, build_snow_filter!, bulk_flux_snow_percolation!, update_state_snow_percolation!, calc_and_apply_aerosol_fluxes!, post_percolation_adjust_layer_thicknesses!, sum_flux_add_snow_percolation!, snow capping (4 fns), snow_compaction! | HOST | variable-`snl` snow-layer loops, loop-carried across layers, per-species aerosol scatter | hardest mechanical cluster (combine/divide/compaction). |
| A8 | **Hydrology/runoff**: set_soil_water_fractions!, set_floodc!, saturated_excess_runoff!, set_qflx_inputs!, infiltration_excess_runoff!, route_infiltration_excess!, infiltration!, total_surface_runoff!, compute_effec_rootfrac_and_vert_tran_sink!, water_table!, renew_condensation! | mixed (water_table! is the most complex, 12+ nested loops) | per-col + sequential water-table search | water_table! is the long pole here. |
| A9 | **Urban path**: urban_radiation!, urban_fluxes! | HOST | per-landunit/per-patch urban canyon | defer if urban not a near-term priority. |
| A10 | **Smaller/init**: biogeophys_pre_flux_calcs! (partial), calc_root_moist_stress! (partial), calc_soilevap_resis!, calc_ozone_stress!, calculate_surface_humidity!, set_actual_roughness_lengths!, dust_dry_dep!, alt_calc!, update_daylength!, clm_drv_init! | mostly HOST | small per-col/per-patch loops | quick wins; finish the partials first. |

## Phase B — BGC / `use_cn` path (deferred biogeochemistry)

**Scope: ~562 host loops across 42 files (~25.5k LOC); only ~98 kernel refs in
11 files → ~83% still host.** Far larger than Phase A. Sub-order by module:

| Module | Files | ~Loops | Kernelized | Priority/notes |
|--------|-------|--------|-----------|----------------|
| Carbon/Nitrogen state updates | 8 | ~100 | NONE | core cascade; do first (element-wise, high reuse) |
| Decomposition & soil cycling (bgc/mimics/competition) | 6 | ~200 | PARTIAL ~20% | biggest; multiple methods |
| Methane (CH4) | 1 | ~80 | MINIMAL <5% | loop-heavy, low coverage |
| Phenology | 2 | ~25 | PARTIAL ~16% | crop/seasonal |
| Fire & disturbance | 2 | ~31 | PARTIAL ~25% | patch→col redistribution |
| Nitrogen cycling (n_dynamics/nitrif_denitrif) | 3 | ~20 | PARTIAL ~60% | most progress already |
| Other (allocation, litter transport, mortality, isotopes) | ~20 | ~106 | MINIMAL | tail |

## Phase C — Reverse-mode AD on GPU (item #3) — HARDWARE-GATED

- Forward-mode (ForwardDiff) already verified on Metal through kernelized ops.
- Reverse-mode (Enzyme) works on CPU; does NOT run on Apple Metal (no Apple-GPU
  Enzyme support → hangs compiling the adjoint). Needs **CUDA/AMD hardware**.
- Also blocked on CPU full-driver by a Julia-1.12 codegen bug ("instruction does
  not dominate all uses") → revisit on a Julia LTS with the fix.
- **This session: write unvalidated scaffolding** (`scripts/gpu_ad_reverse_validate.jl`)
  ready to run when a CUDA box + fixed Julia are available. Cannot be validated here.

---

### Rough effort
Phase A: ~15-20 driver functions, each a focused session (A2/A3 fast via canopy
reuse; A7/A8 hardest). Phase B: many sessions (it's ~half the total loop volume).
Phase C: blocked — scaffolding only until hardware.

### Recommended immediate sequence
A1 (photosynthesis, completes canopy) → A2 (bareground, max pattern reuse) →
A3 (lake_fluxes) → then continue down Phase A. Phase C scaffolding can be written
anytime (no validation feedback). Phase B after Phase A biogeophys is GPU-complete.
