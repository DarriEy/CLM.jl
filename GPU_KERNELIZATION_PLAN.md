# CLM.jl GPU Kernelization — Remaining Work Inventory & Plan

Status as of 2026-06-01. Goal: full `clm_drv!` timestep on GPU + AD-on-GPU.

**Done (whole-function GPU, Metal-validated at parity):** `soil_temperature!`,
`soil_water!` (ZD09 + moisture_form), `lake_temperature!`, `canopy_fluxes!`,
`photosynthesis!` (A1), `bareground_fluxes!` (A2), `lake_fluxes!` (A3),
`soil_fluxes!` (A4), `surface_radiation!` + `canopy_sun_shade_fracs!` (A5),
`canopy_interception_and_throughfall!` (A6), and the **A8 hydrology cluster**:
`water_table!`, `set_soil_water_fractions!`, `set_floodc!`, `set_qflx_inputs!`,
`route_infiltration_excess!`, `infiltration!`, `total_surface_runoff!`,
`renew_condensation!`, `saturated_excess_runoff!`, `infiltration_excess_runoff!`,
`compute_effec_rootfrac_and_vert_tran_sink!`, and the **A7 snow suite**:
`combine_snow_layers!`, `divide_snow_layers!`, `snow_compaction!`, the
snow-percolation chain (`bulk_flux_snow_percolation!`, `update_state_snow_percolation!`,
`calc_and_apply_aerosol_fluxes!`, `post_percolation_adjust_layer_thicknesses!`,
`sum_flux_add_snow_percolation!`), `build_snow_filter!`, the 5 snow-capping fns,
and `handle_new_snow!`'s inline loops.
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
| A1 | **photosynthesis! Ci solver** | ✅ DONE (Metal whole-fn, 14 fields 0.0) | per-patch + per-canopy-level loops; Ci Brent/hybrid core already GPU-callable | Completes the canopy chain (item #2). Self-contained. |
| A2 | **bareground_fluxes!** | ✅ DONE (Metal whole-fn, 24/25 fields 0.0, kbm1 6.6e-08 F32) | stability iteration — SAME shape as canopy_fluxes Newton (active-mask rework + eltype literals + device-view structs) | Also fixed a latent AD bug: grnd_ch4_cond_patch needed its own Vg type param (Float64 kwarg vs Dual state). |
| A3 | **lake_fluxes!** | ✅ DONE (Metal whole-fn, 19 outputs ≤3.7e-05) | per-column Newton-Raphson | mirrors lake_temperature approach; LakeFluxDV device-view struct. |
| A4 | **soil_fluxes!** | ✅ DONE (Metal whole-fn, 17 fields exact 0.0) | per-col/per-patch energy partitioning | feeds soil_temperature outputs; 3 device-view bundles + _scatter_add! errsoi average. |
| A5 | **surface_radiation! + canopy_sun_shade_fracs!** | ✅ DONE (Metal whole-fn, 43 fields exact 0.0) | per-patch/per-band/per-layer, branchy | radiation pair; 5 device-view bundles, both snow/no-snow + urban/veg/bare + near-local-noon branches exercised. |
| A6 | **canopy_interception_and_throughfall!** | ✅ DONE (Metal whole-fn, 15 fields, 12 at 0.0) | finished patch→column scatter (`_scatter_add!`) | scatter exercised by a 2-patch→1-col harness. |
| A7 | **Snow suite**: handle_new_snow!, build_snow_filter!, bulk_flux_snow_percolation!, update_state_snow_percolation!, calc_and_apply_aerosol_fluxes!, post_percolation_adjust_layer_thicknesses!, sum_flux_add_snow_percolation!, snow capping (5 fns), snow_compaction!, combine/divide_snow_layers! | ✅ DONE (all whole-fn on Metal at 0.0; re-validated post-merge across snl=0/-1/-2/-3) | variable-`snl` snow-layer loops, loop-carried across layers, per-species aerosol scatter | hardest mechanical cluster — done via per-column kernels w/ in-thread sequential layer restructuring (combine/divide mutate snl in padded arrays); 4 agents, disjoint line-regions. **Follow-up:** handle_new_snow! orchestrator still calls host loops in snow_cover_fraction.jl (bulkdiag_new_snow_diagnostics!, update_snow_depth_and_frac!) — finish in A10. |
| A8 | **Hydrology/runoff**: set_soil_water_fractions!, set_floodc!, saturated_excess_runoff!, set_qflx_inputs!, infiltration_excess_runoff!, route_infiltration_excess!, infiltration!, total_surface_runoff!, compute_effec_rootfrac_and_vert_tran_sink!, water_table!, renew_condensation! | ✅ DONE (all 11 whole-fn on Metal at 0.0/round-off; re-validated post-merge) | per-col + sequential water-table search | water_table! long pole done: 6 host loops (3 sequential nlevsoi searches w/ loop-carried state) → one per-column kernel, 2 device-view bundles. **Tail also DONE** (separate batch): drainage! + clm_vic_map!, perched_water_table! + theta_based_water_table!, perched/subsurface_lateral_flow!, update_urban_ponding! + irrigation withdrawals — all whole-fn on Metal at 0.0/round-off, validated post-merge. Still HOST: hillslope-connectivity path of lateral flow (kept on CPU; device asserts non-hillslope). |
| A9 | **Urban path**: urban_radiation!, urban_fluxes! | ✅ DONE (both whole-fn on Metal; re-validated post-merge) | per-landunit/per-patch urban canyon | urban_fluxes! iterative canyon-air solve done as ONE thread per urban landunit running the full 3-iter loop-carried stability solve in-thread + member-patch/column gather over contiguous coli:colf/patchi:patchf ranges — each landunit owns its accumulation, so NO cross-thread reduction needed. Day worst 1.3e-7, low-wind worst 1.8e-5 (F32). |
| A10 | **Smaller/init**: biogeophys_pre_flux_calcs!, calc_root_moist_stress!, calc_soilevap_resis!, calc_ozone_stress!, calculate_surface_humidity!, set_actual_roughness_lengths!, dust_dry_dep!, alt_calc!, update_daylength!, clm_drv_init! | ✅ DONE (all whole-fn on Metal at 0.0/round-off; re-validated post-merge via 11 e2e harnesses) | small per-col/per-patch loops | also kernelized the handle_new_snow! tail (bulkdiag_new_snow_diagnostics! + update_snow_depth_and_frac! + add_newsnow_to_intsnow! in snow_cover_fraction.jl), closing A7's follow-up. |

## Phase B — BGC / `use_cn` path (deferred biogeochemistry)

**Scope: ~562 host loops across 42 files (~25.5k LOC).** Far larger than Phase A.
**Module 1 of 7 (C/N state-update cascade, ~100 loops) is now COMPLETE** — every
state-update fn runs whole-function on Metal at parity. Remaining ~460 loops across
decomposition / methane / phenology / fire / N-cycling / tail. Sub-order by module:

| Module | Files | ~Loops | Kernelized | Priority/notes |
|--------|-------|--------|-----------|----------------|
| Carbon/Nitrogen state updates | 8 | ~100 | ✅ **DONE** — entire C/N state-update cascade whole-fn on Metal at 0.0/round-off parity | core cascade complete (element-wise, high reuse). Pattern-setter `c_state_update1!` done by hand; the rest done in one parallel workflow batch (7 disjoint-file agents + per-file review, cherry-picked onto main): `c_state_update0!`/`c_state_update_dyn_patch!`, `c_state_update2!/2h!/2g!`, `c_state_update3!`, `n_state_update1!`/`n_state_update_dyn_patch!`, `n_state_update2!/2h!/2g!`, `n_state_update3!`/`n_state_update_leaching!`, `cn_precision_control!` + `truncate_*`. Every fn: per-col/per-patch kernels (in-thread j/k/l loops, own-index — NO scatter), device-view `@adapt_structure` bundles for field-heavy patch loops, eltype-converted literals, `MetalF32` Float32-down-convert harness. Each gated: full suite 15681/15681 + `--check-bounds` + its own `scripts/gpu_validate_*_e2e.jl` Metal harness. **BGC-harness gotchas that ONLY the Metal harness caught (CPU-green missed both):** (1) crit-threshold scalars passed Float64 → invalid Metal IR; convert to `eltype(state)` at the launch site. (2) a host-built Int index vector (`filter_truncatep`) moved to device via `adapt(typeof(float_ref), …)` coerces indices to Float32 → can't index device arrays; use `similar(ref, Int, n)`+`copyto!` to preserve Int eltype. |
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
~~A1 → A2 → A3 → A4 → A5 → A6 → A7 → A8~~ DONE (all merged to main at 15681/15681).
Phase A COMPLETE — the entire biogeophys/hydrology driver (A1–A10) runs whole-function
on Metal at 0.0/round-off parity, full suite 15681/15681 byte-identical, every function
gated on its own gpu_validate_*_e2e.jl harness. Last piece (urban_fluxes! iterative
canyon-air solve) landed. Remaining HOST in the driver: only the hillslope-connectivity
path of lateral flow (intentionally CPU; device asserts non-hillslope). Next: Phase B
(BGC / use_cn) — the larger remaining body of work. The soil_hydrology.jl tail
(drainage! + clm_vic_map!, perched/subsurface lateral flow, perched/theta water-table
variants, irrigation withdrawals + urban ponding) is now DONE — all whole-fn on Metal
at 0.0/round-off, validated post-merge. Phase C scaffolding can be written anytime
(no validation feedback). Phase B after Phase A biogeophys is GPU-complete.
