# ==========================================================================
# Productionized clm_drv! reverse-mode phases for the compositional engine.
#
# `compositional_reverse!` (canopy_fluxes_reverse.jl) reverses an ordered list of
# `(phase_fn, const_args)` entries, one Enzyme.autodiff call each, on a checkpointed
# bundle. The canopy energy-balance/photosynthesis solve is exposed there as
# `cf_rev_phases`. This module exposes the subsequent clm_drv! timestep phases that
# reverse-differentiate cleanly — soil_temperature!, soil_water! (Zeng-Decker 2009),
# water_table! — as the SAME kind of entries, operating on a unified bundle
# `(; inst, scratch)` so a multi-module / whole-timestep reverse is one call:
#
#     b   = driver_rev_bundle(inst)
#     ph  = driver_rev_phases(bounds, filt, config; canopy_aux = caux, n_canopy = 12)
#     db  = compositional_reverse!(ph, b, (db,b) -> (db.inst.<field> .= 2 .* b.inst.<field>))
#
# Each wrapper mutates `b.inst` in place and calls the exact production physics
# function; the per-phase aux (filters/bounds/config flags/configs) is passed as a
# Const NamedTuple — NEVER closure-captured (Enzyme analyzes a closure's captured
# environment for activity and trips on nested structs/BitVectors).
#
# Validated (Julia 1.10 LTS + Enzyme, real Bow cold-start inst) by
# scripts/enzyme_driver_reverse_full.jl (canopy+soil_temp+soil_water) and
# scripts/enzyme_driver_reverse_hydro.jl (+ water_table!). The forward wrappers are
# guarded WITHOUT Enzyme in test/test_driver_reverse.jl (runs in the 1.12 suite).
# ==========================================================================

# Unified bundle: the differentiated inst + the canopy scratch arrays (so the canopy
# block and the hydrology phases share one checkpointed object).
driver_rev_bundle(inst) =
    (; inst, scratch = cf_rev_scratch(Float64, length(inst.patch.column)))

# --------------------------------------------------------------------------
# Canopy block — runs the N canopy sub-phases (cf_rev_phases) as ONE chain entry
# on the bundle's cf view. Reversed coarsely (one Enzyme call over the block) in the
# chain; for the tight per-sub-phase canopy gradient use canopy_rev_gradient!.
# --------------------------------------------------------------------------
function canopy_rev_block!(b, aux, n::Int)
    cv = cf_rev_bundle(b.inst.canopystate, b.inst.energyflux, b.inst.frictionvel,
        b.inst.temperature, b.inst.solarabs, b.inst.soilstate, b.inst.water.waterfluxbulk_inst,
        b.inst.water.waterstatebulk_inst, b.inst.water.waterdiagnosticbulk_inst,
        b.inst.photosyns, b.scratch)
    for (f, cargs) in cf_rev_phases(aux, n)
        f(cv, cargs...)
    end
    return nothing
end

# --------------------------------------------------------------------------
# soil_temperature!
# --------------------------------------------------------------------------
soiltemp_rev_aux(bounds, filt; dtime = 1800.0) = (;
    urbantv = fill(323.15, bounds.endl), dtime = dtime,
    bc_col = bounds.begc:bounds.endc, bc_lun = bounds.begl:bounds.endl,
    bc_patch = bounds.begp:bounds.endp,
    nolakec = filt.nolakec, nolakep = filt.nolakep, urbanl = filt.urbanl, urbanc = filt.urbanc)

function soiltemp_rev_phase!(b, aux)
    i = b.inst
    soil_temperature!(i.column, i.landunit, i.patch, i.temperature, i.energyflux,
        i.soilstate, i.water.waterstatebulk_inst, i.water.waterdiagnosticbulk_inst,
        i.water.waterfluxbulk_inst, i.solarabs, i.canopystate, i.urbanparams,
        aux.urbantv, i.atm2lnd.forc_lwrad_downscaled_col, aux.nolakec, aux.nolakep,
        aux.urbanl, aux.urbanc, aux.bc_col, aux.bc_lun, aux.bc_patch, aux.dtime)
    return nothing
end

# --------------------------------------------------------------------------
# soil_water! (Zeng-Decker 2009 default; the retention curve is an empty singleton
# and the movement config a small mutable struct — both pass cleanly as Const).
# --------------------------------------------------------------------------
function soilwater_rev_aux(filt, config; dtime = 1800.0)
    cfg = config.use_aquifer_layer ? SoilWaterMovementConfig() :
        SoilWaterMovementConfig(soilwater_movement_method = MOISTURE_FORM,
                                lower_boundary_condition = BC_ZERO_FLUX)
    return (; hydrologyc = filt.hydrologyc, urbanc = filt.urbanc,
              swrc = SoilWaterRetentionCurveClappHornberg1978(), cfg = cfg, dtime = dtime)
end

function soilwater_rev_phase!(b, aux)
    i = b.inst
    soil_water!(i.column, aux.hydrologyc, aux.urbanc,
        i.soilhydrology, i.soilstate, i.water.waterfluxbulk_inst, i.water.waterstatebulk_inst,
        i.temperature, i.canopystate, i.energyflux, aux.swrc, aux.cfg, aux.dtime)
    return nothing
end

# --------------------------------------------------------------------------
# water_table! (aquifer / ThetaBasedWaterTable). The single Bool kwarg is resolved
# inside the wrapper (primal-time), so Enzyme differentiates the body without a
# kwarg-NamedTuple thunk.
# --------------------------------------------------------------------------
watertable_rev_aux(bounds, filt, config; dtime = 1800.0) = (;
    hydrologyc = filt.hydrologyc, bc_col = bounds.begc:bounds.endc,
    nlevsoi = varpar.nlevsoi, dtime = dtime,
    recompute_frost_table = config.use_aquifer_layer)

function watertable_rev_phase!(b, aux)
    i = b.inst
    water_table!(i.soilhydrology, i.soilstate,
        i.temperature.t_soisno_col, i.water.waterstatebulk_inst.ws, i.water.waterfluxbulk_inst,
        i.column.dz, i.column.z, i.column.zi,
        aux.hydrologyc, aux.bc_col, aux.nlevsoi, aux.dtime;
        recompute_frost_table = aux.recompute_frost_table)
    return nothing
end

# --------------------------------------------------------------------------
# hydrology_no_drainage! — the namesake HydrologyNoDrainage end-of-step
# diagnostics (recomputes h2osoi_vol, the 10 cm / total water columns, etc. from
# the updated h2osoi_liq/ice). Runs after water_table! (the intervening snow
# capping/compaction/layer combine-divide ops are discrete and no-ops at snl=0).
# --------------------------------------------------------------------------
hydnodrain_rev_aux(bounds, filt; dtime = 1800.0) = (;
    nolakec = filt.nolakec, hydrologyc = filt.hydrologyc, urbanc = filt.urbanc,
    snowc = filt.snowc, nosnowc = filt.nosnowc, bc_col = bounds.begc:bounds.endc,
    dtime = dtime, nlevsno = varpar.nlevsno, nlevsoi = varpar.nlevsoi,
    nlevgrnd = varpar.nlevgrnd, nlevurb = varpar.nlevurb)

function hydnodrain_rev_phase!(b, aux)
    i = b.inst
    hydrology_no_drainage!(i.temperature, i.soilstate, i.water.waterstatebulk_inst,
        i.water.waterdiagnosticbulk_inst, i.column, i.landunit,
        aux.nolakec, aux.hydrologyc, aux.urbanc, aux.snowc, aux.nosnowc,
        aux.bc_col, aux.dtime, aux.nlevsno, aux.nlevsoi, aux.nlevgrnd, aux.nlevurb)
    return nothing
end

# --------------------------------------------------------------------------
# Pre-soil_water! SURFACE-HYDROLOGY block (HydrologyNoDrainage). The full ordered
# sequence that runs AFTER soil_temperature! and BEFORE soil_water!, partitioning the
# surface water that becomes soil_water!'s top boundary:
#   set_soil_water_fractions! → set_floodc! → saturated_excess_runoff! → set_qflx_inputs!
#   → infiltration_excess_runoff! → route_infiltration_excess! → update_h2osfc!
#   → infiltration! → total_surface_runoff! → compute_effec_rootfrac_and_vert_tran_sink!
# All share one Const aux. saturated_excess_runoff!'s Bool kwargs + update_h2osfc!'s
# dtime/h2osfcflag + compute_effec_rootfrac's use_hydrstress default inside the wrappers;
# qflx_floodg is a zero gridcell vector (no river flood input), captured as Const.
function surfhydro_rev_aux(bounds, filt; dtime = 1800.0)
    return (; hydrologyc = filt.hydrologyc, nolakec = filt.nolakec, urbanc = filt.urbanc,
              bc_col = bounds.begc:bounds.endc, nlevsoi = varpar.nlevsoi, nlevsno = varpar.nlevsno,
              dtime = dtime, qflx_floodg = zeros(Float64, bounds.endg))
end

function setsoilfrac_rev_phase!(b, aux)
    i = b.inst
    set_soil_water_fractions!(i.soilhydrology, i.soilstate, i.water.waterstatebulk_inst.ws,
        i.column.dz, aux.hydrologyc, aux.bc_col, aux.nlevsoi, aux.nlevsno)
    return nothing
end
function setfloodc_rev_phase!(b, aux)
    i = b.inst
    set_floodc!(i.water.waterfluxbulk_inst.wf.qflx_floodc_col, aux.qflx_floodg,
        i.column.gridcell, i.column.itype, aux.nolakec, aux.bc_col)
    return nothing
end
function satexcess_rev_phase!(b, aux)
    i = b.inst
    saturated_excess_runoff!(i.sat_excess_runoff, aux.hydrologyc, aux.bc_col,
        i.column, i.landunit, i.soilhydrology, i.soilstate, i.water.waterfluxbulk_inst)
    return nothing
end
function setqflx_rev_phase!(b, aux)
    i = b.inst
    set_qflx_inputs!(i.water.waterfluxbulk_inst, i.water.waterdiagnosticbulk_inst,
        i.column.snl, aux.hydrologyc, aux.bc_col)
    return nothing
end
function inflexcess_rev_phase!(b, aux)
    i = b.inst
    infiltration_excess_runoff!(i.infilt_excess_runoff, i.soilhydrology, i.soilstate,
        i.sat_excess_runoff.fsat_col, i.water.waterfluxbulk_inst, i.water.waterdiagnosticbulk_inst,
        aux.hydrologyc, aux.bc_col)
    return nothing
end
function routeinfl_rev_phase!(b, aux)
    i = b.inst
    route_infiltration_excess!(i.water.waterfluxbulk_inst, i.soilhydrology,
        i.column.landunit, i.landunit.itype, aux.hydrologyc, aux.bc_col)
    return nothing
end
function updateh2osfc_rev_phase!(b, aux)
    i = b.inst
    update_h2osfc!(i.column, i.soilhydrology, i.energyflux, i.water.waterfluxbulk_inst,
        i.water.waterstatebulk_inst, i.water.waterdiagnosticbulk_inst,
        aux.hydrologyc, aux.bc_col; dtime = aux.dtime)
    return nothing
end
function infil_rev_phase!(b, aux)
    infiltration!(b.inst.water.waterfluxbulk_inst, aux.hydrologyc, aux.bc_col)
    return nothing
end
function totalrunoff_rev_phase!(b, aux)
    i = b.inst
    total_surface_runoff!(i.water.waterfluxbulk_inst, i.soilhydrology, i.water.waterstatebulk_inst.ws,
        i.column.snl, i.column.itype, i.column.landunit, i.landunit.urbpoi,
        aux.hydrologyc, aux.urbanc, aux.bc_col, aux.dtime)
    return nothing
end
function rootsink_rev_phase!(b, aux)
    i = b.inst
    compute_effec_rootfrac_and_vert_tran_sink!(aux.bc_col, aux.nlevsoi, aux.hydrologyc,
        i.soilstate, i.canopystate, i.water.waterfluxbulk_inst, i.energyflux,
        i.column, i.landunit, i.patch)
    return nothing
end

# Ordered pre-soil_water surface-hydrology phase entries.
function surface_hydrology_rev_phases(bounds, filt; dtime = 1800.0)
    a = surfhydro_rev_aux(bounds, filt; dtime)
    return Any[(setsoilfrac_rev_phase!, (a,)), (setfloodc_rev_phase!, (a,)),
               (satexcess_rev_phase!, (a,)),   (setqflx_rev_phase!, (a,)),
               (inflexcess_rev_phase!, (a,)),  (routeinfl_rev_phase!, (a,)),
               (updateh2osfc_rev_phase!, (a,)), (infil_rev_phase!, (a,)),
               (totalrunoff_rev_phase!, (a,)),  (rootsink_rev_phase!, (a,))]
end

# --------------------------------------------------------------------------
# Assembler: the ordered (phase_fn, const_args) list for a clm_drv! reverse, in
# forward order: [canopy] → soil_temp → <surface hydrology block> → soil_water →
# water_table → hydrology_no_drainage. `canopy_aux === nothing` skips the canopy block.
# `include_surface=false` reverts to the soil_temp→soil_water direct jump.
# --------------------------------------------------------------------------
function driver_rev_phases(bounds, filt, config; canopy_aux = nothing,
                           n_canopy::Int = 12, dtime = 1800.0, include_surface::Bool = true)
    phases = Any[]
    canopy_aux === nothing || push!(phases, (canopy_rev_block!, (canopy_aux, n_canopy)))
    push!(phases, (soiltemp_rev_phase!, (soiltemp_rev_aux(bounds, filt; dtime),)))
    include_surface && append!(phases, surface_hydrology_rev_phases(bounds, filt; dtime))
    push!(phases, (soilwater_rev_phase!,  (soilwater_rev_aux(filt, config; dtime),)))
    push!(phases, (watertable_rev_phase!, (watertable_rev_aux(bounds, filt, config; dtime),)))
    push!(phases, (hydnodrain_rev_phase!, (hydnodrain_rev_aux(bounds, filt; dtime),)))
    return phases
end
