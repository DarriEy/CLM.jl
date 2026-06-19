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
# Assembler: the ordered (phase_fn, const_args) list for a clm_drv! reverse.
# `canopy_aux === nothing` skips the canopy block (e.g. a hydrology-only reverse).
# Append further phases (hydrology diagnostics, fluxes, BGC) here as they are
# validated — each is one more entry, same template.
# --------------------------------------------------------------------------
function driver_rev_phases(bounds, filt, config; canopy_aux = nothing,
                           n_canopy::Int = 12, dtime = 1800.0)
    phases = Any[]
    canopy_aux === nothing || push!(phases, (canopy_rev_block!, (canopy_aux, n_canopy)))
    push!(phases, (soiltemp_rev_phase!,   (soiltemp_rev_aux(bounds, filt; dtime),)))
    push!(phases, (soilwater_rev_phase!,  (soilwater_rev_aux(filt, config; dtime),)))
    push!(phases, (watertable_rev_phase!, (watertable_rev_aux(bounds, filt, config; dtime),)))
    push!(phases, (hydnodrain_rev_phase!, (hydnodrain_rev_aux(bounds, filt; dtime),)))
    return phases
end
