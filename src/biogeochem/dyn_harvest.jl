# ==========================================================================
# Ported from: src/dyn_subgrid/dynHarvestMod.F90
#
# Handle reading of transient wood-harvest data and the C/N mortality fluxes
# (to product pools + litter) that result from harvest, for the non-FATES CN
# path.
#
# Public functions (Fortran names preserved):
#   dynHarvest_init                        — open the harvest file, set up the 5
#                                            HARVEST_* DynVarTimeUninterp readers.
#   dynHarvest_interp!                     — read current-year harvest rates,
#                                            summing the 5 types into `harvest`.
#   dynHarvest_interp_resolve_harvesttypes — read current-year rates keeping the
#                                            5 types distinct (FATES path).
#   cn_harvest!            (= CNHarvest)            — patch-level harvest mortality
#                                            C & N fluxes given the harvest fraction.
#   cn_harvest_pft_to_column! (= CNHarvestPftToColumn) — aggregate patch-level
#                                            harvest fluxes to the column level.
#
# Module-private singletons in the Fortran (dynHarvest_file, harvest_inst,
# harvest, do_harvest, harvest_units) are gathered into the explicit mutable
# struct `DynHarvestState`, which the caller passes in (NOT added to CLMInstances
# or any dual-copied struct).
# ==========================================================================

# --- harvest variable names + unit strings (dynHarvestMod parameters) ---
const num_harvest_inst = 5
const harvest_varnames = ["HARVEST_VH1", "HARVEST_VH2", "HARVEST_SH1",
                          "HARVEST_SH2", "HARVEST_SH3"]
const harvest_string_not_set = "not_set"
const harvest_mass_units     = "gC/m2/yr"
const harvest_unitless_units = "unitless"

"""
    DynHarvestState

Holds the module-private state of `dynHarvestMod` (the Fortran module-level
singletons): the open harvest file, the 5 per-type harvest-rate readers, the
gridcell-level summed harvest rate, the `do_harvest` switch and the units string.

The caller owns one instance of this and threads it through `dynHarvest_init`,
`dynHarvest_interp!` and `cn_harvest!`. It must NOT be placed in any
dual-copied / AD instance struct.

Ported from the module-level variables in `dynHarvestMod.F90`.
"""
Base.@kwdef mutable struct DynHarvestState
    dynHarvest_file::Union{DynFile,Nothing} = nothing            # file with harvest data
    harvest_inst::Vector{DynVarTimeUninterp} = DynVarTimeUninterp[]  # 5 per-type readers
    harvest::Vector{Float64} = Float64[]                         # summed harvest rate (per gridcell)
    do_harvest::Bool = false                                     # are we in a harvest period?
    harvest_units::String = harvest_string_not_set              # units from the file
end

# ---------------------------------------------------------------------------
# dynHarvest_init — open the file, set up the 5 harvest-rate variables.
# ---------------------------------------------------------------------------

"""
    dynHarvest_init(begg, endg, harvest_filename; current_year=nothing,
                    use_fates=false, dim1name=grlnd)

Initialize data structures for harvest information. Opens the harvest file as a
`DynFile` (year taken from the START of the timestep) and constructs the 5
`HARVEST_*` `DynVarTimeUninterp` readers. Validates and records the units.

`begg:endg` is the proc-level gridcell range (1-based here). Returns a fully
populated `DynHarvestState`.

Ported from `dynHarvest_init` in `dynHarvestMod.F90`.
"""
function dynHarvest_init(begg::Int, endg::Int, harvest_filename::String;
                         current_year::Union{Int,Nothing} = nothing,
                         use_fates::Bool = false,
                         dim1name::String = grlnd)

    state = DynHarvestState()

    # We only need the summary `harvest` array in the non-FATES CN path.
    if !use_fates
        state.harvest = zeros(Float64, endg - begg + 1)
    end

    # Get the year from the START of the timestep, for consistency with the
    # other dyn-file machinery.
    state.dynHarvest_file = dyn_file_open(harvest_filename,
                                          YEAR_POSITION_START_OF_TIMESTEP;
                                          current_year = current_year)

    num_points = endg - begg + 1
    state.harvest_inst = Vector{DynVarTimeUninterp}(undef, num_harvest_inst)

    for varnum in 1:num_harvest_inst
        state.harvest_inst[varnum] = dyn_var_time_uninterp(
            state.dynHarvest_file, harvest_varnames[varnum], dim1name,
            1.0, false, [num_points])

        # Read + validate the units attribute.
        units = _harvest_get_units(state.dynHarvest_file, harvest_varnames[varnum])
        if units == harvest_string_not_set
            units = harvest_unitless_units
        elseif units == harvest_unitless_units
            # ok
        elseif units != harvest_mass_units
            error("dynHarvest_init: bad units read in from file = $units")
        end
        if varnum > 1 && units != state.harvest_units
            error("dynHarvest_init: harvest units are inconsistent on file = " *
                  harvest_filename)
        end
        state.harvest_units = units
    end

    return state
end

# Read the "units" attribute of a variable on the harvest file. Mirrors the
# Fortran get_att("units", ...) call, which leaves `units` at its initial
# "not_set" value when the attribute is absent.
function _harvest_get_units(df::DynFile, varname::String)
    ds = df.ds
    (ds === nothing || !haskey(ds, varname)) && return harvest_string_not_set
    ncvar = ds[varname]
    return haskey(ncvar.attrib, "units") ? String(ncvar.attrib["units"]) :
           harvest_string_not_set
end

# ---------------------------------------------------------------------------
# dynHarvest_interp! — read current-year harvest rates (summed over types).
# ---------------------------------------------------------------------------

"""
    dynHarvest_interp!(state, begg, endg, year)

Get harvest data for the current model `year`, when needed. Harvest data are
stored as rates (not weights), so no time interpolation is performed — the rate
is held constant through the year.

Sets `state.harvest[begg:endg]` to the sum of the 5 harvest types and sets
`state.do_harvest`. Before the start of the time series, harvest is turned off;
past the end, the last year's rates are maintained.

Ported from `dynHarvest_interp` in `dynHarvestMod.F90`.
"""
function dynHarvest_interp!(state::DynHarvestState, begg::Int, endg::Int, year::Int)
    set_current_year!(state.dynHarvest_file, year)

    # Get total harvest for this time step.
    fill!(state.harvest, 0.0)

    if is_before_time_series(state.dynHarvest_file.time_info)
        # Turn off harvest before the start of the harvest time series.
        state.do_harvest = false
    else
        # do_harvest stays true even past the end of the time series: harvest
        # rates are maintained at the last year's rate for all later years.
        state.do_harvest = true
        for varnum in 1:num_harvest_inst
            this_data = get_current_data_1d(state.harvest_inst[varnum])
            @inbounds for i in 1:(endg - begg + 1)
                state.harvest[i] += this_data[i]
            end
        end
    end

    return nothing
end

# ---------------------------------------------------------------------------
# dynHarvest_interp_resolve_harvesttypes — keep the 5 types distinct (FATES).
# ---------------------------------------------------------------------------

"""
    dynHarvest_interp_resolve_harvesttypes(state, begg, endg, year)
        -> (harvest_rates, after_start_of_harvest_ts)

Like [`dynHarvest_interp!`](@ref) but keeps the different forcing sets distinct
(e.g. for passing to FATES, which has distinct primary and secondary lands), so
the result is a `(num_points, num_harvest_inst)` matrix of per-type rates.

Returns `(harvest_rates, after_start_of_harvest_ts)`.

Ported from `dynHarvest_interp_resolve_harvesttypes` in `dynHarvestMod.F90`.
"""
function dynHarvest_interp_resolve_harvesttypes(state::DynHarvestState,
                                                begg::Int, endg::Int, year::Int)
    set_current_year!(state.dynHarvest_file, year)
    num_points = endg - begg + 1
    harvest_rates = zeros(Float64, num_points, num_harvest_inst)

    if is_before_time_series(state.dynHarvest_file.time_info)
        after_start_of_harvest_ts = false
        # harvest_rates already zeroed.
    else
        after_start_of_harvest_ts = true
        for varnum in 1:num_harvest_inst
            this_data = get_current_data_1d(state.harvest_inst[varnum])
            @inbounds for i in 1:num_points
                harvest_rates[i, varnum] = this_data[i]
            end
        end
    end

    return harvest_rates, after_start_of_harvest_ts
end

# ---------------------------------------------------------------------------
# cn_harvest! (= CNHarvest) — patch-level harvest mortality C & N fluxes.
# ---------------------------------------------------------------------------

"""
    cn_harvest!(state, mask_soilp, patch, pftcon, soilbgc_state,
                cnveg_cs, cnveg_ns, cnveg_cf, cnveg_nf;
                dt, is_beg_curr_year, noveg=noveg, nbrdlf_evr_shrub=nbrdlf_evr_shrub,
                nlevdecomp=1, i_litr_min=1, i_litr_max=3, i_met_lit=1)

Harvest mortality routine for the coupled carbon-nitrogen code (CN), non-FATES,
non-matrix path. Computes the patch-level harvest mortality C and N fluxes from
each vegetation pool (leaf/froot/livestem/deadstem/livecroot/deadcroot + their
storage/transfer pools, gresp, retransn, xsmrpool), routing dead stem to the
product pool (`wood_harvest*`) and everything else to litter, then calls
[`cn_harvest_pft_to_column!`](@ref).

The annual harvest "mortality" rate `am` is derived from `state.harvest[g]`
(the per-gridcell summed harvest rate). For `gC/m2/yr` units it is converted to
a fraction (`min(0.98, cm/thistreec)`); for unitless units it is the rate
directly. All harvest is applied at the start of the year (`is_beg_curr_year`).

Ported from `CNHarvest` in `dynHarvestMod.F90` (non-matrix branch).
"""
function cn_harvest!(state::DynHarvestState,
                     mask_soilp::AbstractVector{Bool},
                     patch::PatchData,
                     pftcon,
                     soilbgc_state::SoilBiogeochemStateData,
                     cnveg_cs::CNVegCarbonStateData,
                     cnveg_ns::CNVegNitrogenStateData,
                     cnveg_cf::CNVegCarbonFluxData,
                     cnveg_nf::CNVegNitrogenFluxData;
                     dt::Real = 1800.0,
                     is_beg_curr_year::Bool = true,
                     noveg::Int = noveg,
                     nbrdlf_evr_shrub::Int = nbrdlf_evr_shrub,
                     nlevdecomp::Int = 1,
                     i_litr_min::Int = 1,
                     i_litr_max::Int = 3,
                     i_met_lit::Int = 1)

    ivt = patch.itype

    # --- C state inputs ---
    leafc              = cnveg_cs.leafc_patch
    frootc             = cnveg_cs.frootc_patch
    livestemc          = cnveg_cs.livestemc_patch
    deadstemc          = cnveg_cs.deadstemc_patch
    livecrootc         = cnveg_cs.livecrootc_patch
    deadcrootc         = cnveg_cs.deadcrootc_patch
    xsmrpool           = cnveg_cs.xsmrpool_patch
    leafc_storage      = cnveg_cs.leafc_storage_patch
    frootc_storage     = cnveg_cs.frootc_storage_patch
    livestemc_storage  = cnveg_cs.livestemc_storage_patch
    deadstemc_storage  = cnveg_cs.deadstemc_storage_patch
    livecrootc_storage = cnveg_cs.livecrootc_storage_patch
    deadcrootc_storage = cnveg_cs.deadcrootc_storage_patch
    gresp_storage      = cnveg_cs.gresp_storage_patch
    leafc_xfer         = cnveg_cs.leafc_xfer_patch
    frootc_xfer        = cnveg_cs.frootc_xfer_patch
    livestemc_xfer     = cnveg_cs.livestemc_xfer_patch
    deadstemc_xfer     = cnveg_cs.deadstemc_xfer_patch
    livecrootc_xfer    = cnveg_cs.livecrootc_xfer_patch
    deadcrootc_xfer    = cnveg_cs.deadcrootc_xfer_patch
    gresp_xfer         = cnveg_cs.gresp_xfer_patch

    # --- N state inputs ---
    leafn              = cnveg_ns.leafn_patch
    frootn             = cnveg_ns.frootn_patch
    livestemn          = cnveg_ns.livestemn_patch
    deadstemn          = cnveg_ns.deadstemn_patch
    livecrootn         = cnveg_ns.livecrootn_patch
    deadcrootn         = cnveg_ns.deadcrootn_patch
    retransn           = cnveg_ns.retransn_patch
    leafn_storage      = cnveg_ns.leafn_storage_patch
    frootn_storage     = cnveg_ns.frootn_storage_patch
    livestemn_storage  = cnveg_ns.livestemn_storage_patch
    deadstemn_storage  = cnveg_ns.deadstemn_storage_patch
    livecrootn_storage = cnveg_ns.livecrootn_storage_patch
    deadcrootn_storage = cnveg_ns.deadcrootn_storage_patch
    leafn_xfer         = cnveg_ns.leafn_xfer_patch
    frootn_xfer        = cnveg_ns.frootn_xfer_patch
    livestemn_xfer     = cnveg_ns.livestemn_xfer_patch
    deadstemn_xfer     = cnveg_ns.deadstemn_xfer_patch
    livecrootn_xfer    = cnveg_ns.livecrootn_xfer_patch
    deadcrootn_xfer    = cnveg_ns.deadcrootn_xfer_patch

    # --- C flux outputs ---
    hrv_leafc_to_litter              = cnveg_cf.hrv_leafc_to_litter_patch
    hrv_frootc_to_litter             = cnveg_cf.hrv_frootc_to_litter_patch
    hrv_livestemc_to_litter          = cnveg_cf.hrv_livestemc_to_litter_patch
    wood_harvestc                    = cnveg_cf.wood_harvestc_patch
    hrv_livecrootc_to_litter         = cnveg_cf.hrv_livecrootc_to_litter_patch
    hrv_deadcrootc_to_litter         = cnveg_cf.hrv_deadcrootc_to_litter_patch
    hrv_xsmrpool_to_atm              = cnveg_cf.hrv_xsmrpool_to_atm_patch
    hrv_leafc_storage_to_litter      = cnveg_cf.hrv_leafc_storage_to_litter_patch
    hrv_frootc_storage_to_litter     = cnveg_cf.hrv_frootc_storage_to_litter_patch
    hrv_livestemc_storage_to_litter  = cnveg_cf.hrv_livestemc_storage_to_litter_patch
    hrv_deadstemc_storage_to_litter  = cnveg_cf.hrv_deadstemc_storage_to_litter_patch
    hrv_livecrootc_storage_to_litter = cnveg_cf.hrv_livecrootc_storage_to_litter_patch
    hrv_deadcrootc_storage_to_litter = cnveg_cf.hrv_deadcrootc_storage_to_litter_patch
    hrv_gresp_storage_to_litter      = cnveg_cf.hrv_gresp_storage_to_litter_patch
    hrv_leafc_xfer_to_litter         = cnveg_cf.hrv_leafc_xfer_to_litter_patch
    hrv_frootc_xfer_to_litter        = cnveg_cf.hrv_frootc_xfer_to_litter_patch
    hrv_livestemc_xfer_to_litter     = cnveg_cf.hrv_livestemc_xfer_to_litter_patch
    hrv_deadstemc_xfer_to_litter     = cnveg_cf.hrv_deadstemc_xfer_to_litter_patch
    hrv_livecrootc_xfer_to_litter    = cnveg_cf.hrv_livecrootc_xfer_to_litter_patch
    hrv_deadcrootc_xfer_to_litter    = cnveg_cf.hrv_deadcrootc_xfer_to_litter_patch
    hrv_gresp_xfer_to_litter         = cnveg_cf.hrv_gresp_xfer_to_litter_patch

    # --- N flux outputs ---
    hrv_leafn_to_litter              = cnveg_nf.hrv_leafn_to_litter_patch
    hrv_frootn_to_litter             = cnveg_nf.hrv_frootn_to_litter_patch
    hrv_livestemn_to_litter          = cnveg_nf.hrv_livestemn_to_litter_patch
    wood_harvestn                    = cnveg_nf.wood_harvestn_patch
    hrv_livecrootn_to_litter         = cnveg_nf.hrv_livecrootn_to_litter_patch
    hrv_deadcrootn_to_litter         = cnveg_nf.hrv_deadcrootn_to_litter_patch
    hrv_retransn_to_litter           = cnveg_nf.hrv_retransn_to_litter_patch
    hrv_leafn_storage_to_litter      = cnveg_nf.hrv_leafn_storage_to_litter_patch
    hrv_frootn_storage_to_litter     = cnveg_nf.hrv_frootn_storage_to_litter_patch
    hrv_livestemn_storage_to_litter  = cnveg_nf.hrv_livestemn_storage_to_litter_patch
    hrv_deadstemn_storage_to_litter  = cnveg_nf.hrv_deadstemn_storage_to_litter_patch
    hrv_livecrootn_storage_to_litter = cnveg_nf.hrv_livecrootn_storage_to_litter_patch
    hrv_deadcrootn_storage_to_litter = cnveg_nf.hrv_deadcrootn_storage_to_litter_patch
    hrv_leafn_xfer_to_litter         = cnveg_nf.hrv_leafn_xfer_to_litter_patch
    hrv_frootn_xfer_to_litter        = cnveg_nf.hrv_frootn_xfer_to_litter_patch
    hrv_livestemn_xfer_to_litter     = cnveg_nf.hrv_livestemn_xfer_to_litter_patch
    hrv_deadstemn_xfer_to_litter     = cnveg_nf.hrv_deadstemn_xfer_to_litter_patch
    hrv_livecrootn_xfer_to_litter    = cnveg_nf.hrv_livecrootn_xfer_to_litter_patch
    hrv_deadcrootn_xfer_to_litter    = cnveg_nf.hrv_deadcrootn_xfer_to_litter_patch

    harvest     = state.harvest
    do_harvest  = state.do_harvest
    is_mass     = (state.harvest_units == harvest_mass_units)

    @inbounds for p in eachindex(mask_soilp)
        mask_soilp[p] || continue
        g = patch.gridcell[p]

        # Only tree PFTs (Fortran: ivt > noveg .and. ivt < nbrdlf_evr_shrub).
        if ivt[p] > noveg && ivt[p] < nbrdlf_evr_shrub

            if do_harvest
                if is_mass
                    thistreec = leafc[p] + frootc[p] + livestemc[p] +
                                deadstemc[p] + livecrootc[p] + deadcrootc[p] +
                                xsmrpool[p]
                    cm = harvest[g]
                    if thistreec > 0.0
                        am = min(0.98, cm / thistreec)  # harvest up to 98% so regrowth is possible
                    else
                        am = 0.0
                    end
                else
                    am = harvest[g]
                end

                # Apply all harvest at the start of the year.
                m = is_beg_curr_year ? am / dt : 0.0
            else
                m = 0.0
            end

            # ---- patch-level harvest carbon fluxes ----
            # displayed pools
            hrv_leafc_to_litter[p]              = leafc[p]              * m
            hrv_frootc_to_litter[p]             = frootc[p]             * m
            hrv_livestemc_to_litter[p]          = livestemc[p]          * m
            wood_harvestc[p]                    = deadstemc[p]          * m
            hrv_livecrootc_to_litter[p]         = livecrootc[p]         * m
            hrv_deadcrootc_to_litter[p]         = deadcrootc[p]         * m
            hrv_xsmrpool_to_atm[p]              = xsmrpool[p]           * m

            # storage pools
            hrv_leafc_storage_to_litter[p]      = leafc_storage[p]      * m
            hrv_frootc_storage_to_litter[p]     = frootc_storage[p]     * m
            hrv_livestemc_storage_to_litter[p]  = livestemc_storage[p]  * m
            hrv_deadstemc_storage_to_litter[p]  = deadstemc_storage[p]  * m
            hrv_livecrootc_storage_to_litter[p] = livecrootc_storage[p] * m
            hrv_deadcrootc_storage_to_litter[p] = deadcrootc_storage[p] * m
            hrv_gresp_storage_to_litter[p]      = gresp_storage[p]      * m

            # transfer pools
            hrv_leafc_xfer_to_litter[p]         = leafc_xfer[p]         * m
            hrv_frootc_xfer_to_litter[p]        = frootc_xfer[p]        * m
            hrv_livestemc_xfer_to_litter[p]     = livestemc_xfer[p]     * m
            hrv_deadstemc_xfer_to_litter[p]     = deadstemc_xfer[p]     * m
            hrv_livecrootc_xfer_to_litter[p]    = livecrootc_xfer[p]    * m
            hrv_deadcrootc_xfer_to_litter[p]    = deadcrootc_xfer[p]    * m
            hrv_gresp_xfer_to_litter[p]         = gresp_xfer[p]         * m

            # ---- patch-level harvest mortality nitrogen fluxes ----
            # displayed pools
            hrv_leafn_to_litter[p]              = leafn[p]              * m
            hrv_frootn_to_litter[p]             = frootn[p]             * m
            hrv_livestemn_to_litter[p]          = livestemn[p]          * m
            wood_harvestn[p]                    = deadstemn[p]          * m
            hrv_livecrootn_to_litter[p]         = livecrootn[p]         * m
            hrv_deadcrootn_to_litter[p]         = deadcrootn[p]         * m
            hrv_retransn_to_litter[p]           = retransn[p]           * m

            # storage pools
            hrv_leafn_storage_to_litter[p]      = leafn_storage[p]      * m
            hrv_frootn_storage_to_litter[p]     = frootn_storage[p]     * m
            hrv_livestemn_storage_to_litter[p]  = livestemn_storage[p]  * m
            hrv_deadstemn_storage_to_litter[p]  = deadstemn_storage[p]  * m
            hrv_livecrootn_storage_to_litter[p] = livecrootn_storage[p] * m
            hrv_deadcrootn_storage_to_litter[p] = deadcrootn_storage[p] * m

            # transfer pools
            hrv_leafn_xfer_to_litter[p]         = leafn_xfer[p]         * m
            hrv_frootn_xfer_to_litter[p]        = frootn_xfer[p]        * m
            hrv_livestemn_xfer_to_litter[p]     = livestemn_xfer[p]     * m
            hrv_deadstemn_xfer_to_litter[p]     = deadstemn_xfer[p]     * m
            hrv_livecrootn_xfer_to_litter[p]    = livecrootn_xfer[p]    * m
            hrv_deadcrootn_xfer_to_litter[p]    = deadcrootn_xfer[p]    * m
        end  # end tree block
    end  # end patch loop

    # Gather all patch-level harvest litterfall fluxes to the column.
    cn_harvest_pft_to_column!(mask_soilp, patch, pftcon, soilbgc_state,
                              cnveg_cf, cnveg_nf;
                              nlevdecomp = nlevdecomp,
                              i_litr_min = i_litr_min, i_litr_max = i_litr_max,
                              i_met_lit = i_met_lit)

    return nothing
end

# ---------------------------------------------------------------------------
# cn_harvest_pft_to_column! (= CNHarvestPftToColumn)
# ---------------------------------------------------------------------------

"""
    cn_harvest_pft_to_column!(mask_soilp, patch, pftcon, soilbgc_state,
                              cnveg_cf, cnveg_nf;
                              nlevdecomp=1, i_litr_min=1, i_litr_max=3, i_met_lit=1)

Gather all patch-level harvest litterfall fluxes to the column level and assign
them to the litter / coarse-woody-debris (CWD) pools (using the leaf/froot/croot/
stem vertical profiles and the per-PFT litter fractions `lf_f`/`fr_f`), and the
product-pool (`wood_harvest*_col`) fluxes.

`soilbgc_state` supplies the vertical profiles
(`leaf_prof_patch`/`froot_prof_patch`/`croot_prof_patch`/`stem_prof_patch`);
`pftcon` supplies the per-PFT litter-fraction matrices (`lf_f`/`fr_f`). The
column accumulators are NOT zeroed here (the caller manages that), matching the
Fortran `+=` semantics.

Ported from `CNHarvestPftToColumn` in `dynHarvestMod.F90`.
"""
function cn_harvest_pft_to_column!(mask_soilp::AbstractVector{Bool},
                                   patch::PatchData,
                                   pftcon,
                                   soilbgc_state::SoilBiogeochemStateData,
                                   cnveg_cf::CNVegCarbonFluxData,
                                   cnveg_nf::CNVegNitrogenFluxData;
                                   nlevdecomp::Int = 1,
                                   i_litr_min::Int = 1,
                                   i_litr_max::Int = 3,
                                   i_met_lit::Int = 1)

    ivt   = patch.itype
    wtcol = patch.wtcol

    lf_f = pftcon.lf_f
    fr_f = pftcon.fr_f

    leaf_prof  = soilbgc_state.leaf_prof_patch
    froot_prof = soilbgc_state.froot_prof_patch
    croot_prof = soilbgc_state.croot_prof_patch
    stem_prof  = soilbgc_state.stem_prof_patch

    # --- patch-level C fluxes (inputs) ---
    hrv_leafc_to_litter              = cnveg_cf.hrv_leafc_to_litter_patch
    hrv_frootc_to_litter             = cnveg_cf.hrv_frootc_to_litter_patch
    hrv_livestemc_to_litter          = cnveg_cf.hrv_livestemc_to_litter_patch
    pwood_harvestc                   = cnveg_cf.wood_harvestc_patch
    hrv_livecrootc_to_litter         = cnveg_cf.hrv_livecrootc_to_litter_patch
    hrv_deadcrootc_to_litter         = cnveg_cf.hrv_deadcrootc_to_litter_patch
    hrv_leafc_storage_to_litter      = cnveg_cf.hrv_leafc_storage_to_litter_patch
    hrv_frootc_storage_to_litter     = cnveg_cf.hrv_frootc_storage_to_litter_patch
    hrv_livestemc_storage_to_litter  = cnveg_cf.hrv_livestemc_storage_to_litter_patch
    hrv_deadstemc_storage_to_litter  = cnveg_cf.hrv_deadstemc_storage_to_litter_patch
    hrv_livecrootc_storage_to_litter = cnveg_cf.hrv_livecrootc_storage_to_litter_patch
    hrv_deadcrootc_storage_to_litter = cnveg_cf.hrv_deadcrootc_storage_to_litter_patch
    hrv_gresp_storage_to_litter      = cnveg_cf.hrv_gresp_storage_to_litter_patch
    hrv_leafc_xfer_to_litter         = cnveg_cf.hrv_leafc_xfer_to_litter_patch
    hrv_frootc_xfer_to_litter        = cnveg_cf.hrv_frootc_xfer_to_litter_patch
    hrv_livestemc_xfer_to_litter     = cnveg_cf.hrv_livestemc_xfer_to_litter_patch
    hrv_deadstemc_xfer_to_litter     = cnveg_cf.hrv_deadstemc_xfer_to_litter_patch
    hrv_livecrootc_xfer_to_litter    = cnveg_cf.hrv_livecrootc_xfer_to_litter_patch
    hrv_deadcrootc_xfer_to_litter    = cnveg_cf.hrv_deadcrootc_xfer_to_litter_patch
    hrv_gresp_xfer_to_litter         = cnveg_cf.hrv_gresp_xfer_to_litter_patch

    cwood_harvestc      = cnveg_cf.wood_harvestc_col
    harvest_c_to_litr_c = cnveg_cf.harvest_c_to_litr_c_col
    harvest_c_to_cwdc   = cnveg_cf.harvest_c_to_cwdc_col

    # --- patch-level N fluxes (inputs) ---
    hrv_leafn_to_litter              = cnveg_nf.hrv_leafn_to_litter_patch
    hrv_frootn_to_litter             = cnveg_nf.hrv_frootn_to_litter_patch
    hrv_livestemn_to_litter          = cnveg_nf.hrv_livestemn_to_litter_patch
    pwood_harvestn                   = cnveg_nf.wood_harvestn_patch
    hrv_livecrootn_to_litter         = cnveg_nf.hrv_livecrootn_to_litter_patch
    hrv_deadcrootn_to_litter         = cnveg_nf.hrv_deadcrootn_to_litter_patch
    hrv_retransn_to_litter           = cnveg_nf.hrv_retransn_to_litter_patch
    hrv_leafn_storage_to_litter      = cnveg_nf.hrv_leafn_storage_to_litter_patch
    hrv_frootn_storage_to_litter     = cnveg_nf.hrv_frootn_storage_to_litter_patch
    hrv_livestemn_storage_to_litter  = cnveg_nf.hrv_livestemn_storage_to_litter_patch
    hrv_deadstemn_storage_to_litter  = cnveg_nf.hrv_deadstemn_storage_to_litter_patch
    hrv_livecrootn_storage_to_litter = cnveg_nf.hrv_livecrootn_storage_to_litter_patch
    hrv_deadcrootn_storage_to_litter = cnveg_nf.hrv_deadcrootn_storage_to_litter_patch
    hrv_leafn_xfer_to_litter         = cnveg_nf.hrv_leafn_xfer_to_litter_patch
    hrv_frootn_xfer_to_litter        = cnveg_nf.hrv_frootn_xfer_to_litter_patch
    hrv_livestemn_xfer_to_litter     = cnveg_nf.hrv_livestemn_xfer_to_litter_patch
    hrv_deadstemn_xfer_to_litter     = cnveg_nf.hrv_deadstemn_xfer_to_litter_patch
    hrv_livecrootn_xfer_to_litter    = cnveg_nf.hrv_livecrootn_xfer_to_litter_patch
    hrv_deadcrootn_xfer_to_litter    = cnveg_nf.hrv_deadcrootn_xfer_to_litter_patch

    cwood_harvestn      = cnveg_nf.wood_harvestn_col
    harvest_n_to_litr_n = cnveg_nf.harvest_n_to_litr_n_col
    harvest_n_to_cwdn   = cnveg_nf.harvest_n_to_cwdn_col

    @inbounds for j in 1:nlevdecomp
        for p in eachindex(mask_soilp)
            mask_soilp[p] || continue
            c    = patch.column[p]
            ivtp = ivt[p] + 1   # 0-based Fortran -> 1-based Julia
            wt   = wtcol[p]
            lp   = leaf_prof[p, j]
            fpr  = froot_prof[p, j]
            cp   = croot_prof[p, j]
            sp   = stem_prof[p, j]

            for i in i_litr_min:i_litr_max
                # leaf harvest mortality carbon fluxes
                harvest_c_to_litr_c[c, j, i] += hrv_leafc_to_litter[p]  * lf_f[ivtp, i] * wt * lp
                # fine root harvest mortality carbon fluxes
                harvest_c_to_litr_c[c, j, i] += hrv_frootc_to_litter[p] * fr_f[ivtp, i] * wt * fpr
            end

            # wood harvest mortality carbon fluxes (to CWD)
            harvest_c_to_cwdc[c, j] += hrv_livestemc_to_litter[p]  * wt * sp
            harvest_c_to_cwdc[c, j] += hrv_livecrootc_to_litter[p] * wt * cp
            harvest_c_to_cwdc[c, j] += hrv_deadcrootc_to_litter[p] * wt * cp

            # storage + transfer harvest mortality carbon fluxes (to metabolic litter)
            harvest_c_to_litr_c[c, j, i_met_lit] +=
                hrv_leafc_storage_to_litter[p]      * wt * lp +
                hrv_frootc_storage_to_litter[p]     * wt * fpr +
                hrv_livestemc_storage_to_litter[p]  * wt * sp +
                hrv_deadstemc_storage_to_litter[p]  * wt * sp +
                hrv_livecrootc_storage_to_litter[p] * wt * cp +
                hrv_deadcrootc_storage_to_litter[p] * wt * cp +
                hrv_gresp_storage_to_litter[p]      * wt * lp +
                hrv_leafc_xfer_to_litter[p]         * wt * lp +
                hrv_frootc_xfer_to_litter[p]        * wt * fpr +
                hrv_livestemc_xfer_to_litter[p]     * wt * sp +
                hrv_deadstemc_xfer_to_litter[p]     * wt * sp +
                hrv_livecrootc_xfer_to_litter[p]    * wt * cp +
                hrv_deadcrootc_xfer_to_litter[p]    * wt * cp +
                hrv_gresp_xfer_to_litter[p]         * wt * lp

            for i in i_litr_min:i_litr_max
                # leaf + fine root harvest mortality nitrogen fluxes
                harvest_n_to_litr_n[c, j, i] +=
                    hrv_leafn_to_litter[p]  * lf_f[ivtp, i] * wt * lp +
                    hrv_frootn_to_litter[p] * fr_f[ivtp, i] * wt * fpr
            end

            # wood harvest mortality nitrogen fluxes (to CWD)
            harvest_n_to_cwdn[c, j] += hrv_livestemn_to_litter[p]  * wt * sp
            harvest_n_to_cwdn[c, j] += hrv_livecrootn_to_litter[p] * wt * cp
            harvest_n_to_cwdn[c, j] += hrv_deadcrootn_to_litter[p] * wt * cp

            # retransn + storage + transfer harvest mortality N (to metabolic litter)
            harvest_n_to_litr_n[c, j, i_met_lit] +=
                hrv_retransn_to_litter[p]           * wt * lp +
                hrv_leafn_storage_to_litter[p]      * wt * lp +
                hrv_frootn_storage_to_litter[p]     * wt * fpr +
                hrv_livestemn_storage_to_litter[p]  * wt * sp +
                hrv_deadstemn_storage_to_litter[p]  * wt * sp +
                hrv_livecrootn_storage_to_litter[p] * wt * cp +
                hrv_deadcrootn_storage_to_litter[p] * wt * cp +
                hrv_leafn_xfer_to_litter[p]         * wt * lp +
                hrv_frootn_xfer_to_litter[p]        * wt * fpr +
                hrv_livestemn_xfer_to_litter[p]     * wt * sp +
                hrv_deadstemn_xfer_to_litter[p]     * wt * sp +
                hrv_livecrootn_xfer_to_litter[p]    * wt * cp +
                hrv_deadcrootn_xfer_to_litter[p]    * wt * cp
        end
    end

    # wood harvest mortality C & N fluxes to the product pools
    @inbounds for p in eachindex(mask_soilp)
        mask_soilp[p] || continue
        c = patch.column[p]
        cwood_harvestc[c] += pwood_harvestc[p] * wtcol[p]
        cwood_harvestn[c] += pwood_harvestn[p] * wtcol[p]
    end

    return nothing
end
