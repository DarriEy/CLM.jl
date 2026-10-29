# fates_driver_landuse.jl — host (CTSM) side of FATES land use: LUH2 reader,
# land-use x PFT static data, and the harvest-mode switches.
#
# Ports:
#   * src/dyn_subgrid/dynFATESLandUseChangeMod.F90 (dynFatesLandUseInit /
#     dynFatesLandUseInterp): the LUH2 state, transition and harvest time series,
#     read uninterpolated (one record per year, YEAR_POSITION_END_OF_TIMESTEP) on the
#     ported dyn_file_io.jl machinery.
#   * clmfates_interfaceMod.F90 GetLandusePFTData: the static land-use x PFT
#     fractions (frac_primr/frac_secnd/frac_pastr/frac_range) and frac_brgnd.
#   * clmfates_interfaceMod.F90 CLMFatesGlobals1 (lines 518-561): fates_harvest_mode
#     -> use_logging / use_lu_harvest / num_lu_harvest_cats.
#   * clmfates_interfaceMod.F90 dynamics_driv / init_coldstart: the per-site pack of
#     hlm_harvest_* and hlm_luh_* into bc_in.
#
# In Fortran the LUH arrays are module globals of dynFATESLandUseChangeMod; here
# they live in a `DynFatesLandUse` held by the FATES interface (`inst.fates.luh`),
# so nothing leaks between instances.

# Number of landuse transition, state and harvest variables (dynFATESLandUseChangeMod).
const num_landuse_transition_vars = 108
const num_landuse_state_vars = 12
const num_landuse_harvest_vars = 5
# clmfates_interfaceMod num_landuse_pft_vars (frac_primr, frac_secnd, frac_pastr, frac_range)
const num_landuse_pft_vars = 4

# fates_harvest_mode namelist values (dynFATESLandUseChangeMod)
const fates_harvest_no_logging   = "no_harvest"
const fates_harvest_logging_only = "event_code"
const fates_harvest_clmlanduse   = "landuse_timeseries"
const fates_harvest_luh_area     = "luhdata_area"
const fates_harvest_luh_mass     = "luhdata_mass"
const FATES_HARVEST_MODES = (fates_harvest_no_logging, fates_harvest_logging_only,
                             fates_harvest_clmlanduse, fates_harvest_luh_area,
                             fates_harvest_luh_mass)

# landuse harvest unit integer representation
const landuse_harvest_area_units = 1
const landuse_harvest_mass_units = 2

# LUH2 raw wood harvest area fraction / biomass carbon variable names
const landuse_harvest_area_varnames =
    ["primf_harv", "primn_harv", "secmf_harv", "secyf_harv", "secnf_harv"]
const landuse_harvest_mass_varnames =
    ["primf_bioh", "primn_bioh", "secmf_bioh", "secyf_bioh", "secnf_bioh"]

const landuse_state_varnames =
    ["primf", "primn", "secdf", "secdn", "pastr", "range",
     "urban", "c3ann", "c4ann", "c3per", "c4per", "c3nfx"]

const landuse_transition_varnames = [
    "primf_to_secdn", "primf_to_pastr", "primf_to_range", "primf_to_urban",
    "primf_to_c3ann", "primf_to_c4ann", "primf_to_c3per", "primf_to_c4per", "primf_to_c3nfx",
    "primn_to_secdf", "primn_to_pastr", "primn_to_range", "primn_to_urban",
    "primn_to_c3ann", "primn_to_c4ann", "primn_to_c3per", "primn_to_c4per", "primn_to_c3nfx",
    "secdf_to_secdn", "secdf_to_pastr", "secdf_to_range", "secdf_to_urban",
    "secdf_to_c3ann", "secdf_to_c4ann", "secdf_to_c3per", "secdf_to_c4per", "secdf_to_c3nfx",
    "secdn_to_secdf", "secdn_to_pastr", "secdn_to_range", "secdn_to_urban",
    "secdn_to_c3ann", "secdn_to_c4ann", "secdn_to_c3per", "secdn_to_c4per", "secdn_to_c3nfx",
    "pastr_to_secdf", "pastr_to_secdn", "pastr_to_range", "pastr_to_urban",
    "pastr_to_c3ann", "pastr_to_c4ann", "pastr_to_c3per", "pastr_to_c4per", "pastr_to_c3nfx",
    "range_to_secdf", "range_to_secdn", "range_to_pastr", "range_to_urban",
    "range_to_c3ann", "range_to_c4ann", "range_to_c3per", "range_to_c4per", "range_to_c3nfx",
    "urban_to_secdf", "urban_to_secdn", "urban_to_pastr", "urban_to_range",
    "urban_to_c3ann", "urban_to_c4ann", "urban_to_c3per", "urban_to_c4per", "urban_to_c3nfx",
    "c3ann_to_c4ann", "c3ann_to_c3per", "c3ann_to_c4per", "c3ann_to_c3nfx",
    "c3ann_to_secdf", "c3ann_to_secdn", "c3ann_to_pastr", "c3ann_to_range", "c3ann_to_urban",
    "c4ann_to_c3ann", "c4ann_to_c3per", "c4ann_to_c4per", "c4ann_to_c3nfx",
    "c4ann_to_secdf", "c4ann_to_secdn", "c4ann_to_pastr", "c4ann_to_range", "c4ann_to_urban",
    "c3per_to_c3ann", "c3per_to_c4ann", "c3per_to_c4per", "c3per_to_c3nfx",
    "c3per_to_secdf", "c3per_to_secdn", "c3per_to_pastr", "c3per_to_range", "c3per_to_urban",
    "c4per_to_c3ann", "c4per_to_c4ann", "c4per_to_c3per", "c4per_to_c3nfx",
    "c4per_to_secdf", "c4per_to_secdn", "c4per_to_pastr", "c4per_to_range", "c4per_to_urban",
    "c3nfx_to_c3ann", "c3nfx_to_c4ann", "c3nfx_to_c3per", "c3nfx_to_c4per",
    "c3nfx_to_secdf", "c3nfx_to_secdn", "c3nfx_to_pastr", "c3nfx_to_range", "c3nfx_to_urban"]

# GetLandusePFTData variable names and natpft size
const landuse_pft_map_varnames = ["frac_primr", "frac_secnd", "frac_pastr", "frac_range"]
const dim_landuse_pft = 14

@assert length(landuse_transition_varnames) == num_landuse_transition_vars
@assert length(landuse_state_varnames) == num_landuse_state_vars

"""
    fates_harvest_ctrl(harvest_mode) -> (pass_logging, pass_lu_harvest, pass_num_lu_harvest_cats)

The FATES logging / harvest control values for a CTSM `fates_harvest_mode`
(CLMFatesGlobals1, clmfates_interfaceMod.F90:518-539): any mode other than
`no_harvest` turns on logging; the LUH2 modes also turn on land-use harvest with
`num_landuse_harvest_vars` categories. `landuse_timeseries` (harvest from the CLM
flanduse_timeseries via dynHarvestMod) is not connected to FATES in the port.
"""
function fates_harvest_ctrl(harvest_mode::AbstractString)
    harvest_mode in FATES_HARVEST_MODES ||
        error("fates_harvest_mode = '$harvest_mode' is not one of $(FATES_HARVEST_MODES)")
    harvest_mode == fates_harvest_no_logging && return (ifalse, ifalse, 0)
    if harvest_mode == fates_harvest_clmlanduse
        error("fates_harvest_mode = 'landuse_timeseries' (CLM flanduse_timeseries " *
              "harvest via dynHarvestMod) is not connected to FATES in CLM.jl; use " *
              "'event_code', 'luhdata_area' or 'luhdata_mass'")
    elseif harvest_mode == fates_harvest_luh_area || harvest_mode == fates_harvest_luh_mass
        return (itrue, itrue, num_landuse_harvest_vars)
    end
    return (itrue, ifalse, 0)   # event_code: time-driven logging only
end

"""
Host LUH2 land-use data (the dynFATESLandUseChangeMod module state).

* `landuse_states[nstate, g]`, `landuse_transitions[ntrans, g]`, `landuse_harvest[nharv, g]`
  (fractions of gridcell; harvest in area fraction or kgC, per `landuse_harvest_units`);
* `file` / `*_vars` — the open LUH2 time series (absent under potential vegetation);
* `landuse_pft_map[g, natpft, 4]` / `landuse_bareground[g]` — GetLandusePFTData
  (only read under `use_fates_fixed_biogeog .and. use_fates_luh`).
"""
Base.@kwdef mutable struct DynFatesLandUse
    harvest_mode::String = fates_harvest_no_logging
    use_potentialveg::Bool = false
    landuse_states::Matrix{Float64} = zeros(num_landuse_state_vars, 0)
    landuse_transitions::Matrix{Float64} = zeros(num_landuse_transition_vars, 0)
    landuse_harvest::Matrix{Float64} = zeros(num_landuse_harvest_vars, 0)
    landuse_harvest_units::Int = 0
    landuse_harvest_varnames::Vector{String} = String[]
    file::Union{DynFile,Nothing} = nothing
    transition_vars::Vector{DynVarTimeUninterp} = DynVarTimeUninterp[]
    state_vars::Vector{DynVarTimeUninterp} = DynVarTimeUninterp[]
    harvest_vars::Vector{DynVarTimeUninterp} = DynVarTimeUninterp[]
    landuse_pft_map::Array{Float64,3} = zeros(0, dim_landuse_pft, num_landuse_pft_vars)
    landuse_bareground::Vector{Float64} = Float64[]
end

_luh_harvest_on(luh::DynFatesLandUse) =
    luh.harvest_mode == fates_harvest_luh_area || luh.harvest_mode == fates_harvest_luh_mass

"""
    dyn_fates_landuse_init(ng, landuse_filename; harvest_mode, use_potentialveg,
                           current_year) -> DynFatesLandUse

`dynFatesLandUseInit`: allocate the state / transition / harvest arrays (zero), and
unless under potential vegetation open the LUH2 time series (YEAR_POSITION_END_OF_
TIMESTEP), register the 12 state and 108 transition variables (plus the 5 harvest
variables of the `luhdata_area` / `luhdata_mass` harvest mode) and read the current
year with `init_state = true`. `current_year` is the host's `get_curr_date` year.
(The Fortran `use_cn` early return does not apply: FATES and CN are exclusive.)
"""
function dyn_fates_landuse_init(ng::Int, landuse_filename::AbstractString;
                                harvest_mode::AbstractString = fates_harvest_no_logging,
                                use_potentialveg::Bool = false,
                                current_year::Int)
    luh = DynFatesLandUse(harvest_mode = String(harvest_mode),
                          use_potentialveg = use_potentialveg,
                          landuse_states = zeros(num_landuse_state_vars, ng),
                          landuse_transitions = zeros(num_landuse_transition_vars, ng),
                          landuse_harvest = zeros(num_landuse_harvest_vars, ng))
    # Avoid initializing the landuse timeseries file in potential vegetation mode
    use_potentialveg && return luh

    isempty(landuse_filename) &&
        error("dyn_fates_landuse_init: fluh_timeseries is required when use_fates_luh " *
              "is set and use_fates_potentialveg is false")
    df = dyn_file_open(String(landuse_filename), YEAR_POSITION_END_OF_TIMESTEP;
                       current_year = current_year)
    luh.file = df
    shape = [ng]
    mk(v) = dyn_var_time_uninterp(df, v, "lndgrid", 1.0, false, shape)
    luh.transition_vars = [mk(v) for v in landuse_transition_varnames]
    luh.state_vars = [mk(v) for v in landuse_state_varnames]
    if _luh_harvest_on(luh)
        if harvest_mode == fates_harvest_luh_area
            luh.landuse_harvest_varnames = copy(landuse_harvest_area_varnames)
            luh.landuse_harvest_units = landuse_harvest_area_units
        else
            luh.landuse_harvest_varnames = copy(landuse_harvest_mass_varnames)
            luh.landuse_harvest_units = landuse_harvest_mass_units
        end
        luh.harvest_vars = [mk(v) for v in luh.landuse_harvest_varnames]
    end
    # FATES needs state data during initialization
    dyn_fates_landuse_interp!(luh, current_year; init_state = true)
    return luh
end

"""
    dyn_fates_landuse_interp!(luh, year; init_state=false)

`dynFatesLandUseInterp`: set the file's current year (the host `get_curr_date`
year) and copy that year's states, transitions and (LUH harvest modes) harvest into
the gridcell arrays. Before the start of the time series (and not at init) all three
are reset to zero. No interpolation: the values hold for the whole year.
"""
function dyn_fates_landuse_interp!(luh::DynFatesLandUse, year::Integer;
                                   init_state::Bool = false)
    df = luh.file
    df === nothing && return luh   # potential vegetation: arrays stay zero
    set_current_year!(df, Int(year))
    if is_before_time_series(df.time_info) && !init_state
        fill!(luh.landuse_transitions, 0.0)
        fill!(luh.landuse_states, 0.0)
        fill!(luh.landuse_harvest, 0.0)
    else
        for (k, v) in enumerate(luh.transition_vars)
            luh.landuse_transitions[k, :] .= get_current_data_1d(v)
        end
        for (k, v) in enumerate(luh.state_vars)
            luh.landuse_states[k, :] .= get_current_data_1d(v)
        end
        if _luh_harvest_on(luh)
            for (k, v) in enumerate(luh.harvest_vars)
                luh.landuse_harvest[k, :] .= get_current_data_1d(v)
            end
        end
    end
    return luh
end

"""
    get_landuse_pft_data!(luh, landuse_pft_file, ng)

`GetLandusePFTData`: read the static land-use x PFT fractions (`frac_primr`,
`frac_secnd`, `frac_pastr`, `frac_range`, each natpft = 14) and `frac_brgnd` into
`luh.landuse_pft_map[g, pft, k]` / `luh.landuse_bareground[g]`.
"""
function get_landuse_pft_data!(luh::DynFatesLandUse, landuse_pft_file::AbstractString,
                               ng::Int)
    isempty(landuse_pft_file) &&
        error("flandusepftdat is required when use_fates_luh and use_fates_fixed_biogeog are set")
    isfile(landuse_pft_file) || error("landuse x pft file not found: $landuse_pft_file")
    luh.landuse_pft_map = fill(fates_unset_r8, ng, dim_landuse_pft, num_landuse_pft_vars)
    luh.landuse_bareground = fill(fates_unset_r8, ng)
    NCDataset(String(landuse_pft_file), "r") do ds
        haskey(ds.dim, "natpft") && ds.dim["natpft"] == dim_landuse_pft ||
            error("natpft dimensions on the landuse x pft file do not match target array size")
        for (k, v) in enumerate(landuse_pft_map_varnames)
            haskey(ds, v) || error("ERROR: $v NOT on landuse x pft file")
            a = Float64.(replace(Array(ds[v]), missing => NaN))
            # file (natpft, lat, lon) -> Julia (lon, lat, natpft): gridcell axis first
            luh.landuse_pft_map[:, :, k] .= reshape(a, ng, dim_landuse_pft)
        end
        haskey(ds, "frac_brgnd") || error("ERROR: frac_brgnd NOT on landuse x pft file")
        luh.landuse_bareground .= vec(Float64.(replace(Array(ds["frac_brgnd"]), missing => NaN)))
    end
    return luh
end

"""
    fates_pack_bcin_landuse!(bc, luh, g)

Pack gridcell `g`'s land-use data into a site's `bc_in` (clmfates_interfaceMod
`dynamics_driv` and `init_coldstart`): under the LUH2 harvest modes the five harvest
rates, their names and units; under `use_fates_luh` the LUH2 states and transitions
with their names.
"""
function fates_pack_bcin_landuse!(bc::bc_in_type, luh::DynFatesLandUse, g::Int)
    if _luh_harvest_on(luh)
        bc.hlm_harvest_rates .= @view luh.landuse_harvest[:, g]
        bc.hlm_harvest_catnames .= luh.landuse_harvest_varnames
        bc.hlm_harvest_units = luh.landuse_harvest_units
    end
    if hlm_use_luh[] == itrue
        bc.hlm_luh_states .= @view luh.landuse_states[:, g]
        bc.hlm_luh_state_names .= landuse_state_varnames
        bc.hlm_luh_transitions .= @view luh.landuse_transitions[:, g]
        bc.hlm_luh_transition_names .= landuse_transition_varnames
    end
    return bc
end
