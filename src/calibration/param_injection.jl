# ==========================================================================
# Full 29-parameter injection for CLM.jl calibration
#
# Matches SYMFLUENCE's CLMWorker parameter application:
# - Namelist params: baseflow_scalar, int_snow_max (via init_* calls)
# - params.nc: snow params, PFT params, hydrology scalars (via NetCDF edit)
# - surfdata: soil multipliers, fmax, organic_max (via NetCDF edit)
# - routing: route_k (post-processing, not written to files)
# ==========================================================================

using NCDatasets

const PARAM_BOUNDS = Dict(
    "baseflow_scalar"        => (0.001, 0.1),
    "fff"                    => (0.02, 1.0),
    "wimp"                   => (0.01, 0.1),
    "ksatdecay"              => (0.1, 10.0),
    "n_baseflow"             => (0.5, 5.0),
    "e_ice"                  => (1.0, 6.0),
    "perched_baseflow_scalar"=> (1e-7, 1e-3),
    "interception_fraction"  => (0.2, 1.0),
    "max_leaf_wetted_frac"   => (0.01, 0.2),
    "fmax"                   => (0.0, 1.0),
    "bsw_mult"               => (0.5, 2.0),
    "sucsat_mult"            => (0.5, 2.0),
    "watsat_mult"            => (0.8, 1.2),
    "hksat_mult"             => (0.1, 10.0),
    "organic_max"            => (0.0, 130.0),
    "fresh_snw_rds_max"      => (50.0, 200.0),
    "snw_aging_bst"          => (0.0, 200.0),
    "SNO_Z0MV"               => (0.0001, 0.01),
    "accum_factor"           => (-0.1, 0.1),
    "SNOW_DENSITY_MAX"       => (250.0, 550.0),
    "SNOW_DENSITY_MIN"       => (50.0, 200.0),
    "n_melt_coef"            => (50.0, 500.0),
    "int_snow_max"           => (500.0, 5000.0),
    "medlynslope"            => (2.0, 12.0),
    "slatop"                 => (0.005, 0.06),
    "flnr"                   => (0.05, 0.25),
    "froot_leaf"             => (0.5, 3.0),
    "stem_leaf"              => (0.5, 3.0),
    "route_k"                => (1.0, 40.0),
)

const PARAM_NAMES = sort(collect(keys(PARAM_BOUNDS)))

# params.nc variables (scalar or PFT-indexed)
const PARAMS_NC_MAP = Dict(
    "fff"                     => "fff",
    "wimp"                    => "wimp",
    "ksatdecay"               => "pc",
    "n_baseflow"              => "n_baseflow",
    "e_ice"                   => "e_ice",
    "perched_baseflow_scalar" => "perched_baseflow_scalar",
    "interception_fraction"   => "interception_fraction",
    "max_leaf_wetted_frac"    => "maximum_leaf_wetted_fraction",
    "fresh_snw_rds_max"       => "fresh_snw_rds_max",
    "snw_aging_bst"           => "snw_aging_bst",  # also read as xdrdt fallback
    "SNO_Z0MV"                => "SNO_Z0MV",
    "accum_factor"            => "accum_factor",
    "SNOW_DENSITY_MAX"        => "SNOW_DENSITY_MAX",
    "SNOW_DENSITY_MIN"        => "SNOW_DENSITY_MIN",
    "n_melt_coef"             => "n_melt_coef",
    "medlynslope"             => "medlynslope",
    "slatop"                  => "slatop",
    "flnr"                    => "flnr",
    "froot_leaf"              => "froot_leaf",
    "stem_leaf"               => "stem_leaf",
)

const PFT_PARAMS = Set(["medlynslope", "slatop", "flnr", "froot_leaf", "stem_leaf"])

const SURFDATA_MULT_MAP = Dict(
    "bsw_mult"    => "bsw",
    "sucsat_mult" => "sucsat",
    "watsat_mult" => "watsat",
    "hksat_mult"  => "hksat",
)

"""
    apply_all_params!(paramfile, surfdata_file, params_dict;
                      base_paramfile=nothing, base_surfdata=nothing)

Apply all 29 calibration parameters to CLM input files.
Modifies paramfile and surfdata_file in-place.
If base files are provided, copies them first (fresh start each iteration).
"""
function apply_all_params!(paramfile::String, surfdata_file::String,
                           params::Dict{String, Float64};
                           base_paramfile::String="",
                           base_surfdata::String="")
    # Copy base files if provided (fresh each iteration)
    if !isempty(base_paramfile) && isfile(base_paramfile) && base_paramfile != paramfile
        cp(base_paramfile, paramfile, force=true)
    end
    if !isempty(base_surfdata) && isfile(base_surfdata) && base_surfdata != surfdata_file
        cp(base_surfdata, surfdata_file, force=true)
    end

    # 1. Update params.nc
    _update_params_nc!(paramfile, params)

    # 2. Update surfdata
    _update_surfdata!(surfdata_file, params)

    return nothing
end

function _update_params_nc!(filepath::String, params::Dict{String, Float64})
    ds = NCDataset(filepath, "a")
    try
        # Find active PFTs from surfdata (if available)
        active_pfts = Int[]

        for (pname, pval) in params
            haskey(PARAMS_NC_MAP, pname) || continue
            nc_var = PARAMS_NC_MAP[pname]
            if !haskey(ds, nc_var)
                if pname in PFT_PARAMS
                    continue  # PFT vars need pft dimension, skip if missing
                end
                defVar(ds, nc_var, Float64, ())  # create scalar variable
            end

            if pname in PFT_PARAMS
                # PFT-indexed: apply to all PFTs (CLM.jl uses all)
                v = ds[nc_var]
                dims = dimnames(v)
                if "pft" in dims
                    data = Array(v[:])
                    # Apply to PFTs 1-16 (skip bare ground index 0 for some)
                    for i in eachindex(data)
                        if data[i] != 0.0  # only modify non-zero PFTs
                            data[i] = pval
                        end
                    end
                    v[:] = data
                else
                    v[:] = fill(pval, size(v))
                end
            else
                ds[nc_var][:] = pval
            end
        end

        # Constraint: SNOW_DENSITY_MIN < SNOW_DENSITY_MAX
        if haskey(ds, "SNOW_DENSITY_MIN") && haskey(ds, "SNOW_DENSITY_MAX")
            dmin = Float64(ds["SNOW_DENSITY_MIN"][1])
            dmax = Float64(ds["SNOW_DENSITY_MAX"][1])
            if dmin >= dmax
                ds["SNOW_DENSITY_MIN"][:] = dmax * 0.5
            end
        end
    finally
        close(ds)
    end
end

function _update_surfdata!(filepath::String, params::Dict{String, Float64})
    ds = NCDataset(filepath, "a")
    try
        # FMAX
        if haskey(params, "fmax") && haskey(ds, "FMAX")
            ds["FMAX"][:] = params["fmax"]
        end

        # Organic max
        if haskey(params, "organic_max") && haskey(ds, "ORGANIC")
            org = Array(ds["ORGANIC"][:])
            org .= min.(org, params["organic_max"])
            ds["ORGANIC"][:] = org
        end

        # Soil multipliers (apply to base values)
        for (mult_name, nc_var) in SURFDATA_MULT_MAP
            haskey(params, mult_name) || continue
            haskey(ds, nc_var) || continue
            base_vals = Array(ds[nc_var][:])
            new_vals = base_vals .* params[mult_name]
            if mult_name == "watsat_mult"
                new_vals = clamp.(new_vals, 0.01, 0.95)
            elseif mult_name == "hksat_mult"
                new_vals = max.(new_vals, 1e-10)
            end
            ds[nc_var][:] = new_vals
        end
    finally
        close(ds)
    end
end

"""
    linear_reservoir_route(q_in, k)

Apply linear reservoir routing: Q_out(t) = (1-1/K)*Q_out(t-1) + (1/K)*Q_in(t)
"""
function linear_reservoir_route(q_in::Vector{Float64}, k::Float64)
    k <= 1.0 && return copy(q_in)
    α = 1.0 / k
    q_out = similar(q_in)
    q_out[1] = q_in[1]
    for i in 2:length(q_out)
        q_out[i] = (1.0 - α) * q_out[i-1] + α * q_in[i]
    end
    return q_out
end
