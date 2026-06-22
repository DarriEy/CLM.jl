# FatesParametersInterface.jl
# Julia port of FATES src/fates/main/FatesParametersInterface.F90
#
# The parameter-reader interface: parameter registration, metadata, and
# retrieval/storage. This is part of the FATES<->HLM interface; the actual file
# read (`Read_interface`) is host-supplied and stubbed here (no-op abstract).
# Type-bound procedures -> functions dispatching on fates_parameters_type, so
# names like Init!/FindIndex coexist with the IO modules' versions.
# Deps: FatesConstantsMod (fates_r8), FatesGlobals (fates_log, fates_endrun).

const max_params = 250
const max_dimensions = 2
const max_used_dimensions = 25
const param_string_length = 40

# Values returned from netcdf when inquiring about number of dimensions.
const dimension_shape_scalar = 0
const dimension_shape_1d = 1
const dimension_shape_2d = 2

# ---------------------------------------------------------------------------
# Dimension names in the fates namespace
# ---------------------------------------------------------------------------
const dimension_name_scalar = ""
const dimension_name_pft = "fates_pft"
const dimension_name_segment = "fates_segment"
const dimension_name_cwd = "fates_NCWD"
const dimension_name_lsc = "fates_litterclass"
const dimension_name_fsc = "fates_litterclass"
const dimension_name_allpfts = "fates_allpfts"
const dimension_name_variants = "fates_variants"
const dimension_name_hydr_organs = "fates_hydr_organs"
const dimension_name_prt_organs = "fates_plant_organs"
const dimension_name_leaf_age = "fates_leafage_class"
const dimension_name_history_size_bins = "fates_history_size_bins"
const dimension_name_history_age_bins = "fates_history_age_bins"
const dimension_name_history_height_bins = "fates_history_height_bins"
const dimension_name_history_coage_bins = "fates_history_coage_bins"
const dimension_name_hlm_pftno = "fates_hlm_pftno"
const dimension_name_history_damage_bins = "fates_history_damage_bins"
const dimension_name_damage = "fates_damage_class"
const dimension_name_landuse = "fates_landuseclass"

# Dimensions in the host namespace
const dimension_name_host_allpfts = "allpfts"

# ---------------------------------------------------------------------------
# parameter_type — one registered parameter.
# ---------------------------------------------------------------------------
Base.@kwdef mutable struct parameter_type
    name::String = ""
    sync_with_host::Bool = false
    dimension_shape::Int = 0
    dimension_sizes::Vector{Int} = zeros(Int, max_dimensions)
    dimension_names::Vector{String} = fill("", max_dimensions)
    dimension_lower_bound::Vector{Int} = ones(Int, max_dimensions)
    data::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
end

# ---------------------------------------------------------------------------
# fates_parameters_type — the registry of parameters.
# ---------------------------------------------------------------------------
Base.@kwdef mutable struct fates_parameters_type
    num_parameters::Int = 0
    parameters::Vector{parameter_type} = [parameter_type() for _ in 1:max_params]
end

# ---------------------------------------------------------------------------
# fates_param_reader_type — abstract reader (host land models implement Read!).
# In Julia this is an abstract supertype; the actual file read is stubbed.
# ---------------------------------------------------------------------------
abstract type fates_param_reader_type end

"""
    Read!(this::fates_param_reader_type, fates_params::fates_parameters_type)

Read 'fates_params' parameters from (HLM-provided) storage. Abstract/stubbed —
host land models override this. The foundation layer keeps registration logic
only; the file-read body is a no-op.
"""
function Read!(this::fates_param_reader_type, fates_params::fates_parameters_type)
    return nothing
end

# =====================================================================================

"""
    Init!(this::fates_parameters_type)

Reset the parameter count to zero.
"""
function Init!(this::fates_parameters_type)
    this.num_parameters = 0
    return nothing
end

# =====================================================================================

"""
    Destroy!(this::fates_parameters_type)

Release per-parameter data arrays.
"""
function Destroy!(this::fates_parameters_type)
    for n in 1:this.num_parameters
        this.parameters[n].data = Matrix{Float64}(undef, 0, 0)
    end
    return nothing
end

# =====================================================================================

"""
    RegisterParameter!(this, name, dimension_shape, dimension_names; sync_with_host=false, lower_bounds=nothing)

Register a new parameter with the given dimension shape and names.
"""
function RegisterParameter!(this::fates_parameters_type, name::AbstractString,
                            dimension_shape::Integer,
                            dimension_names::AbstractVector{<:AbstractString};
                            sync_with_host::Bool=false,
                            lower_bounds::Union{AbstractVector{<:Integer},Nothing}=nothing)

    this.num_parameters += 1
    i = this.num_parameters
    p = this.parameters[i]
    p.name = name
    p.dimension_shape = dimension_shape
    p.dimension_sizes = zeros(Int, max_dimensions)

    num_names = min(max_dimensions, length(dimension_names))
    p.dimension_names = fill("", max_dimensions)
    for n in 1:num_names
        p.dimension_names[n] = dimension_names[n]
    end

    p.sync_with_host = sync_with_host

    # Default 1-based bounds unless otherwise specified.
    p.dimension_lower_bound = ones(Int, max_dimensions)
    if lower_bounds !== nothing
        num_bounds = min(max_dimensions, length(lower_bounds))
        for n in 1:num_bounds
            p.dimension_lower_bound[n] = lower_bounds[n]
        end
    end
    return nothing
end

# =====================================================================================
# Retrieval (Fortran generic RetrieveParameter -> RetrieveParameter{Scalar,1D,2D})
# =====================================================================================

"""
    RetrieveParameterScalar(this, name) -> Float64
"""
function RetrieveParameterScalar(this::fates_parameters_type, name::AbstractString)
    i = FindIndex(this, name)
    return this.parameters[i].data[1, 1]
end

"""
    RetrieveParameter1D(this, name, data)

Copy a 1-D parameter into `data` (length must match the stored size).
"""
function RetrieveParameter1D(this::fates_parameters_type, name::AbstractString,
                             data::AbstractVector)
    i = FindIndex(this, name)
    if length(data) != size(this.parameters[i].data, 1)
        @warn "RetrieveParameter1d : $name size inconsistent. expected=$(length(data)) " *
              "received=$(size(this.parameters[i].data, 1))"
        fates_endrun("size error retrieving 1d parameter.")
    end
    data .= this.parameters[i].data[:, 1]
    return nothing
end

"""
    RetrieveParameter2D(this, name, data)

Copy a 2-D parameter into `data`.
"""
function RetrieveParameter2D(this::fates_parameters_type, name::AbstractString,
                             data::AbstractMatrix)
    i = FindIndex(this, name)
    if size(data, 1) != size(this.parameters[i].data, 1) &&
       size(data, 2) != size(this.parameters[i].data, 2)
        @warn "RetrieveParameter2d : $name size inconsistent."
        fates_endrun("size error retrieving 2d parameter.")
    end
    data .= this.parameters[i].data
    return nothing
end

"""
    RetrieveParameter1DAllocate(this, name) -> Vector

Return a freshly allocated 1-D parameter honoring its lower-bound offset (the
returned vector is 1-based but sized to the stored extent).
"""
function RetrieveParameter1DAllocate(this::fates_parameters_type, name::AbstractString)
    i = FindIndex(this, name)
    return copy(this.parameters[i].data[:, 1])
end

"""
    RetrieveParameter2DAllocate(this, name) -> Matrix
"""
function RetrieveParameter2DAllocate(this::fates_parameters_type, name::AbstractString)
    i = FindIndex(this, name)
    return copy(this.parameters[i].data)
end

# =====================================================================================

"""
    FindIndex(this::fates_parameters_type, name) -> i

Return the index of the parameter named `name`. Returns `num_parameters+1` if
not found (mirroring the Fortran loop-fall-through behaviour).
"""
function FindIndex(this::fates_parameters_type, name::AbstractString)
    i = 1
    while i <= this.num_parameters
        if strip(this.parameters[i].name) == strip(name)
            break
        end
        i += 1
    end
    # if i > num_parameters: parameter name not found (Fortran leaves this silent)
    return i
end

# =====================================================================================

num_params(this::fates_parameters_type) = this.num_parameters

# =====================================================================================

"""
    GetUsedDimensions(this, is_host_file) -> (num_used_dimensions, used_dimensions)

Build the list of unique dimension names used by the parameters matching the
`is_host_file` selector.
"""
function GetUsedDimensions(this::fates_parameters_type, is_host_file::Bool)
    used_dimensions = fill("", max_used_dimensions)
    num_used_dimensions = 0

    for p in 1:this.num_parameters
        if is_host_file == this.parameters[p].sync_with_host
            for d in 1:max_dimensions
                dim_name = this.parameters[p].dimension_names[d]
                if length(strip(dim_name)) != 0
                    # Check if it's already in the list.
                    found = false
                    for ii in 1:num_used_dimensions
                        if used_dimensions[ii] == dim_name
                            found = true
                            break
                        end
                    end
                    if !found
                        num_used_dimensions += 1
                        used_dimensions[num_used_dimensions] = dim_name
                    end
                end
            end
        end
    end
    return num_used_dimensions, used_dimensions
end

# =====================================================================================

"""
    SetDimensionSizes!(this, is_host_file, num_used_dimensions, dimension_names, dimension_sizes)

Set each parameter's per-dimension size from the host-supplied name->size list.
"""
function SetDimensionSizes!(this::fates_parameters_type, is_host_file::Bool,
                            num_used_dimensions::Integer,
                            dimension_names::AbstractVector{<:AbstractString},
                            dimension_sizes::AbstractVector{<:Integer})
    for p in 1:this.num_parameters
        if is_host_file == this.parameters[p].sync_with_host
            for d in 1:max_dimensions
                dim_name = this.parameters[p].dimension_names[d]
                if length(strip(dim_name)) != 0
                    for ii in 1:num_used_dimensions
                        if strip(dimension_names[ii]) == strip(dim_name)
                            this.parameters[p].dimension_sizes[d] = dimension_sizes[ii]
                            break
                        end
                    end
                end
            end
        end
    end
    return nothing
end

# =====================================================================================

"""
    GetMetaData(this, index) -> (name, dimension_shape, dimension_sizes, dimension_names, is_host_param)
"""
function GetMetaData(this::fates_parameters_type, index::Integer)
    p = this.parameters[index]
    return p.name, p.dimension_shape, copy(p.dimension_sizes),
           copy(p.dimension_names), p.sync_with_host
end

# =====================================================================================

"""
    GetMaxDimensionSize(this) -> max_dim_size
"""
function GetMaxDimensionSize(this::fates_parameters_type)
    max_dim_size = 0
    for p in 1:num_params(this)
        for d in 1:max_dimensions
            max_dim_size = max(max_dim_size, this.parameters[p].dimension_sizes[d])
        end
    end
    return max_dim_size
end

# =====================================================================================
# Data storage (Fortran generic SetData -> SetData{Scalar,1D,2D})
# =====================================================================================

"""
    SetDataScalar!(this, index, data)
"""
function SetDataScalar!(this::fates_parameters_type, index::Integer, data::Real)
    this.parameters[index].data = reshape([Float64(data)], 1, 1)
    return nothing
end

"""
    SetData1D!(this, index, data)
"""
function SetData1D!(this::fates_parameters_type, index::Integer, data::AbstractVector)
    size_dim_1 = this.parameters[index].dimension_sizes[1]
    if length(data) != size_dim_1
        @warn "setdata1d : $(this.parameters[index].name) size inconsistent. " *
              "expected=$(length(data)) received=$size_dim_1"
        fates_endrun("size error setting 1d parameter.")
    end
    m = Matrix{Float64}(undef, size_dim_1, 1)
    m[:, 1] .= data
    this.parameters[index].data = m
    return nothing
end

"""
    SetData2D!(this, index, data)
"""
function SetData2D!(this::fates_parameters_type, index::Integer, data::AbstractMatrix)
    this.parameters[index].data = Matrix{Float64}(data)
    return nothing
end
