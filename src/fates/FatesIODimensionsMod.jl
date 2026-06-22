# FatesIODimensionsMod.jl
# Julia port of FATES src/fates/main/FatesIODimensionsMod.F90
#
# IO dimension descriptors: dimension-name string constants, a bounds type
# carrying begin/end for every FATES IO dimension, and a per-dimension
# descriptor type with per-thread (clump) bounds.
# Deps: FatesConstantsMod (fates_short_string_length).

# ---------------------------------------------------------------------------
# Dimension name constants (must match histFileMod / clm_varcon on the HLM side)
# ---------------------------------------------------------------------------
const dimname_levcapf = "fates_levcapf"
const dimname_levcacls = "fates_levcacls"
const dimname_cohort = "cohort"
const dimname_column = "column"
const dimname_levsoil = "levsoi"
const dimname_levscag = "fates_levscag"
const dimname_levscagpft = "fates_levscagpf"
const dimname_levagepft = "fates_levagepft"
const dimname_levscpf = "fates_levscpf"
const dimname_levscls = "fates_levscls"
const dimname_levpft = "fates_levpft"
const dimname_levage = "fates_levage"
const dimname_levheight = "fates_levheight"
const dimname_levfuel = "fates_levfuel"
const dimname_levcwdsc = "fates_levcwdsc"
const dimname_levcan = "fates_levcan"
const dimname_levcnlf = "fates_levcnlf"
const dimname_levcnlfpft = "fates_levcnlfpf"
const dimname_levclscpf = "fates_levclscpf"
const dimname_levcdsc = "fates_levcdsc"
const dimname_levcdpf = "fates_levcdpf"
const dimname_levcdam = "fates_levcdam"
const dimname_levagefuel = "fates_levagefuel"
const dimname_levelem = "fates_levelem"
const dimname_levelpft = "fates_levelpft"
const dimname_levelcwd = "fates_levelcwd"
const dimname_levelage = "fates_levelage"
const dimname_levlanduse = "fates_levlanduse"
const dimname_levlulu = "fates_levlulu"
const dimname_levlupft = "fates_levlupft"

# ---------------------------------------------------------------------------
# fates_bounds_type — begin/end for every IO dimension. Default 0.
# ---------------------------------------------------------------------------
Base.@kwdef mutable struct fates_bounds_type
    cohort_begin::Int = 0
    cohort_end::Int = 0
    column_begin::Int = 0          # FATES has no "column"; we call this a "site"
    column_end::Int = 0
    soil_begin::Int = 0
    soil_end::Int = 0
    sizeage_class_begin::Int = 0
    sizeage_class_end::Int = 0
    sizeagepft_class_begin::Int = 0
    sizeagepft_class_end::Int = 0
    agepft_class_begin::Int = 0
    agepft_class_end::Int = 0
    sizepft_class_begin::Int = 0
    sizepft_class_end::Int = 0
    coagepf_class_begin::Int = 0
    coagepf_class_end::Int = 0
    size_class_begin::Int = 0
    size_class_end::Int = 0
    coage_class_begin::Int = 0
    coage_class_end::Int = 0
    pft_class_begin::Int = 0
    pft_class_end::Int = 0
    age_class_begin::Int = 0
    age_class_end::Int = 0
    height_begin::Int = 0
    height_end::Int = 0
    fuel_begin::Int = 0
    fuel_end::Int = 0
    cwdsc_begin::Int = 0
    cwdsc_end::Int = 0
    can_begin::Int = 0
    can_end::Int = 0
    cnlf_begin::Int = 0
    cnlf_end::Int = 0
    cnlfpft_begin::Int = 0
    cnlfpft_end::Int = 0
    cdsc_begin::Int = 0
    cdsc_end::Int = 0
    cdpf_begin::Int = 0
    cdpf_end::Int = 0
    cdam_begin::Int = 0
    cdam_end::Int = 0
    elem_begin::Int = 0
    elem_end::Int = 0
    elpft_begin::Int = 0
    elpft_end::Int = 0
    elcwd_begin::Int = 0
    elcwd_end::Int = 0
    elage_begin::Int = 0
    elage_end::Int = 0
    agefuel_begin::Int = 0
    agefuel_end::Int = 0
    clscpf_begin::Int = 0
    clscpf_end::Int = 0
    landuse_begin::Int = 0
    landuse_end::Int = 0
    lulu_begin::Int = 0
    lulu_end::Int = 0
    lupft_begin::Int = 0
    lupft_end::Int = 0
end

# ---------------------------------------------------------------------------
# fates_io_dimension_type — one dimension descriptor, with per-thread bounds.
# ---------------------------------------------------------------------------
Base.@kwdef mutable struct fates_io_dimension_type
    name::String = ""
    lower_bound::Int = 0
    upper_bound::Int = 0
    clump_lower_bound::Vector{Int} = Int[]  # lower bound of thread's HIO array portion
    clump_upper_bound::Vector{Int} = Int[]  # upper bound of thread's HIO array portion
end

# =====================================================================================

"""
    Init!(this::fates_io_dimension_type, name, num_threads, lower_bound, upper_bound)

Initialize a dimension descriptor and allocate its per-thread (clump) bound
arrays, filled with -1.
"""
function Init!(this::fates_io_dimension_type, name::AbstractString,
               num_threads::Integer, lower_bound::Integer, upper_bound::Integer)
    this.name = strip(name)
    this.lower_bound = lower_bound
    this.upper_bound = upper_bound

    this.clump_lower_bound = fill(-1, num_threads)
    this.clump_upper_bound = fill(-1, num_threads)
    return nothing
end

# =====================================================================================

"""
    SetThreadBounds!(this::fates_io_dimension_type, thread_index, lower_bound, upper_bound)

Set the HIO-array bounds for one thread (1-based `thread_index`).
"""
function SetThreadBounds!(this::fates_io_dimension_type, thread_index::Integer,
                          lower_bound::Integer, upper_bound::Integer)
    this.clump_lower_bound[thread_index] = lower_bound
    this.clump_upper_bound[thread_index] = upper_bound
    return nothing
end
