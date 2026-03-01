# ==========================================================================
# Ported from: src/biogeophys/LakeStateType.F90
# Lake state data type allocation and initialization
# ==========================================================================

"""
    LakeStateData

Lake state data structure. Holds all lake state variables at the column
and patch levels, including ice fractions, eddy conductivity, friction
velocity, aerodynamic resistance, and surface data fields.

Ported from `lakestate_type` in `LakeStateType.F90`.
"""
Base.@kwdef mutable struct LakeStateData
    # --- Time constant variables (from surface data) ---
    lakefetch_col        ::Vector{Float64} = Float64[]   # col lake fetch from surface data (m)
    etal_col             ::Vector{Float64} = Float64[]   # col lake extinction coefficient from surface data (1/m)

    # --- Time varying variables ---
    lake_raw_col         ::Vector{Float64} = Float64[]   # col aerodynamic resistance for moisture (s/m)
    ks_col               ::Vector{Float64} = Float64[]   # col coefficient for calculation of decay of eddy diffusivity with depth
    ws_col               ::Vector{Float64} = Float64[]   # col surface friction velocity (m/s)
    ust_lake_col         ::Vector{Float64} = Float64[]   # col friction velocity (m/s)
    betaprime_col        ::Vector{Float64} = Float64[]   # col effective beta: sabg_lyr(p,jtop) for snow layers, beta otherwise
    savedtke1_col        ::Vector{Float64} = Float64[]   # col top level eddy conductivity from previous timestep (W/mK)
    lake_icefrac_col     ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col mass fraction of lake layer that is frozen (ncols, nlevlak)
    lake_icefracsurf_col ::Vector{Float64} = Float64[]   # col mass fraction of surface lake layer that is frozen
    lake_icethick_col    ::Vector{Float64} = Float64[]   # col ice thickness (m) (integrated if lakepuddling)
    lakeresist_col       ::Vector{Float64} = Float64[]   # col [s/m] (Needed for calc. of grnd_ch4_cond)
    ram1_lake_patch      ::Vector{Float64} = Float64[]   # patch aerodynamical resistance (s/m)
end

"""
    lakestate_init!(ls::LakeStateData, nc::Int, np::Int)

Allocate and initialize all fields of a `LakeStateData` instance for
`nc` columns and `np` patches. Column-level fields are sized `nc`,
patch-level fields are sized `np`. The 2D `lake_icefrac_col` is
sized `(nc, nlevlak)`.

Most fields are initialized to `NaN`, except `savedtke1_col` and
`ust_lake_col` which are initialized to `SPVAL` (matching Fortran).

Ported from `lakestate_type%InitAllocate` in `LakeStateType.F90`.
"""
function lakestate_init!(ls::LakeStateData, nc::Int, np::Int)
    nlevlak = varpar.nlevlak

    # Time constant variables
    ls.lakefetch_col        = fill(NaN, nc)
    ls.etal_col             = fill(NaN, nc)

    # Time varying variables — column level
    ls.lake_raw_col         = fill(NaN, nc)
    ls.ks_col               = fill(NaN, nc)
    ls.ws_col               = fill(NaN, nc)
    ls.ust_lake_col         = fill(SPVAL, nc)
    ls.betaprime_col        = fill(NaN, nc)
    ls.savedtke1_col        = fill(SPVAL, nc)
    ls.lake_icefrac_col     = fill(NaN, nc, nlevlak)
    ls.lake_icefracsurf_col = fill(NaN, nc)
    ls.lake_icethick_col    = fill(NaN, nc)
    ls.lakeresist_col       = fill(NaN, nc)

    # Time varying variables — patch level
    ls.ram1_lake_patch      = fill(NaN, np)

    return nothing
end

"""
    lakestate_clean!(ls::LakeStateData)

Deallocate (reset to empty) all fields of a `LakeStateData` instance.
"""
function lakestate_clean!(ls::LakeStateData)
    ls.lakefetch_col        = Float64[]
    ls.etal_col             = Float64[]
    ls.lake_raw_col         = Float64[]
    ls.ks_col               = Float64[]
    ls.ws_col               = Float64[]
    ls.ust_lake_col         = Float64[]
    ls.betaprime_col        = Float64[]
    ls.savedtke1_col        = Float64[]
    ls.lake_icefrac_col     = Matrix{Float64}(undef, 0, 0)
    ls.lake_icefracsurf_col = Float64[]
    ls.lake_icethick_col    = Float64[]
    ls.lakeresist_col       = Float64[]
    ls.ram1_lake_patch      = Float64[]

    return nothing
end

"""
    lakestate_init_cold!(ls::LakeStateData, bounds_col::UnitRange{Int};
                          mask_lake::Union{BitVector,Nothing}=nothing)

Initialize cold-start conditions for lake state variables.
Sets lake ice fraction to zero, savedtke1 to `TKWAT`, and
ust_lake to 0.1 for lake columns.

If `mask_lake` is provided, only columns where `mask_lake[c]` is `true`
are initialized. Otherwise, all columns in `bounds_col` are initialized.

Ported from `lakestate_type%InitCold` in `LakeStateType.F90`.
"""
function lakestate_init_cold!(ls::LakeStateData, bounds_col::UnitRange{Int};
                               mask_lake::Union{BitVector,Nothing} = nothing)
    nlevlak = varpar.nlevlak

    for c in bounds_col
        if mask_lake !== nothing
            mask_lake[c] || continue
        end

        # Initialize with no ice
        for j in 1:nlevlak
            ls.lake_icefrac_col[c, j] = 0.0
        end

        # Set top eddy conductivity from previous timestep
        ls.savedtke1_col[c] = TKWAT

        # Set column friction velocity
        ls.ust_lake_col[c] = 0.1
    end

    return nothing
end

"""
    lakestate_init_history!(ls::LakeStateData, bounds_col::UnitRange{Int},
                             bounds_patch::UnitRange{Int})

Register lake state fields for history file output.

Ported from `lakestate_type%InitHistory` in `LakeStateType.F90`.
Requires history infrastructure (histFileMod) — stub until that module is ported.
"""
function lakestate_init_history!(ls::LakeStateData,
                                  bounds_col::UnitRange{Int},
                                  bounds_patch::UnitRange{Int})
    # Stub: history field registration will be added when histFileMod is ported.
    # Fields that would be registered:
    #   LAKEICEFRAC (2D, levlak), LAKEICEFRAC_SURF, LAKEICETHICK,
    #   TKE1, RAM_LAKE, UST_LAKE

    # Set SPVAL for history averaging (matching Fortran InitHistory)
    for c in bounds_col
        for j in 1:varpar.nlevlak
            ls.lake_icefrac_col[c, j] = SPVAL
        end
        ls.lake_icefracsurf_col[c] = SPVAL
        ls.lake_icethick_col[c]    = SPVAL
        ls.savedtke1_col[c]        = SPVAL
        ls.ust_lake_col[c]         = SPVAL
    end
    for p in bounds_patch
        ls.ram1_lake_patch[p] = SPVAL
    end

    return nothing
end

"""
    lakestate_restart!(ls::LakeStateData, bounds_col::UnitRange{Int};
                        flag::String="read")

Read/write lake state from/to restart file.

Ported from `lakestate_type%Restart` in `LakeStateType.F90`.
Requires NetCDF/restart infrastructure — stub until that module is ported.
"""
function lakestate_restart!(ls::LakeStateData,
                             bounds_col::UnitRange{Int};
                             flag::String = "read")
    # Stub: restart variable I/O will be added when restUtilMod/ncdio_pio is ported.
    # Variables that would be read/written:
    #   LAKE_ICEFRAC (2D, levlak), SAVEDTKE1, USTLAKE
    return nothing
end
