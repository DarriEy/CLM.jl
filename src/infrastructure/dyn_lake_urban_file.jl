# dyn_lake_urban_file.jl — Transient lake & urban land-use data readers.
#
# Julia port of two Fortran modules from CLM/CTSM's src/dyn_subgrid/:
#   - dynlakeFileMod.F90   -> dynlake_init  / dynlake_interp   (PCT_LAKE time series)
#   - dynurbanFileMod.F90  -> dynurban_init / dynurban_interp  (PCT_URBAN[numurbl])
#
# Both read a time series of landunit percent-cover from a transient landuse
# (landuse.timeseries) NetCDF file and, each year, snap the corresponding
# landunit weight(s) to the new value (uninterpolated — jumps at the year
# boundary). They reuse the DynFile / DynVarTimeUninterp framework from
# src/infrastructure/dyn_file_io.jl and write the updated landunit weights via
# set_landunit_weight! (mirroring subgridWeightsMod :: set_landunit_weight).
#
# Fortran subroutine / variable names are preserved for traceability. The
# Fortran module-private singletons (dynlake_file/wtlake, dynurban_file/wturban)
# become explicit @kwdef structs returned by *_init and threaded into *_interp.
# This module is standalone (not wired into the driver, not added to any
# CLMInstances / ForwardDiff-dual-copied struct).

# grlnd: name of the lndgrid spatial dimension (clm_varcon :: grlnd). Used only
# as the dyn-var dim1name metadata (the Julia reader slices by position).
const grlnd = "lndgrid"

# ======================================================================
# set_landunit_weight!  (subgridWeightsMod :: set_landunit_weight)
# ======================================================================
"""
    set_landunit_weight!(grc, lun, g, ltype, weight)

Set the subgrid weight of a given landunit type `ltype` on a single grid cell
`g`. Mirrors Fortran `subgridWeightsMod :: set_landunit_weight`.

If a landunit of type `ltype` exists on gridcell `g`, its `wtgcell` is set to
`weight`. If no such landunit exists, this is a no-op for `weight == 0`, but an
attempt to assign a non-zero weight to a non-existent landunit is an error.
"""
function set_landunit_weight!(grc::GridcellData, lun::LandunitData,
                              g::Int, ltype::Int, weight::Real)
    l = grc.landunit_indices[ltype, g]
    if l != ISPVAL
        lun.wtgcell[l] = weight
    elseif weight > 0.0
        error("set_landunit_weight! ERROR: Attempt to assign non-zero weight to " *
              "a non-existent landunit: g, l, ltype, weight = $g, $l, $ltype, $weight")
    end
    return nothing
end

# ======================================================================
# dynlakeFileMod
# ======================================================================
const lake_varname = "PCT_LAKE"

"""
    DynLakeFile

Holds the transient-lake reader state (mirrors the Fortran module-private
singletons `dynlake_file` + `wtlake`): the open `DynFile` and the
`DynVarTimeUninterp` for the lake landunit weight (PCT_LAKE).
"""
Base.@kwdef mutable struct DynLakeFile
    dynlake_file::DynFile                  # file with transient lake data
    wtlake::DynVarTimeUninterp             # weight of the lake landunit
end

"""
    dynlake_init(dynlake_filename, begg, endg; current_year=nothing) -> DynLakeFile

Initialize the dataset containing transient lake info (position it to the time
samples that bound the initial model year). Mirrors Fortran `dynlake_init`.

`begg`/`endg` are the proc-level gridcell bounds; `num_points = endg-begg+1` is
the spatial size. The PCT_LAKE field is read with `conversion_factor=100`
(percent -> fraction) and `do_check_sums_equal_1=false`, exactly as in Fortran.
Reading the year from the START of the timestep (`YEAR_POSITION_START_OF_TIMESTEP`)
matches the Fortran, so lake areas update starting after the year boundary.
"""
function dynlake_init(dynlake_filename::String, begg::Int, endg::Int;
                      current_year::Union{Int,Nothing} = nothing)
    # Get the year from the START of the timestep (consistent with glacier updates).
    dynlake_file = dyn_file_open(dynlake_filename, YEAR_POSITION_START_OF_TIMESTEP;
                                 current_year = current_year)

    # read data PCT_LAKE
    num_points = (endg - begg + 1)
    wtlake = dyn_var_time_uninterp(dynlake_file, lake_varname,
                                   grlnd, 100.0, false, [num_points])

    return DynLakeFile(dynlake_file = dynlake_file, wtlake = wtlake)
end

"""
    dynlake_interp(dl, grc, lun, begg, endg, year)

Get lake cover for the current model `year`, when needed, and set the lake
landunit weight (`lun.wtgcell`) for the deep-lake (`ISTDLAK`) landunit of each
gridcell. Mirrors Fortran `dynlake_interp`.

Lake cover jumps to its new value at the start of the year (uninterpolated).
"""
function dynlake_interp(dl::DynLakeFile, grc::GridcellData, lun::LandunitData,
                        begg::Int, endg::Int, year::Int)
    set_current_year!(dl.dynlake_file, year)

    # Set new landunit area.
    wtlake_cur = get_current_data_1d(dl.wtlake)   # indexed 1:num_points (begg..endg)
    for g in begg:endg
        set_landunit_weight!(grc, lun, g, ISTDLAK, wtlake_cur[g - begg + 1])
    end

    return nothing
end

# ======================================================================
# dynurbanFileMod
# ======================================================================
const urban_varname = "PCT_URBAN"

"""
    DynUrbanFile

Holds the transient-urban reader state (mirrors the Fortran module-private
singletons `dynurban_file` + `wturban`): the open `DynFile` and the 2-d
`DynVarTimeUninterp` for the per-density-class urban landunit weights
(PCT_URBAN[num_points, numurbl]).
"""
Base.@kwdef mutable struct DynUrbanFile
    dynurban_file::DynFile                 # file with transient urban data
    wturban::DynVarTimeUninterp            # weight of the urban landunits (2-d)
end

"""
    dynurban_init(dynurban_filename, begg, endg; current_year=nothing) -> DynUrbanFile

Initialize the dataset containing transient urban info (position it to the time
samples that bound the initial model year). Mirrors Fortran `dynurban_init`.

Checks that the file's `numurbl` dimension equals the model `NUMURBL` (Fortran
`check_dim_size`). The PCT_URBAN field is a 2-d variable shaped
`[num_points, NUMURBL]`, read with `conversion_factor=100` (percent -> fraction)
and `do_check_sums_equal_1=false`, exactly as in Fortran.
"""
function dynurban_init(dynurban_filename::String, begg::Int, endg::Int;
                       current_year::Union{Int,Nothing} = nothing)
    # Get the year from the START of the timestep (consistent with glacier updates).
    dynurban_file = dyn_file_open(dynurban_filename, YEAR_POSITION_START_OF_TIMESTEP;
                                  current_year = current_year)

    # check_dim_size(dynurban_file, 'numurbl', numurbl)
    ds = dynurban_file.ds
    if haskey(ds.dim, "numurbl")
        nfile = ds.dim["numurbl"]
        nfile == NUMURBL || error("dynurban_init ERROR: dimension numurbl on file " *
                                  "= $nfile, expected $NUMURBL")
    end

    # read data PCT_URBAN
    num_points = (endg - begg + 1)
    wturb_shape = [num_points, NUMURBL]
    wturban = dyn_var_time_uninterp(dynurban_file, urban_varname,
                                    grlnd, 100.0, false, wturb_shape)

    return DynUrbanFile(dynurban_file = dynurban_file, wturban = wturban)
end

"""
    dynurban_interp(du, grc, lun, begg, endg, year; urbinp=nothing, caller="dynurban_interp")

Get urban cover for the current model `year`, when needed, and set the urban
landunit weights (`lun.wtgcell`) for the three urban density classes
(`ISTURB_TBD`, `ISTURB_HD`, `ISTURB_MD` — columns 1,2,3 of PCT_URBAN) of each
gridcell. Mirrors Fortran `dynurban_interp`.

Urban cover jumps to its new value at the start of the year (uninterpolated).

If `urbinp` (an `UrbanInputData`) is supplied, the urban-data consistency check
`check_urban` is run on the current weights (mirroring the Fortran
`CheckUrban`); pass `nothing` to skip it (e.g. standalone usage without urban
parameters loaded).
"""
function dynurban_interp(du::DynUrbanFile, grc::GridcellData, lun::LandunitData,
                         begg::Int, endg::Int, year::Int;
                         urbinp = nothing, caller::String = "dynurban_interp")
    set_current_year!(du.dynurban_file, year)

    # Set new landunit areas. wturban_cur is [num_points, NUMURBL].
    wturban_cur = get_current_data_2d(du.wturban)
    for g in begg:endg
        gi = g - begg + 1
        set_landunit_weight!(grc, lun, g, ISTURB_TBD, wturban_cur[gi, 1])
        set_landunit_weight!(grc, lun, g, ISTURB_HD,  wturban_cur[gi, 2])
        set_landunit_weight!(grc, lun, g, ISTURB_MD,  wturban_cur[gi, 3])
    end

    # Check that urban data is valid (Fortran CheckUrban).
    if urbinp !== nothing
        check_urban(urbinp, wturban_cur, 1:(endg - begg + 1); caller = caller)
    end

    return nothing
end
