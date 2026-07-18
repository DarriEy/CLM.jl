# dyn_pft_crop_file.jl — Transient natural-PFT and crop land-use data readers.
#
# Julia port of two tightly-coupled Fortran modules from CLM/CTSM's
# src/dyn_subgrid/:
#   - dynpftFileMod.F90   -> dynpft_init / dynpft_interp   (PCT_NAT_PFT weights)
#   - dyncropFileMod.F90  -> dyncrop_init / dyncrop_interp (PCT_CROP / PCT_CFT /
#                                                           FERTNITRO_CFT)
#
# Both modules read transient (year-by-year) areal weights from the
# flanduse_timeseries file and, each year, push the new weights into the subgrid
# weight arrays:
#   - dynpft_interp sets patch%wtcol for natural-veg (istsoil) patches.
#   - dyncrop_interp sets lun%wtgcell for the crop landunit (via set_landunit_weight)
#     and col%wtlunit for each crop column, plus the annual crop fertilizer.
#
# File reading uses the DynFile / DynVarTimeUninterp framework from dyn_file_io.jl
# (which this file is included AFTER). The variables are NOT interpolated in time:
# they snap to each year's value at the year boundary, exactly as in the Fortran
# (where the comment notes that switching wtpatch / wtcrop / wtcft to a
# dyn_var_time_interp_type is "all you need to do" to enable interpolation).
#
# Fortran subroutine / variable names are preserved for traceability. The Fortran
# module-private singletons (dynpft_file/wtpatch and dyncrop_file/wtcrop/wtcft/
# fertcft) become explicit mutable @kwdef state structs that the caller owns and
# passes back into the interp routines. The subgrid types (grc/lun/col/pch) and
# bounds are passed explicitly rather than read from module globals. None of these
# structs are added to CLMInstances or any ForwardDiff-dual-copied struct.
#
# The namelist-only pieces of the Fortran (dynpft_read_consistency_nl + MPI bcast)
# are not ported; the consistency check instead takes a `check_dynpft_consistency`
# flag (defaulting to true, the Fortran default) directly.

# ======================================================================
# dynpftFileMod :: PCT_NAT_PFT (transient natural-veg PFT weights)
# ======================================================================

# Names of variables on file (dynpftFileMod parameters).
const DYNPFT_VARNAME = "PCT_NAT_PFT"

"""
    DynpftState

Holds the per-run state for the transient natural-PFT reader, mirroring the
Fortran `dynpftFileMod` module-private singletons:

  - `dynpft_file::DynFile`            — the flanduse_timeseries file handle
  - `wtpatch::DynVarTimeUninterp`     — weight of each PFT relative to the natural
                                        veg landunit (PCT_NAT_PFT), as a 2-d
                                        (natpft, gridcell) variable

The data_shape is `[natpft_size, ngridcells]` — i.e. the natpft dimension is
FIRST. This is required by the dyn_file_io `do_check_sums_equal_1` machinery,
which checks that the sum over data_shape dimension 1 equals 1 (here: each grid
cell's PFT weights sum to 1). The on-file `PCT_NAT_PFT` is therefore read with
natpft as its leading spatial dimension.

The caller owns this object: `dynpft_init` returns it and each `dynpft_interp`
call takes it back. NOT added to any dual-copied / CLMInstances struct.
"""
Base.@kwdef mutable struct DynpftState
    dynpft_file::DynFile
    wtpatch::DynVarTimeUninterp
end

"""
    dynpft_init(dynpft_filename; ngridcells, natpft_size, current_year,
                wt_nat_patch=nothing, check_dynpft_consistency=true,
                tol=1e-6) -> DynpftState

Initialize the dynamic-PFT dataset (Fortran `dynpft_init`).

Opens `dynpft_filename` via `dyn_file_open` (positioned to the year from the START
of the timestep, `YEAR_POSITION_START_OF_TIMESTEP`), checks the `natpft` dimension
size, optionally runs the surface-dataset consistency check, then constructs the
`PCT_NAT_PFT` `DynVarTimeUninterp` variable (converted from percent to weight via
`conversion_factor = 100`, with `do_check_sums_equal_1 = true`).

`ngridcells` is the number of grid cells (Fortran `endg - begg + 1`); the wtpatch
data shape is `[ngridcells, natpft_size]`.

If `check_dynpft_consistency` is true and `wt_nat_patch` (the surface-dataset
natural-PFT weights, shape `(natpft_size, ngridcells)`, already as fractions) is
provided, `dynpft_check_consistency` verifies that the file's first time slice
agrees with the surface dataset to within `tol`.
"""
function dynpft_init(dynpft_filename::String;
                     ngridcells::Int,
                     natpft_size::Int,
                     current_year::Int,
                     wt_nat_patch::Union{AbstractMatrix{<:Real},Nothing} = nothing,
                     check_dynpft_consistency::Bool = true,
                     tol::Float64 = 1.0e-6)

    # Get the year from the START of the timestep; this way, we'll update PFT areas
    # starting after the year boundary.
    dynpft_file = dyn_file_open(dynpft_filename, YEAR_POSITION_START_OF_TIMESTEP;
                                current_year = current_year)

    # Consistency check: the file's 'natpft' dimension must equal natpft_size.
    _dynpft_check_dim_size(dynpft_file, "natpft", natpft_size)

    # Consistency check against the surface dataset.
    if check_dynpft_consistency && wt_nat_patch !== nothing
        dynpft_check_consistency(dynpft_file, wt_nat_patch; tol = tol)
    end

    # read data PCT_NAT_PFT corresponding to correct year.
    #
    # Note (matching Fortran): to interpolate rather than jump at Jan 1, change
    # wtpatch to a dyn_var_time_interp_type — that's all that is needed.
    wtpatch_shape = [natpft_size, ngridcells]
    wtpatch = dyn_var_time_uninterp(dynpft_file, DYNPFT_VARNAME, GRLND, 100.0,
                                    true, wtpatch_shape)

    return DynpftState(dynpft_file = dynpft_file, wtpatch = wtpatch)
end

# check_dim_size (ncdio_pio): verify that a file dimension equals an expected size.
function _dynpft_check_dim_size(df::DynFile, dimname::String, expected::Int)
    ds = df.ds
    haskey(ds.dim, dimname) ||
        error("check_dim_size ERROR: dimension '$dimname' not on file $(df.filename)")
    actual = ds.dim[dimname]
    actual == expected ||
        error("check_dim_size ERROR: dimension '$dimname' = $actual on file, " *
              "expected $expected")
    return nothing
end

"""
    dynpft_check_consistency(dynpft_file::DynFile, wt_nat_patch; tol=1e-6)

Check consistency between the dynpft file and the surface dataset (Fortran
`dynpft_check_consistency`).

Reads the first time slice of `PCT_NAT_PFT`, converts from percent to weight
(`/ 100`), and compares each grid cell's PFT weights with `wt_nat_patch` (the
surface-dataset weights, shape `(natpft, ngridcells)`, already as fractions).
Throws (Fortran `endrun`) if any grid cell differs by more than `tol`.
"""
function dynpft_check_consistency(dynpft_file::DynFile,
                                  wt_nat_patch::AbstractMatrix{<:Real};
                                  tol::Float64 = 1.0e-6)
    ds = dynpft_file.ds
    haskey(ds, DYNPFT_VARNAME) ||
        error("dynpft_check_consistency ERROR: $(DYNPFT_VARNAME) NOT on " *
              "landuse_timeseries file")

    # Read first time slice (nt = 1) and convert from PCT to weight on grid cell.
    # The target shape comes from the surface dataset, so real CTSM files (which
    # carry the gridcell axis as an (lsmlon, lsmlat) PAIR) fold correctly here too.
    ncvar = ds[DYNPFT_VARNAME]
    shape = (size(wt_nat_patch, 1), size(wt_nat_patch, 2))
    rawslice = _read_time_slice(ncvar, 1, 2, shape)       # (natpft, ngridcells)
    wtpatch_time1 = Float64.(replace(rawslice, missing => NaN)) ./ 100.0

    ngridcells = size(wtpatch_time1, 2)
    size(wt_nat_patch, 2) == ngridcells ||
        error("dynpft_check_consistency ERROR: wt_nat_patch gridcell dim " *
              "$(size(wt_nat_patch, 2)) != file $ngridcells")

    for g in 1:ngridcells
        if any(abs.(wtpatch_time1[:, g] .- wt_nat_patch[:, g]) .> tol)
            error("dynpft_check_consistency: mismatch between PCT_NAT_PFT at " *
                  "initial time and that obtained from surface dataset at " *
                  "gridcell g=$g.\n" *
                  "On landuse_timeseries file: $(wtpatch_time1[:, g])\n" *
                  "On surface dataset: $(wt_nat_patch[:, g])\n" *
                  "Confirm that the year of your surface dataset corresponds to " *
                  "the first year of your landuse_timeseries file, or bypass this " *
                  "check by passing check_dynpft_consistency=false.")
        end
    end
    return nothing
end

"""
    dynpft_interp!(state::DynpftState, bounds, pch, lun; year=nothing)

Get PFT weights for the current model time and push them into `pch.wtcol`
(Fortran `dynpft_interp`).

Advances the file's `time_info` to the current year (if `year` is given,
`set_current_year!` is called; otherwise the year already set on the file is
used), reads the current `PCT_NAT_PFT` weights, then for every patch on a soil
(istsoil) landunit sets `pch.wtcol[p] = wtpatch_cur[g, m]`, where `m` is the
patch's vegetation type index (`pch.itype[p]`) on the natural-veg axis.

Note (matching Fortran): PFT weights jump to their new value at the start of the
year. This assumes each landunit has only one column.
"""
function dynpft_interp!(state::DynpftState, bounds::BoundsType,
                        pch::PatchData, lun::LandunitData;
                        year::Union{Int,Nothing} = nothing)
    if year !== nothing
        set_current_year!(state.dynpft_file, year)
    end

    # Get pft weights for this time step: (natpft_size, ngridcells).
    wtpatch_cur = get_current_data_2d(state.wtpatch)

    for p in bounds.begp:bounds.endp
        g = pch.gridcell[p]
        l = pch.landunit[p]

        if lun.itype[l] == ISTSOIL
            # Patch vegetation type on the natural-veg axis. The Fortran patch
            # itype is 0-based on natpft (natpft_lb..natpft_ub); the 1-based
            # natpft data index here is m = itype + 1.
            m = pch.itype[p] + 1

            # Note: assumes all patches share a single column.
            pch.wtcol[p] = wtpatch_cur[m, g]
        end
    end

    return nothing
end


# ======================================================================
# dyncropFileMod :: PCT_CROP / PCT_CFT / FERTNITRO_CFT (transient crops)
# ======================================================================

# Names of variables on file (dyncropFileMod parameters).
const DYNCROP_VARNAME = "PCT_CROP"
const DYNCFT_VARNAME  = "PCT_CFT"
const DYNFERT_VARNAME = "FERTNITRO_CFT"

"""
    DyncropState

Holds the per-run state for the transient-crop reader, mirroring the Fortran
`dyncropFileMod` module-private singletons:

  - `dyncrop_file::DynFile`           — the flanduse_timeseries file handle
  - `wtcrop::DynVarTimeUninterp`      — weight of the crop landunit (PCT_CROP), a
                                        1-d (gridcell) variable
  - `wtcft::DynVarTimeUninterp`       — weight of each CFT relative to the crop
                                        landunit (PCT_CFT), a 2-d (cft, gridcell)
                                        variable
  - `fertcft::DynVarTimeUninterp`     — fertilizer of each CFT (FERTNITRO_CFT),
                                        a 2-d (cft, gridcell) variable
                                        (allow_nodata = true)

As with the natural-PFT reader, the CFT dimension is FIRST in `data_shape`
(`[cft_size, ngridcells]`) so the dyn_file_io sum-to-1 check (over data_shape
dimension 1) verifies that each grid cell's CFT weights sum to 1.

The caller owns this object. NOT added to any dual-copied / CLMInstances struct.
"""
Base.@kwdef mutable struct DyncropState
    dyncrop_file::DynFile
    wtcrop::DynVarTimeUninterp
    wtcft::DynVarTimeUninterp
    fertcft::DynVarTimeUninterp
end

"""
    dyncrop_init(dyncrop_filename; ngridcells, cft_size, current_year) -> DyncropState

Initialize the transient-crop dataset (Fortran `dyncrop_init`).

Opens `dyncrop_filename` (positioned to the START-of-timestep year), checks the
`cft` dimension size, and constructs the three `DynVarTimeUninterp` variables:

  - `PCT_CROP`       — 1-d, conversion 100, no sum check
  - `PCT_CFT`        — 2-d `[cft_size, ngridcells]`, conversion 100, sum check
  - `FERTNITRO_CFT`  — 2-d `[cft_size, ngridcells]`, conversion 1, no sum check,
                       `allow_nodata = true`

`ngridcells` is `endg - begg + 1`.
"""
function dyncrop_init(dyncrop_filename::String;
                      ngridcells::Int,
                      cft_size::Int,
                      current_year::Int)

    # Get the year from the START of the timestep.
    dyncrop_file = dyn_file_open(dyncrop_filename, YEAR_POSITION_START_OF_TIMESTEP;
                                 current_year = current_year)
    _dynpft_check_dim_size(dyncrop_file, "cft", cft_size)

    num_points = ngridcells
    wtcrop = dyn_var_time_uninterp(dyncrop_file, DYNCROP_VARNAME, GRLND, 100.0,
                                   false, [num_points])
    wtcft_shape = [cft_size, num_points]
    wtcft = dyn_var_time_uninterp(dyncrop_file, DYNCFT_VARNAME, GRLND, 100.0,
                                  true, wtcft_shape)
    fertcft_shape = [cft_size, num_points]
    fertcft = dyn_var_time_uninterp(dyncrop_file, DYNFERT_VARNAME, GRLND, 1.0,
                                    false, fertcft_shape; allow_nodata = true)

    return DyncropState(dyncrop_file = dyncrop_file, wtcrop = wtcrop,
                        wtcft = wtcft, fertcft = fertcft)
end

# set_landunit_weight! is defined canonically in dyn_landunit_area.jl (loaded
# earlier); the transient readers here reuse that single definition.

"""
    dyncrop_interp!(state::DyncropState, bounds, grc, lun, col, pch;
                    year=nothing, use_crop=false, crop_inst=nothing,
                    collapse_crops=false)

Get crop cover for the current model time and push it into the subgrid weights
(Fortran `dyncrop_interp`).

Sets `lun.wtgcell` for crop landunits (via `set_landunit_weight!`) and
`col.wtlunit` for each crop column. With `use_crop = true` and a `crop_inst`
(`CropData`) provided, also sets `crop_inst.fertnitro_patch[p]` from the current
fertilizer data.

If `collapse_crops` is true, `collapse_crop_types!` is applied to the CFT weights
and fertilizer before they are assigned (this requires the global `pftcon` /
`varctl` state to be initialized; it is off by default so the reader can be used
standalone). The default (false) assigns the raw file CFT weights directly.

The error-check via `col_set` ensures each crop column holds a single CFT.
"""
function dyncrop_interp!(state::DyncropState, bounds::BoundsType,
                         grc::GridcellData, lun::LandunitData,
                         col::ColumnData, pch::PatchData;
                         year::Union{Int,Nothing} = nothing,
                         use_crop::Bool = false,
                         crop_inst = nothing,
                         collapse_crops::Bool = false)
    if year !== nothing
        set_current_year!(state.dyncrop_file, year)
    end

    # Set new landunit area.
    wtcrop_cur = get_current_data_1d(state.wtcrop)             # (ngridcells,)
    for g in bounds.begg:bounds.endg
        set_landunit_weight!(grc, lun, g, ISTCROP, wtcrop_cur[g])
    end

    # Set new CFT weights. Assumes each crop is on its own column.
    wtcft_cur   = get_current_data_2d(state.wtcft)            # (cft, ngridcells)
    fertcft_cur = get_current_data_2d(state.fertcft)          # (cft, ngridcells)

    # Collapse crop types as needed (matches the surfrd_veg_all collapse call).
    # collapse_crop_types! expects a (gridcell, cft) layout, so transpose in/out.
    if collapse_crops
        wt_gc   = permutedims(wtcft_cur, (2, 1))             # (ngridcells, cft)
        fert_gc = permutedims(fertcft_cur, (2, 1))
        collapse_crop_types!(wt_gc, fert_gc, bounds.begg, bounds.endg)
        wtcft_cur   = permutedims(wt_gc, (2, 1))
        fertcft_cur = permutedims(fert_gc, (2, 1))
    end

    # whether we have set the weight for each column (single-CFT-per-column guard).
    col_set = falses(bounds.endc)

    for p in bounds.begp:bounds.endp
        g = pch.gridcell[p]
        l = pch.landunit[p]
        c = pch.column[p]

        if lun.itype[l] == ISTCROP
            # CFT index on the crop axis. Fortran patch itype is on the cft axis
            # (cft_lb..cft_ub); the 2-d data column index is 1-based here. We map
            # the patch itype to a 1-based column via cft_lb.
            m = pch.itype[p] - varpar.cft_lb + 1

            # The following assumes a single CFT on each crop column.
            if col_set[c]
                error("dyncrop_interp! ERROR: attempt to set a column that has " *
                      "already been set (c=$c). This may happen if there are " *
                      "multiple crops on a single column.")
            end

            col.wtlunit[c] = wtcft_cur[m, g]
            if use_crop && crop_inst !== nothing
                crop_inst.fertnitro_patch[p] = fertcft_cur[m, g]
            end
            col_set[c] = true
        end
    end

    return nothing
end

# Convenience: close the underlying files when done.
dynpft_close!(state::DynpftState)  = dyn_file_close!(state.dynpft_file)
dyncrop_close!(state::DyncropState) = dyn_file_close!(state.dyncrop_file)
