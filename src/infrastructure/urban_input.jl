# ==========================================================================
# Ported from: src/biogeophys/UrbanParamsType.F90 (UrbanInput / readUrbanInput)
# Reads the urban morphology / thermal / radiative fields from the surface dataset
# into an UrbanInputData. Without this the urban (isturb) land-unit path runs on
# fill(NaN) morphology -> degenerate. Wired into clm_initialize! after the landunit
# subgrid is built and urban_read_nml! has run.
#
# Surfdata dim layout (CTSM convention):
#   2D fields  (CANYON_HWR, WTLUNIT_ROOF, ...) : (numurbl, lat, lon)
#   albedo     (ALB_*_DIR/DIF)                 : (numrad,  numurbl, lat, lon)
#   thermal    (TK_*, CV_*)                    : (nlevurb, numurbl, lat, lon)
# Julia UrbanInputData wants (ng, numurbl) and (ng, numurbl, numrad|nlevurb).
# ==========================================================================

# Read a (data, numurbl, spatial...) urban field into out[ng, numurbl, ndata].
# Robust to dim order: locates the numurbl axis, the remaining non-spatial (data)
# axis, and the spatial axes by dimension NAME, then permutes to (spatial.., urbl, data).
function _read_urban_layered3!(out::AbstractArray{<:Real,3}, ds::NCDataset, varname::String,
                               ng::Int, numurbl::Int, ndata::Int)
    haskey(ds, varname) || return out
    v = ds[varname]
    dn = collect(NCDatasets.dimnames(v))
    data = Array(v)
    nd = ndims(data)
    urb_ax = findfirst(==("numurbl"), dn)
    urb_ax === nothing && return out
    dat_ax = findfirst(i -> !(dn[i] in _SPATIAL_DIMS) && i != urb_ax, 1:nd)
    dat_ax === nothing && return out
    spat = [i for i in 1:nd if dn[i] in _SPATIAL_DIMS]
    pd = permutedims(data, vcat(spat, [urb_ax, dat_ax]))         # (spatial.., numurbl, ndata)
    flat = reshape(pd, :, size(pd, ndims(pd) - 1), size(pd, ndims(pd)))  # (ng_file, numurbl, ndata)
    @inbounds for g in 1:min(ng, size(flat, 1)), u in 1:min(numurbl, size(flat, 2)),
                  k in 1:min(ndata, size(flat, 3))
        out[g, u, k] = _to_f64(flat[g, u, k])
    end
    return out
end

"""
    read_urban_input!(urbinp, fsurdat, ng, numurbl, numrad, nlevurb)

Populate `urbinp` (UrbanInputData) from the surface dataset `fsurdat`. Fields absent
from the file are left at their `urbinp_init!` defaults (NaN / 0). T_BUILDING_MAX is
not in UrbanInputData (only T_BUILDING_MIN is read).
"""
function read_urban_input!(urbinp::UrbanInputData, fsurdat::String,
                           ng::Int, numurbl::Int, numrad::Int, nlevurb::Int)
    ds = NCDataset(fsurdat, "r")
    try
        # 2D (gridcell, numurbl) — _read_spatial_layered treats numurbl as the layer dim.
        urbinp.canyon_hwr      .= _read_spatial_layered(ds, "CANYON_HWR",      ng, numurbl)
        urbinp.wtlunit_roof    .= _read_spatial_layered(ds, "WTLUNIT_ROOF",    ng, numurbl)
        urbinp.wtroad_perv     .= _read_spatial_layered(ds, "WTROAD_PERV",     ng, numurbl)
        urbinp.em_roof         .= _read_spatial_layered(ds, "EM_ROOF",         ng, numurbl)
        urbinp.em_improad      .= _read_spatial_layered(ds, "EM_IMPROAD",      ng, numurbl)
        urbinp.em_perroad      .= _read_spatial_layered(ds, "EM_PERROAD",      ng, numurbl)
        urbinp.em_wall         .= _read_spatial_layered(ds, "EM_WALL",         ng, numurbl)
        urbinp.ht_roof         .= _read_spatial_layered(ds, "HT_ROOF",         ng, numurbl)
        urbinp.wind_hgt_canyon .= _read_spatial_layered(ds, "WIND_HGT_CANYON", ng, numurbl)
        urbinp.thick_wall      .= _read_spatial_layered(ds, "THICK_WALL",      ng, numurbl)
        urbinp.thick_roof      .= _read_spatial_layered(ds, "THICK_ROOF",      ng, numurbl)
        urbinp.t_building_min  .= _read_spatial_layered(ds, "T_BUILDING_MIN",  ng, numurbl)
        urbinp.nlev_improad    .= round.(Int, _read_spatial_layered(ds, "NLEV_IMPROAD", ng, numurbl))

        # 3D albedo (gridcell, numurbl, numrad)
        for (nm, fld) in (("ALB_ROOF_DIR", urbinp.alb_roof_dir), ("ALB_ROOF_DIF", urbinp.alb_roof_dif),
                          ("ALB_IMPROAD_DIR", urbinp.alb_improad_dir), ("ALB_IMPROAD_DIF", urbinp.alb_improad_dif),
                          ("ALB_PERROAD_DIR", urbinp.alb_perroad_dir), ("ALB_PERROAD_DIF", urbinp.alb_perroad_dif),
                          ("ALB_WALL_DIR", urbinp.alb_wall_dir), ("ALB_WALL_DIF", urbinp.alb_wall_dif))
            _read_urban_layered3!(fld, ds, nm, ng, numurbl, numrad)
        end
        # 3D thermal (gridcell, numurbl, nlevurb)
        for (nm, fld) in (("TK_ROOF", urbinp.tk_roof), ("TK_WALL", urbinp.tk_wall), ("TK_IMPROAD", urbinp.tk_improad),
                          ("CV_ROOF", urbinp.cv_roof), ("CV_WALL", urbinp.cv_wall), ("CV_IMPROAD", urbinp.cv_improad))
            _read_urban_layered3!(fld, ds, nm, ng, numurbl, nlevurb)
        end
    finally
        close(ds)
    end
    return urbinp
end
