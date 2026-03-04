# ==========================================================================
# Ported from: src/main/surfrdMod.F90 + clm_varsur.F90
# Surface data reader — reads NetCDF surface dataset and produces
# landunit/PFT weights plus soil properties.
# ==========================================================================

"""
    SurfaceInputData

Container for all data read from the surface dataset (surfdata_clm.nc).
Replaces Fortran module-level variables in clm_instur/clm_varsur.
"""
Base.@kwdef mutable struct SurfaceInputData
    # Landunit weights [ng, MAX_LUNIT]
    wt_lunit::Matrix{Float64}      = Matrix{Float64}(undef, 0, 0)

    # Natural PFT weights [ng, natpft_size]
    wt_nat_patch::Matrix{Float64}  = Matrix{Float64}(undef, 0, 0)

    # Crop functional type weights [ng, cft_size]
    wt_cft::Matrix{Float64}        = Matrix{Float64}(undef, 0, 0)

    # Urban properties
    urban_valid::Vector{Bool}      = Bool[]
    pct_urban_max::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)

    # Fertilizer for CFTs
    fert_cft::Matrix{Float64}      = Matrix{Float64}(undef, 0, 0)

    # Soil properties [ng, nlevsoi]
    pct_sand::Matrix{Float64}      = Matrix{Float64}(undef, 0, 0)
    pct_clay::Matrix{Float64}      = Matrix{Float64}(undef, 0, 0)
    organic::Matrix{Float64}       = Matrix{Float64}(undef, 0, 0)
    soil_color::Vector{Int}        = Int[]
    zbedrock::Vector{Float64}      = Float64[]

    # Lake/topography
    lakedepth::Vector{Float64}     = Float64[]
    slope::Vector{Float64}         = Float64[]
    std_elev::Vector{Float64}      = Float64[]
    fmax::Vector{Float64}          = Float64[]    # max fractional saturated area (TOPMODEL)

    # Monthly phenology [ng, npft, 12]
    monthly_lai::Array{Float64,3}  = Array{Float64}(undef, 0, 0, 0)
    monthly_sai::Array{Float64,3}  = Array{Float64}(undef, 0, 0, 0)
    monthly_htop::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)
    monthly_hbot::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)

    # Glacier MEC weights [ng, maxpatch_glc]
    wt_glc_mec::Matrix{Float64}    = Matrix{Float64}(undef, 0, 0)
end

# --------------------------------------------------------------------------
# NetCDF helpers: surface data files have spatial dims (lsmlon, lsmlat)
# which may be 2D. These helpers read and flatten to 1D gridcell arrays.
# --------------------------------------------------------------------------

# Spatial dimension names used in CLM surface files
const _SPATIAL_DIMS = Set(["lsmlon", "lsmlat", "gridcell"])

# Convert a value that might be Missing to Float64
_to_f64(val) = ismissing(val) ? 0.0 : Float64(val)

"""
    _read_spatial_scalar(ds, varname, ng, default) -> Vector{Float64}

Read a spatially-indexed scalar variable (dims ⊆ {lsmlon, lsmlat}).
Returns a 1D vector of length ng by flattening the spatial dimensions.
"""
function _read_spatial_scalar(ds::NCDataset, varname::String, ng::Int, default::Float64)
    if !haskey(ds, varname)
        return fill(default, ng)
    end
    raw = Array(ds[varname])
    flat = vec(raw)
    result = fill(default, ng)
    for i in 1:min(ng, length(flat))
        result[i] = ismissing(flat[i]) ? default : Float64(flat[i])
    end
    return result
end

"""
    _read_spatial_layered(ds, varname, ng, nlev) -> Matrix{Float64} [ng, nlev]

Read a variable with a level/type dimension plus spatial dims.
Uses NCDataset dimension names to identify which axis is spatial vs data.
The data dimension (natpft, nlevsoi, etc.) may differ from the model's `nlev`;
we read min(file_nlev, nlev) levels and leave the rest as zero.
Returns [ng, nlev] matrix.
"""
function _read_spatial_layered(ds::NCDataset, varname::String, ng::Int, nlev::Int)
    result = zeros(ng, nlev)
    if !haskey(ds, varname)
        return result
    end

    ncvar = ds[varname]
    dimnames = NCDatasets.dimnames(ncvar)
    data = Array(ncvar)
    dims = size(data)
    nd = ndims(data)

    # Identify which dimension index is the data (non-spatial) dimension
    data_dim_idx = 0
    for i in 1:nd
        if !(dimnames[i] in _SPATIAL_DIMS)
            data_dim_idx = i
            break
        end
    end

    if data_dim_idx == 0
        # All dims are spatial — shouldn't happen for layered data
        flat = vec(data)
        for g in 1:min(ng, length(flat))
            result[g, 1] = _to_f64(flat[g])
        end
        return result
    end

    file_nlev = dims[data_dim_idx]
    nlev_read = min(file_nlev, nlev)

    if nd == 1
        # Only data dim, no spatial (unusual)
        for j in 1:nlev_read
            result[1, j] = _to_f64(data[j])
        end
    elseif nd == 2
        if data_dim_idx == 1
            # (nlev_file, spatial)
            for g in 1:min(ng, dims[2]), j in 1:nlev_read
                result[g, j] = _to_f64(data[j, g])
            end
        else
            # (spatial, nlev_file)
            for g in 1:min(ng, dims[1]), j in 1:nlev_read
                result[g, j] = _to_f64(data[g, j])
            end
        end
    elseif nd == 3
        # Two spatial dims + one data dim
        if data_dim_idx == 1
            # (nlev_file, lsmlon, lsmlat)
            nspat = dims[2] * dims[3]
            reshaped = reshape(data, file_nlev, nspat)
            for g in 1:min(ng, nspat), j in 1:nlev_read
                result[g, j] = _to_f64(reshaped[j, g])
            end
        elseif data_dim_idx == 2
            # (lsmlon, nlev_file, lsmlat)
            for g_spat in 1:min(ng, dims[1] * dims[3]), j in 1:nlev_read
                gi = ((g_spat - 1) % dims[1]) + 1
                gj = ((g_spat - 1) ÷ dims[1]) + 1
                result[g_spat, j] = _to_f64(data[gi, j, gj])
            end
        else
            # (lsmlon, lsmlat, nlev_file)
            nspat = dims[1] * dims[2]
            reshaped = reshape(data, nspat, file_nlev)
            for g in 1:min(ng, nspat), j in 1:nlev_read
                result[g, j] = _to_f64(reshaped[g, j])
            end
        end
    end
    return result
end

"""
    surfrd_get_num_patches(fsurdat) -> (numpft, numcft)

Read surface dataset metadata to determine PFT/CFT dimensions.
"""
function surfrd_get_num_patches(fsurdat::String)
    ds = NCDataset(fsurdat, "r")
    try
        numpft = haskey(ds.dim, "natpft") ? ds.dim["natpft"] : 15
        numcft = haskey(ds.dim, "cft") ? ds.dim["cft"] : 2
        return (numpft, numcft)
    finally
        close(ds)
    end
end

"""
    surfrd_get_nlevurb(fsurdat) -> Int

Read urban soil levels from surface dataset.
"""
function surfrd_get_nlevurb(fsurdat::String)
    ds = NCDataset(fsurdat, "r")
    try
        return haskey(ds.dim, "nlevurb") ? ds.dim["nlevurb"] : 5
    finally
        close(ds)
    end
end

"""
    surfrd_get_data!(surf, begg, endg, fsurdat)

Main reader: populate SurfaceInputData from surface NetCDF file.
"""
function surfrd_get_data!(surf::SurfaceInputData, begg::Int, endg::Int, fsurdat::String)
    ng = endg - begg + 1

    ds = NCDataset(fsurdat, "r")
    try
        # Allocate landunit weights
        surf.wt_lunit = zeros(ng, MAX_LUNIT)

        # Read special landunits
        surfrd_special!(surf, ds, ng)

        # Read vegetation/crop weights
        surfrd_veg_all!(surf, ds, ng)

        # Read soil properties
        surfrd_soil!(surf, ds, ng)

        # Read topography
        surfrd_topo!(surf, ds, ng)

        # Read monthly phenology
        surfrd_phenology!(surf, ds, ng)

    finally
        close(ds)
    end

    # Apply ocean-to-land conversion if requested
    if varctl.convert_ocean_to_land
        apply_convert_ocean_to_land!(surf.wt_lunit)
    end

    # Apply landunit collapses
    if any(x -> x > 0.0, [varctl.toosmall_soil, varctl.toosmall_crop,
                            varctl.toosmall_glacier, varctl.toosmall_lake,
                            varctl.toosmall_wetland, varctl.toosmall_urban])
        collapse_individual_lunits!(surf.wt_lunit;
            toosmall_soil=varctl.toosmall_soil,
            toosmall_crop=varctl.toosmall_crop,
            toosmall_glacier=varctl.toosmall_glacier,
            toosmall_lake=varctl.toosmall_lake,
            toosmall_wetland=varctl.toosmall_wetland,
            toosmall_urban=varctl.toosmall_urban)
    end

    # Collapse to dominant landunits if requested
    if varctl.n_dom_landunits > 0
        collapse_to_dominant!(surf.wt_lunit, varctl.n_dom_landunits)
    end

    # Normalize landunit weights
    renormalize!(surf.wt_lunit; target=1.0)

    nothing
end

"""
    surfrd_special!(surf, ds, ng)

Read special landunit fractions: PCT_LAKE, PCT_WETLAND, PCT_GLACIER, PCT_URBAN.
"""
function surfrd_special!(surf::SurfaceInputData, ds::NCDataset, ng::Int)
    pct_lake = _read_spatial_scalar(ds, "PCT_LAKE", ng, 0.0)
    pct_wetland = _read_spatial_scalar(ds, "PCT_WETLAND", ng, 0.0)
    pct_glacier = _read_spatial_scalar(ds, "PCT_GLACIER", ng, 0.0)

    # Urban: may have density classes (lsmlon, lsmlat, numurbl) or (numurbl, lsmlon, lsmlat)
    pct_urban = zeros(ng, NUMURBL)
    if haskey(ds, "PCT_URBAN")
        data = Array(ds["PCT_URBAN"])
        dims = size(data)
        nd = ndims(data)
        if nd == 1
            for g in 1:min(ng, length(data))
                pct_urban[g, 1] = _to_f64(data[g])
            end
        elseif nd == 2
            # (lsmlon*lsmlat, numurbl) or (numurbl, lsmlon*lsmlat)
            if dims[2] == NUMURBL || (dims[1] != NUMURBL && dims[2] <= NUMURBL)
                for n in 1:min(NUMURBL, dims[2]), g in 1:min(ng, dims[1])
                    pct_urban[g, n] = _to_f64(data[g, n])
                end
            else
                for n in 1:min(NUMURBL, dims[1]), g in 1:min(ng, dims[2])
                    pct_urban[g, n] = _to_f64(data[n, g])
                end
            end
        elseif nd == 3
            # (numurbl, lsmlon, lsmlat) or (lsmlon, lsmlat, numurbl)
            if dims[1] == NUMURBL
                nspat = dims[2] * dims[3]
                reshaped = reshape(data, NUMURBL, nspat)
                for n in 1:NUMURBL, g in 1:min(ng, nspat)
                    pct_urban[g, n] = _to_f64(reshaped[n, g])
                end
            elseif dims[3] == NUMURBL
                nspat = dims[1] * dims[2]
                reshaped = reshape(data, nspat, NUMURBL)
                for n in 1:NUMURBL, g in 1:min(ng, nspat)
                    pct_urban[g, n] = _to_f64(reshaped[g, n])
                end
            end
        end
    end
    surf.pct_urban_max = copy(pct_urban)

    # Convert percentages to fractions
    for g in 1:ng
        surf.wt_lunit[g, ISTDLAK] = pct_lake[g] / 100.0
        surf.wt_lunit[g, ISTWET] = pct_wetland[g] / 100.0
        surf.wt_lunit[g, ISTICE] = pct_glacier[g] / 100.0
        for n in 1:NUMURBL
            surf.wt_lunit[g, ISTURB_MIN + n - 1] = pct_urban[g, n] / 100.0
        end
    end

    # Read PCT_OCEAN if available, store in ISTOCN
    pct_ocean = _read_spatial_scalar(ds, "PCT_OCEAN", ng, 0.0)
    if haskey(ds, "PCT_OCEAN")
        for g in 1:ng
            surf.wt_lunit[g, ISTOCN] = pct_ocean[g] / 100.0
        end
    end

    nothing
end

"""
    surfrd_veg_all!(surf, ds, ng)

Read natural vegetation and crop PFT weights from surface dataset.
"""
function surfrd_veg_all!(surf::SurfaceInputData, ds::NCDataset, ng::Int)
    # Read PCT_NATVEG (landunit fraction for soil)
    # Note: CLM surface files use "PCT_NATVEG" (not "PCT_NAT_VEG")
    pct_nat_veg = _read_spatial_scalar(ds, "PCT_NATVEG", ng, 0.0)

    # Read PCT_CROP (landunit fraction for crop)
    pct_crop = _read_spatial_scalar(ds, "PCT_CROP", ng, 0.0)

    # Set soil and crop landunit weights
    for g in 1:ng
        surf.wt_lunit[g, ISTSOIL] = pct_nat_veg[g] / 100.0
        if varctl.create_crop_landunit
            surf.wt_lunit[g, ISTCROP] = pct_crop[g] / 100.0
        else
            surf.wt_lunit[g, ISTSOIL] += pct_crop[g] / 100.0
        end
    end

    # Read PCT_NAT_PFT: dims are (natpft, lsmlon, lsmlat) in the file
    natpft_size = varpar.natpft_ub - varpar.natpft_lb + 1
    surf.wt_nat_patch = _read_spatial_layered(ds, "PCT_NAT_PFT", ng, natpft_size)

    # Convert from percentages to fractions and normalize
    surf.wt_nat_patch ./= 100.0
    for g in 1:ng
        s = sum(@view surf.wt_nat_patch[g, :])
        if s > 0.0
            surf.wt_nat_patch[g, :] ./= s
        else
            # Default: all bare ground (PFT 0 = index 1)
            surf.wt_nat_patch[g, 1] = 1.0
        end
    end

    # Read PCT_CFT [ng, cft_size]
    cft_size = varpar.cft_size
    if cft_size > 0
        surf.wt_cft = _read_spatial_layered(ds, "PCT_CFT", ng, cft_size)
        surf.wt_cft ./= 100.0
        # Normalize rows
        for g in 1:ng
            s = sum(@view surf.wt_cft[g, :])
            if s > 0.0
                surf.wt_cft[g, :] ./= s
            elseif cft_size >= 1
                surf.wt_cft[g, 1] = 1.0
            end
        end
        # Read fertilizer
        surf.fert_cft = _read_spatial_layered(ds, "FERTNITRO_CFT", ng, cft_size)
    else
        surf.wt_cft = zeros(ng, 0)
        surf.fert_cft = zeros(ng, 0)
    end

    # Glacier MEC weights (simplified: single column)
    maxpatch_glc = varpar.maxpatch_glc
    surf.wt_glc_mec = zeros(ng, maxpatch_glc)
    if maxpatch_glc > 0
        # Default: all glacier weight in class 1
        for g in 1:ng
            if surf.wt_lunit[g, ISTICE] > 0.0
                surf.wt_glc_mec[g, 1] = 1.0
            end
        end
    end

    nothing
end

"""
    surfrd_soil!(surf, ds, ng)

Read soil properties from surface dataset.
"""
function surfrd_soil!(surf::SurfaceInputData, ds::NCDataset, ng::Int)
    nlevsoi = varpar.nlevsoi

    # PCT_SAND, PCT_CLAY, ORGANIC: dims (nlevsoi, lsmlon, lsmlat) in file
    surf.pct_sand = _read_spatial_layered(ds, "PCT_SAND", ng, nlevsoi)
    surf.pct_clay = _read_spatial_layered(ds, "PCT_CLAY", ng, nlevsoi)
    surf.organic = _read_spatial_layered(ds, "ORGANIC", ng, nlevsoi)

    # SOIL_COLOR — scalar spatial
    surf.soil_color = fill(1, ng)
    sc = _read_spatial_scalar(ds, "SOIL_COLOR", ng, 1.0)
    for g in 1:ng
        surf.soil_color[g] = max(1, Int(round(sc[g])))
    end

    # zbedrock — scalar spatial
    surf.zbedrock = fill(zisoi[][varpar.nlevsoi], ng)
    zb = _read_spatial_scalar(ds, "zbedrock", ng, NaN)
    for g in 1:ng
        if !isnan(zb[g])
            surf.zbedrock[g] = zb[g]
        end
    end

    nothing
end

"""
    surfrd_topo!(surf, ds, ng)

Read topographic data from surface dataset.
"""
function surfrd_topo!(surf::SurfaceInputData, ds::NCDataset, ng::Int)
    surf.lakedepth = _read_spatial_scalar(ds, "LAKEDEPTH", ng, SPVAL)
    surf.slope = _read_spatial_scalar(ds, "SLOPE", ng, 0.2)
    # Enforce minimum slope
    for g in 1:ng
        surf.slope[g] = max(surf.slope[g], 0.2)
    end
    surf.std_elev = _read_spatial_scalar(ds, "STD_ELEV", ng, 0.0)
    surf.fmax = _read_spatial_scalar(ds, "FMAX", ng, 0.0)
    nothing
end

"""
    surfrd_phenology!(surf, ds, ng)

Read monthly LAI/SAI/HTOP/HBOT from surface dataset.
The surface file stores these as (lsmlon, lsmlat, lsmpft, time).
`lsmpft` includes ALL PFTs (bare+natural+crop), typically 17.
We read the first `npft` natural PFTs.
"""
function surfrd_phenology!(surf::SurfaceInputData, ds::NCDataset, ng::Int)
    natpft_size = varpar.natpft_ub - varpar.natpft_lb + 1
    npft = natpft_size

    surf.monthly_lai = zeros(ng, npft, 12)
    surf.monthly_sai = zeros(ng, npft, 12)
    surf.monthly_htop = zeros(ng, npft, 12)
    surf.monthly_hbot = zeros(ng, npft, 12)

    for (varname, target) in [("MONTHLY_LAI", surf.monthly_lai),
                               ("MONTHLY_SAI", surf.monthly_sai),
                               ("MONTHLY_HEIGHT_TOP", surf.monthly_htop),
                               ("MONTHLY_HEIGHT_BOT", surf.monthly_hbot)]
        if !haskey(ds, varname)
            continue
        end

        ncvar = ds[varname]
        dimnames = NCDatasets.dimnames(ncvar)
        data = Array(ncvar)
        dims = size(data)
        nd = ndims(data)

        # Identify dimension indices by name
        pft_idx = 0
        time_idx = 0
        spatial_idxs = Int[]
        for i in 1:nd
            dn = dimnames[i]
            if dn == "time" || dn == "nmonth"
                time_idx = i
            elseif dn in ("lsmpft", "natpft", "pft")
                pft_idx = i
            elseif dn in _SPATIAL_DIMS
                push!(spatial_idxs, i)
            else
                # Unknown dim: try to infer from size
                if dims[i] == 12 && time_idx == 0
                    time_idx = i
                elseif pft_idx == 0
                    pft_idx = i
                end
            end
        end

        # Fallback: if we couldn't identify time dim by name, use size
        if time_idx == 0
            for i in 1:nd
                if dims[i] == 12 && i != pft_idx && !(i in spatial_idxs)
                    time_idx = i
                    break
                end
            end
        end

        if pft_idx == 0 || time_idx == 0
            continue  # can't identify dimensions
        end

        npft_file = dims[pft_idx]
        npft_read = min(npft, npft_file)
        ntime = min(12, dims[time_idx])

        # Compute spatial size
        nspat = 1
        for si in spatial_idxs
            nspat *= dims[si]
        end
        ng_read = min(ng, nspat)

        # General N-D copy using computed indices
        if nd == 3
            for g in 1:ng_read, m in 1:npft_read, mon in 1:ntime
                idx = ntuple(nd) do i
                    i == pft_idx ? m : i == time_idx ? mon : g
                end
                target[g, m, mon] = _to_f64(data[idx...])
            end
        elseif nd == 4 && length(spatial_idxs) == 2
            sp1 = spatial_idxs[1]
            sp2 = spatial_idxs[2]
            nsp1 = dims[sp1]
            for g in 1:ng_read, m in 1:npft_read, mon in 1:ntime
                g1 = ((g - 1) % nsp1) + 1
                g2 = ((g - 1) ÷ nsp1) + 1
                idx = ntuple(nd) do i
                    i == pft_idx ? m : i == time_idx ? mon : i == sp1 ? g1 : g2
                end
                target[g, m, mon] = _to_f64(data[idx...])
            end
        end
    end

    nothing
end

"""
    count_subgrid_elements(surf, ng) -> (nl, nc, np)

Count landunits, columns, and patches from surface data weights.
Determines the total number of each subgrid element needed.
"""
function count_subgrid_elements(surf::SurfaceInputData, ng::Int)
    nl = 0  # landunits
    nc = 0  # columns
    np = 0  # patches

    natpft_size = size(surf.wt_nat_patch, 2)
    cft_size = size(surf.wt_cft, 2)

    for g in 1:ng
        # Natural vegetation landunit (ISTSOIL)
        if surf.wt_lunit[g, ISTSOIL] > 0.0 || true  # always create soil landunit
            nl += 1
            nc += 1  # one column for natural veg
            # Count natural PFTs with non-zero weight
            n_nat = 0
            for m in 1:natpft_size
                if surf.wt_nat_patch[g, m] > 0.0
                    n_nat += 1
                end
            end
            np += max(n_nat, 1)  # at least bare ground
        end

        # Crop landunit (ISTCROP)
        if varctl.create_crop_landunit && surf.wt_lunit[g, ISTCROP] > 0.0
            nl += 1
            # Each crop type gets its own column+patch
            n_cft = 0
            for m in 1:cft_size
                if surf.wt_cft[g, m] > 0.0
                    n_cft += 1
                end
            end
            if n_cft > 0
                nc += n_cft
                np += n_cft
            else
                nc += 1
                np += 1
            end
        end

        # Urban landunits (TBD, HD, MD)
        for u in [ISTURB_TBD, ISTURB_HD, ISTURB_MD]
            if surf.wt_lunit[g, u] > 0.0 || varctl.run_zero_weight_urban
                nl += 1
                nc += 5   # roof, sunwall, shadewall, imperv road, perv road
                np += 5   # one patch per urban column
            end
        end

        # Lake (ISTDLAK)
        if surf.wt_lunit[g, ISTDLAK] > 0.0 || true  # always create lake
            nl += 1
            nc += 1
            np += 1
        end

        # Wetland (ISTWET)
        if surf.wt_lunit[g, ISTWET] > 0.0
            nl += 1
            nc += 1
            np += 1
        end

        # Glacier (ISTICE)
        if surf.wt_lunit[g, ISTICE] > 0.0
            nl += 1
            # One column per elevation class that has weight
            n_glc = 0
            if size(surf.wt_glc_mec, 2) > 0
                for m in 1:size(surf.wt_glc_mec, 2)
                    if surf.wt_glc_mec[g, m] > 0.0
                        n_glc += 1
                    end
                end
            end
            n_glc = max(n_glc, 1)
            nc += n_glc
            np += n_glc
        end
    end

    return (nl, nc, np)
end
