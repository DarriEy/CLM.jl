# ==========================================================================
# Ported from: src/main/initVerticalMod.F90
# Initialize vertical components of column datatype
# ==========================================================================

"""
    initVertical!(bounds, grc, lun, col, surf;
                  thick_wall=0.5, thick_roof=0.5)

Initialize vertical layer structure for all columns. Sets col.dz, col.z,
col.zi for soil/urban/lake columns, sets bedrock index, lake layers,
topography, and layer classification.
"""
function initVertical!(bounds::BoundsType, grc::GridcellData,
                        lun::LandunitData, col::ColumnData,
                        surf::SurfaceInputData;
                        thick_wall::Float64=0.5, thick_roof::Float64=0.5)
    begc, endc = bounds.begc, bounds.endc
    begl, endl = bounds.begl, bounds.endl
    begg, endg = bounds.begg, bounds.endg

    nlevsno = varpar.nlevsno
    nlevgrnd = varpar.nlevgrnd
    nlevlak = varpar.nlevlak
    nlevurb = varpar.nlevurb
    nlevmaxurbgrnd = varpar.nlevmaxurbgrnd
    joff = nlevsno  # offset: Fortran j → Julia j+joff

    # Get global soil coordinate arrays
    zsoi_g  = zsoi[]
    dzsoi_g = dzsoi[]
    zisoi_g = zisoi[]
    zlak_g  = zlak[]
    dzlak_g = dzlak[]

    # --- Urban wall/roof layer geometry (per-landunit) ---
    if nlevurb > 0
        zurb_wall  = zeros(endl, nlevurb)
        zurb_roof  = zeros(endl, nlevurb)
        dzurb_wall = zeros(endl, nlevurb)
        dzurb_roof = zeros(endl, nlevurb)
        ziurb_wall = zeros(endl, nlevurb + 1)  # 0:nlevurb → Julia 1:nlevurb+1
        ziurb_roof = zeros(endl, nlevurb + 1)

        for l in begl:endl
            if lun.urbpoi[l]
                # Node depths
                for j in 1:nlevurb
                    zurb_wall[l, j] = (j - 0.5) * (thick_wall / nlevurb)
                    zurb_roof[l, j] = (j - 0.5) * (thick_roof / nlevurb)
                end

                # Thicknesses
                if nlevurb >= 2
                    dzurb_roof[l, 1] = 0.5 * (zurb_roof[l, 1] + zurb_roof[l, 2])
                    dzurb_wall[l, 1] = 0.5 * (zurb_wall[l, 1] + zurb_wall[l, 2])
                    for j in 2:nlevurb-1
                        dzurb_roof[l, j] = 0.5 * (zurb_roof[l, j+1] - zurb_roof[l, j-1])
                        dzurb_wall[l, j] = 0.5 * (zurb_wall[l, j+1] - zurb_wall[l, j-1])
                    end
                    dzurb_roof[l, nlevurb] = zurb_roof[l, nlevurb] - zurb_roof[l, nlevurb-1]
                    dzurb_wall[l, nlevurb] = zurb_wall[l, nlevurb] - zurb_wall[l, nlevurb-1]
                else
                    dzurb_roof[l, 1] = thick_roof
                    dzurb_wall[l, 1] = thick_wall
                end

                # Interfaces (index 1 = surface = 0)
                ziurb_wall[l, 1] = 0.0
                ziurb_roof[l, 1] = 0.0
                for j in 1:nlevurb-1
                    ziurb_wall[l, j+1] = 0.5 * (zurb_wall[l, j] + zurb_wall[l, j+1])
                    ziurb_roof[l, j+1] = 0.5 * (zurb_roof[l, j] + zurb_roof[l, j+1])
                end
                ziurb_wall[l, nlevurb+1] = zurb_wall[l, nlevurb] + 0.5 * dzurb_wall[l, nlevurb]
                ziurb_roof[l, nlevurb+1] = zurb_roof[l, nlevurb] + 0.5 * dzurb_roof[l, nlevurb]
            end
        end
    end

    # --- Set soil/urban vertical structure per column ---
    for c in begc:endc
        l = col.landunit[c]

        if lun.urbpoi[l]
            if col.itype[c] == ICOL_SUNWALL || col.itype[c] == ICOL_SHADEWALL
                for j in 1:nlevurb
                    col.z[c, j + joff]  = zurb_wall[l, j]
                    col.dz[c, j + joff] = dzurb_wall[l, j]
                end
                # zi: interface 0 (surface) at index joff+1, then 1:nlevurb
                col.zi[c, joff + 1] = ziurb_wall[l, 1]  # = 0
                for j in 1:nlevurb
                    col.zi[c, j + joff + 1] = ziurb_wall[l, j + 1]
                end
                # Mark unused layers as SPVAL
                for j in (nlevurb + 1):nlevmaxurbgrnd
                    col.z[c, j + joff]  = SPVAL
                    col.dz[c, j + joff] = SPVAL
                    col.zi[c, j + joff + 1] = SPVAL
                end

            elseif col.itype[c] == ICOL_ROOF
                for j in 1:nlevurb
                    col.z[c, j + joff]  = zurb_roof[l, j]
                    col.dz[c, j + joff] = dzurb_roof[l, j]
                end
                col.zi[c, joff + 1] = ziurb_roof[l, 1]
                for j in 1:nlevurb
                    col.zi[c, j + joff + 1] = ziurb_roof[l, j + 1]
                end
                for j in (nlevurb + 1):nlevmaxurbgrnd
                    col.z[c, j + joff]  = SPVAL
                    col.dz[c, j + joff] = SPVAL
                    col.zi[c, j + joff + 1] = SPVAL
                end

            else
                # Urban road: use standard soil layers
                _set_standard_soil!(col, c, joff, nlevgrnd, nlevmaxurbgrnd, zsoi_g, dzsoi_g, zisoi_g)
            end
        elseif lun.itype[l] != ISTDLAK
            # Non-lake, non-urban: standard soil
            _set_standard_soil!(col, c, joff, nlevgrnd, nlevmaxurbgrnd, zsoi_g, dzsoi_g, zisoi_g)
        end
    end

    # --- Set bedrock index ---
    init_bedrock!(bounds, grc, col, surf)

    # --- Set lake layers ---
    init_lake_layers!(bounds, col, lun, surf, zlak_g, dzlak_g,
                       zsoi_g, dzsoi_g, zisoi_g, joff, nlevgrnd, nlevlak, nlevmaxurbgrnd)

    # --- Set layer classification ---
    setSoilLayerClass!(bounds, col, lun)

    # --- Set topography ---
    init_topography!(bounds, col, surf)

    nothing
end

"""
    _set_standard_soil!(col, c, joff, nlevgrnd, nlevmaxurbgrnd, zsoi_g, dzsoi_g, zisoi_g)

Copy global soil coordinates into column c.
"""
function _set_standard_soil!(col::ColumnData, c::Int, joff::Int,
                              nlevgrnd::Int, nlevmaxurbgrnd::Int,
                              zsoi_g::Vector{Float64}, dzsoi_g::Vector{Float64},
                              zisoi_g::Vector{Float64})
    for j in 1:nlevgrnd
        col.z[c, j + joff]  = zsoi_g[j]
        col.dz[c, j + joff] = dzsoi_g[j]
    end
    # Interface: zisoi_g[1]=0 (surface), zisoi_g[j+1]=interface below layer j
    col.zi[c, joff + 1] = zisoi_g[1]  # = 0
    for j in 1:nlevgrnd
        col.zi[c, j + joff + 1] = zisoi_g[j + 1]
    end
    # Mark any extra levels as SPVAL
    for j in (nlevgrnd + 1):nlevmaxurbgrnd
        col.z[c, j + joff]  = SPVAL
        col.dz[c, j + joff] = SPVAL
        col.zi[c, j + joff + 1] = SPVAL
    end
    nothing
end

"""
    init_bedrock!(bounds, grc, col, surf)

Set gridcell and column bedrock index from surface data zbedrock.
"""
function init_bedrock!(bounds::BoundsType, grc::GridcellData,
                        col::ColumnData, surf::SurfaceInputData)
    nlevsoi = varpar.nlevsoi
    zisoi_g = zisoi[]
    zbedrock_in = copy(surf.zbedrock)

    # Apply parameter overrides
    if _initvert_zbedrock[] >= 0.0
        fill!(zbedrock_in, _initvert_zbedrock[])
    end
    if _initvert_zbedrock_sf[] != 1.0
        zbedrock_in .*= _initvert_zbedrock_sf[]
    end

    # If not using bedrock, set to bottom of soil column
    if !varctl.use_bedrock
        fill!(zbedrock_in, zisoi_g[nlevsoi + 1])  # zisoi[nlevsoi] in Fortran
    end

    # Determine minimum bedrock layer
    jmin_bedrock = 3
    for j in 3:nlevsoi
        if zisoi_g[j] < ZMIN_BEDROCK && zisoi_g[j + 1] >= ZMIN_BEDROCK
            jmin_bedrock = j
        end
    end

    # Set gridcell bedrock
    for g in bounds.begg:bounds.endg
        grc.nbedrock[g] = nlevsoi
        gi = g - bounds.begg + 1  # index into zbedrock_in
        for j in jmin_bedrock:nlevsoi
            if zisoi_g[j] < zbedrock_in[gi] && zisoi_g[j + 1] >= zbedrock_in[gi]
                grc.nbedrock[g] = j
            end
        end
    end

    # Copy to columns
    for c in bounds.begc:bounds.endc
        g = col.gridcell[c]
        col.nbedrock[c] = grc.nbedrock[g]
    end

    nothing
end

"""
    init_lake_layers!(bounds, col, lun, surf, zlak_g, dzlak_g,
                       zsoi_g, dzsoi_g, zisoi_g, joff, nlevgrnd, nlevlak, nlevmaxurbgrnd)

Set lake depth, z_lake, dz_lake for lake columns. Also set soil layers below lake.
"""
function init_lake_layers!(bounds::BoundsType, col::ColumnData,
                            lun::LandunitData, surf::SurfaceInputData,
                            zlak_g::Vector{Float64}, dzlak_g::Vector{Float64},
                            zsoi_g::Vector{Float64}, dzsoi_g::Vector{Float64},
                            zisoi_g::Vector{Float64},
                            joff::Int, nlevgrnd::Int, nlevlak::Int,
                            nlevmaxurbgrnd::Int)
    # Copy lakedepth from surface data
    for c in bounds.begc:bounds.endc
        g = col.gridcell[c]
        gi = g - bounds.begg + 1
        col.lakedepth[c] = surf.lakedepth[gi]
    end

    # Set lake layers and soil layers below lake
    for c in bounds.begc:bounds.endc
        l = col.landunit[c]
        lun.itype[l] == ISTDLAK || continue

        # Default lake depth if missing
        if col.lakedepth[c] == SPVAL || col.lakedepth[c] <= 0.0
            col.lakedepth[c] = zlak_g[nlevlak] + 0.5 * dzlak_g[nlevlak]
            for k in 1:nlevlak
                col.z_lake[c, k] = zlak_g[k]
                col.dz_lake[c, k] = dzlak_g[k]
            end
        elseif col.lakedepth[c] > 1.0 && col.lakedepth[c] < 5000.0
            depthratio = col.lakedepth[c] / (zlak_g[nlevlak] + 0.5 * dzlak_g[nlevlak])
            col.z_lake[c, 1]  = zlak_g[1]
            col.dz_lake[c, 1] = dzlak_g[1]
            for k in 2:nlevlak-1
                col.dz_lake[c, k] = dzlak_g[k] * depthratio
            end
            col.dz_lake[c, nlevlak] = dzlak_g[nlevlak] * depthratio -
                                       (col.dz_lake[c, 1] - dzlak_g[1] * depthratio)
            for k in 2:nlevlak
                col.z_lake[c, k] = col.z_lake[c, k-1] +
                    0.5 * (col.dz_lake[c, k-1] + col.dz_lake[c, k])
            end
        elseif col.lakedepth[c] > 0.0 && col.lakedepth[c] <= 1.0
            dz_uniform = col.lakedepth[c] / nlevlak
            for k in 1:nlevlak
                col.dz_lake[c, k] = dz_uniform
            end
            col.z_lake[c, 1] = dz_uniform / 2.0
            for k in 2:nlevlak
                col.z_lake[c, k] = col.z_lake[c, k-1] +
                    0.5 * (col.dz_lake[c, k-1] + col.dz_lake[c, k])
            end
        end

        # Lake columns also have soil layers below the lake bed
        _set_standard_soil!(col, c, joff, nlevgrnd, nlevmaxurbgrnd,
                            zsoi_g, dzsoi_g, zisoi_g)
    end

    nothing
end

"""
    setSoilLayerClass!(bounds, col, lun)

Set levgrnd_class for each column layer based on bedrock depth.
"""
function setSoilLayerClass!(bounds::BoundsType, col::ColumnData,
                             lun::LandunitData)
    nlevsoi = varpar.nlevsoi
    nlevmaxurbgrnd = varpar.nlevmaxurbgrnd
    nlevsno = varpar.nlevsno
    joff = nlevsno

    LEVGRND_CLASS_STANDARD        = 1
    LEVGRND_CLASS_DEEP_BEDROCK    = 2
    LEVGRND_CLASS_SHALLOW_BEDROCK = 3

    for c in bounds.begc:bounds.endc
        l = col.landunit[c]
        if hasBedrock(col.itype[c], lun.itype[l])
            nb = col.nbedrock[c]
            for j in 1:min(nb, nlevmaxurbgrnd)
                col.levgrnd_class[c, j] = LEVGRND_CLASS_STANDARD
            end
            if nb < nlevsoi
                for j in (nb + 1):nlevsoi
                    col.levgrnd_class[c, j] = LEVGRND_CLASS_SHALLOW_BEDROCK
                end
            end
            for j in (nlevsoi + 1):nlevmaxurbgrnd
                col.levgrnd_class[c, j] = LEVGRND_CLASS_DEEP_BEDROCK
            end
        else
            for j in 1:nlevmaxurbgrnd
                col.levgrnd_class[c, j] = LEVGRND_CLASS_STANDARD
            end
        end
    end

    # Mark unused layers as ISPVAL
    for j in 1:nlevmaxurbgrnd
        for c in bounds.begc:bounds.endc
            if col.z[c, j + joff] == SPVAL
                col.levgrnd_class[c, j] = ISPVAL
            end
        end
    end

    nothing
end

"""
    init_topography!(bounds, col, surf)

Set topographic slope, std_elev, and micro_sigma for each column.
"""
function init_topography!(bounds::BoundsType, col::ColumnData,
                           surf::SurfaceInputData)
    slopebeta = _initvert_slopebeta[]
    slopemax  = _initvert_slopemax[]

    for c in bounds.begc:bounds.endc
        g = col.gridcell[c]
        gi = g - bounds.begg + 1
        col.topo_slope[c] = surf.slope[gi]
        col.topo_std[c] = surf.std_elev[gi]
    end

    # Compute micro_sigma
    slope0 = slopemax > 0.0 && slopebeta > 0.0 ?
        slopemax^(1.0 / slopebeta) : 0.0

    for c in bounds.begc:bounds.endc
        if col.is_hillslope_column[c]
            col.micro_sigma[c] = (atan(col.hill_slope[c]) + slope0)^slopebeta
        else
            col.micro_sigma[c] = (col.topo_slope[c] + slope0)^slopebeta
        end
    end

    nothing
end

"""
    hasBedrock(col_itype, lun_itype) -> Bool

Returns true if the given column type includes bedrock layers.
"""
function hasBedrock(col_itype::Int, lun_itype::Int)
    if lun_itype == ISTICE
        return false
    elseif ISTURB_MIN <= lun_itype <= ISTURB_MAX
        return col_itype == ICOL_ROAD_PERV
    else
        return true
    end
end

"""
    find_soil_layer_containing_depth(depth) -> Int

Find the soil layer that contains the given depth (m).
"""
function find_soil_layer_containing_depth(depth::Float64)
    zisoi_g = zisoi[]
    nlevgrnd = varpar.nlevgrnd

    depth > zisoi_g[1] || error("find_soil_layer_containing_depth: depth $depth above top of soil")

    for i in 1:nlevgrnd
        if depth <= zisoi_g[i + 1]
            return i
        end
    end

    error("find_soil_layer_containing_depth: depth $depth below bottom of soil")
end
