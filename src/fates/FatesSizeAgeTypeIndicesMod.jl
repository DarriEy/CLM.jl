# =====================================================================================
# FatesSizeAgeTypeIndicesMod — Julia port of
#   src/fates/main/FatesSizeAgeTypeIndicesMod.F90
#
# Pure index-mapping helpers for the FATES history/diagnostic multiplexed
# dimensions: bin a continuous quantity (dbh, patch age, cohort age, height) into
# a 1-based class index, and combine (size, pft), (age), (height), (fuel),
# (damage), (canopy layer) etc. into a single flat 1-D index (and back, in the
# analysis-code comments preserved from the Fortran source).
#
# Dependencies (all from the merged FATES foundation + Batch-1):
#   - FatesInterfaceTypesMod : nlevsclass, nlevage, nlevheight, nlevcoage,
#                              nlevdamage  (Ref{Int}, dereferenced with `[]`)
#   - EDParamsMod            : nclmax (plain const) and the runtime bin-edge
#                              vectors exposed on `ed_params()`:
#                                ED_val_history_sizeclass_bin_edges
#                                ED_val_history_ageclass_bin_edges
#                                ED_val_history_height_bin_edges
#                                ED_val_history_coageclass_bin_edges
#                                ED_val_history_damage_bin_edges
#
# All modules are `include`d directly into the top-level `CLM` module, so the
# above names are referenced unqualified. Fortran `fates_r8` -> Float64.
#
# NOTE on the bin lookup: Fortran's
#     count(value - bin_edges .ge. 0.0_r8)
# counts how many edge values are <= `value`. With monotonically increasing
# edges whose first edge is the minimum, this yields a 1-based class index:
# below the first edge -> 0, in [edge[k], edge[k+1]) -> k, at/above the last
# edge -> nedges. We reproduce that semantics exactly (no clamping).
# =====================================================================================

"""
    _count_le(value, edges)

Count how many `edges` satisfy `value - edge >= 0`, i.e. how many edges are
`<= value`. Direct translation of the Fortran `count(value - edges .ge. 0)`
idiom used for every size/age/height/coage bin lookup. Returns a 1-based class
index for monotonically increasing edge arrays.
"""
@inline function _count_le(value::Float64, edges::AbstractVector{Float64})
    n = 0
    @inbounds for e in edges
        if value - e >= 0.0
            n += 1
        end
    end
    return n
end

# =====================================================================================

"""
    get_age_class_index(age) -> Int

Patch-age class index from continuous patch age (`get_age_class_index`).
"""
function get_age_class_index(age::Float64)::Int
    return _count_le(age, ed_params().ED_val_history_ageclass_bin_edges)
end

# =====================================================================================

"""
    get_size_class_index(dbh) -> Int

Cohort size (dbh) class index (`get_size_class_index`).
"""
function get_size_class_index(dbh::Float64)::Int
    return _count_le(dbh, ed_params().ED_val_history_sizeclass_bin_edges)
end

# =====================================================================================

"""
    get_coage_class_index(coage) -> Int

Cohort-age class index (`get_coage_class_index`).
"""
function get_coage_class_index(coage::Float64)::Int
    return _count_le(coage, ed_params().ED_val_history_coageclass_bin_edges)
end

# =====================================================================================

"""
    get_height_index(height) -> Int

Cohort height class index (`get_height_index`).
"""
function get_height_index(height::Float64)::Int
    return _count_le(height, ed_params().ED_val_history_height_bin_edges)
end

# =====================================================================================

"""
    get_layersizetype_class_index(layer, dbh, pft) -> Int

Flat 1-D index for a (canopy layer, size class, pft) triplet
(`get_layersizetype_class_index`).

Reverse mapping (1-based), preserved from the Fortran comment:
    pft        = ceil(index / (nlevsclass*nclmax))
    size_class = ceil((index - (pft-1)*nlevsclass*nclmax) / nclmax)
    layer      = index - ((pft-1)*nlevsclass*nclmax + (size_class-1)*nclmax)
"""
function get_layersizetype_class_index(layer::Int, dbh::Float64, pft::Int)::Int
    size_class = get_size_class_index(dbh)
    return (pft - 1) * nlevsclass[] * nclmax + (size_class - 1) * nclmax + layer
end

# =====================================================================================

"""
    get_sizeage_class_index(dbh, age) -> Int

Flat 1-D index for a (size class, age class) pair (`get_sizeage_class_index`).
"""
function get_sizeage_class_index(dbh::Float64, age::Float64)::Int
    size_class = get_size_class_index(dbh)
    age_class  = get_age_class_index(age)
    return (age_class - 1) * nlevsclass[] + size_class
end

# =====================================================================================

"""
    get_cdamagesize_class_index(dbh, cdamage) -> Int

Flat 1-D index for a (crown-damage class, size class) pair
(`get_cdamagesize_class_index`). `cdamage` is already an integer class.
"""
function get_cdamagesize_class_index(dbh::Float64, cdamage::Int)::Int
    size_class = get_size_class_index(dbh)
    return (cdamage - 1) * nlevsclass[] + size_class
end

# =====================================================================================

"""
    sizetype_class_index(dbh, pft) -> (size_class, size_by_pft_class)

Fortran subroutine `sizetype_class_index` with two `intent(out)` args; returned
here as a tuple `(size_class, size_by_pft_class)`.
"""
function sizetype_class_index(dbh::Float64, pft::Int)
    size_class        = get_size_class_index(dbh)
    size_by_pft_class = (pft - 1) * nlevsclass[] + size_class
    return size_class, size_by_pft_class
end

# =====================================================================================

"""
    coagetype_class_index(coage, pft) -> (coage_class, coage_by_pft_class)

Fortran subroutine `coagetype_class_index` with two `intent(out)` args; returned
here as a tuple `(coage_class, coage_by_pft_class)`.
"""
function coagetype_class_index(coage::Float64, pft::Int)
    coage_class        = get_coage_class_index(coage)
    coage_by_pft_class = (pft - 1) * nlevcoage[] + coage_class
    return coage_class, coage_by_pft_class
end

# =====================================================================================

"""
    get_sizeagepft_class_index(dbh, age, pft) -> Int

Flat 1-D index for a (size class, age class, pft) triplet
(`get_sizeagepft_class_index`).
"""
function get_sizeagepft_class_index(dbh::Float64, age::Float64, pft::Int)::Int
    size_class = get_size_class_index(dbh)
    age_class  = get_age_class_index(age)
    return (age_class - 1) * nlevsclass[] + size_class +
           (pft - 1) * nlevsclass[] * nlevage[]
end

# =====================================================================================

"""
    get_cdamagesizepft_class_index(dbh, cdamage, pft) -> Int

Flat 1-D index for a (size class, crown-damage class, pft) triplet
(`get_cdamagesizepft_class_index`).
"""
function get_cdamagesizepft_class_index(dbh::Float64, cdamage::Int, pft::Int)::Int
    size_class = get_size_class_index(dbh)
    return (cdamage - 1) * nlevsclass[] + size_class +
           (pft - 1) * nlevsclass[] * nlevdamage[]
end

# =====================================================================================

"""
    get_agepft_class_index(age, pft) -> Int

Flat 1-D index for an (age class, pft) pair (`get_agepft_class_index`).
"""
function get_agepft_class_index(age::Float64, pft::Int)::Int
    age_class = get_age_class_index(age)
    return age_class + (pft - 1) * nlevage[]
end

# =====================================================================================

"""
    get_agefuel_class_index(age, fuel) -> Int

Flat 1-D index for an (age class, fuel class) pair (`get_agefuel_class_index`).
`fuel` is already an integer class.
"""
function get_agefuel_class_index(age::Float64, fuel::Int)::Int
    age_class = get_age_class_index(age)
    return age_class + (fuel - 1) * nlevage[]
end
