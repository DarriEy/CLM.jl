# ==========================================================================
# Ported from: src/biogeochem/CNProductsMod.F90
# Carbon and Nitrogen product pools (wood products, crop products)
# ==========================================================================

"""
    CNProductsData

Carbon or nitrogen product pool data at the gridcell level.
Holds wood and crop product pool states and fluxes.

Ported from `cn_products_type` in `CNProductsMod.F90`.
"""
Base.@kwdef mutable struct CNProductsData{FT<:Real}
    cropprod1_grc      ::Vector{FT} = Float64[]  # (g/m2) crop product pool (1-year)
    tot_woodprod_grc   ::Vector{FT} = Float64[]  # (g/m2) total wood product pools
    product_loss_grc   ::Vector{FT} = Float64[]  # (g/m2/s) losses from wood & crop products
end

"""
    cn_products_init!(prod::CNProductsData, ng::Int)

Allocate and initialize CNProductsData for `ng` gridcells.
"""
function cn_products_init!(prod::CNProductsData{FT}, ng::Int) where {FT}
    prod.cropprod1_grc    = fill(zero(FT), ng)
    prod.tot_woodprod_grc = fill(zero(FT), ng)
    prod.product_loss_grc = fill(zero(FT), ng)
    return nothing
end
