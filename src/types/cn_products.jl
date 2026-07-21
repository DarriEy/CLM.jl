# ==========================================================================
# Ported from: src/biogeochem/CNProductsMod.F90
# Carbon and Nitrogen product pools (wood products, crop products)
# ==========================================================================

"""
    AbstractCNProducts

Common supertype for the two product-pool data structures:
[`CNProductsData`](@ref) (the compact three-field gridcell view consumed by the
CN balance check) and `CNProductsFullData` (the full Fortran `cn_products_type`
surface, defined in `biogeochem/cn_products_mod.jl`, that `cn_products_update!`
advances). Both expose `cropprod1_grc`, `tot_woodprod_grc`, and
`product_loss_grc`, so the balance-check routines dispatch on this supertype.
"""
abstract type AbstractCNProducts end

"""
    CNProductsData

Carbon or nitrogen product pool data at the gridcell level.
Holds wood and crop product pool states and fluxes.

Ported from `cn_products_type` in `CNProductsMod.F90`.
"""
Base.@kwdef mutable struct CNProductsData{FT<:Real,
                              V<:AbstractVector{FT}} <: AbstractCNProducts
    cropprod1_grc      ::V = Float64[]  # (g/m2) crop product pool (1-year)
    tot_woodprod_grc   ::V = Float64[]  # (g/m2) total wood product pools
    product_loss_grc   ::V = Float64[]  # (g/m2/s) losses from wood & crop products
end

CNProductsData{FT}(; kwargs...) where {FT<:Real} =
    CNProductsData{FT, Vector{FT}}(; kwargs...)
Adapt.@adapt_structure CNProductsData


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
