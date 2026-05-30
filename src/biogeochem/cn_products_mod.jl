# ==========================================================================
# Ported from: src/biogeochem/CNProductsMod.F90
# Calculate loss fluxes from wood products pools, and update product pool
# state variables.
#
# The type definition (CNProductsData) lives in src/types/cn_products.jl.
# This file contains all the science/infrastructure functions that operate
# on that type.
#
# Public functions:
#   cn_products_full_init!        -- Allocate the extended (full) product data
#   cn_products_set_values!       -- Zero gain-flux arrays each timestep
#   cn_products_update!           -- Top-level driver: partition + update pools
#   cn_products_compute_summary!  -- Compute summary variables (totals)
# ==========================================================================

# --------------------------------------------------------------------------
# Extended data type that holds the full set of product-pool fields
# (private states, patch-level fluxes, loss fluxes, gains).
# The compact CNProductsData in types/cn_products.jl only carries the
# three publicly-visible gridcell arrays.  For the module functions we
# need the complete Fortran type surface.
# --------------------------------------------------------------------------

"""
    CNProductsFullData

Extended product-pool data structure that holds *all* fields from the
Fortran `cn_products_type` (states, gain fluxes, loss fluxes, both at
patch and gridcell level).

The compact `CNProductsData` (in `src/types/cn_products.jl`) only carries
the three publicly-visible gridcell-level fields.

Ported from `cn_products_type` in `CNProductsMod.F90`.
"""
Base.@kwdef mutable struct CNProductsFullData
    # ---- Public states (same as CNProductsData) --------------------------
    cropprod1_grc      ::Vector{Float64} = Float64[]  # (g/m2) crop product pool, 1-year lifespan
    prod10_grc         ::Vector{Float64} = Float64[]  # (g/m2) wood product pool, 10-year lifespan
    prod100_grc        ::Vector{Float64} = Float64[]  # (g/m2) wood product pool, 100-year lifespan
    tot_woodprod_grc   ::Vector{Float64} = Float64[]  # (g/m2) total wood product pool
    product_loss_grc   ::Vector{Float64} = Float64[]  # (g/m2/s) total decomposition loss

    # ---- Gain fluxes (gridcell level) ------------------------------------
    dwt_prod10_gain_grc            ::Vector{Float64} = Float64[]
    dwt_prod100_gain_grc           ::Vector{Float64} = Float64[]
    dwt_woodprod_gain_grc          ::Vector{Float64} = Float64[]
    dwt_cropprod1_gain_grc         ::Vector{Float64} = Float64[]
    gru_prod10_gain_grc            ::Vector{Float64} = Float64[]
    gru_prod100_gain_grc           ::Vector{Float64} = Float64[]
    gru_woodprod_gain_grc          ::Vector{Float64} = Float64[]
    hrv_deadstem_to_prod10_grc     ::Vector{Float64} = Float64[]
    hrv_deadstem_to_prod100_grc    ::Vector{Float64} = Float64[]
    crop_harvest_to_cropprod1_grc  ::Vector{Float64} = Float64[]

    # ---- Gain fluxes (patch level) ---------------------------------------
    gru_prod10_gain_patch              ::Vector{Float64} = Float64[]
    gru_prod100_gain_patch             ::Vector{Float64} = Float64[]
    hrv_deadstem_to_prod10_patch       ::Vector{Float64} = Float64[]
    hrv_deadstem_to_prod100_patch      ::Vector{Float64} = Float64[]
    crop_harvest_to_cropprod1_patch    ::Vector{Float64} = Float64[]

    # ---- Loss fluxes (gridcell level) ------------------------------------
    cropprod1_loss_grc       ::Vector{Float64} = Float64[]
    prod10_loss_grc          ::Vector{Float64} = Float64[]
    prod100_loss_grc         ::Vector{Float64} = Float64[]
    tot_woodprod_loss_grc    ::Vector{Float64} = Float64[]
end

# --------------------------------------------------------------------------
# cn_products_full_init!
# --------------------------------------------------------------------------

"""
    cn_products_full_init!(prod, ng, np)

Allocate and zero-initialize all fields in `CNProductsFullData` for
`ng` gridcells and `np` patches.

Ported from `InitAllocate` + `InitCold` in `CNProductsMod.F90`.
"""
function cn_products_full_init!(prod::CNProductsFullData, ng::Int, np::Int)
    # ---- Gridcell-level states -------------------------------------------
    prod.cropprod1_grc      = zeros(ng)
    prod.prod10_grc         = zeros(ng)
    prod.prod100_grc        = zeros(ng)
    prod.tot_woodprod_grc   = zeros(ng)
    prod.product_loss_grc   = zeros(ng)

    # ---- Gridcell-level gain fluxes --------------------------------------
    prod.dwt_prod10_gain_grc           = zeros(ng)
    prod.dwt_prod100_gain_grc          = zeros(ng)
    prod.dwt_woodprod_gain_grc         = zeros(ng)
    prod.dwt_cropprod1_gain_grc        = zeros(ng)
    prod.gru_prod10_gain_grc           = zeros(ng)
    prod.gru_prod100_gain_grc          = zeros(ng)
    prod.gru_woodprod_gain_grc         = zeros(ng)
    prod.hrv_deadstem_to_prod10_grc    = zeros(ng)
    prod.hrv_deadstem_to_prod100_grc   = zeros(ng)
    prod.crop_harvest_to_cropprod1_grc = zeros(ng)

    # ---- Patch-level gain fluxes -----------------------------------------
    prod.gru_prod10_gain_patch             = zeros(np)
    prod.gru_prod100_gain_patch            = zeros(np)
    prod.hrv_deadstem_to_prod10_patch      = zeros(np)
    prod.hrv_deadstem_to_prod100_patch     = zeros(np)
    prod.crop_harvest_to_cropprod1_patch   = zeros(np)

    # ---- Gridcell-level loss fluxes --------------------------------------
    prod.cropprod1_loss_grc    = zeros(ng)
    prod.prod10_loss_grc       = zeros(ng)
    prod.prod100_loss_grc      = zeros(ng)
    prod.tot_woodprod_loss_grc = zeros(ng)

    return nothing
end

# --------------------------------------------------------------------------
# cn_products_set_values!
# --------------------------------------------------------------------------

# GPU/CPU kernel: zero (or set) the six per-gridcell gain-flux arrays. Each
# gridcell writes only its own index — fully independent.
@kernel function _cnprod_set_values_kernel!(dwt_prod10, dwt_prod100, dwt_cropprod1,
                                            crop_harv_cropprod1, hrv_prod10,
                                            hrv_prod100, setval)
    g = @index(Global)
    @inbounds begin
        dwt_prod10[g]          = setval
        dwt_prod100[g]         = setval
        dwt_cropprod1[g]       = setval
        crop_harv_cropprod1[g] = setval
        hrv_prod10[g]          = setval
        hrv_prod100[g]         = setval
    end
end

cnprod_set_values!(dwt_prod10, dwt_prod100, dwt_cropprod1, crop_harv_cropprod1,
                   hrv_prod10, hrv_prod100, setval) =
    _launch!(_cnprod_set_values_kernel!, dwt_prod10, dwt_prod100, dwt_cropprod1,
             crop_harv_cropprod1, hrv_prod10, hrv_prod100, setval)

"""
    cn_products_set_values!(prod, bounds, setval)

Set the gain-flux arrays that are incremented each timestep to `setval`.
This is typically called with `setval = 0.0` at the start of each timestep.

Ported from `SetValues` in `CNProductsMod.F90`.
"""
function cn_products_set_values!(prod::CNProductsFullData, bounds::BoundsType,
                                  setval::Real)
    # Assumes a single clump with begg == 1 (ndrange = length of gridcell arrays).
    cnprod_set_values!(prod.dwt_prod10_gain_grc,
                       prod.dwt_prod100_gain_grc,
                       prod.dwt_cropprod1_gain_grc,
                       prod.crop_harvest_to_cropprod1_grc,
                       prod.hrv_deadstem_to_prod10_grc,
                       prod.hrv_deadstem_to_prod100_grc,
                       Float64(setval))
    return nothing
end

# --------------------------------------------------------------------------
# cn_products_update! (top-level driver)
# --------------------------------------------------------------------------

"""
    cn_products_update!(prod, bounds, num_soilp, filter_soilp,
                        dwt_wood_product_gain_patch,
                        gru_wood_product_gain_patch,
                        wood_harvest_patch,
                        dwt_crop_product_gain_patch,
                        crop_harvest_to_cropprod_patch;
                        pprod10, pprod100, pprodharv10,
                        patch_gridcell, patch_itype,
                        pch, col, lun)

Update all loss fluxes from wood and crop product pools, and update product
pool state variables for both loss and gain terms.

This calls the internal partition routines and then
`cn_products_compute_product_summary!` to advance the state.

**Keyword arguments for subgrid data** (`pch`, `col`, `lun`) are only
required when `p2g` averaging is desired; when they are omitted a
simplified direct-summation path is used (suitable for testing / single-
gridcell setups where every patch maps 1-to-1 to its gridcell).

Ported from `UpdateProducts` in `CNProductsMod.F90`.
"""
function cn_products_update!(prod::CNProductsFullData,
                              bounds::BoundsType,
                              num_soilp::Int,
                              filter_soilp::Vector{Int},
                              dwt_wood_product_gain_patch::Vector{<:Real},
                              gru_wood_product_gain_patch::Vector{<:Real},
                              wood_harvest_patch::Vector{<:Real},
                              dwt_crop_product_gain_patch::Vector{<:Real},
                              crop_harvest_to_cropprod_patch::Vector{<:Real};
                              pprod10::Vector{<:Real},
                              pprod100::Vector{<:Real},
                              pprodharv10::Vector{<:Real},
                              patch_gridcell::Vector{Int},
                              patch_itype::Vector{Int},
                              pch::PatchData,
                              col::ColumnData,
                              lun::LandunitData)

    cn_products_partition_wood_fluxes!(prod, bounds,
        num_soilp, filter_soilp,
        dwt_wood_product_gain_patch,
        gru_wood_product_gain_patch,
        wood_harvest_patch;
        pprod10 = pprod10,
        pprod100 = pprod100,
        pprodharv10 = pprodharv10,
        patch_gridcell = patch_gridcell,
        patch_itype = patch_itype,
        pch = pch, col = col, lun = lun)

    cn_products_partition_crop_fluxes!(prod, bounds,
        num_soilp, filter_soilp,
        dwt_crop_product_gain_patch,
        crop_harvest_to_cropprod_patch;
        patch_gridcell = patch_gridcell,
        pch = pch, col = col, lun = lun)

    return nothing
end

# --------------------------------------------------------------------------
# cn_products_partition_wood_fluxes! (private helper)
# --------------------------------------------------------------------------

"""
    cn_products_partition_wood_fluxes!(prod, bounds, ...)

Partition input wood fluxes into 10-year and 100-year product pools.

Ported from `PartitionWoodFluxes` in `CNProductsMod.F90`.
"""
function cn_products_partition_wood_fluxes!(prod::CNProductsFullData,
        bounds::BoundsType,
        num_soilp::Int,
        filter_soilp::Vector{Int},
        dwt_wood_product_gain_patch::Vector{<:Real},
        gru_wood_product_gain_patch::Vector{<:Real},
        wood_harvest_patch::Vector{<:Real};
        pprod10::Vector{<:Real},
        pprod100::Vector{<:Real},
        pprodharv10::Vector{<:Real},
        patch_gridcell::Vector{Int},
        patch_itype::Vector{Int},
        pch::PatchData,
        col::ColumnData,
        lun::LandunitData)

    # ---- Partition patch-level gross unrepresented fluxes ----
    for fp in 1:num_soilp
        p = filter_soilp[fp]
        ivt = patch_itype[p] + 1  # 0-based Fortran → 1-based Julia

        pp10  = pprod10[ivt]
        pp100 = pprod100[ivt]
        pp_tot = pp10 + pp100

        if pp_tot > 0.0
            pprod10_frac  = pp10  / pp_tot
            pprod100_frac = pp100 / pp_tot
        else
            pprod10_frac  = 0.0
            pprod100_frac = 0.0
        end

        prod.gru_prod10_gain_patch[p]  = gru_wood_product_gain_patch[p] * pprod10_frac
        prod.gru_prod100_gain_patch[p] = gru_wood_product_gain_patch[p] * pprod100_frac
    end

    # ---- Average gross unrepresented fluxes from patch to gridcell ----
    p2g_1d!(prod.gru_prod10_gain_grc,
            prod.gru_prod10_gain_patch,
            bounds, "unity", "unity", "unity",
            pch, col, lun)
    p2g_1d!(prod.gru_prod100_gain_grc,
            prod.gru_prod100_gain_patch,
            bounds, "unity", "unity", "unity",
            pch, col, lun)

    # ---- Partition patch-level harvest fluxes ----
    for fp in 1:num_soilp
        p = filter_soilp[fp]
        ivt = patch_itype[p] + 1  # 0-based Fortran → 1-based Julia

        prod.hrv_deadstem_to_prod10_patch[p]  =
            wood_harvest_patch[p] * pprodharv10[ivt]
        prod.hrv_deadstem_to_prod100_patch[p] =
            wood_harvest_patch[p] * (1.0 - pprodharv10[ivt])
    end

    # ---- Average harvest fluxes from patch to gridcell ----
    p2g_1d!(prod.hrv_deadstem_to_prod10_grc,
            prod.hrv_deadstem_to_prod10_patch,
            bounds, "unity", "unity", "unity",
            pch, col, lun)
    p2g_1d!(prod.hrv_deadstem_to_prod100_grc,
            prod.hrv_deadstem_to_prod100_patch,
            bounds, "unity", "unity", "unity",
            pch, col, lun)

    # ---- Partition dynamic land cover fluxes (dwt) ----
    # Note: dwt fluxes are already expressed per unit gridcell area, so we
    # simply sum the patch contributions to the owning gridcell.
    for p in bounds.begp:bounds.endp
        g   = patch_gridcell[p]
        ivt = patch_itype[p] + 1  # 0-based Fortran → 1-based Julia

        pp10  = pprod10[ivt]
        pp100 = pprod100[ivt]
        pp_tot = pp10 + pp100

        if pp_tot > 0.0
            pprod10_frac  = pp10  / pp_tot
            pprod100_frac = pp100 / pp_tot

            prod.dwt_prod10_gain_grc[g]  += dwt_wood_product_gain_patch[p] * pprod10_frac
            prod.dwt_prod100_gain_grc[g] += dwt_wood_product_gain_patch[p] * pprod100_frac
        elseif dwt_wood_product_gain_patch[p] > 0.0
            error("cn_products_partition_wood_fluxes!: " *
                  "dwt_wood_product_gain_patch[p] > 0 but pprod_tot == 0 at p=$p")
        end
    end

    return nothing
end

# --------------------------------------------------------------------------
# cn_products_partition_crop_fluxes! (private helper)
# --------------------------------------------------------------------------

"""
    cn_products_partition_crop_fluxes!(prod, bounds, ...)

Partition input crop fluxes into crop product pools. Currently only a
single 1-year crop product pool exists.

Ported from `PartitionCropFluxes` in `CNProductsMod.F90`.
"""
function cn_products_partition_crop_fluxes!(prod::CNProductsFullData,
        bounds::BoundsType,
        num_soilp::Int,
        filter_soilp::Vector{Int},
        dwt_crop_product_gain_patch::Vector{<:Real},
        crop_harvest_to_cropprod_patch::Vector{<:Real};
        patch_gridcell::Vector{Int},
        pch::PatchData,
        col::ColumnData,
        lun::LandunitData)

    # ---- Determine gains from crop harvest ----
    for fp in 1:num_soilp
        p = filter_soilp[fp]
        # All crop product goes into the 1-year pool
        prod.crop_harvest_to_cropprod1_patch[p] = crop_harvest_to_cropprod_patch[p]
    end

    # ---- Average crop harvest from patch to gridcell ----
    p2g_1d!(prod.crop_harvest_to_cropprod1_grc,
            prod.crop_harvest_to_cropprod1_patch,
            bounds, "unity", "unity", "unity",
            pch, col, lun)

    # ---- Determine gains from dynamic landcover ----
    # dwt fluxes are per-gridcell area -> simply sum patch contributions
    for p in bounds.begp:bounds.endp
        g = patch_gridcell[p]
        prod.dwt_cropprod1_gain_grc[g] += dwt_crop_product_gain_patch[p]
    end

    return nothing
end

# --------------------------------------------------------------------------
# cn_products_compute_product_summary!
# --------------------------------------------------------------------------

"""
    cn_products_compute_product_summary!(prod, bounds, dt)

Calculate losses from product pools, then update pool state variables by
adding gains and subtracting losses over one timestep of length `dt` seconds.

Decay constants (1/s) give ~90% loss of initial state over 1, 10, and
100 years using a discrete-time fractional decay algorithm.

Ported from `ComputeProductSummaryVars` in `CNProductsMod.F90`.
"""
# GPU/CPU kernel: per-gridcell decay losses + pool state update. Each gridcell
# touches only its own index in every array (loss computed then applied within
# the same thread), so iterations are fully independent.
@kernel function _cnprod_summary_kernel!(cropprod1_grc, prod10_grc, prod100_grc,
                                         cropprod1_loss, prod10_loss, prod100_loss,
                                         @Const(dwt_cropprod1_gain),
                                         @Const(dwt_prod10_gain),
                                         @Const(dwt_prod100_gain),
                                         @Const(gru_prod10_gain),
                                         @Const(gru_prod100_gain),
                                         @Const(crop_harvest_to_cropprod1),
                                         @Const(hrv_deadstem_to_prod10),
                                         @Const(hrv_deadstem_to_prod100),
                                         kprod1, kprod10, kprod100, dt)
    g = @index(Global)
    @inbounds begin
        # ---- Calculate loss fluxes (g/m2/s) ----
        cropprod1_loss[g] = cropprod1_grc[g] * kprod1
        prod10_loss[g]    = prod10_grc[g]    * kprod10
        prod100_loss[g]   = prod100_grc[g]   * kprod100

        # ---- Gains from dynamic landcover change ----
        cropprod1_grc[g] += dwt_cropprod1_gain[g] * dt
        prod10_grc[g]    += dwt_prod10_gain[g]    * dt
        prod100_grc[g]   += dwt_prod100_gain[g]   * dt

        # ---- Gains from gross unrepresented landcover change ----
        prod10_grc[g]  += gru_prod10_gain[g]  * dt
        prod100_grc[g] += gru_prod100_gain[g] * dt

        # ---- Gains from harvest ----
        cropprod1_grc[g] += crop_harvest_to_cropprod1[g] * dt
        prod10_grc[g]    += hrv_deadstem_to_prod10[g]    * dt
        prod100_grc[g]   += hrv_deadstem_to_prod100[g]   * dt

        # ---- Losses from decomposition ----
        cropprod1_grc[g] -= cropprod1_loss[g] * dt
        prod10_grc[g]    -= prod10_loss[g]    * dt
        prod100_grc[g]   -= prod100_loss[g]   * dt
    end
end

function cnprod_compute_product_summary!(cropprod1_grc, prod10_grc, prod100_grc,
                                         cropprod1_loss, prod10_loss, prod100_loss,
                                         dwt_cropprod1_gain, dwt_prod10_gain,
                                         dwt_prod100_gain, gru_prod10_gain,
                                         gru_prod100_gain, crop_harvest_to_cropprod1,
                                         hrv_deadstem_to_prod10,
                                         hrv_deadstem_to_prod100,
                                         kprod1, kprod10, kprod100, dt)
    _launch!(_cnprod_summary_kernel!, cropprod1_grc, prod10_grc, prod100_grc,
             cropprod1_loss, prod10_loss, prod100_loss, dwt_cropprod1_gain,
             dwt_prod10_gain, dwt_prod100_gain, gru_prod10_gain, gru_prod100_gain,
             crop_harvest_to_cropprod1, hrv_deadstem_to_prod10,
             hrv_deadstem_to_prod100, kprod1, kprod10, kprod100, dt)
end

function cn_products_compute_product_summary!(prod::CNProductsFullData,
                                               bounds::BoundsType,
                                               dt::Real)
    # Decay constants (1/s) for 1-yr, 10-yr, 100-yr pools
    kprod1   = 7.2e-8
    kprod10  = 7.2e-9
    kprod100 = 7.2e-10

    # Assumes a single clump with begg == 1 (ndrange = length of gridcell arrays).
    cnprod_compute_product_summary!(prod.cropprod1_grc, prod.prod10_grc,
        prod.prod100_grc, prod.cropprod1_loss_grc, prod.prod10_loss_grc,
        prod.prod100_loss_grc, prod.dwt_cropprod1_gain_grc,
        prod.dwt_prod10_gain_grc, prod.dwt_prod100_gain_grc,
        prod.gru_prod10_gain_grc, prod.gru_prod100_gain_grc,
        prod.crop_harvest_to_cropprod1_grc, prod.hrv_deadstem_to_prod10_grc,
        prod.hrv_deadstem_to_prod100_grc, kprod1, kprod10, kprod100, Float64(dt))

    return nothing
end

# --------------------------------------------------------------------------
# cn_products_compute_summary!
# --------------------------------------------------------------------------

"""
    cn_products_compute_summary!(prod, bounds)

Compute summary variables: sums across multiple product pools.

Ported from `ComputeSummaryVars` in `CNProductsMod.F90`.
"""
# GPU/CPU kernel: per-gridcell summary sums. Each gridcell writes only its own
# index from read-only inputs — fully independent.
@kernel function _cnprod_compute_summary_kernel!(tot_woodprod, tot_woodprod_loss,
                                                 product_loss, dwt_woodprod_gain,
                                                 gru_woodprod_gain,
                                                 @Const(prod10_grc),
                                                 @Const(prod100_grc),
                                                 @Const(prod10_loss),
                                                 @Const(prod100_loss),
                                                 @Const(cropprod1_loss),
                                                 @Const(dwt_prod10_gain),
                                                 @Const(dwt_prod100_gain),
                                                 @Const(gru_prod10_gain),
                                                 @Const(gru_prod100_gain))
    g = @index(Global)
    @inbounds begin
        tot_woodprod[g]      = prod10_grc[g] + prod100_grc[g]
        tot_woodprod_loss[g] = prod10_loss[g] + prod100_loss[g]
        product_loss[g]      = cropprod1_loss[g] + prod10_loss[g] + prod100_loss[g]
        dwt_woodprod_gain[g] = dwt_prod10_gain[g] + dwt_prod100_gain[g]
        gru_woodprod_gain[g] = gru_prod10_gain[g] + gru_prod100_gain[g]
    end
end

function cnprod_compute_summary!(tot_woodprod, tot_woodprod_loss, product_loss,
                                 dwt_woodprod_gain, gru_woodprod_gain, prod10_grc,
                                 prod100_grc, prod10_loss, prod100_loss,
                                 cropprod1_loss, dwt_prod10_gain, dwt_prod100_gain,
                                 gru_prod10_gain, gru_prod100_gain)
    _launch!(_cnprod_compute_summary_kernel!, tot_woodprod, tot_woodprod_loss,
             product_loss, dwt_woodprod_gain, gru_woodprod_gain, prod10_grc,
             prod100_grc, prod10_loss, prod100_loss, cropprod1_loss,
             dwt_prod10_gain, dwt_prod100_gain, gru_prod10_gain, gru_prod100_gain)
end

function cn_products_compute_summary!(prod::CNProductsFullData,
                                       bounds::BoundsType)
    # Assumes a single clump with begg == 1 (ndrange = length of gridcell arrays).
    cnprod_compute_summary!(prod.tot_woodprod_grc, prod.tot_woodprod_loss_grc,
        prod.product_loss_grc, prod.dwt_woodprod_gain_grc,
        prod.gru_woodprod_gain_grc, prod.prod10_grc, prod.prod100_grc,
        prod.prod10_loss_grc, prod.prod100_loss_grc, prod.cropprod1_loss_grc,
        prod.dwt_prod10_gain_grc, prod.dwt_prod100_gain_grc,
        prod.gru_prod10_gain_grc, prod.gru_prod100_gain_grc)

    return nothing
end
