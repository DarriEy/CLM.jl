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

"""
    cn_products_set_values!(prod, bounds, setval)

Set the gain-flux arrays that are incremented each timestep to `setval`.
This is typically called with `setval = 0.0` at the start of each timestep.

Ported from `SetValues` in `CNProductsMod.F90`.
"""
function cn_products_set_values!(prod::CNProductsFullData, bounds::BoundsType,
                                  setval::Float64)
    for g in bounds.begg:bounds.endg
        prod.dwt_prod10_gain_grc[g]           = setval
        prod.dwt_prod100_gain_grc[g]          = setval
        prod.dwt_cropprod1_gain_grc[g]        = setval
        prod.crop_harvest_to_cropprod1_grc[g] = setval
        prod.hrv_deadstem_to_prod10_grc[g]    = setval
        prod.hrv_deadstem_to_prod100_grc[g]   = setval
    end
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
                              dwt_wood_product_gain_patch::Vector{Float64},
                              gru_wood_product_gain_patch::Vector{Float64},
                              wood_harvest_patch::Vector{Float64},
                              dwt_crop_product_gain_patch::Vector{Float64},
                              crop_harvest_to_cropprod_patch::Vector{Float64};
                              pprod10::Vector{Float64},
                              pprod100::Vector{Float64},
                              pprodharv10::Vector{Float64},
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
        dwt_wood_product_gain_patch::Vector{Float64},
        gru_wood_product_gain_patch::Vector{Float64},
        wood_harvest_patch::Vector{Float64};
        pprod10::Vector{Float64},
        pprod100::Vector{Float64},
        pprodharv10::Vector{Float64},
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
        dwt_crop_product_gain_patch::Vector{Float64},
        crop_harvest_to_cropprod_patch::Vector{Float64};
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
function cn_products_compute_product_summary!(prod::CNProductsFullData,
                                               bounds::BoundsType,
                                               dt::Float64)
    # Decay constants (1/s) for 1-yr, 10-yr, 100-yr pools
    kprod1   = 7.2e-8
    kprod10  = 7.2e-9
    kprod100 = 7.2e-10

    for g in bounds.begg:bounds.endg
        # ---- Calculate loss fluxes (g/m2/s) ----
        prod.cropprod1_loss_grc[g] = prod.cropprod1_grc[g] * kprod1
        prod.prod10_loss_grc[g]    = prod.prod10_grc[g]    * kprod10
        prod.prod100_loss_grc[g]   = prod.prod100_grc[g]   * kprod100
    end

    for g in bounds.begg:bounds.endg
        # ---- Gains from dynamic landcover change ----
        prod.cropprod1_grc[g] += prod.dwt_cropprod1_gain_grc[g] * dt
        prod.prod10_grc[g]    += prod.dwt_prod10_gain_grc[g]    * dt
        prod.prod100_grc[g]   += prod.dwt_prod100_gain_grc[g]   * dt

        # ---- Gains from gross unrepresented landcover change ----
        prod.prod10_grc[g]  += prod.gru_prod10_gain_grc[g]  * dt
        prod.prod100_grc[g] += prod.gru_prod100_gain_grc[g] * dt

        # ---- Gains from harvest ----
        prod.cropprod1_grc[g] += prod.crop_harvest_to_cropprod1_grc[g] * dt
        prod.prod10_grc[g]    += prod.hrv_deadstem_to_prod10_grc[g]    * dt
        prod.prod100_grc[g]   += prod.hrv_deadstem_to_prod100_grc[g]   * dt

        # ---- Losses from decomposition ----
        prod.cropprod1_grc[g] -= prod.cropprod1_loss_grc[g] * dt
        prod.prod10_grc[g]    -= prod.prod10_loss_grc[g]    * dt
        prod.prod100_grc[g]   -= prod.prod100_loss_grc[g]   * dt
    end

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
function cn_products_compute_summary!(prod::CNProductsFullData,
                                       bounds::BoundsType)
    for g in bounds.begg:bounds.endg
        # Total wood products
        prod.tot_woodprod_grc[g] =
            prod.prod10_grc[g] + prod.prod100_grc[g]

        # Total loss from wood products
        prod.tot_woodprod_loss_grc[g] =
            prod.prod10_loss_grc[g] + prod.prod100_loss_grc[g]

        # Total loss from ALL products
        prod.product_loss_grc[g] =
            prod.cropprod1_loss_grc[g] +
            prod.prod10_loss_grc[g] +
            prod.prod100_loss_grc[g]

        # Summary gain fluxes
        prod.dwt_woodprod_gain_grc[g] =
            prod.dwt_prod10_gain_grc[g] + prod.dwt_prod100_gain_grc[g]

        prod.gru_woodprod_gain_grc[g] =
            prod.gru_prod10_gain_grc[g] + prod.gru_prod100_gain_grc[g]
    end

    return nothing
end
