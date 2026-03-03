# ==========================================================================
# Ported from: src/biogeochem/CNBalanceCheckMod.F90
# Carbon/nitrogen mass balance checking
# ==========================================================================

"""
    CNBalanceData

Data structure for carbon/nitrogen mass balance checking.
Stores beginning and end-of-timestep masses at column and gridcell levels,
plus warning/error thresholds.

Ported from `cn_balance_type` in `CNBalanceCheckMod.F90`.
"""
Base.@kwdef mutable struct CNBalanceData
    begcb_col ::Vector{Float64} = Float64[]  # (gC/m2) column carbon mass, beginning of time step
    endcb_col ::Vector{Float64} = Float64[]  # (gC/m2) column carbon mass, end of time step
    begnb_col ::Vector{Float64} = Float64[]  # (gN/m2) column nitrogen mass, beginning of time step
    endnb_col ::Vector{Float64} = Float64[]  # (gN/m2) column nitrogen mass, end of time step
    begcb_grc ::Vector{Float64} = Float64[]  # (gC/m2) gridcell carbon mass, beginning of time step
    endcb_grc ::Vector{Float64} = Float64[]  # (gC/m2) gridcell carbon mass, end of time step
    begnb_grc ::Vector{Float64} = Float64[]  # (gN/m2) gridcell nitrogen mass, beginning of time step
    endnb_grc ::Vector{Float64} = Float64[]  # (gN/m2) gridcell nitrogen mass, end of time step
    cwarning  ::Float64 = 1.0e-8   # (gC/m2) carbon balance warning threshold
    nwarning  ::Float64 = 1.0e-7   # (gN/m2) nitrogen balance warning threshold
    cerror    ::Float64 = 1.0e-7   # (gC/m2) carbon balance error threshold
    nerror    ::Float64 = 1.0e-3   # (gN/m2) nitrogen balance error threshold
end

"""
    cn_balance_init!(bal::CNBalanceData, nc::Int, ng::Int)

Allocate and initialize a `CNBalanceData` instance for `nc` columns and `ng` gridcells.

Ported from `Init` and `InitAllocate` in `CNBalanceCheckMod.F90`.
"""
function cn_balance_init!(bal::CNBalanceData, nc::Int, ng::Int)
    bal.begcb_col = fill(NaN, nc)
    bal.endcb_col = fill(NaN, nc)
    bal.begnb_col = fill(NaN, nc)
    bal.endnb_col = fill(NaN, nc)
    bal.begcb_grc = fill(NaN, ng)
    bal.endcb_grc = fill(NaN, ng)
    bal.begnb_grc = fill(NaN, ng)
    bal.endnb_grc = fill(NaN, ng)

    bal.cwarning = 1.0e-8
    bal.nwarning = 1.0e-7
    bal.nerror   = 1.0e-3
    bal.cerror   = 1.0e-7
    return nothing
end

# --------------------------------------------------------------------------
# Column-to-gridcell aggregation (unity weighting)
# --------------------------------------------------------------------------

"""
    c2g_unity!(garr, carr, col_gridcell, col_wtgcell, bounds_c, bounds_g)

Column-to-gridcell area-weighted average with unity scaling.
Simple reimplementation of `c2g` from `subgridAveMod.F90` used by balance checks.
"""
function c2g_unity!(
    garr::Vector{Float64},
    carr::Vector{Float64},
    col_gridcell::Vector{Int},
    col_wtgcell::Vector{Float64},
    bounds_c::UnitRange{Int},
    bounds_g::UnitRange{Int}
)
    for g in bounds_g
        garr[g] = 0.0
    end
    for c in bounds_c
        g = col_gridcell[c]
        if g in bounds_g
            garr[g] += carr[c] * col_wtgcell[c]
        end
    end
    return nothing
end

# --------------------------------------------------------------------------
# BeginCNGridcellBalance
# --------------------------------------------------------------------------

"""
    begin_cn_gridcell_balance!(bal, soilbgc_cstate, soilbgc_nstate,
        c_products, n_products, bounds_g;
        use_fates_bgc, hrv_xsmrpool_amount_left, gru_conv_cflux_amount_left,
        dwt_conv_cflux_amount_left)

Calculate beginning gridcell-level carbon/nitrogen balance for mass conservation check.

Should be called after CN state summaries have been computed and before
the dynamic landunit area updates.

Ported from `BeginCNGridcellBalance` in `CNBalanceCheckMod.F90`.

The dribbler amounts (`hrv_xsmrpool_amount_left`, `gru_conv_cflux_amount_left`,
`dwt_conv_cflux_amount_left`) must be pre-computed by the caller when
`use_fates_bgc` is false. When `use_fates_bgc` is true, they are ignored
and treated as zero.
"""
function begin_cn_gridcell_balance!(
    bal::CNBalanceData,
    soilbgc_cstate::SoilBiogeochemCarbonStateData,
    soilbgc_nstate::SoilBiogeochemNitrogenStateData,
    c_products::CNProductsData,
    n_products::CNProductsData,
    bounds_g::UnitRange{Int};
    use_fates_bgc::Bool=false,
    hrv_xsmrpool_amount_left::Vector{Float64}=Float64[],
    gru_conv_cflux_amount_left::Vector{Float64}=Float64[],
    dwt_conv_cflux_amount_left::Vector{Float64}=Float64[]
)
    totc = soilbgc_cstate.totc_grc
    totn = soilbgc_nstate.totn_grc
    c_cropprod1    = c_products.cropprod1_grc
    n_cropprod1    = n_products.cropprod1_grc
    c_tot_woodprod = c_products.tot_woodprod_grc
    n_tot_woodprod = n_products.tot_woodprod_grc

    for g in bounds_g
        if use_fates_bgc
            bal.begcb_grc[g] = totc[g] + c_tot_woodprod[g] + c_cropprod1[g]
        else
            bal.begcb_grc[g] = totc[g] + c_tot_woodprod[g] + c_cropprod1[g] +
                               hrv_xsmrpool_amount_left[g] +
                               gru_conv_cflux_amount_left[g] +
                               dwt_conv_cflux_amount_left[g]
        end
        bal.begnb_grc[g] = totn[g] + n_tot_woodprod[g] + n_cropprod1[g]
    end

    return nothing
end

# --------------------------------------------------------------------------
# BeginCNColumnBalance
# --------------------------------------------------------------------------

"""
    begin_cn_column_balance!(bal, soilbgc_cstate, soilbgc_nstate,
        mask_soil, bounds_c)

Calculate beginning column-level carbon/nitrogen balance for mass conservation check.

Should be called after CN state summaries have been recomputed for this time step
(after dynamic landunit area updates and associated filter updates).

Ported from `BeginCNColumnBalance` in `CNBalanceCheckMod.F90`.
"""
function begin_cn_column_balance!(
    bal::CNBalanceData,
    soilbgc_cstate::SoilBiogeochemCarbonStateData,
    soilbgc_nstate::SoilBiogeochemNitrogenStateData,
    mask_soil::BitVector,
    bounds_c::UnitRange{Int}
)
    totcolc = soilbgc_cstate.totc_col
    totcoln = soilbgc_nstate.totn_col

    for c in bounds_c
        mask_soil[c] || continue
        bal.begcb_col[c] = totcolc[c]
        bal.begnb_col[c] = totcoln[c]
    end

    return nothing
end

# --------------------------------------------------------------------------
# CBalanceCheck
# --------------------------------------------------------------------------

"""
    c_balance_check!(bal, soilbgc_cflux, soilbgc_cstate, cnveg_cflux,
        c_products, col_data, grc_data,
        mask_soil, bounds_c, bounds_g, dt;
        is_fates_col, use_fates_bgc,
        hrv_xsmrpool_amount_left, gru_conv_cflux_amount_left,
        dwt_conv_cflux_amount_left)

Perform carbon mass conservation check for column and gridcell.

Ported from `CBalanceCheck` in `CNBalanceCheckMod.F90`.

# Arguments
- `bal`: CNBalanceData with beginning-of-step masses already set
- `soilbgc_cflux`: SoilBiogeochemCarbonFluxData
- `soilbgc_cstate`: SoilBiogeochemCarbonStateData
- `cnveg_cflux`: CNVegCarbonFluxData
- `c_products`: CNProductsData for carbon
- `col_data`: ColumnData for subgrid mapping
- `grc_data`: GridcellData for lat/lon
- `mask_soil`: BitVector mask for soil columns
- `bounds_c`: column index range
- `bounds_g`: gridcell index range
- `dt`: timestep size (seconds)
- `is_fates_col`: per-column vector indicating FATES columns
- `use_fates_bgc`: global flag for FATES BGC
- `hrv_xsmrpool_amount_left`: dribbler amount (end-of-step, gridcell)
- `gru_conv_cflux_amount_left`: dribbler amount (end-of-step, gridcell)
- `dwt_conv_cflux_amount_left`: dribbler amount (end-of-step, gridcell)
"""
function c_balance_check!(
    bal::CNBalanceData,
    soilbgc_cflux::SoilBiogeochemCarbonFluxData,
    soilbgc_cstate::SoilBiogeochemCarbonStateData,
    cnveg_cflux::CNVegCarbonFluxData,
    c_products::CNProductsData,
    col_data::ColumnData,
    grc_data::GridcellData,
    mask_soil::BitVector,
    bounds_c::UnitRange{Int},
    bounds_g::UnitRange{Int},
    dt::Float64;
    is_fates_col::Vector{Bool}=Bool[],
    use_fates_bgc::Bool=false,
    hrv_xsmrpool_amount_left::Vector{Float64}=Float64[],
    gru_conv_cflux_amount_left::Vector{Float64}=Float64[],
    dwt_conv_cflux_amount_left::Vector{Float64}=Float64[]
)
    # Local aliases for column-level fields
    col_begcb = bal.begcb_col
    col_endcb = bal.endcb_col
    grc_begcb = bal.begcb_grc
    grc_endcb = bal.endcb_grc
    totcolc   = soilbgc_cstate.totc_col

    gpp                     = cnveg_cflux.gpp_col
    er                      = cnveg_cflux.er_col
    col_fire_closs          = cnveg_cflux.fire_closs_col
    col_hrv_xsmrpool_to_atm = cnveg_cflux.hrv_xsmrpool_to_atm_col
    col_xsmrpool_to_atm     = cnveg_cflux.xsmrpool_to_atm_col
    wood_harvestc           = cnveg_cflux.wood_harvestc_col
    gru_conv_cflux          = cnveg_cflux.gru_conv_cflux_col
    gru_wood_productc_gain  = cnveg_cflux.gru_wood_productc_gain_col
    crop_harvestc_to_cropprodc = cnveg_cflux.crop_harvestc_to_cropprodc_col

    som_c_leached    = soilbgc_cflux.som_c_leached_col
    hr_col           = soilbgc_cflux.hr_col
    fates_litter_flux = soilbgc_cflux.fates_litter_flux

    nbp_grc                   = cnveg_cflux.nbp_grc
    dwt_seedc_to_leaf_grc     = cnveg_cflux.dwt_seedc_to_leaf_grc
    dwt_seedc_to_deadstem_grc = cnveg_cflux.dwt_seedc_to_deadstem_grc

    cropprod1_grc    = c_products.cropprod1_grc
    tot_woodprod_grc = c_products.tot_woodprod_grc

    # Column-level allocation for error tracking
    nc = length(bounds_c) > 0 ? last(bounds_c) : 0
    col_errcb = zeros(nc)

    err_found = false
    err_index = 0

    # --- Column-level balance check ---
    for c in bounds_c
        mask_soil[c] || continue

        col_endcb[c] = totcolc[c]

        if !isempty(is_fates_col) && is_fates_col[c]
            col_cinputs = fates_litter_flux[c]
            col_coutputs = hr_col[c]
        else
            col_cinputs = gpp[c]
            col_coutputs = er[c] + col_fire_closs[c] + col_hrv_xsmrpool_to_atm[c] +
                           col_xsmrpool_to_atm[c] + gru_conv_cflux[c]
            col_coutputs += wood_harvestc[c] +
                            gru_wood_productc_gain[c] +
                            crop_harvestc_to_cropprodc[c]
        end

        # subtract leaching flux
        col_coutputs -= som_c_leached[c]

        # calculate the total column-level carbon balance error for this time step
        col_errcb[c] = (col_cinputs - col_coutputs) * dt -
                        (col_endcb[c] - col_begcb[c])

        if abs(col_errcb[c]) > bal.cerror
            err_found = true
            err_index = c
        end
        if abs(col_errcb[c]) > bal.cwarning
            @warn "cbalance warning at c = $c" col_errcb=col_errcb[c] col_endcb=col_endcb[c]
        end
    end

    if err_found
        c = err_index
        is_fates = !isempty(is_fates_col) && is_fates_col[c]
        g = col_data.gridcell[c]
        msg = string(
            "column cbalance error    = ", col_errcb[c], " c=", c, "\n",
            "is fates column?         = ", is_fates, "\n",
            "Latdeg,Londeg            = ", grc_data.latdeg[g], ",", grc_data.londeg[g], "\n",
            "begcb                    = ", col_begcb[c], "\n",
            "endcb                    = ", col_endcb[c], "\n",
            "delta store              = ", col_endcb[c] - col_begcb[c], "\n",
            "--- Inputs ---\n",
            is_fates ? "fates litter_flux        = $(fates_litter_flux[c]*dt)\n" :
                       "gpp                      = $(gpp[c]*dt)\n",
            "--- Outputs ---\n",
            !is_fates ? string(
                "er                       = ", er[c]*dt, "\n",
                "col_fire_closs           = ", col_fire_closs[c]*dt, "\n",
                "col_hrv_xsmrpool_to_atm  = ", col_hrv_xsmrpool_to_atm[c]*dt, "\n",
                "col_xsmrpool_to_atm      = ", col_xsmrpool_to_atm[c]*dt, "\n",
                "wood_harvestc            = ", wood_harvestc[c]*dt, "\n",
                "crop_harvestc_to_cropprodc = ", crop_harvestc_to_cropprodc[c]*dt, "\n"
            ) : string("hr                       = ", hr_col[c]*dt, "\n"),
            "-1*som_c_leached         = ", som_c_leached[c]*dt
        )
        error("CNBalanceCheck carbon error:\n$msg")
    end

    # --- Gridcell-level balance check ---

    # Column-to-gridcell aggregation
    totgrcc = soilbgc_cstate.totc_grc
    c2g_unity!(totgrcc, totcolc, col_data.gridcell, col_data.wtgcell, bounds_c, bounds_g)

    som_c_leached_grc = zeros(length(grc_data.latdeg))
    c2g_unity!(som_c_leached_grc, som_c_leached, col_data.gridcell, col_data.wtgcell, bounds_c, bounds_g)

    ng = length(grc_data.latdeg)
    grc_errcb = zeros(ng)

    err_found = false
    err_index = 0

    for g in bounds_g
        if !use_fates_bgc
            grc_endcb[g] = totgrcc[g] + tot_woodprod_grc[g] + cropprod1_grc[g]
            if !isempty(hrv_xsmrpool_amount_left)
                grc_endcb[g] += hrv_xsmrpool_amount_left[g] +
                                gru_conv_cflux_amount_left[g] +
                                dwt_conv_cflux_amount_left[g]
            end

            grc_cinputs = nbp_grc[g] +
                          dwt_seedc_to_leaf_grc[g] + dwt_seedc_to_deadstem_grc[g]
            grc_coutputs = -som_c_leached_grc[g]

            grc_errcb[g] = (grc_cinputs - grc_coutputs) * dt -
                            (grc_endcb[g] - grc_begcb[g])
        else
            grc_endcb[g] = grc_begcb[g]
            grc_errcb[g] = 0.0
        end

        if abs(grc_errcb[g]) > bal.cerror
            err_found = true
            err_index = g
        end
        if abs(grc_errcb[g]) > bal.cwarning
            @warn "cbalance warning at g = $g" grc_errcb=grc_errcb[g] grc_endcb=grc_endcb[g]
        end
    end

    if err_found
        g = err_index
        msg = string(
            "gridcell cbalance error  = ", grc_errcb[g], " g=", g, "\n",
            "latdeg, londeg           = ", grc_data.latdeg[g], ",", grc_data.londeg[g], "\n",
            "begcb                    = ", grc_begcb[g], "\n",
            "endcb                    = ", grc_endcb[g], "\n",
            "delta store              = ", grc_endcb[g] - grc_begcb[g], "\n",
            "--- Inputs ---\n",
            "nbp_grc                  = ", nbp_grc[g] * dt, "\n",
            "dwt_seedc_to_leaf_grc    = ", dwt_seedc_to_leaf_grc[g] * dt, "\n",
            "dwt_seedc_to_deadstem_grc = ", dwt_seedc_to_deadstem_grc[g] * dt, "\n",
            "--- Outputs ---\n",
            "-1*som_c_leached_grc     = ", som_c_leached_grc[g] * dt
        )
        error("CNBalanceCheck gridcell carbon error:\n$msg")
    end

    return nothing
end

# --------------------------------------------------------------------------
# NBalanceCheck
# --------------------------------------------------------------------------

"""
    n_balance_check!(bal, soilbgc_nflux, soilbgc_nstate, cnveg_nflux,
        n_products, col_data, grc_data,
        mask_soil, bounds_c, bounds_g, dt;
        is_fates_col, use_fates_bgc, use_nitrif_denitrif,
        use_crop, use_fun)

Perform nitrogen mass conservation check for column and gridcell.

Ported from `NBalanceCheck` in `CNBalanceCheckMod.F90`.

# Arguments
- `bal`: CNBalanceData with beginning-of-step masses already set
- `soilbgc_nflux`: SoilBiogeochemNitrogenFluxData
- `soilbgc_nstate`: SoilBiogeochemNitrogenStateData
- `cnveg_nflux`: CNVegNitrogenFluxData
- `n_products`: CNProductsData for nitrogen
- `col_data`: ColumnData for subgrid mapping
- `grc_data`: GridcellData for lat/lon
- `mask_soil`: BitVector mask for soil columns
- `bounds_c`: column index range
- `bounds_g`: gridcell index range
- `dt`: timestep size (seconds)
- `is_fates_col`: per-column flag for FATES columns
- `use_fates_bgc`: global FATES BGC flag
- `use_nitrif_denitrif`: flag for nitrification/denitrification
- `use_crop`: flag for crop model
- `use_fun`: flag for FUN model
"""
function n_balance_check!(
    bal::CNBalanceData,
    soilbgc_nflux::SoilBiogeochemNitrogenFluxData,
    soilbgc_nstate::SoilBiogeochemNitrogenStateData,
    cnveg_nflux::CNVegNitrogenFluxData,
    n_products::CNProductsData,
    col_data::ColumnData,
    grc_data::GridcellData,
    mask_soil::BitVector,
    bounds_c::UnitRange{Int},
    bounds_g::UnitRange{Int},
    dt::Float64;
    is_fates_col::Vector{Bool}=Bool[],
    use_fates_bgc::Bool=false,
    use_nitrif_denitrif::Bool=true,
    use_crop::Bool=false,
    use_fun::Bool=false
)
    # Local aliases
    col_begnb = bal.begnb_col
    col_endnb = bal.endnb_col
    grc_begnb = bal.begnb_grc
    grc_endnb = bal.endnb_grc

    totcoln           = soilbgc_nstate.totn_col
    ndep_to_sminn     = soilbgc_nflux.ndep_to_sminn_col
    nfix_to_sminn     = soilbgc_nflux.nfix_to_sminn_col
    ffix_to_sminn     = soilbgc_nflux.ffix_to_sminn_col
    fert_to_sminn     = soilbgc_nflux.fert_to_sminn_col
    soyfixn_to_sminn  = soilbgc_nflux.soyfixn_to_sminn_col
    supplement_to_sminn = soilbgc_nflux.supplement_to_sminn_col
    denit             = soilbgc_nflux.denit_col
    sminn_leached     = soilbgc_nflux.sminn_leached_col
    smin_no3_leached  = soilbgc_nflux.smin_no3_leached_col
    smin_no3_runoff   = soilbgc_nflux.smin_no3_runoff_col
    f_n2o_nit         = soilbgc_nflux.f_n2o_nit_col
    som_n_leached     = soilbgc_nflux.som_n_leached_col
    sminn_to_plant    = soilbgc_nflux.sminn_to_plant_col
    fates_litter_flux = soilbgc_nflux.fates_litter_flux

    col_fire_nloss    = cnveg_nflux.fire_nloss_col
    wood_harvestn     = cnveg_nflux.wood_harvestn_col
    gru_conv_nflux    = cnveg_nflux.gru_conv_nflux_col
    gru_wood_productn_gain = cnveg_nflux.gru_wood_productn_gain_col
    crop_harvestn_to_cropprodn = cnveg_nflux.crop_harvestn_to_cropprodn_col

    dwt_seedn_to_leaf_grc     = cnveg_nflux.dwt_seedn_to_leaf_grc
    dwt_seedn_to_deadstem_grc = cnveg_nflux.dwt_seedn_to_deadstem_grc
    dwt_conv_nflux_grc        = cnveg_nflux.dwt_conv_nflux_grc
    gru_conv_nflux_grc        = cnveg_nflux.gru_conv_nflux_grc
    gru_wood_productn_gain_grc = cnveg_nflux.gru_wood_productn_gain_grc

    cropprod1_grc    = n_products.cropprod1_grc
    tot_woodprod_grc = n_products.tot_woodprod_grc
    product_loss_grc = n_products.product_loss_grc

    # Column-level allocations
    nc = length(bounds_c) > 0 ? last(bounds_c) : 0
    col_ninputs         = zeros(nc)
    col_noutputs        = zeros(nc)
    col_errnb           = zeros(nc)
    col_ninputs_partial  = zeros(nc)
    col_noutputs_partial = zeros(nc)

    err_found = false
    err_index = 0

    for c in bounds_c
        mask_soil[c] || continue

        # calculate total column-level nitrogen storage
        col_endnb[c] = totcoln[c]

        # calculate total column-level inputs
        col_ninputs[c] = ndep_to_sminn[c] + nfix_to_sminn[c] + supplement_to_sminn[c]

        if !isempty(is_fates_col) && is_fates_col[c]
            col_ninputs[c] += fates_litter_flux[c]
        end

        if use_fun
            col_ninputs[c] += ffix_to_sminn[c]
        end

        if use_crop
            col_ninputs[c] += fert_to_sminn[c] + soyfixn_to_sminn[c]
        end

        col_ninputs_partial[c] = col_ninputs[c]

        # calculate total column-level outputs
        col_noutputs[c] = denit[c]

        if isempty(is_fates_col) || !is_fates_col[c]
            col_noutputs[c] += col_fire_nloss[c] + gru_conv_nflux[c]
            col_noutputs[c] += wood_harvestn[c] +
                               gru_wood_productn_gain[c] +
                               crop_harvestn_to_cropprodn[c]
        else
            col_noutputs[c] += sminn_to_plant[c]
        end

        if !use_nitrif_denitrif
            col_noutputs[c] += sminn_leached[c]
        else
            col_noutputs[c] += f_n2o_nit[c]
            col_noutputs[c] += smin_no3_leached[c] + smin_no3_runoff[c]
        end

        col_noutputs[c] -= som_n_leached[c]

        col_noutputs_partial[c] = col_noutputs[c]

        if isempty(is_fates_col) || !is_fates_col[c]
            col_noutputs_partial[c] -= wood_harvestn[c] -
                                       crop_harvestn_to_cropprodn[c]
        end

        # calculate column-level nitrogen balance error
        col_errnb[c] = (col_ninputs[c] - col_noutputs[c]) * dt -
                        (col_endnb[c] - col_begnb[c])

        if abs(col_errnb[c]) > bal.nerror
            err_found = true
            err_index = c
        end

        if abs(col_errnb[c]) > bal.nwarning
            @warn "nbalance warning at c = $c" col_errnb=col_errnb[c] col_endnb=col_endnb[c] inputs_ffix=ffix_to_sminn[c]*dt inputs_nfix=nfix_to_sminn[c]*dt inputs_ndep=ndep_to_sminn[c]*dt outputs_lch=smin_no3_leached[c]*dt outputs_roff=smin_no3_runoff[c]*dt outputs_dnit=f_n2o_nit[c]*dt
        end
    end

    if err_found
        c = err_index
        is_fates = !isempty(is_fates_col) && is_fates_col[c]
        g = col_data.gridcell[c]
        msg = string(
            "column nbalance error    = ", col_errnb[c], " c=", c, "\n",
            "Latdeg,Londeg            = ", grc_data.latdeg[g], ",", grc_data.londeg[g], "\n",
            "begnb                    = ", col_begnb[c], "\n",
            "endnb                    = ", col_endnb[c], "\n",
            "delta store              = ", col_endnb[c] - col_begnb[c], "\n",
            "input mass               = ", col_ninputs[c] * dt, "\n",
            "output mass              = ", col_noutputs[c] * dt, "\n",
            "net flux                 = ", (col_ninputs[c] - col_noutputs[c]) * dt, "\n",
            is_fates ?
                "inputs,ndep,nfix,suppn   = $(ndep_to_sminn[c]*dt),$(nfix_to_sminn[c]*dt),$(supplement_to_sminn[c]*dt)\n" :
                "inputs,ffix,nfix,ndep    = $(ffix_to_sminn[c]*dt),$(nfix_to_sminn[c]*dt),$(ndep_to_sminn[c]*dt)\n",
            is_fates ?
                "outputs,lch,roff,dnit,plnt = $(smin_no3_leached[c]*dt),$(smin_no3_runoff[c]*dt),$(f_n2o_nit[c]*dt),$(sminn_to_plant[c]*dt)\n" :
                "outputs,lch,roff,dnit    = $(smin_no3_leached[c]*dt),$(smin_no3_runoff[c]*dt),$(f_n2o_nit[c]*dt)\n"
        )
        error("CNBalanceCheck nitrogen error:\n$msg")
    end

    # --- Gridcell-level balance check (only when not using FATES BGC) ---
    if !use_fates_bgc
        totgrcn = soilbgc_nstate.totn_grc
        c2g_unity!(totgrcn, totcoln, col_data.gridcell, col_data.wtgcell, bounds_c, bounds_g)

        ng = length(grc_data.latdeg)
        grc_ninputs_partial  = zeros(ng)
        grc_noutputs_partial = zeros(ng)

        c2g_unity!(grc_ninputs_partial, col_ninputs_partial,
                   col_data.gridcell, col_data.wtgcell, bounds_c, bounds_g)
        c2g_unity!(grc_noutputs_partial, col_noutputs_partial,
                   col_data.gridcell, col_data.wtgcell, bounds_c, bounds_g)

        grc_errnb = zeros(ng)
        err_found = false
        err_index = 0

        for g in bounds_g
            grc_endnb[g] = totgrcn[g] + tot_woodprod_grc[g] + cropprod1_grc[g]

            grc_ninputs = grc_ninputs_partial[g] +
                          dwt_seedn_to_leaf_grc[g] +
                          dwt_seedn_to_deadstem_grc[g]

            grc_noutputs = grc_noutputs_partial[g] +
                           dwt_conv_nflux_grc[g] +
                           product_loss_grc[g] -
                           gru_wood_productn_gain_grc[g]

            grc_errnb[g] = (grc_ninputs - grc_noutputs) * dt -
                            (grc_endnb[g] - grc_begnb[g])

            if abs(grc_errnb[g]) > bal.nerror
                err_found = true
                err_index = g
            end
            if abs(grc_errnb[g]) > bal.nwarning
                @warn "nbalance warning at g = $g" grc_errnb=grc_errnb[g] grc_endnb=grc_endnb[g]
            end
        end

        if err_found
            g = err_index
            msg = string(
                "gridcell nbalance error  = ", grc_errnb[g], " g=", g, "\n",
                "latdeg, londeg           = ", grc_data.latdeg[g], ",", grc_data.londeg[g], "\n",
                "begnb                    = ", grc_begnb[g], "\n",
                "endnb                    = ", grc_endnb[g], "\n",
                "delta store              = ", grc_endnb[g] - grc_begnb[g], "\n",
                "input mass               = ", grc_ninputs_partial[g] * dt, "\n",
                "output mass              = ", grc_noutputs_partial[g] * dt, "\n",
                "net flux                 = ", (grc_ninputs_partial[g] - grc_noutputs_partial[g]) * dt, "\n",
                "--- Inputs ---\n",
                "grc_ninputs_partial      = ", grc_ninputs_partial[g] * dt, "\n",
                "dwt_seedn_to_leaf_grc    = ", dwt_seedn_to_leaf_grc[g] * dt, "\n",
                "dwt_seedn_to_deadstem_grc = ", dwt_seedn_to_deadstem_grc[g] * dt, "\n",
                "--- Outputs ---\n",
                "grc_noutputs_partial     = ", grc_noutputs_partial[g] * dt, "\n",
                "dwt_conv_nflux_grc       = ", dwt_conv_nflux_grc[g] * dt, "\n",
                "-gru_wood_productn_gain_grc = ", -gru_wood_productn_gain_grc[g] * dt, "\n",
                "product_loss_grc         = ", product_loss_grc[g] * dt
            )
            error("CNBalanceCheck gridcell nitrogen error:\n$msg")
        end
    end

    return nothing
end
