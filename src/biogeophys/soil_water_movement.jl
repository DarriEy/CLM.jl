# ==========================================================================
# Ported from: src/biogeophys/SoilWaterMovementMod.F90 (~2210 lines)
# Soil water movement — coupling soil and root water interactions.
#
# Public functions:
#   soil_water!                        — main entry point (dispatches to method)
#   init_soilwater_movement            — create config with defaults
#   use_aquifer_layer                  — check if aquifer layer is used
#   baseflow_sink!                     — vertically distributed baseflow sink
#
# Private helper functions:
#   soilwater_zengdecker2009!          — Zeng-Decker 2009 method
#   soilwater_moisture_form!           — moisture-based Richards eqn
#   compute_hydraulic_properties!      — hk, smp and derivatives
#   compute_moisture_fluxes_and_derivs!— flux at layer boundaries
#   compute_RHS_moisture_form!         — RHS of moisture-based Richards eqn
#   compute_LHS_moisture_form!         — LHS tridiagonal matrix
#   compute_qcharge!                   — aquifer recharge
#   ice_impedance                      — ice impedance factor
#   tridiagonal_col!                   — single-column tridiagonal solver
# ==========================================================================

# ---- Solution method constants ----
const ZENGDECKER_2009 = 0
const MOISTURE_FORM   = 1
const MIXED_FORM      = 2
const HEAD_FORM       = 3

# ---- Boundary condition constants ----
const BC_HEAD        = 0
const BC_FLUX        = 1
const BC_ZERO_FLUX   = 2
const BC_WATERTABLE  = 3
const BC_AQUIFER     = 4

# ---- Unit conversion ----
const M_TO_MM = 1.0e3

"""
    SoilWaterMovementConfig

Configuration for the soil water movement module.
Holds method selection, boundary conditions, adaptive time stepping
parameters, and the ice impedance parameter e_ice.

Ported from module-level variables in `SoilWaterMovementMod.F90`.
"""
Base.@kwdef mutable struct SoilWaterMovementConfig
    # Solution method
    soilwater_movement_method::Int = ZENGDECKER_2009
    upper_boundary_condition::Int  = BC_FLUX
    lower_boundary_condition::Int  = BC_AQUIFER

    # Adaptive time stepping parameters
    dtmin::Float64         = 60.0       # minimum time step length (seconds)
    verySmall::Float64     = 1.0e-8     # check for sub step completion
    xTolerUpper::Float64   = 1.0e-1     # tolerance to halve length of substep
    xTolerLower::Float64   = 1.0e-2     # tolerance to double length of substep
    expensive::Int         = 42
    inexpensive::Int       = 1
    flux_calculation::Int  = 1          # default = inexpensive

    # Soil ice impedance factor (unitless) — from params file
    e_ice::Float64 = 6.0
end

"""
    init_soilwater_movement(; kwargs...) -> SoilWaterMovementConfig

Create and return a `SoilWaterMovementConfig` with the given keyword
arguments. Validates consistency of method and boundary conditions.

Ported from `init_soilwater_movement` in `SoilWaterMovementMod.F90`.
"""
function init_soilwater_movement(; kwargs...)
    cfg = SoilWaterMovementConfig(; kwargs...)

    # Validate consistency
    if cfg.soilwater_movement_method == ZENGDECKER_2009 &&
       cfg.lower_boundary_condition != BC_AQUIFER
        error("init_soilwater_movement: ZD09 must use bc_aquifer lower boundary condition")
    end

    return cfg
end

"""
    use_aquifer_layer(cfg::SoilWaterMovementConfig) -> Bool

Return true if an aquifer layer is used (lower boundary is aquifer or water table).

Ported from `use_aquifer_layer` in `SoilWaterMovementMod.F90`.
"""
function use_aquifer_layer(cfg::SoilWaterMovementConfig)
    return cfg.lower_boundary_condition == BC_AQUIFER ||
           cfg.lower_boundary_condition == BC_WATERTABLE
end

# ===========================================================================
# ice_impedance — compute hydraulic conductivity reduction due to ice
# ===========================================================================
"""
    ice_impedance(icefrac, e_ice) -> Float64

Compute the hydraulic conductivity reduction factor due to ice in pore space.

Ported from `IceImpedance` in `SoilWaterMovementMod.F90`.
"""
function ice_impedance(icefrac::Float64, e_ice::Float64)
    return 10.0^(-e_ice * icefrac)
end

# ===========================================================================
# tridiagonal_col! — single-column tridiagonal solver
# ===========================================================================
"""
    tridiagonal_col!(u, a, b, c, r, jtop, lbj, ubj)

Tridiagonal matrix solution for a single column.

Ported from `TridiagonalCol` in `SoilWaterMovementMod.F90`.
"""
function tridiagonal_col!(u::AbstractVector{Float64},
                          a::AbstractVector{Float64},
                          b::AbstractVector{Float64},
                          c::AbstractVector{Float64},
                          r::AbstractVector{Float64},
                          jtop::Int, lbj::Int, ubj::Int)
    n = ubj - lbj + 1
    gam = zeros(n)

    bet = b[jtop]

    for j in lbj:ubj
        if j >= jtop
            if j == jtop
                u[j] = r[j] / bet
            else
                gam[j] = c[j-1] / bet
                bet = b[j] - a[j] * gam[j]
                u[j] = (r[j] - a[j] * u[j-1]) / bet
            end
        end
    end

    for j in (ubj-1):-1:lbj
        if j >= jtop
            u[j] = u[j] - gam[j+1] * u[j+1]
        end
    end

    nothing
end

# ===========================================================================
# baseflow_sink! — vertically distributed baseflow (placeholder)
# ===========================================================================
"""
    baseflow_sink!(baseflow_sink, mask_hydrology, nlevsoi)

Apply baseflow as a vertically distributed sink condition.
Currently a placeholder that sets baseflow_sink to zero.

Ported from `BaseflowSink` in `SoilWaterMovementMod.F90`.
"""
function baseflow_sink!(baseflow_sink::Matrix{Float64},
                        mask_hydrology::BitVector,
                        nlevsoi::Int)
    # Placeholder — just zero out
    fill!(baseflow_sink, 0.0)
    nothing
end

# ===========================================================================
# compute_hydraulic_properties! — hk, smp and derivatives for moisture form
# ===========================================================================
"""
    compute_hydraulic_properties!(c, nlayers, soilhydrology, soilstate,
        swrc, cfg, vwc_liq, hk, smp, dhkdw, dsmpdw, imped_out)

Calculate soil hydraulic conductivity, matric potential and their derivatives.

Ported from `compute_hydraulic_properties` in `SoilWaterMovementMod.F90`.
"""
function compute_hydraulic_properties!(c::Int, nlayers::Int,
        soilhydrology::SoilHydrologyData, soilstate::SoilStateData,
        swrc::SoilWaterRetentionCurve, cfg::SoilWaterMovementConfig,
        vwc_liq::AbstractVector{Float64},
        hk::AbstractVector{Float64}, smp_out::AbstractVector{Float64},
        dhkdw::AbstractVector{Float64}, dsmpdw::AbstractVector{Float64},
        imped_out::AbstractVector{Float64})

    icefrac = soilhydrology.icefrac_col
    watsat  = soilstate.watsat_col
    smp_l   = soilstate.smp_l_col
    hk_l    = soilstate.hk_l_col

    fill!(hk, 0.0)
    fill!(smp_out, 0.0)
    fill!(dhkdw, 0.0)
    fill!(dsmpdw, 0.0)
    fill!(imped_out, 0.0)

    # Compute relative saturation at each layer node
    s2 = zeros(nlayers)
    for j in 1:nlayers
        s2[j] = vwc_liq[j] / watsat[c, j]
        s2[j] = min(s2[j], 1.0)
        s2[j] = max(0.01, s2[j])
    end

    for j in 1:nlayers
        # s1 is interface value, s2 is node value
        if j == nlayers
            s1 = s2[j]
            imped_out[j] = ice_impedance(icefrac[c, j], cfg.e_ice)
        else
            s1 = 0.5 * (s2[j] + s2[j+1])
            imped_out[j] = ice_impedance(0.5 * (icefrac[c, j] + icefrac[c, j+1]), cfg.e_ice)
        end

        # Impose constraints on relative saturation at interface
        s1 = min(s1, 1.0)
        s1 = max(0.01, s1)

        # Compute hydraulic conductivity
        hk_val, dhkds = soil_hk!(swrc, c, j, s1, imped_out[j], soilstate)
        hk[j] = hk_val

        # Compute matric potential
        smp_val, dsmpds = soil_suction!(swrc, c, j, s2[j], soilstate)
        smp_out[j] = smp_val

        # Save derivative of hk w.r.t. relative saturation at interface
        dhkdw[j] = dhkds

        # Compute derivative w.r.t. volumetric liquid water content
        dsmpdw[j] = dsmpds / watsat[c, j]

        # Save for output
        smp_l[c, j] = smp_out[j]
        hk_l[c, j] = hk[j]
    end

    nothing
end

# ===========================================================================
# compute_moisture_fluxes_and_derivs! — fluxes at layer boundaries
# ===========================================================================
"""
    compute_moisture_fluxes_and_derivs!(c, nlayers, soilhydrology, soilstate,
        temperature, waterfluxbulk, swrc, cfg, vwc_liq, hk, smp, dhkdw,
        dsmpdw, imped, qin, qout, dqidw0, dqidw1, dqodw1, dqodw2)

Calculate fluxes at the boundary of each layer and derivatives.

Ported from `compute_moisture_fluxes_and_derivs` in `SoilWaterMovementMod.F90`.
"""
function compute_moisture_fluxes_and_derivs!(c::Int, nlayers::Int,
        soilhydrology::SoilHydrologyData, soilstate::SoilStateData,
        temperature::TemperatureData, waterfluxbulk::WaterFluxBulkData,
        swrc::SoilWaterRetentionCurve, cfg::SoilWaterMovementConfig,
        vwc_liq::AbstractVector{Float64},
        hk::AbstractVector{Float64}, smp::AbstractVector{Float64},
        dhkdw::AbstractVector{Float64}, dsmpdw::AbstractVector{Float64},
        imped::AbstractVector{Float64},
        qin::AbstractVector{Float64}, qout::AbstractVector{Float64},
        dqidw0::AbstractVector{Float64}, dqidw1::AbstractVector{Float64},
        dqodw1::AbstractVector{Float64}, dqodw2::AbstractVector{Float64})

    col_z      = soilhydrology.zwt_col  # use just for dispatch; actual z from col
    nlevsoi    = varpar.nlevsoi
    watsat     = soilstate.watsat_col
    qflx_infl  = waterfluxbulk.wf.qflx_infl_col
    zwt        = soilhydrology.zwt_col

    fill!(qin, 0.0)
    fill!(qout, 0.0)
    fill!(dqidw0, 0.0)
    fill!(dqidw1, 0.0)
    fill!(dqodw1, 0.0)
    fill!(dqodw2, 0.0)

    # We need column z, zi, dz — these come from a ColumnData struct
    # But in the moisture form, z/zi/dz are accessed via col%z etc.
    # For this port, we pass z/zi/dz separately or access from a ColumnData
    # Here we assume the caller provides col_data with z, zi, dz fields.
    # This function is called from soilwater_moisture_form! which has access to col.

    # NOTE: This function is an internal helper called from soilwater_moisture_form!
    # The actual col z/zi/dz access happens through the col argument passed by the caller.
    # For this port, we'll accept z_col, zi_col as additional parameters.
    error("compute_moisture_fluxes_and_derivs! requires col data — use the full-argument version")
end

"""
    compute_moisture_fluxes_and_derivs!(c, nlayers, z_col, zi_col, dz_col,
        soilhydrology, soilstate, waterfluxbulk, swrc, cfg, vwc_liq,
        hk, smp, dhkdw, dsmpdw, imped,
        qin, qout, dqidw0, dqidw1, dqodw1, dqodw2)

Full-argument version with explicit column geometry arrays.
"""
function compute_moisture_fluxes_and_derivs!(c::Int, nlayers::Int,
        z_col::Matrix{Float64}, zi_col::Matrix{Float64}, dz_col::Matrix{Float64},
        soilhydrology::SoilHydrologyData, soilstate::SoilStateData,
        waterfluxbulk::WaterFluxBulkData,
        swrc::SoilWaterRetentionCurve, cfg::SoilWaterMovementConfig,
        vwc_liq::AbstractVector{Float64},
        hk::AbstractVector{Float64}, smp::AbstractVector{Float64},
        dhkdw::AbstractVector{Float64}, dsmpdw::AbstractVector{Float64},
        imped::AbstractVector{Float64},
        qin::AbstractVector{Float64}, qout::AbstractVector{Float64},
        dqidw0::AbstractVector{Float64}, dqidw1::AbstractVector{Float64},
        dqodw1::AbstractVector{Float64}, dqodw2::AbstractVector{Float64})

    nlevsoi    = varpar.nlevsoi
    watsat     = soilstate.watsat_col
    qflx_infl  = waterfluxbulk.wf.qflx_infl_col
    zwt        = soilhydrology.zwt_col

    fill!(qin, 0.0)
    fill!(qout, 0.0)
    fill!(dqidw0, 0.0)
    fill!(dqidw1, 0.0)
    fill!(dqodw1, 0.0)
    fill!(dqodw2, 0.0)

    zdflag = 0  # matching Fortran default

    # ---- Top layer (j=1) ----
    j = 1

    # Upper boundary condition
    if cfg.upper_boundary_condition == BC_FLUX
        # Specify infiltration
        qin[j] = qflx_infl[c]
        dqidw1[j] = 0.0

    elseif cfg.upper_boundary_condition == BC_HEAD
        # Head boundary condition (not commonly used)
        error("compute_moisture_fluxes_and_derivs!: BC_HEAD upper boundary not fully implemented")

    else
        error("compute_moisture_fluxes_and_derivs!: upper boundary condition must be specified!")
    end

    # Compute derivatives in hk at interface w.r.t. vwc in layers above and below
    dhkds1 = 0.5 * dhkdw[j] / watsat[c, j]
    dhkds2 = 0.5 * dhkdw[j] / watsat[c, j+1]

    if zdflag == 1
        dhkds1 = dhkdw[j] / (watsat[c, j] + watsat[c, min(nlevsoi, j+1)])
        dhkds2 = dhkds1
    end

    # Compute flux at bottom of j-th layer
    num = smp[j+1] - smp[j]
    den = M_TO_MM * (z_col[c, j+1] - z_col[c, j])
    qout[j] = -hk[j] * num / den + hk[j]

    # Compute flux derivatives
    dqodw1[j] = (hk[j] * dsmpdw[j] - dhkds1 * num) / den + dhkds1
    dqodw2[j] = (-hk[j] * dsmpdw[j+1] - dhkds2 * num) / den + dhkds2

    # ---- Interior nodes (j=2 to nlayers-1) ----
    for j in 2:(nlayers-1)
        qin[j] = qout[j-1]
        dqidw0[j] = dqodw1[j-1]
        dqidw1[j] = dqodw2[j-1]

        dhkds1 = 0.5 * dhkdw[j] / watsat[c, j]
        dhkds2 = 0.5 * dhkdw[j] / watsat[c, j+1]

        if zdflag == 1
            dhkds1 = dhkdw[j] / (watsat[c, j] + watsat[c, min(nlevsoi, j+1)])
            dhkds2 = dhkds1
        end

        num = smp[j+1] - smp[j]
        den = M_TO_MM * (z_col[c, j+1] - z_col[c, j])
        qout[j] = -hk[j] * num / den + hk[j]

        dqodw1[j] = (hk[j] * dsmpdw[j] - dhkds1 * num) / den + dhkds1
        dqodw2[j] = (-hk[j] * dsmpdw[j+1] - dhkds2 * num) / den + dhkds2
    end

    # ---- Bottom node (j=nlayers) ----
    j = nlayers

    qin[j] = qout[j-1]
    dqidw0[j] = dqodw1[j-1]
    dqidw1[j] = dqodw2[j-1]

    if cfg.lower_boundary_condition == BC_WATERTABLE
        jwt = nlevsoi
        for jj in 1:nlevsoi
            if zwt[c] <= zi_col[c, jj]
                jwt = jj - 1
                break
            end
        end

        if j > jwt
            # Water table within soil column — no drainage
            qout[j] = 0.0
            dqodw1[j] = 0.0
            dqodw2[j] = 0.0
        else
            # Water table below soil column
            dhkds1_local = dhkdw[j] / watsat[c, j]

            if zdflag == 1
                dhkds1_local = dhkdw[j] / (watsat[c, j] + watsat[c, min(nlevsoi, j+1)])
            end

            num = -smp[j]  # assume saturation at water table depth (smp=0)
            den = M_TO_MM * (zwt[c] - z_col[c, j])
            qout[j] = -hk[j] * num / den + hk[j]

            dqodw1[j] = (hk[j] * dsmpdw[j] - dhkds1_local * num) / den + dhkds1_local
            dqodw2[j] = 0.0
        end

    elseif cfg.lower_boundary_condition == BC_FLUX
        # Free drainage
        qout[j] = hk[j]
        dqodw1[j] = dhkdw[j] / watsat[c, j]

    elseif cfg.lower_boundary_condition == BC_ZERO_FLUX
        qout[j] = 0.0
        dqodw1[j] = 0.0

    elseif cfg.lower_boundary_condition == BC_HEAD
        error("compute_moisture_fluxes_and_derivs!: BC_HEAD lower boundary not fully implemented/tested")

    else
        error("compute_moisture_fluxes_and_derivs!: lower boundary condition must be specified!")
    end

    nothing
end

# ===========================================================================
# compute_RHS_moisture_form! — RHS of moisture-based Richards equation
# ===========================================================================
"""
    compute_RHS_moisture_form!(c, nlayers, vert_trans_sink, vwc_liq,
        qin, qout, dt_dz, rmx)

Calculate RHS of moisture-based form of Richards equation.

Ported from `compute_RHS_moisture_form` in `SoilWaterMovementMod.F90`.
"""
function compute_RHS_moisture_form!(c::Int, nlayers::Int,
        vert_trans_sink::AbstractVector{Float64},
        vwc_liq::AbstractVector{Float64},
        qin::AbstractVector{Float64}, qout::AbstractVector{Float64},
        dt_dz::AbstractVector{Float64},
        rmx::AbstractVector{Float64})

    fill!(rmx, 0.0)

    for j in 1:nlayers
        fluxNet = qin[j] - qout[j] - vert_trans_sink[j]
        # Non-iterative solution
        rmx[j] = -fluxNet * dt_dz[j]
    end

    nothing
end

# ===========================================================================
# compute_LHS_moisture_form! — LHS tridiagonal matrix
# ===========================================================================
"""
    compute_LHS_moisture_form!(c, nlayers, dt_dz, dqidw0, dqidw1, dqodw1,
        dqodw2, amx, bmx, cmx)

Calculate LHS of moisture-based form of Richards equation (tridiagonal matrix).

Ported from `compute_LHS_moisture_form` in `SoilWaterMovementMod.F90`.
"""
function compute_LHS_moisture_form!(c::Int, nlayers::Int,
        dt_dz::AbstractVector{Float64},
        dqidw0::AbstractVector{Float64}, dqidw1::AbstractVector{Float64},
        dqodw1::AbstractVector{Float64}, dqodw2::AbstractVector{Float64},
        amx::AbstractVector{Float64}, bmx::AbstractVector{Float64},
        cmx::AbstractVector{Float64})

    fill!(amx, 0.0)
    fill!(bmx, 0.0)
    fill!(cmx, 0.0)

    # Top soil layer
    j = 1
    amx[j] = 0.0
    bmx[j] = -1.0 - (-dqidw1[j] + dqodw1[j]) * dt_dz[j]
    cmx[j] = -dqodw2[j] * dt_dz[j]

    # Interior soil layers
    for j in 2:(nlayers-1)
        amx[j] = dqidw0[j] * dt_dz[j]
        bmx[j] = -1.0 - (-dqidw1[j] + dqodw1[j]) * dt_dz[j]
        cmx[j] = -dqodw2[j] * dt_dz[j]
    end

    # Bottom soil layer
    j = nlayers
    amx[j] = dqidw0[j] * dt_dz[j]
    bmx[j] = -1.0 - (-dqidw1[j] + dqodw1[j]) * dt_dz[j]
    cmx[j] = 0.0

    nothing
end

# ===========================================================================
# compute_qcharge! — aquifer recharge for moisture form
# ===========================================================================
"""
    compute_qcharge!(col_data, mask_hydrology, soilhydrology, soilstate,
        waterstatebulk, swrc, cfg, dwat, smp, imped, vwc_liq, dtime, bounds)

Calculate additional aquifer recharge terms for the moisture form method.

Ported from `compute_qcharge` in `SoilWaterMovementMod.F90`.
"""
function compute_qcharge!(col_data::ColumnData,
        mask_hydrology::BitVector,
        soilhydrology::SoilHydrologyData, soilstate::SoilStateData,
        waterstatebulk::WaterStateBulkData,
        swrc::SoilWaterRetentionCurve, cfg::SoilWaterMovementConfig,
        dwat::Matrix{Float64}, smp::Matrix{Float64},
        imped::Matrix{Float64}, vwc_liq::Matrix{Float64},
        dtime::Float64, bounds::UnitRange{Int})

    nlevsoi = varpar.nlevsoi
    qcharge = soilhydrology.qcharge_col
    zwt     = soilhydrology.zwt_col
    sucsat  = soilstate.sucsat_col
    watsat  = soilstate.watsat_col
    smpmin  = soilstate.smpmin_col
    z       = col_data.z
    zi      = col_data.zi

    for c in eachindex(mask_hydrology)
        mask_hydrology[c] || continue

        # Locate index of layer above water table
        jwt = nlevsoi
        for j in 1:nlevsoi
            if zwt[c] <= zi[c, j]
                jwt = j - 1
                break
            end
        end

        # Calculate qcharge for case jwt < nlevsoi
        if jwt < nlevsoi
            wh_zwt = -sucsat[c, jwt+1] - zwt[c] * M_TO_MM

            # Recharge rate to groundwater (positive to aquifer)
            s1 = max(vwc_liq[c, jwt+1] / watsat[c, jwt+1], 0.01)
            s1 = min(1.0, s1)

            # Unsaturated hydraulic conductivity
            ka, _ = soil_hk!(swrc, c, jwt+1, s1, imped[c, jwt+1], soilstate)

            smp1 = max(smpmin[c], smp[c, max(1, jwt)])
            wh = smp1 - z[c, max(1, jwt)] * M_TO_MM

            if jwt == 0
                qcharge[c] = -ka * (wh_zwt - wh) / ((zwt[c] + 1.0e-3) * M_TO_MM)
            else
                qcharge[c] = -ka * (wh_zwt - wh) / ((zwt[c] - z[c, jwt]) * M_TO_MM * 2.0)
            end

            # Limit qcharge
            qcharge[c] = max(-10.0 / dtime, qcharge[c])
            qcharge[c] = min(10.0 / dtime, qcharge[c])
        end
    end

    nothing
end

# ===========================================================================
# soilwater_moisture_form! — moisture-based Richards equation solver
# ===========================================================================
"""
    soilwater_moisture_form!(col_data, mask_hydrology, mask_urban,
        soilhydrology, soilstate, waterfluxbulk, waterstatebulk,
        temperature, canopystate, energyflux, swrc, cfg, dtime)

Solve Richards equation using the moisture-based form with adaptive
time stepping.

Ported from `soilwater_moisture_form` in `SoilWaterMovementMod.F90`.
"""
function soilwater_moisture_form!(col_data::ColumnData,
        mask_hydrology::BitVector, mask_urban::BitVector,
        soilhydrology::SoilHydrologyData, soilstate::SoilStateData,
        waterfluxbulk::WaterFluxBulkData, waterstatebulk::WaterStateBulkData,
        temperature::TemperatureData, canopystate::CanopyStateData,
        energyflux::EnergyFluxData,
        swrc::SoilWaterRetentionCurve, cfg::SoilWaterMovementConfig,
        dtime::Float64)

    nlevsoi = varpar.nlevsoi
    nbedrock = col_data.nbedrock
    z  = col_data.z
    zi = col_data.zi
    dz = col_data.dz

    nsubsteps_col = soilhydrology.num_substeps_col
    qcharge       = soilhydrology.qcharge_col
    h2osoi_liq    = waterstatebulk.ws.h2osoi_liq_col
    qflx_rootsoi  = waterfluxbulk.qflx_rootsoi_col

    for c in eachindex(mask_hydrology)
        mask_hydrology[c] || continue

        nlayers = nbedrock[c]

        # Local arrays
        hk_loc     = zeros(nlayers)
        smp_loc    = zeros(nlayers)
        dhkdw_loc  = zeros(nlayers)
        dsmpdw_loc = zeros(nlayers)
        imped_loc  = zeros(nlayers)
        vwc_liq_loc = zeros(nlayers)
        dt_dz_loc  = zeros(nlayers)
        qin_loc    = zeros(nlayers)
        qout_loc   = zeros(nlayers)
        dqidw0_loc = zeros(nlayers)
        dqidw1_loc = zeros(nlayers)
        dqodw1_loc = zeros(nlayers)
        dqodw2_loc = zeros(nlayers)
        dwat_loc   = zeros(nlayers)
        amx_loc    = zeros(nlayers)
        bmx_loc    = zeros(nlayers)
        cmx_loc    = zeros(nlayers)
        rmx_loc    = zeros(nlayers)

        # Initialize adaptive substeps
        nsubstep = 0
        dtsub = dtime
        dtdone = 0.0
        qcharge[c] = 0.0

        # Adaptive sub-step loop
        while true
            nsubstep += 1

            # Calculate commonly used variables
            for j in 1:nlayers
                vwc_liq_loc[j] = max(h2osoi_liq[c, j], 1.0e-6) / (dz[c, j] * DENH2O)
                dt_dz_loc[j] = dtsub / (M_TO_MM * dz[c, j])
            end

            # Hydraulic conductivity and soil matric potential
            compute_hydraulic_properties!(c, nlayers,
                soilhydrology, soilstate, swrc, cfg,
                vwc_liq_loc, hk_loc, smp_loc, dhkdw_loc, dsmpdw_loc, imped_loc)

            # Soil moisture fluxes and derivatives
            compute_moisture_fluxes_and_derivs!(c, nlayers, z, zi, dz,
                soilhydrology, soilstate, waterfluxbulk, swrc, cfg,
                vwc_liq_loc, hk_loc, smp_loc, dhkdw_loc, dsmpdw_loc, imped_loc,
                qin_loc, qout_loc, dqidw0_loc, dqidw1_loc, dqodw1_loc, dqodw2_loc)

            # RHS
            rootsoi_slice = view(qflx_rootsoi, c, 1:nlayers)
            compute_RHS_moisture_form!(c, nlayers,
                rootsoi_slice,
                vwc_liq_loc, qin_loc, qout_loc, dt_dz_loc, rmx_loc)

            # LHS
            compute_LHS_moisture_form!(c, nlayers,
                dt_dz_loc, dqidw0_loc, dqidw1_loc, dqodw1_loc, dqodw2_loc,
                amx_loc, bmx_loc, cmx_loc)

            # Solve tridiagonal system
            tridiagonal_col!(dwat_loc, amx_loc, bmx_loc, cmx_loc, rmx_loc,
                1, 1, nlayers)

            # Error estimation
            fluxNet0 = zeros(nlayers)
            fluxNet1 = zeros(nlayers)

            for j in 1:nlayers
                if cfg.flux_calculation == cfg.expensive
                    # Expensive option
                    if j == 1
                        qin_test = qin_loc[j] + dqidw1_loc[j] * dwat_loc[j]
                    else
                        qin_test = qin_loc[j] + dqidw0_loc[j] * dwat_loc[j-1] + dqidw1_loc[j] * dwat_loc[j]
                    end

                    if j == nlayers
                        qout_test = qout_loc[j] + dqodw1_loc[j] * dwat_loc[j]
                    else
                        qout_test = qout_loc[j] + dqodw1_loc[j] * dwat_loc[j] + dqodw2_loc[j] * dwat_loc[j+1]
                    end

                    fluxNet0[j] = qin_test - qout_test - qflx_rootsoi[c, j]
                else
                    # Inexpensive: convert iteration increment to net flux
                    fluxNet0[j] = dwat_loc[j] / dt_dz_loc[j]
                end

                fluxNet1[j] = qin_loc[j] - qout_loc[j] - qflx_rootsoi[c, j]
            end

            # Compute absolute errors
            errorMax = 0.0
            for j in 1:nlayers
                err_val = abs(fluxNet1[j] - fluxNet0[j]) * dtsub * 0.5
                errorMax = max(errorMax, err_val)
            end

            # Check if error is above upper tolerance
            if errorMax > cfg.xTolerUpper && dtsub > cfg.dtmin
                dtsub = max(dtsub / 2.0, cfg.dtmin)
                continue  # substep rejected; try again
            end

            # Renew the mass of liquid water
            for j in 1:nlayers
                h2osoi_liq[c, j] = h2osoi_liq[c, j] + dwat_loc[j] * (M_TO_MM * dz[c, j])
            end

            # Compute drainage from bottom of soil column
            if cfg.lower_boundary_condition == BC_FLUX
                qcTemp = hk_loc[nlayers] + dhkdw_loc[nlayers] * dwat_loc[nlayers]
            elseif cfg.lower_boundary_condition == BC_ZERO_FLUX
                qcTemp = 0.0
            elseif cfg.lower_boundary_condition == BC_WATERTABLE
                qcTemp = qout_loc[nlayers] + dqodw1_loc[nlayers] * dwat_loc[nlayers]
            else
                error("soilwater_moisture_form!: lower boundary condition must be specified!")
            end

            # Increment qcharge flux
            qcharge[c] = qcharge[c] + qcTemp * (dtsub / dtime)

            # Increment substep and check for completion
            dtdone = dtdone + dtsub
            if abs(dtime - dtdone) < cfg.verySmall
                break  # time step completed
            end

            # Check if error is below lower tolerance
            if errorMax < cfg.xTolerLower
                dtsub = dtsub * 2.0
            end

            # Ensure substep does not exceed time remaining
            dtsub = min(dtsub, dtime - dtdone)
        end  # substep loop

        # Save number of adaptive substeps
        nsubsteps_col[c] = Float64(nsubstep)
    end  # spatial loop

    # Calculate qcharge when water table is in soil column (bc_watertable)
    if use_aquifer_layer(cfg)
        # Build temporary 2D arrays for compute_qcharge!
        nc = length(mask_hydrology)
        dwat_2d  = zeros(nc, nlevsoi)
        smp_2d   = zeros(nc, nlevsoi)
        imped_2d = zeros(nc, nlevsoi)
        vwc_2d   = zeros(nc, nlevsoi)

        for c in eachindex(mask_hydrology)
            mask_hydrology[c] || continue
            nlayers = nbedrock[c]
            for j in 1:min(nlayers, nlevsoi)
                vwc_2d[c, j] = max(h2osoi_liq[c, j], 1.0e-6) / (dz[c, j] * DENH2O)
            end
        end

        compute_qcharge!(col_data, mask_hydrology, soilhydrology, soilstate,
            waterstatebulk, swrc, cfg, dwat_2d, smp_2d, imped_2d, vwc_2d,
            dtime, 1:nc)
    end

    nothing
end

# ===========================================================================
# soilwater_zengdecker2009! — Zeng-Decker 2009 method
# ===========================================================================
"""
    soilwater_zengdecker2009!(col_data, mask_hydrology, mask_urban,
        soilhydrology, soilstate, waterfluxbulk, waterstatebulk,
        temperature, canopystate, energyflux, swrc, cfg, dtime)

Soil hydrology using the Zeng & Decker (2009) method with equilibrium
matric potential and coupled aquifer layer.

Ported from `soilwater_zengdecker2009` in `SoilWaterMovementMod.F90`.
"""
function soilwater_zengdecker2009!(col_data::ColumnData,
        mask_hydrology::BitVector, mask_urban::BitVector,
        soilhydrology::SoilHydrologyData, soilstate::SoilStateData,
        waterfluxbulk::WaterFluxBulkData, waterstatebulk::WaterStateBulkData,
        temperature::TemperatureData, canopystate::CanopyStateData,
        energyflux::EnergyFluxData,
        swrc::SoilWaterRetentionCurve, cfg::SoilWaterMovementConfig,
        dtime::Float64)

    nlevsoi  = varpar.nlevsoi
    nlevgrnd = varpar.nlevgrnd
    nlevsno  = varpar.nlevsno

    # Offsets for snow+soil arrays (z, dz, h2osoi_liq/ice, t_soisno)
    # Soil layer j in Fortran maps to index j + joff in Julia 1-based arrays
    joff    = nlevsno          # for z, dz, h2osoi_liq/ice, t_soisno
    joff_zi = nlevsno + 1     # for zi (has extra element at top)

    z  = col_data.z
    zi = col_data.zi
    dz = col_data.dz

    qcharge       = soilhydrology.qcharge_col
    zwt           = soilhydrology.zwt_col
    icefrac       = soilhydrology.icefrac_col
    hkdepth       = soilhydrology.hkdepth_col

    smpmin        = soilstate.smpmin_col
    watsat        = soilstate.watsat_col
    hksat         = soilstate.hksat_col
    bsw           = soilstate.bsw_col
    sucsat        = soilstate.sucsat_col
    eff_porosity  = soilstate.eff_porosity_col
    smp_l         = soilstate.smp_l_col
    hk_l          = soilstate.hk_l_col

    h2osoi_ice    = waterstatebulk.ws.h2osoi_ice_col
    h2osoi_liq    = waterstatebulk.ws.h2osoi_liq_col
    h2osoi_vol    = waterstatebulk.ws.h2osoi_vol_col

    qflx_deficit  = waterfluxbulk.qflx_deficit_col
    qflx_infl     = waterfluxbulk.wf.qflx_infl_col
    qflx_rootsoi  = waterfluxbulk.qflx_rootsoi_col

    t_soisno      = temperature.t_soisno_col

    nc = length(mask_hydrology)

    # Local arrays
    hk      = zeros(nc, nlevsoi)
    dhkdw   = zeros(nc, nlevsoi)
    smp_arr = zeros(nc, nlevsoi)
    dsmpdw  = zeros(nc, nlevsoi)
    imped_arr = zeros(nc, nlevsoi)
    vol_ice = zeros(nc, nlevsoi)
    vwc_liq = zeros(nc, nlevsoi + 1)
    zmm     = zeros(nc, nlevsoi + 1)
    dzmm    = zeros(nc, nlevsoi + 1)
    zimm    = zeros(nc, 0:nlevsoi)  # Can't do 0-indexed; use offset
    # Use 1-indexed with offset: zimm_arr[c, j+1] for Fortran zimm(c,j) where j=0:nlevsoi
    zimm_arr = zeros(nc, nlevsoi + 1)  # index 1 = Fortran j=0, index nlevsoi+1 = Fortran j=nlevsoi
    zwtmm   = zeros(nc)
    jwt     = zeros(Int, nc)
    vwc_zwt = zeros(nc)
    vol_eq  = zeros(nc, nlevsoi + 1)
    zq      = zeros(nc, nlevsoi + 1)
    qin     = zeros(nc, nlevsoi + 1)
    qout_arr = zeros(nc, nlevsoi + 1)
    dqidw0  = zeros(nc, nlevsoi + 1)
    dqidw1  = zeros(nc, nlevsoi + 1)
    dqodw1  = zeros(nc, nlevsoi + 1)
    dqodw2  = zeros(nc, nlevsoi + 1)
    amx     = zeros(nc, nlevsoi + 1)
    bmx     = zeros(nc, nlevsoi + 1)
    cmx     = zeros(nc, nlevsoi + 1)
    rmx     = zeros(nc, nlevsoi + 1)
    dwat2   = zeros(nc, nlevsoi + 1)
    smp_grad = zeros(nc, nlevsoi + 1)

    sdamp = 0.0

    # Convert to mm and compute ice fractions
    for j in 1:nlevsoi
        for c in eachindex(mask_hydrology)
            mask_hydrology[c] || continue
            zmm[c, j]  = z[c, joff + j] * 1.0e3
            dzmm[c, j] = dz[c, joff + j] * 1.0e3
            zimm_arr[c, j+1] = zi[c, joff_zi + j] * 1.0e3  # j+1 because index 1 = Fortran j=0

            vol_ice[c, j] = min(watsat[c, j], h2osoi_ice[c, joff + j] / (dz[c, joff + j] * DENICE))
            icefrac[c, j] = min(1.0, vol_ice[c, j] / watsat[c, j])
            vwc_liq[c, j] = max(h2osoi_liq[c, joff + j], 1.0e-6) / (dz[c, joff + j] * DENH2O)
        end
    end

    for c in eachindex(mask_hydrology)
        mask_hydrology[c] || continue
        zimm_arr[c, 1] = 0.0  # Fortran zimm(c,0) = 0
        zwtmm[c] = zwt[c] * 1.0e3
    end

    # Compute jwt index — layer right above water table
    for c in eachindex(mask_hydrology)
        mask_hydrology[c] || continue
        jwt[c] = nlevsoi
        for j in 1:nlevsoi
            if zwt[c] <= zi[c, joff_zi + j]
                jwt[c] = j - 1
                break
            end
        end

        # Compute vwc at water table depth
        vwc_zwt[c] = watsat[c, nlevsoi]
        if t_soisno[c, joff + jwt[c]+1] < TFRZ
            vwc_zwt[c] = vwc_liq[c, nlevsoi]
            for j in nlevsoi:nlevgrnd
                if zwt[c] <= zi[c, joff_zi + j]
                    smp1_val = HFUS * (TFRZ - t_soisno[c, joff + j]) / (GRAV * t_soisno[c, joff + j]) * 1000.0
                    smp1_val = max(sucsat[c, nlevsoi], smp1_val)
                    vwc_zwt[c] = watsat[c, nlevsoi] * (smp1_val / sucsat[c, nlevsoi])^(-1.0 / bsw[c, nlevsoi])
                    vwc_zwt[c] = min(vwc_zwt[c], 0.5 * (watsat[c, nlevsoi] + h2osoi_vol[c, nlevsoi]))
                    break
                end
            end
        end
    end

    # Calculate equilibrium water content based on water table depth
    for j in 1:nlevsoi
        for c in eachindex(mask_hydrology)
            mask_hydrology[c] || continue
            # zimm_arr index: j maps to Fortran zimm(c,j-1), j+1 maps to Fortran zimm(c,j)
            zimm_jm1 = zimm_arr[c, j]      # Fortran zimm(c,j-1)
            zimm_j   = zimm_arr[c, j+1]    # Fortran zimm(c,j)

            if zwtmm[c] <= zimm_jm1
                vol_eq[c, j] = watsat[c, j]
            elseif zwtmm[c] < zimm_j && zwtmm[c] > zimm_jm1
                tempi = 1.0
                temp0 = ((sucsat[c, j] + zwtmm[c] - zimm_jm1) / sucsat[c, j])^(1.0 - 1.0 / bsw[c, j])
                voleq1 = -sucsat[c, j] * watsat[c, j] / (1.0 - 1.0 / bsw[c, j]) / (zwtmm[c] - zimm_jm1) * (tempi - temp0)
                vol_eq[c, j] = (voleq1 * (zwtmm[c] - zimm_jm1) + watsat[c, j] * (zimm_j - zwtmm[c])) / (zimm_j - zimm_jm1)
                vol_eq[c, j] = min(watsat[c, j], vol_eq[c, j])
                vol_eq[c, j] = max(vol_eq[c, j], 0.0)
            else
                tempi = ((sucsat[c, j] + zwtmm[c] - zimm_j) / sucsat[c, j])^(1.0 - 1.0 / bsw[c, j])
                temp0 = ((sucsat[c, j] + zwtmm[c] - zimm_jm1) / sucsat[c, j])^(1.0 - 1.0 / bsw[c, j])
                vol_eq[c, j] = -sucsat[c, j] * watsat[c, j] / (1.0 - 1.0 / bsw[c, j]) / (zimm_j - zimm_jm1) * (tempi - temp0)
                vol_eq[c, j] = max(vol_eq[c, j], 0.0)
                vol_eq[c, j] = min(watsat[c, j], vol_eq[c, j])
            end
            zq[c, j] = -sucsat[c, j] * (max(vol_eq[c, j] / watsat[c, j], 0.01))^(-bsw[c, j])
            zq[c, j] = max(smpmin[c], zq[c, j])
        end
    end

    # If water table is below soil column, calculate zq for the 11th layer
    j = nlevsoi
    for c in eachindex(mask_hydrology)
        mask_hydrology[c] || continue
        if jwt[c] == nlevsoi
            zimm_j = zimm_arr[c, j+1]  # Fortran zimm(c,j)
            tempi = 1.0
            temp0 = ((sucsat[c, j] + zwtmm[c] - zimm_j) / sucsat[c, j])^(1.0 - 1.0 / bsw[c, j])
            vol_eq[c, j+1] = -sucsat[c, j] * watsat[c, j] / (1.0 - 1.0 / bsw[c, j]) / (zwtmm[c] - zimm_j) * (tempi - temp0)
            vol_eq[c, j+1] = max(vol_eq[c, j+1], 0.0)
            vol_eq[c, j+1] = min(watsat[c, j], vol_eq[c, j+1])
            zq[c, j+1] = -sucsat[c, j] * (max(vol_eq[c, j+1] / watsat[c, j], 0.01))^(-bsw[c, j])
            zq[c, j+1] = max(smpmin[c], zq[c, j+1])
        end
    end

    # Hydraulic conductivity and soil matric potential
    for j in 1:nlevsoi
        for c in eachindex(mask_hydrology)
            mask_hydrology[c] || continue

            s1 = 0.5 * (vwc_liq[c, j] + vwc_liq[c, min(nlevsoi, j+1)]) /
                 (0.5 * (watsat[c, j] + watsat[c, min(nlevsoi, j+1)]))
            s1 = min(1.0, s1)
            s2 = hksat[c, j] * s1^(2.0 * bsw[c, j] + 2.0)

            imped_arr[c, j] = 10.0^(-cfg.e_ice * (0.5 * (icefrac[c, j] + icefrac[c, min(nlevsoi, j+1)])))

            hk[c, j] = imped_arr[c, j] * s1 * s2
            dhkdw[c, j] = imped_arr[c, j] * (2.0 * bsw[c, j] + 3.0) * s2 *
                 (1.0 / (watsat[c, j] + watsat[c, min(nlevsoi, j+1)]))

            # Matric potential and derivative
            s_node = max(vwc_liq[c, j] / watsat[c, j], 0.01)
            s_node = min(1.0, s_node)

            smp_arr[c, j] = -sucsat[c, j] * s_node^(-bsw[c, j])
            smp_arr[c, j] = max(smpmin[c], smp_arr[c, j])

            dsmpdw[c, j] = -bsw[c, j] * smp_arr[c, j] / vwc_liq[c, j]

            smp_l[c, j] = smp_arr[c, j]
            hk_l[c, j] = hk[c, j]
        end
    end

    # Aquifer (11th) layer
    for c in eachindex(mask_hydrology)
        mask_hydrology[c] || continue
        zmm[c, nlevsoi+1] = 0.5 * (1.0e3 * zwt[c] + zmm[c, nlevsoi])
        if jwt[c] < nlevsoi
            dzmm[c, nlevsoi+1] = dzmm[c, nlevsoi]
        else
            dzmm[c, nlevsoi+1] = (1.0e3 * zwt[c] - zmm[c, nlevsoi])
        end
    end

    # ---- Set up tridiagonal system ----

    # Node j=1 (top)
    j = 1
    for c in eachindex(mask_hydrology)
        mask_hydrology[c] || continue
        qin[c, j] = qflx_infl[c]
        den = zmm[c, j+1] - zmm[c, j]
        dzq = zq[c, j+1] - zq[c, j]
        num = (smp_arr[c, j+1] - smp_arr[c, j]) - dzq
        qout_arr[c, j] = -hk[c, j] * num / den
        dqodw1[c, j] = -(-hk[c, j] * dsmpdw[c, j] + num * dhkdw[c, j]) / den
        dqodw2[c, j] = -(hk[c, j] * dsmpdw[c, j+1] + num * dhkdw[c, j]) / den
        rmx[c, j] = qin[c, j] - qout_arr[c, j] - qflx_rootsoi[c, j]
        amx[c, j] = 0.0
        bmx[c, j] = dzmm[c, j] * (sdamp + 1.0 / dtime) + dqodw1[c, j]
        cmx[c, j] = dqodw2[c, j]
    end

    # Nodes j=2 to j=nlevsoi-1
    for j in 2:(nlevsoi-1)
        for c in eachindex(mask_hydrology)
            mask_hydrology[c] || continue
            den = zmm[c, j] - zmm[c, j-1]
            dzq = zq[c, j] - zq[c, j-1]
            num = (smp_arr[c, j] - smp_arr[c, j-1]) - dzq
            qin[c, j] = -hk[c, j-1] * num / den
            dqidw0[c, j] = -(-hk[c, j-1] * dsmpdw[c, j-1] + num * dhkdw[c, j-1]) / den
            dqidw1[c, j] = -(hk[c, j-1] * dsmpdw[c, j] + num * dhkdw[c, j-1]) / den
            den = zmm[c, j+1] - zmm[c, j]
            dzq = zq[c, j+1] - zq[c, j]
            num = (smp_arr[c, j+1] - smp_arr[c, j]) - dzq
            qout_arr[c, j] = -hk[c, j] * num / den
            dqodw1[c, j] = -(-hk[c, j] * dsmpdw[c, j] + num * dhkdw[c, j]) / den
            dqodw2[c, j] = -(hk[c, j] * dsmpdw[c, j+1] + num * dhkdw[c, j]) / den
            rmx[c, j] = qin[c, j] - qout_arr[c, j] - qflx_rootsoi[c, j]
            amx[c, j] = -dqidw0[c, j]
            bmx[c, j] = dzmm[c, j] / dtime - dqidw1[c, j] + dqodw1[c, j]
            cmx[c, j] = dqodw2[c, j]
        end
    end

    # Node j=nlevsoi (bottom)
    j = nlevsoi
    for c in eachindex(mask_hydrology)
        mask_hydrology[c] || continue
        if j > jwt[c]  # water table is in soil column
            den = zmm[c, j] - zmm[c, j-1]
            dzq = zq[c, j] - zq[c, j-1]
            num = (smp_arr[c, j] - smp_arr[c, j-1]) - dzq
            qin[c, j] = -hk[c, j-1] * num / den
            dqidw0[c, j] = -(-hk[c, j-1] * dsmpdw[c, j-1] + num * dhkdw[c, j-1]) / den
            dqidw1[c, j] = -(hk[c, j-1] * dsmpdw[c, j] + num * dhkdw[c, j-1]) / den
            qout_arr[c, j] = 0.0
            dqodw1[c, j] = 0.0
            rmx[c, j] = qin[c, j] - qout_arr[c, j] - qflx_rootsoi[c, j]
            amx[c, j] = -dqidw0[c, j]
            bmx[c, j] = dzmm[c, j] / dtime - dqidw1[c, j] + dqodw1[c, j]
            cmx[c, j] = 0.0

            # Aquifer layer — hydrologically inactive
            rmx[c, j+1] = 0.0
            amx[c, j+1] = 0.0
            bmx[c, j+1] = dzmm[c, j+1] / dtime
            cmx[c, j+1] = 0.0
        else  # water table is below soil column
            # Compute aquifer soil moisture
            s_node = max(0.5 * ((vwc_zwt[c] + vwc_liq[c, j]) / watsat[c, j]), 0.01)
            s_node = min(1.0, s_node)

            # Compute smp for aquifer layer
            smp1 = -sucsat[c, j] * s_node^(-bsw[c, j])
            smp1 = max(smpmin[c], smp1)

            # Compute dsmpdw for aquifer layer
            dsmpdw1 = -bsw[c, j] * smp1 / (s_node * watsat[c, j])

            # Bottom layer of soil column
            den = zmm[c, j] - zmm[c, j-1]
            dzq = zq[c, j] - zq[c, j-1]
            num = (smp_arr[c, j] - smp_arr[c, j-1]) - dzq
            qin[c, j] = -hk[c, j-1] * num / den
            dqidw0[c, j] = -(-hk[c, j-1] * dsmpdw[c, j-1] + num * dhkdw[c, j-1]) / den
            dqidw1[c, j] = -(hk[c, j-1] * dsmpdw[c, j] + num * dhkdw[c, j-1]) / den
            den = zmm[c, j+1] - zmm[c, j]
            dzq = zq[c, j+1] - zq[c, j]
            num = (smp1 - smp_arr[c, j]) - dzq
            qout_arr[c, j] = -hk[c, j] * num / den
            dqodw1[c, j] = -(-hk[c, j] * dsmpdw[c, j] + num * dhkdw[c, j]) / den
            dqodw2[c, j] = -(hk[c, j] * dsmpdw1 + num * dhkdw[c, j]) / den

            rmx[c, j] = qin[c, j] - qout_arr[c, j] - qflx_rootsoi[c, j]
            amx[c, j] = -dqidw0[c, j]
            bmx[c, j] = dzmm[c, j] / dtime - dqidw1[c, j] + dqodw1[c, j]
            cmx[c, j] = dqodw2[c, j]

            # Aquifer layer
            qin[c, j+1] = qout_arr[c, j]
            dqidw0[c, j+1] = -(-hk[c, j] * dsmpdw[c, j] + num * dhkdw[c, j]) / den
            dqidw1[c, j+1] = -(hk[c, j] * dsmpdw1 + num * dhkdw[c, j]) / den
            qout_arr[c, j+1] = 0.0
            dqodw1[c, j+1] = 0.0
            rmx[c, j+1] = qin[c, j+1] - qout_arr[c, j+1]
            amx[c, j+1] = -dqidw0[c, j+1]
            bmx[c, j+1] = dzmm[c, j+1] / dtime - dqidw1[c, j+1] + dqodw1[c, j+1]
            cmx[c, j+1] = 0.0
        end
    end

    # Solve for dwat using multi-column tridiagonal solver
    jtop_arr = fill(1, nc)
    mask_jtop = mask_hydrology
    tridiagonal_multi!(dwat2, amx, bmx, cmx, rmx,
        jtop_arr, mask_jtop, nc, nlevsoi + 1)

    # Renew the mass of liquid water and compute qcharge
    for c in eachindex(mask_hydrology)
        mask_hydrology[c] || continue

        for j in 1:nlevsoi
            h2osoi_liq[c, joff + j] = h2osoi_liq[c, joff + j] + dwat2[c, j] * dzmm[c, j]
        end

        # Calculate qcharge
        if jwt[c] < nlevsoi
            wh_zwt = 0.0

            s_node = max(h2osoi_vol[c, jwt[c]+1] / watsat[c, jwt[c]+1], 0.01)
            s1 = min(1.0, s_node)

            ka = imped_arr[c, jwt[c]+1] * hksat[c, jwt[c]+1] *
                 s1^(2.0 * bsw[c, jwt[c]+1] + 3.0)

            smp1 = max(smpmin[c], smp_arr[c, max(1, jwt[c])])
            wh = smp1 - zq[c, max(1, jwt[c])]

            if jwt[c] == 0
                qcharge[c] = -ka * (wh_zwt - wh) / ((zwt[c] + 1.0e-3) * 1000.0)
            else
                qcharge[c] = -ka * (wh_zwt - wh) / ((zwt[c] - z[c, joff + jwt[c]]) * 1000.0 * 2.0)
            end

            # Limit qcharge
            qcharge[c] = max(-10.0 / dtime, qcharge[c])
            qcharge[c] = min(10.0 / dtime, qcharge[c])
        else
            # Water table below soil column
            qcharge[c] = dwat2[c, nlevsoi+1] * dzmm[c, nlevsoi+1] / dtime
        end
    end

    # Compute water deficit and reset negative liquid water
    for c in eachindex(mask_hydrology)
        mask_hydrology[c] || continue
        qflx_deficit[c] = 0.0
        for j in 1:nlevsoi
            if h2osoi_liq[c, joff + j] < 0.0
                qflx_deficit[c] = qflx_deficit[c] - h2osoi_liq[c, joff + j]
            end
        end
    end

    nothing
end

# ===========================================================================
# soil_water! — main entry point
# ===========================================================================
"""
    soil_water!(col_data, mask_hydrology, mask_urban,
        soilhydrology, soilstate, waterfluxbulk, waterstatebulk,
        temperature, canopystate, energyflux, swrc, cfg, dtime;
        use_flexibleCN=false)

Main entry point for soil water movement. Dispatches to the appropriate
solver method based on `cfg.soilwater_movement_method`.

Ported from `SoilWater` in `SoilWaterMovementMod.F90`.
"""
function soil_water!(col_data::ColumnData,
        mask_hydrology::BitVector, mask_urban::BitVector,
        soilhydrology::SoilHydrologyData, soilstate::SoilStateData,
        waterfluxbulk::WaterFluxBulkData, waterstatebulk::WaterStateBulkData,
        temperature::TemperatureData, canopystate::CanopyStateData,
        energyflux::EnergyFluxData,
        swrc::SoilWaterRetentionCurve, cfg::SoilWaterMovementConfig,
        dtime::Float64;
        use_flexibleCN::Bool = false)

    nlevsoi = varpar.nlevsoi

    wa          = waterstatebulk.ws.wa_col
    dz          = col_data.dz
    h2osoi_ice  = waterstatebulk.ws.h2osoi_ice_col
    h2osoi_vol  = waterstatebulk.ws.h2osoi_vol_col
    h2osoi_liq  = waterstatebulk.ws.h2osoi_liq_col

    if cfg.soilwater_movement_method == ZENGDECKER_2009
        soilwater_zengdecker2009!(col_data, mask_hydrology, mask_urban,
            soilhydrology, soilstate, waterfluxbulk, waterstatebulk,
            temperature, canopystate, energyflux, swrc, cfg, dtime)

    elseif cfg.soilwater_movement_method == MOISTURE_FORM
        soilwater_moisture_form!(col_data, mask_hydrology, mask_urban,
            soilhydrology, soilstate, waterfluxbulk, waterstatebulk,
            temperature, canopystate, energyflux, swrc, cfg, dtime)

    else
        error("soil_water!: a SoilWater implementation must be specified!")
    end

    # FlexibleCN workaround for negative liquid water
    if use_flexibleCN
        nlevsno_  = varpar.nlevsno
        joff_     = nlevsno_
        watmin = 0.001

        for j in 1:(nlevsoi-1)
            for c in eachindex(mask_hydrology)
                mask_hydrology[c] || continue
                if h2osoi_liq[c, joff_ + j] < 0.0
                    xs = watmin - h2osoi_liq[c, joff_ + j]
                else
                    xs = 0.0
                end
                h2osoi_liq[c, joff_ + j]   = h2osoi_liq[c, joff_ + j] + xs
                h2osoi_liq[c, joff_ + j+1] = h2osoi_liq[c, joff_ + j+1] - xs
            end
        end

        j = nlevsoi
        for c in eachindex(mask_hydrology)
            mask_hydrology[c] || continue
            if h2osoi_liq[c, joff_ + j] < watmin
                xs = watmin - h2osoi_liq[c, joff_ + j]
            else
                xs = 0.0
            end
            h2osoi_liq[c, joff_ + j] = h2osoi_liq[c, joff_ + j] + xs
            wa[c] = wa[c] - xs
        end

        # Update volumetric soil moisture
        for j in 1:nlevsoi
            for c in eachindex(mask_hydrology)
                mask_hydrology[c] || continue
                h2osoi_vol[c, j] = h2osoi_liq[c, joff_ + j] / (dz[c, joff_ + j] * DENH2O) +
                                   h2osoi_ice[c, joff_ + j] / (dz[c, joff_ + j] * DENICE)
            end
        end
    end

    nothing
end
