# ==========================================================================
# UrbanFluxes — ported from UrbanFluxesMod.F90
#
# Calculate turbulent fluxes for urban landunit (roof, sunwall, shadewall,
# pervious and impervious road).
# ==========================================================================

# ---------------------------------------------------------------------------
# Module-level parameters (read from parameter file in Fortran)
# ---------------------------------------------------------------------------

Base.@kwdef mutable struct UrbanFluxesParamsData
    wind_min::Float64 = 0.0  # Minimum wind speed at atmospheric forcing height (m/s)
end

const urban_fluxes_params = UrbanFluxesParamsData()

"""
    urban_fluxes_read_params!(; wind_min)

Set UrbanFluxes module parameters (replaces Fortran readParams).
"""
function urban_fluxes_read_params!(; wind_min::Real = urban_fluxes_params.wind_min)
    urban_fluxes_params.wind_min = wind_min
    return nothing
end

# ---------------------------------------------------------------------------
# simple_wasteheatfromac — Calculate waste heat from AC (simple CLM4.5 method)
# ---------------------------------------------------------------------------

"""
    simple_wasteheatfromac(eflx_urban_ac, eflx_urban_heat)
        -> (eflx_wasteheat, eflx_heat_from_ac)

Calculate waste heat from air-conditioning with the simpler method (CLM4.5).

Ported from: UrbanFluxesMod.F90 :: simple_wasteheatfromac
"""
function simple_wasteheatfromac(eflx_urban_ac::Real, eflx_urban_heat::Real)
    # wasteheat from heating/cooling
    if urban_ctrl.urban_hac == URBAN_WASTEHEAT_ON
        eflx_wasteheat = AC_WASTEHEAT_FACTOR * eflx_urban_ac +
                          HT_WASTEHEAT_FACTOR * eflx_urban_heat
    else
        eflx_wasteheat = 0.0
    end

    # If air conditioning on, always replace heat removed with heat into canyon
    if urban_ctrl.urban_hac == URBAN_HAC_ON || urban_ctrl.urban_hac == URBAN_WASTEHEAT_ON
        eflx_heat_from_ac = smooth_abs(eflx_urban_ac)
    else
        eflx_heat_from_ac = 0.0
    end

    return (eflx_wasteheat, eflx_heat_from_ac)
end

# ---------------------------------------------------------------------------
# wasteheat! — Calculate wasteheat from urban heating/cooling
# ---------------------------------------------------------------------------

"""
    wasteheat!(energyflux, lun_data, num_urbanl, filter_urbanl,
               eflx_wasteheat_roof, eflx_wasteheat_sunwall,
               eflx_wasteheat_shadewall, eflx_heat_from_ac_roof,
               eflx_heat_from_ac_sunwall, eflx_heat_from_ac_shadewall)

Calculate the wasteheat flux from urban heating or air-conditioning.

Ported from: UrbanFluxesMod.F90 :: wasteheat
"""
function wasteheat!(
        energyflux       ::EnergyFluxData,
        lun_data         ::LandunitData,
        num_urbanl       ::Int,
        filter_urbanl    ::Vector{Int},
        eflx_wasteheat_roof      ::Vector{<:Real},
        eflx_wasteheat_sunwall   ::Vector{<:Real},
        eflx_wasteheat_shadewall ::Vector{<:Real},
        eflx_heat_from_ac_roof      ::Vector{<:Real},
        eflx_heat_from_ac_sunwall   ::Vector{<:Real},
        eflx_heat_from_ac_shadewall ::Vector{<:Real})

    for fl in 1:num_urbanl
        l = filter_urbanl[fl]
        g = lun_data.gridcell[l]

        if is_simple_build_temp()
            # Total waste heat and heat from AC is sum of heat for walls and roofs
            # accounting for different surface areas
            energyflux.eflx_wasteheat_lun[l] = lun_data.wtlunit_roof[l] * eflx_wasteheat_roof[l] +
                (1.0 - lun_data.wtlunit_roof[l]) * (lun_data.canyon_hwr[l] * (eflx_wasteheat_sunwall[l] +
                eflx_wasteheat_shadewall[l]))

        elseif is_prog_build_temp()
            # wasteheat from heating/cooling
            if urban_ctrl.urban_hac == URBAN_WASTEHEAT_ON
                energyflux.eflx_wasteheat_lun[l] = AC_WASTEHEAT_FACTOR * energyflux.eflx_urban_ac_lun[l] +
                    HT_WASTEHEAT_FACTOR * energyflux.eflx_urban_heat_lun[l]
            else
                energyflux.eflx_wasteheat_lun[l] = 0.0
            end
        end

        # Limit wasteheat to ensure no unrealistically strong positive feedbacks
        energyflux.eflx_wasteheat_lun[l] = smooth_min(energyflux.eflx_wasteheat_lun[l], WASTEHEAT_LIMIT)

        if is_simple_build_temp()
            energyflux.eflx_heat_from_ac_lun[l] = lun_data.wtlunit_roof[l] * eflx_heat_from_ac_roof[l] +
                (1.0 - lun_data.wtlunit_roof[l]) * (lun_data.canyon_hwr[l] * (eflx_heat_from_ac_sunwall[l] +
                eflx_heat_from_ac_shadewall[l]))

        elseif is_prog_build_temp()
            # If air conditioning on, always replace heat removed with heat into canyon
            if urban_ctrl.urban_hac == URBAN_HAC_ON || urban_ctrl.urban_hac == URBAN_WASTEHEAT_ON
                energyflux.eflx_heat_from_ac_lun[l] = smooth_abs(energyflux.eflx_urban_ac_lun[l])
            else
                energyflux.eflx_heat_from_ac_lun[l] = 0.0
            end
        end
    end

    return nothing
end

# ---------------------------------------------------------------------------
# calc_simple_internal_building_temp! — Simple building temp (CLM4.5)
# ---------------------------------------------------------------------------

"""
    calc_simple_internal_building_temp!(temperature, col_data, lun_data,
        num_urbanc, filter_urbanc, num_urbanl, filter_urbanl)

Calculate the internal building temperature by simpler method (CLM4.5).

Ported from: UrbanFluxesMod.F90 :: calc_simple_internal_building_temp
"""
function calc_simple_internal_building_temp!(
        temperature   ::TemperatureData,
        col_data      ::ColumnData,
        lun_data      ::LandunitData,
        num_urbanc    ::Int,
        filter_urbanc ::Vector{Int},
        num_urbanl    ::Int,
        filter_urbanl ::Vector{Int})

    nlevurb = varpar.nlevurb

    nl = length(lun_data.gridcell)
    FT = eltype(temperature.t_soisno_col)
    t_sunwall_innerl   = zeros(FT, nl)
    t_shadewall_innerl = zeros(FT, nl)
    t_roof_innerl      = zeros(FT, nl)

    # Gather inner layer temperatures from urban columns
    for fc in 1:num_urbanc
        c = filter_urbanc[fc]
        l = col_data.landunit[c]

        if col_data.itype[c] == ICOL_ROOF
            t_roof_innerl[l] = temperature.t_soisno_col[c, nlevurb]
        elseif col_data.itype[c] == ICOL_SUNWALL
            t_sunwall_innerl[l] = temperature.t_soisno_col[c, nlevurb]
        elseif col_data.itype[c] == ICOL_SHADEWALL
            t_shadewall_innerl[l] = temperature.t_soisno_col[c, nlevurb]
        end
    end

    # Calculate internal building temperature
    for fl in 1:num_urbanl
        l = filter_urbanl[fl]

        lngth_roof = (lun_data.ht_roof[l] / lun_data.canyon_hwr[l]) *
            lun_data.wtlunit_roof[l] / (1.0 - lun_data.wtlunit_roof[l])
        temperature.t_building_lun[l] = (lun_data.ht_roof[l] * (t_shadewall_innerl[l] + t_sunwall_innerl[l]) +
            lngth_roof * t_roof_innerl[l]) / (2.0 * lun_data.ht_roof[l] + lngth_roof)
    end

    return nothing
end

# ---------------------------------------------------------------------------
# GPU kernels for the column-INDEPENDENT, element-wise passes of urban_fluxes!.
#
# These passes have no cross-index coupling (each landunit/column/patch writes
# only its own slot) and run WHOLE on the device. The iterative canyon-air solve
# (canyontop wind, stability iteration via friction_velocity!'s landunit path,
# the urban-column gather/reduction into taf/qaf, canyon-surface fluxes and the
# sum()-based error checks) is loop-carried + reduction-coupled and stays HOST.
#
# Masks (Vector{Bool} over the full index range) replace the integer filters for
# the kernelized passes so MtlArray/device-Bool dispatch works; the host core
# keeps using the integer filters. Literals are eltype-converted (zero(T)/T(...))
# so no Float64 reaches a Float32-only backend; on Float64 this is byte-identical.
# ---------------------------------------------------------------------------

# Pass 1: restart fields for NON-urban landunits (taf/qaf = SPVAL).
@kernel function _uf_nourban_restart_kernel!(taf_lun, qaf_lun, @Const(mask), spval)
    l = @index(Global)
    @inbounds if mask[l]
        taf_lun[l] = spval
        qaf_lun[l] = spval
    end
end

# Pass 2 (set-constants beta/zii) is NOT kernelized: beta_arr/zii_arr are host-only
# local work arrays consumed exclusively by the host-resident canyon-air solve.

# Pass 3: urban roots — only pervious road has roots, all others zero.
# Per-patch; the level loop is sequential-independent and lives inside the thread.
@kernel function _uf_roots_kernel!(rootr_patch, @Const(mask), @Const(column),
                                   @Const(itype), @Const(rootr_road_perv_col),
                                   nlevgrnd::Int)
    p = @index(Global)
    @inbounds if mask[p]
        T = eltype(rootr_patch)
        c = column[p]
        if itype[c] == ICOL_ROAD_PERV
            for j in 1:nlevgrnd
                rootr_patch[p, j] = rootr_road_perv_col[c, j]
            end
        else
            for j in 1:nlevgrnd
                rootr_patch[p, j] = zero(T)
            end
        end
    end
end

# Pass 4: 2-m temperature, humidity, RH and t_veg diagnostics (per urban patch).
# Bundle keeps the kernel under Metal's ~31-arg cap; field names mirror the
# original variable paths so the body reads like the scalar loop.
Base.@kwdef struct _UfDiagDV{V,VI}
    t_ref2m_u_patch::V
    t_veg_patch::V
    q_ref2m_patch::V
    rh_ref2m_patch::V
    rh_ref2m_u_patch::V
    taf_lun::V
    qaf_lun::V
    forc_t_grc::V
    forc_pbot_grc::V
    gridcell::VI
    landunit::VI
end
Adapt.@adapt_structure _UfDiagDV

# `t_ref2m_patch` is the standalone first arg (drives backend detection in _launch!);
# the rest of the per-patch fields ride in the bundle (keeps under Metal's arg cap).
@kernel function _uf_diag_kernel!(t_ref2m_patch, @Const(mask), o::_UfDiagDV)
    p = @index(Global)
    @inbounds if mask[p]
        T = eltype(t_ref2m_patch)
        g = o.gridcell[p]
        l = o.landunit[p]

        t_ref2m_patch[p]     = o.taf_lun[l]
        o.q_ref2m_patch[p]   = o.qaf_lun[l]
        o.t_ref2m_u_patch[p] = o.taf_lun[l]

        (qsat_ref2m, _, _, _) = qsat(t_ref2m_patch[p], o.forc_pbot_grc[g])
        o.rh_ref2m_patch[p] = smooth_min(T(100.0),
            o.q_ref2m_patch[p] / qsat_ref2m * T(100.0))
        o.rh_ref2m_u_patch[p] = o.rh_ref2m_patch[p]

        o.t_veg_patch[p] = o.forc_t_grc[g]
    end
end

# Move a host forcing vector onto the backend (and eltype) of a state prototype so
# a kernel reading it dispatches on-device. CPU prototype -> converted host Vector
# (byte-identical to the Float64 source on Float64 CPU); device prototype -> device
# array at the prototype's working precision.
function _uf_to_backend(proto::AbstractArray, v::AbstractVector{<:Real})
    T = eltype(proto)
    h = convert(Vector{T}, collect(v))
    d = similar(proto, T, length(h))
    copyto!(d, h)
    return d
end

# Build a length-n Vector{Bool} mask from an integer filter, then move it to the
# backend of `proto` (a state array) so kernels dispatch onto the same device.
# CPU prototype -> plain Vector{Bool}; device prototype -> device Bool array.
function _uf_device_mask(proto::AbstractArray, filter::AbstractVector{<:Integer},
                         fn::Int, n::Int)
    h = Vector{Bool}(undef, n); fill!(h, false)
    @inbounds for f in 1:fn
        h[filter[f]] = true
    end
    d = similar(proto, Bool, n)
    copyto!(d, h)
    return d
end

"""
    urban_fluxes_diagnostics!(temperature, waterdiagbulk, soilstate, patch_data,
        col_data, m_nourbanl, m_urbanp, forc_t_grc, forc_pbot_grc, nlevgrnd) -> nothing

Run the three column-INDEPENDENT, device-resident passes of urban_fluxes! as a
single whole-function unit: (1) the non-urban taf/qaf restart, (2) the urban-roots
fill (pervious-road only), and (3) the 2-m temperature/humidity/RH/t_veg
diagnostics. Every loop is a KernelAbstractions kernel, so this runs WHOLE on the
device (Metal) when the state structs and masks live there, and is CPU
byte-identical on the host path. The iterative canyon-air solve in urban_fluxes!
is loop-carried + reduction-coupled and remains on the host (see that function).
"""
function urban_fluxes_diagnostics!(
        temperature, waterdiagbulk, soilstate, patch_data, col_data,
        m_nourbanl::AbstractVector{Bool}, m_urbanp::AbstractVector{Bool},
        forc_t_grc::AbstractVector{<:Real}, forc_pbot_grc::AbstractVector{<:Real},
        nlevgrnd::Int)

    FT = eltype(temperature.taf_lun)
    spval = convert(FT, SPVAL)
    _launch!(_uf_nourban_restart_kernel!, temperature.taf_lun,
             waterdiagbulk.qaf_lun, m_nourbanl, spval)

    _launch!(_uf_roots_kernel!, soilstate.rootr_patch, m_urbanp, patch_data.column,
             col_data.itype, soilstate.rootr_road_perv_col, nlevgrnd;
             ndrange=size(soilstate.rootr_patch, 1))

    proto_p = temperature.t_ref2m_patch
    fpbot_d = _uf_to_backend(proto_p, forc_pbot_grc)
    ft_d    = _uf_to_backend(proto_p, forc_t_grc)
    diagdv = _UfDiagDV(
        t_ref2m_u_patch  = temperature.t_ref2m_u_patch,
        t_veg_patch      = temperature.t_veg_patch,
        q_ref2m_patch    = waterdiagbulk.q_ref2m_patch,
        rh_ref2m_patch   = waterdiagbulk.rh_ref2m_patch,
        rh_ref2m_u_patch = waterdiagbulk.rh_ref2m_u_patch,
        taf_lun          = temperature.taf_lun,
        qaf_lun          = waterdiagbulk.qaf_lun,
        forc_t_grc       = ft_d,
        forc_pbot_grc    = fpbot_d,
        gridcell         = patch_data.gridcell,
        landunit         = patch_data.landunit)
    _launch!(_uf_diag_kernel!, temperature.t_ref2m_patch, m_urbanp, diagdv;
             ndrange=length(temperature.t_ref2m_patch))
    return nothing
end

# ---------------------------------------------------------------------------
# Canyon-surface flux pass — per-patch, column-INDEPENDENT given the SOLVED
# per-landunit canopy-air state (taf/qaf, the resistances ramu/zeta and the
# conductance weights wt{u,s}*). Each urban patch reads only its own landunit's
# already-solved scalars + its own column state + gridcell forcing and writes only
# its own patch slot, so the whole block (original lines "Determine fluxes from
# canyon surfaces") runs WHOLE on the device as a single KernelAbstractions kernel.
#
# `dth`/`dqh` in the original were per-landunit scratch reused inside the per-patch
# loop; here they are PATCH-LOCAL temporaries (read+written within the same thread),
# which removes the false landunit coupling without changing results.
#
# The per-landunit weight arrays are bundled into @adapt_structure device-view
# structs (Metal caps total kernel args ~31). Literals are eltype-converted so no
# Float64 reaches a Float32-only backend; on Float64 CPU this is byte-identical.
# ---------------------------------------------------------------------------

# Per-landunit conductance weights, split into sensible (S) and latent (Q) bundles
# plus the "_unscl" companions used to back out each surface's own contribution.
Base.@kwdef struct _UfWtsDV{V}
    wtas::V
    wtus_roof::V
    wtus_road_perv::V
    wtus_road_imperv::V
    wtus_sunwall::V
    wtus_shadewall::V
    wtus_roof_unscl::V
    wtus_road_perv_unscl::V
    wtus_road_imperv_unscl::V
    wtus_sunwall_unscl::V
    wtus_shadewall_unscl::V
    wts_sum::V
end
Adapt.@adapt_structure _UfWtsDV

Base.@kwdef struct _UfWtqDV{V}
    wtaq::V
    wtuq_roof::V
    wtuq_road_perv::V
    wtuq_road_imperv::V
    wtuq_sunwall::V
    wtuq_shadewall::V
    wtuq_roof_unscl::V
    wtuq_road_perv_unscl::V
    wtuq_road_imperv_unscl::V
    wtq_sum::V
end
Adapt.@adapt_structure _UfWtqDV

# Patch-level energyflux/waterflux outputs + the per-landunit solved scalars + the
# per-column state and gridcell forcing the kernel reads. `cgrnd_patch` is the
# standalone first kernel arg (drives backend detection in _launch!).
Base.@kwdef struct _UfCanyonDV{V,VI}
    cgrnds_patch::V
    cgrndl_patch::V
    taux_patch::V
    tauy_patch::V
    ulrad_patch::V
    dlrad_patch::V
    eflx_sh_grnd_patch::V
    eflx_sh_snow_patch::V
    eflx_sh_soil_patch::V
    eflx_sh_h2osfc_patch::V
    eflx_sh_tot_patch::V
    eflx_sh_tot_u_patch::V
    qflx_evap_soi_patch::V
    qflx_tran_veg_patch::V
    qflx_evap_veg_patch::V
    ram1_patch::V
    zeta_patch::V
    eflx_sh_grnd_scale::V
    qflx_evap_soi_scale::V
    # per-landunit solved scalars
    ramu::V
    zeta_lunit::V
    taf_lun::V
    qaf_lun::V
    # per-column state
    t_grnd_col::V
    qg_col::V
    htvp_col::V
    dqgdT_col::V
    frac_sno_col::V
    soilalpha_u_col::V
    wtus_col::V
    wtuq_col::V
    itype::VI
    # gridcell forcing
    forc_rho_grc::V
    forc_u_grc::V
    forc_v_grc::V
    # index maps
    column::VI
    gridcell::VI
    landunit::VI
end
Adapt.@adapt_structure _UfCanyonDV

@kernel function _uf_canyon_surface_kernel!(cgrnd_patch, @Const(mask),
        o::_UfCanyonDV, s::_UfWtsDV, q::_UfWtqDV, CPAIR_)
    p = @index(Global)
    @inbounds if mask[p]
        T = eltype(cgrnd_patch)
        c = o.column[p]
        g = o.gridcell[p]
        l = o.landunit[p]
        cpair = T(CPAIR_)

        o.ram1_patch[p] = o.ramu[l]
        o.zeta_patch[p] = o.zeta_lunit[l]

        # Upward and downward canopy longwave are zero
        o.ulrad_patch[p] = zero(T)
        o.dlrad_patch[p] = zero(T)

        ct = o.itype[c]

        # Derivative of sensible and latent heat fluxes wrt ground temperature
        if ct == ICOL_ROOF
            o.cgrnds_patch[p] = o.forc_rho_grc[g] * cpair *
                (s.wtas[l] + s.wtus_road_perv[l] + s.wtus_road_imperv[l] + s.wtus_sunwall[l] + s.wtus_shadewall[l]) *
                (s.wtus_roof_unscl[l] / s.wts_sum[l])
            o.cgrndl_patch[p] = o.forc_rho_grc[g] *
                (q.wtaq[l] + q.wtuq_road_perv[l] + q.wtuq_road_imperv[l] + q.wtuq_sunwall[l] + q.wtuq_shadewall[l]) *
                (q.wtuq_roof_unscl[l] / q.wtq_sum[l]) * o.dqgdT_col[c]
        elseif ct == ICOL_ROAD_PERV
            o.cgrnds_patch[p] = o.forc_rho_grc[g] * cpair *
                (s.wtas[l] + s.wtus_roof[l] + s.wtus_road_imperv[l] + s.wtus_sunwall[l] + s.wtus_shadewall[l]) *
                (s.wtus_road_perv_unscl[l] / s.wts_sum[l])
            o.cgrndl_patch[p] = o.forc_rho_grc[g] *
                (q.wtaq[l] + q.wtuq_roof[l] + q.wtuq_road_imperv[l] + q.wtuq_sunwall[l] + q.wtuq_shadewall[l]) *
                (q.wtuq_road_perv_unscl[l] / q.wtq_sum[l]) * o.dqgdT_col[c]
        elseif ct == ICOL_ROAD_IMPERV
            o.cgrnds_patch[p] = o.forc_rho_grc[g] * cpair *
                (s.wtas[l] + s.wtus_roof[l] + s.wtus_road_perv[l] + s.wtus_sunwall[l] + s.wtus_shadewall[l]) *
                (s.wtus_road_imperv_unscl[l] / s.wts_sum[l])
            o.cgrndl_patch[p] = o.forc_rho_grc[g] *
                (q.wtaq[l] + q.wtuq_roof[l] + q.wtuq_road_perv[l] + q.wtuq_sunwall[l] + q.wtuq_shadewall[l]) *
                (q.wtuq_road_imperv_unscl[l] / q.wtq_sum[l]) * o.dqgdT_col[c]
        elseif ct == ICOL_SUNWALL
            o.cgrnds_patch[p] = o.forc_rho_grc[g] * cpair *
                (s.wtas[l] + s.wtus_roof[l] + s.wtus_road_perv[l] + s.wtus_road_imperv[l] + s.wtus_shadewall[l]) *
                (s.wtus_sunwall_unscl[l] / s.wts_sum[l])
            o.cgrndl_patch[p] = zero(T)
        elseif ct == ICOL_SHADEWALL
            o.cgrnds_patch[p] = o.forc_rho_grc[g] * cpair *
                (s.wtas[l] + s.wtus_roof[l] + s.wtus_road_perv[l] + s.wtus_road_imperv[l] + s.wtus_sunwall[l]) *
                (s.wtus_shadewall_unscl[l] / s.wts_sum[l])
            o.cgrndl_patch[p] = zero(T)
        end
        cgrnd_patch[p] = o.cgrnds_patch[p] + o.cgrndl_patch[p] * o.htvp_col[c]

        # Surface fluxes of momentum
        o.taux_patch[p] = -o.forc_rho_grc[g] * o.forc_u_grc[g] / o.ramu[l]
        o.tauy_patch[p] = -o.forc_rho_grc[g] * o.forc_v_grc[g] / o.ramu[l]

        # Sensible heat flux from ground using new canopy air temperature
        dth = o.taf_lun[l] - o.t_grnd_col[c]

        if ct == ICOL_ROOF
            o.eflx_sh_grnd_patch[p] = -o.forc_rho_grc[g] * cpair * s.wtus_roof_unscl[l] * dth
        elseif ct == ICOL_ROAD_PERV
            o.eflx_sh_grnd_patch[p] = -o.forc_rho_grc[g] * cpair * s.wtus_road_perv_unscl[l] * dth
        elseif ct == ICOL_ROAD_IMPERV
            o.eflx_sh_grnd_patch[p] = -o.forc_rho_grc[g] * cpair * s.wtus_road_imperv_unscl[l] * dth
        elseif ct == ICOL_SUNWALL
            o.eflx_sh_grnd_patch[p] = -o.forc_rho_grc[g] * cpair * s.wtus_sunwall_unscl[l] * dth
        elseif ct == ICOL_SHADEWALL
            o.eflx_sh_grnd_patch[p] = -o.forc_rho_grc[g] * cpair * s.wtus_shadewall_unscl[l] * dth
        end
        o.eflx_sh_snow_patch[p]   = zero(T)
        o.eflx_sh_soil_patch[p]   = zero(T)
        o.eflx_sh_h2osfc_patch[p] = zero(T)

        o.eflx_sh_tot_patch[p]   = o.eflx_sh_grnd_patch[p]
        o.eflx_sh_tot_u_patch[p] = o.eflx_sh_tot_patch[p]

        # Latent heat flux
        dqh = o.qaf_lun[l] - o.qg_col[c]

        if ct == ICOL_ROOF
            o.qflx_evap_soi_patch[p] = -o.forc_rho_grc[g] * q.wtuq_roof_unscl[l] * dqh
        elseif ct == ICOL_ROAD_PERV
            if dqh > zero(T) || o.frac_sno_col[c] > zero(T) || o.soilalpha_u_col[c] <= zero(T)
                o.qflx_evap_soi_patch[p] = -o.forc_rho_grc[g] * q.wtuq_road_perv_unscl[l] * dqh
                o.qflx_tran_veg_patch[p] = zero(T)
            else
                o.qflx_evap_soi_patch[p] = zero(T)
                o.qflx_tran_veg_patch[p] = -o.forc_rho_grc[g] * q.wtuq_road_perv_unscl[l] * dqh
            end
            o.qflx_evap_veg_patch[p] = o.qflx_tran_veg_patch[p]
        elseif ct == ICOL_ROAD_IMPERV
            o.qflx_evap_soi_patch[p] = -o.forc_rho_grc[g] * q.wtuq_road_imperv_unscl[l] * dqh
        elseif ct == ICOL_SUNWALL
            o.qflx_evap_soi_patch[p] = zero(T)
        elseif ct == ICOL_SHADEWALL
            o.qflx_evap_soi_patch[p] = zero(T)
        end

        # SCALED sensible and latent heat flux for error check
        o.eflx_sh_grnd_scale[p]  = -o.forc_rho_grc[g] * cpair * o.wtus_col[c] * dth
        o.qflx_evap_soi_scale[p] = -o.forc_rho_grc[g] * o.wtuq_col[c] * dqh
    end
end

"""
    urban_fluxes_canyon_surface!(energyflux, frictionvel, waterfluxbulk,
        waterdiagbulk, soilstate, temperature, patch_data, col_data, m_urbanp,
        wts::_UfWtsDV, wtq::_UfWtqDV, ramu, zeta_lunit, wtus_col, wtuq_col,
        forc_rho_grc, forc_u_grc, forc_v_grc) -> nothing

Run the canyon-surface flux pass (per urban patch) WHOLE on the device as a single
kernel, given the SOLVED per-landunit canopy-air state. Writes every patch-level
energyflux/waterflux/frictionvel output of urban_fluxes!' "Determine fluxes from
canyon surfaces" block, plus the scaled-flux scratch used by the error check.
CPU byte-identical on the host path.
"""
function urban_fluxes_canyon_surface!(
        energyflux, frictionvel, waterfluxbulk, waterdiagbulk, soilstate,
        temperature, patch_data, col_data,
        m_urbanp::AbstractVector{Bool},
        wts::_UfWtsDV, wtq::_UfWtqDV,
        ramu, zeta_lunit, wtus_col, wtuq_col,
        eflx_sh_grnd_scale, qflx_evap_soi_scale,
        forc_rho_grc::AbstractVector{<:Real},
        forc_u_grc::AbstractVector{<:Real},
        forc_v_grc::AbstractVector{<:Real})

    o = _UfCanyonDV(
        cgrnds_patch         = energyflux.cgrnds_patch,
        cgrndl_patch         = energyflux.cgrndl_patch,
        taux_patch           = energyflux.taux_patch,
        tauy_patch           = energyflux.tauy_patch,
        ulrad_patch          = energyflux.ulrad_patch,
        dlrad_patch          = energyflux.dlrad_patch,
        eflx_sh_grnd_patch   = energyflux.eflx_sh_grnd_patch,
        eflx_sh_snow_patch   = energyflux.eflx_sh_snow_patch,
        eflx_sh_soil_patch   = energyflux.eflx_sh_soil_patch,
        eflx_sh_h2osfc_patch = energyflux.eflx_sh_h2osfc_patch,
        eflx_sh_tot_patch    = energyflux.eflx_sh_tot_patch,
        eflx_sh_tot_u_patch  = energyflux.eflx_sh_tot_u_patch,
        qflx_evap_soi_patch  = waterfluxbulk.wf.qflx_evap_soi_patch,
        qflx_tran_veg_patch  = waterfluxbulk.wf.qflx_tran_veg_patch,
        qflx_evap_veg_patch  = waterfluxbulk.wf.qflx_evap_veg_patch,
        ram1_patch           = frictionvel.ram1_patch,
        zeta_patch           = frictionvel.zeta_patch,
        eflx_sh_grnd_scale   = eflx_sh_grnd_scale,
        qflx_evap_soi_scale  = qflx_evap_soi_scale,
        ramu                 = ramu,
        zeta_lunit           = zeta_lunit,
        taf_lun              = temperature.taf_lun,
        qaf_lun              = waterdiagbulk.qaf_lun,
        t_grnd_col           = temperature.t_grnd_col,
        qg_col               = waterdiagbulk.qg_col,
        htvp_col             = energyflux.htvp_col,
        dqgdT_col            = waterdiagbulk.dqgdT_col,
        frac_sno_col         = waterdiagbulk.frac_sno_col,
        soilalpha_u_col      = soilstate.soilalpha_u_col,
        wtus_col             = wtus_col,
        wtuq_col             = wtuq_col,
        itype                = col_data.itype,
        forc_rho_grc         = forc_rho_grc,
        forc_u_grc           = forc_u_grc,
        forc_v_grc           = forc_v_grc,
        column               = patch_data.column,
        gridcell             = patch_data.gridcell,
        landunit             = patch_data.landunit)

    FT = eltype(energyflux.cgrnd_patch)
    _launch!(_uf_canyon_surface_kernel!, energyflux.cgrnd_patch, m_urbanp, o, wts, wtq,
             convert(FT, CPAIR); ndrange=length(energyflux.cgrnd_patch))
    return nothing
end

# ---------------------------------------------------------------------------
# Main subroutine: urban_fluxes!
# ---------------------------------------------------------------------------

"""
    urban_fluxes!(...)

Compute turbulent and momentum fluxes from urban canyon (roof, sunwall,
shadewall, pervious and impervious road).

Ported from: UrbanFluxesMod.F90 :: UrbanFluxes
"""
function urban_fluxes!(
        # Data structures
        energyflux       ::EnergyFluxData,
        frictionvel      ::FrictionVelocityData,
        temperature      ::TemperatureData,
        soilstate        ::SoilStateData,
        urbanparams      ::UrbanParamsData,
        waterfluxbulk    ::WaterFluxBulkData,
        waterstatebulk   ::WaterStateBulkData,
        waterdiagbulk    ::WaterDiagnosticBulkData,
        patch_data       ::PatchData,
        col_data         ::ColumnData,
        lun_data         ::LandunitData,
        # Filters
        num_nourbanl     ::Int,
        filter_nourbanl  ::Vector{Int},
        num_urbanl       ::Int,
        filter_urbanl    ::Vector{Int},
        num_urbanc       ::Int,
        filter_urbanc    ::Vector{Int},
        num_urbanp       ::Int,
        filter_urbanp    ::Vector{Int},
        # Bounds
        bounds_lun       ::UnitRange{Int},
        bounds_col       ::UnitRange{Int},
        bounds_patch     ::UnitRange{Int},
        # Atmospheric forcing (gridcell-level)
        forc_t_grc       ::Vector{<:Real},
        forc_th_grc      ::Vector{<:Real},
        forc_rho_grc     ::Vector{<:Real},
        forc_q_grc       ::Vector{<:Real},
        forc_pbot_grc    ::Vector{<:Real},
        forc_u_grc       ::Vector{<:Real},
        forc_v_grc       ::Vector{<:Real};
        # Optional masks for the device-resident element-wise passes. When omitted
        # they are built from the integer filters on the host (CPU byte-identical).
        # The iterative canyon-air solve always runs on the host (CPU state).
        mask_nourbanl    ::Union{AbstractVector{Bool},Nothing} = nothing,
        mask_urbanp      ::Union{AbstractVector{Bool},Nothing} = nothing,
        # Optional time step info
        dtime            ::Real = 3600.0,
        nstep            ::Int     = 1)

    params = urban_fluxes_params
    nlevgrnd = varpar.nlevgrnd

    begl = first(bounds_lun)
    endl = last(bounds_lun)
    begc = first(bounds_col)
    endc = last(bounds_col)
    begp = first(bounds_patch)
    endp = last(bounds_patch)

    lapse_rate = 0.0098   # Dry adiabatic lapse rate (K/m)
    niters = 3            # maximum number of iterations for surface temperature

    # --- Local work arrays (landunit-indexed) ---
    nl = endl
    nc = endc
    np_local = endp

    FT = eltype(forc_t_grc)
    canyontop_wind      = zeros(FT, nl)
    canyon_u_wind       = zeros(FT, nl)
    canyon_wind         = zeros(FT, nl)
    canyon_resistance   = zeros(FT, nl)
    ur                  = zeros(FT, nl)
    ustar_loc           = zeros(FT, nl)
    ramu                = zeros(FT, nl)
    rahu                = zeros(FT, nl)
    rawu                = zeros(FT, nl)
    temp1               = zeros(FT, nl)
    temp12m             = zeros(FT, nl)
    temp2               = zeros(FT, nl)
    temp22m             = zeros(FT, nl)
    thm_g               = zeros(FT, nl)
    thv_g               = zeros(FT, nl)
    dth                 = zeros(FT, nl)
    dqh                 = zeros(FT, nl)
    zldis_arr           = zeros(FT, nl)
    zeta_lunit          = zeros(FT, nl)
    um                  = zeros(FT, nl)
    obu                 = zeros(FT, nl)
    taf_numer           = zeros(FT, nl)
    taf_denom           = zeros(FT, nl)
    qaf_numer           = zeros(FT, nl)
    qaf_denom           = zeros(FT, nl)
    wtas                = zeros(FT, nl)
    wtaq                = zeros(FT, nl)
    wts_sum             = zeros(FT, nl)
    wtq_sum             = zeros(FT, nl)
    beta_arr            = zeros(FT, nl)
    zii_arr             = zeros(FT, nl)
    fm_arr              = zeros(FT, nl)

    wtus_col            = zeros(FT, nc)
    wtuq_col            = zeros(FT, nc)

    wtus_roof           = zeros(FT, nl)
    wtuq_roof           = zeros(FT, nl)
    wtus_road_perv      = zeros(FT, nl)
    wtuq_road_perv      = zeros(FT, nl)
    wtus_road_imperv    = zeros(FT, nl)
    wtuq_road_imperv    = zeros(FT, nl)
    wtus_sunwall        = zeros(FT, nl)
    wtuq_sunwall        = zeros(FT, nl)
    wtus_shadewall      = zeros(FT, nl)
    wtuq_shadewall      = zeros(FT, nl)

    wtus_roof_unscl        = zeros(FT, nl)
    wtuq_roof_unscl        = zeros(FT, nl)
    wtus_road_perv_unscl   = zeros(FT, nl)
    wtuq_road_perv_unscl   = zeros(FT, nl)
    wtus_road_imperv_unscl = zeros(FT, nl)
    wtuq_road_imperv_unscl = zeros(FT, nl)
    wtus_sunwall_unscl     = zeros(FT, nl)
    wtuq_sunwall_unscl     = zeros(FT, nl)
    wtus_shadewall_unscl   = zeros(FT, nl)
    wtuq_shadewall_unscl   = zeros(FT, nl)

    eflx_sh_grnd_scale     = zeros(FT, np_local)
    qflx_evap_soi_scale    = zeros(FT, np_local)

    eflx_wasteheat_roof      = zeros(FT, nl)
    eflx_wasteheat_sunwall   = zeros(FT, nl)
    eflx_wasteheat_shadewall = zeros(FT, nl)
    eflx_heat_from_ac_roof      = zeros(FT, nl)
    eflx_heat_from_ac_sunwall   = zeros(FT, nl)
    eflx_heat_from_ac_shadewall = zeros(FT, nl)

    eflx_arr       = zeros(FT, nl)
    qflx_arr       = zeros(FT, nl)
    eflx_scale_arr = zeros(FT, nl)
    qflx_scale_arr = zeros(FT, nl)
    eflx_err       = zeros(FT, nl)
    qflx_err       = zeros(FT, nl)

    # =========================================================================
    # Set restart fields for non-urban landunits
    # =========================================================================
    # Element-wise passes that write STATE fields run as device kernels (mask-based).
    # Masks default-built from the integer filters on the host (CPU byte-identical).
    m_nourbanl = mask_nourbanl === nothing ?
        _uf_device_mask(temperature.taf_lun, filter_nourbanl, num_nourbanl, nl) :
        mask_nourbanl
    m_urbanp = mask_urbanp === nothing ?
        _uf_device_mask(temperature.t_ref2m_patch, filter_urbanp, num_urbanp, np_local) :
        mask_urbanp

    spval = convert(FT, SPVAL)
    _launch!(_uf_nourban_restart_kernel!, temperature.taf_lun,
             waterdiagbulk.qaf_lun, m_nourbanl, spval)

    # =========================================================================
    # Set constants
    # =========================================================================
    for l in begl:endl
        beta_arr[l] = 1.0
        zii_arr[l]  = 1000.0
    end

    # =========================================================================
    # Compute canyontop wind using Masson (2000)
    # =========================================================================
    for fl in 1:num_urbanl
        l = filter_urbanl[fl]
        g = lun_data.gridcell[l]

        # Error checks
        if lun_data.ht_roof[l] - lun_data.z_d_town[l] <= lun_data.z_0_town[l]
            if _is_ad_type(eltype(frictionvel.forc_hgt_u_patch))
                @warn "urban_fluxes! aerodynamic parameter error (AD mode, continuing)" maxlog=1
            else
                error("aerodynamic parameter error in urban_fluxes!: ht_roof - z_d_town <= z_0_town " *
                      "ht_roof=$(lun_data.ht_roof[l]) z_d_town=$(lun_data.z_d_town[l]) z_0_town=$(lun_data.z_0_town[l])")
            end
        end
        if frictionvel.forc_hgt_u_patch[lun_data.patchi[l]] - lun_data.z_d_town[l] <= lun_data.z_0_town[l]
            if _is_ad_type(eltype(frictionvel.forc_hgt_u_patch))
                @warn "urban_fluxes! aerodynamic parameter error (AD mode, continuing)" maxlog=1
            else
                error("aerodynamic parameter error in urban_fluxes!: forc_hgt_u - z_d_town <= z_0_town " *
                      "forc_hgt_u=$(frictionvel.forc_hgt_u_patch[lun_data.patchi[l]]) z_d_town=$(lun_data.z_d_town[l]) z_0_town=$(lun_data.z_0_town[l])")
            end
        end

        # Magnitude of atmospheric wind
        ur[l] = smooth_max(params.wind_min, sqrt(forc_u_grc[g]^2 + forc_v_grc[g]^2))

        # Canyon top wind
        canyontop_wind[l] = ur[l] *
            log((lun_data.ht_roof[l] - lun_data.z_d_town[l]) / lun_data.z_0_town[l]) /
            log((frictionvel.forc_hgt_u_patch[lun_data.patchi[l]] - lun_data.z_d_town[l]) / lun_data.z_0_town[l])

        # U component of canyon wind
        if lun_data.canyon_hwr[l] < 0.5  # isolated roughness flow
            canyon_u_wind[l] = canyontop_wind[l] * exp(-0.5 * lun_data.canyon_hwr[l] *
                (1.0 - (urbanparams.wind_hgt_canyon[l] / lun_data.ht_roof[l])))
        elseif lun_data.canyon_hwr[l] < 1.0  # wake interference flow
            canyon_u_wind[l] = canyontop_wind[l] * (1.0 + 2.0 * (2.0 / RPI - 1.0) *
                (lun_data.ht_roof[l] / (lun_data.ht_roof[l] / lun_data.canyon_hwr[l]) - 0.5)) *
                exp(-0.5 * lun_data.canyon_hwr[l] * (1.0 - (urbanparams.wind_hgt_canyon[l] / lun_data.ht_roof[l])))
        else  # skimming flow
            canyon_u_wind[l] = canyontop_wind[l] * (2.0 / RPI) *
                exp(-0.5 * lun_data.canyon_hwr[l] * (1.0 - (urbanparams.wind_hgt_canyon[l] / lun_data.ht_roof[l])))
        end
    end

    # =========================================================================
    # Compute fluxes — follows CLM approach for bare soils
    # =========================================================================
    for fl in 1:num_urbanl
        l = filter_urbanl[fl]
        g = lun_data.gridcell[l]

        thm_g[l] = forc_t_grc[g] + lapse_rate * frictionvel.forc_hgt_t_patch[lun_data.patchi[l]]
        thv_g[l] = forc_th_grc[g] * (1.0 + 0.61 * forc_q_grc[g])
        dth[l]   = thm_g[l] - temperature.taf_lun[l]
        dqh[l]   = forc_q_grc[g] - waterdiagbulk.qaf_lun[l]
        dthv     = dth[l] * (1.0 + 0.61 * forc_q_grc[g]) + 0.61 * forc_th_grc[g] * dqh[l]
        zldis_arr[l] = frictionvel.forc_hgt_u_patch[lun_data.patchi[l]] - lun_data.z_d_town[l]

        # Initialize Monin-Obukhov length and wind speed
        (um_val, obu_val) = monin_obuk_ini(frictionvel.zetamaxstable,
            ur[l], thv_g[l], dthv, zldis_arr[l], lun_data.z_0_town[l])
        um[l]  = um_val
        obu[l] = obu_val
    end

    # =========================================================================
    # Stability iteration
    # =========================================================================
    for iter in 1:niters

        # Get friction velocity
        if num_urbanl > 0
            z_d_vec  = lun_data.z_d_town
            z_0_vec  = lun_data.z_0_town

            friction_velocity!(frictionvel, num_urbanl, filter_urbanl,
                z_d_vec, z_0_vec, z_0_vec, z_0_vec,
                obu, iter, ur, um, ustar_loc,
                temp1, temp2, temp12m, temp22m, fm_arr;
                landunit_index=true,
                lun_gridcell=lun_data.gridcell,
                lun_patchi=lun_data.patchi,
                lun_patchf=lun_data.patchf)
        end

        for fl in 1:num_urbanl
            l = filter_urbanl[fl]
            g = lun_data.gridcell[l]

            # Aerodynamic resistance from urban canopy air to atmosphere
            ramu[l] = 1.0 / (ustar_loc[l]^2 / um[l])
            rahu[l] = 1.0 / (temp1[l] * ustar_loc[l])
            rawu[l] = 1.0 / (temp2[l] * ustar_loc[l])

            # Canyon wind magnitude
            canyon_wind[l] = sqrt(canyon_u_wind[l]^2.0 + ustar_loc[l]^2.0)

            # Canyon resistance (Masson 2000)
            canyon_resistance[l] = CPAIR * forc_rho_grc[g] / (11.8 + 4.2 * canyon_wind[l])
        end

        # First term in taf/qaf equations
        for fl in 1:num_urbanl
            l = filter_urbanl[fl]
            g = lun_data.gridcell[l]

            taf_numer[l] = thm_g[l] / rahu[l]
            taf_denom[l] = 1.0 / rahu[l]
            qaf_numer[l] = forc_q_grc[g] / rawu[l]
            qaf_denom[l] = 1.0 / rawu[l]

            wtas[l] = 1.0 / rahu[l]
            wtaq[l] = 1.0 / rawu[l]
        end

        # Gather terms from urban columns
        for fc in 1:num_urbanc
            c = filter_urbanc[fc]
            l = col_data.landunit[c]

            if col_data.itype[c] == ICOL_ROOF
                # Scaled sensible heat conductance
                wtus_col[c] = lun_data.wtlunit_roof[l] / canyon_resistance[l]
                wtus_roof[l] = wtus_col[c]
                wtus_roof_unscl[l] = 1.0 / canyon_resistance[l]

                # Wetness fraction for roof
                if waterdiagbulk.snow_depth_col[c] > 0.0
                    fwet_roof = smooth_min(waterdiagbulk.snow_depth_col[c] / 0.05, 1.0)
                else
                    fwet_roof = (smooth_max(0.0, waterstatebulk.ws.h2osoi_liq_col[c, 1] +
                        waterstatebulk.ws.h2osoi_ice_col[c, 1]) / PONDMX_URBAN)^0.666666666666
                    fwet_roof = smooth_min(fwet_roof, 1.0)
                end
                if waterdiagbulk.qaf_lun[l] > waterdiagbulk.qg_col[c]
                    fwet_roof = 1.0
                end

                # Scaled latent heat conductance
                wtuq_col[c] = fwet_roof * (lun_data.wtlunit_roof[l] / canyon_resistance[l])
                wtuq_roof[l] = wtuq_col[c]
                wtuq_roof_unscl[l] = fwet_roof * (1.0 / canyon_resistance[l])

                if is_simple_build_temp()
                    (wh, hfac) = simple_wasteheatfromac(
                        energyflux.eflx_urban_ac_col[c], energyflux.eflx_urban_heat_col[c])
                    eflx_wasteheat_roof[l] = wh
                    eflx_heat_from_ac_roof[l] = hfac
                end

            elseif col_data.itype[c] == ICOL_ROAD_PERV
                # Scaled sensible heat conductance
                wtus_col[c] = lun_data.wtroad_perv[l] * (1.0 - lun_data.wtlunit_roof[l]) / canyon_resistance[l]
                wtus_road_perv[l] = wtus_col[c]
                wtus_road_perv_unscl[l] = 1.0 / canyon_resistance[l]

                # Scaled latent heat conductance
                wtuq_col[c] = lun_data.wtroad_perv[l] * (1.0 - lun_data.wtlunit_roof[l]) / canyon_resistance[l]
                wtuq_road_perv[l] = wtuq_col[c]
                wtuq_road_perv_unscl[l] = 1.0 / canyon_resistance[l]

            elseif col_data.itype[c] == ICOL_ROAD_IMPERV
                # Scaled sensible heat conductance
                wtus_col[c] = (1.0 - lun_data.wtroad_perv[l]) * (1.0 - lun_data.wtlunit_roof[l]) / canyon_resistance[l]
                wtus_road_imperv[l] = wtus_col[c]
                wtus_road_imperv_unscl[l] = 1.0 / canyon_resistance[l]

                # Wetness fraction for impervious road
                if waterdiagbulk.snow_depth_col[c] > 0.0
                    fwet_road_imperv = smooth_min(waterdiagbulk.snow_depth_col[c] / 0.05, 1.0)
                else
                    fwet_road_imperv = (smooth_max(0.0, waterstatebulk.ws.h2osoi_liq_col[c, 1] +
                        waterstatebulk.ws.h2osoi_ice_col[c, 1]) / PONDMX_URBAN)^0.666666666666
                    fwet_road_imperv = smooth_min(fwet_road_imperv, 1.0)
                end
                if waterdiagbulk.qaf_lun[l] > waterdiagbulk.qg_col[c]
                    fwet_road_imperv = 1.0
                end

                # Scaled latent heat conductance
                wtuq_col[c] = fwet_road_imperv * (1.0 - lun_data.wtroad_perv[l]) *
                    (1.0 - lun_data.wtlunit_roof[l]) / canyon_resistance[l]
                wtuq_road_imperv[l] = wtuq_col[c]
                wtuq_road_imperv_unscl[l] = fwet_road_imperv * (1.0 / canyon_resistance[l])

            elseif col_data.itype[c] == ICOL_SUNWALL
                # Scaled sensible heat conductance
                wtus_col[c] = lun_data.canyon_hwr[l] * (1.0 - lun_data.wtlunit_roof[l]) / canyon_resistance[l]
                wtus_sunwall[l] = wtus_col[c]
                wtus_sunwall_unscl[l] = 1.0 / canyon_resistance[l]

                # Walls have zero latent heat conductance
                wtuq_col[c] = 0.0
                wtuq_sunwall[l] = 0.0
                wtuq_sunwall_unscl[l] = 0.0

                if is_simple_build_temp()
                    (wh, hfac) = simple_wasteheatfromac(
                        energyflux.eflx_urban_ac_col[c], energyflux.eflx_urban_heat_col[c])
                    eflx_wasteheat_sunwall[l] = wh
                    eflx_heat_from_ac_sunwall[l] = hfac
                end

            elseif col_data.itype[c] == ICOL_SHADEWALL
                # Scaled sensible heat conductance
                wtus_col[c] = lun_data.canyon_hwr[l] * (1.0 - lun_data.wtlunit_roof[l]) / canyon_resistance[l]
                wtus_shadewall[l] = wtus_col[c]
                wtus_shadewall_unscl[l] = 1.0 / canyon_resistance[l]

                # Walls have zero latent heat conductance
                wtuq_col[c] = 0.0
                wtuq_shadewall[l] = 0.0
                wtuq_shadewall_unscl[l] = 0.0

                if is_simple_build_temp()
                    (wh, hfac) = simple_wasteheatfromac(
                        energyflux.eflx_urban_ac_col[c], energyflux.eflx_urban_heat_col[c])
                    eflx_wasteheat_shadewall[l] = wh
                    eflx_heat_from_ac_shadewall[l] = hfac
                end

            else
                if _is_ad_type(eltype(energyflux.eflx_sh_tot_patch))
                    @warn "urban_fluxes! ctype out of range (AD mode, continuing)" maxlog=1
                else
                    error("ERROR: ctype out of range in urban_fluxes! c=$c ctype=$(col_data.itype[c])")
                end
            end

            taf_numer[l] = taf_numer[l] + temperature.t_grnd_col[c] * wtus_col[c]
            taf_denom[l] = taf_denom[l] + wtus_col[c]
            qaf_numer[l] = qaf_numer[l] + waterdiagbulk.qg_col[c] * wtuq_col[c]
            qaf_denom[l] = qaf_denom[l] + wtuq_col[c]
        end

        # Calculate new urban canopy air temperature and specific humidity
        wasteheat!(energyflux, lun_data, num_urbanl, filter_urbanl,
            eflx_wasteheat_roof, eflx_wasteheat_sunwall, eflx_wasteheat_shadewall,
            eflx_heat_from_ac_roof, eflx_heat_from_ac_sunwall, eflx_heat_from_ac_shadewall)

        for fl in 1:num_urbanl
            l = filter_urbanl[fl]
            g = lun_data.gridcell[l]

            # Calculate traffic heat flux (only from impervious road)
            energyflux.eflx_traffic_lun[l] = (1.0 - lun_data.wtlunit_roof[l]) *
                (1.0 - lun_data.wtroad_perv[l]) * urbanparams.eflx_traffic_factor[l]

            temperature.taf_lun[l] = taf_numer[l] / taf_denom[l]
            waterdiagbulk.qaf_lun[l] = qaf_numer[l] / qaf_denom[l]

            wts_sum[l] = wtas[l] + wtus_roof[l] + wtus_road_perv[l] +
                wtus_road_imperv[l] + wtus_sunwall[l] + wtus_shadewall[l]

            wtq_sum[l] = wtaq[l] + wtuq_roof[l] + wtuq_road_perv[l] +
                wtuq_road_imperv[l] + wtuq_sunwall[l] + wtuq_shadewall[l]
        end

        # Determine stability using new taf and qaf
        for fl in 1:num_urbanl
            l = filter_urbanl[fl]
            g = lun_data.gridcell[l]

            dth[l] = thm_g[l] - temperature.taf_lun[l]
            dqh[l] = forc_q_grc[g] - waterdiagbulk.qaf_lun[l]
            tstar = temp1[l] * dth[l]
            qstar = temp2[l] * dqh[l]
            thvstar = tstar * (1.0 + 0.61 * forc_q_grc[g]) + 0.61 * forc_th_grc[g] * qstar
            zeta_lunit[l] = zldis_arr[l] * VKC * GRAV * thvstar / (ustar_loc[l]^2 * thv_g[l])

            if zeta_lunit[l] >= 0.0  # stable
                zeta_lunit[l] = smooth_min(frictionvel.zetamaxstable, smooth_max(zeta_lunit[l], 0.01))
                um[l] = smooth_max(ur[l], 0.1)
            else  # unstable
                zeta_lunit[l] = smooth_max(-100.0, smooth_min(zeta_lunit[l], -0.01))
                wc = beta_arr[l] * (-GRAV * ustar_loc[l] * thvstar * zii_arr[l] / thv_g[l])^0.333
                um[l] = sqrt(ur[l]^2 + wc^2)
            end

            obu[l] = zldis_arr[l] / zeta_lunit[l]
        end
    end  # end stability iteration

    # =========================================================================
    # Determine fluxes from canyon surfaces
    # =========================================================================
    # Per-patch + column-INDEPENDENT given the SOLVED canopy-air state: runs WHOLE
    # on the device as a single kernel (see urban_fluxes_canyon_surface!). The
    # per-landunit work arrays are host-local here; _uf_to_backend moves them onto
    # the state structs' backend (identity / byte-identical on Float64 CPU). The
    # scaled-flux scratch is initialised to zero on-device via fill! over the full
    # patch range (matches the original begp:endp init), then overwritten for urban
    # patches by the kernel.
    proto_pp = energyflux.cgrnd_patch
    fill!(eflx_sh_grnd_scale, zero(eltype(eflx_sh_grnd_scale)))
    fill!(qflx_evap_soi_scale, zero(eltype(qflx_evap_soi_scale)))

    b(v) = _uf_to_backend(proto_pp, v)
    wts_dv = _UfWtsDV(
        wtas = b(wtas),
        wtus_roof = b(wtus_roof), wtus_road_perv = b(wtus_road_perv),
        wtus_road_imperv = b(wtus_road_imperv), wtus_sunwall = b(wtus_sunwall),
        wtus_shadewall = b(wtus_shadewall),
        wtus_roof_unscl = b(wtus_roof_unscl), wtus_road_perv_unscl = b(wtus_road_perv_unscl),
        wtus_road_imperv_unscl = b(wtus_road_imperv_unscl), wtus_sunwall_unscl = b(wtus_sunwall_unscl),
        wtus_shadewall_unscl = b(wtus_shadewall_unscl), wts_sum = b(wts_sum))
    wtq_dv = _UfWtqDV(
        wtaq = b(wtaq),
        wtuq_roof = b(wtuq_roof), wtuq_road_perv = b(wtuq_road_perv),
        wtuq_road_imperv = b(wtuq_road_imperv), wtuq_sunwall = b(wtuq_sunwall),
        wtuq_shadewall = b(wtuq_shadewall),
        wtuq_roof_unscl = b(wtuq_roof_unscl), wtuq_road_perv_unscl = b(wtuq_road_perv_unscl),
        wtuq_road_imperv_unscl = b(wtuq_road_imperv_unscl), wtq_sum = b(wtq_sum))

    urban_fluxes_canyon_surface!(
        energyflux, frictionvel, waterfluxbulk, waterdiagbulk, soilstate,
        temperature, patch_data, col_data, m_urbanp, wts_dv, wtq_dv,
        b(ramu), b(zeta_lunit), b(wtus_col), b(wtuq_col),
        eflx_sh_grnd_scale, qflx_evap_soi_scale,
        b(forc_rho_grc), b(forc_u_grc), b(forc_v_grc))

    # Host copies of the scaled-flux scratch for the sum()-based error check below
    # (the check is a reduction over each landunit's patch range and stays HOST).
    eflx_sh_grnd_scale = collect(eflx_sh_grnd_scale)
    qflx_evap_soi_scale = collect(qflx_evap_soi_scale)

    # =========================================================================
    # Error checking: total fluxes should equal sum of scaled fluxes
    # =========================================================================
    for fl in 1:num_urbanl
        l = filter_urbanl[fl]
        g = lun_data.gridcell[l]
        eflx_arr[l] = -(forc_rho_grc[g] * CPAIR / rahu[l]) * (thm_g[l] - temperature.taf_lun[l])
        qflx_arr[l] = -(forc_rho_grc[g] / rawu[l]) * (forc_q_grc[g] - waterdiagbulk.qaf_lun[l])
        eflx_scale_arr[l] = sum(eflx_sh_grnd_scale[lun_data.patchi[l]:lun_data.patchf[l]])
        qflx_scale_arr[l] = sum(qflx_evap_soi_scale[lun_data.patchi[l]:lun_data.patchf[l]])
        eflx_err[l] = eflx_scale_arr[l] - eflx_arr[l]
        qflx_err[l] = qflx_scale_arr[l] - qflx_arr[l]
    end

    # Check sensible heat flux error
    for fl in 1:num_urbanl
        l = filter_urbanl[fl]
        if abs(eflx_err[l]) > 0.01
            @warn "Total sensible heat does not equal sum of scaled heat fluxes for urban columns" nstep l eflx_err=eflx_err[l]
            if abs(eflx_err[l]) > 0.01
                if _is_ad_type(eltype(eflx_err))
                    @warn "urban_fluxes! sensible heat flux error (AD mode, continuing)" maxlog=1
                else
                    error("urban_fluxes! sensible heat flux error > 0.01 W/m**2: " *
                          "eflx_scale=$(eflx_scale_arr[l]) eflx=$(eflx_arr[l]) err=$(eflx_err[l])")
                end
            end
        end
    end

    # Check water vapor flux error
    for fl in 1:num_urbanl
        l = filter_urbanl[fl]
        if abs(qflx_err[l]) > 4.0e-9
            @warn "Total water vapor flux does not equal sum of scaled fluxes for urban columns" nstep l qflx_err=qflx_err[l]
            if abs(qflx_err[l]) > 4.0e-9
                if _is_ad_type(eltype(qflx_err))
                    @warn "urban_fluxes! water vapor flux error (AD mode, continuing)" maxlog=1
                else
                    error("urban_fluxes! water vapor flux error > 4e-9 kg/m**2/s: " *
                          "qflx_scale=$(qflx_scale_arr[l]) qflx=$(qflx_arr[l]) err=$(qflx_err[l])")
                end
            end
        end
    end

    # =========================================================================
    # Internal building temperature (simple method)
    # =========================================================================
    if is_simple_build_temp()
        calc_simple_internal_building_temp!(temperature, col_data, lun_data,
            num_urbanc, filter_urbanc, num_urbanl, filter_urbanl)
    end

    # =========================================================================
    # Roots for urban (pervious-road only) + 2-m temperature/humidity/RH/t_veg
    # diagnostics — both column-independent device kernels. Re-runs the (idempotent)
    # non-urban taf/qaf SPVAL restart too; the non-urban landunits are untouched by
    # the host canyon-air solve, so this is byte-identical to setting it once.
    # =========================================================================
    urban_fluxes_diagnostics!(temperature, waterdiagbulk, soilstate, patch_data,
        col_data, m_nourbanl, m_urbanp, forc_t_grc, forc_pbot_grc, nlevgrnd)

    return nothing
end
