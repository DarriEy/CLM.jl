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
# Device-capable internal building temperature — ONE THREAD PER URBAN LANDUNIT.
# Each thread gathers its member columns' (coli:colf) inner-layer temperatures
# (t_soisno at nlevurb) for roof/sunwall/shadewall, then computes t_building_lun.
# Runs WHOLE on the device; CPU byte-identical to calc_simple_internal_building_temp!.
# (urban_fluxes! always has coli/colf set; the standalone scalar function above is
# kept for callers that build only col_data.landunit.)
# ---------------------------------------------------------------------------
@kernel function _uf_building_temp_kernel!(t_building_lun, @Const(mask),
        @Const(t_soisno_col), @Const(col_itype), @Const(coli), @Const(colf),
        @Const(ht_roof), @Const(canyon_hwr), @Const(wtlunit_roof), nlevurb::Int)
    fl = @index(Global)
    @inbounds if mask[fl]
        l = fl
        T = eltype(t_building_lun)
        t_roof = zero(T); t_sun = zero(T); t_sha = zero(T)
        for c in coli[l]:colf[l]
            ct = col_itype[c]
            if ct == ICOL_ROOF
                t_roof = t_soisno_col[c, nlevurb]
            elseif ct == ICOL_SUNWALL
                t_sun = t_soisno_col[c, nlevurb]
            elseif ct == ICOL_SHADEWALL
                t_sha = t_soisno_col[c, nlevurb]
            end
        end
        ht = ht_roof[l]
        lngth_roof = (ht / canyon_hwr[l]) * wtlunit_roof[l] / (one(T) - wtlunit_roof[l])
        t_building_lun[l] = (ht * (t_sha + t_sun) + lngth_roof * t_roof) /
            (T(2.0) * ht + lngth_roof)
    end
end

function calc_simple_internal_building_temp_dev!(temperature, col_data, lun_data,
        m_urbanl::AbstractVector{Bool})
    nlevurb = varpar.nlevurb
    _launch!(_uf_building_temp_kernel!, temperature.t_building_lun, m_urbanl,
        temperature.t_soisno_col, col_data.itype, lun_data.coli, lun_data.colf,
        lun_data.ht_roof, lun_data.canyon_hwr, lun_data.wtlunit_roof, nlevurb;
        ndrange = length(temperature.t_building_lun))
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
# ITERATIVE CANYON-AIR SOLVE — the last HOST piece, now a single device kernel.
#
# ONE THREAD PER URBAN LANDUNIT. Each thread runs the whole solve for its
# landunit in-thread: canyontop-wind setup (Masson 2000), monin_obuk_ini init,
# the niters=3 stability iteration loop with the FrictionVelocity landunit-path
# math INLINED (urban uses z0m=z0h=z0q=z_0_town and displa=z_d_town), the gather
# over the landunit's member COLUMNS (coli[l]:colf[l], contiguous — same FP
# accumulation order as the host filter_urbanc loop) into taf/qaf, the per-member
# PATCH friction-velocity side-effects (u10/u10_clm/va/fv/vds, over patchi:patchf),
# wasteheat / heat-from-AC / traffic, and the final stability update.
#
# Because each landunit owns its own accumulation, taf/qaf and every reduction are
# in-thread scalars — NO cross-thread reduction is needed. The aerodynamic-parameter
# guards (which `error()`/`@warn`) stay on the host (validation only, no state).
#
# Literals are eltype-converted (T(...)/zero/one) so no Float64 reaches a
# Float32-only backend; on Float64 CPU this is byte-identical to the scalar loops.
# Feature flags (simple/prog build temp, urban_hac) and the AC/HT wasteheat
# factors are module globals -> resolved on the host and passed as scalars.
# ---------------------------------------------------------------------------

# Per-landunit work + lun geometry + gridcell forcing (everything the iteration
# reads/writes at landunit granularity). Split into two bundles to stay well under
# Metal's ~31-arg kernel cap.
Base.@kwdef struct _UfIterLunDV{V,VI}
    # outputs consumed downstream (canyon-surface kernel + error check)
    ramu::V
    rahu::V
    rawu::V
    zeta_lunit::V
    taf_lun::V
    qaf_lun::V
    wtas::V
    wtaq::V
    wts_sum::V
    wtq_sum::V
    wtus_roof::V
    wtuq_roof::V
    wtus_road_perv::V
    wtuq_road_perv::V
    wtus_road_imperv::V
    wtuq_road_imperv::V
    wtus_sunwall::V
    wtuq_sunwall::V
    wtus_shadewall::V
    wtuq_shadewall::V
    wtus_roof_unscl::V
    wtuq_roof_unscl::V
    wtus_road_perv_unscl::V
    wtuq_road_perv_unscl::V
    wtus_road_imperv_unscl::V
    wtuq_road_imperv_unscl::V
    wtus_sunwall_unscl::V
    wtus_shadewall_unscl::V
    # energyflux landunit-level outputs
    eflx_wasteheat_lun::V
    eflx_heat_from_ac_lun::V
    eflx_traffic_lun::V
    # lun geometry / indexing
    lun_gridcell::VI
    coli::VI
    colf::VI
    patchi::VI
    patchf::VI
    ht_roof::V
    z_d_town::V
    z_0_town::V
    canyon_hwr::V
    wtlunit_roof::V
    wtroad_perv::V
    wind_hgt_canyon::V
    eflx_traffic_factor::V
    # gridcell forcing
    forc_u_grc::V
    forc_v_grc::V
    forc_t_grc::V
    forc_th_grc::V
    forc_q_grc::V
    forc_rho_grc::V
    # patch forcing heights (indexed by lun.patchi[l])
    forc_hgt_u_patch::V
    forc_hgt_t_patch::V
end
Adapt.@adapt_structure _UfIterLunDV

# Per-column + per-member-patch state the iteration reads/writes, plus the friction-
# velocity per-patch side-effect outputs.
Base.@kwdef struct _UfIterColDV{V,M,VI}
    # per-column outputs / state
    wtus_col::V
    wtuq_col::V
    t_grnd_col::V
    qg_col::V
    snow_depth_col::V
    h2osoi_liq_col::M   # (column, layer) — only layer 1 is read
    h2osoi_ice_col::M
    eflx_urban_ac_col::V
    eflx_urban_heat_col::V
    col_itype::VI
    # friction-velocity per-PATCH side effects (written for every member patch)
    u10_clm_patch::V
    va_patch::V
    u10_patch::V
    fv_patch::V
    vds_patch::V
end
Adapt.@adapt_structure _UfIterColDV

# The whole iterative canyon-air solve for one urban landunit.
# `taf_lun` is the standalone first arg (drives backend detection in _launch!).
@kernel function _uf_iter_kernel!(taf_lun_out, @Const(mask),
        o::_UfIterLunDV, c_::_UfIterColDV,
        niters::Int, zetamaxstable, wind_min,
        VKC_, GRAV_, CPAIR_, RPI_, PONDMX_URBAN_, lapse_rate_,
        is_simple::Bool, urban_hac_on::Bool, wasteheat_on::Bool,
        ac_wh_factor, ht_wh_factor, wasteheat_limit)
    fl = @index(Global)
    @inbounds if mask[fl]
        l = fl
        T = eltype(taf_lun_out)
        vkc   = T(VKC_)
        grav  = T(GRAV_)
        cpair = T(CPAIR_)
        rpi   = T(RPI_)
        pondmx = T(PONDMX_URBAN_)
        lapse = T(lapse_rate_)
        zetamax = T(zetamaxstable)
        wmin  = T(wind_min)
        zetam = T(1.574)
        zetat = T(0.465)

        g  = o.lun_gridcell[l]
        pi_ = o.patchi[l]                # representative patch (forcing heights)
        ci  = o.coli[l]; cf = o.colf[l]
        ppi = o.patchi[l]; ppf = o.patchf[l]

        z_d   = o.z_d_town[l]
        z_0   = o.z_0_town[l]
        ht    = o.ht_roof[l]
        hwr   = o.canyon_hwr[l]
        wtroof = o.wtlunit_roof[l]
        wtrp   = o.wtroad_perv[l]

        # --- canyontop wind (Masson 2000) ---
        ur = max(wmin, sqrt(o.forc_u_grc[g]^2 + o.forc_v_grc[g]^2))
        canyontop_wind = ur *
            log((ht - z_d) / z_0) /
            log((o.forc_hgt_u_patch[pi_] - z_d) / z_0)

        whc = o.wind_hgt_canyon[l]
        if hwr < T(0.5)
            canyon_u_wind = canyontop_wind * exp(-T(0.5) * hwr * (one(T) - (whc / ht)))
        elseif hwr < one(T)
            canyon_u_wind = canyontop_wind * (one(T) + T(2.0) * (T(2.0) / rpi - one(T)) *
                (ht / (ht / hwr) - T(0.5))) *
                exp(-T(0.5) * hwr * (one(T) - (whc / ht)))
        else
            canyon_u_wind = canyontop_wind * (T(2.0) / rpi) *
                exp(-T(0.5) * hwr * (one(T) - (whc / ht)))
        end

        # --- initialize Monin-Obukhov length and wind speed ---
        thm_g = o.forc_t_grc[g] + lapse * o.forc_hgt_t_patch[pi_]
        thv_g = o.forc_th_grc[g] * (one(T) + T(0.61) * o.forc_q_grc[g])
        dth = thm_g - taf_lun_out[l]
        dqh = o.forc_q_grc[g] - o.qaf_lun[l]
        dthv = dth * (one(T) + T(0.61) * o.forc_q_grc[g]) + T(0.61) * o.forc_th_grc[g] * dqh
        zldis = o.forc_hgt_u_patch[pi_] - z_d

        # monin_obuk_ini (inlined; unused local `ustar=0.06` dropped)
        wc = T(0.5)
        um = dthv >= zero(T) ? max(ur, T(0.1)) : sqrt(ur * ur + wc * wc)
        rib = grav * zldis * dthv / (thv_g * um * um)
        if rib >= zero(T)
            zeta = rib * log(zldis / z_0) / (one(T) - T(5.0) * min(rib, T(0.19)))
            zeta = min(zetamax, max(zeta, T(0.01)))
        else
            zeta = rib * log(zldis / z_0)
            zeta = max(T(-100.0), min(zeta, T(-0.01)))
        end
        obu = zldis / zeta

        # scratch carried across iterations
        ustar = zero(T); temp1 = zero(T); temp2 = zero(T)
        temp12m = zero(T); temp22m = zero(T); fm = zero(T)
        rahu = zero(T); rawu = zero(T)
        canyon_resistance = zero(T)

        # ============== stability iteration ==============
        for iter in 1:niters
            # ---- FrictionVelocity (landunit path), z0m=z0h=z0q=z_0, displa=z_d ----
            # Wind profile -> ustar
            zldisw = o.forc_hgt_u_patch[pi_] - z_d
            zeta_w = zldisw / obu
            if zeta_w < -zetam
                ustar = vkc * um / (log(-zetam * obu / z_0) -
                    stability_func1(-zetam) + stability_func1(z_0 / obu) +
                    T(1.14) * ((-zeta_w)^T(0.333) - (zetam)^T(0.333)))
            elseif zeta_w < zero(T)
                ustar = vkc * um / (log(zldisw / z_0) -
                    stability_func1(zeta_w) + stability_func1(z_0 / obu))
            elseif zeta_w <= one(T)
                ustar = vkc * um / (log(zldisw / z_0) + T(5.0) * zeta_w - T(5.0) * z_0 / obu)
            else
                ustar = vkc * um / (log(obu / z_0) + T(5.0) - T(5.0) * z_0 / obu +
                    (T(5.0) * log(zeta_w) + zeta_w - one(T)))
            end

            # Deposition velocity (per member patch)
            if zeta_w < zero(T)
                vds_tmp = T(2.0e-3) * ustar * (one(T) + (T(300.0) / max(-obu, T(1.0e-10)))^T(0.666))
            else
                vds_tmp = T(2.0e-3) * ustar
            end
            for pp in ppi:ppf
                c_.vds_patch[pp] = vds_tmp
            end

            # 10-m wind (CLM) per member patch
            for pp in ppi:ppf
                if zldisw - z_0 <= T(10.0)
                    c_.u10_clm_patch[pp] = um
                else
                    if zeta_w < -zetam
                        c_.u10_clm_patch[pp] = um - (ustar / vkc * (log(-zetam * obu / (T(10.0) + z_0)) -
                            stability_func1(-zetam) + stability_func1((T(10.0) + z_0) / obu) +
                            T(1.14) * ((-zeta_w)^T(0.333) - (zetam)^T(0.333))))
                    elseif zeta_w < zero(T)
                        c_.u10_clm_patch[pp] = um - (ustar / vkc * (log(zldisw / (T(10.0) + z_0)) -
                            stability_func1(zeta_w) + stability_func1((T(10.0) + z_0) / obu)))
                    elseif zeta_w <= one(T)
                        c_.u10_clm_patch[pp] = um - (ustar / vkc * (log(zldisw / (T(10.0) + z_0)) +
                            T(5.0) * zeta_w - T(5.0) * (T(10.0) + z_0) / obu))
                    else
                        c_.u10_clm_patch[pp] = um - (ustar / vkc * (log(obu / (T(10.0) + z_0)) +
                            T(5.0) - T(5.0) * (T(10.0) + z_0) / obu +
                            (T(5.0) * log(zeta_w) + zeta_w - one(T))))
                    end
                end
                c_.va_patch[pp] = um
            end

            # Temperature profile -> temp1
            zldist = o.forc_hgt_t_patch[pi_] - z_d
            zeta_t = zldist / obu
            if zeta_t < -zetat
                temp1 = vkc / (log(-zetat * obu / z_0) -
                    stability_func2(-zetat) + stability_func2(z_0 / obu) +
                    T(0.8) * ((zetat)^(-T(0.333)) - (-zeta_t)^(-T(0.333))))
            elseif zeta_t < zero(T)
                temp1 = vkc / (log(zldist / z_0) -
                    stability_func2(zeta_t) + stability_func2(z_0 / obu))
            elseif zeta_t <= one(T)
                temp1 = vkc / (log(zldist / z_0) + T(5.0) * zeta_t - T(5.0) * z_0 / obu)
            else
                temp1 = vkc / (log(obu / z_0) + T(5.0) - T(5.0) * z_0 / obu +
                    (T(5.0) * log(zeta_t) + zeta_t - one(T)))
            end

            # Humidity profile: hgt_q==hgt_t and z0q==z0h -> temp2 = temp1
            temp2 = temp1

            # Temperature profile at 2m -> temp12m
            zldis2 = T(2.0) + z_0
            zeta2 = zldis2 / obu
            if zeta2 < -zetat
                temp12m = vkc / (log(-zetat * obu / z_0) -
                    stability_func2(-zetat) + stability_func2(z_0 / obu) +
                    T(0.8) * ((zetat)^(-T(0.333)) - (-zeta2)^(-T(0.333))))
            elseif zeta2 < zero(T)
                temp12m = vkc / (log(zldis2 / z_0) -
                    stability_func2(zeta2) + stability_func2(z_0 / obu))
            elseif zeta2 <= one(T)
                temp12m = vkc / (log(zldis2 / z_0) + T(5.0) * zeta2 - T(5.0) * z_0 / obu)
            else
                temp12m = vkc / (log(obu / z_0) + T(5.0) - T(5.0) * z_0 / obu +
                    (T(5.0) * log(zeta2) + zeta2 - one(T)))
            end
            # Humidity profile at 2m: z0q==z0h -> temp22m = temp12m
            temp22m = temp12m

            # 10-m wind for dust model (per member patch), needs fm
            zldisd = o.forc_hgt_u_patch[pi_] - z_d
            zeta_d = zldisd / obu
            if min(zeta_d, one(T)) < zero(T)
                tmp1 = (one(T) - T(16.0) * min(zeta_d, one(T)))^T(0.25)
                tmp2 = log((one(T) + tmp1 * tmp1) / T(2.0))
                tmp3 = log((one(T) + tmp1) / T(2.0))
                fmnew = T(2.0) * tmp3 + tmp2 - T(2.0) * atan(tmp1) + T(1.5707963)
            else
                fmnew = -T(5.0) * min(zeta_d, one(T))
            end
            if iter == 1
                fm = fmnew
            else
                fm = T(0.5) * (fm + fmnew)
            end

            zeta10 = min(T(10.0) / obu, one(T))
            if zeta_d == zero(T)
                zeta10 = zero(T)
            end
            if zeta10 < zero(T)
                tmp1 = (one(T) - T(16.0) * zeta10)^T(0.25)
                tmp2 = log((one(T) + tmp1 * tmp1) / T(2.0))
                tmp3 = log((one(T) + tmp1) / T(2.0))
                fm10 = T(2.0) * tmp3 + tmp2 - T(2.0) * atan(tmp1) + T(1.5707963)
            else
                fm10 = -T(5.0) * zeta10
            end
            tmp4 = log(max(one(T), o.forc_hgt_u_patch[pi_] / T(10.0)))
            for pp in ppi:ppf
                c_.u10_patch[pp] = ur - ustar / vkc * (tmp4 - fm + fm10)
                c_.fv_patch[pp] = ustar
            end

            # ---- aerodynamic resistances ----
            ramu = one(T) / (ustar^2 / um)
            rahu = one(T) / (temp1 * ustar)
            rawu = one(T) / (temp2 * ustar)

            canyon_wind = sqrt(canyon_u_wind^T(2.0) + ustar^T(2.0))
            canyon_resistance = cpair * o.forc_rho_grc[g] / (T(11.8) + T(4.2) * canyon_wind)

            # ---- first term in taf/qaf equations ----
            taf_numer = thm_g / rahu
            taf_denom = one(T) / rahu
            qaf_numer = o.forc_q_grc[g] / rawu
            qaf_denom = one(T) / rawu
            wtas = one(T) / rahu
            wtaq = one(T) / rawu

            # reset per-iteration the per-landunit weight accumulators
            wtus_roof = zero(T); wtuq_roof = zero(T)
            wtus_road_perv = zero(T); wtuq_road_perv = zero(T)
            wtus_road_imperv = zero(T); wtuq_road_imperv = zero(T)
            wtus_sunwall = zero(T); wtuq_sunwall = zero(T)
            wtus_shadewall = zero(T); wtuq_shadewall = zero(T)
            wtus_roof_unscl = zero(T); wtuq_roof_unscl = zero(T)
            wtus_road_perv_unscl = zero(T); wtuq_road_perv_unscl = zero(T)
            wtus_road_imperv_unscl = zero(T); wtuq_road_imperv_unscl = zero(T)
            wtus_sunwall_unscl = zero(T)
            wtus_shadewall_unscl = zero(T)

            # ---- gather over member columns (coli:colf = column-index order) ----
            for c in ci:cf
                ct = c_.col_itype[c]
                wtus_c = zero(T); wtuq_c = zero(T)
                if ct == ICOL_ROOF
                    wtus_c = wtroof / canyon_resistance
                    wtus_roof = wtus_c
                    wtus_roof_unscl = one(T) / canyon_resistance
                    if c_.snow_depth_col[c] > zero(T)
                        fwet = min(c_.snow_depth_col[c] / T(0.05), one(T))
                    else
                        fwet = (max(zero(T), c_.h2osoi_liq_col[c, 1] + c_.h2osoi_ice_col[c, 1]) / pondmx)^T(0.666666666666)
                        fwet = min(fwet, one(T))
                    end
                    if o.qaf_lun[l] > c_.qg_col[c]
                        fwet = one(T)
                    end
                    wtuq_c = fwet * (wtroof / canyon_resistance)
                    wtuq_roof = wtuq_c
                    wtuq_roof_unscl = fwet * (one(T) / canyon_resistance)
                elseif ct == ICOL_ROAD_PERV
                    wtus_c = wtrp * (one(T) - wtroof) / canyon_resistance
                    wtus_road_perv = wtus_c
                    wtus_road_perv_unscl = one(T) / canyon_resistance
                    wtuq_c = wtrp * (one(T) - wtroof) / canyon_resistance
                    wtuq_road_perv = wtuq_c
                    wtuq_road_perv_unscl = one(T) / canyon_resistance
                elseif ct == ICOL_ROAD_IMPERV
                    wtus_c = (one(T) - wtrp) * (one(T) - wtroof) / canyon_resistance
                    wtus_road_imperv = wtus_c
                    wtus_road_imperv_unscl = one(T) / canyon_resistance
                    if c_.snow_depth_col[c] > zero(T)
                        fwet = min(c_.snow_depth_col[c] / T(0.05), one(T))
                    else
                        fwet = (max(zero(T), c_.h2osoi_liq_col[c, 1] + c_.h2osoi_ice_col[c, 1]) / pondmx)^T(0.666666666666)
                        fwet = min(fwet, one(T))
                    end
                    if o.qaf_lun[l] > c_.qg_col[c]
                        fwet = one(T)
                    end
                    wtuq_c = fwet * (one(T) - wtrp) * (one(T) - wtroof) / canyon_resistance
                    wtuq_road_imperv = wtuq_c
                    wtuq_road_imperv_unscl = fwet * (one(T) / canyon_resistance)
                elseif ct == ICOL_SUNWALL
                    wtus_c = hwr * (one(T) - wtroof) / canyon_resistance
                    wtus_sunwall = wtus_c
                    wtus_sunwall_unscl = one(T) / canyon_resistance
                    wtuq_c = zero(T)
                    wtuq_sunwall = zero(T)
                elseif ct == ICOL_SHADEWALL
                    wtus_c = hwr * (one(T) - wtroof) / canyon_resistance
                    wtus_shadewall = wtus_c
                    wtus_shadewall_unscl = one(T) / canyon_resistance
                    wtuq_c = zero(T)
                    wtuq_shadewall = zero(T)
                end
                c_.wtus_col[c] = wtus_c
                c_.wtuq_col[c] = wtuq_c

                taf_numer = taf_numer + c_.t_grnd_col[c] * wtus_c
                taf_denom = taf_denom + wtus_c
                qaf_numer = qaf_numer + c_.qg_col[c] * wtuq_c
                qaf_denom = qaf_denom + wtuq_c
            end

            # ---- wasteheat / heat-from-AC at landunit level (simple build) ----
            if is_simple
                # recompute per-surface wasteheat the same way the host wasteheat! does
                # (roof + (1-wtroof)*(hwr*(sunwall+shadewall))). Recompute each surface's
                # wh/hfac from its column's AC/heat (walls share canyon_hwr weighting).
                wh_roof = zero(T); hf_roof = zero(T)
                wh_sun = zero(T);  hf_sun = zero(T)
                wh_sha = zero(T);  hf_sha = zero(T)
                for c in ci:cf
                    ct = c_.col_itype[c]
                    wh = wasteheat_on ?
                        ac_wh_factor * c_.eflx_urban_ac_col[c] + ht_wh_factor * c_.eflx_urban_heat_col[c] :
                        zero(T)
                    hf = (urban_hac_on || wasteheat_on) ? abs(c_.eflx_urban_ac_col[c]) : zero(T)
                    if ct == ICOL_ROOF
                        wh_roof = wh; hf_roof = hf
                    elseif ct == ICOL_SUNWALL
                        wh_sun = wh; hf_sun = hf
                    elseif ct == ICOL_SHADEWALL
                        wh_sha = wh; hf_sha = hf
                    end
                end
                eflx_wh_lun = wtroof * wh_roof + (one(T) - wtroof) * (hwr * (wh_sun + wh_sha))
                eflx_wh_lun = min(eflx_wh_lun, wasteheat_limit)
                eflx_hfac_lun = wtroof * hf_roof + (one(T) - wtroof) * (hwr * (hf_sun + hf_sha))
                o.eflx_wasteheat_lun[l] = eflx_wh_lun
                o.eflx_heat_from_ac_lun[l] = eflx_hfac_lun
            else
                # prog build temp: wasteheat from landunit-level AC/heat (handled on host
                # via energyflux.eflx_urban_*_lun; not exercised here). Leave untouched.
            end

            # traffic heat flux (impervious road only)
            o.eflx_traffic_lun[l] = (one(T) - wtroof) * (one(T) - wtrp) * o.eflx_traffic_factor[l]

            # ---- new canopy-air temperature/humidity ----
            taf_lun_out[l] = taf_numer / taf_denom
            o.qaf_lun[l] = qaf_numer / qaf_denom

            wts_sum = wtas + wtus_roof + wtus_road_perv + wtus_road_imperv + wtus_sunwall + wtus_shadewall
            wtq_sum = wtaq + wtuq_roof + wtuq_road_perv + wtuq_road_imperv + wtuq_sunwall + wtuq_shadewall

            # ---- stability update using new taf/qaf ----
            dth2 = thm_g - taf_lun_out[l]
            dqh2 = o.forc_q_grc[g] - o.qaf_lun[l]
            tstar = temp1 * dth2
            qstar = temp2 * dqh2
            thvstar = tstar * (one(T) + T(0.61) * o.forc_q_grc[g]) + T(0.61) * o.forc_th_grc[g] * qstar
            zeta_l = zldis * vkc * grav * thvstar / (ustar^2 * thv_g)
            if zeta_l >= zero(T)
                zeta_l = min(zetamax, max(zeta_l, T(0.01)))
                um = max(ur, T(0.1))
            else
                zeta_l = max(T(-100.0), min(zeta_l, T(-0.01)))
                wc2 = one(T) * (-grav * ustar * thvstar * T(1000.0) / thv_g)^T(0.333)
                um = sqrt(ur^2 + wc2^2)
            end
            obu = zldis / zeta_l

            # commit per-landunit outputs every iteration (last iteration wins,
            # matching the host which leaves these at their final-iteration values)
            o.ramu[l] = ramu
            o.rahu[l] = rahu
            o.rawu[l] = rawu
            o.zeta_lunit[l] = zeta_l
            o.wtas[l] = wtas; o.wtaq[l] = wtaq
            o.wts_sum[l] = wts_sum; o.wtq_sum[l] = wtq_sum
            o.wtus_roof[l] = wtus_roof; o.wtuq_roof[l] = wtuq_roof
            o.wtus_road_perv[l] = wtus_road_perv; o.wtuq_road_perv[l] = wtuq_road_perv
            o.wtus_road_imperv[l] = wtus_road_imperv; o.wtuq_road_imperv[l] = wtuq_road_imperv
            o.wtus_sunwall[l] = wtus_sunwall; o.wtuq_sunwall[l] = wtuq_sunwall
            o.wtus_shadewall[l] = wtus_shadewall; o.wtuq_shadewall[l] = wtuq_shadewall
            o.wtus_roof_unscl[l] = wtus_roof_unscl; o.wtuq_roof_unscl[l] = wtuq_roof_unscl
            o.wtus_road_perv_unscl[l] = wtus_road_perv_unscl; o.wtuq_road_perv_unscl[l] = wtuq_road_perv_unscl
            o.wtus_road_imperv_unscl[l] = wtus_road_imperv_unscl; o.wtuq_road_imperv_unscl[l] = wtuq_road_imperv_unscl
            o.wtus_sunwall_unscl[l] = wtus_sunwall_unscl
            o.wtus_shadewall_unscl[l] = wtus_shadewall_unscl
        end
    end
end

# ---------------------------------------------------------------------------
# Error-check kernel: per-landunit, in-thread sum over the landunit's member
# patches of the scaled fluxes vs the bulk fluxes (matches the host sum()-based
# reduction). Writes eflx_err/qflx_err which the host reads back for the warn/error.
# ---------------------------------------------------------------------------
@kernel function _uf_errcheck_kernel!(eflx_err, @Const(mask), @Const(qflx_err_in),
        @Const(eflx_sh_grnd_scale), @Const(qflx_evap_soi_scale),
        @Const(patchi), @Const(patchf), @Const(lun_gridcell),
        @Const(forc_rho_grc), @Const(forc_q_grc), @Const(thm_g_lun),
        @Const(taf_lun), @Const(qaf_lun), @Const(rahu), @Const(rawu),
        qflx_err, eflx_arr_out, qflx_arr_out, eflx_scale_out, qflx_scale_out,
        CPAIR_)
    fl = @index(Global)
    @inbounds if mask[fl]
        l = fl
        T = eltype(eflx_err)
        cpair = T(CPAIR_)
        g = lun_gridcell[l]
        eflx_arr = -(forc_rho_grc[g] * cpair / rahu[l]) * (thm_g_lun[l] - taf_lun[l])
        qflx_arr = -(forc_rho_grc[g] / rawu[l]) * (forc_q_grc[g] - qaf_lun[l])
        es = zero(T); qs = zero(T)
        for p in patchi[l]:patchf[l]
            es += eflx_sh_grnd_scale[p]
            qs += qflx_evap_soi_scale[p]
        end
        eflx_arr_out[l] = eflx_arr
        qflx_arr_out[l] = qflx_arr
        eflx_scale_out[l] = es
        qflx_scale_out[l] = qs
        eflx_err[l] = es - eflx_arr
        qflx_err[l] = qs - qflx_arr
    end
end

"""
    urban_fluxes_iterate!(energyflux, frictionvel, temperature, soilstate,
        urbanparams, waterstatebulk, waterdiagbulk, lun_data, col_data,
        m_urbanl, num_urbanl, filter_urbanl,
        ramu, zeta_lunit, wtus_col, wtuq_col, rahu, rawu, thm_g_lun, wts, wtq,
        forc_*, niters) -> nothing

Run the WHOLE iterative canyon-air solve (canyontop wind, monin_obuk_ini init, the
niters stability loop with the FrictionVelocity landunit math inlined, the
per-column gather into taf/qaf, wasteheat/traffic, and the stability update) as a
single KernelAbstractions kernel — ONE THREAD PER URBAN LANDUNIT. Runs WHOLE on the
device (Metal) when the state structs + masks live there, and is CPU byte-identical
on the host path. Writes the per-landunit solved scalars / weights, the per-column
wtus/wtuq, and the per-member-patch friction-velocity side effects.
"""
function urban_fluxes_iterate!(
        energyflux, frictionvel, temperature, soilstate, urbanparams,
        waterstatebulk, waterdiagbulk, lun_data, col_data,
        m_urbanl::AbstractVector{Bool}, num_urbanl::Int, filter_urbanl::Vector{Int},
        ramu, zeta_lunit, wtus_col, wtuq_col, rahu, rawu, thm_g_lun,
        wts::_UfWtsDV, wtq::_UfWtqDV,
        forc_t_grc, forc_th_grc, forc_rho_grc, forc_q_grc, forc_u_grc, forc_v_grc,
        niters::Int; lapse_rate::Real = 0.0098)

    proto = temperature.taf_lun
    FT = eltype(proto)
    b(v) = _uf_to_backend(proto, v)

    # Host-only feature flags + wasteheat factors resolved to scalars / Bools.
    is_simple    = is_simple_build_temp()
    urban_hac_on = (urban_ctrl.urban_hac == URBAN_HAC_ON || urban_ctrl.urban_hac == URBAN_WASTEHEAT_ON)
    wasteheat_on = (urban_ctrl.urban_hac == URBAN_WASTEHEAT_ON)

    o = _UfIterLunDV(
        ramu = ramu, rahu = rahu, rawu = rawu, zeta_lunit = zeta_lunit,
        taf_lun = temperature.taf_lun, qaf_lun = waterdiagbulk.qaf_lun,
        wtas = wts.wtas, wtaq = wtq.wtaq, wts_sum = wts.wts_sum, wtq_sum = wtq.wtq_sum,
        wtus_roof = wts.wtus_roof, wtuq_roof = wtq.wtuq_roof,
        wtus_road_perv = wts.wtus_road_perv, wtuq_road_perv = wtq.wtuq_road_perv,
        wtus_road_imperv = wts.wtus_road_imperv, wtuq_road_imperv = wtq.wtuq_road_imperv,
        wtus_sunwall = wts.wtus_sunwall, wtuq_sunwall = wtq.wtuq_sunwall,
        wtus_shadewall = wts.wtus_shadewall, wtuq_shadewall = wtq.wtuq_shadewall,
        wtus_roof_unscl = wts.wtus_roof_unscl, wtuq_roof_unscl = wtq.wtuq_roof_unscl,
        wtus_road_perv_unscl = wts.wtus_road_perv_unscl, wtuq_road_perv_unscl = wtq.wtuq_road_perv_unscl,
        wtus_road_imperv_unscl = wts.wtus_road_imperv_unscl, wtuq_road_imperv_unscl = wtq.wtuq_road_imperv_unscl,
        wtus_sunwall_unscl = wts.wtus_sunwall_unscl,
        wtus_shadewall_unscl = wts.wtus_shadewall_unscl,
        eflx_wasteheat_lun = energyflux.eflx_wasteheat_lun,
        eflx_heat_from_ac_lun = energyflux.eflx_heat_from_ac_lun,
        eflx_traffic_lun = energyflux.eflx_traffic_lun,
        lun_gridcell = lun_data.gridcell, coli = lun_data.coli, colf = lun_data.colf,
        patchi = lun_data.patchi, patchf = lun_data.patchf,
        ht_roof = lun_data.ht_roof, z_d_town = lun_data.z_d_town, z_0_town = lun_data.z_0_town,
        canyon_hwr = lun_data.canyon_hwr, wtlunit_roof = lun_data.wtlunit_roof,
        wtroad_perv = lun_data.wtroad_perv,
        wind_hgt_canyon = urbanparams.wind_hgt_canyon,
        eflx_traffic_factor = urbanparams.eflx_traffic_factor,
        forc_u_grc = b(forc_u_grc), forc_v_grc = b(forc_v_grc),
        forc_t_grc = b(forc_t_grc), forc_th_grc = b(forc_th_grc),
        forc_q_grc = b(forc_q_grc), forc_rho_grc = b(forc_rho_grc),
        forc_hgt_u_patch = frictionvel.forc_hgt_u_patch,
        forc_hgt_t_patch = frictionvel.forc_hgt_t_patch)

    c_ = _UfIterColDV(
        wtus_col = wtus_col, wtuq_col = wtuq_col,
        t_grnd_col = temperature.t_grnd_col, qg_col = waterdiagbulk.qg_col,
        snow_depth_col = waterdiagbulk.snow_depth_col,
        h2osoi_liq_col = waterstatebulk.ws.h2osoi_liq_col,
        h2osoi_ice_col = waterstatebulk.ws.h2osoi_ice_col,
        eflx_urban_ac_col = energyflux.eflx_urban_ac_col,
        eflx_urban_heat_col = energyflux.eflx_urban_heat_col,
        col_itype = col_data.itype,
        u10_clm_patch = frictionvel.u10_clm_patch, va_patch = frictionvel.va_patch,
        u10_patch = frictionvel.u10_patch, fv_patch = frictionvel.fv_patch,
        vds_patch = frictionvel.vds_patch)

    # thm_g per landunit (needed by the error check downstream); iteration-invariant,
    # so computed once here (host) rather than threaded out of the kernel. Host copies
    # of the forcing + forc_hgt avoid scalar-indexing a device array. thm_g_lun itself
    # may be device-resident, so it is written via a host staging vector + copyto!.
    forc_t_h    = Array(forc_t_grc)
    forc_hgt_h  = Array(frictionvel.forc_hgt_t_patch)
    gridcell_h  = Array(lun_data.gridcell)
    patchi_h    = Array(lun_data.patchi)
    thm_g_h     = zeros(FT, length(thm_g_lun))
    @inbounds for fl in 1:num_urbanl
        l = filter_urbanl[fl]
        g = gridcell_h[l]
        thm_g_h[l] = forc_t_h[g] + FT(lapse_rate) * forc_hgt_h[patchi_h[l]]
    end
    copyto!(thm_g_lun, thm_g_h)

    _launch!(_uf_iter_kernel!, temperature.taf_lun, m_urbanl, o, c_,
        niters, convert(FT, frictionvel.zetamaxstable), convert(FT, urban_fluxes_params.wind_min),
        convert(FT, VKC), convert(FT, GRAV), convert(FT, CPAIR), convert(FT, RPI),
        convert(FT, PONDMX_URBAN), convert(FT, lapse_rate),
        is_simple, urban_hac_on, wasteheat_on,
        convert(FT, AC_WASTEHEAT_FACTOR), convert(FT, HT_WASTEHEAT_FACTOR),
        convert(FT, WASTEHEAT_LIMIT);
        ndrange = length(temperature.taf_lun))
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
        # Atmospheric forcing (gridcell-level). Widened to AbstractVector so device
        # (MtlArray) forcing dispatches; plain-Vector callsites are unaffected.
        forc_t_grc       ::AbstractVector{<:Real},
        forc_th_grc      ::AbstractVector{<:Real},
        forc_rho_grc     ::AbstractVector{<:Real},
        forc_q_grc       ::AbstractVector{<:Real},
        forc_pbot_grc    ::AbstractVector{<:Real},
        forc_u_grc       ::AbstractVector{<:Real},
        forc_v_grc       ::AbstractVector{<:Real};
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

    FT = eltype(temperature.taf_lun)

    # Per-landunit / per-column solved outputs of the iterative canyon-air solve.
    # These are allocated on the STATE backend (device-resident on Metal) so the
    # iteration kernel and the downstream canyon-surface kernel both read/write them
    # in place — no host<->device round trips. `_uf_dev_zeros` mirrors the working
    # precision + backend of the state. Bundled (wts/wtq) so the iteration and the
    # canyon-surface pass share the same device arrays.
    proto_l = temperature.taf_lun
    dz(n) = (a = similar(proto_l, FT, n); fill!(a, zero(FT)); a)
    proto_c = temperature.t_grnd_col
    dzc(n) = (a = similar(proto_c, FT, n); fill!(a, zero(FT)); a)

    ramu       = dz(nl)
    rahu       = dz(nl)
    rawu       = dz(nl)
    zeta_lunit = dz(nl)
    thm_g_lun  = dz(nl)
    wtus_col   = dzc(nc)
    wtuq_col   = dzc(nc)

    wts_dv = _UfWtsDV(
        wtas = dz(nl),
        wtus_roof = dz(nl), wtus_road_perv = dz(nl),
        wtus_road_imperv = dz(nl), wtus_sunwall = dz(nl),
        wtus_shadewall = dz(nl),
        wtus_roof_unscl = dz(nl), wtus_road_perv_unscl = dz(nl),
        wtus_road_imperv_unscl = dz(nl), wtus_sunwall_unscl = dz(nl),
        wtus_shadewall_unscl = dz(nl), wts_sum = dz(nl))
    wtq_dv = _UfWtqDV(
        wtaq = dz(nl),
        wtuq_roof = dz(nl), wtuq_road_perv = dz(nl),
        wtuq_road_imperv = dz(nl), wtuq_sunwall = dz(nl),
        wtuq_shadewall = dz(nl),
        wtuq_roof_unscl = dz(nl), wtuq_road_perv_unscl = dz(nl),
        wtuq_road_imperv_unscl = dz(nl), wtq_sum = dz(nl))

    proto_pp = energyflux.cgrnd_patch
    eflx_sh_grnd_scale  = (a = similar(proto_pp, FT, np_local); fill!(a, zero(FT)); a)
    qflx_evap_soi_scale = (a = similar(proto_pp, FT, np_local); fill!(a, zero(FT)); a)

    eflx_err       = dz(nl)
    qflx_err       = dz(nl)
    eflx_arr       = dz(nl)
    qflx_arr       = dz(nl)
    eflx_scale_arr = dz(nl)
    qflx_scale_arr = dz(nl)

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
    # Per-landunit URBAN mask (over the full landunit range) for the iteration +
    # error-check kernels; one thread per landunit, urban ones do work.
    m_urbanl = _uf_device_mask(temperature.taf_lun, filter_urbanl, num_urbanl, nl)

    spval = convert(FT, SPVAL)
    _launch!(_uf_nourban_restart_kernel!, temperature.taf_lun,
             waterdiagbulk.qaf_lun, m_nourbanl, spval)

    # =========================================================================
    # Aerodynamic-parameter guards (HOST, validation only — no state change). The
    # iteration kernel cannot `error()`/`@warn`; these checks read host-resident
    # lun geometry + forcing heights and gate exactly the original error paths.
    # Host copies (single small read-back) so a device-resident state doesn't trip
    # GPU scalar-indexing; the validation is identical to the original.
    # =========================================================================
    htr_h   = Array(lun_data.ht_roof)
    zdt_h   = Array(lun_data.z_d_town)
    z0t_h   = Array(lun_data.z_0_town)
    hgtu_h  = Array(frictionvel.forc_hgt_u_patch)
    lpatchi = Array(lun_data.patchi)
    for fl in 1:num_urbanl
        l = filter_urbanl[fl]
        if htr_h[l] - zdt_h[l] <= z0t_h[l]
            if _is_ad_type(eltype(frictionvel.forc_hgt_u_patch))
                @warn "urban_fluxes! aerodynamic parameter error (AD mode, continuing)" maxlog=1
            else
                error("aerodynamic parameter error in urban_fluxes!: ht_roof - z_d_town <= z_0_town " *
                      "ht_roof=$(htr_h[l]) z_d_town=$(zdt_h[l]) z_0_town=$(z0t_h[l])")
            end
        end
        if hgtu_h[lpatchi[l]] - zdt_h[l] <= z0t_h[l]
            if _is_ad_type(eltype(frictionvel.forc_hgt_u_patch))
                @warn "urban_fluxes! aerodynamic parameter error (AD mode, continuing)" maxlog=1
            else
                error("aerodynamic parameter error in urban_fluxes!: forc_hgt_u - z_d_town <= z_0_town " *
                      "forc_hgt_u=$(hgtu_h[lpatchi[l]]) z_d_town=$(zdt_h[l]) z_0_town=$(z0t_h[l])")
            end
        end
    end

    # =========================================================================
    # ITERATIVE CANYON-AIR SOLVE — WHOLE on device, one thread per urban landunit.
    # (canyontop wind + monin_obuk_ini + niters stability loop + column gather +
    #  wasteheat/traffic + stability update). See urban_fluxes_iterate!.
    # =========================================================================
    urban_fluxes_iterate!(
        energyflux, frictionvel, temperature, soilstate, urbanparams,
        waterstatebulk, waterdiagbulk, lun_data, col_data,
        m_urbanl, num_urbanl, filter_urbanl,
        ramu, zeta_lunit, wtus_col, wtuq_col, rahu, rawu, thm_g_lun,
        wts_dv, wtq_dv,
        forc_t_grc, forc_th_grc, forc_rho_grc, forc_q_grc, forc_u_grc, forc_v_grc,
        niters; lapse_rate = lapse_rate)

    # Prog build temp: wasteheat is computed from landunit-level AC/heat (untouched
    # by the simple-build kernel path). Run the host wasteheat! to fill it. (Simple
    # build temp is handled inside the iteration kernel.)
    if is_prog_build_temp()
        empty_l = FT[]
        wasteheat!(energyflux, lun_data, num_urbanl, filter_urbanl,
            empty_l, empty_l, empty_l, empty_l, empty_l, empty_l)
    end

    # =========================================================================
    # Determine fluxes from canyon surfaces
    # =========================================================================
    # Per-patch + column-INDEPENDENT given the SOLVED canopy-air state: runs WHOLE
    # on the device as a single kernel (see urban_fluxes_canyon_surface!). The
    # per-landunit weight arrays (wts_dv/wtq_dv) and ramu/zeta_lunit/wtus_col/
    # wtuq_col are ALREADY device-resident (written by the iteration kernel above),
    # so they pass straight through. Only the gridcell forcing vectors are still
    # host-resident; _uf_to_backend moves them onto the state backend (identity /
    # byte-identical on Float64 CPU). The scaled-flux scratch was zero-initialised
    # on-device (matches the begp:endp init), then overwritten for urban patches.
    bf(v) = _uf_to_backend(energyflux.cgrnd_patch, v)
    urban_fluxes_canyon_surface!(
        energyflux, frictionvel, waterfluxbulk, waterdiagbulk, soilstate,
        temperature, patch_data, col_data, m_urbanp, wts_dv, wtq_dv,
        ramu, zeta_lunit, wtus_col, wtuq_col,
        eflx_sh_grnd_scale, qflx_evap_soi_scale,
        bf(forc_rho_grc), bf(forc_u_grc), bf(forc_v_grc))

    # =========================================================================
    # Error checking: total fluxes should equal sum of scaled fluxes — per-landunit
    # in-thread reduction over the landunit's member patches (one thread per
    # landunit), so the whole check runs on-device. The warn/error gating reads the
    # small per-landunit err arrays back on the host (a kernel cannot error()/@warn).
    # =========================================================================
    _launch!(_uf_errcheck_kernel!, eflx_err, m_urbanl, qflx_err,
        eflx_sh_grnd_scale, qflx_evap_soi_scale,
        lun_data.patchi, lun_data.patchf, lun_data.gridcell,
        bf(forc_rho_grc), bf(forc_q_grc), thm_g_lun,
        temperature.taf_lun, waterdiagbulk.qaf_lun, rahu, rawu,
        qflx_err, eflx_arr, qflx_arr, eflx_scale_arr, qflx_scale_arr,
        convert(FT, CPAIR); ndrange = nl)

    # Read back the small per-landunit error arrays for the host warn/error gating.
    eflx_err_h = collect(eflx_err); qflx_err_h = collect(qflx_err)
    eflx_arr_h = collect(eflx_arr); qflx_arr_h = collect(qflx_arr)
    eflx_scale_h = collect(eflx_scale_arr); qflx_scale_h = collect(qflx_scale_arr)

    # Check sensible heat flux error
    for fl in 1:num_urbanl
        l = filter_urbanl[fl]
        if abs(eflx_err_h[l]) > 0.01
            @warn "Total sensible heat does not equal sum of scaled heat fluxes for urban columns" nstep l eflx_err=eflx_err_h[l]
            if abs(eflx_err_h[l]) > 0.01
                if _is_ad_type(eltype(eflx_err))
                    @warn "urban_fluxes! sensible heat flux error (AD mode, continuing)" maxlog=1
                else
                    error("urban_fluxes! sensible heat flux error > 0.01 W/m**2: " *
                          "eflx_scale=$(eflx_scale_h[l]) eflx=$(eflx_arr_h[l]) err=$(eflx_err_h[l])")
                end
            end
        end
    end

    # Check water vapor flux error
    for fl in 1:num_urbanl
        l = filter_urbanl[fl]
        if abs(qflx_err_h[l]) > 4.0e-9
            @warn "Total water vapor flux does not equal sum of scaled fluxes for urban columns" nstep l qflx_err=qflx_err_h[l]
            if abs(qflx_err_h[l]) > 4.0e-9
                if _is_ad_type(eltype(qflx_err))
                    @warn "urban_fluxes! water vapor flux error (AD mode, continuing)" maxlog=1
                else
                    error("urban_fluxes! water vapor flux error > 4e-9 kg/m**2/s: " *
                          "qflx_scale=$(qflx_scale_h[l]) qflx=$(qflx_arr_h[l]) err=$(qflx_err_h[l])")
                end
            end
        end
    end

    # =========================================================================
    # Internal building temperature (simple method)
    # =========================================================================
    if is_simple_build_temp()
        calc_simple_internal_building_temp_dev!(temperature, col_data, lun_data, m_urbanl)
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
