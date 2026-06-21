# ==========================================================================
# Ported from: src/biogeophys/UrbBuildTempOleson2015Mod.F90 (~960 lines)
#
# Prognostic (CLM5.0 / Oleson 2015) interior building-temperature model.
#
# Solves for five unknowns at each urban landunit:
#   t_roof_inner, t_sunw_inner, t_shdw_inner, t_floor, t_building
# from a 5-equation linear system (energy balance at each inner surface plus
# the building-air heat budget). The original Fortran solves the 5x5 system
# with LAPACK dgesv; here we use a self-contained partial-pivoted Gaussian
# elimination (`_solve5x5!`) so the routine is dependency-free and AD-friendly.
#
# Public function:
#   building_temperature!  — prognostic building energy balance (CLM5)
#
# This complements the simple (CLM4.5) `building_hac!` in soil_temperature.jl.
# Wired into soil_temperature! via the is_prog_build_temp() branch.
#
# Conventions (CLAUDE.md): Fortran variable names preserved, SoA layout,
# mask-based loops, bang function, Float64 (parametric working type for AD).
# ==========================================================================

# Urban interior building-temperature physical constants (clm_varcon.F90).
# (Defined here next to their sole consumer; values from CLM5.)
const EM_ROOF_INT  = 0.9      # emissivity of interior surface of roof (Bueno et al. 2012, GMD)
const EM_SUNW_INT  = 0.9      # emissivity of interior surface of sunwall
const EM_SHDW_INT  = 0.9      # emissivity of interior surface of shadewall
const EM_FLOOR_INT = 0.9      # emissivity of interior surface of floor
const HCV_ROOF          = 0.948   # interior convective heat transfer coeff, roof (W m-2 K-1)
const HCV_ROOF_ENHANCED = 4.040   # enhanced (t_roof_int <= t_room) roof coeff (W m-2 K-1)
const HCV_FLOOR          = 0.948  # interior convective heat transfer coeff, floor (W m-2 K-1)
const HCV_FLOOR_ENHANCED = 4.040  # enhanced (t_floor_int >= t_room) floor coeff (W m-2 K-1)
const HCV_SUNW = 3.076        # interior convective heat transfer coeff, sunwall (W m-2 K-1)
const HCV_SHDW = 3.076        # interior convective heat transfer coeff, shadewall (W m-2 K-1)
const DZ_FLOOR   = 0.1        # floor thickness - concrete (Salamanca et al. 2010) (m)
const DENS_FLOOR = 2.35e3     # density of floor - concrete (kg m-3)
const SH_FLOOR   = 880.0      # specific heat of floor - concrete (J kg-1 K-1)
const CP_FLOOR   = DENS_FLOOR * SH_FLOOR  # volumetric heat capacity of floor (J m-3 K-1)
const VENT_ACH   = 0.3        # ventilation rate (air exchanges per hour)

# --------------------------------------------------------------------------
# Dense 5x5 linear solve A x = b by partial-pivoted Gaussian elimination,
# in-place on a (mutable) 5x5 matrix `a` and length-5 vector `b`. On return
# `b` holds the solution x. Carries the working element type (Float64 or a
# ForwardDiff.Dual) so AD flows through. Mirrors LAPACK dgesv for n=5.
# --------------------------------------------------------------------------
@inline function _solve5x5!(a::AbstractMatrix{T}, b::AbstractVector{T}) where {T}
    n = 5
    @inbounds for k in 1:n
        # Partial pivot: find row p>=k with max |a[p,k]|
        p = k
        amax = abs(a[k, k])
        for i in (k+1):n
            ai = abs(a[i, k])
            if ai > amax
                amax = ai
                p = i
            end
        end
        if p != k
            for j in 1:n
                a[k, j], a[p, j] = a[p, j], a[k, j]
            end
            b[k], b[p] = b[p], b[k]
        end
        akk = a[k, k]
        # Eliminate below the pivot
        for i in (k+1):n
            f = a[i, k] / akk
            for j in (k+1):n
                a[i, j] -= f * a[k, j]
            end
            b[i] -= f * b[k]
        end
    end
    # Back substitution
    @inbounds for k in n:-1:1
        s = b[k]
        for j in (k+1):n
            s -= a[k, j] * b[j]
        end
        b[k] = s / a[k, k]
    end
    return b
end

"""
    building_temperature!(col, lun, temperature, energyflux, urbanparams,
                          urbantv_t_building_max, urbantv_p_ac,
                          tk, mask_urbanl, mask_nolakec, bounds_lun, dtime)

Prognostic (CLM5.0) interior building-air temperature.

Solves the coupled energy balance at the roof/sunwall/shadewall/floor inner
surfaces and the building air, updating `t_roof_inner_lun`, `t_sunw_inner_lun`,
`t_shdw_inner_lun`, `t_floor_lun`, `t_building_lun`, and the building energy
fluxes (`eflx_building_lun`, `eflx_urban_ac_lun`, `eflx_urban_heat_lun`,
`eflx_ventilation_lun`).

`tk` is the interface thermal conductivity matrix (col, nlevsno+level) as
produced inside `soil_temperature!`.

Ported from `BuildingTemperature` in `UrbBuildTempOleson2015Mod.F90`.
Host (CPU) implementation; runs per urban landunit (mask-based loop).
"""
function building_temperature!(col::ColumnData, lun::LandunitData,
                               temperature::TemperatureData,
                               energyflux::EnergyFluxData,
                               urbanparams::UrbanParamsData,
                               urbantv_t_building_max::AbstractVector{<:Real},
                               urbantv_p_ac::AbstractVector{<:Real},
                               tk::AbstractMatrix{<:Real},
                               mask_urbanl::AbstractVector{Bool},
                               mask_nolakec::AbstractVector{Bool},
                               bounds_lun::UnitRange{Int},
                               dtime::Real)

    nlevsno = varpar.nlevsno
    nlevurb = varpar.nlevurb
    joff = nlevsno

    FT = eltype(temperature.t_soisno_col)

    # Constants at the working type
    sb     = FT(SB)
    rair   = FT(RAIR)
    pstd   = FT(PSTD)
    cpair  = FT(CPAIR)
    dt     = FT(dtime)

    em_roofi_c  = FT(EM_ROOF_INT)
    em_sunwi_c  = FT(EM_SUNW_INT)
    em_shdwi_c  = FT(EM_SHDW_INT)
    em_floori_c = FT(EM_FLOOR_INT)
    dz_floori   = FT(DZ_FLOOR)
    cp_floori   = FT(CP_FLOOR)
    vent_ach    = FT(VENT_ACH)

    # State aliases
    t_soisno     = temperature.t_soisno_col
    tssbef       = temperature.t_ssbef_col
    taf          = temperature.taf_lun
    t_roof_inner = temperature.t_roof_inner_lun
    t_sunw_inner = temperature.t_sunw_inner_lun
    t_shdw_inner = temperature.t_shdw_inner_lun
    t_floor      = temperature.t_floor_lun
    t_building   = temperature.t_building_lun

    canyon_hwr   = lun.canyon_hwr
    wtlunit_roof = lun.wtlunit_roof
    ht_roof      = lun.ht_roof
    urbpoi       = lun.urbpoi

    t_building_min = urbanparams.t_building_min
    t_building_max = urbantv_t_building_max
    p_ac           = urbantv_p_ac

    eflx_building    = energyflux.eflx_building_lun
    eflx_urban_ac    = energyflux.eflx_urban_ac_lun
    eflx_urban_heat  = energyflux.eflx_urban_heat_lun
    eflx_ventilation = energyflux.eflx_ventilation_lun

    # HAC mode / explicit AC flags (resolved on host)
    urban_hac         = urban_ctrl.urban_hac
    urban_explicit_ac = urban_ctrl.urban_explicit_ac
    hac_active = (urban_hac == URBAN_HAC_ON) || (urban_hac == URBAN_WASTEHEAT_ON)

    # Inner-node (nlevurb) Julia level index = nlevurb + joff. The corresponding
    # lower interface is zi[c, nlevurb + joff + 1] (matches _fn1_kernel!).
    nlevurb_jj = nlevurb + joff

    # ----------------------------------------------------------------------
    # Gather, per urban landunit, the inner-node (nlevurb) geometry/temps/conductivity
    # from each of the three wall/roof columns of the landunit. (Fortran loops
    # over num_nolakec; here we sweep masked columns and dispatch by column type.)
    # ----------------------------------------------------------------------
    nl = length(bounds_lun)
    zi_roof_innerl = fill(FT(NaN), nl); z_roof_innerl = fill(FT(NaN), nl)
    zi_sunw_innerl = fill(FT(NaN), nl); z_sunw_innerl = fill(FT(NaN), nl)
    zi_shdw_innerl = fill(FT(NaN), nl); z_shdw_innerl = fill(FT(NaN), nl)
    t_roof_innerl_bef = fill(FT(NaN), nl); t_roof_innerl = fill(FT(NaN), nl)
    t_sunw_innerl_bef = fill(FT(NaN), nl); t_sunw_innerl = fill(FT(NaN), nl)
    t_shdw_innerl_bef = fill(FT(NaN), nl); t_shdw_innerl = fill(FT(NaN), nl)
    tk_roof_innerl = fill(FT(NaN), nl)
    tk_sunw_innerl = fill(FT(NaN), nl)
    tk_shdw_innerl = fill(FT(NaN), nl)

    @inbounds for c in eachindex(mask_nolakec)
        mask_nolakec[c] || continue
        l = col.landunit[c]
        urbpoi[l] || continue
        it = col.itype[c]
        if it == ICOL_ROOF
            zi_roof_innerl[l]    = col.zi[c, nlevurb_jj + 1]
            z_roof_innerl[l]     = col.z[c, nlevurb_jj]
            t_roof_innerl_bef[l] = tssbef[c, nlevurb_jj]
            t_roof_innerl[l]     = t_soisno[c, nlevurb_jj]
            tk_roof_innerl[l]    = tk[c, nlevurb_jj]
        elseif it == ICOL_SUNWALL
            zi_sunw_innerl[l]    = col.zi[c, nlevurb_jj + 1]
            z_sunw_innerl[l]     = col.z[c, nlevurb_jj]
            t_sunw_innerl_bef[l] = tssbef[c, nlevurb_jj]
            t_sunw_innerl[l]     = t_soisno[c, nlevurb_jj]
            tk_sunw_innerl[l]    = tk[c, nlevurb_jj]
        elseif it == ICOL_SHADEWALL
            zi_shdw_innerl[l]    = col.zi[c, nlevurb_jj + 1]
            z_shdw_innerl[l]     = col.z[c, nlevurb_jj]
            t_shdw_innerl_bef[l] = tssbef[c, nlevurb_jj]
            t_shdw_innerl[l]     = t_soisno[c, nlevurb_jj]
            tk_shdw_innerl[l]    = tk[c, nlevurb_jj]
        end
    end

    # ----------------------------------------------------------------------
    # Per urban landunit: assemble and solve the 5x5 system, then HAC clamp.
    # ----------------------------------------------------------------------
    a = Matrix{FT}(undef, 5, 5)
    res = Vector{FT}(undef, 5)

    @inbounds for l in bounds_lun
        mask_urbanl[l] || continue
        urbpoi[l]      || continue

        # 1. Save previous-step temperatures
        t_roof_inner_bef = t_roof_inner[l]
        t_sunw_inner_bef = t_sunw_inner[l]
        t_shdw_inner_bef = t_shdw_inner[l]
        t_floor_bef      = t_floor[l]
        t_building_bef   = t_building[l]

        # 2. Convective heat transfer coefficients (Bueno et al. 2012)
        hcv_roofi  = t_roof_inner_bef <= t_building_bef ? FT(HCV_ROOF_ENHANCED)  : FT(HCV_ROOF)
        hcv_floori = t_floor_bef      >= t_building_bef ? FT(HCV_FLOOR_ENHANCED) : FT(HCV_FLOOR)
        hcv_sunwi  = FT(HCV_SUNW)
        hcv_shdwi  = FT(HCV_SHDW)

        em_roofi  = em_roofi_c
        em_sunwi  = em_sunwi_c
        em_shdwi  = em_shdwi_c
        em_floori = em_floori_c

        # Intermediate calc for concrete floor (W m-2 K-1)
        cv_floori = (dz_floori * cp_floori) / dt
        # Density of dry air at standard pressure and t_building
        rho_dair = pstd / (rair * t_building_bef)
        # Building height to building width ratio
        building_hwr = canyon_hwr[l] * (one(FT) - wtlunit_roof[l]) / wtlunit_roof[l]

        # 3. View factors (inner box)
        vf_rf = sqrt(one(FT) + building_hwr^2) - building_hwr
        vf_fr = vf_rf
        vf_wf = FT(0.5) * (one(FT) - vf_rf)   # per unit wall -> per unit floor
        vf_fw = vf_wf
        vf_rw = vf_fw
        vf_wr = vf_wf                          # per unit wall -> per unit roof
        vf_ww = one(FT) - vf_rw - vf_fw

        # Geometry / conductivity at the inner (nlevurb) node
        zir = zi_roof_innerl[l]; zr = z_roof_innerl[l]; tkr = tk_roof_innerl[l]
        zis = zi_sunw_innerl[l]; zs = z_sunw_innerl[l]; tks = tk_sunw_innerl[l]
        zih = zi_shdw_innerl[l]; zh = z_shdw_innerl[l]; tkh = tk_shdw_innerl[l]
        tril = t_roof_innerl[l];  trilb = t_roof_innerl_bef[l]
        tsil = t_sunw_innerl[l];  tsilb = t_sunw_innerl_bef[l]
        thil = t_shdw_innerl[l];  thilb = t_shdw_innerl_bef[l]

        trib = t_roof_inner_bef
        tsib = t_sunw_inner_bef
        thib = t_shdw_inner_bef
        tfb  = t_floor_bef
        tbb  = t_building_bef

        # ---- ROOF (equation 1) ----
        a[1,1] =   FT(0.5)*hcv_roofi +
                   FT(0.5)*tkr/(zir - zr) +
                   FT(4)*em_roofi*sb*trib^3 -
                   FT(4)*em_roofi*sb*trib^3*vf_rw*(one(FT)-em_sunwi)*vf_wr -
                   FT(4)*em_roofi*sb*trib^3*vf_rw*(one(FT)-em_shdwi)*vf_wr -
                   FT(4)*em_roofi*sb*trib^3*vf_rf*(one(FT)-em_floori)*vf_fr
        a[1,2] = - FT(4)*em_roofi*em_sunwi*sb*tsib^3*vf_wr -
                   FT(4)*em_sunwi*sb*tsib^3*vf_ww*(one(FT)-em_shdwi)*vf_wr -
                   FT(4)*em_sunwi*sb*tsib^3*vf_wf*(one(FT)-em_floori)*vf_fr
        a[1,3] = - FT(4)*em_roofi*em_shdwi*sb*thib^3*vf_wr -
                   FT(4)*em_shdwi*sb*thib^3*vf_ww*(one(FT)-em_sunwi)*vf_wr -
                   FT(4)*em_shdwi*sb*thib^3*vf_wf*(one(FT)-em_floori)*vf_fr
        a[1,4] = - FT(4)*em_roofi*em_floori*sb*tfb^3*vf_fr -
                   FT(4)*em_floori*sb*tfb^3*vf_fw*(one(FT)-em_sunwi)*vf_wr -
                   FT(4)*em_floori*sb*tfb^3*vf_fw*(one(FT)-em_shdwi)*vf_wr
        a[1,5] = - FT(0.5)*hcv_roofi
        res[1] =   FT(0.5)*tkr*tril/(zir - zr) -
                   FT(0.5)*tkr*(trib - trilb)/(zir - zr) -
                   FT(3)*em_roofi*em_sunwi*sb*tsib^4*vf_wr -
                   FT(3)*em_roofi*em_shdwi*sb*thib^4*vf_wr -
                   FT(3)*em_roofi*em_floori*sb*tfb^4*vf_fr +
                   FT(3)*em_roofi*sb*trib^4 -
                   FT(3)*em_roofi*sb*trib^4*vf_rw*(one(FT)-em_sunwi)*vf_wr -
                   FT(3)*em_roofi*sb*trib^4*vf_rw*(one(FT)-em_shdwi)*vf_wr -
                   FT(3)*em_roofi*sb*trib^4*vf_rf*(one(FT)-em_floori)*vf_fr -
                   FT(3)*em_sunwi*sb*tsib^4*vf_ww*(one(FT)-em_shdwi)*vf_wr -
                   FT(3)*em_sunwi*sb*tsib^4*vf_wf*(one(FT)-em_floori)*vf_fr -
                   FT(3)*em_shdwi*sb*thib^4*vf_ww*(one(FT)-em_sunwi)*vf_wr -
                   FT(3)*em_shdwi*sb*thib^4*vf_wf*(one(FT)-em_floori)*vf_fr -
                   FT(3)*em_floori*sb*tfb^4*vf_fw*(one(FT)-em_sunwi)*vf_wr -
                   FT(3)*em_floori*sb*tfb^4*vf_fw*(one(FT)-em_shdwi)*vf_wr -
                   FT(0.5)*hcv_roofi*(trib - tbb)

        # ---- SUNWALL (equation 2) ----
        a[2,1] = - FT(4)*em_sunwi*em_roofi*sb*trib^3*vf_rw -
                   FT(4)*em_roofi*sb*trib^3*vf_rw*(one(FT)-em_shdwi)*vf_ww -
                   FT(4)*em_roofi*sb*trib^3*vf_rf*(one(FT)-em_floori)*vf_fw
        a[2,2] =   FT(0.5)*hcv_sunwi*building_hwr +
                   FT(0.5)*tks/(zis - zs)*building_hwr +
                   FT(4)*em_sunwi*sb*tsib^3 -
                   FT(4)*em_sunwi*sb*tsib^3*vf_wr*(one(FT)-em_roofi)*vf_rw -
                   FT(4)*em_sunwi*sb*tsib^3*vf_ww*(one(FT)-em_shdwi)*vf_ww -
                   FT(4)*em_sunwi*sb*tsib^3*vf_wf*(one(FT)-em_floori)*vf_fw
        a[2,3] = - FT(4)*em_sunwi*em_shdwi*sb*thib^3*vf_ww -
                   FT(4)*em_shdwi*sb*thib^3*vf_wf*(one(FT)-em_floori)*vf_fw -
                   FT(4)*em_shdwi*sb*thib^3*vf_wr*(one(FT)-em_roofi)*vf_rw
        a[2,4] = - FT(4)*em_sunwi*em_floori*sb*tfb^3*vf_fw -
                   FT(4)*em_floori*sb*tfb^3*vf_fr*(one(FT)-em_roofi)*vf_rw -
                   FT(4)*em_floori*sb*tfb^3*vf_fw*(one(FT)-em_shdwi)*vf_ww
        a[2,5] = - FT(0.5)*hcv_sunwi*building_hwr
        res[2] =   FT(0.5)*tks*tsil/(zis - zs)*building_hwr -
                   FT(0.5)*tks*(tsib - tsilb)/(zis - zs)*building_hwr -
                   FT(3)*em_sunwi*em_roofi*sb*trib^4*vf_rw -
                   FT(3)*em_sunwi*em_shdwi*sb*thib^4*vf_ww -
                   FT(3)*em_sunwi*em_floori*sb*tfb^4*vf_fw +
                   FT(3)*em_sunwi*sb*tsib^4 -
                   FT(3)*em_sunwi*sb*tsib^4*vf_wr*(one(FT)-em_roofi)*vf_rw -
                   FT(3)*em_sunwi*sb*tsib^4*vf_ww*(one(FT)-em_shdwi)*vf_ww -
                   FT(3)*em_sunwi*sb*tsib^4*vf_wf*(one(FT)-em_floori)*vf_fw -
                   FT(3)*em_shdwi*sb*thib^4*vf_wf*(one(FT)-em_floori)*vf_fw -
                   FT(3)*em_shdwi*sb*thib^4*vf_wr*(one(FT)-em_roofi)*vf_rw -
                   FT(3)*em_roofi*sb*trib^4*vf_rw*(one(FT)-em_shdwi)*vf_ww -
                   FT(3)*em_roofi*sb*trib^4*vf_rf*(one(FT)-em_floori)*vf_fw -
                   FT(3)*em_floori*sb*tfb^4*vf_fr*(one(FT)-em_roofi)*vf_rw -
                   FT(3)*em_floori*sb*tfb^4*vf_fw*(one(FT)-em_shdwi)*vf_ww -
                   FT(0.5)*hcv_sunwi*(tsib - tbb)*building_hwr

        # ---- SHADEWALL (equation 3) ----
        a[3,1] = - FT(4)*em_shdwi*em_roofi*sb*trib^3*vf_rw -
                   FT(4)*em_roofi*sb*trib^3*vf_rw*(one(FT)-em_sunwi)*vf_ww -
                   FT(4)*em_roofi*sb*trib^3*vf_rf*(one(FT)-em_floori)*vf_fw
        a[3,2] = - FT(4)*em_shdwi*em_sunwi*sb*tsib^3*vf_ww -
                   FT(4)*em_sunwi*sb*tsib^3*vf_wf*(one(FT)-em_floori)*vf_fw -
                   FT(4)*em_sunwi*sb*tsib^3*vf_wr*(one(FT)-em_roofi)*vf_rw
        a[3,3] =   FT(0.5)*hcv_shdwi*building_hwr +
                   FT(0.5)*tkh/(zih - zh)*building_hwr +
                   FT(4)*em_shdwi*sb*thib^3 -
                   FT(4)*em_shdwi*sb*thib^3*vf_wr*(one(FT)-em_roofi)*vf_rw -
                   FT(4)*em_shdwi*sb*thib^3*vf_ww*(one(FT)-em_sunwi)*vf_ww -
                   FT(4)*em_shdwi*sb*thib^3*vf_wf*(one(FT)-em_floori)*vf_fw
        a[3,4] = - FT(4)*em_shdwi*em_floori*sb*tfb^3*vf_fw -
                   FT(4)*em_floori*sb*tfb^3*vf_fr*(one(FT)-em_roofi)*vf_rw -
                   FT(4)*em_floori*sb*tfb^3*vf_fw*(one(FT)-em_sunwi)*vf_ww
        a[3,5] = - FT(0.5)*hcv_shdwi*building_hwr
        res[3] =   FT(0.5)*tkh*thil/(zih - zh)*building_hwr -
                   FT(0.5)*tkh*(thib - thilb)/(zih - zh)*building_hwr -
                   FT(3)*em_shdwi*em_roofi*sb*trib^4*vf_rw -
                   FT(3)*em_shdwi*em_sunwi*sb*tsib^4*vf_ww -
                   FT(3)*em_shdwi*em_floori*sb*tfb^4*vf_fw +
                   FT(3)*em_shdwi*sb*thib^4 -
                   FT(3)*em_shdwi*sb*thib^4*vf_wr*(one(FT)-em_roofi)*vf_rw -
                   FT(3)*em_shdwi*sb*thib^4*vf_ww*(one(FT)-em_sunwi)*vf_ww -
                   FT(3)*em_shdwi*sb*thib^4*vf_wf*(one(FT)-em_floori)*vf_fw -
                   FT(3)*em_sunwi*sb*tsib^4*vf_wf*(one(FT)-em_floori)*vf_fw -
                   FT(3)*em_sunwi*sb*tsib^4*vf_wr*(one(FT)-em_roofi)*vf_rw -
                   FT(3)*em_roofi*sb*trib^4*vf_rw*(one(FT)-em_sunwi)*vf_ww -
                   FT(3)*em_roofi*sb*trib^4*vf_rf*(one(FT)-em_floori)*vf_fw -
                   FT(3)*em_floori*sb*tfb^4*vf_fr*(one(FT)-em_roofi)*vf_rw -
                   FT(3)*em_floori*sb*tfb^4*vf_fw*(one(FT)-em_sunwi)*vf_ww -
                   FT(0.5)*hcv_shdwi*(thib - tbb)*building_hwr

        # ---- FLOOR (equation 4) ----
        a[4,1] = - FT(4)*em_floori*em_roofi*sb*trib^3*vf_rf -
                   FT(4)*em_roofi*sb*trib^3*vf_rw*(one(FT)-em_sunwi)*vf_wf -
                   FT(4)*em_roofi*sb*trib^3*vf_rw*(one(FT)-em_shdwi)*vf_wf
        a[4,2] = - FT(4)*em_floori*em_sunwi*sb*tsib^3*vf_wf -
                   FT(4)*em_sunwi*sb*tsib^3*vf_ww*(one(FT)-em_shdwi)*vf_wf -
                   FT(4)*em_sunwi*sb*tsib^3*vf_wr*(one(FT)-em_roofi)*vf_rf
        a[4,3] = - FT(4)*em_floori*em_shdwi*sb*thib^3*vf_wf -
                   FT(4)*em_shdwi*sb*thib^3*vf_wr*(one(FT)-em_roofi)*vf_rf -
                   FT(4)*em_shdwi*sb*thib^3*vf_ww*(one(FT)-em_sunwi)*vf_wf
        a[4,4] =   (cv_floori + FT(0.5)*hcv_floori) +
                   FT(4)*em_floori*sb*tfb^3 -
                   FT(4)*em_floori*sb*tfb^3*vf_fr*(one(FT)-em_roofi)*vf_rf -
                   FT(4)*em_floori*sb*tfb^3*vf_fw*(one(FT)-em_sunwi)*vf_wf -
                   FT(4)*em_floori*sb*tfb^3*vf_fw*(one(FT)-em_shdwi)*vf_wf
        a[4,5] = - FT(0.5)*hcv_floori
        res[4] =   cv_floori*tfb -
                   FT(3)*em_floori*em_roofi*sb*trib^4*vf_rf -
                   FT(3)*em_floori*em_sunwi*sb*tsib^4*vf_wf -
                   FT(3)*em_floori*em_shdwi*sb*thib^4*vf_wf +
                   FT(3)*em_floori*sb*tfb^4 -
                   FT(3)*em_floori*sb*tfb^4*vf_fr*(one(FT)-em_roofi)*vf_rf -
                   FT(3)*em_floori*sb*tfb^4*vf_fw*(one(FT)-em_sunwi)*vf_wf -
                   FT(3)*em_floori*sb*tfb^4*vf_fw*(one(FT)-em_shdwi)*vf_wf -
                   FT(3)*em_sunwi*sb*tsib^4*vf_ww*(one(FT)-em_shdwi)*vf_wf -
                   FT(3)*em_sunwi*sb*tsib^4*vf_wr*(one(FT)-em_roofi)*vf_rf -
                   FT(3)*em_shdwi*sb*thib^4*vf_wr*(one(FT)-em_roofi)*vf_rf -
                   FT(3)*em_shdwi*sb*thib^4*vf_ww*(one(FT)-em_sunwi)*vf_wf -
                   FT(3)*em_roofi*sb*trib^4*vf_rw*(one(FT)-em_sunwi)*vf_wf -
                   FT(3)*em_roofi*sb*trib^4*vf_rw*(one(FT)-em_shdwi)*vf_wf -
                   FT(0.5)*hcv_floori*(tfb - tbb)

        # ---- BUILDING AIR (equation 5) ----
        a[5,1] = - FT(0.5)*hcv_roofi
        a[5,2] = - FT(0.5)*hcv_sunwi*building_hwr
        a[5,3] = - FT(0.5)*hcv_shdwi*building_hwr
        a[5,4] = - FT(0.5)*hcv_floori
        a[5,5] =   ((ht_roof[l]*rho_dair*cpair)/dt) +
                   ((ht_roof[l]*vent_ach)/FT(3600))*rho_dair*cpair +
                   FT(0.5)*hcv_roofi +
                   FT(0.5)*hcv_sunwi*building_hwr +
                   FT(0.5)*hcv_shdwi*building_hwr +
                   FT(0.5)*hcv_floori
        res[5] =   (ht_roof[l]*rho_dair*cpair/dt)*tbb +
                   ((ht_roof[l]*vent_ach)/FT(3600))*rho_dair*cpair*taf[l] +
                   FT(0.5)*hcv_roofi*(trib - tbb) +
                   FT(0.5)*hcv_sunwi*(tsib - tbb)*building_hwr +
                   FT(0.5)*hcv_shdwi*(thib - tbb)*building_hwr +
                   FT(0.5)*hcv_floori*(tfb - tbb)

        # Solve the 5x5 system (res holds the solution on return)
        _solve5x5!(a, res)

        t_roof_inner[l] = res[1]
        t_sunw_inner[l] = res[2]
        t_shdw_inner[l] = res[3]
        t_floor[l]      = res[4]
        t_building[l]   = res[5]

        # Sensible heat flux from ventilation (W/m2 urban area via wtlunit_roof)
        eflx_ventilation[l] = wtlunit_roof[l] * ( - ht_roof[l]*(vent_ach/FT(3600)) *
                              rho_dair * cpair * (taf[l] - t_building[l]) )

        # ------------------------------------------------------------------
        # Restrict building air temperature to [t_building_min, t_building_max]
        # and compute heating / air-conditioning flux.
        # ------------------------------------------------------------------
        if hac_active
            t_building_bef_hac = t_building[l]
            if t_building_bef_hac > t_building_max[l]
                if urban_explicit_ac
                    eflx_urban_ac_sat = wtlunit_roof[l] * abs(
                        (ht_roof[l]*rho_dair*cpair/dt)*t_building_max[l] -
                        (ht_roof[l]*rho_dair*cpair/dt)*t_building_bef_hac )
                    t_building[l] = t_building_max[l] + (one(FT) - p_ac[l]) * eflx_urban_ac_sat *
                                    dt / (ht_roof[l]*rho_dair*cpair*wtlunit_roof[l])
                    eflx_urban_ac[l] = p_ac[l] * eflx_urban_ac_sat
                else
                    t_building[l] = t_building_max[l]
                    eflx_urban_ac[l] = wtlunit_roof[l] * abs(
                        (ht_roof[l]*rho_dair*cpair/dt)*t_building[l] -
                        (ht_roof[l]*rho_dair*cpair/dt)*t_building_bef_hac )
                end
            elseif t_building_bef_hac < t_building_min[l]
                t_building[l] = t_building_min[l]
                eflx_urban_heat[l] = wtlunit_roof[l] * abs(
                    (ht_roof[l]*rho_dair*cpair/dt)*t_building[l] -
                    (ht_roof[l]*rho_dair*cpair/dt)*t_building_bef_hac )
            else
                eflx_urban_ac[l]   = zero(FT)
                eflx_urban_heat[l] = zero(FT)
            end
        else
            eflx_urban_ac[l]   = zero(FT)
            eflx_urban_heat[l] = zero(FT)
        end

        eflx_building[l] = wtlunit_roof[l] * (ht_roof[l]*rho_dair*cpair/dt) *
                           (t_building[l] - t_building_bef)
    end

    return nothing
end
