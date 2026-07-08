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

# ==========================================================================
# GPU KERNELIZATION (KernelAbstractions, runs on host + Metal at parity).
#
# The routine has two data-parallel phases:
#   1. GATHER (one thread per COLUMN): each urban roof/sunwall/shadewall column
#      writes the inner-node (nlevurb) geometry/temps/conductivity into DISTINCT
#      per-landunit slots (roof- OR sunwall- OR shadewall-fields, selected by
#      col.itype). Each urban landunit has exactly one column of each type, so the
#      three columns touch disjoint fields of the same landunit — race-free, no atomics.
#   2. SOLVE (one thread per LANDUNIT): assemble the 5×5 energy-balance system into
#      per-thread scratch slices `a_scr[l,:,:]`/`res_scr[l,:]` (a backend array,
#      allocation-free in-kernel — same device-scratch trick as the batched band
#      solver), run partial-pivoted Gaussian elimination, then the HAC clamp + fluxes.
#
# Field bundles keep the launch under Metal's ~31-buffer limit; a `_BTScal{FT}` isbits
# struct carries the constants so no `Float64` scalar reaches a Metal kernel. The KA
# CPU backend runs the identical kernels on `Array`s (byte-identical) and on
# `ForwardDiff.Dual` arrays (the AD path flows through the pivoted GE).
# ==========================================================================

# Per-landunit gather arrays (inner-node values from the 3 wall/roof columns).
struct _BTGath{V}
    zi_roof::V; z_roof::V; tk_roof::V; t_roof::V; tbef_roof::V
    zi_sunw::V; z_sunw::V; tk_sunw::V; t_sunw::V; tbef_sunw::V
    zi_shdw::V; z_shdw::V; tk_shdw::V; t_shdw::V; tbef_shdw::V
end
Adapt.@adapt_structure _BTGath

# Landunit-indexed inputs to the solve (V float vectors, VB the Bool urbpoi mask).
struct _BTLunIn{V,VB}
    canyon_hwr::V; wtlunit_roof::V; ht_roof::V; taf::V
    t_building_min::V; t_building_max::V; p_ac::V
    urbpoi::VB
end
Adapt.@adapt_structure _BTLunIn

# Landunit-indexed outputs mutated by the solve.
struct _BTOut{V}
    t_roof_inner::V; t_sunw_inner::V; t_shdw_inner::V; t_floor::V; t_building::V
    eflx_building::V; eflx_urban_ac::V; eflx_urban_heat::V; eflx_ventilation::V
end
Adapt.@adapt_structure _BTOut

# Working-type scalar constants (isbits — passed by value, no Float64 on Metal).
struct _BTScal{FT}
    sb::FT; rair::FT; pstd::FT; cpair::FT; dt::FT
    em_roofi::FT; em_sunwi::FT; em_shdwi::FT; em_floori::FT
    dz_floori::FT; cp_floori::FT; vent_ach::FT
    hcv_roof::FT; hcv_roof_enh::FT; hcv_floor::FT; hcv_floor_enh::FT
    hcv_sunw::FT; hcv_shdw::FT
end

# --------------------------------------------------------------------------
# GATHER kernel — one thread per column; dispatch by column type.
# --------------------------------------------------------------------------
@kernel function _bt_gather_kernel!(g, @Const(col_zi), @Const(col_z), @Const(t_soisno),
        @Const(tssbef), @Const(tk), @Const(col_landunit), @Const(col_itype),
        @Const(urbpoi), @Const(mask_nolakec), nlevurb_jj::Int)
    c = @index(Global)
    @inbounds if mask_nolakec[c]
        l = col_landunit[c]
        if urbpoi[l]
            it = col_itype[c]
            if it == ICOL_ROOF
                g.zi_roof[l]   = col_zi[c, nlevurb_jj + 1]
                g.z_roof[l]    = col_z[c, nlevurb_jj]
                g.tbef_roof[l] = tssbef[c, nlevurb_jj]
                g.t_roof[l]    = t_soisno[c, nlevurb_jj]
                g.tk_roof[l]   = tk[c, nlevurb_jj]
            elseif it == ICOL_SUNWALL
                g.zi_sunw[l]   = col_zi[c, nlevurb_jj + 1]
                g.z_sunw[l]    = col_z[c, nlevurb_jj]
                g.tbef_sunw[l] = tssbef[c, nlevurb_jj]
                g.t_sunw[l]    = t_soisno[c, nlevurb_jj]
                g.tk_sunw[l]   = tk[c, nlevurb_jj]
            elseif it == ICOL_SHADEWALL
                g.zi_shdw[l]   = col_zi[c, nlevurb_jj + 1]
                g.z_shdw[l]    = col_z[c, nlevurb_jj]
                g.tbef_shdw[l] = tssbef[c, nlevurb_jj]
                g.t_shdw[l]    = t_soisno[c, nlevurb_jj]
                g.tk_shdw[l]   = tk[c, nlevurb_jj]
            end
        end
    end
end

# --------------------------------------------------------------------------
# SOLVE kernel — one thread per landunit; assemble + pivoted GE + HAC clamp.
# --------------------------------------------------------------------------
@kernel function _bt_solve_kernel!(o, lin, g, a, res, @Const(mask_urbanl), s,
        hac_active::Bool, urban_explicit_ac::Bool)
    l = @index(Global)
    @inbounds if mask_urbanl[l] && lin.urbpoi[l]
        FT = eltype(res)

        # 1. Save previous-step temperatures
        t_roof_inner_bef = o.t_roof_inner[l]
        t_sunw_inner_bef = o.t_sunw_inner[l]
        t_shdw_inner_bef = o.t_shdw_inner[l]
        t_floor_bef      = o.t_floor[l]
        t_building_bef   = o.t_building[l]

        # 2. Convective heat transfer coefficients (Bueno et al. 2012)
        hcv_roofi  = t_roof_inner_bef <= t_building_bef ? s.hcv_roof_enh  : s.hcv_roof
        hcv_floori = t_floor_bef      >= t_building_bef ? s.hcv_floor_enh : s.hcv_floor
        hcv_sunwi  = s.hcv_sunw
        hcv_shdwi  = s.hcv_shdw

        sb = s.sb; cpair = s.cpair; dt = s.dt
        em_roofi  = s.em_roofi
        em_sunwi  = s.em_sunwi
        em_shdwi  = s.em_shdwi
        em_floori = s.em_floori

        # Intermediate calc for concrete floor (W m-2 K-1)
        cv_floori = (s.dz_floori * s.cp_floori) / dt
        # Density of dry air at standard pressure and t_building
        rho_dair = s.pstd / (s.rair * t_building_bef)
        # Building height to building width ratio
        ht_roof_l    = lin.ht_roof[l]
        wtlunit_roof = lin.wtlunit_roof[l]
        building_hwr = lin.canyon_hwr[l] * (one(FT) - wtlunit_roof) / wtlunit_roof

        # 3. View factors (inner box)
        vf_rf = sqrt(one(FT) + building_hwr^2) - building_hwr
        vf_fr = vf_rf
        vf_wf = FT(0.5) * (one(FT) - vf_rf)   # per unit wall -> per unit floor
        vf_fw = vf_wf
        vf_rw = vf_fw
        vf_wr = vf_wf                          # per unit wall -> per unit roof
        vf_ww = one(FT) - vf_rw - vf_fw

        # Geometry / conductivity at the inner (nlevurb) node
        zir = g.zi_roof[l]; zr = g.z_roof[l]; tkr = g.tk_roof[l]
        zis = g.zi_sunw[l]; zs = g.z_sunw[l]; tks = g.tk_sunw[l]
        zih = g.zi_shdw[l]; zh = g.z_shdw[l]; tkh = g.tk_shdw[l]
        tril = g.t_roof[l];  trilb = g.tbef_roof[l]
        tsil = g.t_sunw[l];  tsilb = g.tbef_sunw[l]
        thil = g.t_shdw[l];  thilb = g.tbef_shdw[l]

        trib = t_roof_inner_bef
        tsib = t_sunw_inner_bef
        thib = t_shdw_inner_bef
        tfb  = t_floor_bef
        tbb  = t_building_bef

        # ---- ROOF (equation 1) ----
        a[l,1,1] =   FT(0.5)*hcv_roofi +
                   FT(0.5)*tkr/(zir - zr) +
                   FT(4)*em_roofi*sb*trib^3 -
                   FT(4)*em_roofi*sb*trib^3*vf_rw*(one(FT)-em_sunwi)*vf_wr -
                   FT(4)*em_roofi*sb*trib^3*vf_rw*(one(FT)-em_shdwi)*vf_wr -
                   FT(4)*em_roofi*sb*trib^3*vf_rf*(one(FT)-em_floori)*vf_fr
        a[l,1,2] = - FT(4)*em_roofi*em_sunwi*sb*tsib^3*vf_wr -
                   FT(4)*em_sunwi*sb*tsib^3*vf_ww*(one(FT)-em_shdwi)*vf_wr -
                   FT(4)*em_sunwi*sb*tsib^3*vf_wf*(one(FT)-em_floori)*vf_fr
        a[l,1,3] = - FT(4)*em_roofi*em_shdwi*sb*thib^3*vf_wr -
                   FT(4)*em_shdwi*sb*thib^3*vf_ww*(one(FT)-em_sunwi)*vf_wr -
                   FT(4)*em_shdwi*sb*thib^3*vf_wf*(one(FT)-em_floori)*vf_fr
        a[l,1,4] = - FT(4)*em_roofi*em_floori*sb*tfb^3*vf_fr -
                   FT(4)*em_floori*sb*tfb^3*vf_fw*(one(FT)-em_sunwi)*vf_wr -
                   FT(4)*em_floori*sb*tfb^3*vf_fw*(one(FT)-em_shdwi)*vf_wr
        a[l,1,5] = - FT(0.5)*hcv_roofi
        res[l,1] =   FT(0.5)*tkr*tril/(zir - zr) -
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
        a[l,2,1] = - FT(4)*em_sunwi*em_roofi*sb*trib^3*vf_rw -
                   FT(4)*em_roofi*sb*trib^3*vf_rw*(one(FT)-em_shdwi)*vf_ww -
                   FT(4)*em_roofi*sb*trib^3*vf_rf*(one(FT)-em_floori)*vf_fw
        a[l,2,2] =   FT(0.5)*hcv_sunwi*building_hwr +
                   FT(0.5)*tks/(zis - zs)*building_hwr +
                   FT(4)*em_sunwi*sb*tsib^3 -
                   FT(4)*em_sunwi*sb*tsib^3*vf_wr*(one(FT)-em_roofi)*vf_rw -
                   FT(4)*em_sunwi*sb*tsib^3*vf_ww*(one(FT)-em_shdwi)*vf_ww -
                   FT(4)*em_sunwi*sb*tsib^3*vf_wf*(one(FT)-em_floori)*vf_fw
        a[l,2,3] = - FT(4)*em_sunwi*em_shdwi*sb*thib^3*vf_ww -
                   FT(4)*em_shdwi*sb*thib^3*vf_wf*(one(FT)-em_floori)*vf_fw -
                   FT(4)*em_shdwi*sb*thib^3*vf_wr*(one(FT)-em_roofi)*vf_rw
        a[l,2,4] = - FT(4)*em_sunwi*em_floori*sb*tfb^3*vf_fw -
                   FT(4)*em_floori*sb*tfb^3*vf_fr*(one(FT)-em_roofi)*vf_rw -
                   FT(4)*em_floori*sb*tfb^3*vf_fw*(one(FT)-em_shdwi)*vf_ww
        a[l,2,5] = - FT(0.5)*hcv_sunwi*building_hwr
        res[l,2] =   FT(0.5)*tks*tsil/(zis - zs)*building_hwr -
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
        a[l,3,1] = - FT(4)*em_shdwi*em_roofi*sb*trib^3*vf_rw -
                   FT(4)*em_roofi*sb*trib^3*vf_rw*(one(FT)-em_sunwi)*vf_ww -
                   FT(4)*em_roofi*sb*trib^3*vf_rf*(one(FT)-em_floori)*vf_fw
        a[l,3,2] = - FT(4)*em_shdwi*em_sunwi*sb*tsib^3*vf_ww -
                   FT(4)*em_sunwi*sb*tsib^3*vf_wf*(one(FT)-em_floori)*vf_fw -
                   FT(4)*em_sunwi*sb*tsib^3*vf_wr*(one(FT)-em_roofi)*vf_rw
        a[l,3,3] =   FT(0.5)*hcv_shdwi*building_hwr +
                   FT(0.5)*tkh/(zih - zh)*building_hwr +
                   FT(4)*em_shdwi*sb*thib^3 -
                   FT(4)*em_shdwi*sb*thib^3*vf_wr*(one(FT)-em_roofi)*vf_rw -
                   FT(4)*em_shdwi*sb*thib^3*vf_ww*(one(FT)-em_sunwi)*vf_ww -
                   FT(4)*em_shdwi*sb*thib^3*vf_wf*(one(FT)-em_floori)*vf_fw
        a[l,3,4] = - FT(4)*em_shdwi*em_floori*sb*tfb^3*vf_fw -
                   FT(4)*em_floori*sb*tfb^3*vf_fr*(one(FT)-em_roofi)*vf_rw -
                   FT(4)*em_floori*sb*tfb^3*vf_fw*(one(FT)-em_sunwi)*vf_ww
        a[l,3,5] = - FT(0.5)*hcv_shdwi*building_hwr
        res[l,3] =   FT(0.5)*tkh*thil/(zih - zh)*building_hwr -
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
        a[l,4,1] = - FT(4)*em_floori*em_roofi*sb*trib^3*vf_rf -
                   FT(4)*em_roofi*sb*trib^3*vf_rw*(one(FT)-em_sunwi)*vf_wf -
                   FT(4)*em_roofi*sb*trib^3*vf_rw*(one(FT)-em_shdwi)*vf_wf
        a[l,4,2] = - FT(4)*em_floori*em_sunwi*sb*tsib^3*vf_wf -
                   FT(4)*em_sunwi*sb*tsib^3*vf_ww*(one(FT)-em_shdwi)*vf_wf -
                   FT(4)*em_sunwi*sb*tsib^3*vf_wr*(one(FT)-em_roofi)*vf_rf
        a[l,4,3] = - FT(4)*em_floori*em_shdwi*sb*thib^3*vf_wf -
                   FT(4)*em_shdwi*sb*thib^3*vf_wr*(one(FT)-em_roofi)*vf_rf -
                   FT(4)*em_shdwi*sb*thib^3*vf_ww*(one(FT)-em_sunwi)*vf_wf
        a[l,4,4] =   (cv_floori + FT(0.5)*hcv_floori) +
                   FT(4)*em_floori*sb*tfb^3 -
                   FT(4)*em_floori*sb*tfb^3*vf_fr*(one(FT)-em_roofi)*vf_rf -
                   FT(4)*em_floori*sb*tfb^3*vf_fw*(one(FT)-em_sunwi)*vf_wf -
                   FT(4)*em_floori*sb*tfb^3*vf_fw*(one(FT)-em_shdwi)*vf_wf
        a[l,4,5] = - FT(0.5)*hcv_floori
        res[l,4] =   cv_floori*tfb -
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
        a[l,5,1] = - FT(0.5)*hcv_roofi
        a[l,5,2] = - FT(0.5)*hcv_sunwi*building_hwr
        a[l,5,3] = - FT(0.5)*hcv_shdwi*building_hwr
        a[l,5,4] = - FT(0.5)*hcv_floori
        a[l,5,5] =   ((ht_roof_l*rho_dair*cpair)/dt) +
                   ((ht_roof_l*s.vent_ach)/FT(3600))*rho_dair*cpair +
                   FT(0.5)*hcv_roofi +
                   FT(0.5)*hcv_sunwi*building_hwr +
                   FT(0.5)*hcv_shdwi*building_hwr +
                   FT(0.5)*hcv_floori
        res[l,5] =   (ht_roof_l*rho_dair*cpair/dt)*tbb +
                   ((ht_roof_l*s.vent_ach)/FT(3600))*rho_dair*cpair*lin.taf[l] +
                   FT(0.5)*hcv_roofi*(trib - tbb) +
                   FT(0.5)*hcv_sunwi*(tsib - tbb)*building_hwr +
                   FT(0.5)*hcv_shdwi*(thib - tbb)*building_hwr +
                   FT(0.5)*hcv_floori*(tfb - tbb)

        # Solve the 5x5 system in the per-landunit scratch slice by partial-pivoted
        # Gaussian elimination (mirrors _solve5x5! / the batched band solver). NOTE:
        # `akk` is captured as a[l,piv,k] BEFORE the row swap (== the post-swap
        # diagonal a[l,k,k]); reading a[l,k,k] after the swap miscompiles under the
        # KA CPU backend (the pivot-search read is cached across the swap). Urban
        # roof/wall rows have genuinely different magnitudes, so pivots DO swap here.
        for k in 1:5
            piv = k
            amax = abs(a[l,k,k])
            for i in (k+1):5
                ai = abs(a[l,i,k])
                if ai > amax
                    amax = ai
                    piv = i
                end
            end
            akk = a[l,piv,k]
            if piv != k
                for j in 1:5
                    akj = a[l,k,j]
                    apj = a[l,piv,j]
                    a[l,k,j] = apj
                    a[l,piv,j] = akj
                end
                rk = res[l,k]; rp = res[l,piv]
                res[l,k] = rp; res[l,piv] = rk
            end
            for i in (k+1):5
                f = a[l,i,k] / akk
                for j in (k+1):5
                    a[l,i,j] -= f * a[l,k,j]
                end
                res[l,i] -= f * res[l,k]
            end
        end
        for k in 5:-1:1
            sm = res[l,k]
            for j in (k+1):5
                sm -= a[l,k,j] * res[l,j]
            end
            res[l,k] = sm / a[l,k,k]
        end

        o.t_roof_inner[l] = res[l,1]
        o.t_sunw_inner[l] = res[l,2]
        o.t_shdw_inner[l] = res[l,3]
        o.t_floor[l]      = res[l,4]
        o.t_building[l]   = res[l,5]

        # Sensible heat flux from ventilation (W/m2 urban area via wtlunit_roof)
        o.eflx_ventilation[l] = wtlunit_roof * ( - ht_roof_l*(s.vent_ach/FT(3600)) *
                              rho_dair * cpair * (lin.taf[l] - o.t_building[l]) )

        # ------------------------------------------------------------------
        # Restrict building air temperature to [t_building_min, t_building_max]
        # and compute heating / air-conditioning flux.
        # ------------------------------------------------------------------
        if hac_active
            t_building_bef_hac = o.t_building[l]
            if t_building_bef_hac > lin.t_building_max[l]
                if urban_explicit_ac
                    eflx_urban_ac_sat = wtlunit_roof * abs(
                        (ht_roof_l*rho_dair*cpair/dt)*lin.t_building_max[l] -
                        (ht_roof_l*rho_dair*cpair/dt)*t_building_bef_hac )
                    o.t_building[l] = lin.t_building_max[l] + (one(FT) - lin.p_ac[l]) * eflx_urban_ac_sat *
                                    dt / (ht_roof_l*rho_dair*cpair*wtlunit_roof)
                    o.eflx_urban_ac[l] = lin.p_ac[l] * eflx_urban_ac_sat
                else
                    o.t_building[l] = lin.t_building_max[l]
                    o.eflx_urban_ac[l] = wtlunit_roof * abs(
                        (ht_roof_l*rho_dair*cpair/dt)*o.t_building[l] -
                        (ht_roof_l*rho_dair*cpair/dt)*t_building_bef_hac )
                end
            elseif t_building_bef_hac < lin.t_building_min[l]
                o.t_building[l] = lin.t_building_min[l]
                o.eflx_urban_heat[l] = wtlunit_roof * abs(
                    (ht_roof_l*rho_dair*cpair/dt)*o.t_building[l] -
                    (ht_roof_l*rho_dair*cpair/dt)*t_building_bef_hac )
            else
                o.eflx_urban_ac[l]   = zero(FT)
                o.eflx_urban_heat[l] = zero(FT)
            end
        else
            o.eflx_urban_ac[l]   = zero(FT)
            o.eflx_urban_heat[l] = zero(FT)
        end

        o.eflx_building[l] = wtlunit_roof * (ht_roof_l*rho_dair*cpair/dt) *
                           (o.t_building[l] - t_building_bef)
    end
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

Ported from `BuildingTemperature` in `UrbBuildTempOleson2015Mod.F90`. Runs as two
KernelAbstractions kernels (gather over columns, solve over landunits) — backend-
agnostic: the KA CPU backend is byte-identical to the old scalar loop and carries
`ForwardDiff.Dual` for AD; device arrays run on the GPU (Metal validated).
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

    # State aliases
    t_soisno     = temperature.t_soisno_col
    tssbef       = temperature.t_ssbef_col

    # HAC mode / explicit AC flags (resolved on host — passed as kernel scalars)
    urban_hac         = urban_ctrl.urban_hac
    urban_explicit_ac = urban_ctrl.urban_explicit_ac
    hac_active = (urban_hac == URBAN_HAC_ON) || (urban_hac == URBAN_WASTEHEAT_ON)

    # Inner-node (nlevurb) Julia level index = nlevurb + joff. The corresponding
    # lower interface is zi[c, nlevurb + joff + 1] (matches _fn1_kernel!).
    nlevurb_jj = nlevurb + joff

    # Working-type scalar constants bundle.
    scal = _BTScal{FT}(FT(SB), FT(RAIR), FT(PSTD), FT(CPAIR), FT(dtime),
        FT(EM_ROOF_INT), FT(EM_SUNW_INT), FT(EM_SHDW_INT), FT(EM_FLOOR_INT),
        FT(DZ_FLOOR), FT(CP_FLOOR), FT(VENT_ACH),
        FT(HCV_ROOF), FT(HCV_ROOF_ENHANCED), FT(HCV_FLOOR), FT(HCV_FLOOR_ENHANCED),
        FT(HCV_SUNW), FT(HCV_SHDW))

    # ----------------------------------------------------------------------
    # GATHER (one thread per column). Per-landunit inner-node arrays allocated on
    # the state backend so the kernel writes device-resident memory. Sized to the
    # landunit count; NaN-initialized (unused non-urban slots never read).
    # ----------------------------------------------------------------------
    nl = length(mask_urbanl)
    mk() = fill!(similar(t_soisno, FT, nl), FT(NaN))
    gath = _BTGath(mk(), mk(), mk(), mk(), mk(), mk(), mk(), mk(),
                   mk(), mk(), mk(), mk(), mk(), mk(), mk())

    be = _kernel_backend(t_soisno)
    _bt_gather_kernel!(be)(gath, col.zi, col.z, t_soisno, tssbef, tk,
        col.landunit, col.itype, lun.urbpoi, mask_nolakec, nlevurb_jj;
        ndrange = length(mask_nolakec))
    KA.synchronize(be)

    # ----------------------------------------------------------------------
    # SOLVE (one thread per landunit). Per-thread 5x5 scratch slices on the backend.
    # ----------------------------------------------------------------------
    a_scr   = similar(t_soisno, FT, nl, 5, 5)
    res_scr = similar(t_soisno, FT, nl, 5)

    lin = _BTLunIn(lun.canyon_hwr, lun.wtlunit_roof, lun.ht_roof, temperature.taf_lun,
        urbanparams.t_building_min, urbantv_t_building_max, urbantv_p_ac, lun.urbpoi)
    out = _BTOut(temperature.t_roof_inner_lun, temperature.t_sunw_inner_lun,
        temperature.t_shdw_inner_lun, temperature.t_floor_lun, temperature.t_building_lun,
        energyflux.eflx_building_lun, energyflux.eflx_urban_ac_lun,
        energyflux.eflx_urban_heat_lun, energyflux.eflx_ventilation_lun)

    _bt_solve_kernel!(be)(out, lin, gath, a_scr, res_scr, mask_urbanl, scal,
        hac_active, urban_explicit_ac; ndrange = nl)
    KA.synchronize(be)

    return nothing
end
