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

# ===========================================================================
# GPU kernels (KernelAbstractions) for fully per-element soil-water loops.
# Each (column, layer) element is computed from inputs at its OWN index
# (with constant snow offsets); no neighbor-layer reads, no accumulation,
# no tridiagonal coupling. See src/infrastructure/kernels.jl for the API.
# ===========================================================================

# --------------------------------------------------------------------------
# Convert geometry to mm and compute per-layer ice/liquid volume fractions
# (Zeng-Decker 2009 prep). 2D over (column, soil layer); masked in-kernel.
# Writes zmm/dzmm at [c,j], zimm_arr at [c,j+1] (unique per j), and
# vol_ice/icefrac/vwc_liq at [c,j]. All reads are at the element's own index.
# --------------------------------------------------------------------------
@kernel function _soilwm_mm_icefrac_kernel!(zmm, dzmm, zimm_arr, vol_ice, icefrac,
        vwc_liq, @Const(mask), @Const(z), @Const(dz), @Const(zi),
        @Const(watsat), @Const(h2osoi_ice), @Const(h2osoi_liq),
        joff::Int, joff_zi::Int, denice, denh2o)
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        T = eltype(zmm)
        zmm[c, j]  = z[c, joff + j] * T(1.0e3)
        dzmm[c, j] = dz[c, joff + j] * T(1.0e3)
        zimm_arr[c, j+1] = zi[c, joff_zi + j] * T(1.0e3)  # j+1: index 1 = Fortran j=0

        vol_ice[c, j] = smooth_min(watsat[c, j], h2osoi_ice[c, joff + j] / (dz[c, joff + j] * denice))
        icefrac[c, j] = smooth_min(one(T), vol_ice[c, j] / watsat[c, j])
        vwc_liq[c, j] = smooth_max(h2osoi_liq[c, joff + j], T(1.0e-6)) / (dz[c, joff + j] * denh2o)
    end
end

"""
    soilwm_mm_icefrac!(zmm, dzmm, zimm_arr, vol_ice, icefrac, vwc_liq, mask,
        z, dz, zi, watsat, h2osoi_ice, h2osoi_liq, joff, joff_zi, nlevsoi; denice, denh2o)

Convert column geometry to mm and compute per-layer ice/liquid volume
fractions over all active (column, soil layer) pairs. Backend-agnostic 2D
kernel; replaces the inline double loop in `soilwater_zengdecker2009!`.
"""
function soilwm_mm_icefrac!(zmm, dzmm, zimm_arr, vol_ice, icefrac, vwc_liq,
        mask, z, dz, zi, watsat, h2osoi_ice, h2osoi_liq,
        joff::Int, joff_zi::Int, nlevsoi::Int; denice, denh2o)
    T = eltype(zmm)
    _launch!(_soilwm_mm_icefrac_kernel!, zmm, dzmm, zimm_arr, vol_ice, icefrac,
             vwc_liq, mask, z, dz, zi, watsat, h2osoi_ice, h2osoi_liq,
             joff, joff_zi, convert(T, denice), convert(T, denh2o);
             ndrange = (length(mask), nlevsoi))
end

# --------------------------------------------------------------------------
# Water-table depth to mm and the top zimm entry (Zeng-Decker 2009 prep).
# 1D over columns; masked in-kernel. Each column writes only its own entries.
# --------------------------------------------------------------------------
@kernel function _soilwm_zwtmm_kernel!(zwtmm, zimm_arr, @Const(mask), @Const(zwt))
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(zwtmm)
        zimm_arr[c, 1] = zero(T)       # Fortran zimm(c,0) = 0
        zwtmm[c] = zwt[c] * T(1.0e3)
    end
end

"""
    soilwm_zwtmm!(zwtmm, zimm_arr, mask, zwt)

Set water-table depth in mm and the top zimm entry for each active column.
Backend-agnostic 1D kernel; replaces an inline per-column loop in
`soilwater_zengdecker2009!`.
"""
soilwm_zwtmm!(zwtmm, zimm_arr, mask, zwt) =
    _launch!(_soilwm_zwtmm_kernel!, zwtmm, zimm_arr, mask, zwt)

# --------------------------------------------------------------------------
# Equilibrium volumetric water content vol_eq and equilibrium matric
# potential zq based on water-table depth (Zeng-Decker 2009). 2D over
# (column, soil layer); masked in-kernel. Each (c,j) reads its own layer plus
# the precomputed interface depths zimm_arr[c,j] (Fortran zimm(c,j-1)) and
# zimm_arr[c,j+1] (Fortran zimm(c,j)) — both read-only — so writes are
# disjoint and order-independent. Literals are eltype-converted so no Float64
# reaches a Float32-only backend (byte-identical on Float64).
# --------------------------------------------------------------------------
@kernel function _soilwm_voleq_zq_kernel!(vol_eq, zq, @Const(mask), @Const(zwtmm),
        @Const(zimm_arr), @Const(watsat), @Const(sucsat), @Const(bsw), @Const(smpmin))
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        T = eltype(vol_eq)
        zimm_jm1 = zimm_arr[c, j]      # Fortran zimm(c,j-1)
        zimm_j   = zimm_arr[c, j+1]    # Fortran zimm(c,j)
        ws   = watsat[c, j]
        suc  = sucsat[c, j]
        b    = bsw[c, j]
        zwt  = zwtmm[c]
        onem1b = one(T) - one(T) / b   # (1 - 1/bsw); matches Fortran exactly on Float64
        if zwt <= zimm_jm1
            vol_eq[c, j] = ws
        elseif zwt < zimm_j && zwt > zimm_jm1
            tempi  = one(T)
            temp0  = ((suc + zwt - zimm_jm1) / suc)^onem1b
            voleq1 = -suc * ws / onem1b / (zwt - zimm_jm1) * (tempi - temp0)
            v = (voleq1 * (zwt - zimm_jm1) + ws * (zimm_j - zwt)) / (zimm_j - zimm_jm1)
            v = smooth_min(ws, v)
            vol_eq[c, j] = smooth_max(v, zero(T))
        else
            tempi = ((suc + zwt - zimm_j) / suc)^onem1b
            temp0 = ((suc + zwt - zimm_jm1) / suc)^onem1b
            v = -suc * ws / onem1b / (zimm_j - zimm_jm1) * (tempi - temp0)
            v = smooth_max(v, zero(T))
            vol_eq[c, j] = smooth_min(ws, v)
        end
        zqv = -suc * (smooth_max(vol_eq[c, j] / ws, T(0.01)))^(-b)
        zq[c, j] = smooth_max(smpmin[c], zqv)
    end
end

"""
    soilwm_voleq_zq!(vol_eq, zq, mask, zwtmm, zimm_arr, watsat, sucsat, bsw, smpmin, nlevsoi)

Equilibrium water content and matric potential over all active (column, soil
layer) pairs (Zeng-Decker 2009). Backend-agnostic 2D kernel; replaces the inline
double loop in `soilwater_zengdecker2009!`.
"""
soilwm_voleq_zq!(vol_eq, zq, mask, zwtmm, zimm_arr, watsat, sucsat, bsw, smpmin,
        nlevsoi::Int) =
    _launch!(_soilwm_voleq_zq_kernel!, vol_eq, zq, mask, zwtmm, zimm_arr,
             watsat, sucsat, bsw, smpmin; ndrange = (length(mask), nlevsoi))

# --------------------------------------------------------------------------
# Hydraulic conductivity hk, matric potential smp and their derivatives
# (Zeng-Decker 2009 interface/node values). 2D over (column, soil layer);
# masked in-kernel. The interface value reads layer j and j+1 (clamped to
# nlevsoi) — read-only — so the per-(c,j) writes are disjoint. Writes hk,
# dhkdw, imped, smp, dsmpdw and mirrors smp/hk to the per-layer state outputs.
# --------------------------------------------------------------------------
@kernel function _soilwm_hk_smp_kernel!(hk, dhkdw, imped_arr, smp_arr, dsmpdw,
        smp_l, hk_l, @Const(mask), @Const(vwc_liq), @Const(watsat), @Const(hksat),
        @Const(bsw), @Const(icefrac), @Const(sucsat), @Const(smpmin),
        nlevsoi::Int, e_ice)
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        T   = eltype(hk)
        jp1 = min(nlevsoi, j + 1)
        b   = bsw[c, j]

        s1 = T(0.5) * (vwc_liq[c, j] + vwc_liq[c, jp1]) /
             (T(0.5) * (watsat[c, j] + watsat[c, jp1]))
        s1 = smooth_min(one(T), s1)
        s2 = hksat[c, j] * s1^(T(2.0) * b + T(2.0))

        imp = T(10.0)^(-e_ice * (T(0.5) * (icefrac[c, j] + icefrac[c, jp1])))
        imped_arr[c, j] = imp

        hk[c, j] = imp * s1 * s2
        dhkdw[c, j] = imp * (T(2.0) * b + T(3.0)) * s2 *
                      (one(T) / (watsat[c, j] + watsat[c, jp1]))

        s_node = smooth_max(vwc_liq[c, j] / watsat[c, j], T(0.01))
        s_node = smooth_min(one(T), s_node)

        smpv = -sucsat[c, j] * s_node^(-b)
        smpv = smooth_max(smpmin[c], smpv)
        smp_arr[c, j] = smpv

        dsmpdw[c, j] = -b * smpv / vwc_liq[c, j]

        smp_l[c, j] = smpv
        hk_l[c, j]  = hk[c, j]
    end
end

"""
    soilwm_hk_smp!(hk, dhkdw, imped_arr, smp_arr, dsmpdw, smp_l, hk_l, mask,
        vwc_liq, watsat, hksat, bsw, icefrac, sucsat, smpmin, nlevsoi, e_ice)

Hydraulic conductivity, matric potential and derivatives over all active
(column, soil layer) pairs (Zeng-Decker 2009). Backend-agnostic 2D kernel;
replaces the inline double loop in `soilwater_zengdecker2009!`.
"""
function soilwm_hk_smp!(hk, dhkdw, imped_arr, smp_arr, dsmpdw, smp_l, hk_l,
        mask, vwc_liq, watsat, hksat, bsw, icefrac, sucsat, smpmin,
        nlevsoi::Int, e_ice)
    ee = convert(eltype(hk), e_ice)
    _launch!(_soilwm_hk_smp_kernel!, hk, dhkdw, imped_arr, smp_arr, dsmpdw,
             smp_l, hk_l, mask, vwc_liq, watsat, hksat, bsw, icefrac, sucsat,
             smpmin, nlevsoi, ee; ndrange = (length(mask), nlevsoi))
end

# --------------------------------------------------------------------------
# Locate the layer index jwt right above the water table, and the volumetric
# water content vwc_zwt at the water-table depth (Zeng-Decker 2009). 1D per
# column; masked in-kernel. Each column runs two short sequential search loops
# (with break) — the per-column-kernel-with-internal-loops pattern — and writes
# only its own jwt[c] (Int) and vwc_zwt[c]. Physical constants are
# eltype-converted so no Float64 reaches a Float32-only backend.
# --------------------------------------------------------------------------
@kernel function _soilwm_jwt_vwczwt_kernel!(vwc_zwt, jwt, @Const(mask), @Const(zwt),
        @Const(zi), @Const(watsat), @Const(vwc_liq), @Const(t_soisno),
        @Const(sucsat), @Const(bsw), @Const(h2osoi_vol),
        joff::Int, joff_zi::Int, nlevsoi::Int, nlevgrnd::Int)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(vwc_zwt)
        # layer right above the water table
        jw = nlevsoi
        for j in 1:nlevsoi
            if zwt[c] <= zi[c, joff_zi + j]
                jw = j - 1
                break
            end
        end
        jwt[c] = jw

        # vwc at water-table depth
        vwc_zwt[c] = watsat[c, nlevsoi]
        if t_soisno[c, joff + jw + 1] < T(TFRZ)
            v = vwc_liq[c, nlevsoi]
            for j in nlevsoi:nlevgrnd
                if zwt[c] <= zi[c, joff_zi + j]
                    smp1 = T(HFUS) * (T(TFRZ) - t_soisno[c, joff + j]) /
                           (T(GRAV) * t_soisno[c, joff + j]) * T(1000.0)
                    smp1 = smooth_max(sucsat[c, nlevsoi], smp1)
                    v = watsat[c, nlevsoi] *
                        (smp1 / sucsat[c, nlevsoi])^(-one(T) / bsw[c, nlevsoi])
                    v = smooth_min(v, T(0.5) * (watsat[c, nlevsoi] + h2osoi_vol[c, nlevsoi]))
                    break
                end
            end
            vwc_zwt[c] = v
        end
    end
end

"""
    soilwm_jwt_vwczwt!(vwc_zwt, jwt, mask, zwt, zi, watsat, vwc_liq, t_soisno,
        sucsat, bsw, h2osoi_vol, joff, joff_zi, nlevsoi, nlevgrnd)

Water-table layer index jwt and water-table-depth water content vwc_zwt for each
active column. Backend-agnostic 1D kernel; replaces an inline per-column loop in
`soilwater_zengdecker2009!`.
"""
soilwm_jwt_vwczwt!(vwc_zwt, jwt, mask, zwt, zi, watsat, vwc_liq, t_soisno,
        sucsat, bsw, h2osoi_vol, joff::Int, joff_zi::Int, nlevsoi::Int,
        nlevgrnd::Int) =
    _launch!(_soilwm_jwt_vwczwt_kernel!, vwc_zwt, jwt, mask, zwt, zi, watsat,
             vwc_liq, t_soisno, sucsat, bsw, h2osoi_vol, joff, joff_zi,
             nlevsoi, nlevgrnd; ndrange = length(mask))

# --------------------------------------------------------------------------
# Equilibrium water content / matric potential for the aquifer (11th) layer
# when the water table is below the soil column (jwt == nlevsoi). 1D per column;
# masked in-kernel; writes only its own [c, nlevsoi+1] entries.
# --------------------------------------------------------------------------
@kernel function _soilwm_voleq_zq_aqu_kernel!(vol_eq, zq, @Const(mask), @Const(jwt),
        @Const(zwtmm), @Const(zimm_arr), @Const(watsat), @Const(sucsat),
        @Const(bsw), @Const(smpmin), nlevsoi::Int)
    c = @index(Global)
    @inbounds if mask[c] && jwt[c] == nlevsoi
        T = eltype(vol_eq)
        j = nlevsoi
        suc = sucsat[c, j]; ws = watsat[c, j]; b = bsw[c, j]
        zwt = zwtmm[c]
        zimm_j = zimm_arr[c, j+1]
        onem1b = one(T) - one(T) / b
        tempi  = one(T)
        temp0  = ((suc + zwt - zimm_j) / suc)^onem1b
        v = -suc * ws / onem1b / (zwt - zimm_j) * (tempi - temp0)
        v = smooth_max(v, zero(T))
        v = smooth_min(ws, v)
        vol_eq[c, j+1] = v
        zqv = -suc * (smooth_max(v / ws, T(0.01)))^(-b)
        zq[c, j+1] = smooth_max(smpmin[c], zqv)
    end
end

"""
    soilwm_voleq_zq_aqu!(vol_eq, zq, mask, jwt, zwtmm, zimm_arr, watsat, sucsat,
        bsw, smpmin, nlevsoi)

Aquifer-layer equilibrium water content/matric potential for columns whose water
table is below the soil column. Backend-agnostic 1D kernel.
"""
soilwm_voleq_zq_aqu!(vol_eq, zq, mask, jwt, zwtmm, zimm_arr, watsat, sucsat,
        bsw, smpmin, nlevsoi::Int) =
    _launch!(_soilwm_voleq_zq_aqu_kernel!, vol_eq, zq, mask, jwt, zwtmm,
             zimm_arr, watsat, sucsat, bsw, smpmin, nlevsoi;
             ndrange = length(mask))

# --------------------------------------------------------------------------
# Aquifer (11th) layer geometry: mid-depth zmm and thickness dzmm. 1D per
# column; masked in-kernel; writes only its own [c, nlevsoi+1] entries.
# --------------------------------------------------------------------------
@kernel function _soilwm_aqu_zmm_kernel!(zmm, dzmm, @Const(mask), @Const(jwt),
        @Const(zwt), nlevsoi::Int)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(zmm)
        zmm[c, nlevsoi+1] = T(0.5) * (T(1.0e3) * zwt[c] + zmm[c, nlevsoi])
        if jwt[c] < nlevsoi
            dzmm[c, nlevsoi+1] = dzmm[c, nlevsoi]
        else
            dzmm[c, nlevsoi+1] = T(1.0e3) * zwt[c] - zmm[c, nlevsoi]
        end
    end
end

"""
    soilwm_aqu_zmm!(zmm, dzmm, mask, jwt, zwt, nlevsoi)

Aquifer-layer mid-depth and thickness for each active column. Backend-agnostic
1D kernel.
"""
soilwm_aqu_zmm!(zmm, dzmm, mask, jwt, zwt, nlevsoi::Int) =
    _launch!(_soilwm_aqu_zmm_kernel!, zmm, dzmm, mask, jwt, zwt, nlevsoi;
             ndrange = length(mask))

# ==========================================================================
# Zeng-Decker 2009 tridiagonal-system assembly (top + interior + bottom/aquifer
# nodes). The loop touches ~25 distinct arrays — past the per-array kernel-arg
# limit — so they are grouped into two immutable device-view structs (the
# soil_temperature! grouped-struct template). The three original loops write
# DISJOINT rows (row 1; rows 2..nlevsoi-1; rows nlevsoi & nlevsoi+1) with no
# cross-row output reads, so they fuse into ONE per-column kernel that runs the
# interior as an internal j-loop — one launch, fully parallel over columns.
# ==========================================================================
Base.@kwdef struct SwmTriOut{M}   # tridiagonal outputs (all nc×(nlevsoi+1) matrices)
    qin::M; qout::M; dqidw0::M; dqidw1::M; dqodw1::M; dqodw2::M
    rmx::M; amx::M; bmx::M; cmx::M
end
Base.@kwdef struct SwmTriIn{M,V}  # read-only inputs (M = per-layer matrix, V = per-column vector)
    zmm::M; dzmm::M; zq::M; smp::M; hk::M; dhkdw::M; dsmpdw::M
    qflx_rootsoi::M; watsat::M; sucsat::M; bsw::M; vwc_liq::M
    qflx_infl::V; vwc_zwt::V; smpmin::V
end
Adapt.@adapt_structure SwmTriOut
Adapt.@adapt_structure SwmTriIn

@kernel function _soilwm_tridiag_assemble_kernel!(out, tin, @Const(mask), @Const(jwt),
        sdamp, dtime, nlevsoi::Int)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(out.amx)

        # ---- Node j=1 (top) ----
        j = 1
        out.qin[c, j] = tin.qflx_infl[c]
        den = tin.zmm[c, j+1] - tin.zmm[c, j]
        dzq = tin.zq[c, j+1] - tin.zq[c, j]
        num = (tin.smp[c, j+1] - tin.smp[c, j]) - dzq
        out.qout[c, j] = -tin.hk[c, j] * num / den
        out.dqodw1[c, j] = -(-tin.hk[c, j] * tin.dsmpdw[c, j] + num * tin.dhkdw[c, j]) / den
        out.dqodw2[c, j] = -(tin.hk[c, j] * tin.dsmpdw[c, j+1] + num * tin.dhkdw[c, j]) / den
        out.rmx[c, j] = out.qin[c, j] - out.qout[c, j] - tin.qflx_rootsoi[c, j]
        out.amx[c, j] = zero(T)
        # NOTE: top uses (sdamp + 1/dtime) [reciprocal], interior/bottom use /dtime
        # [division] — kept as distinct forms to stay byte-identical on Float64.
        out.bmx[c, j] = tin.dzmm[c, j] * (sdamp + one(T) / dtime) + out.dqodw1[c, j]
        out.cmx[c, j] = out.dqodw2[c, j]

        # ---- Interior nodes j=2..nlevsoi-1 ----
        for j in 2:(nlevsoi - 1)
            den = tin.zmm[c, j] - tin.zmm[c, j-1]
            dzq = tin.zq[c, j] - tin.zq[c, j-1]
            num = (tin.smp[c, j] - tin.smp[c, j-1]) - dzq
            out.qin[c, j] = -tin.hk[c, j-1] * num / den
            out.dqidw0[c, j] = -(-tin.hk[c, j-1] * tin.dsmpdw[c, j-1] + num * tin.dhkdw[c, j-1]) / den
            out.dqidw1[c, j] = -(tin.hk[c, j-1] * tin.dsmpdw[c, j] + num * tin.dhkdw[c, j-1]) / den
            den = tin.zmm[c, j+1] - tin.zmm[c, j]
            dzq = tin.zq[c, j+1] - tin.zq[c, j]
            num = (tin.smp[c, j+1] - tin.smp[c, j]) - dzq
            out.qout[c, j] = -tin.hk[c, j] * num / den
            out.dqodw1[c, j] = -(-tin.hk[c, j] * tin.dsmpdw[c, j] + num * tin.dhkdw[c, j]) / den
            out.dqodw2[c, j] = -(tin.hk[c, j] * tin.dsmpdw[c, j+1] + num * tin.dhkdw[c, j]) / den
            out.rmx[c, j] = out.qin[c, j] - out.qout[c, j] - tin.qflx_rootsoi[c, j]
            out.amx[c, j] = -out.dqidw0[c, j]
            out.bmx[c, j] = tin.dzmm[c, j] / dtime - out.dqidw1[c, j] + out.dqodw1[c, j]
            out.cmx[c, j] = out.dqodw2[c, j]
        end

        # ---- Node j=nlevsoi (bottom) ----
        j = nlevsoi
        if j > jwt[c]   # water table is in soil column
            den = tin.zmm[c, j] - tin.zmm[c, j-1]
            dzq = tin.zq[c, j] - tin.zq[c, j-1]
            num = (tin.smp[c, j] - tin.smp[c, j-1]) - dzq
            out.qin[c, j] = -tin.hk[c, j-1] * num / den
            out.dqidw0[c, j] = -(-tin.hk[c, j-1] * tin.dsmpdw[c, j-1] + num * tin.dhkdw[c, j-1]) / den
            out.dqidw1[c, j] = -(tin.hk[c, j-1] * tin.dsmpdw[c, j] + num * tin.dhkdw[c, j-1]) / den
            out.qout[c, j] = zero(T)
            out.dqodw1[c, j] = zero(T)
            out.rmx[c, j] = out.qin[c, j] - out.qout[c, j] - tin.qflx_rootsoi[c, j]
            out.amx[c, j] = -out.dqidw0[c, j]
            out.bmx[c, j] = tin.dzmm[c, j] / dtime - out.dqidw1[c, j] + out.dqodw1[c, j]
            out.cmx[c, j] = zero(T)

            # Aquifer layer — hydrologically inactive
            out.rmx[c, j+1] = zero(T)
            out.amx[c, j+1] = zero(T)
            out.bmx[c, j+1] = tin.dzmm[c, j+1] / dtime
            out.cmx[c, j+1] = zero(T)
        else            # water table is below soil column
            s_node = smooth_max(T(0.5) * ((tin.vwc_zwt[c] + tin.vwc_liq[c, j]) / tin.watsat[c, j]), T(0.01))
            s_node = smooth_min(one(T), s_node)
            smp1 = -tin.sucsat[c, j] * s_node^(-tin.bsw[c, j])
            smp1 = smooth_max(tin.smpmin[c], smp1)
            dsmpdw1 = -tin.bsw[c, j] * smp1 / (s_node * tin.watsat[c, j])

            den = tin.zmm[c, j] - tin.zmm[c, j-1]
            dzq = tin.zq[c, j] - tin.zq[c, j-1]
            num = (tin.smp[c, j] - tin.smp[c, j-1]) - dzq
            out.qin[c, j] = -tin.hk[c, j-1] * num / den
            out.dqidw0[c, j] = -(-tin.hk[c, j-1] * tin.dsmpdw[c, j-1] + num * tin.dhkdw[c, j-1]) / den
            out.dqidw1[c, j] = -(tin.hk[c, j-1] * tin.dsmpdw[c, j] + num * tin.dhkdw[c, j-1]) / den
            den = tin.zmm[c, j+1] - tin.zmm[c, j]
            dzq = tin.zq[c, j+1] - tin.zq[c, j]
            num = (smp1 - tin.smp[c, j]) - dzq
            out.qout[c, j] = -tin.hk[c, j] * num / den
            out.dqodw1[c, j] = -(-tin.hk[c, j] * tin.dsmpdw[c, j] + num * tin.dhkdw[c, j]) / den
            out.dqodw2[c, j] = -(tin.hk[c, j] * dsmpdw1 + num * tin.dhkdw[c, j]) / den
            out.rmx[c, j] = out.qin[c, j] - out.qout[c, j] - tin.qflx_rootsoi[c, j]
            out.amx[c, j] = -out.dqidw0[c, j]
            out.bmx[c, j] = tin.dzmm[c, j] / dtime - out.dqidw1[c, j] + out.dqodw1[c, j]
            out.cmx[c, j] = out.dqodw2[c, j]

            # Aquifer layer
            out.qin[c, j+1] = out.qout[c, j]
            out.dqidw0[c, j+1] = -(-tin.hk[c, j] * tin.dsmpdw[c, j] + num * tin.dhkdw[c, j]) / den
            out.dqidw1[c, j+1] = -(tin.hk[c, j] * dsmpdw1 + num * tin.dhkdw[c, j]) / den
            out.qout[c, j+1] = zero(T)
            out.dqodw1[c, j+1] = zero(T)
            out.rmx[c, j+1] = out.qin[c, j+1] - out.qout[c, j+1]
            out.amx[c, j+1] = -out.dqidw0[c, j+1]
            out.bmx[c, j+1] = tin.dzmm[c, j+1] / dtime - out.dqidw1[c, j+1] + out.dqodw1[c, j+1]
            out.cmx[c, j+1] = zero(T)
        end
    end
end

"""
    soilwm_tridiag_assemble!(out::SwmTriOut, tin::SwmTriIn, mask, jwt, sdamp, dtime, nlevsoi)

Assemble the Zeng-Decker 2009 tridiagonal system (a/b/c/r + flux/derivative
arrays) for all active columns. One per-column kernel; scalars are
eltype-converted so no Float64 reaches a Float32-only backend.
"""
function soilwm_tridiag_assemble!(out::SwmTriOut, tin::SwmTriIn, mask, jwt,
        sdamp, dtime, nlevsoi::Int)
    nc = length(mask)
    nc == 0 && return nothing
    T = eltype(out.amx)
    backend = _kernel_backend(out.amx)
    _soilwm_tridiag_assemble_kernel!(backend)(out, tin, mask, jwt,
        convert(T, sdamp), convert(T, dtime), nlevsoi; ndrange = nc)
    KA.synchronize(backend)
    return nothing
end

# --------------------------------------------------------------------------
# Post-solve: renew the layer liquid water mass from the tridiagonal solution
# dwat2 and compute aquifer recharge qcharge. 1D per column; masked in-kernel.
# Each column updates its own h2osoi_liq[joff+1..joff+nlevsoi] (internal j-loop)
# and writes qcharge[c], branching on whether the water table is in the soil
# column (jwt < nlevsoi) or below it. ~15 array args (under the grouped-struct
# threshold). Literals eltype-converted; /dtime kept as division (byte-identical).
# --------------------------------------------------------------------------
@kernel function _soilwm_renew_qcharge_kernel!(h2osoi_liq, qcharge, @Const(mask),
        @Const(jwt), @Const(dwat2), @Const(dzmm), @Const(h2osoi_vol), @Const(watsat),
        @Const(imped_arr), @Const(hksat), @Const(bsw), @Const(smpmin), @Const(smp_arr),
        @Const(zq), @Const(zwt), @Const(z), joff::Int, nlevsoi::Int, dtime)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(h2osoi_liq)
        for j in 1:nlevsoi
            h2osoi_liq[c, joff + j] = h2osoi_liq[c, joff + j] + dwat2[c, j] * dzmm[c, j]
        end

        jw = jwt[c]
        if jw < nlevsoi
            wh_zwt = zero(T)
            s_node = smooth_max(h2osoi_vol[c, jw+1] / watsat[c, jw+1], T(0.01))
            s1 = smooth_min(one(T), s_node)
            ka = imped_arr[c, jw+1] * hksat[c, jw+1] * s1^(T(2.0) * bsw[c, jw+1] + T(3.0))
            smp1 = smooth_max(smpmin[c], smp_arr[c, max(1, jw)])
            wh = smp1 - zq[c, max(1, jw)]
            if jw == 0
                qcharge[c] = -ka * (wh_zwt - wh) / ((zwt[c] + T(1.0e-3)) * T(1000.0))
            else
                qcharge[c] = -ka * (wh_zwt - wh) / ((zwt[c] - z[c, joff + jw]) * T(1000.0) * T(2.0))
            end
            qcharge[c] = smooth_max(-T(10.0) / dtime, qcharge[c])
            qcharge[c] = smooth_min(T(10.0) / dtime, qcharge[c])
        else
            qcharge[c] = dwat2[c, nlevsoi+1] * dzmm[c, nlevsoi+1] / dtime
        end
    end
end

"""
    soilwm_renew_qcharge!(h2osoi_liq, qcharge, mask, jwt, dwat2, dzmm, h2osoi_vol,
        watsat, imped_arr, hksat, bsw, smpmin, smp_arr, zq, zwt, z, joff, nlevsoi, dtime)

Update layer liquid water mass and aquifer recharge after the tridiagonal solve,
for each active column. Backend-agnostic 1D kernel.
"""
function soilwm_renew_qcharge!(h2osoi_liq, qcharge, mask, jwt, dwat2, dzmm,
        h2osoi_vol, watsat, imped_arr, hksat, bsw, smpmin, smp_arr, zq, zwt, z,
        joff::Int, nlevsoi::Int, dtime)
    dt = convert(eltype(h2osoi_liq), dtime)
    _launch!(_soilwm_renew_qcharge_kernel!, h2osoi_liq, qcharge, mask, jwt, dwat2,
        dzmm, h2osoi_vol, watsat, imped_arr, hksat, bsw, smpmin, smp_arr, zq, zwt,
        z, joff, nlevsoi, dt; ndrange = length(mask))
end

# --------------------------------------------------------------------------
# Water deficit: per-column sum of any negative layer liquid water (as a
# positive deficit) via smooth_ifelse. 1D per column; internal j-reduction in
# ascending order (matches the scalar accumulation).
# --------------------------------------------------------------------------
@kernel function _soilwm_deficit_kernel!(qflx_deficit, @Const(mask),
        @Const(h2osoi_liq), joff::Int, nlevsoi::Int)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(qflx_deficit)
        acc = zero(T)
        for j in 1:nlevsoi
            acc = acc + smooth_ifelse(h2osoi_liq[c, joff + j], zero(T),
                                      -h2osoi_liq[c, joff + j])
        end
        qflx_deficit[c] = acc
    end
end

"""
    soilwm_deficit!(qflx_deficit, mask, h2osoi_liq, joff, nlevsoi)

Per-column water deficit (sum of negative layer liquid water as positive)
for each active column. Backend-agnostic 1D kernel.
"""
soilwm_deficit!(qflx_deficit, mask, h2osoi_liq, joff::Int, nlevsoi::Int) =
    _launch!(_soilwm_deficit_kernel!, qflx_deficit, mask, h2osoi_liq, joff, nlevsoi;
             ndrange = length(mask))

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
function ice_impedance(icefrac::Real, e_ice::Real)
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
function tridiagonal_col!(u::AbstractVector{<:Real},
                          a::AbstractVector{<:Real},
                          b::AbstractVector{<:Real},
                          c::AbstractVector{<:Real},
                          r::AbstractVector{<:Real},
                          jtop::Int, lbj::Int, ubj::Int)
    T = promote_type(eltype(a), eltype(b), eltype(c), eltype(r))
    n = ubj - lbj + 1
    gam = zeros(T, n)

    bet = b[jtop]
    if abs(bet) < 1.0e-30
        bet = 1.0e-30
    end

    for j in lbj:ubj
        if j >= jtop
            if j == jtop
                u[j] = r[j] / bet
            else
                gam[j] = c[j-1] / bet
                bet = b[j] - a[j] * gam[j]
                if abs(bet) < 1.0e-30
                    bet = copysign(1.0e-30, bet == 0.0 ? one(bet) : bet)
                end
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
function baseflow_sink!(baseflow_sink::AbstractMatrix{<:Real},
                        mask_hydrology::AbstractVector{Bool},
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
    FT = eltype(vwc_liq)
    s2 = zeros(FT, nlayers)
    for j in 1:nlayers
        s2[j] = vwc_liq[j] / watsat[c, j]
        s2[j] = smooth_min(s2[j], 1.0)
        s2[j] = smooth_max(0.01, s2[j])
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
        s1 = smooth_min(s1, 1.0)
        s1 = smooth_max(0.01, s1)

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
        z_col::Matrix{<:Real}, zi_col::Matrix{<:Real}, dz_col::Matrix{<:Real},
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
    nlevsno    = varpar.nlevsno
    joff       = nlevsno          # offset for z, dz arrays (snow+soil combined)
    joff_zi    = nlevsno + 1     # offset for zi array
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
    den = M_TO_MM * (z_col[c, joff + j+1] - z_col[c, joff + j])
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
        den = M_TO_MM * (z_col[c, joff + j+1] - z_col[c, joff + j])
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
            if zwt[c] <= zi_col[c, joff_zi + jj]
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
            den = M_TO_MM * (zwt[c] - z_col[c, joff + j])
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
        mask_hydrology::AbstractVector{Bool},
        soilhydrology::SoilHydrologyData, soilstate::SoilStateData,
        waterstatebulk::WaterStateBulkData,
        swrc::SoilWaterRetentionCurve, cfg::SoilWaterMovementConfig,
        dwat::Matrix{<:Real}, smp::Matrix{<:Real},
        imped::Matrix{<:Real}, vwc_liq::Matrix{<:Real},
        dtime::Real, bounds::UnitRange{Int})

    nlevsoi = varpar.nlevsoi
    nlevsno = varpar.nlevsno
    joff    = nlevsno          # offset for combined snow+soil arrays (z)
    joff_zi = nlevsno + 1     # offset for zi
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
            if zwt[c] <= zi[c, joff_zi + j]
                jwt = j - 1
                break
            end
        end

        # Calculate qcharge for case jwt < nlevsoi
        if jwt < nlevsoi
            wh_zwt = -sucsat[c, jwt+1] - zwt[c] * M_TO_MM

            # Recharge rate to groundwater (positive to aquifer)
            s1 = smooth_max(vwc_liq[c, jwt+1] / watsat[c, jwt+1], 0.01)
            s1 = smooth_min(1.0, s1)

            # Unsaturated hydraulic conductivity
            ka, _ = soil_hk!(swrc, c, jwt+1, s1, imped[c, jwt+1], soilstate)

            smp1 = smooth_max(smpmin[c], smp[c, max(1, jwt)])
            wh = smp1 - z[c, joff + max(1, jwt)] * M_TO_MM

            if jwt == 0
                qcharge[c] = -ka * (wh_zwt - wh) / ((zwt[c] + 1.0e-3) * M_TO_MM)
            else
                qcharge[c] = -ka * (wh_zwt - wh) / ((zwt[c] - z[c, joff + jwt]) * M_TO_MM * 2.0)
            end

            # Limit qcharge
            qcharge[c] = smooth_max(-10.0 / dtime, qcharge[c])
            qcharge[c] = smooth_min(10.0 / dtime, qcharge[c])
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
        mask_hydrology::AbstractVector{Bool}, mask_urban::AbstractVector{Bool},
        soilhydrology::SoilHydrologyData, soilstate::SoilStateData,
        waterfluxbulk::WaterFluxBulkData, waterstatebulk::WaterStateBulkData,
        temperature::TemperatureData, canopystate::CanopyStateData,
        energyflux::EnergyFluxData,
        swrc::SoilWaterRetentionCurve, cfg::SoilWaterMovementConfig,
        dtime::Real)

    nlevsoi = varpar.nlevsoi
    nlevsno = varpar.nlevsno
    joff    = nlevsno          # offset for combined snow+soil arrays (z, dz, h2osoi_liq/ice)
    nbedrock = col_data.nbedrock
    z  = col_data.z
    zi = col_data.zi
    dz = col_data.dz

    nsubsteps_col = soilhydrology.num_substeps_col
    qcharge       = soilhydrology.qcharge_col
    h2osoi_liq    = waterstatebulk.ws.h2osoi_liq_col
    qflx_rootsoi  = waterfluxbulk.qflx_rootsoi_col
    FT = eltype(h2osoi_liq)
    hk_loc     = zeros(FT, nlevsoi)
    smp_loc    = zeros(FT, nlevsoi)
    dhkdw_loc  = zeros(FT, nlevsoi)
    dsmpdw_loc = zeros(FT, nlevsoi)
    imped_loc  = zeros(FT, nlevsoi)
    vwc_liq_loc = zeros(FT, nlevsoi)
    dt_dz_loc  = zeros(FT, nlevsoi)
    qin_loc    = zeros(FT, nlevsoi)
    qout_loc   = zeros(FT, nlevsoi)
    dqidw0_loc = zeros(FT, nlevsoi)
    dqidw1_loc = zeros(FT, nlevsoi)
    dqodw1_loc = zeros(FT, nlevsoi)
    dqodw2_loc = zeros(FT, nlevsoi)
    dwat_loc   = zeros(FT, nlevsoi)
    amx_loc    = zeros(FT, nlevsoi)
    bmx_loc    = zeros(FT, nlevsoi)
    cmx_loc    = zeros(FT, nlevsoi)
    rmx_loc    = zeros(FT, nlevsoi)
    fluxNet0   = zeros(FT, nlevsoi)
    fluxNet1   = zeros(FT, nlevsoi)

    for c in eachindex(mask_hydrology)
        mask_hydrology[c] || continue

        nlayers = nbedrock[c]

        # Zero workspace for this column's active layers
        for arr in (hk_loc, smp_loc, dhkdw_loc, dsmpdw_loc, imped_loc,
                    vwc_liq_loc, dt_dz_loc, qin_loc, qout_loc,
                    dqidw0_loc, dqidw1_loc, dqodw1_loc, dqodw2_loc,
                    dwat_loc, amx_loc, bmx_loc, cmx_loc, rmx_loc,
                    fluxNet0, fluxNet1)
            @inbounds for j in 1:nlayers
                arr[j] = 0.0
            end
        end

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
                vwc_liq_loc[j] = smooth_max(h2osoi_liq[c, joff + j], 1.0e-6) / (dz[c, joff + j] * DENH2O)
                dt_dz_loc[j] = dtsub / (M_TO_MM * dz[c, joff + j])
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
            @inbounds for j in 1:nlayers
                fluxNet0[j] = 0.0
                fluxNet1[j] = 0.0
            end

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
                h2osoi_liq[c, joff + j] = h2osoi_liq[c, joff + j] + dwat_loc[j] * (M_TO_MM * dz[c, joff + j])
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
        FT = eltype(h2osoi_liq)
        dwat_2d  = zeros(FT, nc, nlevsoi)
        smp_2d   = zeros(FT, nc, nlevsoi)
        imped_2d = zeros(FT, nc, nlevsoi)
        vwc_2d   = zeros(FT, nc, nlevsoi)

        for c in eachindex(mask_hydrology)
            mask_hydrology[c] || continue
            nlayers = nbedrock[c]
            for j in 1:min(nlayers, nlevsoi)
                vwc_2d[c, j] = smooth_max(h2osoi_liq[c, joff + j], 1.0e-6) / (dz[c, joff + j] * DENH2O)
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
        mask_hydrology::AbstractVector{Bool}, mask_urban::AbstractVector{Bool},
        soilhydrology::SoilHydrologyData, soilstate::SoilStateData,
        waterfluxbulk::WaterFluxBulkData, waterstatebulk::WaterStateBulkData,
        temperature::TemperatureData, canopystate::CanopyStateData,
        energyflux::EnergyFluxData,
        swrc::SoilWaterRetentionCurve, cfg::SoilWaterMovementConfig,
        dtime::Real)

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

    # Local scratch — allocated on the backend of the state via similar()+fill!
    # (NOT zeros(), which is always host memory) so the kernels below run
    # device-resident on GPU. zerocol2(m) → nc×m float matrix of zeros.
    FT = eltype(h2osoi_liq)
    zerocol2(m::Int) = fill!(similar(h2osoi_liq, FT, nc, m), zero(FT))
    zerocol1()       = fill!(similar(h2osoi_liq, FT, nc), zero(FT))
    hk      = zerocol2(nlevsoi)
    dhkdw   = zerocol2(nlevsoi)
    smp_arr = zerocol2(nlevsoi)
    dsmpdw  = zerocol2(nlevsoi)
    imped_arr = zerocol2(nlevsoi)
    vol_ice = zerocol2(nlevsoi)
    vwc_liq = zerocol2(nlevsoi + 1)
    zmm     = zerocol2(nlevsoi + 1)
    dzmm    = zerocol2(nlevsoi + 1)
    # 1-indexed with offset: zimm_arr[c, j+1] for Fortran zimm(c,j), j=0:nlevsoi
    # (index 1 = Fortran j=0, index nlevsoi+1 = Fortran j=nlevsoi)
    zimm_arr = zerocol2(nlevsoi + 1)
    zwtmm   = zerocol1()
    jwt     = fill!(similar(h2osoi_liq, Int, nc), 0)
    vwc_zwt = zerocol1()
    vol_eq  = zerocol2(nlevsoi + 1)
    zq      = zerocol2(nlevsoi + 1)
    qin     = zerocol2(nlevsoi + 1)
    qout_arr = zerocol2(nlevsoi + 1)
    dqidw0  = zerocol2(nlevsoi + 1)
    dqidw1  = zerocol2(nlevsoi + 1)
    dqodw1  = zerocol2(nlevsoi + 1)
    dqodw2  = zerocol2(nlevsoi + 1)
    amx     = zerocol2(nlevsoi + 1)
    bmx     = zerocol2(nlevsoi + 1)
    cmx     = zerocol2(nlevsoi + 1)
    rmx     = zerocol2(nlevsoi + 1)
    dwat2   = zerocol2(nlevsoi + 1)

    sdamp = 0.0

    # Convert to mm and compute ice fractions (per-(column,layer) kernel)
    soilwm_mm_icefrac!(zmm, dzmm, zimm_arr, vol_ice, icefrac, vwc_liq,
        mask_hydrology, z, dz, zi, watsat, h2osoi_ice, h2osoi_liq,
        joff, joff_zi, nlevsoi; denice = DENICE, denh2o = DENH2O)

    soilwm_zwtmm!(zwtmm, zimm_arr, mask_hydrology, zwt)

    # Compute jwt index (layer above water table) and vwc at water-table depth
    # (per-column kernel; each column runs the two short search loops internally).
    soilwm_jwt_vwczwt!(vwc_zwt, jwt, mask_hydrology, zwt, zi, watsat, vwc_liq,
        t_soisno, sucsat, bsw, h2osoi_vol, joff, joff_zi, nlevsoi, nlevgrnd)

    # Calculate equilibrium water content based on water table depth
    # (per-(column, soil layer) kernel; reads its own layer + interface depths).
    soilwm_voleq_zq!(vol_eq, zq, mask_hydrology, zwtmm, zimm_arr,
        watsat, sucsat, bsw, smpmin, nlevsoi)

    # If water table is below soil column, calculate vol_eq/zq for the 11th
    # (aquifer) layer (per-column kernel; only columns with jwt==nlevsoi write).
    soilwm_voleq_zq_aqu!(vol_eq, zq, mask_hydrology, jwt, zwtmm, zimm_arr,
        watsat, sucsat, bsw, smpmin, nlevsoi)

    # Hydraulic conductivity and soil matric potential
    # (per-(column, soil layer) kernel; interface value reads layer j and j+1).
    soilwm_hk_smp!(hk, dhkdw, imped_arr, smp_arr, dsmpdw, smp_l, hk_l,
        mask_hydrology, vwc_liq, watsat, hksat, bsw, icefrac, sucsat, smpmin,
        nlevsoi, cfg.e_ice)

    # Aquifer (11th) layer geometry (per-column kernel).
    soilwm_aqu_zmm!(zmm, dzmm, mask_hydrology, jwt, zwt, nlevsoi)

    # ---- Set up tridiagonal system (single per-column kernel: top node + an
    # internal interior j-loop + bottom/aquifer node; writes disjoint rows). The
    # ~25 arrays are grouped into two device-view structs (grouped-struct template).
    tri_out = SwmTriOut(; qin = qin, qout = qout_arr, dqidw0 = dqidw0, dqidw1 = dqidw1,
        dqodw1 = dqodw1, dqodw2 = dqodw2, rmx = rmx, amx = amx, bmx = bmx, cmx = cmx)
    tri_in = SwmTriIn(; zmm = zmm, dzmm = dzmm, zq = zq, smp = smp_arr, hk = hk,
        dhkdw = dhkdw, dsmpdw = dsmpdw, qflx_rootsoi = qflx_rootsoi, watsat = watsat,
        sucsat = sucsat, bsw = bsw, vwc_liq = vwc_liq, qflx_infl = qflx_infl,
        vwc_zwt = vwc_zwt, smpmin = smpmin)
    soilwm_tridiag_assemble!(tri_out, tri_in, mask_hydrology, jwt, sdamp, dtime, nlevsoi)

    # Solve for dwat using multi-column tridiagonal solver
    jtop_arr = fill!(similar(h2osoi_liq, Int, nc), 1)
    mask_jtop = mask_hydrology
    tridiagonal_multi!(dwat2, amx, bmx, cmx, rmx,
        jtop_arr, mask_jtop, nc, nlevsoi + 1)

    # Renew the layer liquid water mass and compute qcharge (per-column kernel;
    # internal j-loop for the mass update, then the jwt-branched recharge).
    soilwm_renew_qcharge!(h2osoi_liq, qcharge, mask_hydrology, jwt, dwat2, dzmm,
        h2osoi_vol, watsat, imped_arr, hksat, bsw, smpmin, smp_arr, zq, zwt, z,
        joff, nlevsoi, dtime)

    # Compute water deficit (per-column reduction over the soil layers).
    soilwm_deficit!(qflx_deficit, mask_hydrology, h2osoi_liq, joff, nlevsoi)

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
        mask_hydrology::AbstractVector{Bool}, mask_urban::AbstractVector{Bool},
        soilhydrology::SoilHydrologyData, soilstate::SoilStateData,
        waterfluxbulk::WaterFluxBulkData, waterstatebulk::WaterStateBulkData,
        temperature::TemperatureData, canopystate::CanopyStateData,
        energyflux::EnergyFluxData,
        swrc::SoilWaterRetentionCurve, cfg::SoilWaterMovementConfig,
        dtime::Real;
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
                xs = smooth_ifelse(h2osoi_liq[c, joff_ + j], zero(eltype(h2osoi_liq)), watmin - h2osoi_liq[c, joff_ + j])
                h2osoi_liq[c, joff_ + j]   = h2osoi_liq[c, joff_ + j] + xs
                h2osoi_liq[c, joff_ + j+1] = h2osoi_liq[c, joff_ + j+1] - xs
            end
        end

        j = nlevsoi
        for c in eachindex(mask_hydrology)
            mask_hydrology[c] || continue
            xs = smooth_ifelse(h2osoi_liq[c, joff_ + j] - watmin, zero(eltype(h2osoi_liq)), watmin - h2osoi_liq[c, joff_ + j])
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
