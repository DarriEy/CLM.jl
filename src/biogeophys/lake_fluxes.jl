# ==========================================================================
# Ported from: src/biogeophys/LakeFluxesMod.F90 (781 lines)
# Lake surface flux calculations
#
# Calculates surface fluxes and surface temperature for lakes using
# Monin-Obukhov similarity theory with stability iteration.
#
# Handles: open water, frozen lakes (with/without snow), variable
# depth, fetch-dependent roughness, Newton-Raphson temperature solving.
#
# Public functions:
#   lake_fluxes!  — Main lake flux calculation
#
# Lake constants (from LakeCon):
#   emg_lake, betavis, minz0lake, tdmax
#
# GPU kernelization (Phase A, A3): the per-patch driver loop is moved into a
# single per-patch KernelAbstractions kernel (_lake_fluxes_kernel!), one thread
# per filtered patch. Each patch gathers its column/gridcell/landunit indices and
# runs the whole stability iteration sequentially in-thread (the iteration is
# loop-carried: obu/ustar/z0* of one pass feed the next), then writes only its own
# [p]/[c] outputs — no reductions, no atomics. The many state arrays are bundled
# into ONE @adapt_structure device-view struct (LakeFluxDV) so KA adapts MtlArray
# fields to device arrays at launch and the arg count stays under Metal's ~31 cap;
# struct fields are aliased back to locals inside the kernel so the body stays
# verbatim. All Float64 literals are eltype-converted (T(...)) so no Float64
# reaches a Float32-only backend (Metal); on Float64 this is byte-identical.
# ==========================================================================

# Lake-specific constants (normally from LakeCon.F90)
const EMG_LAKE = 0.97          # lake surface emissivity
const BETAVIS_LAKE = 0.0       # LakeCon.F90:45 betavis (was 0.4, which CTSM never uses)
# LakeCon.F90:143-158 sets fcrit and minz0lake as a PAIR, chosen by the
# `lake_use_old_fcrit_minz0` namelist flag (default .false., LakeCon.F90:98):
#   .true.  (Subin et al. 2011)  -> fcrit = 22.0,  minz0lake = 1.0e-5
#   .false. (Vickers&Mahrt 1997) -> fcrit = 100.0, minz0lake = 1.0e-10
# This file previously took FCRIT from the new branch and MINZ0LAKE from the old
# one, which is not a configuration CTSM can ever be in. `lake_con.jl` already
# defaults `lake_use_old_fcrit_minz0 = false`, so the new branch is the one the
# port is on: minz0lake is 1e-10, "set low so it is only to avoid floating point
# exceptions" — NOT a physical floor. At 1e-5 it is a 1e5x floor waiting to bite
# any case where the Charnock roughness legitimately drops below 1e-5 m.
const MINZ0LAKE = 1.0e-10      # minimum lake roughness length [m] (LakeCon minz0lake)
const TDMAX_LAKE = 277.0       # temperature of maximum water density [K]

# Aerodynamic parameters
const PRN_AIR = 0.713          # Prandtl number for air
const SCH_WATER = 0.66         # Schmidt number for water in air
const KVA0 = 1.51e-5           # kinematic viscosity of air at 20°C [m²/s]
const CUS = 0.1                # smooth flow regime coefficient
const CUR0 = 0.01              # min Charnock parameter (LakeCon cur0)
const CURM = 0.1               # max Charnock parameter (LakeCon curm; was wrongly 0 → no
                               # fetch-limited enhancement, leaving z0mg far too small)
const FCRIT = 100.0            # critical dimensionless fetch for the Charnock parameter
const BETA1 = 1.0              # coefficient of convective velocity
const ZETAMAX_LAKE = 0.5       # max zeta under stable conditions
const KVA0TEMP = 293.15        # 20 C + tfrz: reference temperature for KVA0 (LakeFluxesMod kva0temp)
const KVA0PRES = 1.013e5       # pressure for KVA0 [Pa] (LakeFluxesMod.F90:192)
const NU_PARAM_LAKE = 1.5e-5   # clm_varcon nu_param: viscosity in the ZengWang2007 z0hg
const ZII_LAKE = 1000.0        # convective boundary height [m] (LakeFluxesMod.F90:136)

# --------------------------------------------------------------------------
# Device-view struct: every array touched by the per-patch loop, bundled into
# ONE struct so the kernel takes a single struct arg (keeps the arg count under
# Metal's ~31 cap) and `Adapt.@adapt_structure` moves all MtlArray fields to the
# device at launch. Array fields are shared refs from the real state structs, so
# writes inside the kernel flow straight back to the caller's state. Type params:
# Vi=Int vec, V=float vec, M=float matrix (so a Float32 device build and a Float64
# host build both specialize cleanly).
# --------------------------------------------------------------------------
Base.@kwdef struct LakeFluxDV{Vi,V,M}
    # patch -> column/gridcell/landunit maps
    p_column::Vi; p_gridcell::Vi; p_landunit::Vi
    # column-level columns map (for ks_col latitude — simplified, uses col gridcell)
    c_gridcell::Vi
    # column scalars
    snl::Vi; lakedepth::V; savedtke1::V; ws_col::V; ks_col::V; htvp_col::V; ust_lake::V
    # column matrices
    dz::M; dz_lake::M; t_soisno::M; t_lake::M
    # column vectors
    t_grnd::V
    # patch solar (+ NIR split and SNICAR layer absorption for betaprime)
    sabg::V
    fsds_nir_d::V; fsds_nir_i::V; fsr_nir_d::V; fsr_nir_i::V
    sabg_lyr::M
    # column roughness lengths carried between steps (LakeFluxesMod.F90:327,764-766)
    z0mg_col::V; z0hg_col::V; z0qg_col::V
    # water state (column, layer)
    h2osoi_liq::M; h2osoi_ice::M
    # forcings (column-indexed)
    forc_t::V; forc_th::V; forc_q::V; forc_pbot::V; forc_rho::V; forc_lwrad::V
    # forcings (gridcell-indexed)
    forc_u::V; forc_v::V; forc_hgt_u::V; forc_hgt_t::V; forc_hgt_q::V
    # gridcell latitude [radians] — LakeFluxesMod.F90:755 needs grc%lat(g) for ks
    grc_lat::V
    # outputs (patch)
    eflx_sh_tot::V; eflx_sh_grnd::V; eflx_lh_tot::V; eflx_lh_grnd::V
    eflx_lwrad_out::V; eflx_lwrad_net::V; eflx_soil_grnd::V; eflx_gnet::V; taux::V; tauy::V
    qflx_evap_soi::V; qflx_evap_tot::V; ustar::V; z0mg::V; z0hg::V; z0qg::V
    # 2-m diagnostics (LakeFluxesMod.F90:708-717). These were NEVER written on a
    # lake patch: `t_ref2m_patch` is touched only by bareground_fluxes.jl and
    # urban_fluxes.jl, neither of which runs on a lake, so TSA on a 100%-lake
    # gridcell sat at its cold-start 283.000 K for the whole run while Fortran's
    # TSA ranged over 255-270 K. A dead write, not a physics error.
    t_ref2m::V; q_ref2m::V; rh_ref2m::V
    # Converged Monin-Obukhov length + stability parameter. Fortran carries both
    # as patch state; the port computed and discarded them, so nothing outside
    # the kernel could tell WHICH stability regime a lake step was running in.
    obu_out::V; zeta_out::V
end
Adapt.@adapt_structure LakeFluxDV

# Coerce the gridcell latitude to the SAME array type as the other gridcell-indexed
# forcings. `LakeFluxDV{Vi,V,M}` requires every float-vector field to share one type
# V, so handing it a host `Vector{Float64}` while the state is on the device (or a
# Float64 lat against a Float32 device build) would fail to construct. `similar` +
# `copyto!` keeps the device path working with one small gridcell-length allocation,
# and the common case (same type already) is a no-op passthrough.
function _lake_lat_like(proto, grc_lat)
    grc_lat === nothing && return zero(proto)
    typeof(grc_lat) === typeof(proto) && return grc_lat
    return copyto!(similar(proto), grc_lat)
end

# Build a LakeFluxDV from the caller's state structs (array fields are shared refs).
function _lake_flux_dv(temperature, energyflux, frictionvel, solarabs, lakestate,
                       waterstatebulk, waterdiagbulk, waterfluxbulk, col_data, patch_data,
                       forc_t, forc_th, forc_q, forc_pbot, forc_rho, forc_lwrad,
                       forc_u, forc_v, forc_hgt_u_grc, forc_hgt_t_grc, forc_hgt_q_grc,
                       grc_lat)
    return LakeFluxDV(;
        p_column = patch_data.column, p_gridcell = patch_data.gridcell,
        p_landunit = patch_data.landunit, c_gridcell = col_data.gridcell,
        snl = col_data.snl, lakedepth = col_data.lakedepth,
        savedtke1 = lakestate.savedtke1_col, ws_col = lakestate.ws_col,
        ks_col = lakestate.ks_col, htvp_col = energyflux.htvp_col,
        ust_lake = lakestate.ust_lake_col,
        dz = col_data.dz, dz_lake = col_data.dz_lake,
        t_soisno = temperature.t_soisno_col, t_lake = temperature.t_lake_col,
        t_grnd = temperature.t_grnd_col, sabg = solarabs.sabg_patch,
        fsds_nir_d = solarabs.fsds_nir_d_patch, fsds_nir_i = solarabs.fsds_nir_i_patch,
        fsr_nir_d = solarabs.fsr_nir_d_patch, fsr_nir_i = solarabs.fsr_nir_i_patch,
        sabg_lyr = solarabs.sabg_lyr_patch,
        z0mg_col = frictionvel.z0mg_col, z0hg_col = frictionvel.z0hg_col,
        z0qg_col = frictionvel.z0qg_col,
        h2osoi_liq = waterstatebulk.ws.h2osoi_liq_col,
        h2osoi_ice = waterstatebulk.ws.h2osoi_ice_col,
        forc_t = forc_t, forc_th = forc_th, forc_q = forc_q, forc_pbot = forc_pbot,
        forc_rho = forc_rho, forc_lwrad = forc_lwrad,
        forc_u = forc_u, forc_v = forc_v, forc_hgt_u = forc_hgt_u_grc,
        forc_hgt_t = forc_hgt_t_grc, forc_hgt_q = forc_hgt_q_grc, grc_lat = grc_lat,
        eflx_sh_tot = energyflux.eflx_sh_tot_patch,
        eflx_sh_grnd = energyflux.eflx_sh_grnd_patch,
        eflx_lh_tot = energyflux.eflx_lh_tot_patch,
        eflx_lh_grnd = energyflux.eflx_lh_grnd_patch,
        eflx_lwrad_out = energyflux.eflx_lwrad_out_patch,
        eflx_lwrad_net = energyflux.eflx_lwrad_net_patch,
        eflx_soil_grnd = energyflux.eflx_soil_grnd_patch,
        eflx_gnet = energyflux.eflx_gnet_patch,
        taux = energyflux.taux_patch, tauy = energyflux.tauy_patch,
        qflx_evap_soi = waterfluxbulk.wf.qflx_evap_soi_patch,
        qflx_evap_tot = waterfluxbulk.wf.qflx_evap_tot_patch,
        ustar = frictionvel.ustar_patch, z0mg = frictionvel.z0mg_patch,
        z0hg = frictionvel.z0hg_patch, z0qg = frictionvel.z0qg_patch,
        t_ref2m = temperature.t_ref2m_patch,
        q_ref2m = waterdiagbulk.q_ref2m_patch,
        rh_ref2m = waterdiagbulk.rh_ref2m_patch,
        obu_out = frictionvel.obu_patch, zeta_out = frictionvel.zeta_patch)
end

# --------------------------------------------------------------------------
# Per-patch kernel: ONE THREAD PER FILTERED PATCH. p = @index(Global) over
# ndrange = length(mask_lakep). The internal stability iteration (1:niters) runs
# sequentially in-thread because it is loop-carried (obu/ustar/z0* carry between
# passes). Struct fields are aliased to locals so the loop body matches the
# original scalar code; all Float64 literals are converted to the working element
# type T so no Float64 reaches a Float32-only backend.
# --------------------------------------------------------------------------
@kernel function _lake_fluxes_kernel!(d, @Const(mask_lakep), nlevsno_val::Int,
                                      niters::Int, dtime, a_coef, a_exp, zsno)
    p = @index(Global)
    @inbounds if mask_lakep[p]
        # ----- alias struct fields to locals (body stays verbatim) -----
        p_column = d.p_column; p_gridcell = d.p_gridcell; p_landunit = d.p_landunit
        c_gridcell = d.c_gridcell
        snl_arr = d.snl; lakedepth = d.lakedepth; savedtke1 = d.savedtke1
        ws_col = d.ws_col; ks_col = d.ks_col; htvp_col = d.htvp_col; ust_lake = d.ust_lake
        dz = d.dz; dz_lake = d.dz_lake; t_soisno = d.t_soisno; t_lake = d.t_lake
        t_grnd = d.t_grnd; sabg_arr = d.sabg
        fsds_nir_d = d.fsds_nir_d; fsds_nir_i = d.fsds_nir_i
        fsr_nir_d = d.fsr_nir_d; fsr_nir_i = d.fsr_nir_i; sabg_lyr = d.sabg_lyr
        z0mg_col = d.z0mg_col; z0hg_col = d.z0hg_col; z0qg_col = d.z0qg_col
        h2osoi_liq = d.h2osoi_liq; h2osoi_ice = d.h2osoi_ice
        grc_lat = d.grc_lat
        forc_t = d.forc_t; forc_th = d.forc_th; forc_q = d.forc_q
        forc_pbot = d.forc_pbot; forc_rho = d.forc_rho; forc_lwrad = d.forc_lwrad
        forc_u = d.forc_u; forc_v = d.forc_v
        forc_hgt_u_grc = d.forc_hgt_u; forc_hgt_t_grc = d.forc_hgt_t
        forc_hgt_q_grc = d.forc_hgt_q

        T = eltype(t_grnd)

        c = p_column[p]
        g = p_gridcell[p]
        l = p_landunit[p]

        snl = snl_arr[c]

        # Wind minimum (CLM params_inst%wind_min = 1.0; Fortran LakeFluxes floors ur
        # at this, not 0.1 — the low-wind ur feeds ustar).
        wind_min = T(1.0)

        # ============================================================
        # Phase 1: Initialization
        # ============================================================

        # Fetch estimate
        lakedepth_c = lakedepth[c]
        if lakedepth_c < T(4.0)
            fetch = T(100.0)
        else
            fetch = T(25.0) * lakedepth_c
        end

        # Top layer index (Julia 1-based)
        jtop = snl + 1 + nlevsno_val

        # Solar absorbed at surface (LakeFluxesMod.F90:376-397). NOT clamped:
        # Fortran uses sabg(p) as-is.
        sabg = sabg_arr[p]

        # betaprime = fraction of sabg absorbed at the skin (LakeFluxesMod.F90:383-393):
        #   snow layers -> SNICAR top-layer absorption sabg_lyr(jtop)/sabg;
        #   no snow     -> the NIR fraction of absorbed solar, + (1-NIR)*betavis.
        # betavis is 0 in CTSM (LakeCon.F90:45), so betaprime is the NIR fraction --
        # the SAME `beta` LakeTemperature uses to put the remaining (visible) solar
        # into the water column. The port had a "simplified" constant 0.4 here,
        # inconsistent with its own lake_temperature.jl (which already computes the
        # NIR fraction), so the skin and the column double-booked daytime solar.
        if snl < 0
            betaprime = sabg_lyr[p, jtop] / max(T(1.0e-5), sabg)
            dzsur = dz[c, jtop] / T(2.0)
        else
            sabg_nir = fsds_nir_d[p] + fsds_nir_i[p] - fsr_nir_d[p] - fsr_nir_i[p]
            sabg_nir = min(sabg_nir, sabg)
            betaprime = sabg_nir / max(T(1.0e-5), sabg)
            betaprime = betaprime + (one(T) - betaprime) * T(BETAVIS_LAKE)
            dzsur = dz_lake[c, 1] / T(2.0)
        end

        # Saturated specific humidity + derivative at the surface -- CLM QSat
        # (QSatMod.F90: Flatau polynomials, ICE saturation below tfrz), F90:402.
        # The port used a Tetens water-only formula throughout this kernel, i.e.
        # liquid-water saturation over a frozen lake skin.
        (qsatg, _, qsatgdT, _) = qsat(t_grnd[c], forc_pbot[c])

        # Initial roughness lengths (LakeFluxesMod.F90:327-367, z0param_method =
        # 'ZengWang2007'; use_z0m_snowmelt = .false.):
        #   unfrozen -> the carried column z0mg_col, classic z0hg/z0qg from ust_lake;
        #   frozen, no snow -> z0frzlake + ZengWang z0hg with nu_param;
        #   snow -> zsno + ZengWang z0hg with nu_param.
        if t_grnd[c] > T(TFRZ)
            z0mg = z0mg_col[c]
            kva = T(KVA0) * (t_grnd[c] / T(KVA0TEMP))^T(1.5) * T(KVA0PRES) / forc_pbot[c]
            sqre0 = (max(z0mg * ust_lake[c] / kva, T(0.1)))^T(0.5)
            z0hg = z0mg * exp(-T(VKC) / T(PRN_AIR) * (T(4.0) * sqre0 - T(3.2)))
            z0qg = z0mg * exp(-T(VKC) / T(SCH_WATER) * (T(4.0) * sqre0 - T(4.2)))
            z0qg = max(z0qg, T(MINZ0LAKE))
            z0hg = max(z0hg, T(MINZ0LAKE))
        elseif snl == 0
            z0mg = T(z0frzlake)
            z0hg = z0mg / exp(a_coef * (ust_lake[c] * z0mg / T(NU_PARAM_LAKE))^a_exp)
            z0qg = z0hg
        else
            z0mg = zsno
            z0hg = z0mg / exp(a_coef * (ust_lake[c] * z0mg / T(NU_PARAM_LAKE))^a_exp)
            z0qg = z0hg
        end

        # Forcing heights for the roughness lengths (F90:370-372). Fortran adds z0 to
        # the atmospheric reference heights (forc_hgt_*_patch) and refreshes them at
        # the end of every pass; thm and zldis are taken ONCE, here.
        hgt_u = forc_hgt_u_grc[g] + z0mg
        hgt_t = forc_hgt_t_grc[g] + z0hg
        hgt_q = forc_hgt_q_grc[g] + z0qg
        thm = forc_t[c] + T(0.0098) * hgt_t
        thv = forc_th[c] * (one(T) + T(0.61) * forc_q[c])

        # Latent heat by the GROUND temperature (F90:426-430), re-evaluated after
        # every Newton solve and after the freeze/convection reset. The port keyed it
        # once per step on the top-layer liquid/ice content.
        htvp = t_grnd[c] > T(TFRZ) ? T(HVAP) : T(HSUB)

        # Wind speed
        ur = max(sqrt(forc_u[g]^2 + forc_v[g]^2), wind_min)
        dth = thm - t_grnd[c]
        dqh = forc_q[c] - qsatg
        dthv = dth * (one(T) + T(0.61) * forc_q[c]) + T(0.61) * forc_th[c] * dqh
        zldis = hgt_u - zero(T)

        # Displacement height (zero for lakes)
        displa = zero(T)

        # MoninObukIni (FrictionVelocityMod), bulk-Richardson form.
        ustar = T(0.06)   # overwritten by FrictionVelocity on pass 1
        if dthv >= zero(T)
            um = max(ur, T(0.1))
        else
            um = sqrt(ur^2 + T(0.25))   # convective velocity wc = 0.5 m/s
        end
        rib = T(GRAV) * zldis * dthv / (thv * um^2)
        if rib >= zero(T)
            zeta = rib * log(zldis / z0mg) / (one(T) - T(5.0) * min(rib, T(0.19)))
            zeta = min(T(ZETAMAX_LAKE), max(zeta, T(0.01)))
        else
            zeta = rib * log(zldis / z0mg)
            zeta = max(T(-100.0), min(zeta, T(-0.01)))
        end
        obu = zldis / zeta

        # ============================================================
        # Phase 2: Stability iteration (F90:439-631)
        # ============================================================
        # Fortran re-linearises the surface Newton solve EVERY pass about the latest
        # ground temperature: `tgbef(c) = t_grnd(c)` at the top of the pass, stftg3
        # from that tgbef, and qsatg/qsatgdT refreshed by QSat(t_grnd) at the end of
        # the previous pass (F90:532). The port linearised every pass about the
        # step's INPUT t_grnd with qsat frozen there, so a skin that moved several K
        # within the step carried a stale qsat/derivative into its LH and Newton terms.
        t_grnd_new = t_grnd[c]
        tgbef = t_grnd_new
        rah_c = one(T); raw_c = one(T); ram_c = one(T)
        temp1_c = one(T); temp2_c = one(T); temp12m_c = one(T); temp22m_c = one(T)
        dth_c = dth; dqh_c = dqh
        zeta_conv = zeta
        eflx_sh_grnd_p = zero(T); qflx_evap_soi_p = zero(T)
        obuold = zero(T); nmozsgn = 0

        for iter in 1:niters
            # ---- FrictionVelocity (FrictionVelocityMod.F90), all four regimes ----
            zldis_u = hgt_u - displa
            zldis_t = hgt_t - displa
            zldis_q = hgt_q - displa
            zeta_u = zldis_u / obu
            zeta_t = zldis_t / obu
            zeta_q = zldis_q / obu
            ustar = T(VKC) * um / mo_profile_denom_m(zeta_u, zldis_u, z0mg, obu)
            temp1 = T(VKC) / mo_profile_denom_h(zeta_t, zldis_t, z0hg, obu)
            temp2 = T(VKC) / mo_profile_denom_h(zeta_q, zldis_q, z0qg, obu)
            zldis_2h = T(2.0) + z0hg
            temp12m = T(VKC) / mo_profile_denom_h(zldis_2h / obu, zldis_2h, z0hg, obu)
            if z0qg == z0hg
                temp22m = temp12m
            else
                zldis_2q = T(2.0) + z0qg
                temp22m = T(VKC) / mo_profile_denom_h(zldis_2q / obu, zldis_2q, z0qg, obu)
            end

            # ---- per-pass surface state (F90:460-487) ----
            tgbef = t_grnd_new
            if t_grnd_new > T(TFRZ) && t_lake[c, 1] > T(TFRZ) && snl == 0
                tksur = savedtke1[c]
                tsur = t_lake[c, 1]
            elseif snl == 0
                tksur = T(TKICE)
                tsur = t_lake[c, 1]
            else
                bw = (h2osoi_ice[c, jtop] + h2osoi_liq[c, jtop]) / dz[c, jtop]
                tksur = T(TKAIR) + (T(7.75e-5) * bw + T(1.105e-6) * bw * bw) * (T(TKICE) - T(TKAIR))
                tsur = t_soisno[c, jtop]
            end

            # Aerodynamic resistances (F90:489-491)
            ram = one(T) / (ustar * ustar / um)
            rah = one(T) / (temp1 * ustar)
            raw = one(T) / (temp2 * ustar)

            # Newton solve for the skin temperature (F90:500-514)
            stftg3 = T(EMG_LAKE) * T(SB) * tgbef * tgbef * tgbef
            ax = betaprime * sabg + T(EMG_LAKE) * forc_lwrad[c] + T(3.0) * stftg3 * tgbef +
                 forc_rho[c] * T(CPAIR) / rah * thm -
                 htvp * forc_rho[c] / raw * (qsatg - qsatgdT * tgbef - forc_q[c]) +
                 tksur * tsur / dzsur
            bx = T(4.0) * stftg3 + forc_rho[c] * T(CPAIR) / rah +
                 htvp * forc_rho[c] / raw * qsatgdT + tksur / dzsur
            t_grnd_new = ax / bx

            htvp = t_grnd_new > T(TFRZ) ? T(HVAP) : T(HSUB)

            # Surface fluxes, linearised about this pass's tgbef (F90:526-527)
            eflx_sh_grnd_p = forc_rho[c] * T(CPAIR) * (t_grnd_new - thm) / rah
            qflx_evap_soi_p = forc_rho[c] * (qsatg + qsatgdT * (t_grnd_new - tgbef) - forc_q[c]) / raw

            # Refresh qsat at the new skin temperature (F90:532-536)
            (qsatg, _, qsatgdT, _) = qsat(t_grnd_new, forc_pbot[c])
            dth = thm - t_grnd_new
            dqh = forc_q[c] - qsatg

            tstar = temp1 * dth
            qstar = temp2 * dqh
            thvstar = tstar * (one(T) + T(0.61) * forc_q[c]) + T(0.61) * forc_th[c] * qstar
            zeta = zldis * T(VKC) * T(GRAV) * thvstar / (ustar^2 * thv)

            if zeta >= zero(T)
                zeta = min(T(ZETAMAX_LAKE), max(zeta, T(0.01)))
                um = max(ur, T(0.1))
            else
                zeta = max(T(-100.0), min(zeta, T(-0.01)))
                wc = T(BETA1) * (-T(GRAV) * ustar * thvstar * T(ZII_LAKE) / thv)^T(0.333)
                um = sqrt(ur * ur + wc * wc)
            end
            obu = zldis / zeta
            obuold * obu < zero(T) && (nmozsgn += 1)
            obuold = obu

            zeta_conv = zeta
            rah_c = rah; raw_c = raw; ram_c = ram
            temp1_c = temp1; temp2_c = temp2; temp12m_c = temp12m; temp22m_c = temp22m
            dth_c = dth; dqh_c = dqh

            # Roughness update from the just-solved skin temperature (F90:558-617)
            if t_grnd_new > T(TFRZ) && snl == 0
                cur = T(CUR0) + T(CURM) * exp(max(
                          -(fetch * T(GRAV) / ustar / ustar)^(one(T) / T(3.0)) / T(FCRIT),
                          -(lakedepth_c * T(GRAV) / ur / ur)^T(0.5)))
                kva = T(KVA0) * (t_grnd_new / T(KVA0TEMP))^T(1.5) * T(KVA0PRES) / forc_pbot[c]
                z0mg = max(T(CUS) * kva / max(ustar, T(1.0e-4)), cur * ustar * ustar / T(GRAV))
                z0mg = max(z0mg, T(MINZ0LAKE))
                sqre0 = (max(z0mg * ustar / kva, T(0.1)))^T(0.5)
                z0hg = z0mg * exp(-T(VKC) / T(PRN_AIR) * (T(4.0) * sqre0 - T(3.2)))
                z0qg = z0mg * exp(-T(VKC) / T(SCH_WATER) * (T(4.0) * sqre0 - T(4.2)))
                z0qg = max(z0qg, T(MINZ0LAKE))
                z0hg = max(z0hg, T(MINZ0LAKE))
            elseif snl == 0
                z0mg = T(z0frzlake)
                z0hg = z0mg / exp(a_coef * (ustar * z0mg / T(NU_PARAM_LAKE))^a_exp)
                z0qg = z0hg
            else
                # snow layers: z0mg keeps its pre-loop zsno (use_z0m_snowmelt=.false.)
                z0hg = z0mg / exp(a_coef * (ustar * z0mg / T(NU_PARAM_LAKE))^a_exp)
                z0qg = z0hg
            end
            hgt_u = forc_hgt_u_grc[g] + z0mg
            hgt_t = forc_hgt_t_grc[g] + z0hg
            hgt_q = forc_hgt_q_grc[g] + z0qg

            # Fortran drops a patch from further passes once obu has changed sign
            # three times (F90:625-631).
            nmozsgn >= 3 && break
        end

        # ============================================================
        # Phase 3: freeze / convective-mixing resets (F90:654-674). The fluxes are
        # re-evaluated ONLY when t_grnd is reset, linearised about the pre-reset skin
        # temperature with the qsat refreshed there; otherwise the last pass's
        # fluxes stand. (The port recomputed both fluxes from a nonlinear qsat at the
        # final temperature on every step.)
        # ============================================================
        if (snl < 0 || t_lake[c, 1] <= T(TFRZ)) && t_grnd_new > T(TFRZ)
            t_grnd_temp = t_grnd_new
            t_grnd_new = T(TFRZ)
            eflx_sh_grnd_p = forc_rho[c] * T(CPAIR) * (t_grnd_new - thm) / rah_c
            qflx_evap_soi_p = forc_rho[c] * (qsatg + qsatgdT * (t_grnd_new - t_grnd_temp) - forc_q[c]) / raw_c
        elseif (t_lake[c, 1] > t_grnd_new && t_grnd_new > T(TDMAX_LAKE)) ||
               (t_lake[c, 1] < t_grnd_new && t_lake[c, 1] > T(TFRZ) && t_grnd_new < T(TDMAX_LAKE))
            t_grnd_temp = t_grnd_new
            t_grnd_new = t_lake[c, 1]
            eflx_sh_grnd_p = forc_rho[c] * T(CPAIR) * (t_grnd_new - thm) / rah_c
            qflx_evap_soi_p = forc_rho[c] * (qsatg + qsatgdT * (t_grnd_new - t_grnd_temp) - forc_q[c]) / raw_c
        end
        htvp = t_grnd_new > T(TFRZ) ? T(HVAP) : T(HSUB)

        t_grnd[c] = t_grnd_new

        # ============================================================
        # Phase 4: final fluxes (F90:685-705)
        # ============================================================
        eflx_lwrad_out_p = (one(T) - T(EMG_LAKE)) * forc_lwrad[c] + T(EMG_LAKE) * T(SB) * t_grnd_new^4

        # eflx_soil_grnd uses the FULL sabg (F90:691) so the energy-balance check
        # closes; eflx_gnet (F90:743) is the skin flux into the lake, betaprime*sabg.
        eflx_soil_grnd_p = sabg + forc_lwrad[c] - eflx_lwrad_out_p -
                           eflx_sh_grnd_p - htvp * qflx_evap_soi_p

        d.eflx_sh_tot[p] = eflx_sh_grnd_p
        d.eflx_sh_grnd[p] = eflx_sh_grnd_p
        d.eflx_lh_tot[p] = htvp * qflx_evap_soi_p
        d.eflx_lh_grnd[p] = htvp * qflx_evap_soi_p
        d.eflx_lwrad_out[p] = eflx_lwrad_out_p
        # net longwave is upward-positive [+ = to atm] = out - in (LakeFluxesMod.F90:772).
        d.eflx_lwrad_net[p] = eflx_lwrad_out_p - forc_lwrad[c]
        d.eflx_soil_grnd[p] = eflx_soil_grnd_p
        # eflx_gnet: lake_temperature reads it as `fin`, its top-layer BC.
        d.eflx_gnet[p] = betaprime * sabg + forc_lwrad[c] -
                         (eflx_lwrad_out_p + eflx_sh_grnd_p + htvp * qflx_evap_soi_p)

        d.qflx_evap_soi[p] = qflx_evap_soi_p
        d.qflx_evap_tot[p] = qflx_evap_soi_p

        # ---- 2-m diagnostics (LakeFluxesMod.F90:707-717) ----
        # DEAD WRITE until now: nothing wrote t_ref2m_patch on a lake patch, so TSA
        # on a 100%-lake gridcell reported the untouched cold-start 283.000 K at
        # every step. Writing it here is lake-landunit-gated by construction — this
        # kernel only ever runs over `mask_lakep`, and no other landunit's 2-m
        # diagnostic reaches this code path.
        d.t_ref2m[p] = thm + temp1_c * dth_c * (one(T) / temp12m_c - one(T) / temp1_c)
        q_ref2m_p = forc_q[c] + temp2_c * dqh_c * (one(T) / temp22m_c - one(T) / temp2_c)
        d.q_ref2m[p] = q_ref2m_p
        (qsat_ref2m, _, _, _) = qsat(d.t_ref2m[p], forc_pbot[c])
        d.rh_ref2m[p] = min(T(100.0), q_ref2m_p / qsat_ref2m * T(100.0))

        d.ustar[p] = ustar
        d.obu_out[p] = obu
        d.zeta_out[p] = zeta_conv
        d.z0mg[p] = z0mg
        d.z0hg[p] = z0hg
        d.z0qg[p] = z0qg

        # Momentum stress
        ram_final = ram_c
        d.taux[p] = -forc_rho[c] * forc_u[g] / ram_final
        d.tauy[p] = -forc_rho[c] * forc_v[g] / ram_final

        # 2m wind speed for mixing parameters
        u2m = max(T(0.1), ustar / T(VKC) * log(T(2.0) / z0mg))
        ust_lake[c] = ustar     # carry the friction velocity for next step's z0hg (LakeFluxesMod:761)
        ws_col[c] = T(1.2e-3) * u2m
        # Wave-driven eddy-extinction coefficient (LakeFluxesMod.F90:755)
        #   ks = 6.6 * sqrt(|sin(lat)|) * u2m^-1.84      [1/m], lat in RADIANS
        # `ks` is the DEPTH-DECAY RATE of the wind-driven eddy diffusivity in
        # LakeTemperature: ke = vkc*ws*z/p0 * exp(-ks*z) / (1+37*ri^2). It was
        # hardwired to zero here because the gridcell latitude had never been
        # wired into the lake kernel — and `exp(-0*z) == 1` means the surface
        # eddy diffusivity NEVER DECAYS WITH DEPTH: `ke` then grows linearly all
        # the way to the lake bed, coupling the whole column to the skin. At Bow
        # (lat 51.2 deg, u2m ~1.6 m/s) the true ks is ~2.4 /m, so at z = 5 m the
        # port's diffusivity was exp(2.4*5) = 1.6e5 times too large. That is why
        # the surface could never cool to freezing: every W/m2 the surface lost
        # was replenished from 50 m of 277 K water instead of from the top ~0.3 m.
        #
        # AD note: the previous literal zero was justified by sqrt(0) having an
        # Inf ForwardDiff derivative. That hazard is real but only at |sin(lat)|
        # == 0 exactly (the equator), so branch on it rather than delete the term.
        sinlat = abs(sin(grc_lat[g]))
        ks_col[c] = sinlat > zero(T) ? T(6.6) * sqrt(sinlat) * u2m^T(-1.84) : zero(T)

        htvp_col[c] = htvp
        # Carry the column roughness lengths to the next step (F90:764-766).
        z0mg_col[c] = z0mg
        z0hg_col[c] = z0hg
        z0qg_col[c] = z0qg
    end
end

"""
    lake_fluxes!(temperature, energyflux, frictionvel, solarabs,
                  lakestate, waterstatebulk, waterdiagbulk,
                  waterfluxbulk, col_data, patch_data, lun_data,
                  forc_t, forc_th, forc_q, forc_pbot,
                  forc_rho, forc_lwrad, forc_u, forc_v,
                  forc_hgt_u_grc, forc_hgt_t_grc, forc_hgt_q_grc,
                  mask_lakec, mask_lakep,
                  bounds_col, bounds_patch;
                  dtime)

Calculate surface fluxes and surface temperature for lake columns.

Uses iterative Monin-Obukhov similarity theory with Newton-Raphson
surface temperature solution.

Ported from `LakeFluxes` in `LakeFluxesMod.F90`. The per-patch loop runs as a
single KernelAbstractions kernel (one thread per filtered patch, internal
sequential stability iteration), so the whole function executes on the GPU when
the state lives on the device.
"""
function lake_fluxes!(temperature::TemperatureData,
                       energyflux::EnergyFluxData,
                       frictionvel::FrictionVelocityData,
                       solarabs::SolarAbsorbedData,
                       lakestate::LakeStateData,
                       waterstatebulk::WaterStateBulkData,
                       waterdiagbulk::WaterDiagnosticBulkData,
                       waterfluxbulk::WaterFluxBulkData,
                       col_data::ColumnData,
                       patch_data::PatchData,
                       lun_data::LandunitData,
                       forc_t::AbstractVector{<:Real},
                       forc_th::AbstractVector{<:Real},
                       forc_q::AbstractVector{<:Real},
                       forc_pbot::AbstractVector{<:Real},
                       forc_rho::AbstractVector{<:Real},
                       forc_lwrad::AbstractVector{<:Real},
                       forc_u::AbstractVector{<:Real},
                       forc_v::AbstractVector{<:Real},
                       forc_hgt_u_grc::AbstractVector{<:Real},
                       forc_hgt_t_grc::AbstractVector{<:Real},
                       forc_hgt_q_grc::AbstractVector{<:Real},
                       mask_lakec::AbstractVector{Bool},
                       mask_lakep::AbstractVector{Bool},
                       bounds_col::UnitRange{Int},
                       bounds_patch::UnitRange{Int};
                       dtime::Real = 1800.0,
                       grc_lat::Union{Nothing,AbstractVector{<:Real}} = nothing)

    nlevsno_val = varpar.nlevsno
    niters = 4  # stability iterations

    np = length(mask_lakep)
    np == 0 && return nothing

    dv = _lake_flux_dv(temperature, energyflux, frictionvel, solarabs, lakestate,
                       waterstatebulk, waterdiagbulk, waterfluxbulk, col_data, patch_data,
                       forc_t, forc_th, forc_q, forc_pbot, forc_rho, forc_lwrad,
                       forc_u, forc_v, forc_hgt_u_grc, forc_hgt_t_grc, forc_hgt_q_grc,
                       _lake_lat_like(forc_hgt_u_grc, grc_lat))

    FT = eltype(temperature.t_grnd_col)
    # Struct-first kernel (the DV is passed whole as the first arg): pick the
    # backend off a real state array, launch directly, then synchronize — matches
    # the photosynthesis/_psn_dv struct-first launch convention.
    be = _kernel_backend(temperature.t_grnd_col)
    # a_coef/a_exp = the BareGround/ZengWang2007 drag params (clm5_params.nc), shared with the
    # canopy roughness formula. Used for the FROZEN-lake heat-roughness z0hg (LakeFluxesMod 'ZengWang2007').
    _lake_fluxes_kernel!(be)(dv, mask_lakep, nlevsno_val, niters, FT(dtime),
                             FT(canopy_fluxes_params.a_coef), FT(canopy_fluxes_params.a_exp),
                             FT(frictionvel.zsno);
                             ndrange = np)
    KA.synchronize(be)

    return nothing
end

# Helper: saturation specific humidity for liquid water (eltype-generic so it
# lowers to valid Float32 Metal IR; byte-identical to the Float64 literals on CPU).
function qsat_water(t::Real, p::Real)
    T = promote_type(typeof(t), typeof(p))
    tt = T(t); pp = T(p)
    # Tetens formula for saturation vapor pressure over water
    es = T(611.2) * exp(T(17.67) * (tt - T(TFRZ)) / (tt - T(TFRZ) + T(243.5)))
    return T(0.622) * es / (pp - T(0.378) * es)
end

# Helper: derivative of qsat w.r.t. temperature
function qsat_water_dT(t::Real, p::Real)
    T = promote_type(typeof(t), typeof(p))
    tt = T(t); pp = T(p)
    es = T(611.2) * exp(T(17.67) * (tt - T(TFRZ)) / (tt - T(TFRZ) + T(243.5)))
    desdT = es * T(17.67) * T(243.5) / (tt - T(TFRZ) + T(243.5))^2
    qs = T(0.622) * es / (pp - T(0.378) * es)
    return T(0.622) * pp * desdT / (pp - T(0.378) * es)^2
end
