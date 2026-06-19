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
const BETAVIS_LAKE = 0.4       # fraction of visible radiation absorbed at surface
const MINZ0LAKE = 1.0e-5       # minimum lake roughness length [m]
const TDMAX_LAKE = 277.0       # temperature of maximum water density [K]

# Aerodynamic parameters
const PRN_AIR = 0.713          # Prandtl number for air
const SCH_WATER = 0.66         # Schmidt number for water in air
const KVA0 = 1.51e-5           # kinematic viscosity of air at 20°C [m²/s]
const CUS = 0.1                # smooth flow regime coefficient
const CUR0 = 0.01              # base Charnock constant
const CURM = 0.0                # modified Charnock coefficient
const BETA1 = 1.0              # coefficient of convective velocity
const ZETAMAX_LAKE = 0.5       # max zeta under stable conditions

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
    snl::Vi; lakedepth::V; savedtke1::V; ws_col::V; ks_col::V; htvp_col::V
    # column matrices
    dz::M; dz_lake::M; t_soisno::M; t_lake::M
    # column vectors
    t_grnd::V
    # patch solar
    sabg::V
    # water state (column, layer)
    h2osoi_liq::M; h2osoi_ice::M
    # forcings (column-indexed)
    forc_t::V; forc_th::V; forc_q::V; forc_pbot::V; forc_rho::V; forc_lwrad::V
    # forcings (gridcell-indexed)
    forc_u::V; forc_v::V; forc_hgt_u::V; forc_hgt_t::V; forc_hgt_q::V
    # outputs (patch)
    eflx_sh_tot::V; eflx_sh_grnd::V; eflx_lh_tot::V; eflx_lh_grnd::V
    eflx_lwrad_out::V; eflx_lwrad_net::V; eflx_soil_grnd::V; taux::V; tauy::V
    qflx_evap_soi::V; qflx_evap_tot::V; ustar::V; z0mg::V; z0hg::V; z0qg::V
end
Adapt.@adapt_structure LakeFluxDV

# Build a LakeFluxDV from the caller's state structs (array fields are shared refs).
function _lake_flux_dv(temperature, energyflux, frictionvel, solarabs, lakestate,
                       waterstatebulk, waterfluxbulk, col_data, patch_data,
                       forc_t, forc_th, forc_q, forc_pbot, forc_rho, forc_lwrad,
                       forc_u, forc_v, forc_hgt_u_grc, forc_hgt_t_grc, forc_hgt_q_grc)
    return LakeFluxDV(;
        p_column = patch_data.column, p_gridcell = patch_data.gridcell,
        p_landunit = patch_data.landunit, c_gridcell = col_data.gridcell,
        snl = col_data.snl, lakedepth = col_data.lakedepth,
        savedtke1 = lakestate.savedtke1_col, ws_col = lakestate.ws_col,
        ks_col = lakestate.ks_col, htvp_col = energyflux.htvp_col,
        dz = col_data.dz, dz_lake = col_data.dz_lake,
        t_soisno = temperature.t_soisno_col, t_lake = temperature.t_lake_col,
        t_grnd = temperature.t_grnd_col, sabg = solarabs.sabg_patch,
        h2osoi_liq = waterstatebulk.ws.h2osoi_liq_col,
        h2osoi_ice = waterstatebulk.ws.h2osoi_ice_col,
        forc_t = forc_t, forc_th = forc_th, forc_q = forc_q, forc_pbot = forc_pbot,
        forc_rho = forc_rho, forc_lwrad = forc_lwrad,
        forc_u = forc_u, forc_v = forc_v, forc_hgt_u = forc_hgt_u_grc,
        forc_hgt_t = forc_hgt_t_grc, forc_hgt_q = forc_hgt_q_grc,
        eflx_sh_tot = energyflux.eflx_sh_tot_patch,
        eflx_sh_grnd = energyflux.eflx_sh_grnd_patch,
        eflx_lh_tot = energyflux.eflx_lh_tot_patch,
        eflx_lh_grnd = energyflux.eflx_lh_grnd_patch,
        eflx_lwrad_out = energyflux.eflx_lwrad_out_patch,
        eflx_lwrad_net = energyflux.eflx_lwrad_net_patch,
        eflx_soil_grnd = energyflux.eflx_soil_grnd_patch,
        taux = energyflux.taux_patch, tauy = energyflux.tauy_patch,
        qflx_evap_soi = waterfluxbulk.wf.qflx_evap_soi_patch,
        qflx_evap_tot = waterfluxbulk.wf.qflx_evap_tot_patch,
        ustar = frictionvel.ustar_patch, z0mg = frictionvel.z0mg_patch,
        z0hg = frictionvel.z0hg_patch, z0qg = frictionvel.z0qg_patch)
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
                                      niters::Int, dtime)
    p = @index(Global)
    @inbounds if mask_lakep[p]
        # ----- alias struct fields to locals (body stays verbatim) -----
        p_column = d.p_column; p_gridcell = d.p_gridcell; p_landunit = d.p_landunit
        c_gridcell = d.c_gridcell
        snl_arr = d.snl; lakedepth = d.lakedepth; savedtke1 = d.savedtke1
        ws_col = d.ws_col; ks_col = d.ks_col; htvp_col = d.htvp_col
        dz = d.dz; dz_lake = d.dz_lake; t_soisno = d.t_soisno; t_lake = d.t_lake
        t_grnd = d.t_grnd; sabg_arr = d.sabg
        h2osoi_liq = d.h2osoi_liq; h2osoi_ice = d.h2osoi_ice
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
        if snl < 0
            dzsur = dz[c, jtop] / T(2.0)
        else
            dzsur = dz_lake[c, 1] / T(2.0)
        end

        # Thermal conductivity of surface layer
        tksur = savedtke1[c]

        # Temperature of subsurface layer
        if snl < 0
            # Snow present: subsurface is next layer
            tsur = t_soisno[c, jtop + 1]
        else
            # No snow: second lake layer
            if size(t_lake, 2) >= 2
                tsur = t_lake[c, 2]
            else
                tsur = t_lake[c, 1]
            end
        end

        # Saved ground temperature
        tgbef = t_grnd[c]

        # Solar absorbed at surface
        sabg = sabg_arr[p]
        sabg = smooth_max(sabg, zero(T))

        # Betaprime: fraction of solar absorbed at surface
        if snl < 0
            betaprime = one(T)  # all absorbed in snow
        else
            betaprime = T(BETAVIS_LAKE)  # simplified
        end

        # Virtual potential temperature
        thv = forc_th[c] * (one(T) + T(0.61) * forc_q[c])

        # Wind speed
        ur = smooth_max(sqrt(forc_u[g]^2 + forc_v[g]^2), wind_min)

        # Forcing heights
        forc_hgt_u = forc_hgt_u_grc[g]
        forc_hgt_t = forc_hgt_t_grc[g]
        forc_hgt_q = forc_hgt_q_grc[g]

        # Intermediate temperature
        thm = forc_t[c] + T(0.0098) * forc_hgt_t

        # Kinematic viscosity
        kva = T(KVA0) * (t_grnd[c] / T(293.15))^T(1.5) * (T(1013.25e2) / forc_pbot[c])

        # Latent heat
        h2o_liq = h2osoi_liq[c, jtop]
        h2o_ice = h2osoi_ice[c, jtop]
        htvp = (h2o_liq <= zero(T) && h2o_ice > zero(T)) ? T(HSUB) : T(HVAP)

        # ============================================================
        # Initial roughness lengths
        # ============================================================
        if tgbef > T(TFRZ)  # Unfrozen
            z0mg = smooth_max(T(MINZ0LAKE), T(CUS) * kva / smooth_max(ur * T(0.1), T(1e-4)))
        else  # Frozen
            if snl < 0
                z0mg = T(0.00085)  # snow roughness
            else
                z0mg = T(0.001)  # ice roughness
            end
        end

        sqre0 = sqrt(smooth_max(z0mg * ur * T(0.1) / kva, T(0.1)))
        z0hg = z0mg * exp(-T(VKC) / T(PRN_AIR) * (T(4.0) * sqre0 - T(3.2)))
        z0qg = z0mg * exp(-T(VKC) / T(SCH_WATER) * (T(4.0) * sqre0 - T(4.2)))
        z0mg = smooth_max(z0mg, T(1.0e-10))
        z0hg = smooth_max(z0hg, T(1.0e-10))
        z0qg = smooth_max(z0qg, T(1.0e-10))

        # Reference displacement height (zero for lakes)
        displa = zero(T)

        # Monin-Obukhov initialization
        zldis = forc_hgt_u - displa
        zldis = smooth_max(zldis, z0mg + T(0.01))
        dth = thm - tgbef
        dqh = forc_q[c] - qsat_water(tgbef, forc_pbot[c])
        dthv = dth * (one(T) + T(0.61) * forc_q[c]) + T(0.61) * forc_th[c] * dqh

        # Initial Monin-Obukhov length + wind speed — CLM FrictionVelocityMod
        # MoninObukIni (bulk-Richardson form). The previous initialization used
        # obu = -ustar^3*thv/(vkc*grav*thvstar), which gives the WRONG SIGN for
        # UNSTABLE conditions (it returned positive obu / the stable branch when the
        # lake surface is warmer than the air, dthv<0). That collapsed ustar to ~0.01
        # (rah ~ 1e3-1e4), drove the turbulent fluxes to ~0, and let the surface
        # radiatively over-cool (TG 250 K vs Fortran 271 K). The Richardson form gives
        # obu<0 + a convective wind speed for unstable, matching Fortran.
        ustar = smooth_max(T(VKC) * ur / log(zldis / z0mg), T(0.001))  # overwritten in the iteration
        if dthv >= zero(T)
            um = smooth_max(ur, T(0.1))
        else
            um = sqrt(ur^2 + T(0.25))   # convective velocity wc = 0.5 m/s
        end
        rib = T(GRAV) * zldis * dthv / (thv * um^2)
        if rib >= zero(T)
            zeta = rib * log(zldis / z0mg) / (one(T) - T(5.0) * min(rib, T(0.19)))
            zeta = clamp(zeta, T(0.01), T(ZETAMAX_LAKE))
        else
            zeta = rib * log(zldis / z0mg)
            zeta = clamp(zeta, T(-100.0), T(-0.01))
        end
        obu = zldis / zeta

        # ============================================================
        # Phase 2: Stability iteration
        # ============================================================
        t_grnd_new = tgbef
        rah_c = one(T); raw_c = one(T); ram_c = one(T)   # converged resistances (set in loop)

        for iter in 1:niters
            # Aerodynamic resistances
            zldis_u = forc_hgt_u - displa
            zldis_t = forc_hgt_t - displa
            zldis_q = forc_hgt_q - displa
            zldis_u = smooth_max(zldis_u, z0mg + T(0.01))
            zldis_t = smooth_max(zldis_t, z0hg + T(0.01))
            zldis_q = smooth_max(zldis_q, z0qg + T(0.01))

            # Stability functions
            zeta_u = zldis_u / obu
            zeta_t = zldis_t / obu
            zeta_q = zldis_q / obu
            zeta0m = z0mg / obu
            zeta0h = z0hg / obu
            zeta0q = z0qg / obu

            # Clamp zeta values
            zeta_u = clamp(zeta_u, T(-100.0), T(ZETAMAX_LAKE))
            zeta_t = clamp(zeta_t, T(-100.0), T(ZETAMAX_LAKE))
            zeta_q = clamp(zeta_q, T(-100.0), T(ZETAMAX_LAKE))
            zeta0m = clamp(zeta0m, T(-100.0), T(ZETAMAX_LAKE))
            zeta0h = clamp(zeta0h, T(-100.0), T(ZETAMAX_LAKE))
            zeta0q = clamp(zeta0q, T(-100.0), T(ZETAMAX_LAKE))

            # Friction velocity
            ustar = T(VKC) * um / (log(zldis_u / z0mg) - stability_func1(zeta_u) + stability_func1(zeta0m))
            ustar = smooth_max(ustar, T(0.001))

            # Transfer coefficients
            temp1 = T(VKC) / (log(zldis_t / z0hg) - stability_func2(zeta_t) + stability_func2(zeta0h))
            temp2 = T(VKC) / (log(zldis_q / z0qg) - stability_func2(zeta_q) + stability_func2(zeta0q))

            # Aerodynamic resistances
            ram = one(T) / (ustar * T(VKC) / (log(zldis_u / z0mg) - stability_func1(zeta_u) + stability_func1(zeta0m)))
            rah = one(T) / (temp1 * ustar)
            raw = one(T) / (temp2 * ustar)
            ram = smooth_max(ram, one(T))
            rah = smooth_max(rah, one(T))
            raw = smooth_max(raw, one(T))

            # Saturation specific humidity at surface
            qsatg = qsat_water(tgbef, forc_pbot[c])
            qsatgdT = qsat_water_dT(tgbef, forc_pbot[c])

            # Newton-Raphson for surface temperature
            stftg3 = T(EMG_LAKE) * T(SB) * tgbef^3

            ax = betaprime * sabg +
                 T(EMG_LAKE) * forc_lwrad[c] +
                 T(3.0) * stftg3 * tgbef +
                 forc_rho[c] * T(CPAIR) / rah * thm +
                 tksur * tsur / dzsur -
                 htvp * forc_rho[c] / raw * (qsatg - qsatgdT * tgbef - forc_q[c])

            bx = T(4.0) * stftg3 +
                 forc_rho[c] * T(CPAIR) / rah +
                 htvp * forc_rho[c] / raw * qsatgdT +
                 tksur / dzsur

            bx = smooth_max(bx, T(1.0e-10))
            t_grnd_new = ax / bx

            # Surface fluxes
            eflx_sh = forc_rho[c] * T(CPAIR) * (t_grnd_new - thm) / rah
            qflx_evap = forc_rho[c] * (qsatg + qsatgdT * (t_grnd_new - tgbef) - forc_q[c]) / raw

            # Update Obukhov length
            dth_new = thm - t_grnd_new
            dqh_new = forc_q[c] - (qsatg + qsatgdT * (t_grnd_new - tgbef))
            tstar = temp1 * dth_new
            qstar = temp2 * dqh_new
            thvstar = tstar * (one(T) + T(0.61) * forc_q[c]) + T(0.61) * forc_th[c] * qstar

            if abs(thvstar) > T(1e-10)
                zeta_val = zldis_u * T(VKC) * T(GRAV) * thvstar / (ustar^2 * thv)
            else
                zeta_val = zero(T)
            end

            if zeta_val >= zero(T)
                zeta_val = clamp(zeta_val, T(0.01), T(ZETAMAX_LAKE))
                um = smooth_max(ur, wind_min)
            else
                zeta_val = clamp(zeta_val, T(-100.0), T(-0.01))
                wc = T(BETA1) * (-T(GRAV) * ustar * thvstar * T(1000.0) / thv)^(one(T)/T(3.0))
                um = sqrt(ur^2 + wc^2)
            end
            obu = zldis_u / zeta_val
            obu = clamp(obu, T(-1e4), T(1e4))
            rah_c = rah; raw_c = raw; ram_c = ram   # carry converged (log + stability) resistances

            # Update roughness for unfrozen lakes
            if tgbef > T(TFRZ)
                z0mg = smooth_max(smooth_max(T(MINZ0LAKE), T(CUS) * kva / smooth_max(ustar, T(1e-6))),
                           T(CUR0) * ustar^2 / T(GRAV))
                sqre0 = sqrt(smooth_max(z0mg * ustar / kva, T(0.1)))
                z0hg = z0mg * exp(-T(VKC) / T(PRN_AIR) * (T(4.0) * sqre0 - T(3.2)))
                z0qg = z0mg * exp(-T(VKC) / T(SCH_WATER) * (T(4.0) * sqre0 - T(4.2)))
                z0mg = smooth_max(z0mg, T(1.0e-10))
                z0hg = smooth_max(z0hg, T(1.0e-10))
                z0qg = smooth_max(z0qg, T(1.0e-10))
            end
        end

        # ============================================================
        # Phase 3: Temperature corrections
        # ============================================================

        # Prevent unfreezing when lake/snow is frozen
        if (snl < 0 || t_lake[c, 1] <= T(TFRZ)) && t_grnd_new > T(TFRZ)
            t_grnd_new = T(TFRZ)
        end

        # Convective mixing correction
        if t_lake[c, 1] > t_grnd_new && t_grnd_new > T(TDMAX_LAKE)
            t_grnd_new = t_lake[c, 1]
        elseif t_lake[c, 1] < t_grnd_new &&
               t_lake[c, 1] > T(TFRZ) && t_grnd_new < T(TDMAX_LAKE)
            t_grnd_new = t_lake[c, 1]
        end

        # Update ground temperature
        t_grnd[c] = t_grnd_new

        # ============================================================
        # Phase 4: Final flux calculations
        # ============================================================

        # Recalculate fluxes at final temperature
        qsatg_final = qsat_water(t_grnd_new, forc_pbot[c])
        thm_local = forc_t[c] + T(0.0098) * forc_hgt_t_grc[g]

        # Use the converged resistances from the stability iteration (log profile +
        # stability corrections). The previous rah_final = zldis_t/(vkc*ustar) dropped
        # the log(z/z0) factor and the stability functions, inflating the resistance.
        rah_final = rah_c
        raw_final = raw_c

        eflx_sh_grnd_p = forc_rho[c] * T(CPAIR) * (t_grnd_new - thm_local) / rah_final
        qflx_evap_soi_p = forc_rho[c] * (qsatg_final - forc_q[c]) / raw_final

        # Outgoing longwave
        eflx_lwrad_out_p = (one(T) - T(EMG_LAKE)) * forc_lwrad[c] + T(EMG_LAKE) * T(SB) * t_grnd_new^4

        # Ground heat flux (energy balance residual)
        eflx_soil_grnd_p = betaprime * sabg + forc_lwrad[c] - eflx_lwrad_out_p -
                         eflx_sh_grnd_p - htvp * qflx_evap_soi_p

        # Store results
        d.eflx_sh_tot[p] = eflx_sh_grnd_p
        d.eflx_sh_grnd[p] = eflx_sh_grnd_p
        d.eflx_lh_tot[p] = htvp * qflx_evap_soi_p
        d.eflx_lh_grnd[p] = htvp * qflx_evap_soi_p
        d.eflx_lwrad_out[p] = eflx_lwrad_out_p
        d.eflx_lwrad_net[p] = forc_lwrad[c] - eflx_lwrad_out_p
        d.eflx_soil_grnd[p] = eflx_soil_grnd_p

        d.qflx_evap_soi[p] = qflx_evap_soi_p
        d.qflx_evap_tot[p] = qflx_evap_soi_p

        d.ustar[p] = ustar
        d.z0mg[p] = z0mg
        d.z0hg[p] = z0hg
        d.z0qg[p] = z0qg

        # Momentum stress
        ram_final = ram_c
        d.taux[p] = -forc_rho[c] * forc_u[g] / ram_final
        d.tauy[p] = -forc_rho[c] * forc_v[g] / ram_final

        # 2m wind speed for mixing parameters
        u2m = smooth_max(T(0.1), ustar / T(VKC) * log(T(2.0) / z0mg))
        ws_col[c] = T(1.2e-3) * u2m
        # Wave-driven eddy-extinction coefficient ks = 6.6·sqrt(|sin(lat)|)·u2m^-1.84.
        # The gridcell latitude is not yet wired in (simplified to 0), so the sqrt(|sin 0|)
        # factor — and hence the whole term — is exactly 0. It MUST be written as a literal
        # zero, not sqrt(smooth_abs(sin(zero(T))))·…: sqrt(0) has an Inf ForwardDiff
        # derivative, so even with a zero-valued input the AD partial is Inf·0 = NaN, which
        # poisons ks → the lake eddy diffusivity kme → the entire coupled lake temperature
        # solve under AD (finite values, NaN gradients). Value is unchanged (ks = 0).
        ks_col[c] = zero(T)

        htvp_col[c] = htvp
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
                       dtime::Real = 1800.0)

    nlevsno_val = varpar.nlevsno
    niters = 4  # stability iterations

    np = length(mask_lakep)
    np == 0 && return nothing

    dv = _lake_flux_dv(temperature, energyflux, frictionvel, solarabs, lakestate,
                       waterstatebulk, waterfluxbulk, col_data, patch_data,
                       forc_t, forc_th, forc_q, forc_pbot, forc_rho, forc_lwrad,
                       forc_u, forc_v, forc_hgt_u_grc, forc_hgt_t_grc, forc_hgt_q_grc)

    FT = eltype(temperature.t_grnd_col)
    # Struct-first kernel (the DV is passed whole as the first arg): pick the
    # backend off a real state array, launch directly, then synchronize — matches
    # the photosynthesis/_psn_dv struct-first launch convention.
    be = _kernel_backend(temperature.t_grnd_col)
    _lake_fluxes_kernel!(be)(dv, mask_lakep, nlevsno_val, niters, FT(dtime);
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
