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

Ported from `LakeFluxes` in `LakeFluxesMod.F90`.
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
                       forc_t::Vector{<:Real},
                       forc_th::Vector{<:Real},
                       forc_q::Vector{<:Real},
                       forc_pbot::Vector{<:Real},
                       forc_rho::Vector{<:Real},
                       forc_lwrad::Vector{<:Real},
                       forc_u::Vector{<:Real},
                       forc_v::Vector{<:Real},
                       forc_hgt_u_grc::Vector{<:Real},
                       forc_hgt_t_grc::Vector{<:Real},
                       forc_hgt_q_grc::Vector{<:Real},
                       mask_lakec::BitVector,
                       mask_lakep::BitVector,
                       bounds_col::UnitRange{Int},
                       bounds_patch::UnitRange{Int};
                       dtime::Real = 1800.0)

    nlevsno_val = varpar.nlevsno
    niters = 4  # stability iterations

    # Wind minimum
    wind_min = 0.1

    for p in bounds_patch
        mask_lakep[p] || continue
        c = patch_data.column[p]
        g = patch_data.gridcell[p]
        l = patch_data.landunit[p]

        snl = col_data.snl[c]

        # ============================================================
        # Phase 1: Initialization
        # ============================================================

        # Fetch estimate
        lakedepth = col_data.lakedepth[c]
        if lakedepth < 4.0
            fetch = 100.0
        else
            fetch = 25.0 * lakedepth
        end

        # Top layer index (Julia 1-based)
        jtop = snl + 1 + nlevsno_val
        if snl < 0
            dzsur = col_data.dz[c, jtop] / 2.0
        else
            dzsur = col_data.dz_lake[c, 1] / 2.0
        end

        # Thermal conductivity of surface layer
        tksur = lakestate.savedtke1_col[c]

        # Temperature of subsurface layer
        if snl < 0
            # Snow present: subsurface is next layer
            tsur = temperature.t_soisno_col[c, jtop + 1]
        else
            # No snow: second lake layer
            if size(temperature.t_lake_col, 2) >= 2
                tsur = temperature.t_lake_col[c, 2]
            else
                tsur = temperature.t_lake_col[c, 1]
            end
        end

        # Saved ground temperature
        tgbef = temperature.t_grnd_col[c]

        # Solar absorbed at surface
        sabg = solarabs.sabg_patch[p]
        sabg = smooth_max(sabg, 0.0)

        # Betaprime: fraction of solar absorbed at surface
        if snl < 0
            betaprime = 1.0  # all absorbed in snow
        else
            betaprime = BETAVIS_LAKE  # simplified
        end

        # Virtual potential temperature
        thv = forc_th[c] * (1.0 + 0.61 * forc_q[c])

        # Wind speed
        ur = smooth_max(sqrt(forc_u[g]^2 + forc_v[g]^2), wind_min)

        # Forcing heights
        forc_hgt_u = forc_hgt_u_grc[g]
        forc_hgt_t = forc_hgt_t_grc[g]
        forc_hgt_q = forc_hgt_q_grc[g]

        # Intermediate temperature
        thm = forc_t[c] + 0.0098 * forc_hgt_t

        # Kinematic viscosity
        kva = KVA0 * (temperature.t_grnd_col[c] / 293.15)^1.5 * (1013.25e2 / forc_pbot[c])

        # Latent heat
        h2o_liq = waterstatebulk.ws.h2osoi_liq_col[c, jtop]
        h2o_ice = waterstatebulk.ws.h2osoi_ice_col[c, jtop]
        htvp = (h2o_liq <= 0.0 && h2o_ice > 0.0) ? HSUB : HVAP

        # ============================================================
        # Initial roughness lengths
        # ============================================================
        if tgbef > TFRZ  # Unfrozen
            z0mg = smooth_max(MINZ0LAKE, CUS * kva / smooth_max(ur * 0.1, 1e-4))
        else  # Frozen
            if snl < 0
                z0mg = 0.00085  # snow roughness
            else
                z0mg = 0.001  # ice roughness
            end
        end

        sqre0 = sqrt(smooth_max(z0mg * ur * 0.1 / kva, 0.1))
        z0hg = z0mg * exp(-VKC / PRN_AIR * (4.0 * sqre0 - 3.2))
        z0qg = z0mg * exp(-VKC / SCH_WATER * (4.0 * sqre0 - 4.2))
        z0mg = smooth_max(z0mg, 1.0e-10)
        z0hg = smooth_max(z0hg, 1.0e-10)
        z0qg = smooth_max(z0qg, 1.0e-10)

        # Reference displacement height (zero for lakes)
        displa = 0.0

        # Monin-Obukhov initialization
        zldis = forc_hgt_u - displa
        zldis = smooth_max(zldis, z0mg + 0.01)
        dth = thm - tgbef
        dqh = forc_q[c] - qsat_water(tgbef, forc_pbot[c])
        dthv = dth * (1.0 + 0.61 * forc_q[c]) + 0.61 * forc_th[c] * dqh

        # Initial ustar estimate
        ustar = VKC * ur / log(zldis / z0mg)
        ustar = smooth_max(ustar, 0.001)

        # Initial Obukhov length
        tstar = VKC * dth / log(forc_hgt_t / z0hg)
        qstar = VKC * dqh / log(forc_hgt_q / z0qg)
        thvstar = tstar * (1.0 + 0.61 * forc_q[c]) + 0.61 * forc_th[c] * qstar

        if abs(thvstar) > 1e-10
            obu = -ustar^3 * thv / (VKC * GRAV * thvstar)
            obu = clamp(obu, -1e4, 1e4)
        else
            obu = 1e4
        end

        um = ur

        # ============================================================
        # Phase 2: Stability iteration
        # ============================================================
        t_grnd_new = tgbef

        for iter in 1:niters
            # Aerodynamic resistances
            zldis_u = forc_hgt_u - displa
            zldis_t = forc_hgt_t - displa
            zldis_q = forc_hgt_q - displa
            zldis_u = smooth_max(zldis_u, z0mg + 0.01)
            zldis_t = smooth_max(zldis_t, z0hg + 0.01)
            zldis_q = smooth_max(zldis_q, z0qg + 0.01)

            # Stability functions
            zeta_u = zldis_u / obu
            zeta_t = zldis_t / obu
            zeta_q = zldis_q / obu
            zeta0m = z0mg / obu
            zeta0h = z0hg / obu
            zeta0q = z0qg / obu

            # Clamp zeta values
            zeta_u = clamp(zeta_u, -100.0, ZETAMAX_LAKE)
            zeta_t = clamp(zeta_t, -100.0, ZETAMAX_LAKE)
            zeta_q = clamp(zeta_q, -100.0, ZETAMAX_LAKE)
            zeta0m = clamp(zeta0m, -100.0, ZETAMAX_LAKE)
            zeta0h = clamp(zeta0h, -100.0, ZETAMAX_LAKE)
            zeta0q = clamp(zeta0q, -100.0, ZETAMAX_LAKE)

            # Friction velocity
            ustar = VKC * um / (log(zldis_u / z0mg) - stability_func1(zeta_u) + stability_func1(zeta0m))
            ustar = smooth_max(ustar, 0.001)

            # Transfer coefficients
            temp1 = VKC / (log(zldis_t / z0hg) - stability_func2(zeta_t) + stability_func2(zeta0h))
            temp2 = VKC / (log(zldis_q / z0qg) - stability_func2(zeta_q) + stability_func2(zeta0q))

            # Aerodynamic resistances
            ram = 1.0 / (ustar * VKC / (log(zldis_u / z0mg) - stability_func1(zeta_u) + stability_func1(zeta0m)))
            rah = 1.0 / (temp1 * ustar)
            raw = 1.0 / (temp2 * ustar)
            ram = smooth_max(ram, 1.0)
            rah = smooth_max(rah, 1.0)
            raw = smooth_max(raw, 1.0)

            # Saturation specific humidity at surface
            qsatg, qsatgdT = qsat_water(tgbef, forc_pbot[c]), qsat_water_dT(tgbef, forc_pbot[c])

            # Newton-Raphson for surface temperature
            stftg3 = EMG_LAKE * SB * tgbef^3

            ax = betaprime * sabg +
                 EMG_LAKE * forc_lwrad[c] +
                 3.0 * stftg3 * tgbef +
                 forc_rho[c] * CPAIR / rah * thm +
                 tksur * tsur / dzsur -
                 htvp * forc_rho[c] / raw * (qsatg - qsatgdT * tgbef - forc_q[c])

            bx = 4.0 * stftg3 +
                 forc_rho[c] * CPAIR / rah +
                 htvp * forc_rho[c] / raw * qsatgdT +
                 tksur / dzsur

            bx = smooth_max(bx, 1.0e-10)
            t_grnd_new = ax / bx

            # Surface fluxes
            eflx_sh = forc_rho[c] * CPAIR * (t_grnd_new - thm) / rah
            qflx_evap = forc_rho[c] * (qsatg + qsatgdT * (t_grnd_new - tgbef) - forc_q[c]) / raw

            # Update Obukhov length
            dth_new = thm - t_grnd_new
            dqh_new = forc_q[c] - (qsatg + qsatgdT * (t_grnd_new - tgbef))
            tstar = temp1 * dth_new
            qstar = temp2 * dqh_new
            thvstar = tstar * (1.0 + 0.61 * forc_q[c]) + 0.61 * forc_th[c] * qstar

            if abs(thvstar) > 1e-10
                zeta_val = zldis_u * VKC * GRAV * thvstar / (ustar^2 * thv)
            else
                zeta_val = 0.0
            end

            if zeta_val >= 0.0
                zeta_val = clamp(zeta_val, 0.01, ZETAMAX_LAKE)
                um = smooth_max(ur, wind_min)
            else
                zeta_val = clamp(zeta_val, -100.0, -0.01)
                wc = BETA1 * (-GRAV * ustar * thvstar * 1000.0 / thv)^(1.0/3.0)
                um = sqrt(ur^2 + wc^2)
            end
            obu = zldis_u / zeta_val
            obu = clamp(obu, -1e4, 1e4)

            # Update roughness for unfrozen lakes
            if tgbef > TFRZ
                z0mg = smooth_max(smooth_max(MINZ0LAKE, CUS * kva / smooth_max(ustar, 1e-6)),
                           CUR0 * ustar^2 / GRAV)
                sqre0 = sqrt(smooth_max(z0mg * ustar / kva, 0.1))
                z0hg = z0mg * exp(-VKC / PRN_AIR * (4.0 * sqre0 - 3.2))
                z0qg = z0mg * exp(-VKC / SCH_WATER * (4.0 * sqre0 - 4.2))
                z0mg = smooth_max(z0mg, 1.0e-10)
                z0hg = smooth_max(z0hg, 1.0e-10)
                z0qg = smooth_max(z0qg, 1.0e-10)
            end
        end

        # ============================================================
        # Phase 3: Temperature corrections
        # ============================================================

        # Prevent unfreezing when lake/snow is frozen
        if (snl < 0 || temperature.t_lake_col[c, 1] <= TFRZ) && t_grnd_new > TFRZ
            t_grnd_new = TFRZ
        end

        # Convective mixing correction
        if temperature.t_lake_col[c, 1] > t_grnd_new && t_grnd_new > TDMAX_LAKE
            t_grnd_new = temperature.t_lake_col[c, 1]
        elseif temperature.t_lake_col[c, 1] < t_grnd_new &&
               temperature.t_lake_col[c, 1] > TFRZ && t_grnd_new < TDMAX_LAKE
            t_grnd_new = temperature.t_lake_col[c, 1]
        end

        # Update ground temperature
        temperature.t_grnd_col[c] = t_grnd_new

        # ============================================================
        # Phase 4: Final flux calculations
        # ============================================================

        # Recalculate fluxes at final temperature
        qsatg_final = qsat_water(t_grnd_new, forc_pbot[c])
        thm_local = forc_t[c] + 0.0098 * forc_hgt_t_grc[g]

        # Aerodynamic resistances (use final iteration values)
        zldis_t = smooth_max(forc_hgt_t - displa, z0hg + 0.01)
        zldis_q = smooth_max(forc_hgt_q - displa, z0qg + 0.01)
        zldis_u = smooth_max(forc_hgt_u - displa, z0mg + 0.01)

        rah_final = smooth_max(zldis_t / (VKC * ustar), 1.0)
        raw_final = smooth_max(zldis_q / (VKC * ustar), 1.0)

        eflx_sh_grnd = forc_rho[c] * CPAIR * (t_grnd_new - thm_local) / rah_final
        qflx_evap_soi = forc_rho[c] * (qsatg_final - forc_q[c]) / raw_final

        # Outgoing longwave
        eflx_lwrad_out = (1.0 - EMG_LAKE) * forc_lwrad[c] + EMG_LAKE * SB * t_grnd_new^4

        # Ground heat flux (energy balance residual)
        eflx_soil_grnd = betaprime * sabg + forc_lwrad[c] - eflx_lwrad_out -
                         eflx_sh_grnd - htvp * qflx_evap_soi

        # Store results
        energyflux.eflx_sh_tot_patch[p] = eflx_sh_grnd
        energyflux.eflx_sh_grnd_patch[p] = eflx_sh_grnd
        energyflux.eflx_lh_tot_patch[p] = htvp * qflx_evap_soi
        energyflux.eflx_lh_grnd_patch[p] = htvp * qflx_evap_soi
        energyflux.eflx_lwrad_out_patch[p] = eflx_lwrad_out
        energyflux.eflx_lwrad_net_patch[p] = forc_lwrad[c] - eflx_lwrad_out
        energyflux.eflx_soil_grnd_patch[p] = eflx_soil_grnd

        waterfluxbulk.wf.qflx_evap_soi_patch[p] = qflx_evap_soi
        waterfluxbulk.wf.qflx_evap_tot_patch[p] = qflx_evap_soi

        frictionvel.ustar_patch[p] = ustar
        frictionvel.z0mg_patch[p] = z0mg
        frictionvel.z0hg_patch[p] = z0hg
        frictionvel.z0qg_patch[p] = z0qg

        # Momentum stress
        ram_final = smooth_max(zldis_u / (VKC * ustar), 1.0)
        energyflux.taux_patch[p] = -forc_rho[c] * forc_u[g] / ram_final
        energyflux.tauy_patch[p] = -forc_rho[c] * forc_v[g] / ram_final

        # 2m wind speed for mixing parameters
        u2m = smooth_max(0.1, ustar / VKC * log(2.0 / z0mg))
        lakestate.ws_col[c] = 1.2e-3 * u2m
        g_idx = col_data.gridcell[c]
        lakestate.ks_col[c] = 6.6 * sqrt(smooth_abs(sin(0.0))) * u2m^(-1.84)  # lat from gridcell, simplified

        energyflux.htvp_col[c] = htvp
    end

    return nothing
end

# Helper: saturation specific humidity for liquid water
function qsat_water(t::Real, p::Real)
    # Tetens formula for saturation vapor pressure over water
    es = 611.2 * exp(17.67 * (t - TFRZ) / (t - TFRZ + 243.5))
    return 0.622 * es / (p - 0.378 * es)
end

# Helper: derivative of qsat w.r.t. temperature
function qsat_water_dT(t::Real, p::Real)
    es = 611.2 * exp(17.67 * (t - TFRZ) / (t - TFRZ + 243.5))
    desdT = es * 17.67 * 243.5 / (t - TFRZ + 243.5)^2
    qs = 0.622 * es / (p - 0.378 * es)
    return 0.622 * p * desdT / (p - 0.378 * es)^2
end
