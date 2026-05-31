# ==========================================================================
# Ported from: src/biogeophys/LakeTemperatureMod.F90 (~1483 lines)
# Calculates surface fluxes and temperature for lakes.
# Created by Zack Subin, 2009
#
# Public subroutines:
#   lake_temperature!      — main driver for lake temperature calculation
#
# Private subroutines:
#   soil_therm_prop_lake!  — thermal conductivities and heat capacities of snow/soil layers
#   phase_change_lake!     — phase change within snow/soil/lake layers
# ==========================================================================

# --- Lake constants (from LakeCon.F90) ---
const BETAVIS      = 0.0       # fraction of visible sunlight absorbed at surface
const ZA_LAKE      = 0.6       # base of surface absorption layer (m)
const N2MIN        = 7.5e-5    # min Brunt-Vaisala freq for enhanced diffusivity (s^-2)
const TDMAX        = 277.0     # temperature of maximum water density (K)
const DEPTHCRIT    = 25.0      # depth beneath which to increase mixing (m)
const MIXFACT      = 10.0      # mixing increase factor for deep lakes
const LAKEPUDDLING = false     # suppress convection when ice present
const LAKE_NO_ED   = false     # suppress enhanced diffusion
const PUDZ         = 0.2       # min total ice thickness for puddling (m)

# =========================================================================
# soil_therm_prop_lake! — Set thermal conductivities and heat capacities
#                         of snow/soil layers for lake columns
# =========================================================================
# ==========================================================================
# GPU kernels for soil_therm_prop_lake! and calculate_total_h2osno!. 2D (column,
# layer) kernels indexed by jj = j + nlevsno (so jj in 1:(nlevsno+nlevgrnd) maps
# Fortran j in -nlevsno+1 : nlevgrnd); per-column kernels for the snow reduction.
# Physical constants are eltype-converted (T(...)/one/zero) so no Float64 reaches
# a Float32-only backend (byte-identical on Float64). Lake soil is saturated
# (satw=1), so the Farouki branch is simpler than soil_temperature's.
# ==========================================================================

# Soil/snow layer thermal conductivity thk (Farouki 1981 soil; Jordan 1991 snow).
@kernel function _lake_soil_tk_kernel!(thk, @Const(mask), @Const(snl), @Const(dz),
        @Const(t_soisno), @Const(h2osoi_liq), @Const(h2osoi_ice),
        @Const(watsat), @Const(tksatu), @Const(tkmg), @Const(tkdry),
        nlevsno::Int, nlevsoi::Int)
    c, jj = @index(Global, NTuple)
    @inbounds if mask[c]
        T = eltype(thk)
        j = jj - nlevsno
        if j >= 1 && j <= nlevsoi
            satw = one(T)
            fl = h2osoi_liq[c, jj] / (h2osoi_ice[c, jj] + h2osoi_liq[c, jj])
            if t_soisno[c, jj] >= T(TFRZ)
                dke = smooth_max(zero(T), log10(satw) + one(T))
                dksat = tksatu[c, j]
            else
                dke = satw
                dksat = tkmg[c, j] * T(0.249)^(fl * watsat[c, j]) * T(2.29)^watsat[c, j]
            end
            thk[c, jj] = dke * dksat + (one(T) - dke) * tkdry[c, j]
            satw = (h2osoi_liq[c, jj] / T(DENH2O) + h2osoi_ice[c, jj] / T(DENICE)) /
                   (dz[c, jj] * watsat[c, j])
            if satw > one(T)
                xicevol = (satw - one(T)) * watsat[c, j]
                thk[c, jj] = (thk[c, jj] + xicevol * T(TKICE)) / (one(T) + xicevol) / (one(T) + xicevol)
            end
        elseif j > nlevsoi
            thk[c, jj] = T(THK_BEDROCK)
        end
        if snl[c] + 1 < 1 && j >= snl[c] + 1 && j <= 0
            bw = (h2osoi_ice[c, jj] + h2osoi_liq[c, jj]) / dz[c, jj]
            thk[c, jj] = T(TKAIR) + (T(7.75e-5) * bw + T(1.105e-6) * bw * bw) * (T(TKICE) - T(TKAIR))
        end
    end
end

# Thermal conductivity at the layer interface + top-soil-layer value.
@kernel function _lake_tk_interface_kernel!(tk, tktopsoillay, @Const(thk), @Const(mask),
        @Const(snl), @Const(z), @Const(zi), nlevsno::Int, nlevgrnd::Int)
    c, jj = @index(Global, NTuple)
    @inbounds if mask[c]
        T = eltype(tk)
        j = jj - nlevsno
        if j >= snl[c] + 1 && j <= nlevgrnd - 1 && j != 0
            jj1 = jj + 1
            tk[c, jj] = thk[c, jj] * thk[c, jj1] * (z[c, jj1] - z[c, jj]) /
                (thk[c, jj] * (z[c, jj1] - zi[c, jj + 1]) + thk[c, jj1] * (zi[c, jj + 1] - z[c, jj]))
        elseif j == 0 && j >= snl[c] + 1
            tk[c, jj] = thk[c, jj]
        elseif j == nlevgrnd
            tk[c, jj] = zero(T)
        end
        if j == 1
            tktopsoillay[c] = thk[c, jj]
        end
    end
end

# Soil heat capacity (de Vries 1963). Indexed by soil layer j in 1:nlevgrnd.
@kernel function _lake_soil_cv_kernel!(cv, @Const(mask), @Const(csol), @Const(watsat),
        @Const(dz), @Const(h2osoi_ice), @Const(h2osoi_liq), nlevsno::Int, nlevsoi::Int)
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        T = eltype(cv)
        jj = j + nlevsno
        cv[c, jj] = csol[c, j] * (one(T) - watsat[c, j]) * dz[c, jj] +
            (h2osoi_ice[c, jj] * T(CPICE) + h2osoi_liq[c, jj] * T(CPLIQ))
        if j > nlevsoi
            cv[c, jj] = T(CSOL_BEDROCK) * dz[c, jj]
        end
    end
end

# Snow layer heat capacity. Indexed by jj in 1:nlevsno (j = jj - nlevsno ≤ 0).
@kernel function _lake_snow_cv_kernel!(cv, @Const(mask), @Const(snl),
        @Const(h2osoi_liq), @Const(h2osoi_ice), nlevsno::Int)
    c, jj = @index(Global, NTuple)
    @inbounds if mask[c]
        T = eltype(cv)
        j = jj - nlevsno
        if snl[c] + 1 < 1 && j >= snl[c] + 1
            cv[c, jj] = T(CPLIQ) * h2osoi_liq[c, jj] + T(CPICE) * h2osoi_ice[c, jj]
        end
    end
end

# Total column snow water = no-layer reservoir + resolved snow ice+liquid.
@kernel function _lake_total_h2osno_kernel!(h2osno_total, @Const(mask), @Const(snl),
        @Const(h2osno_no_layers), @Const(h2osoi_ice), @Const(h2osoi_liq), nlevsno::Int)
    c = @index(Global)
    @inbounds if mask[c]
        tot = h2osno_no_layers[c]
        for jj in 1:nlevsno
            j = jj - nlevsno
            if j >= snl[c] + 1
                tot += h2osoi_ice[c, jj] + h2osoi_liq[c, jj]
            end
        end
        h2osno_total[c] = tot
    end
end

"""
    soil_therm_prop_lake!(col, soilstate, waterstatebulk, temperature,
                          tk, cv, tktopsoillay,
                          mask_lakec, bounds_col)

Calculate thermal conductivities and heat capacities of snow/soil layers
for lake columns.

Ported from `SoilThermProp_Lake` in `LakeTemperatureMod.F90`.
"""
function soil_therm_prop_lake!(col::ColumnData, soilstate::SoilStateData,
                               waterstatebulk::WaterStateBulkData,
                               temperature::TemperatureData,
                               tk::AbstractMatrix{<:Real}, cv::AbstractMatrix{<:Real},
                               tktopsoillay::AbstractVector{<:Real},
                               mask_lakec::AbstractVector{Bool}, bounds_col::UnitRange{Int})
    nlevsno = varpar.nlevsno
    nlevgrnd = varpar.nlevgrnd
    nlevsoi = varpar.nlevsoi

    snl = col.snl
    dz = col.dz
    z = col.z
    zi = col.zi

    watsat = soilstate.watsat_col
    tksatu = soilstate.tksatu_col
    tkmg = soilstate.tkmg_col
    tkdry = soilstate.tkdry_col
    csol = soilstate.csol_col

    h2osoi_liq = waterstatebulk.ws.h2osoi_liq_col
    h2osoi_ice = waterstatebulk.ws.h2osoi_ice_col
    t_soisno = temperature.t_soisno_col

    # Working array for layer thermal conductivity (device-resident via similar()).
    nc = length(mask_lakec)
    nlevtot = nlevsno + varpar.nlevmaxurbgrnd
    FT = eltype(t_soisno)
    thk = fill!(similar(tk, FT, nc, nlevtot), zero(FT))

    # Per-(column, layer) kernels; each launch is ordered (interface tk reads thk).
    _launch!(_lake_soil_tk_kernel!, thk, mask_lakec, snl, dz, t_soisno, h2osoi_liq,
             h2osoi_ice, watsat, tksatu, tkmg, tkdry, nlevsno, nlevsoi;
             ndrange = (nc, nlevsno + nlevgrnd))
    _launch!(_lake_tk_interface_kernel!, tk, tktopsoillay, thk, mask_lakec, snl, z, zi,
             nlevsno, nlevgrnd; ndrange = (nc, nlevsno + nlevgrnd))
    _launch!(_lake_soil_cv_kernel!, cv, mask_lakec, csol, watsat, dz, h2osoi_ice,
             h2osoi_liq, nlevsno, nlevsoi; ndrange = (nc, nlevgrnd))
    _launch!(_lake_snow_cv_kernel!, cv, mask_lakec, snl, h2osoi_liq, h2osoi_ice,
             nlevsno; ndrange = (nc, nlevsno))
    return nothing
end

# =========================================================================
# phase_change_lake! kernels. The 6 scalar loops share the lhabs / qflx_snomelt /
# qflx_snofrz accumulators across loops, so they fuse into TWO per-column kernels
# (lake-water side + snow/soil side) that run the internal j-loops sequentially in
# each thread with local accumulators — no atomics. The water kernel runs first
# (writes the partial accumulators); the soil kernel adds its contribution and
# writes the final eflx_snomelt. Constants eltype-converted; imelt is Int.
# =========================================================================

# Init snow-layer fluxes + snow-without-layers melt + lake-layer phase change.
@kernel function _phase_change_lake_water_kernel!(t_lake, snow_depth, h2osno_no_layers,
        cv_lake, lake_icefrac, qflx_snomelt, qflx_snofrz, qflx_snow_drain, lhabs,
        qflx_snomelt_lyr, qflx_snofrz_lyr, imelt, @Const(mask), @Const(dz_lake),
        nlevsno::Int, nlevlak::Int, dtime)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(t_lake)
        small = T(1.0e-12)
        qsm = zero(T); qsd = zero(T); lh = zero(T)
        for jj in 1:nlevsno
            qflx_snomelt_lyr[c, jj] = zero(T)
            qflx_snofrz_lyr[c, jj] = zero(T)
            imelt[c, jj] = 0
        end
        # snow without snow layers + top lake layer above freezing
        if h2osno_no_layers[c] > zero(T) && t_lake[c, 1] > T(TFRZ)
            heatavail = (t_lake[c, 1] - T(TFRZ)) * cv_lake[c, 1]
            melt = smooth_min(h2osno_no_layers[c], heatavail / T(HFUS))
            heatrem = smooth_max(heatavail - melt * T(HFUS), zero(T))
            t_lake[c, 1] = T(TFRZ) + heatrem / cv_lake[c, 1]
            snow_depth[c] = snow_depth[c] * (one(T) - melt / h2osno_no_layers[c])
            h2osno_no_layers[c] = h2osno_no_layers[c] - melt
            lh += melt * T(HFUS)
            qsm += melt / dtime
            qsd += melt / dtime
            if h2osno_no_layers[c] < small; h2osno_no_layers[c] = zero(T); end
            if snow_depth[c] < small; snow_depth[c] = zero(T); end
        end
        # lake phase change
        for j in 1:nlevlak
            dophase = false
            heatavail = zero(T); melt = zero(T); heatrem = zero(T)
            if t_lake[c, j] > T(TFRZ) && lake_icefrac[c, j] > zero(T)
                dophase = true
                heatavail = (t_lake[c, j] - T(TFRZ)) * cv_lake[c, j]
                melt = smooth_min(lake_icefrac[c, j] * T(DENH2O) * dz_lake[c, j], heatavail / T(HFUS))
                heatrem = smooth_max(heatavail - melt * T(HFUS), zero(T))
            elseif t_lake[c, j] < T(TFRZ) && lake_icefrac[c, j] < one(T)
                dophase = true
                heatavail = (t_lake[c, j] - T(TFRZ)) * cv_lake[c, j]
                melt = smooth_max(-(one(T) - lake_icefrac[c, j]) * T(DENH2O) * dz_lake[c, j], heatavail / T(HFUS))
                heatrem = smooth_min(heatavail - melt * T(HFUS), zero(T))
            end
            if dophase
                lake_icefrac[c, j] = lake_icefrac[c, j] - melt / (T(DENH2O) * dz_lake[c, j])
                lh += melt * T(HFUS)
                cv_lake[c, j] = cv_lake[c, j] + melt * (T(CPLIQ) - T(CPICE))
                t_lake[c, j] = T(TFRZ) + heatrem / cv_lake[c, j]
                if lake_icefrac[c, j] > one(T) - small; lake_icefrac[c, j] = one(T); end
                if lake_icefrac[c, j] < small; lake_icefrac[c, j] = zero(T); end
            end
        end
        qflx_snomelt[c]    = qsm
        qflx_snofrz[c]     = zero(T)
        qflx_snow_drain[c] = qsd
        lhabs[c]           = lh
    end
end

# Snow/soil phase change (adds to lhabs/qflx_snomelt/qflx_snofrz) + final eflx_snomelt.
@kernel function _phase_change_lake_soil_kernel!(t_soisno, h2osoi_ice, h2osoi_liq, cv,
        lhabs, qflx_snomelt, qflx_snofrz, eflx_snomelt, imelt, qflx_snomelt_lyr,
        qflx_snofrz_lyr, @Const(mask), @Const(snl), nlevsno::Int, nlevgrnd::Int, dtime)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(t_soisno)
        small = T(1.0e-12)
        lh = lhabs[c]; qsm = qflx_snomelt[c]; qsf = qflx_snofrz[c]
        for j in (-nlevsno + 1):nlevgrnd
            jj = j + nlevsno
            if j >= snl[c] + 1
                dophase = false
                heatavail = zero(T); melt = zero(T); heatrem = zero(T)
                if t_soisno[c, jj] > T(TFRZ) && h2osoi_ice[c, jj] > zero(T)
                    dophase = true
                    heatavail = (t_soisno[c, jj] - T(TFRZ)) * cv[c, jj]
                    melt = smooth_min(h2osoi_ice[c, jj], heatavail / T(HFUS))
                    heatrem = smooth_max(heatavail - melt * T(HFUS), zero(T))
                    if j <= 0
                        imelt[c, jj] = 1
                        qflx_snomelt_lyr[c, jj] = melt / dtime
                        qsm += qflx_snomelt_lyr[c, jj]
                    end
                elseif t_soisno[c, jj] < T(TFRZ) && h2osoi_liq[c, jj] > zero(T)
                    dophase = true
                    heatavail = (t_soisno[c, jj] - T(TFRZ)) * cv[c, jj]
                    melt = smooth_max(-h2osoi_liq[c, jj], heatavail / T(HFUS))
                    heatrem = smooth_min(heatavail - melt * T(HFUS), zero(T))
                    if j <= 0
                        imelt[c, jj] = 2
                        qflx_snofrz_lyr[c, jj] = -melt / dtime
                        qsf += qflx_snofrz_lyr[c, jj]
                    end
                end
                if dophase
                    h2osoi_ice[c, jj] = h2osoi_ice[c, jj] - melt
                    h2osoi_liq[c, jj] = h2osoi_liq[c, jj] + melt
                    lh += melt * T(HFUS)
                    cv[c, jj] = cv[c, jj] + melt * (T(CPLIQ) - T(CPICE))
                    t_soisno[c, jj] = T(TFRZ) + heatrem / cv[c, jj]
                    if h2osoi_ice[c, jj] < small; h2osoi_ice[c, jj] = zero(T); end
                    if h2osoi_liq[c, jj] < small; h2osoi_liq[c, jj] = zero(T); end
                end
            end
        end
        lhabs[c]        = lh
        qflx_snomelt[c] = qsm
        qflx_snofrz[c]  = qsf
        eflx_snomelt[c] = qsm * T(HFUS)
    end
end

# =========================================================================
# phase_change_lake! — Phase change within snow, soil, and lake layers
# =========================================================================
"""
    phase_change_lake!(col, waterstatebulk, waterdiagbulk, waterfluxbulk,
                       temperature, energyflux, lakestate,
                       cv, cv_lake, lhabs,
                       mask_lakec, bounds_col, dtime)

Calculate the phase change within snow, soil, and lake layers.

Ported from `PhaseChange_Lake` in `LakeTemperatureMod.F90`.
"""
function phase_change_lake!(col::ColumnData,
                            waterstatebulk::WaterStateBulkData,
                            waterdiagbulk::WaterDiagnosticBulkData,
                            waterfluxbulk::WaterFluxBulkData,
                            temperature::TemperatureData,
                            energyflux::EnergyFluxData,
                            lakestate::LakeStateData,
                            cv::AbstractMatrix{<:Real}, cv_lake::AbstractMatrix{<:Real},
                            lhabs::AbstractVector{<:Real},
                            mask_lakec::AbstractVector{Bool}, bounds_col::UnitRange{Int},
                            dtime::Real)
    nlevsno = varpar.nlevsno
    nlevgrnd = varpar.nlevgrnd
    nlevlak = varpar.nlevlak

    dz_lake = col.dz_lake
    snl = col.snl

    snow_depth = waterdiagbulk.snow_depth_col
    h2osno_no_layers = waterstatebulk.ws.h2osno_no_layers_col
    h2osoi_liq = waterstatebulk.ws.h2osoi_liq_col
    h2osoi_ice = waterstatebulk.ws.h2osoi_ice_col

    lake_icefrac = lakestate.lake_icefrac_col

    qflx_snofrz_lyr = waterfluxbulk.wf.qflx_snofrz_lyr_col
    qflx_snomelt_lyr = waterfluxbulk.wf.qflx_snomelt_lyr_col
    qflx_snow_drain = waterfluxbulk.wf.qflx_snow_drain_col
    qflx_snomelt = waterfluxbulk.wf.qflx_snomelt_col
    qflx_snofrz = waterfluxbulk.wf.qflx_snofrz_col

    t_soisno = temperature.t_soisno_col
    t_lake = temperature.t_lake_col
    imelt = temperature.imelt_col

    eflx_snomelt = energyflux.eflx_snomelt_col

    nc = length(mask_lakec)
    dt = convert(eltype(t_lake), dtime)

    # Two fused per-column kernels (water side then snow/soil side); the soil
    # kernel reads the partial lhabs/qflx accumulators the water kernel wrote.
    _launch!(_phase_change_lake_water_kernel!, t_lake, snow_depth, h2osno_no_layers,
        cv_lake, lake_icefrac, qflx_snomelt, qflx_snofrz, qflx_snow_drain, lhabs,
        qflx_snomelt_lyr, qflx_snofrz_lyr, imelt, mask_lakec, dz_lake, nlevsno, nlevlak,
        dt; ndrange = nc)
    _launch!(_phase_change_lake_soil_kernel!, t_soisno, h2osoi_ice, h2osoi_liq, cv,
        lhabs, qflx_snomelt, qflx_snofrz, eflx_snomelt, imelt, qflx_snomelt_lyr,
        qflx_snofrz_lyr, mask_lakec, snl, nlevsno, nlevgrnd, dt; ndrange = nc)
    return nothing
end

# =========================================================================
# lake_temperature! — Main driver for lake temperature calculation
# =========================================================================
# lake_temperature! pre-solve kernels (lt3). Per-column / per-(column,layer) /
# per-patch. LAKEPUDDLING and LAKE_NO_ED are compile-time `const false`, so their
# branches (and the frzn state) constant-fold away exactly as in the scalar loop;
# they are kept referenced for fidelity. Constants + local scalars (km/p0/cwat/
# tkice_eff) eltype-converted so no Float64 reaches a Float32-only backend.
# =========================================================================
@kernel function _lake_init_kernel!(ocvts, ncvts, esum1, esum2, jconvect, jconvectbot,
        lakeresist, bottomconvect, @Const(mask), nlevlak::Int, use_lch4::Bool)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(ocvts)
        ocvts[c] = zero(T); ncvts[c] = zero(T); esum1[c] = zero(T); esum2[c] = zero(T)
        if use_lch4
            jconvect[c] = 0
            jconvectbot[c] = nlevlak + 1
            lakeresist[c] = zero(T)
        end
        bottomconvect[c] = false
    end
end

@kernel function _lake_prior_ice_kernel!(frac_iceold, @Const(mask), @Const(snl),
        @Const(h2osoi_liq), @Const(h2osoi_ice), nlevsno::Int)
    c, jj = @index(Global, NTuple)
    @inbounds if mask[c]
        j = jj - nlevsno
        if j >= snl[c] + 1
            frac_iceold[c, jj] = h2osoi_ice[c, jj] / (h2osoi_liq[c, jj] + h2osoi_ice[c, jj])
        end
    end
end

# Per-patch: fin (= ground net heat flux) + NIR-fraction beta. Lake columns have
# one patch each, so the patch→column writes are 1:1 (no scatter conflict).
@kernel function _lake_solar_nir_kernel!(fin, beta, @Const(maskp), @Const(column),
        @Const(eflx_gnet), @Const(fsds_nir_d), @Const(fsds_nir_i), @Const(fsr_nir_d),
        @Const(fsr_nir_i), @Const(sabg))
    p = @index(Global)
    @inbounds if maskp[p]
        T = eltype(fin)
        c = column[p]
        fin[c] = eflx_gnet[p]
        sabg_nir = fsds_nir_d[p] + fsds_nir_i[p] - fsr_nir_d[p] - fsr_nir_i[p]
        sabg_nir = smooth_min(sabg_nir, sabg[p])
        beta[c] = sabg_nir / smooth_max(T(1.0e-5), sabg[p])
        beta[c] = beta[c] + (one(T) - beta[c]) * T(BETAVIS)
    end
end

@kernel function _lake_density_kernel!(rhow, @Const(mask), @Const(lake_icefrac),
        @Const(t_lake), nlevlak::Int)
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        T = eltype(rhow)
        rhow[c, j] = (one(T) - lake_icefrac[c, j]) *
            T(1000.0) * (one(T) - T(1.9549e-05) * smooth_abs(t_lake[c, j] - T(TDMAX))^T(1.68)) +
            lake_icefrac[c, j] * T(DENICE)
    end
end

# Diffusivity / implied thermal conductivity for interior lake layers j=1:nlevlak-1.
@kernel function _lake_diffusivity_kernel!(kme, tk_lake, frzn, @Const(mask), @Const(rhow),
        @Const(z_lake), @Const(ws), @Const(ks), @Const(t_grnd), @Const(t_lake),
        @Const(snl), @Const(lakedepth), @Const(lake_icefrac),
        nlevlak::Int, km, p0, cwat, tkice_eff)
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        T = eltype(kme)
        drhodz = (rhow[c, j+1] - rhow[c, j]) / (z_lake[c, j+1] - z_lake[c, j])
        n2 = T(GRAV) / rhow[c, j] * drhodz
        num = T(40.0) * n2 * (T(VKC) * z_lake[c, j])^T(2.0)
        den = smooth_max(ws[c]^T(2.0) * exp(-T(2.0) * ks[c] * z_lake[c, j]), T(1.0e-10))
        ri = (-one(T) + sqrt(smooth_max(one(T) + num / den, zero(T)))) / T(20.0)
        if LAKEPUDDLING && j == 1
            frzn[c] = false
        end
        if t_grnd[c] > T(TFRZ) && t_lake[c, 1] > T(TFRZ) && snl[c] == 0 &&
           (!LAKEPUDDLING || (lake_icefrac[c, j] == zero(T) && !frzn[c]))
            ke = T(VKC) * ws[c] * z_lake[c, j] / p0 * exp(-ks[c] * z_lake[c, j]) /
                 (one(T) + T(37.0) * ri * ri)
            kme[c, j] = km + ke
            if !LAKE_NO_ED
                fangkm = T(1.039e-8) * smooth_max(n2, T(N2MIN))^T(-0.43)
                kme[c, j] = kme[c, j] + fangkm
            end
            if lakedepth[c] >= T(DEPTHCRIT)
                kme[c, j] = kme[c, j] * T(MIXFACT)
            end
            tk_lake[c, j] = kme[c, j] * cwat
        else
            kme[c, j] = km
            if !LAKE_NO_ED
                fangkm = T(1.039e-8) * smooth_max(n2, T(N2MIN))^T(-0.43)
                kme[c, j] = kme[c, j] + fangkm
                if lakedepth[c] >= T(DEPTHCRIT)
                    kme[c, j] = kme[c, j] * T(MIXFACT)
                end
                tk_lake[c, j] = kme[c, j] * cwat * tkice_eff /
                    ((one(T) - lake_icefrac[c, j]) * tkice_eff + kme[c, j] * cwat * lake_icefrac[c, j])
            else
                tk_lake[c, j] = T(TKWAT) * tkice_eff /
                    ((one(T) - lake_icefrac[c, j]) * tkice_eff + T(TKWAT) * lake_icefrac[c, j])
            end
            if LAKEPUDDLING
                frzn[c] = true
            end
        end
    end
end

# Bottom lake layer conductivity + savedtke1 + jtop (per column).
@kernel function _lake_diffusivity_bottom_kernel!(kme, tk_lake, savedtke1, jtop,
        @Const(mask), @Const(frzn), @Const(t_grnd), @Const(t_lake), @Const(snl),
        @Const(lake_icefrac), nlevlak::Int, cwat, tkice_eff)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(kme)
        j = nlevlak
        kme[c, nlevlak] = kme[c, nlevlak - 1]
        if t_grnd[c] > T(TFRZ) && t_lake[c, 1] > T(TFRZ) && snl[c] == 0 &&
           (!LAKEPUDDLING || (lake_icefrac[c, j] == zero(T) && !frzn[c]))
            tk_lake[c, j] = tk_lake[c, j - 1]
        else
            if !LAKE_NO_ED
                tk_lake[c, j] = kme[c, j] * cwat * tkice_eff /
                    ((one(T) - lake_icefrac[c, j]) * tkice_eff + kme[c, j] * cwat * lake_icefrac[c, j])
            else
                tk_lake[c, j] = T(TKWAT) * tkice_eff /
                    ((one(T) - lake_icefrac[c, j]) * tkice_eff + T(TKWAT) * lake_icefrac[c, j])
            end
        end
        savedtke1[c] = kme[c, 1] * cwat
        jtop[c] = snl[c] + 1
    end
end

# =========================================================================
# lake_temperature! pre-solve part 2 kernels (lt4).
# =========================================================================
# Solar heat source penetrating the lake (per-patch × lake layer).
@kernel function _lake_heat_source_kernel!(phi, phi_soil, @Const(maskp), @Const(column),
        @Const(etal), @Const(lakedepth), @Const(z_lake), @Const(dz_lake), @Const(t_grnd),
        @Const(t_lake), @Const(snl), @Const(sabg), @Const(beta), @Const(sabg_lyr),
        nlevlak::Int)
    p, j = @index(Global, NTuple)
    @inbounds if maskp[p]
        T = eltype(phi)
        c = column[p]
        eta = etal[c] > zero(T) ? etal[c] : T(1.1925) * smooth_max(lakedepth[c], one(T))^T(-0.424)
        zin  = z_lake[c, j] - T(0.5) * dz_lake[c, j]
        zout = z_lake[c, j] + T(0.5) * dz_lake[c, j]
        rsfin  = exp(-eta * smooth_max(zin - T(ZA_LAKE), zero(T)))
        rsfout = exp(-eta * smooth_max(zout - T(ZA_LAKE), zero(T)))
        phidum = zero(T)
        if t_grnd[c] > T(TFRZ) && t_lake[c, 1] > T(TFRZ) && snl[c] == 0
            phidum = (rsfin - rsfout) * sabg[p] * (one(T) - beta[c])
            if j == nlevlak
                phi_soil[c] = rsfout * sabg[p] * (one(T) - beta[c])
            end
        elseif j == 1 && snl[c] == 0
            phidum = sabg[p] * (one(T) - beta[c])
        elseif j == 1
            phidum = sabg_lyr[p, j]
        else
            phidum = zero(T)
            if j == nlevlak
                phi_soil[c] = zero(T)
            end
        end
        phi[c, j] = phidum
    end
end

# Lake-layer volumetric heat capacity.
@kernel function _lake_cv_lake_kernel!(cv_lake, @Const(mask), @Const(dz_lake),
        @Const(lake_icefrac), nlevlak::Int, cwat, cice_eff)
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        T = eltype(cv_lake)
        cv_lake[c, j] = dz_lake[c, j] * (cwat * (one(T) - lake_icefrac[c, j]) + cice_eff * lake_icefrac[c, j])
    end
end

# Initial lake internal energy sum (per-column, accumulates ocvts) + t_lake_bef.
@kernel function _lake_energy_sum_lake_kernel!(ocvts, t_lake_bef, @Const(mask),
        @Const(cv_lake), @Const(t_lake), @Const(lake_icefrac), @Const(dz_lake),
        nlevlak::Int, cfus)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(ocvts)
        acc = ocvts[c]
        for j in 1:nlevlak
            acc += cv_lake[c, j] * (t_lake[c, j] - T(TFRZ)) + cfus * dz_lake[c, j] * (one(T) - lake_icefrac[c, j])
            t_lake_bef[c, j] = t_lake[c, j]
        end
        ocvts[c] = acc
    end
end

# Initial snow/soil internal energy sum (per-column, accumulates ocvts) + t_soisno_bef.
@kernel function _lake_energy_sum_soil_kernel!(ocvts, t_soisno_bef, @Const(mask),
        @Const(jtop), @Const(cv), @Const(t_soisno), @Const(h2osoi_liq),
        @Const(h2osno_total), nlevsno::Int, nlevgrnd::Int)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(ocvts)
        acc = ocvts[c]
        for j in (-nlevsno + 1):nlevgrnd
            jj = j + nlevsno
            if j >= jtop[c]
                acc += cv[c, jj] * (t_soisno[c, jj] - T(TFRZ)) + T(HFUS) * h2osoi_liq[c, jj]
                if j == 1 && h2osno_total[c] > zero(T) && j == jtop[c]
                    acc -= h2osno_total[c] * T(HFUS)
                end
                t_soisno_bef[c, jj] = t_soisno[c, jj]
            end
        end
        ocvts[c] = acc
    end
end

# Transfer patch sabg / sabg_lyr to column level (per-patch × snow level).
@kernel function _lake_solar_transfer_kernel!(sabg_col, sabg_lyr_col, @Const(maskp),
        @Const(column), @Const(jtop), @Const(sabg), @Const(sabg_lyr), nlevsno::Int)
    p, jidx = @index(Global, NTuple)
    @inbounds if maskp[p]
        c = column[p]
        j = jidx - nlevsno
        if j >= jtop[c]
            if j == jtop[c]
                sabg_col[c] = sabg[p]
            end
            sabg_lyr_col[c, jidx] = sabg_lyr[p, jidx]
        end
    end
end

# Extended-column arrays zx/cvx/phix/tx (per-column; internal ascending j_ext loop
# because the soil branch reads zx[c, nlevlak+nlevsno] written by the lake branch).
@kernel function _lake_extcol_setup_kernel!(zx, cvx, phix, tx, @Const(mask), @Const(jtop),
        @Const(z), @Const(cv), @Const(sabg_lyr_col), @Const(z_lake), @Const(cv_lake),
        @Const(phi), @Const(t_lake), @Const(dz_lake), @Const(phi_soil), @Const(t_soisno),
        nlevsno::Int, nlevlak::Int, nlevgrnd::Int)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(zx)
        extoff = nlevsno
        for j_ext in (-nlevsno + 1):(nlevlak + nlevgrnd)
            ji = j_ext + extoff
            jprime = j_ext - nlevlak
            jj = j_ext + nlevsno
            if j_ext >= jtop[c]
                if j_ext < 1
                    zx[c, ji] = z[c, jj]
                    cvx[c, ji] = cv[c, jj]
                    phix[c, ji] = j_ext == jtop[c] ? zero(T) : sabg_lyr_col[c, j_ext + nlevsno]
                    tx[c, ji] = t_soisno[c, jj]
                elseif j_ext <= nlevlak
                    zx[c, ji] = z_lake[c, j_ext]
                    cvx[c, ji] = cv_lake[c, j_ext]
                    phix[c, ji] = phi[c, j_ext]
                    tx[c, ji] = t_lake[c, j_ext]
                else
                    jprime_jj = jprime + nlevsno
                    zx[c, ji] = zx[c, nlevlak + extoff] + dz_lake[c, nlevlak] / T(2.0) + z[c, jprime_jj]
                    cvx[c, ji] = cv[c, jprime_jj]
                    phix[c, ji] = j_ext == nlevlak + 1 ? phi_soil[c] : zero(T)
                    tx[c, ji] = t_soisno[c, jprime_jj]
                end
            end
        end
    end
end

# =========================================================================
# lake_temperature! tridiagonal-assembly kernels (lt5). All 2D (column, extended
# level ji=j_ext+nlevsno); neighbor reads (zx/tx/tkix/fnx at ji±1) are into arrays
# fully written by earlier kernels, so no intra-kernel race. Constants/dtime
# eltype-converted. The batched solve (tridiagonal_multi!) is already kernelized.
# =========================================================================
@kernel function _lake_tkix_kernel!(tkix, @Const(mask), @Const(jtop), @Const(tk),
        @Const(tk_lake), @Const(z), @Const(z_lake), @Const(dz_lake), @Const(tktopsoillay),
        @Const(zx), nlevsno::Int, nlevlak::Int, nlevgrnd::Int)
    c, ji = @index(Global, NTuple)
    @inbounds if mask[c]
        T = eltype(tkix)
        j_ext = ji - nlevsno
        jprime = j_ext - nlevlak
        jj = j_ext + nlevsno
        if j_ext >= jtop[c]
            if j_ext < 0
                tkix[c, ji] = tk[c, jj]
            elseif j_ext == 0
                dzp = zx[c, ji + 1] - zx[c, ji]
                tkix[c, ji] = tk_lake[c, 1] * tk[c, jj] * dzp /
                    (tk[c, jj] * z_lake[c, 1] + tk_lake[c, 1] * (-z[c, jj]))
            elseif j_ext < nlevlak
                tkix[c, ji] = (tk_lake[c, j_ext] * tk_lake[c, j_ext + 1] * (dz_lake[c, j_ext + 1] + dz_lake[c, j_ext])) /
                    (tk_lake[c, j_ext] * dz_lake[c, j_ext + 1] + tk_lake[c, j_ext + 1] * dz_lake[c, j_ext])
            elseif j_ext == nlevlak
                dzp = zx[c, ji + 1] - zx[c, ji]
                jprime_jj = 1 + nlevsno
                tkix[c, ji] = tktopsoillay[c] * tk_lake[c, j_ext] * dzp /
                    (tktopsoillay[c] * dz_lake[c, j_ext] / T(2.0) + tk_lake[c, j_ext] * z[c, jprime_jj])
            else
                jprime_jj = jprime + nlevsno
                tkix[c, ji] = tk[c, jprime_jj]
            end
        end
    end
end

@kernel function _lake_tridiag_factors_kernel!(factx, fnx, @Const(mask), @Const(jtop),
        @Const(cvx), @Const(tkix), @Const(tx), @Const(zx), nlevsno::Int, nlevlak::Int,
        nlevgrnd::Int, dtime)
    c, ji = @index(Global, NTuple)
    @inbounds if mask[c]
        T = eltype(factx)
        j_ext = ji - nlevsno
        if j_ext >= jtop[c]
            if j_ext < nlevlak + nlevgrnd
                factx[c, ji] = dtime / cvx[c, ji]
                fnx[c, ji] = tkix[c, ji] * (tx[c, ji + 1] - tx[c, ji]) / (zx[c, ji + 1] - zx[c, ji])
            else
                factx[c, ji] = dtime / cvx[c, ji]
                fnx[c, ji] = zero(T)
            end
        end
    end
end

@kernel function _lake_tridiag_coeffs_kernel!(a, b, c1, r, @Const(mask), @Const(jtop),
        @Const(zx), @Const(factx), @Const(tkix), @Const(tx), @Const(fin), @Const(phix),
        @Const(fnx), nlevsno::Int, nlevlak::Int, nlevgrnd::Int, cnfac)
    c, ji = @index(Global, NTuple)
    @inbounds if mask[c]
        T = eltype(a)
        j_ext = ji - nlevsno
        if j_ext >= jtop[c]
            omc = one(T) - cnfac
            if j_ext == jtop[c]
                dzp = zx[c, ji + 1] - zx[c, ji]
                a[c, ji] = zero(T)
                b[c, ji] = one(T) + omc * factx[c, ji] * tkix[c, ji] / dzp
                c1[c, ji] = -omc * factx[c, ji] * tkix[c, ji] / dzp
                r[c, ji] = tx[c, ji] + factx[c, ji] * (fin[c] + phix[c, ji] + cnfac * fnx[c, ji])
            elseif j_ext < nlevlak + nlevgrnd
                dzm = zx[c, ji] - zx[c, ji - 1]
                dzp = zx[c, ji + 1] - zx[c, ji]
                a[c, ji] = -omc * factx[c, ji] * tkix[c, ji - 1] / dzm
                b[c, ji] = one(T) + omc * factx[c, ji] * (tkix[c, ji] / dzp + tkix[c, ji - 1] / dzm)
                c1[c, ji] = -omc * factx[c, ji] * tkix[c, ji] / dzp
                r[c, ji] = tx[c, ji] + cnfac * factx[c, ji] * (fnx[c, ji] - fnx[c, ji - 1]) + factx[c, ji] * phix[c, ji]
            else
                dzm = zx[c, ji] - zx[c, ji - 1]
                a[c, ji] = -omc * factx[c, ji] * tkix[c, ji - 1] / dzm
                b[c, ji] = one(T) + omc * factx[c, ji] * tkix[c, ji - 1] / dzm
                c1[c, ji] = zero(T)
                r[c, ji] = tx[c, ji] - cnfac * factx[c, ji] * fnx[c, ji - 1]
            end
        end
    end
end

@kernel function _lake_unpack_kernel!(t_soisno, t_lake, @Const(mask), @Const(jtop),
        @Const(tx), nlevsno::Int, nlevlak::Int, nlevgrnd::Int)
    c, ji = @index(Global, NTuple)
    @inbounds if mask[c]
        j_ext = ji - nlevsno
        jprime = j_ext - nlevlak
        jj = j_ext + nlevsno
        if j_ext >= jtop[c]
            if j_ext < 1
                t_soisno[c, jj] = tx[c, ji]
            elseif j_ext <= nlevlak
                t_lake[c, j_ext] = tx[c, ji]
            else
                jprime_jj = jprime + nlevsno
                t_soisno[c, jprime_jj] = tx[c, ji]
            end
        end
    end
end

# =========================================================================
"""
    lake_temperature!(col, patch_data, solarabs, soilstate, waterstatebulk,
                      waterdiagbulk, waterfluxbulk, energyflux, temperature,
                      lakestate, grnd_ch4_cond,
                      mask_lakec, mask_lakep, bounds_col, bounds_patch, dtime)

Calculate temperatures in the 25-45 layer column of (possible) snow,
lake water, soil, and bedrock beneath lake.

Uses a Crank-Nicholson tridiagonal system to solve for temperature,
then applies phase change and convective mixing.

Ported from `LakeTemperature` in `LakeTemperatureMod.F90`.
"""
function lake_temperature!(col::ColumnData, patch_data::PatchData,
                           solarabs::SolarAbsorbedData, soilstate::SoilStateData,
                           waterstatebulk::WaterStateBulkData,
                           waterdiagbulk::WaterDiagnosticBulkData,
                           waterfluxbulk::WaterFluxBulkData,
                           energyflux::EnergyFluxData,
                           temperature::TemperatureData,
                           lakestate::LakeStateData,
                           grnd_ch4_cond::Vector{<:Real},
                           mask_lakec::BitVector, mask_lakep::BitVector,
                           bounds_col::UnitRange{Int}, bounds_patch::UnitRange{Int},
                           dtime::Real)
    nlevsno = varpar.nlevsno
    nlevgrnd = varpar.nlevgrnd
    nlevlak = varpar.nlevlak
    joff = nlevsno  # offset: Fortran j → Julia j+joff

    p0 = 1.0  # neutral turbulent Prandtl number

    # Initialize constants
    cwat = CPLIQ * DENH2O       # water heat capacity per unit volume
    cice_eff = CPICE * DENH2O   # effective heat capacity of ice
    cfus = HFUS * DENH2O        # latent heat per unit volume
    tkice_eff = TKICE * DENICE / DENH2O  # effective conductivity
    km = TKWAT / cwat           # molecular diffusivity

    nc = length(bounds_col)

    # Aliases for instance fields
    dz_lake = col.dz_lake
    z_lake = col.z_lake
    dz = col.dz
    z = col.z
    snl = col.snl
    lakedepth = col.lakedepth

    sabg = solarabs.sabg_patch
    sabg_lyr = solarabs.sabg_lyr_patch
    fsds_nir_d = solarabs.fsds_nir_d_patch
    fsds_nir_i = solarabs.fsds_nir_i_patch
    fsr_nir_d = solarabs.fsr_nir_d_patch
    fsr_nir_i = solarabs.fsr_nir_i_patch

    etal = lakestate.etal_col
    ks = lakestate.ks_col
    ws = lakestate.ws_col
    lake_raw = lakestate.lake_raw_col

    h2osoi_liq = waterstatebulk.ws.h2osoi_liq_col
    h2osoi_ice = waterstatebulk.ws.h2osoi_ice_col
    frac_iceold = waterdiagbulk.frac_iceold_col

    t_grnd = temperature.t_grnd_col
    t_soisno = temperature.t_soisno_col
    t_lake = temperature.t_lake_col

    beta = lakestate.betaprime_col
    lake_icefrac = lakestate.lake_icefrac_col
    lake_icefracsurf = lakestate.lake_icefracsurf_col
    lake_icethick = lakestate.lake_icethick_col
    savedtke1 = lakestate.savedtke1_col
    lakeresist = lakestate.lakeresist_col

    eflx_soil_grnd = energyflux.eflx_soil_grnd_patch
    eflx_sh_grnd = energyflux.eflx_sh_grnd_patch
    eflx_sh_tot = energyflux.eflx_sh_tot_patch
    eflx_gnet = energyflux.eflx_gnet_patch
    errsoi = energyflux.errsoi_col

    # Total column levels for tridiagonal: from -nlevsno+1 to nlevlak+nlevgrnd
    ncol_total = nlevsno + nlevlak + nlevgrnd  # total number of levels

    # Local scratch — device-resident via similar()+fill! (referenced off t_lake, a
    # state matrix) so the kernelized loops below run on the state's backend (GPU).
    FT = eltype(temperature.t_grnd_col)
    mF(m::Int)  = fill!(similar(t_lake, FT, nc, m), zero(FT))       # nc×m float matrix
    vF()        = fill!(similar(t_lake, FT, nc), zero(FT))          # nc float vector
    nlevtot = nlevsno + varpar.nlevmaxurbgrnd
    rhow = mF(nlevlak)
    phi = mF(nlevlak)
    kme = mF(nlevlak)
    phi_soil = vF()
    fin = vF()
    ocvts = vF()
    ncvts_arr = vF()
    jtop = fill!(similar(t_lake, Int, nc), 1)
    cv = mF(nlevtot)
    tk = mF(nlevtot)
    cv_lake = mF(nlevlak)
    tk_lake = mF(nlevlak)
    tktopsoillay = vF()
    lhabs = vF()

    # Extended column arrays (from -nlevsno+1 to nlevlak+nlevgrnd, stored 1:ncol_total)
    cvx = mF(ncol_total)
    tkix = mF(ncol_total)
    tx = mF(ncol_total)
    phix = mF(ncol_total)
    zx = mF(ncol_total)
    fnx = mF(ncol_total)
    factx = mF(ncol_total)
    a = mF(ncol_total)
    b = mF(ncol_total)
    c1 = mF(ncol_total)
    r = mF(ncol_total)

    # Other local arrays
    t_lake_bef = mF(nlevlak)
    t_soisno_bef = mF(nlevsno + nlevgrnd)
    sabg_col = vF()
    sabg_lyr_col = mF(nlevsno + 1)  # -nlevsno+1:1
    tav_froz = vF()
    tav_unfr = vF()
    nav = vF()
    iceav = vF()
    qav = vF()
    zsum = vF()
    esum1 = vF()
    esum2 = vF()
    h2osno_total = vF()
    jconvect = fill!(similar(t_lake, Int, nc), 0)
    jconvectbot = fill!(similar(t_lake, Int, nc), nlevlak + 1)
    bottomconvect = fill!(similar(t_lake, Bool, nc), false)
    puddle = fill!(similar(t_lake, Bool, nc), false)
    frzn = fill!(similar(t_lake, Bool, nc), false)
    icesum = vF()

    # Helper: convert extended column index j (Fortran: -nlevsno+1 to nlevlak+nlevgrnd)
    # to Julia 1-based index
    # j_ext → j_ext + nlevsno (since min j_ext = -nlevsno+1, maps to index 1)
    ext_off = nlevsno  # j_ext + ext_off = 1-based index (but -nlevsno+1 + nlevsno = 1)
    # Actually: Fortran index -nlevsno+1 → Julia index 1, so offset = nlevsno
    # But we already have nlevsno = ext_off. Let me be precise:
    # j_ext ranges from (-nlevsno+1) to (nlevlak+nlevgrnd)
    # Julia index = j_ext - (-nlevsno+1) + 1 = j_ext + nlevsno
    # So ext_off = nlevsno (same as joff for snow/soil)
    # BUT for extended column, max j_ext = nlevlak+nlevgrnd
    # Julia index = nlevlak + nlevgrnd + nlevsno = ncol_total ✓

    # 1-3) Init, prior-ice, solar NIR, lake density, diffusivity (per-column /
    # per-(column,layer) / per-patch kernels; launched in order — the diffusivity
    # reads rhow produced by the density kernel; the bottom kernel reads kme).
    use_lch4 = varctl.use_lch4
    _launch!(_lake_init_kernel!, ocvts, ncvts_arr, esum1, esum2, jconvect, jconvectbot,
        lakeresist, bottomconvect, mask_lakec, nlevlak, use_lch4; ndrange = nc)
    _launch!(_lake_prior_ice_kernel!, frac_iceold, mask_lakec, snl, h2osoi_liq, h2osoi_ice,
        nlevsno; ndrange = (nc, nlevsno))
    _launch!(_lake_solar_nir_kernel!, fin, beta, mask_lakep, patch_data.column, eflx_gnet,
        fsds_nir_d, fsds_nir_i, fsr_nir_d, fsr_nir_i, sabg; ndrange = length(mask_lakep))
    _launch!(_lake_density_kernel!, rhow, mask_lakec, lake_icefrac, t_lake, nlevlak;
        ndrange = (nc, nlevlak))
    _launch!(_lake_diffusivity_kernel!, kme, tk_lake, frzn, mask_lakec, rhow, z_lake, ws,
        ks, t_grnd, t_lake, snl, lakedepth, lake_icefrac, nlevlak,
        FT(km), FT(p0), FT(cwat), FT(tkice_eff); ndrange = (nc, nlevlak - 1))
    _launch!(_lake_diffusivity_bottom_kernel!, kme, tk_lake, savedtke1, jtop, mask_lakec,
        frzn, t_grnd, t_lake, snl, lake_icefrac, nlevlak, FT(cwat), FT(tkice_eff);
        ndrange = nc)

    # 4!) Solar heat source penetrating the lake (per-patch × lake layer).
    _launch!(_lake_heat_source_kernel!, phi, phi_soil, mask_lakep, patch_data.column,
        etal, lakedepth, z_lake, dz_lake, t_grnd, t_lake, snl, sabg, beta, sabg_lyr,
        nlevlak; ndrange = (length(mask_lakep), nlevlak))

    # 5!) Thermal properties + initial energy content.
    _launch!(_lake_cv_lake_kernel!, cv_lake, mask_lakec, dz_lake, lake_icefrac, nlevlak,
        FT(cwat), FT(cice_eff); ndrange = (nc, nlevlak))
    soil_therm_prop_lake!(col, soilstate, waterstatebulk, temperature,
                          tk, cv, tktopsoillay, mask_lakec, bounds_col)
    _launch!(_lake_energy_sum_lake_kernel!, ocvts, t_lake_bef, mask_lakec, cv_lake, t_lake,
        lake_icefrac, dz_lake, nlevlak, FT(cfus); ndrange = nc)
    calculate_total_h2osno!(waterstatebulk.ws, snl, h2osno_total, mask_lakec, bounds_col)
    _launch!(_lake_energy_sum_soil_kernel!, ocvts, t_soisno_bef, mask_lakec, jtop, cv,
        t_soisno, h2osoi_liq, h2osno_total, nlevsno, nlevgrnd; ndrange = nc)

    # 6!) Transfer solar to column level + build extended-column arrays.
    _launch!(_lake_solar_transfer_kernel!, sabg_col, sabg_lyr_col, mask_lakep,
        patch_data.column, jtop, sabg, sabg_lyr, nlevsno;
        ndrange = (length(mask_lakep), nlevsno + 1))
    _launch!(_lake_extcol_setup_kernel!, zx, cvx, phix, tx, mask_lakec, jtop, z, cv,
        sabg_lyr_col, z_lake, cv_lake, phi, t_lake, dz_lake, phi_soil, t_soisno,
        nlevsno, nlevlak, nlevgrnd; ndrange = nc)

    # Interface conductivity + Crank-Nicolson tridiagonal assembly (2D kernels).
    _launch!(_lake_tkix_kernel!, tkix, mask_lakec, jtop, tk, tk_lake, z, z_lake, dz_lake,
        tktopsoillay, zx, nlevsno, nlevlak, nlevgrnd; ndrange = (nc, ncol_total))
    _launch!(_lake_tridiag_factors_kernel!, factx, fnx, mask_lakec, jtop, cvx, tkix, tx, zx,
        nlevsno, nlevlak, nlevgrnd, FT(dtime); ndrange = (nc, ncol_total))
    _launch!(_lake_tridiag_coeffs_kernel!, a, b, c1, r, mask_lakec, jtop, zx, factx, tkix,
        tx, fin, phix, fnx, nlevsno, nlevlak, nlevgrnd, FT(CNFAC); ndrange = (nc, ncol_total))

    # 7!) Solve the batched tridiagonal system (already a kernel).
    nlevs_total = ncol_total
    jtop_ext = jtop .+ ext_off   # extended-column jtop (device-friendly broadcast)
    tridiagonal_multi!(tx, a, b, c1, r, jtop_ext, mask_lakec, nc, nlevs_total)

    # Unpack the solution back into t_soisno / t_lake.
    _launch!(_lake_unpack_kernel!, t_soisno, t_lake, mask_lakec, jtop, tx,
        nlevsno, nlevlak, nlevgrnd; ndrange = (nc, ncol_total))

    # 9!) Phase change
    phase_change_lake!(col, waterstatebulk, waterdiagbulk, waterfluxbulk,
                       temperature, energyflux, lakestate,
                       cv, cv_lake, lhabs, mask_lakec, bounds_col, dtime)

    # 10!) Convective mixing
    # Recalculate density
    for j in 1:nlevlak
        for c in bounds_col
            mask_lakec[c] || continue
            rhow[c, j] = (1.0 - lake_icefrac[c, j]) *
                1000.0 * (1.0 - 1.9549e-05 * smooth_abs(t_lake[c, j] - TDMAX)^1.68) +
                lake_icefrac[c, j] * DENICE
        end
    end

    if LAKEPUDDLING
        for j in 1:nlevlak
            for c in bounds_col
                mask_lakec[c] || continue
                if j == 1
                    icesum[c] = 0.0
                    puddle[c] = false
                end
                icesum[c] = icesum[c] + lake_icefrac[c, j] * dz[c, j + joff]
                if j == nlevlak
                    if icesum[c] >= PUDZ
                        puddle[c] = true
                    end
                end
            end
        end
    end

    # Top nlevlak-2 layers convection
    for j in 1:(nlevlak - 2)
        for c in bounds_col
            mask_lakec[c] || continue
            qav[c] = 0.0
            nav[c] = 0.0
            iceav[c] = 0.0
        end

        for i in 1:(j + 1)
            for c in bounds_col
                mask_lakec[c] || continue
                if (!LAKEPUDDLING || !puddle[c]) && (rhow[c, j] > rhow[c, j+1] ||
                    (lake_icefrac[c, j] < 1.0 && lake_icefrac[c, j+1] > 0.0))
                    qav[c] = qav[c] + dz_lake[c, i] * (t_lake[c, i] - TFRZ) *
                        ((1.0 - lake_icefrac[c, i]) * cwat + lake_icefrac[c, i] * cice_eff)
                    iceav[c] = iceav[c] + lake_icefrac[c, i] * dz_lake[c, i]
                    nav[c] = nav[c] + dz_lake[c, i]
                    if varctl.use_lch4
                        jconvect[c] = j + 1
                    end
                end
            end
        end

        for c in bounds_col
            mask_lakec[c] || continue
            if (!LAKEPUDDLING || !puddle[c]) && (rhow[c, j] > rhow[c, j+1] ||
                (lake_icefrac[c, j] < 1.0 && lake_icefrac[c, j+1] > 0.0))
                qav[c] = qav[c] / nav[c]
                iceav[c] = iceav[c] / nav[c]
                if qav[c] > 0.0
                    tav_froz[c] = 0.0
                    tav_unfr[c] = qav[c] / ((1.0 - iceav[c]) * cwat)
                elseif qav[c] < 0.0
                    tav_froz[c] = qav[c] / (iceav[c] * cice_eff)
                    tav_unfr[c] = 0.0
                else
                    tav_froz[c] = 0.0
                    tav_unfr[c] = 0.0
                end
            end
        end

        for i in 1:(j + 1)
            for c in bounds_col
                mask_lakec[c] || continue
                if nav[c] > 0.0
                    if i == 1
                        zsum[c] = 0.0
                    end
                    if (zsum[c] + dz_lake[c, i]) / nav[c] <= iceav[c]
                        lake_icefrac[c, i] = 1.0
                        t_lake[c, i] = tav_froz[c] + TFRZ
                    elseif zsum[c] / nav[c] < iceav[c]
                        lake_icefrac[c, i] = (iceav[c] * nav[c] - zsum[c]) / dz_lake[c, i]
                        t_lake[c, i] = (lake_icefrac[c, i] * tav_froz[c] * cice_eff +
                            (1.0 - lake_icefrac[c, i]) * tav_unfr[c] * cwat) /
                            (lake_icefrac[c, i] * cice_eff + (1 - lake_icefrac[c, i]) * cwat) + TFRZ
                    else
                        lake_icefrac[c, i] = 0.0
                        t_lake[c, i] = tav_unfr[c] + TFRZ
                    end
                    zsum[c] = zsum[c] + dz_lake[c, i]

                    rhow[c, i] = (1.0 - lake_icefrac[c, i]) *
                        1000.0 * (1.0 - 1.9549e-05 * smooth_abs(t_lake[c, i] - TDMAX)^1.68) +
                        lake_icefrac[c, i] * DENICE
                end
            end
        end
    end

    # Check bottom layer for convection
    j = nlevlak - 1
    for c in bounds_col
        mask_lakec[c] || continue
        if (!LAKEPUDDLING || !puddle[c]) && (rhow[c, j] > rhow[c, j+1] ||
            (lake_icefrac[c, j] < 1.0 && lake_icefrac[c, j+1] > 0.0))
            bottomconvect[c] = true
        end
    end

    # Bottom-up convective mixing
    for j in (nlevlak - 1):-1:1
        for c in bounds_col
            mask_lakec[c] || continue
            qav[c] = 0.0
            nav[c] = 0.0
            iceav[c] = 0.0
        end

        for i in j:nlevlak
            for c in bounds_col
                mask_lakec[c] || continue
                if bottomconvect[c] &&
                   (!LAKEPUDDLING || !puddle[c]) && (rhow[c, j] > rhow[c, j+1] ||
                    (lake_icefrac[c, j] < 1.0 && lake_icefrac[c, j+1] > 0.0))
                    qav[c] = qav[c] + dz_lake[c, i] * (t_lake[c, i] - TFRZ) *
                        ((1.0 - lake_icefrac[c, i]) * cwat + lake_icefrac[c, i] * cice_eff)
                    iceav[c] = iceav[c] + lake_icefrac[c, i] * dz_lake[c, i]
                    nav[c] = nav[c] + dz_lake[c, i]
                    if varctl.use_lch4
                        jconvectbot[c] = j
                    end
                end
            end
        end

        for c in bounds_col
            mask_lakec[c] || continue
            if bottomconvect[c] &&
               (!LAKEPUDDLING || !puddle[c]) && (rhow[c, j] > rhow[c, j+1] ||
                (lake_icefrac[c, j] < 1.0 && lake_icefrac[c, j+1] > 0.0))
                qav[c] = qav[c] / nav[c]
                iceav[c] = iceav[c] / nav[c]
                if qav[c] > 0.0
                    tav_froz[c] = 0.0
                    tav_unfr[c] = qav[c] / ((1.0 - iceav[c]) * cwat)
                elseif qav[c] < 0.0
                    tav_froz[c] = qav[c] / (iceav[c] * cice_eff)
                    tav_unfr[c] = 0.0
                else
                    tav_froz[c] = 0.0
                    tav_unfr[c] = 0.0
                end
            end
        end

        for i in j:nlevlak
            for c in bounds_col
                mask_lakec[c] || continue
                if bottomconvect[c] && nav[c] > 0.0
                    if i == j
                        zsum[c] = 0.0
                    end
                    if (zsum[c] + dz_lake[c, i]) / nav[c] <= iceav[c]
                        lake_icefrac[c, i] = 1.0
                        t_lake[c, i] = tav_froz[c] + TFRZ
                    elseif zsum[c] / nav[c] < iceav[c]
                        lake_icefrac[c, i] = (iceav[c] * nav[c] - zsum[c]) / dz_lake[c, i]
                        t_lake[c, i] = (lake_icefrac[c, i] * tav_froz[c] * cice_eff +
                            (1.0 - lake_icefrac[c, i]) * tav_unfr[c] * cwat) /
                            (lake_icefrac[c, i] * cice_eff + (1 - lake_icefrac[c, i]) * cwat) + TFRZ
                    else
                        lake_icefrac[c, i] = 0.0
                        t_lake[c, i] = tav_unfr[c] + TFRZ
                    end
                    zsum[c] = zsum[c] + dz_lake[c, i]

                    rhow[c, i] = (1.0 - lake_icefrac[c, i]) *
                        1000.0 * (1.0 - 1.9549e-05 * smooth_abs(t_lake[c, i] - TDMAX)^1.68) +
                        lake_icefrac[c, i] * DENICE
                end
            end
        end
    end

    # Calculate lakeresist and grnd_ch4_cond for CH4 module
    if varctl.use_lch4
        for j in 1:nlevlak
            for c in bounds_col
                mask_lakec[c] || continue
                if j > jconvect[c] && j < jconvectbot[c]
                    lakeresist[c] = lakeresist[c] + dz_lake[c, j] / kme[c, j]
                end
                if j == nlevlak
                    grnd_ch4_cond[c] = 1.0 / (lakeresist[c] + lake_raw[c])
                    if lake_icefrac[c, 1] > 0.1
                        grnd_ch4_cond[c] = 0.0
                    end
                end
            end
        end
    end

    # 11!) Re-evaluate thermal properties and sum energy content
    # Lake thermal properties
    for j in 1:nlevlak
        for c in bounds_col
            mask_lakec[c] || continue
            cv_lake[c, j] = dz_lake[c, j] * (cwat * (1.0 - lake_icefrac[c, j]) + cice_eff * lake_icefrac[c, j])
        end
    end

    # Snow/soil thermal properties
    soil_therm_prop_lake!(col, soilstate, waterstatebulk, temperature,
                          tk, cv, tktopsoillay, mask_lakec, bounds_col)

    # Sum energy content
    for j in 1:nlevlak
        for c in bounds_col
            mask_lakec[c] || continue
            ncvts_arr[c] = ncvts_arr[c] + cv_lake[c, j] * (t_lake[c, j] - TFRZ) +
                cfus * dz_lake[c, j] * (1.0 - lake_icefrac[c, j])
            fin[c] = fin[c] + phi[c, j]
        end
    end

    calculate_total_h2osno!(waterstatebulk.ws, snl, h2osno_total, mask_lakec, bounds_col)

    for j in (-nlevsno + 1):nlevgrnd
        jj = j + joff
        for c in bounds_col
            mask_lakec[c] || continue
            if j >= jtop[c]
                ncvts_arr[c] = ncvts_arr[c] + cv[c, jj] * (t_soisno[c, jj] - TFRZ) +
                    HFUS * h2osoi_liq[c, jj]
                if j < 1
                    fin[c] = fin[c] + phix[c, j + ext_off]
                end
                if j == 1 && h2osno_total[c] > 0.0 && j == jtop[c]
                    ncvts_arr[c] = ncvts_arr[c] - h2osno_total[c] * HFUS
                end
            end
            if j == 1
                fin[c] = fin[c] + phi_soil[c]
            end
        end
    end

    # Check energy conservation
    for p in bounds_patch
        mask_lakep[p] || continue
        c = patch_data.column[p]
        errsoi[c] = (ncvts_arr[c] - ocvts[c]) / dtime - fin[c]
        if abs(errsoi[c]) < 0.10
            eflx_sh_tot[p] = eflx_sh_tot[p] - errsoi[c]
            eflx_sh_grnd[p] = eflx_sh_grnd[p] - errsoi[c]
            eflx_soil_grnd[p] = eflx_soil_grnd[p] + errsoi[c]
            eflx_gnet[p] = eflx_gnet[p] + errsoi[c]
            if abs(errsoi[c]) > 1.0e-3
                @warn "errsoi incorporated into sensible heat in lake_temperature!: c=$c, errsoi=$(errsoi[c]) W/m^2"
            end
            errsoi[c] = 0.0
        end
    end

    # Lake ice thickness and surface ice fraction diagnostic
    for j in 1:nlevlak
        for c in bounds_col
            mask_lakec[c] || continue
            if j == 1
                lake_icethick[c] = 0.0
                lake_icefracsurf[c] = lake_icefrac[c, 1]
            end
            lake_icethick[c] = lake_icethick[c] + lake_icefrac[c, j] * dz_lake[c, j] * DENH2O / DENICE
        end
    end

    return nothing
end

# =========================================================================
# Helper: calculate total h2osno (simplified version for lake use)
# =========================================================================
"""
    calculate_total_h2osno!(ws, snl, h2osno_total, mask, bounds)

Calculate total snow water (h2osno_no_layers + snow layer ice + liq).
Simplified version matching `CalculateTotalH2osno` for lake columns.
"""
function calculate_total_h2osno!(ws::WaterStateData, snl::AbstractVector{<:Integer},
                                 h2osno_total::AbstractVector{<:Real},
                                 mask::AbstractVector{Bool}, bounds::UnitRange{Int})
    nlevsno = varpar.nlevsno
    _launch!(_lake_total_h2osno_kernel!, h2osno_total, mask, snl,
             ws.h2osno_no_layers_col, ws.h2osoi_ice_col, ws.h2osoi_liq_col, nlevsno;
             ndrange = length(mask))
    return nothing
end
