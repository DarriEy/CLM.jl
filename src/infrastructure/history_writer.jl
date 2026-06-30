# ==========================================================================
# Minimal history output writer — NetCDF time series
# No direct Fortran equivalent (Fortran histFileMod is very complex)
#
# Aggregates patch- and column-level fields to gridcell level using subgrid
# area weights (patch%wtgcell, col%wtgcell), matching Fortran CLM h0 output.
#
# Public functions:
#   history_writer_init!   — Create output file with dimensions
#   history_write_step!    — Write one timestep of output fields
#   history_writer_close!  — Close output file
#   default_hist_fields    — Return list of default output fields
# ==========================================================================

"""
    HistFieldDef

Definition of a single history output field.
"""
struct HistFieldDef
    name::String
    long_name::String
    units::String
    level::String  # "gridcell", "column", or "patch"
    getter::Function  # (inst::CLMInstances) -> AbstractArray
    reduce::Symbol    # :avg (daily mean, default) or :snapshot (end-of-day value)
end

# Convenience constructor: fields default to daily-mean reduction.
HistFieldDef(name, long_name, units, level, getter) =
    HistFieldDef(name, long_name, units, level, getter, :avg)

"""
    HistoryWriter

Stateful NetCDF history output writer. All fields are aggregated to gridcell
level using subgrid area weights before writing.
"""
Base.@kwdef mutable struct HistoryWriter
    ds::Union{NCDataset, Nothing} = nothing
    fields::Vector{HistFieldDef} = HistFieldDef[]
    time_index::Int = 0
    filepath::String = ""
    # Daily averaging accumulators
    accum::Dict{String, Vector{Float64}} = Dict{String, Vector{Float64}}()
    accum_count::Int = 0
end

"""
    default_hist_fields() -> Vector{HistFieldDef}

Return default set of history output fields.
"""
function default_hist_fields()
    return HistFieldDef[
        HistFieldDef("T_GRND", "ground temperature", "K", "column",
            inst -> inst.temperature.t_grnd_col),
        HistFieldDef("QFLX_EVAP_TOT", "total evaporation", "mm/s", "column",
            inst -> inst.water.waterfluxbulk_inst.wf.qflx_evap_tot_col),
        HistFieldDef("EFLX_LH_TOT", "total latent heat flux", "W/m2", "patch",
            inst -> inst.energyflux.eflx_lh_tot_patch),
        HistFieldDef("EFLX_SH_TOT", "total sensible heat flux", "W/m2", "patch",
            inst -> inst.energyflux.eflx_sh_tot_patch),
        HistFieldDef("FSA", "absorbed solar radiation", "W/m2", "patch",
            inst -> inst.solarabs.fsa_patch),
        HistFieldDef("H2OSNO", "snow water equivalent", "mm", "column",
            history_h2osno_total_col),
        HistFieldDef("TSA", "2m air temperature", "K", "patch",
            inst -> inst.temperature.t_ref2m_patch),
        HistFieldDef("RAIN", "atmospheric rain", "mm/s", "column",
            inst -> length(inst.atm2lnd.forc_rain_downscaled_col) > 0 ?
                inst.atm2lnd.forc_rain_downscaled_col :
                zeros(length(inst.column.gridcell))),
        HistFieldDef("SNOW", "atmospheric snow", "mm/s", "column",
            inst -> length(inst.atm2lnd.forc_snow_downscaled_col) > 0 ?
                inst.atm2lnd.forc_snow_downscaled_col :
                zeros(length(inst.column.gridcell))),
        HistFieldDef("QRUNOFF", "total runoff", "mm/s", "column",
            inst -> inst.water.waterfluxbulk_inst.wf.qflx_runoff_col),
        HistFieldDef("ZWT", "water table depth", "m", "column",
            history_zwt_veg_col),
        HistFieldDef("TLAI", "total projected leaf area index", "m2/m2", "patch",
            inst -> inst.canopystate.tlai_patch),
        HistFieldDef("ELAI", "exposed one-sided leaf area index", "m2/m2", "patch",
            inst -> inst.canopystate.elai_patch),
        HistFieldDef("QFLX_DRAIN", "sub-surface drainage", "mm/s", "column",
            inst -> inst.water.waterfluxbulk_inst.wf.qflx_drain_col),
        HistFieldDef("QFLX_INFL", "infiltration", "mm/s", "column",
            inst -> inst.water.waterfluxbulk_inst.wf.qflx_infl_col),
        HistFieldDef("QFLX_TRAN_VEG", "vegetation transpiration", "mm/s", "column",
            inst -> inst.water.waterfluxbulk_inst.wf.qflx_tran_veg_col),
        HistFieldDef("SNOW_DEPTH", "snow depth", "m", "column",
            inst -> inst.water.waterdiagnosticbulk_inst.snow_depth_col),
        HistFieldDef("SNOWDP", "area-averaged snow depth", "m", "column",
            inst -> inst.water.waterdiagnosticbulk_inst.snowdp_col),
        HistFieldDef("FRAC_SNO", "fraction of ground covered by snow", "unitless", "column",
            inst -> inst.water.waterdiagnosticbulk_inst.frac_sno_col),
        # BTRANMN is a daily-MIN, not a daily-mean: snapshot the completed daily-min
        # (btran_min_patch) at the end-of-day flush. energyflux_update_acc_vars! runs
        # before the history flush, so at that moment btran_min_patch holds THIS day's
        # min — writing it as a snapshot puts today's min in today's record (matching
        # the averaged fields), instead of the day-average of a field that holds
        # yesterday's completed min for most of the day (which lagged BTRANMN by 1 day).
        HistFieldDef("BTRAN", "transpiration beta factor", "unitless", "patch",
            history_btran_daily_min_patch, :snapshot),
        # PHS diagnostics (late-season over-transpiration investigation)
        HistFieldDef("VEGWPSUN", "sunlit leaf water potential", "mm", "patch",
            history_vegwpsun_patch),
        HistFieldDef("VEGWPROOT", "root water potential", "mm", "patch",
            history_vegwproot_patch),
        HistFieldDef("KSOILROOT", "soil-to-root conductance (layer sum)", "mm/s", "patch",
            history_ksoilroot_patch),
        HistFieldDef("SMPLMEAN", "soil matric potential (soil-layer mean)", "mm", "column",
            history_smplmean_col),
        HistFieldDef("SOILICE", "soil column ice", "mm", "column",
            history_soilice_col),
        HistFieldDef("BSUNDIAG", "sunlit PHS water-stress factor", "unitless", "patch",
            inst -> isempty(inst.energyflux.bsun_patch) ? Float64[] :
                    Float64.(inst.energyflux.bsun_patch)),
        HistFieldDef("SOILWATER_10CM", "soil water+ice in top 0.1 m", "mm", "column",
            history_soilwater_10cm_col),
        HistFieldDef("SOILRESIS", "soil evaporative resistance (S&L14)", "s/m", "column",
            inst -> isempty(inst.soilstate.soilresis_col) ? Float64[] :
                    Float64.(inst.soilstate.soilresis_col)),
        HistFieldDef("DSL", "dry surface layer thickness", "mm", "column",
            inst -> isempty(inst.soilstate.dsl_col) ? Float64[] :
                    Float64.(inst.soilstate.dsl_col)),
        HistFieldDef("FSR", "reflected solar radiation", "W/m2", "patch",
            inst -> isempty(inst.solarabs.fsr_patch) ? Float64[] :
                    Float64.(inst.solarabs.fsr_patch)),
        HistFieldDef("DHSDTCANOPY", "canopy heat storage rate", "W/m2", "patch",
            inst -> isempty(inst.energyflux.dhsdt_canopy_patch) ? Float64[] :
                    Float64.(inst.energyflux.dhsdt_canopy_patch)),
        HistFieldDef("FSRVD", "reflected direct beam vis solar", "W/m2", "patch",
            inst -> isempty(inst.surfrad.fsr_vis_d_patch) ? Float64[] :
                    Float64.(inst.surfrad.fsr_vis_d_patch)),
        HistFieldDef("FSRND", "reflected direct beam NIR solar", "W/m2", "patch",
            inst -> isempty(inst.solarabs.fsr_nir_d_patch) ? Float64[] :
                    Float64.(inst.solarabs.fsr_nir_d_patch)),
        HistFieldDef("FSRVI", "reflected diffuse vis solar", "W/m2", "patch",
            inst -> isempty(inst.surfrad.fsr_vis_i_patch) ? Float64[] :
                    Float64.(inst.surfrad.fsr_vis_i_patch)),
        HistFieldDef("FSRNI", "reflected diffuse NIR solar", "W/m2", "patch",
            inst -> isempty(inst.solarabs.fsr_nir_i_patch) ? Float64[] :
                    Float64.(inst.solarabs.fsr_nir_i_patch)),
        HistFieldDef("FSDSVD", "incident direct beam vis solar", "W/m2", "patch",
            inst -> isempty(inst.surfrad.fsds_vis_d_patch) ? Float64[] :
                    Float64.(inst.surfrad.fsds_vis_d_patch)),
        HistFieldDef("FSDSND", "incident direct beam NIR solar", "W/m2", "patch",
            inst -> isempty(inst.solarabs.fsds_nir_d_patch) ? Float64[] :
                    Float64.(inst.solarabs.fsds_nir_d_patch)),
        HistFieldDef("FSDSVI", "incident diffuse vis solar", "W/m2", "patch",
            inst -> isempty(inst.surfrad.fsds_vis_i_patch) ? Float64[] :
                    Float64.(inst.surfrad.fsds_vis_i_patch)),
        HistFieldDef("FSDSNI", "incident diffuse NIR solar", "W/m2", "patch",
            inst -> isempty(inst.solarabs.fsds_nir_i_patch) ? Float64[] :
                    Float64.(inst.solarabs.fsds_nir_i_patch)),
        HistFieldDef("FIRA", "net infrared (longwave) radiation", "W/m2", "patch",
            inst -> isempty(inst.energyflux.eflx_lwrad_net_patch) ? Float64[] :
                    Float64.(inst.energyflux.eflx_lwrad_net_patch)),
        HistFieldDef("FGR", "soil heat flux into ground", "W/m2", "patch",
            inst -> isempty(inst.energyflux.eflx_soil_grnd_patch) ? Float64[] :
                    Float64.(inst.energyflux.eflx_soil_grnd_patch)),
        HistFieldDef("TV", "vegetation temperature", "K", "patch",
            inst -> isempty(inst.temperature.t_veg_patch) ? Float64[] :
                    Float64.(inst.temperature.t_veg_patch)),
        HistFieldDef("TSOI_10CM", "soil temperature in top 0.1 m", "K", "column",
            history_tsoi_10cm_col),
        HistFieldDef("FPSN", "gross photosynthesis", "umol/m2/s", "patch",
            history_fpsn_patch),
        HistFieldDef("LAISUN", "sunlit projected LAI", "m2/m2", "patch",
            inst -> isempty(inst.canopystate.laisun_patch) ? Float64[] :
                    Float64.(inst.canopystate.laisun_patch)),
        HistFieldDef("LAISHA", "shaded projected LAI", "m2/m2", "patch",
            inst -> isempty(inst.canopystate.laisha_patch) ? Float64[] :
                    Float64.(inst.canopystate.laisha_patch)),
        HistFieldDef("FSUNDIAG", "sunlit fraction of canopy", "1", "patch",
            inst -> isempty(inst.canopystate.fsun_patch) ? Float64[] :
                    [(@inbounds v=inst.canopystate.fsun_patch[p]; isfinite(v) ? Float64(v) : 0.0)
                     for p in 1:length(inst.canopystate.fsun_patch)]),
        HistFieldDef("PARVEG", "canopy absorbed PAR (sun+sha)", "W/m2", "patch",
            history_parveg_patch),
        # NOTE: this is the LUNA-acclimated capacity (vcmx25_z), which drives
        # photosynthesis. It is NOT comparable to the Fortran h0 'VCMX25T'
        # (luvcmax25top), which on this run is the vestigial vcmax_opt=3 diagnostic
        # (i_vcad+s_vcad*lnc)*dayl_factor, not the LUNA value. Named VCMX25LUNA to
        # avoid the false comparison. Bare patches → NaN so the gridcell aggregate
        # is veg-area-weighted (per-patch already validated vs Fortran: tree 55 vs
        # 52, grass 128 vs 124).
        HistFieldDef("VCMX25LUNA", "canopy-top LUNA vcmax25 (acclimated capacity)", "umol/m2/s", "patch",
            history_vcmx25luna_patch),
        HistFieldDef("JMX25LUNA", "canopy-top LUNA jmax25 (acclimated capacity)", "umol/m2/s", "patch",
            history_jmx25luna_patch),
        # Tree-only photosynthesis (isolates the evergreen tree from the grass-dominated
        # gridcell aggregate; for attributing the winter GPP/transpiration deficit).
        HistFieldDef("VCMX25TREE", "tree canopy-top LUNA vcmax25", "umol/m2/s", "patch",
            history_vcmx25tree),
        HistFieldDef("JMX25TREE", "tree canopy-top LUNA jmax25", "umol/m2/s", "patch",
            history_jmx25tree),
        HistFieldDef("FPSNTREE", "tree photosynthesis (GPP)", "umol/m2/s", "patch",
            history_fpsntree),
        HistFieldDef("PARVEGTREE", "tree canopy absorbed PAR", "W/m2", "patch",
            history_parvegtree),
        HistFieldDef("BTRANTREE", "tree instantaneous btran", "1", "patch",
            history_btrantree),
        HistFieldDef("PSNSUNTREE", "tree sunlit photosynthesis", "umol/m2/s", "patch",
            history_psnsuntree),
        HistFieldDef("ANSUNTREE", "tree sunlit net photosynthesis", "umol/m2/s", "patch",
            history_ansuntree),
        HistFieldDef("KSRTREE", "tree soil-root conductance (sum top 10 layers)", "mm/s", "patch",
            history_ksrtree),
        # Grass-only (winter-grass-dormancy chase)
        HistFieldDef("VCMX25GRASS", "grass canopy-top LUNA vcmax25", "umol/m2/s", "patch",
            history_vcmx25grass),
        HistFieldDef("FPSNGRASS", "grass photosynthesis (GPP)", "umol/m2/s", "patch",
            history_fpsngrass),
        HistFieldDef("BSUNGRASS", "grass sunlit btran (bsun)", "1", "patch",
            history_bsungrass),
        HistFieldDef("BTRANGRASS", "grass instantaneous btran", "1", "patch",
            history_btrangrass),
        HistFieldDef("VEGWPSUNGRASS", "grass sunlit leaf water potential", "mm", "patch",
            history_vegwpsungrass),
        HistFieldDef("ELAIGRASS", "grass exposed LAI", "m2/m2", "patch",
            history_elaigrass),
        HistFieldDef("TLAIGRASS", "grass total LAI", "m2/m2", "patch",
            history_tlaigrass),
        HistFieldDef("QOVER", "surface runoff", "mm/s", "column",
            inst -> inst.water.waterfluxbulk_inst.wf.qflx_surf_col),
        HistFieldDef("FSAT", "saturated fraction", "unitless", "column",
            inst -> inst.sat_excess_runoff.fsat_col),
        HistFieldDef("QFLX_RAIN_PLUS_SNOMELT", "rain plus snowmelt", "mm/s", "column",
            inst -> inst.water.waterfluxbulk_inst.wf.qflx_rain_plus_snomelt_col),
        HistFieldDef("QFLX_SNOMELT", "snowmelt", "mm/s", "column",
            inst -> inst.water.waterfluxbulk_inst.wf.qflx_snomelt_col),
        HistFieldDef("FCTR", "transpiration heat flux", "W/m2", "patch",
            inst -> inst.energyflux.eflx_lh_vegt_patch),
        HistFieldDef("FCEV", "canopy evaporation heat flux", "W/m2", "patch",
            inst -> inst.energyflux.eflx_lh_vege_patch),
        HistFieldDef("FGEV", "ground evaporation heat flux", "W/m2", "patch",
            inst -> inst.energyflux.eflx_lh_grnd_patch),
        HistFieldDef("FSH", "sensible heat flux", "W/m2", "patch",
            inst -> inst.energyflux.eflx_sh_tot_patch),
        HistFieldDef("QVEGT", "canopy transpiration", "mm/s", "column",
            inst -> inst.water.waterfluxbulk_inst.wf.qflx_tran_veg_col),
        HistFieldDef("QVEGE", "canopy evaporation", "mm/s", "column",
            history_qvege_col),
        HistFieldDef("QSOIL", "ground evaporation", "mm/s", "column",
            inst -> inst.water.waterfluxbulk_inst.wf.qflx_evap_soi_col),
    ]
end

"""
    history_h2osno_total_col(inst) -> Vector{Float64}

Return total snow water equivalent (mm) per column, including both unresolved
snow (`h2osno_no_layers_col`) and explicit snow layers.
"""
function history_h2osno_total_col(inst::CLMInstances)
    col = inst.column
    ws = inst.water.waterstatebulk_inst.ws
    nc = length(col.snl)
    h2osno_total = zeros(nc)
    waterstate_calculate_total_h2osno!(ws, trues(nc), 1:nc, col.snl, h2osno_total)
    return h2osno_total
end

"""
    history_btran_daily_min_patch(inst) -> Vector{Float64}

Return daily minimum BTRAN to align with Fortran `BTRANMN`.

Preference order:
1. `btran_min_patch` when finite (COMPLETED-day minimum — constant over the day)
2. `btran_min_inst_patch` when finite (running minimum; only the first day, before
   any completed-day min exists)
3. fallback to instantaneous `btran_patch`

CRITICAL: must use the COMPLETED daily minimum (`btran_min_patch`), not the running
one (`btran_min_inst_patch`). The history writer DAILY-AVERAGES this field; the
completed min is constant over the day so its average equals the daily min (matching
Fortran's `BTRANMN`, which points at `btran_min_patch` with avgflag='A'). The running
min varies over the day (high until the midday btran dip, then drops), so averaging it
overstates the daily min by ~5x (the bug that made annual BTRAN read +150%).
"""
function history_btran_daily_min_patch(inst::CLMInstances)
    ef = inst.energyflux
    pch = inst.patch
    np = length(ef.btran_patch)

    out = fill(NaN, np)
    for p in 1:np
        p <= length(pch.itype) || continue
        # Fortran BTRANMN uses a vegetation mask: exclude bare-ground patches.
        if pch.itype[p] != noveg
            bt = ef.btran_patch[p]
            if length(ef.btran_min_patch) >= p && isfinite(ef.btran_min_patch[p]) &&
               ef.btran_min_patch[p] != SPVAL
                bt = ef.btran_min_patch[p]
            elseif length(ef.btran_min_inst_patch) >= p && isfinite(ef.btran_min_inst_patch[p]) &&
                   ef.btran_min_inst_patch[p] != SPVAL
                bt = ef.btran_min_inst_patch[p]
            end
            out[p] = bt
        end
    end
    return out
end

# --- PHS diagnostic accessors (DBG_PHS_FIELDS) -------------------------------
# Temporary diagnostics for the late-season PHS over-transpiration investigation.
# Expose the chain smp_l → k_soil_root → vegwp → bsun, plus soil ice, so the
# Julia annual run can be compared link-by-link to the Fortran h0 VEGWP/BTRANMN.
_phs_vegwp_seg(inst, seg) = begin
    cs = inst.canopystate
    isempty(cs.vegwp_patch) ? Float64[] : Float64.(@view cs.vegwp_patch[:, seg])
end
history_vegwpsun_patch(inst) = _phs_vegwp_seg(inst, 1)
history_vegwproot_patch(inst) = _phs_vegwp_seg(inst, 4)
function history_ksoilroot_patch(inst::CLMInstances)
    ss = inst.soilstate
    isempty(ss.k_soil_root_patch) && return Float64[]
    np, nl = size(ss.k_soil_root_patch)
    out = zeros(np)
    @inbounds for p in 1:np
        s = 0.0
        for j in 1:nl
            v = ss.k_soil_root_patch[p, j]
            isfinite(v) && (s += v)
        end
        out[p] = s
    end
    return out
end
function history_smplmean_col(inst::CLMInstances)
    ss = inst.soilstate
    isempty(ss.smp_l_col) && return Float64[]
    nc = size(ss.smp_l_col, 1)
    nsoi = min(varpar.nlevsoi, size(ss.smp_l_col, 2))
    out = fill(NaN, nc)
    @inbounds for c in 1:nc
        s = 0.0; n = 0
        for j in 1:nsoi
            v = ss.smp_l_col[c, j]
            isfinite(v) && (s += v; n += 1)
        end
        n > 0 && (out[c] = s / n)
    end
    return out
end
function history_soilwater_10cm_col(inst::CLMInstances)
    ws = inst.water.waterstatebulk_inst.ws
    col = inst.column
    isempty(ws.h2osoi_liq_col) && return Float64[]
    nc = size(ws.h2osoi_liq_col, 1)
    nsno = varpar.nlevsno
    nsoi = varpar.nlevsoi
    out = zeros(nc)
    @inbounds for c in 1:nc
        s = 0.0
        for j in 1:nsoi
            jj = j + nsno
            zi_top = col.zi[c, j + nsno]      # interface above layer j (m)
            zi_top >= 0.1 && break
            # fraction of this layer within top 0.1 m
            zi_bot = col.zi[c, j + nsno + 1]
            frac = zi_bot <= 0.1 ? 1.0 : (0.1 - zi_top) / max(1e-9, zi_bot - zi_top)
            liq = ws.h2osoi_liq_col[c, jj]
            ice = ws.h2osoi_ice_col[c, jj]
            isfinite(liq) && (s += liq * frac)
            isfinite(ice) && (s += ice * frac)
        end
        out[c] = s
    end
    return out
end

# Fortran QVEGE = canopy INTERCEPTION evaporation = qflx_evap_veg - qflx_tran_veg
# (qflx_evap_veg_col is the TOTAL leaf water flux incl. transpiration). Matches the
# energy form FCEV = hvap*(qflx_evap_veg - qflx_tran_veg).
# FPSN = canopy gross photosynthesis (umol CO2 m-2 s-1) = fpsn_patch
# (= psnsun*laisun + psnsha*laisha). NaN/non-exposed → 0 (no photosynthesis), matching
# Fortran where psn=0 at night; the daily mean then includes night zeros.
# Canopy absorbed PAR (W/m2) = sum over radiative layers of parsun_z + parsha_z
# (= Fortran parveg = solad_vis*fabd_vis + solai_vis*fabi_vis, the VIS band). Daily
# mean here (Fortran PARVEGLN is local-noon, so expect Julia daily-mean < that).
function history_parveg_patch(inst::CLMInstances)
    sa = inst.solarabs
    isempty(sa.parsun_z_patch) && return Float64[]
    np = size(sa.parsun_z_patch, 1)
    nz = size(sa.parsun_z_patch, 2)
    nrad = inst.surfalb.nrad_patch
    out = zeros(np)
    @inbounds for p in 1:np
        nr = (p <= length(nrad) && nrad[p] >= 1) ? min(nrad[p], nz) : nz
        s = 0.0
        for z in 1:nr
            vsun = sa.parsun_z_patch[p, z]; vsha = sa.parsha_z_patch[p, z]
            isfinite(vsun) && (s += vsun)
            isfinite(vsha) && (s += vsha)
        end
        out[p] = s
    end
    return out
end

# LUNA-acclimated capacity per patch, top layer; bare → NaN (excluded by the
# gridcell aggregator → veg-area-weighted). The h0 VCMX25T is a different quantity
# (vcmax_opt=3), so compare per-patch, not gridcell-vs-h0.
function _history_lunacap(arr, pch)
    isempty(arr) && return Float64[]
    np = size(arr, 1)
    out = fill(NaN, np)
    @inbounds for p in 1:np
        (p <= length(pch.itype) && pch.itype[p] == noveg) && continue
        v = arr[p, 1]
        isfinite(v) && (out[p] = Float64(v))
    end
    return out
end
history_vcmx25luna_patch(inst) = _history_lunacap(inst.photosyns.vcmx25_z_patch, inst.patch)
history_jmx25luna_patch(inst) = _history_lunacap(inst.photosyns.jmx25_z_patch, inst.patch)

# Tree-only (trees+shrubs, itype 1:11) canopy-top value — isolates the evergreen
# tree from the grass-dominated gridcell aggregate. In winter the dormant grass
# carries a high vcmax25 *capacity* that contaminates the gridcell mean, yet it
# contributes ~0 to the GPP flux; these fields expose the tree's actual winter
# photosynthesis so the winter GPP/transpiration deficit can be attributed.
function _history_tree(arr, pch)
    isempty(arr) && return Float64[]
    np = ndims(arr) == 2 ? size(arr, 1) : length(arr)
    out = fill(NaN, np)
    @inbounds for p in 1:np
        p <= length(pch.itype) || continue
        it = pch.itype[p]
        (1 <= it <= 11) || continue
        v = ndims(arr) == 2 ? arr[p, 1] : arr[p]
        isfinite(v) && (out[p] = Float64(v))
    end
    return out
end
history_vcmx25tree(inst) = _history_tree(inst.photosyns.vcmx25_z_patch, inst.patch)
history_jmx25tree(inst)  = _history_tree(inst.photosyns.jmx25_z_patch, inst.patch)
history_fpsntree(inst)   = _history_tree(inst.photosyns.fpsn_patch, inst.patch)
history_parvegtree(inst) = _history_tree(history_parveg_patch(inst), inst.patch)
history_btrantree(inst)  = _history_tree(inst.energyflux.btran_patch, inst.patch)
# Realized tree sunlit photosynthesis (psnsun) + net (an_sun) for the winter-GPP
# attribution. (vcmax_z is a scratch array cleared after the solve, so it is not
# readable from end-of-step history — use the realized rates + PAR/btran inputs.)
history_psnsuntree(inst) = _history_tree(inst.photosyns.psnsun_patch, inst.patch)
history_ansuntree(inst)  = _history_tree(inst.photosyns.an_sun_patch, inst.patch)
# Tree total soil-root conductance (sum over the top 10 soil layers, to match the
# Fortran FCALCSTRESS dump). In frozen soil this tiny conductance sets the intraday
# vegwp draw-down depth → the BTRAN daily-min. Fortran tree winter sum ~2.2e-12 mm/s.
function history_ksrtree(inst)
    ksr = inst.soilstate.k_soil_root_patch
    isempty(ksr) && return Float64[]
    np = size(ksr, 1); nl = min(10, size(ksr, 2))
    s = [sum(@view ksr[p, 1:nl]) for p in 1:np]
    return _history_tree(s, inst.patch)
end

# Grass-only (itype 12:15) canopy value — to chase the winter grass dormancy
# (Julia's winter grass goes dormant, GPP~0, while Fortran's stays lightly active,
# an_sun~0.08/bsun~0.10). bsun<0.01 is zeroed in calcstress!, so a slightly more
# negative grass vegwp tips Julia's grass below the threshold into dormancy.
function _history_grass(arr, pch)
    isempty(arr) && return Float64[]
    np = ndims(arr) == 2 ? size(arr, 1) : length(arr)
    out = fill(NaN, np)
    @inbounds for p in 1:np
        p <= length(pch.itype) || continue
        it = pch.itype[p]
        (12 <= it <= 15) || continue
        v = ndims(arr) == 2 ? arr[p, 1] : arr[p]
        isfinite(v) && (out[p] = Float64(v))
    end
    return out
end
history_vcmx25grass(inst)   = _history_grass(inst.photosyns.vcmx25_z_patch, inst.patch)
history_fpsngrass(inst)     = _history_grass(inst.photosyns.fpsn_patch, inst.patch)
history_bsungrass(inst)     = _history_grass(inst.energyflux.bsun_patch, inst.patch)
history_btrangrass(inst)    = _history_grass(inst.energyflux.btran_patch, inst.patch)
history_vegwpsungrass(inst) = _history_grass(inst.canopystate.vegwp_patch, inst.patch)
history_elaigrass(inst)     = _history_grass(inst.canopystate.elai_patch, inst.patch)
history_tlaigrass(inst)     = _history_grass(inst.canopystate.tlai_patch, inst.patch)

function history_fpsn_patch(inst::CLMInstances)
    fp = inst.photosyns.fpsn_patch
    isempty(fp) && return Float64[]
    pch = inst.patch
    n = length(fp)
    out = fill(NaN, n)
    @inbounds for p in 1:n
        # Bare-ground patches → NaN so the gridcell aggregator EXCLUDES them
        # (matches Fortran FPSN=SPVAL on bare). Vegetated patches keep their value
        # (a finite 0 at night, included as a real zero). Returning 0 for bare
        # diluted the gridcell mean by the bare-area fraction (~-6% on FPSN).
        (p <= length(pch.itype) && pch.itype[p] == noveg) && continue
        v = fp[p]
        out[p] = isfinite(v) ? Float64(v) : 0.0
    end
    return out
end

function history_qvege_col(inst::CLMInstances)
    wf = inst.water.waterfluxbulk_inst.wf
    ev = wf.qflx_evap_veg_col; tr = wf.qflx_tran_veg_col
    isempty(ev) && return Float64[]
    n = length(ev)
    out = zeros(n)
    @inbounds for c in 1:n
        e = ev[c]; t = c <= length(tr) ? tr[c] : 0.0
        out[c] = (isfinite(e) ? e : 0.0) - (isfinite(t) ? t : 0.0)
    end
    return out
end

function history_tsoi_10cm_col(inst::CLMInstances)
    temp = inst.temperature
    col = inst.column
    isempty(temp.t_soisno_col) && return Float64[]
    nc = size(temp.t_soisno_col, 1)
    nsno = varpar.nlevsno
    nsoi = varpar.nlevsoi
    out = fill(NaN, nc)
    @inbounds for c in 1:nc
        s = 0.0; wsum = 0.0
        for j in 1:nsoi
            zi_top = col.zi[c, j + nsno]
            zi_top >= 0.1 && break
            zi_bot = col.zi[c, j + nsno + 1]
            dz_in = (zi_bot <= 0.1 ? zi_bot : 0.1) - zi_top
            t = temp.t_soisno_col[c, j + nsno]
            isfinite(t) && (s += t * dz_in; wsum += dz_in)
        end
        wsum > 0 && (out[c] = s / wsum)
    end
    return out
end

function history_soilice_col(inst::CLMInstances)
    ws = inst.water.waterstatebulk_inst.ws
    isempty(ws.h2osoi_ice_col) && return Float64[]
    nc = size(ws.h2osoi_ice_col, 1)
    nsno = varpar.nlevsno
    nsoi = varpar.nlevsoi
    out = zeros(nc)
    @inbounds for c in 1:nc
        s = 0.0
        for j in 1:nsoi
            jj = j + nsno
            jj <= size(ws.h2osoi_ice_col, 2) || continue
            v = ws.h2osoi_ice_col[c, jj]
            isfinite(v) && (s += v)
        end
        out[c] = s
    end
    return out
end

"""
    history_zwt_veg_col(inst) -> Vector{Float64}

Return water table depth for vegetated/crop landunits only, matching Fortran
`ZWT` landunit mask semantics.
"""
function history_zwt_veg_col(inst::CLMInstances)
    sh = inst.soilhydrology
    col = inst.column
    nc = length(sh.zwt_col)
    out = fill(NaN, nc)

    # Fortran ZWT uses a vegetation mask at column level.
    # Mark columns with at least one vegetated patch (itype != noveg).
    has_veg_col = falses(nc)
    pch = inst.patch
    for p in eachindex(pch.column)
        p <= length(pch.itype) || continue
        c = pch.column[p]
        c >= 1 && c <= nc || continue
        if pch.itype[p] != noveg
            has_veg_col[c] = true
        end
    end

    for c in 1:nc
        c <= length(col.landunit) || continue
        if has_veg_col[c]
            out[c] = sh.zwt_col[c]
        end
    end
    return out
end

"""
    history_writer_init!(hw, filepath, fields, ng, nc, np)

Create the NetCDF output file with dimensions and variable definitions.
All fields are written on the gridcell ("lndgrid") dimension after area-weighted
aggregation, matching Fortran CLM h0 output format.
"""
function history_writer_init!(hw::HistoryWriter, filepath::String,
                              fields::Vector{HistFieldDef},
                              ng::Int, nc::Int, np::Int)
    hw.filepath = filepath
    hw.fields = fields
    hw.time_index = 0

    hw.ds = NCDataset(filepath, "c")
    ds = hw.ds

    # Define dimensions — all output is on lndgrid (gridcell) level
    defDim(ds, "lndgrid", ng)
    defDim(ds, "time", Inf)  # unlimited

    # Time variable
    defVar(ds, "time", Float64, ("time",);
           attrib = Dict("units" => "days since 2000-01-01",
                         "calendar" => "noleap"))

    # Define each field on the gridcell dimension
    for f in fields
        defVar(ds, f.name, Float64, ("lndgrid", "time");
               attrib = Dict("long_name" => f.long_name,
                             "units" => f.units,
                             "subgrid_level" => f.level))
    end

    return nothing
end

"""
    _aggregate_to_gridcell(f, inst, ng) -> Vector{Float64}

Aggregate a field from its native level (patch/column/gridcell) to gridcell
level using area weights with renormalization.
"""
function _aggregate_to_gridcell(f::HistFieldDef, inst::CLMInstances, ng::Int)
    data = try
        f.getter(inst)
    catch e
        @debug "History field $(f.name) unavailable: $e"
        nothing
    end
    data === nothing && return nothing
    length(data) == 0 && return nothing

    gcell_val = zeros(ng)
    gcell_wt  = zeros(ng)

    pch = inst.patch
    col = inst.column

    if f.level == "patch"
        np = length(data)
        for p in 1:np
            val = Float64(data[p])
            if isfinite(val) && p <= length(pch.gridcell) && p <= length(pch.wtgcell)
                g = pch.gridcell[p]
                if g >= 1 && g <= ng
                    gcell_val[g] += val * pch.wtgcell[p]
                    gcell_wt[g]  += pch.wtgcell[p]
                end
            end
        end
    elseif f.level == "column"
        nc = length(data)
        for c in 1:nc
            val = Float64(data[c])
            if isfinite(val) && c <= length(col.gridcell) && c <= length(col.wtgcell)
                g = col.gridcell[c]
                if g >= 1 && g <= ng
                    gcell_val[g] += val * col.wtgcell[c]
                    gcell_wt[g]  += col.wtgcell[c]
                end
            end
        end
    else
        for g in 1:min(ng, length(data))
            val = Float64(data[g])
            gcell_val[g] = isfinite(val) ? val : -9999.0
            gcell_wt[g]  = isfinite(val) ? 1.0 : 0.0
        end
    end

    for g in 1:ng
        if gcell_wt[g] > 0.0
            gcell_val[g] /= gcell_wt[g]
        else
            gcell_val[g] = -9999.0
        end
    end

    return gcell_val
end

"""
    history_write_step!(hw, inst, current_time; is_end_curr_day=false)

Accumulate one timestep of all registered fields. When `is_end_curr_day` is
true, write the daily average to the output file and reset accumulators.
Matches Fortran CLM h0 daily-average output convention.
"""
# Days since 2000-01-01 on the NOLEAP (365-day) calendar — matching the `calendar=
# "noleap"` attribute on the time variable and CLM's datm convention. The previous
# code wrote Gregorian DateTime day-counts under a noleap label, so noleap-aware
# readers (and the CLM h0 reference) mis-decoded the stamps by one day (the 2000 leap
# day), making an index-by-index comparison off by a day. Drop Feb 29 in leap years.
function _noleap_days_since_2000(dt::DateTime)
    y = Dates.year(dt)
    doy = Dates.dayofyear(dt)
    if Dates.isleapyear(y) && doy > 59   # day 60 = Feb 29; collapse onto the noleap axis
        doy -= 1
    end
    frac = (Dates.hour(dt) * 3600 + Dates.minute(dt) * 60 + Dates.second(dt)) / 86400.0
    return (y - 2000) * 365 + (doy - 1) + frac
end

function history_write_step!(hw::HistoryWriter, inst::CLMInstances,
                             current_time::DateTime;
                             is_end_curr_day::Bool=false)
    ds = hw.ds
    ds === nothing && return nothing

    col = inst.column
    ng = length(col.gridcell) > 0 ? maximum(col.gridcell) : 0

    # Accumulate gridcell-level values for each field
    for f in hw.fields
        gcv = _aggregate_to_gridcell(f, inst, ng)
        gcv === nothing && continue

        if !haskey(hw.accum, f.name)
            hw.accum[f.name] = zeros(ng)
        end
        for g in 1:ng
            if gcv[g] != -9999.0
                hw.accum[f.name][g] += gcv[g]
            end
        end
    end
    hw.accum_count += 1

    # Write daily average at end of day
    if is_end_curr_day && hw.accum_count > 0
        hw.time_index += 1
        ti = hw.time_index

        # Write time = the daily-average interval END (next midnight) on the noleap
        # calendar — matching CLM, which stamps day D's average at D+1 00:00. (CLM.jl's
        # is_end_curr_day fires on the day's last sub-step, so current_time is e.g.
        # 23:30; the interval end is the following midnight.)
        interval_end = DateTime(Dates.year(current_time), Dates.month(current_time),
                                Dates.day(current_time)) + Dates.Day(1)
        ds["time"][ti] = _noleap_days_since_2000(interval_end)

        n = hw.accum_count
        for f in hw.fields
            if f.reduce == :snapshot
                # Daily-min/end-of-day field: write its value AT THE FLUSH (the
                # completed daily statistic), not the day-average of the accumulator.
                snap = _aggregate_to_gridcell(f, inst, ng)
                snap === nothing && continue
                ds[f.name][:, ti] = snap
            else
                haskey(hw.accum, f.name) || continue
                ds[f.name][:, ti] = hw.accum[f.name] ./ n
            end
        end

        # Reset accumulators
        for k in keys(hw.accum)
            fill!(hw.accum[k], 0.0)
        end
        hw.accum_count = 0
    end

    return nothing
end

"""
    history_writer_close!(hw)

Close the history NetCDF file.
"""
function history_writer_close!(hw::HistoryWriter)
    if hw.ds !== nothing
        # Flush partial-day data on close only if no complete day was written yet
        # (avoids spurious partial records at simulation end, but ensures short
        # sub-day runs still produce output)
        if hw.accum_count > 0 && hw.time_index == 0
            hw.time_index += 1
            ti = hw.time_index
            ds = hw.ds
            ds["time"][ti] = 0.0  # placeholder
            n = hw.accum_count
            for f in hw.fields
                haskey(hw.accum, f.name) || continue
                avg = hw.accum[f.name] ./ n
                ds[f.name][:, ti] = avg
            end
            hw.accum_count = 0
        end
        close(hw.ds)
        hw.ds = nothing
    end
    return nothing
end
