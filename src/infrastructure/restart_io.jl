# ==========================================================================
# Restart I/O — write/read the model's prognostic state to a CLM-style
# restart NetCDF file and round-trip it back.
#
# Ported from: src/main/restFileMod.F90 (file-level driver only).
#
# restFileMod's job in Fortran is to open the NetCDF, define the subgrid
# dimensions (gridcell/landunit/column/pft + the level dims levgrnd/levsno/
# levtot/levlak/levdcmp), then call each type's `*_Restart` subroutine which
# in turn issues `restartvar` calls (define on write, read on read). We mirror
# that structure here: `write_restart` / `read_restart!` are the file-level
# driver, and a per-variable registry plays the role of the scattered
# `restartvar` calls — each entry names the CLM variable, its subgrid level and
# vertical dimension, and a getter/setter onto the live `CLMInstances` tree.
#
# Variable names and dimension names follow CLM (T_SOISNO/T_GRND/T_VEG; dims
# column/pft/levtot/levlak/levgrnd/levdcmp/levsoi) where practical, for
# traceability back to the Fortran restart files.
#
# Scope: the prognostic biogeophysical state (temperature, soil/snow water,
# snow-layer geometry, canopy water, lake state, canopy state) plus, when CN
# is active, the soil-BGC and CN-vegetation carbon/nitrogen pools. This is a
# faithful subset — history/diagnostic-only fields are intentionally omitted.
#
# Public functions:
#   write_restart   — create the NetCDF and write the prognostic state
#   read_restart!   — read it back into a (fresh) CLMInstances
# ==========================================================================

const _RESTART_FILL = -9999.0

"""
    RestartVarDef

Definition of a single restart variable. `dims` is the rank of the array
(1=Vector, 2=Matrix, 3=3-D); `level` is the subgrid dimension the first axis
lives on ("column" or "patch"); `levdim` is the NetCDF name for the (optional)
vertical/pool axes — empty for rank-1 variables. `getter`/`setter!` read/write
the corresponding field on the `CLMInstances` tree.
"""
struct RestartVarDef
    name::String
    dims::Int
    level::String          # "column" or "patch"
    levdim::String         # NetCDF level-dim name (or "")
    leveldim2::String      # second level dim for rank-3 vars (or "")
    getter::Function       # (inst) -> AbstractArray
    setter!::Function      # (inst, data) -> nothing
end

# convenience constructors -----------------------------------------------------
_rv1(name, level, get, set) = RestartVarDef(name, 1, level, "", "", get, set)
_rv2(name, level, levdim, get, set) = RestartVarDef(name, 2, level, levdim, "", get, set)
_rv3(name, level, levdim, levdim2, get, set) =
    RestartVarDef(name, 3, level, levdim, levdim2, get, set)

"""
    _restart_registry_biogeophys() -> Vector{RestartVarDef}

Prognostic biogeophysical state: temperature, soil/snow water, snow-layer
geometry, canopy water, lake and canopy state. Always written.
"""
function _restart_registry_biogeophys()
    return RestartVarDef[
        # --- Temperature state ---
        _rv2("T_SOISNO", "column", "levtot",
            i -> i.temperature.t_soisno_col,
            (i, d) -> i.temperature.t_soisno_col .= d),
        _rv1("T_GRND", "column",
            i -> i.temperature.t_grnd_col,
            (i, d) -> i.temperature.t_grnd_col .= d),
        _rv2("T_LAKE", "column", "levlak",
            i -> i.temperature.t_lake_col,
            (i, d) -> i.temperature.t_lake_col .= d),
        _rv1("T_H2OSFC", "column",
            i -> i.temperature.t_h2osfc_col,
            (i, d) -> i.temperature.t_h2osfc_col .= d),
        _rv1("T_VEG", "patch",
            i -> i.temperature.t_veg_patch,
            (i, d) -> i.temperature.t_veg_patch .= d),
        _rv1("T_REF2M", "patch",
            i -> i.temperature.t_ref2m_patch,
            (i, d) -> i.temperature.t_ref2m_patch .= d),

        # --- Water state ---
        _rv2("H2OSOI_LIQ", "column", "levtot",
            i -> i.water.waterstatebulk_inst.ws.h2osoi_liq_col,
            (i, d) -> i.water.waterstatebulk_inst.ws.h2osoi_liq_col .= d),
        _rv2("H2OSOI_ICE", "column", "levtot",
            i -> i.water.waterstatebulk_inst.ws.h2osoi_ice_col,
            (i, d) -> i.water.waterstatebulk_inst.ws.h2osoi_ice_col .= d),
        _rv1("H2OSNO", "column",
            i -> i.water.waterstatebulk_inst.ws.h2osno_no_layers_col,
            (i, d) -> i.water.waterstatebulk_inst.ws.h2osno_no_layers_col .= d),
        _rv1("H2OSFC", "column",
            i -> i.water.waterstatebulk_inst.ws.h2osfc_col,
            (i, d) -> i.water.waterstatebulk_inst.ws.h2osfc_col .= d),
        _rv1("WA", "column",
            i -> i.water.waterstatebulk_inst.ws.wa_col,
            (i, d) -> i.water.waterstatebulk_inst.ws.wa_col .= d),
        _rv1("INT_SNOW", "column",
            i -> i.water.waterstatebulk_inst.int_snow_col,
            (i, d) -> i.water.waterstatebulk_inst.int_snow_col .= d),

        # --- Soil hydrology ---
        _rv1("ZWT", "column",
            i -> i.soilhydrology.zwt_col,
            (i, d) -> i.soilhydrology.zwt_col .= d),
        _rv1("ZWT_PERCH", "column",
            i -> i.soilhydrology.zwt_perched_col,
            (i, d) -> i.soilhydrology.zwt_perched_col .= d),

        # --- Snow state (counts + diagnostics + geometry) ---
        _rv1("SNLSNO", "column",
            i -> Float64.(i.column.snl),
            # snl is an integer count; its natural missing marker is ISPVAL, which
            # equals the float fill (-9999) and is restored to NaN on read — map
            # NaN back to ISPVAL rather than letting round(Int, NaN) throw.
            (i, d) -> i.column.snl .= (x -> isnan(x) ? ISPVAL : round(Int, x)).(d)),
        _rv1("SNOW_DEPTH", "column",
            i -> i.water.waterdiagnosticbulk_inst.snow_depth_col,
            (i, d) -> i.water.waterdiagnosticbulk_inst.snow_depth_col .= d),
        _rv1("frac_sno", "column",
            i -> i.water.waterdiagnosticbulk_inst.frac_sno_col,
            (i, d) -> i.water.waterdiagnosticbulk_inst.frac_sno_col .= d),
        _rv1("frac_sno_eff", "column",
            i -> i.water.waterdiagnosticbulk_inst.frac_sno_eff_col,
            (i, d) -> i.water.waterdiagnosticbulk_inst.frac_sno_eff_col .= d),
        _rv2("DZSNO", "column", "levtot",
            i -> i.column.dz,
            (i, d) -> i.column.dz .= d),
        _rv2("ZSNO", "column", "levtot",
            i -> i.column.z,
            (i, d) -> i.column.z .= d),
        _rv2("ZISNO", "column", "leviface",
            i -> i.column.zi,
            (i, d) -> i.column.zi .= d),

        # --- Canopy water ---
        _rv1("LIQCAN", "patch",
            i -> i.water.waterstatebulk_inst.ws.liqcan_patch,
            (i, d) -> i.water.waterstatebulk_inst.ws.liqcan_patch .= d),
        _rv1("SNOCAN", "patch",
            i -> i.water.waterstatebulk_inst.ws.snocan_patch,
            (i, d) -> i.water.waterstatebulk_inst.ws.snocan_patch .= d),

        # --- Lake state ---
        _rv2("LAKE_ICEFRAC", "column", "levlak",
            i -> i.lakestate.lake_icefrac_col,
            (i, d) -> i.lakestate.lake_icefrac_col .= d),
        _rv1("SAVEDTKE1", "column",
            i -> i.lakestate.savedtke1_col,
            (i, d) -> i.lakestate.savedtke1_col .= d),

        # --- Canopy state ---
        _rv1("TLAI", "patch",
            i -> i.canopystate.tlai_patch,
            (i, d) -> i.canopystate.tlai_patch .= d),
        _rv1("ELAI", "patch",
            i -> i.canopystate.elai_patch,
            (i, d) -> i.canopystate.elai_patch .= d),
        _rv1("HTOP", "patch",
            i -> i.canopystate.htop_patch,
            (i, d) -> i.canopystate.htop_patch .= d),
    ]
end

"""
    _restart_registry_cn() -> Vector{RestartVarDef}

CN prognostic state: soil-BGC vertically-resolved carbon/nitrogen pools and
mineral N, plus CN-vegetation carbon/nitrogen pools. Written only when CN is
active. Names follow CLM (`decomp_cpools_vr`, `sminn_vr`, `leafc`, …).
"""
function _restart_registry_cn()
    return RestartVarDef[
        # --- Soil BGC carbon (vertically + pool resolved) ---
        _rv3("decomp_cpools_vr", "column", "levdcmp", "ndecomp_pools",
            i -> i.soilbiogeochem_carbonstate.decomp_cpools_vr_col,
            (i, d) -> i.soilbiogeochem_carbonstate.decomp_cpools_vr_col .= d),
        _rv2("ctrunc_vr", "column", "levdcmp",
            i -> i.soilbiogeochem_carbonstate.ctrunc_vr_col,
            (i, d) -> i.soilbiogeochem_carbonstate.ctrunc_vr_col .= d),

        # --- Soil BGC nitrogen ---
        _rv3("decomp_npools_vr", "column", "levdcmp", "ndecomp_pools",
            i -> i.soilbiogeochem_nitrogenstate.decomp_npools_vr_col,
            (i, d) -> i.soilbiogeochem_nitrogenstate.decomp_npools_vr_col .= d),
        _rv2("sminn_vr", "column", "levdcmp",
            i -> i.soilbiogeochem_nitrogenstate.sminn_vr_col,
            (i, d) -> i.soilbiogeochem_nitrogenstate.sminn_vr_col .= d),
        _rv2("smin_no3_vr", "column", "levdcmp",
            i -> i.soilbiogeochem_nitrogenstate.smin_no3_vr_col,
            (i, d) -> i.soilbiogeochem_nitrogenstate.smin_no3_vr_col .= d),
        _rv2("smin_nh4_vr", "column", "levdcmp",
            i -> i.soilbiogeochem_nitrogenstate.smin_nh4_vr_col,
            (i, d) -> i.soilbiogeochem_nitrogenstate.smin_nh4_vr_col .= d),

        # --- CN vegetation carbon pools ---
        _rv1("leafc", "patch",
            i -> i.bgc_vegetation.cnveg_carbonstate_inst.leafc_patch,
            (i, d) -> i.bgc_vegetation.cnveg_carbonstate_inst.leafc_patch .= d),
        _rv1("frootc", "patch",
            i -> i.bgc_vegetation.cnveg_carbonstate_inst.frootc_patch,
            (i, d) -> i.bgc_vegetation.cnveg_carbonstate_inst.frootc_patch .= d),
        _rv1("livestemc", "patch",
            i -> i.bgc_vegetation.cnveg_carbonstate_inst.livestemc_patch,
            (i, d) -> i.bgc_vegetation.cnveg_carbonstate_inst.livestemc_patch .= d),
        _rv1("deadstemc", "patch",
            i -> i.bgc_vegetation.cnveg_carbonstate_inst.deadstemc_patch,
            (i, d) -> i.bgc_vegetation.cnveg_carbonstate_inst.deadstemc_patch .= d),
        _rv1("cpool", "patch",
            i -> i.bgc_vegetation.cnveg_carbonstate_inst.cpool_patch,
            (i, d) -> i.bgc_vegetation.cnveg_carbonstate_inst.cpool_patch .= d),

        # --- CN vegetation nitrogen pools ---
        _rv1("leafn", "patch",
            i -> i.bgc_vegetation.cnveg_nitrogenstate_inst.leafn_patch,
            (i, d) -> i.bgc_vegetation.cnveg_nitrogenstate_inst.leafn_patch .= d),
        _rv1("frootn", "patch",
            i -> i.bgc_vegetation.cnveg_nitrogenstate_inst.frootn_patch,
            (i, d) -> i.bgc_vegetation.cnveg_nitrogenstate_inst.frootn_patch .= d),
        _rv1("npool", "patch",
            i -> i.bgc_vegetation.cnveg_nitrogenstate_inst.npool_patch,
            (i, d) -> i.bgc_vegetation.cnveg_nitrogenstate_inst.npool_patch .= d),
    ]
end

# Integer/boolean setters mirror the `snl` handling above: the float round-trip
# carries an integer (or 0/1 boolean) value, but a fill/ISPVAL slot comes back as
# NaN — map NaN back to the type's missing marker rather than letting
# round(Int, NaN) throw. A value too large to represent as a Julia Int (e.g. the
# `typemax(Int)` planting-date sentinel, which a Float64 channel cannot carry
# losslessly) is also mapped to ISPVAL rather than overflowing round(Int, ·).
function _to_int_or_ispval(x)
    (isnan(x) || abs(x) >= 9.007199254740992e15) && return ISPVAL  # > 2^53 ⇒ not exact
    return round(Int, x)
end
_set_int!(field, d) = (field .= _to_int_or_ispval.(d); nothing)
_set_bool!(field, d) =
    (field .= (x -> isnan(x) ? false : (round(Int, clamp(x, 0.0, 1.0)) != 0)).(d); nothing)

"""
    _restart_registry_cn_flux() -> Vector{RestartVarDef}

Additional CN-vegetation prognostic state that must persist across restart: the
storage/transfer (`_storage`/`_xfer`) pools, the growth-respiration and excess-MR
pools, and the retranslocated-N pool. These are the flux/accumulator pools the
Fortran `CNVeg{Carbon,Nitrogen}StateType%Restart` write beyond the bulk pools
already covered by [`_restart_registry_cn`]. Written only when CN is active.
Names follow CLM's `restartvar` varnames (`xsmrpool`, `gresp_storage`, …).
"""
function _restart_registry_cn_flux()
    cs(field) = (i -> getproperty(i.bgc_vegetation.cnveg_carbonstate_inst, field),
                 (i, d) -> getproperty(i.bgc_vegetation.cnveg_carbonstate_inst, field) .= d)
    ns(field) = (i -> getproperty(i.bgc_vegetation.cnveg_nitrogenstate_inst, field),
                 (i, d) -> getproperty(i.bgc_vegetation.cnveg_nitrogenstate_inst, field) .= d)

    reg = RestartVarDef[]

    # --- Carbon storage/transfer pools + special C pools (all patch, 1-D) ---
    for (name, field) in (
        ("leafc_storage",      :leafc_storage_patch),
        ("leafc_xfer",         :leafc_xfer_patch),
        ("frootc_storage",     :frootc_storage_patch),
        ("frootc_xfer",        :frootc_xfer_patch),
        ("livestemc_storage",  :livestemc_storage_patch),
        ("livestemc_xfer",     :livestemc_xfer_patch),
        ("deadstemc_storage",  :deadstemc_storage_patch),
        ("deadstemc_xfer",     :deadstemc_xfer_patch),
        ("livecrootc",         :livecrootc_patch),
        ("livecrootc_storage", :livecrootc_storage_patch),
        ("livecrootc_xfer",    :livecrootc_xfer_patch),
        ("deadcrootc",         :deadcrootc_patch),
        ("deadcrootc_storage", :deadcrootc_storage_patch),
        ("deadcrootc_xfer",    :deadcrootc_xfer_patch),
        ("gresp_storage",      :gresp_storage_patch),
        ("gresp_xfer",         :gresp_xfer_patch),
        ("xsmrpool",           :xsmrpool_patch),
        ("xsmrpool_loss",      :xsmrpool_loss_patch),
        ("storage_cdemand",    :storage_cdemand_patch),
        ("leafc_storage_xfer_acc", :leafc_storage_xfer_acc_patch),
        ("pft_ctrunc",         :ctrunc_patch),
        ("cropseedc_deficit",  :cropseedc_deficit_patch),
    )
        g, s = cs(field)
        push!(reg, _rv1(name, "patch", g, s))
    end

    # --- Reproductive (grain) carbon pools (patch × nrepr, 2-D) ---
    for (name, field) in (
        ("reproductivec",         :reproductivec_patch),
        ("reproductivec_storage", :reproductivec_storage_patch),
        ("reproductivec_xfer",    :reproductivec_xfer_patch),
    )
        g, s = cs(field)
        push!(reg, _rv2(name, "patch", "nrepr", g, s))
    end

    # --- Nitrogen storage/transfer pools + retranslocated/special N pools ---
    for (name, field) in (
        ("leafn_storage",      :leafn_storage_patch),
        ("leafn_xfer",         :leafn_xfer_patch),
        ("frootn_storage",     :frootn_storage_patch),
        ("frootn_xfer",        :frootn_xfer_patch),
        ("livestemn",          :livestemn_patch),
        ("livestemn_storage",  :livestemn_storage_patch),
        ("livestemn_xfer",     :livestemn_xfer_patch),
        ("deadstemn",          :deadstemn_patch),
        ("deadstemn_storage",  :deadstemn_storage_patch),
        ("deadstemn_xfer",     :deadstemn_xfer_patch),
        ("livecrootn",         :livecrootn_patch),
        ("livecrootn_storage", :livecrootn_storage_patch),
        ("livecrootn_xfer",    :livecrootn_xfer_patch),
        ("deadcrootn",         :deadcrootn_patch),
        ("deadcrootn_storage", :deadcrootn_storage_patch),
        ("deadcrootn_xfer",    :deadcrootn_xfer_patch),
        ("retransn",           :retransn_patch),
        ("storage_ndemand",    :storage_ndemand_patch),
        ("pft_ntrunc",         :ntrunc_patch),
        ("cropseedn_deficit",  :cropseedn_deficit_patch),
    )
        g, s = ns(field)
        push!(reg, _rv1(name, "patch", g, s))
    end

    # --- Reproductive (grain) nitrogen pools (patch × nrepr, 2-D) ---
    for (name, field) in (
        ("reproductiven",         :reproductiven_patch),
        ("reproductiven_storage", :reproductiven_storage_patch),
        ("reproductiven_xfer",    :reproductiven_xfer_patch),
    )
        g, s = ns(field)
        push!(reg, _rv2(name, "patch", "nrepr", g, s))
    end

    return reg
end

"""
    _restart_registry_pheno() -> Vector{RestartVarDef}

CN phenology prognostic counters/accumulators (`CNVegStateType%Restart`,
always-on block): the dormancy/onset/offset state-machine flags and counters,
litterfall/transfer rates, T2m running means and the annual GPP/retransn sums.
Written when CN is active. Names follow CLM `restartvar` varnames.
"""
function _restart_registry_pheno()
    vs(field) = (i -> getproperty(i.bgc_vegetation.cnveg_state_inst, field),
                 (i, d) -> getproperty(i.bgc_vegetation.cnveg_state_inst, field) .= d)
    reg = RestartVarDef[]
    for (name, field) in (
        ("dormant_flag",          :dormant_flag_patch),
        ("days_active",           :days_active_patch),
        ("onset_flag",            :onset_flag_patch),
        ("onset_counter",         :onset_counter_patch),
        ("onset_gddflag",         :onset_gddflag_patch),
        ("onset_fdd",             :onset_fdd_patch),
        ("onset_gdd",             :onset_gdd_patch),
        ("onset_swi",             :onset_swi_patch),
        ("offset_flag",           :offset_flag_patch),
        ("offset_counter",        :offset_counter_patch),
        ("offset_fdd",            :offset_fdd_patch),
        ("offset_swi",            :offset_swi_patch),
        ("lgsf",                  :lgsf_patch),
        ("bglfr",                 :bglfr_patch),
        ("bgtr",                  :bgtr_patch),
        ("annavg_t2m",            :annavg_t2m_patch),
        ("tempavg_t2m",           :tempavg_t2m_patch),
        ("c_allometry",           :c_allometry_patch),
        ("n_allometry",           :n_allometry_patch),
        ("tempsum_potential_gpp", :tempsum_potential_gpp_patch),
        ("annsum_potential_gpp",  :annsum_potential_gpp_patch),
        ("tempmax_retransn",      :tempmax_retransn_patch),
        ("annmax_retransn",       :annmax_retransn_patch),
        ("downreg",               :downreg_patch),
    )
        g, s = vs(field)
        push!(reg, _rv1(name, "patch", g, s))
    end
    # Column-level annual accumulators
    push!(reg, _rv1("annsum_counter", "column",
        i -> i.bgc_vegetation.cnveg_state_inst.annsum_counter_col,
        (i, d) -> i.bgc_vegetation.cnveg_state_inst.annsum_counter_col .= d))
    push!(reg, _rv1("cannavg_t2m", "column",
        i -> i.bgc_vegetation.cnveg_state_inst.annavg_t2m_col,
        (i, d) -> i.bgc_vegetation.cnveg_state_inst.annavg_t2m_col .= d))
    return reg
end

"""
    _restart_registry_crop() -> Vector{RestartVarDef}

Crop prognostic phenology/scheduling state, written when `use_crop` is active.
Covers the crop block of `CNVegStateType%Restart` (planting/allocation/vernal-
ization counters) and `CropType%Restart` (crop-alive flag, harvest date, phase,
the per-sowing / per-harvest "thisyr" arrays). Integer/boolean fields use the
NaN↔ISPVAL/false round-trip helpers. Names follow CLM `restartvar` varnames.
"""
function _restart_registry_crop()
    vs(field) = (i -> getproperty(i.bgc_vegetation.cnveg_state_inst, field),
                 (i, d) -> getproperty(i.bgc_vegetation.cnveg_state_inst, field) .= d)
    reg = RestartVarDef[]

    # --- CNVegState crop block (patch, 1-D real) ---
    for (name, field) in (
        ("htmx",        :htmx_patch),
        ("aleaf",       :aleaf_patch),
        ("aleafi",      :aleafi_patch),
        ("astem",       :astem_patch),
        ("astemi",      :astemi_patch),
        ("hdidx",       :hdidx_patch),
        ("cumvd",       :cumvd_patch),
        ("gddmaturity", :gddmaturity_patch),
        ("huileaf",     :huileaf_patch),
        ("huigrain",    :huigrain_patch),
        ("grain_flag",  :grain_flag_patch),
    )
        g, s = vs(field)
        push!(reg, _rv1(name, "patch", g, s))
    end
    # CNVegState crop integer fields (patch)
    push!(reg, _rv1("peaklai", "patch",
        i -> Float64.(i.bgc_vegetation.cnveg_state_inst.peaklai_patch),
        (i, d) -> _set_int!(i.bgc_vegetation.cnveg_state_inst.peaklai_patch, d)))
    push!(reg, _rv1("idop", "patch",
        i -> Float64.(i.bgc_vegetation.cnveg_state_inst.idop_patch),
        (i, d) -> _set_int!(i.bgc_vegetation.cnveg_state_inst.idop_patch, d)))
    # gddmaturity_thisyr (patch × mxharvests, 2-D)
    push!(reg, _rv2("gddmaturity_thisyr", "patch", "mxharvests",
        i -> i.bgc_vegetation.cnveg_state_inst.gddmaturity_thisyr,
        (i, d) -> i.bgc_vegetation.cnveg_state_inst.gddmaturity_thisyr .= d))

    # --- CropType state ---
    cr(field) = (i -> getproperty(i.crop, field),
                 (i, d) -> getproperty(i.crop, field) .= d)
    # Real 1-D (patch)
    for (name, field) in (("vf", :vf_patch), ("cphase", :cphase_patch))
        g, s = cr(field)
        push!(reg, _rv1(name, "patch", g, s))
    end
    # Integer 1-D (patch)
    push!(reg, _rv1("nyrs_crop_active", "patch",
        i -> Float64.(i.crop.nyrs_crop_active_patch),
        (i, d) -> _set_int!(i.crop.nyrs_crop_active_patch, d)))
    push!(reg, _rv1("harvdate", "patch",
        i -> Float64.(i.crop.harvdate_patch),
        (i, d) -> _set_int!(i.crop.harvdate_patch, d)))
    push!(reg, _rv1("sowing_reason_patch", "patch",
        i -> Float64.(i.crop.sowing_reason_patch),
        (i, d) -> _set_int!(i.crop.sowing_reason_patch, d)))
    # Boolean 1-D (patch)
    push!(reg, _rv1("croplive", "patch",
        i -> Float64.(i.crop.croplive_patch),
        (i, d) -> _set_bool!(i.crop.croplive_patch, d)))
    push!(reg, _rv1("sown_in_this_window", "patch",
        i -> Float64.(i.crop.sown_in_this_window),
        (i, d) -> _set_bool!(i.crop.sown_in_this_window, d)))
    # Real 2-D thisyr / perharv arrays (patch × mxsowings | mxharvests)
    for (name, field) in (
        ("sdates_thisyr_patch",         :sdates_thisyr_patch),
        ("swindow_starts_thisyr_patch", :swindow_starts_thisyr_patch),
        ("swindow_ends_thisyr_patch",   :swindow_ends_thisyr_patch),
        ("sowing_reason_thisyr_patch",  :sowing_reason_thisyr_patch),
        ("sdates_perharv_patch",        :sdates_perharv_patch),
        ("syears_perharv_patch",        :syears_perharv_patch),
        ("hdates_thisyr_patch",         :hdates_thisyr_patch),
        ("gddaccum_thisyr_patch",       :gddaccum_thisyr_patch),
        ("hui_thisyr_patch",            :hui_thisyr_patch),
        ("sowing_reason_perharv_patch", :sowing_reason_perharv_patch),
        ("harvest_reason_thisyr_patch", :harvest_reason_thisyr_patch),
    )
        g, s = cr(field)
        push!(reg, _rv2(name, "patch", "mxgrowseason", g, s))
    end
    return reg
end

"""
    _restart_registry_isotope(suffix, getinst) -> Vector{RestartVarDef}

Carbon-isotope vegetation pools for one isotope species. `suffix` is the CLM
restart-var suffix (`_13` or `_14`); `getinst` returns the matching carbon-state
instance from the `CLMInstances` tree. Written only when the corresponding
`use_c13`/`use_c14` flag is active. The pool set mirrors C12 (bulk pools +
storage/transfer/special pools + reproductive pools).
"""
function _restart_registry_isotope(suffix::String, getinst::Function)
    cs(field) = (i -> getproperty(getinst(i), field),
                 (i, d) -> getproperty(getinst(i), field) .= d)
    reg = RestartVarDef[]
    for (base, field) in (
        ("leafc",              :leafc_patch),
        ("leafc_storage",      :leafc_storage_patch),
        ("leafc_xfer",         :leafc_xfer_patch),
        ("frootc",             :frootc_patch),
        ("frootc_storage",     :frootc_storage_patch),
        ("frootc_xfer",        :frootc_xfer_patch),
        ("livestemc",          :livestemc_patch),
        ("livestemc_storage",  :livestemc_storage_patch),
        ("livestemc_xfer",     :livestemc_xfer_patch),
        ("deadstemc",          :deadstemc_patch),
        ("deadstemc_storage",  :deadstemc_storage_patch),
        ("deadstemc_xfer",     :deadstemc_xfer_patch),
        ("livecrootc",         :livecrootc_patch),
        ("livecrootc_storage", :livecrootc_storage_patch),
        ("livecrootc_xfer",    :livecrootc_xfer_patch),
        ("deadcrootc",         :deadcrootc_patch),
        ("deadcrootc_storage", :deadcrootc_storage_patch),
        ("deadcrootc_xfer",    :deadcrootc_xfer_patch),
        ("gresp_storage",      :gresp_storage_patch),
        ("gresp_xfer",         :gresp_xfer_patch),
        ("cpool",              :cpool_patch),
        ("xsmrpool",           :xsmrpool_patch),
        ("pft_ctrunc",         :ctrunc_patch),
    )
        g, s = cs(field)
        push!(reg, _rv1(base * suffix, "patch", g, s))
    end
    for (base, field) in (
        ("reproductivec",         :reproductivec_patch),
        ("reproductivec_storage", :reproductivec_storage_patch),
        ("reproductivec_xfer",    :reproductivec_xfer_patch),
    )
        g, s = cs(field)
        push!(reg, _rv2(base * suffix, "patch", "nrepr", g, s))
    end
    return reg
end

# Replace non-finite entries with the NetCDF fill value (in place, on a copy).
function _clean_for_write!(a::AbstractArray)
    @inbounds for i in eachindex(a)
        isfinite(a[i]) || (a[i] = _RESTART_FILL)
    end
    return a
end

# Restore fill values back to NaN after read (in place).
function _restore_after_read!(a::AbstractArray)
    @inbounds for i in eachindex(a)
        a[i] == _RESTART_FILL && (a[i] = NaN)
    end
    return a
end

# Ensure a NetCDF dim of name `nm` exists with size `n`; if it exists with a
# different size, fall back to a unique per-variable name and define that.
function _ensure_dim!(ds, nm::AbstractString, n::Int, fallback::AbstractString)
    if haskey(ds.dim, nm)
        ds.dim[nm] == n && return nm
        # size clash → use a variable-unique dim name
        if !haskey(ds.dim, fallback)
            defDim(ds, fallback, n)
        end
        return fallback
    end
    defDim(ds, nm, n)
    return nm
end

"""
    _assemble_registry(; use_cn, use_crop, use_c13, use_c14) -> Vector{RestartVarDef}

Build the full ordered restart registry for the active configuration. The
biogeophysical block is always present; CN state/flux + phenology are added when
`use_cn`; the crop block when `use_crop`; and the C13/C14 isotope pools when the
matching isotope flag is set. Used by both `write_restart` and `read_restart!`
so the two stay in lock-step.
"""
function _assemble_registry(; use_cn::Bool, use_crop::Bool,
                            use_c13::Bool, use_c14::Bool)
    registry = _restart_registry_biogeophys()
    if use_cn
        append!(registry, _restart_registry_cn())
        append!(registry, _restart_registry_cn_flux())
        append!(registry, _restart_registry_pheno())
        use_c13 && append!(registry,
            _restart_registry_isotope("_13",
                i -> i.bgc_vegetation.c13_cnveg_carbonstate_inst))
        use_c14 && append!(registry,
            _restart_registry_isotope("_14",
                i -> i.bgc_vegetation.c14_cnveg_carbonstate_inst))
    end
    # Crop phenology/scheduling state is gated on use_crop (its always-on
    # phenology counters live with the CN block, mirroring CNVegStateType%Restart).
    use_crop && append!(registry, _restart_registry_crop())
    return registry
end

"""
    write_restart(inst, filename; bounds, use_cn=false, use_crop=false,
                  use_c13=false, use_c14=false, time=nothing)

Create a CLM-style restart NetCDF at `filename` and write the prognostic state
held in `inst`. Mirrors `restFileMod::restFile_write`: defines the subgrid and
level dimensions, then writes each registered variable (the role of the
per-type `restartvar` calls).

- `bounds`   — `BoundsType` giving `endc` (#columns) and `endp` (#patches).
- `use_cn`   — also write the CN carbon/nitrogen pools (bulk + storage/transfer/
  flux pools) and the phenology counters.
- `use_crop` — also write the crop/phenology prognostic counters.
- `use_c13`/`use_c14` — also write the carbon-isotope vegetation pools.
- `time`     — optional valid time, stored as a scalar `timemgr_rst_curr_date`.
"""
function write_restart(inst::CLMInstances, filename::String;
                       bounds::BoundsType,
                       use_cn::Bool = false,
                       use_crop::Bool = false,
                       use_c13::Bool = false,
                       use_c14::Bool = false,
                       time::Union{DateTime,Nothing} = nothing)
    nc = bounds.endc
    np = bounds.endp

    registry = _assemble_registry(use_cn = use_cn, use_crop = use_crop,
                                  use_c13 = use_c13, use_c14 = use_c14)

    ds = NCDataset(filename, "c")
    try
        # Subgrid dimensions (restFileMod defines column/pft; we add what we use).
        defDim(ds, "column", nc)
        defDim(ds, "pft", np)

        if time !== nothing
            tvar = defVar(ds, "timemgr_rst_curr_date", Float64, ())
            tvar[] = Dates.value(time - DateTime(2000, 1, 1)) / (1000 * 86400)
        end

        for rv in registry
            data = try
                rv.getter(inst)
            catch
                continue
            end
            (data === nothing || length(data) == 0) && continue

            dim1 = rv.level == "column" ? "column" : "pft"

            if rv.dims == 1
                v = defVar(ds, rv.name, Float64, (dim1,))
                clean = _clean_for_write!(collect(Float64, data))
                v[:] = clean
            elseif rv.dims == 2
                n2 = size(data, 2)
                d2 = _ensure_dim!(ds, rv.levdim, n2, "lev_$(rv.name)")
                v = defVar(ds, rv.name, Float64, (dim1, d2))
                clean = _clean_for_write!(collect(Float64, data))
                v[:, :] = clean
            else # rv.dims == 3
                n2 = size(data, 2); n3 = size(data, 3)
                d2 = _ensure_dim!(ds, rv.levdim,  n2, "lev_$(rv.name)")
                d3 = _ensure_dim!(ds, rv.leveldim2, n3, "pool_$(rv.name)")
                v = defVar(ds, rv.name, Float64, (dim1, d2, d3))
                clean = _clean_for_write!(collect(Float64, data))
                v[:, :, :] = clean
            end
        end
    finally
        close(ds)
    end
    return nothing
end

"""
    read_restart!(inst, filename; bounds, use_cn=false, use_crop=false,
                  use_c13=false, use_c14=false)

Read the prognostic state from a restart NetCDF written by [`write_restart`]
back into `inst`, overwriting each registered field. Mirrors
`restFileMod::restFile_read`. Variables absent from the file are left untouched;
fill values are restored to `NaN`. The `use_cn`/`use_crop`/`use_c13`/`use_c14`
flags select the same blocks they do on write — a flag left `false` leaves the
corresponding fields untouched.
"""
function read_restart!(inst::CLMInstances, filename::String;
                       bounds::BoundsType,
                       use_cn::Bool = false,
                       use_crop::Bool = false,
                       use_c13::Bool = false,
                       use_c14::Bool = false)
    isfile(filename) || error("Restart file not found: $filename")

    registry = _assemble_registry(use_cn = use_cn, use_crop = use_crop,
                                  use_c13 = use_c13, use_c14 = use_c14)

    NCDataset(filename, "r") do ds
        for rv in registry
            haskey(ds, rv.name) || continue
            raw = Array(ds[rv.name])
            data = Float64.(replace(raw, missing => NaN))
            _restore_after_read!(data)
            try
                rv.setter!(inst, data)
            catch e
                @warn "Skipping restart variable $(rv.name): $e" maxlog = 1
            end
        end
    end
    return nothing
end

# --- Backwards-compatible aliases (old arg order: filepath, inst, bounds) -----
# The earlier API was `write_restart!(filepath, inst, bounds)`. Keep thin
# shims so existing callers/scripts continue to work.
function write_restart!(filepath::String, inst::CLMInstances, bounds::BoundsType;
                        use_cn::Bool = false, use_crop::Bool = false,
                        use_c13::Bool = false, use_c14::Bool = false,
                        time::Union{DateTime,Nothing} = nothing)
    return write_restart(inst, filepath; bounds = bounds, use_cn = use_cn,
                         use_crop = use_crop, use_c13 = use_c13,
                         use_c14 = use_c14, time = time)
end

function read_restart!(filepath::String, inst::CLMInstances, bounds::BoundsType;
                       use_cn::Bool = false, use_crop::Bool = false,
                       use_c13::Bool = false, use_c14::Bool = false)
    return read_restart!(inst, filepath; bounds = bounds, use_cn = use_cn,
                         use_crop = use_crop, use_c13 = use_c13, use_c14 = use_c14)
end
