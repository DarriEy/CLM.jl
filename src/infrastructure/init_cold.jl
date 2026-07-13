# ==========================================================================
# init_cold.jl — the single dispatch point for Fortran's `InitCold` phase.
#
# WHY THIS FILE EXISTS
# --------------------
# Every CLM Fortran state type follows
#
#     subroutine Init(...)
#        call this%InitAllocate(bounds)   ! allocates + fills with NaN / spval
#        call this%InitHistory(bounds)    ! history registration
#        call this%InitCold(bounds)       ! ZEROES / SEEDS the prognostic fields
#     end subroutine
#
# The Julia port faithfully translated ~20 of those `InitCold` routines into
# `*_init_cold!` functions — and then NEVER CALLED most of them. `*_init!`
# (the InitAllocate analogue) NaN-fills, so the NaN survived into the run. This
# was a bug FACTORY, not a bug: three separate live failures were traced to it
# (waterflux -> `qflx_evap_tot_col` NaN on every urban roof/wall column, a live
# lnd2atm output; aerosol -> `mss_*` NaN -> SNICAR albedo NaN -> `sabg_lyr` NaN
# -> the glacier soil column's band solve NaN'd; activelayer -> altmax never
# tracked). Each was found only when a domain happened to blow up.
#
# This file closes the class. It routes EVERY ported `*_init_cold!` through two
# entry points on the live init path, and `test_init_cold_wiring.jl` asserts
# statically that no `*_init_cold!` defined in `src/` is missing from them — so
# a future dead port is caught at test time, not by a NaN three domains later.
#
# ORDERING (mirrors Fortran `clm_instMod.F90::clm_instInit`)
# ----------------------------------------------------------
# Fortran runs InitCold inside `Init`, i.e. AFTER the subgrid hierarchy is built
# (`initGridcells` precedes `clm_instInit`) and BEFORE any restart read (the
# restart file then overwrites whatever InitCold seeded). In CLM.jl the
# equivalent window is `cold_start_initialize!` (clm_initialize! Step 14): the
# subgrid + vertical grid are ready, and nothing has read a restart yet.
#
#   `init_cold_biogeophys!`  — called at the TOP of `cold_start_initialize!`.
#   `init_cold_biogeochem!`  — called from `clm_initialize!` AFTER the
#                              decomposition cascade is built, because Fortran's
#                              `soilbiogeochem_carbonstate_inst%Init` (whose
#                              InitCold seeds `decomp_cpools_vr` from
#                              `decomp_cascade_con%initial_stock`) runs AFTER
#                              `init_decompcascade_bgc` — clm_instMod.F90:414-425.
#
# WHY "InitCold FIRST, bespoke seeds SECOND"
# ------------------------------------------
# `cold_start_initialize!` already hand-rolls a lot of cold-start state
# (`init_temperatures!`, `init_soil_moisture!`, `init_misc_state!`, …). Several
# of those seeds intentionally DIVERGE from the plain Fortran InitCold because
# they were tuned against Fortran parity dumps (e.g. the 277 K lake branch, the
# physical step-1 `savedtke1`, the fresh-snow `snw_rds` = 54.526). Running
# `init_cold_biogeophys!` FIRST and letting the bespoke seeds overwrite on top
# means:
#   * every field the bespoke code owns keeps its validated value  -> Bow stays
#     BIT-IDENTICAL;
#   * every field NOTHING else touched gets its faithful Fortran InitCold value
#     instead of the allocator's NaN  -> the bug class is closed.
#
# WHAT WE DELIBERATELY LEAVE NON-FINITE
# -------------------------------------
# Fortran uses NaN/`spval` as a genuine MISSING-VALUE FLAG in places, and the
# history buffer masks those out. Zeroing them would be WRONG — a real 0.0 flux
# is not the same as "this column has no glacier mass balance". The ports keep
# them, and we do not second-guess them:
#   * `qflx_glcice*`      — only meaningful inside the `do_smb` filter.
#   * `t_ref2m_r/u`, `eflx_*_r/u` — the rural/urban split diagnostics are spval
#     on the landunit they do not apply to.
#   * `wf_col` / `wf2_col`, `snw_rds_top`, `sno_liq_top` — spval when undefined.
#   * `fsun_patch`        — spval so the accumulator averages correctly.
#   * `latbaset_patch`    — NaN unless `baset_mapping == LATVARY`.
#   * `alphapsnsun/sha`   — spval (PhotosynthesisMod::InitCold sets exactly this).
# ==========================================================================

"""
    init_cold_biogeophys!(inst, bounds)

Run every ported biogeophysics `*_init_cold!` — the `InitCold` phase of Fortran's
`Init` for the temperature / energy-flux / canopy / soil / water / radiation
types.

Called at the TOP of [`cold_start_initialize!`](@ref), so the bespoke Julia
cold-start seeds that follow it overwrite the fields they own (keeping validated
domains bit-identical) while the fields nothing else writes get their faithful
Fortran InitCold value instead of the allocator's NaN.

Requires: the subgrid hierarchy (`initGridCells!`), the vertical grid
(`initVertical!`), `readParameters!`, and — when urban landunits are present —
`urbanparams_populate!`. All of those precede Step 14 in `clm_initialize!`.
"""
function init_cold_biogeophys!(inst::CLMInstances, bounds::BoundsType;
                               use_aquifer_layer::Bool = true)
    col = inst.column
    lun = inst.landunit
    pch = inst.patch
    bc  = bounds.begc:bounds.endc
    bp  = bounds.begp:bounds.endp
    bl  = bounds.begl:bounds.endl
    bg  = bounds.begg:bounds.endg

    # Seed BOTH building-temperature methods' fields.
    #
    # Fortran picks one (`is_simple_buildtemp` / `is_prog_buildtemp`) from a
    # namelist that UrbanParamsType::Init reads BEFORE temperature_inst%Init
    # (clm_instMod.F90:268 vs :315), so its InitCold knows the answer. CLM.jl
    # reads the urban namelist AFTER cold_start_initialize!, and
    # `urban_ctrl.building_temp_method` stays mutable afterwards — scripts/
    # urban_smoke.jl flips it to SIMPLE once clm_initialize! has already
    # returned. Committing to one method here would therefore leave the OTHER
    # method's fields at NaN for any caller that switches.
    #
    # Both branches only ever WRITE (they never read each other's fields), so
    # seeding both is safe: a field the active method never reads is inert at
    # 0.0 and lethal at NaN. This is the one place we deliberately do more than
    # Fortran, and it is strictly a superset.
    is_prog   = true
    is_simple = true

    # ---- TemperatureType::InitCold -----------------------------------------
    # Must precede waterstate_init_cold! (which splits h2osoi_vol into liq/ice by
    # t_soisno) and energyflux_init_cold! (which needs t_grnd_col), exactly as
    # Fortran orders temperature_inst%Init before energyflux_inst%Init
    # (clm_instMod.F90:315 vs :345). `col.snl` must already be 0 — see the
    # init_snow_state! hoist in cold_start_initialize!.
    up = inst.urbanparams
    temperature_init_cold!(inst.temperature, col, lun, pch, bc, bp, bl;
                           em_roof_lun    = up.em_roof,
                           em_wall_lun    = up.em_wall,
                           em_improad_lun = up.em_improad,
                           em_perroad_lun = up.em_perroad,
                           is_prog_buildtemp = is_prog)

    # ---- CanopyStateType::InitCold -----------------------------------------
    # Seeds vegwp = -2.5e4 mm (the PHS plant water potential — a REAL physics
    # input under use_hydrstress, previously NaN on every active patch), the
    # leaf/stem biomass, and the *_hist diagnostics. Runs BEFORE
    # init_satellite_phenology!, which then installs the real LAI/SAI/htop.
    canopystate_init_cold!(inst.canopystate, bp;
                           landunit_patch = pch.landunit,
                           lun_itype      = lun.itype)

    # ---- SoilStateType::InitCold -------------------------------------------
    soilstate_init_cold!(inst.soilstate, bc)

    # ---- EnergyFluxType::InitCold ------------------------------------------
    # Zeroes the urban anthropogenic-flux terms (wasteheat / traffic /
    # ventilation / heat-from-AC) on NON-urban patches and `rresis_patch`; sets
    # the rural/urban split diagnostics to spval on the landunit they do not
    # apply to. Those were NaN on every active patch of every domain.
    energyflux_init_cold!(inst.energyflux, col, lun, pch,
                          inst.temperature.t_grnd_col, bc, bp, bl;
                          is_simple_buildtemp = is_simple,
                          is_prog_buildtemp   = is_prog)

    # ---- SolarAbsorbedType::InitCold ---------------------------------------
    # Zeroes the urban sabs_* absorbed fractions. Runs BEFORE
    # init_surface_albedo_cold!, whose `_init_urban_albedo_cold!` then seeds the
    # physically-bounded (1 - albedo) values for the ACTUAL urban landunits.
    solarabs_init_cold!(inst.solarabs, bl)

    # ---- SurfaceRadiationMod::InitCold -------------------------------------
    # Genuinely a no-op in Fortran too (SurfaceRadiationMod.F90). Called anyway
    # so the wiring assertion holds and the pairing stays explicit.
    surfrad_init_cold!(inst.surfrad, bp)

    # ---- PhotosynthesisMod::InitCold ---------------------------------------
    photosyns_init_cold!(inst.photosyns, pch, lun, bp)

    # ---- WaterStateType / WaterStateBulkType::InitCold ----------------------
    # Cold start has no snow, so h2osno_input / snow_depth_input are 0.
    ws  = inst.water.waterstatebulk_inst.ws
    wsb = inst.water.waterstatebulk_inst
    h2osno_input     = zeros(Float64, bounds.endc)
    snow_depth_input = zeros(Float64, bounds.endc)
    waterstate_init_cold!(ws, bc, bp, bl, bg;
                          h2osno_input_col = h2osno_input,
                          watsat_col       = Matrix(inst.soilstate.watsat_col),
                          t_soisno_col     = Matrix(inst.temperature.t_soisno_col),
                          use_aquifer_layer = use_aquifer_layer,
                          ratio            = 1.0,
                          snl_col          = Vector{Int}(col.snl),
                          dz_col           = Matrix(col.dz),
                          landunit_col     = Vector{Int}(col.landunit),
                          lakpoi           = BitVector(lun.lakpoi),
                          urbpoi           = BitVector(lun.urbpoi),
                          lun_itype        = Vector{Int}(lun.itype),
                          col_itype        = Vector{Int}(col.itype),
                          nbedrock_col     = Vector{Int}(col.nbedrock),
                          gridcell_col     = Vector{Int}(col.gridcell),
                          exice_init_conc_col = zeros(Float64, bounds.endc))
    waterstatebulk_init_cold!(wsb, bc; h2osno_input_col = h2osno_input)

    # ---- WaterDiagnosticBulkType::InitCold ---------------------------------
    waterdiagnosticbulk_init_cold!(inst.water.waterdiagnosticbulk_inst, bc, bp;
                                   snow_depth_input_col = snow_depth_input,
                                   h2osno_input_col     = h2osno_input,
                                   snl_col              = Vector{Int}(col.snl),
                                   landunit_col         = Vector{Int}(col.landunit),
                                   urbpoi               = BitVector(lun.urbpoi))

    # ---- WaterFluxType / WaterFluxBulkType::InitCold ------------------------
    # NOT called here — PR #209 owns this one, via `init_water_flux_cold!` at the
    # top of `cold_start_initialize!` (it threads `use_hillslope_routing`, which
    # this entry point does not receive). Same sweep, different call site; the
    # wiring assertion in test_init_cold_wiring.jl covers both.

    # ---- WaterBalanceType::InitCold ----------------------------------------
    waterbalance_init_cold!(inst.water.waterbalancebulk_inst, bc)

    # ---- AerosolType::InitCold ---------------------------------------------
    # NOT called here — PR #210 owns this one, via `clm_instInit!` (right after
    # aerosol_init!), so the RESTART path is zeroed before the dump overwrites it
    # too, exactly as Fortran's Init does. Zeroing the mss_* snow-aerosol masses
    # is what stopped SNICAR NaN'ing the glacier column's albedo.

    # ---- CropType::InitCold ------------------------------------------------
    # `latbaset_patch` is deliberately left NaN unless baset_mapping == LATVARY;
    # crop_init_cold! encodes that. Safe to call with use_crop=false (the crop
    # fields are simply zeroed / flagged and never read).
    crop_init_cold!(inst.crop, bp;
                    patch_landunit = pch.landunit,
                    lun_itype      = lun.itype,
                    istcrop        = ISTCROP,
                    patch_gridcell = pch.gridcell,
                    patch_itype    = pch.itype,
                    pftcon_baset   = isempty(pftcon.baset) ? nothing : pftcon.baset,
                    grc_latdeg     = isempty(inst.gridcell.latdeg) ? nothing :
                                     Vector{Float64}(inst.gridcell.latdeg))

    return nothing
end

"""
    init_cold_biogeochem!(inst, bounds; nlevdecomp, nlevdecomp_full, ndecomp_pools)

Run every ported soil-BGC / CN-vegetation `*_init_cold!`.

Called from `clm_initialize!` AFTER the decomposition cascade is built, because
Fortran's `soilbiogeochem_carbonstate_inst%Init` (→ `InitCold`, which seeds
`decomp_cpools_vr` from `decomp_cascade_con%initial_stock` and the exponential
depth profile) runs AFTER `init_decompcascade_bgc` — clm_instMod.F90:414-425.

Called UNCONDITIONALLY, not gated on `use_cn`. Fortran gates the whole
`Init` (allocate + InitCold) on `decomp_method /= no_soil_decomp`, so a non-CN
Fortran run has these arrays UNALLOCATED. CLM.jl's `clm_instInit!` allocates
them unconditionally — so to preserve the Fortran invariant "every allocated
array has been InitCold'd", the InitCold must be unconditional too. With
`use_cn=false` the cascade's `initial_stock` is all-zero, so this seeds zeros
into arrays that are never read: inert, and strictly safer than NaN.
"""
function init_cold_biogeochem!(inst::CLMInstances, bounds::BoundsType;
                               nlevdecomp::Int = varpar.nlevdecomp,
                               nlevdecomp_full::Int = varpar.nlevdecomp_full,
                               ndecomp_pools::Int = 7)
    col = inst.column
    lun = inst.landunit
    pch = inst.patch
    bc  = bounds.begc:bounds.endc
    bp  = bounds.begp:bounds.endp

    cfg = inst.bgc_vegetation.config
    use_matrixcn      = getfield_or(cfg, :use_matrixcn, false)
    use_soil_matrixcn = getfield_or(cfg, :use_soil_matrixcn, false)
    use_crop          = getfield_or(cfg, :use_crop, false)
    use_nitrif        = getfield_or(cfg, :use_nitrif_denitrif, false)

    # Special (non-soil/crop) landunits: Fortran's InitCold zeroes their BGC
    # fluxes and skips the soil-column pool seeding.
    mask_special   = falses(bounds.endc)
    mask_soil_crop = falses(bounds.endc)
    for c in bc
        l = col.landunit[c]
        (l >= 1 && l <= length(lun.itype)) || continue
        it = lun.itype[l]
        if it == ISTSOIL || it == ISTCROP
            mask_soil_crop[c] = true
        else
            mask_special[c] = true
        end
    end

    cascade = inst.decomp_cascade
    initial_stock = isempty(cascade.initial_stock) ? zeros(ndecomp_pools) :
                    Vector{Float64}(cascade.initial_stock)
    initial_cn    = isempty(cascade.initial_cn_ratio) ? ones(ndecomp_pools) :
                    Vector{Float64}(cascade.initial_cn_ratio)
    zsoi_vals = length(zsoi[]) >= nlevdecomp ?
                Vector{Float64}(zsoi[][1:nlevdecomp]) : zeros(nlevdecomp)
    cstocks_depth = inst.decomp_bgc_params.bgc_initial_Cstocks_depth
    cstocks_depth = (cstocks_depth isa Real && cstocks_depth > 0) ? Float64(cstocks_depth) : 0.3

    # ---- SoilBiogeochemStateType::InitCold ---------------------------------
    soil_bgc_state_init_cold!(inst.soilbiogeochem_state, bc, nlevdecomp_full;
                              mask_special   = mask_special,
                              mask_soil_crop = mask_soil_crop)

    # ---- SoilBiogeochemCarbonStateType::InitCold (c12 + isotopes) ----------
    scs = inst.soilbiogeochem_carbonstate
    soil_bgc_carbon_state_init_cold!(scs, bc, 1.0;
                                     nlevdecomp = nlevdecomp,
                                     nlevdecomp_full = nlevdecomp_full,
                                     ndecomp_pools = ndecomp_pools,
                                     initial_stock = initial_stock,
                                     initial_stock_soildepth = cstocks_depth,
                                     zsoi_vals = zsoi_vals,
                                     use_soil_matrixcn = use_soil_matrixcn,
                                     mask_soil_crop = mask_soil_crop)
    for (iso_cs, ratio) in ((inst.c13_soilbiogeochem_carbonstate, C13RATIO),
                            (inst.c14_soilbiogeochem_carbonstate, C14RATIO))
        isempty(iso_cs.decomp_cpools_vr_col) && continue
        soil_bgc_carbon_state_init_cold!(iso_cs, bc, ratio;
                                         nlevdecomp = nlevdecomp,
                                         nlevdecomp_full = nlevdecomp_full,
                                         ndecomp_pools = ndecomp_pools,
                                         initial_stock = initial_stock,
                                         initial_stock_soildepth = cstocks_depth,
                                         zsoi_vals = zsoi_vals,
                                         use_soil_matrixcn = use_soil_matrixcn,
                                         mask_soil_crop = mask_soil_crop,
                                         c12_inst = scs)
    end

    # ---- SoilBiogeochemCarbonFluxType::InitCold ----------------------------
    soil_bgc_carbon_flux_init_cold!(inst.soilbiogeochem_carbonflux, bc;
                                    mask_special = mask_special)
    for iso_cf in (inst.c13_soilbiogeochem_carbonflux, inst.c14_soilbiogeochem_carbonflux)
        isempty(iso_cf.decomp_cpools_leached_col) && continue
        soil_bgc_carbon_flux_init_cold!(iso_cf, bc; mask_special = mask_special)
    end

    # ---- SoilBiogeochemNitrogenStateType::InitCold -------------------------
    # Seeds decomp_npools_vr from the C pools just seeded above / initial C:N.
    soil_bgc_nitrogen_state_init_cold!(inst.soilbiogeochem_nitrogenstate, bc;
                                       nlevdecomp = nlevdecomp,
                                       nlevdecomp_full = nlevdecomp_full,
                                       ndecomp_pools = ndecomp_pools,
                                       decomp_cpools_vr_col = Array(scs.decomp_cpools_vr_col),
                                       decomp_cpools_col    = Array(scs.decomp_cpools_col),
                                       decomp_cpools_1m_col = Array(scs.decomp_cpools_1m_col),
                                       initial_cn_ratio = initial_cn,
                                       use_soil_matrixcn = use_soil_matrixcn,
                                       use_nitrif_denitrif = use_nitrif,
                                       mask_soil_crop = mask_soil_crop)

    # ---- SoilBiogeochemNitrogenFluxType::InitCold --------------------------
    soil_bgc_nitrogen_flux_init_cold!(inst.soilbiogeochem_nitrogenflux, bc;
                                      mask_special = mask_special)

    # ---- CNVeg*Type::InitCold ----------------------------------------------
    # THE BIG ONE. On a use_cn=true COLD START the ENTIRE CN carbon/nitrogen
    # vegetation state was left at the allocator's NaN: nothing on the live init
    # path writes leafc / deadstemc / the N pools — only a restart or a
    # Fortran-dump injection does. That is why the CN parity harness never caught
    # it (it injects a dump) and why scripts/clmdrv_cn_fixture.jl had to hand-roll
    # these exact calls to get "a self-consistent finite initial state".
    #
    # `cn_vegetation_init!` sizes these arrays only when the corresponding
    # sub-facade is active, so each call is guarded on its target actually being
    # allocated (an unallocated facade has nothing to InitCold — Fortran's
    # allocate+InitCold pairing, preserved).
    veg = inst.bgc_vegetation
    isempty(veg.cnveg_state_inst.dormant_flag_patch) ||
        cnveg_state_init_cold!(veg.cnveg_state_inst, bc, bp;
                               col_landunit   = col.landunit,
                               patch_landunit = pch.landunit,
                               lun_ifspecial  = Vector{Bool}(lun.ifspecial),
                               lun_itype      = lun.itype)
    # Pass the PFT data so the CN cold start can take Fortran's per-PFT branches
    # (noveg / evergreen / crop / deciduous, and the woody deadstemc seed). Without it
    # every patch was seeded as deciduous: evergreen PFTs cold-started with ZERO
    # displayed leaf and root carbon, and bare-soil patches got storage C.
    # pftcon is filled by readParameters! before the cold start; if a caller has not
    # loaded a parameter file these stay empty, in which case fall back to the legacy
    # (all-deciduous) seeding rather than index an empty array.
    _pft_ok  = !isempty(pftcon.evergreen) && !isempty(pftcon.woody)
    _pft_itype = _pft_ok ? pch.itype : nothing
    isempty(veg.cnveg_carbonstate_inst.leafc_patch) ||
        cnveg_carbon_state_init_cold!(veg.cnveg_carbonstate_inst, bp;
                                      ratio = 1.0,
                                      use_matrixcn = use_matrixcn,
                                      use_crop = use_crop,
                                      patch_itype = _pft_itype,
                                      evergreen   = _pft_ok ? pftcon.evergreen : nothing,
                                      woody       = _pft_ok ? pftcon.woody : nothing)
    # N is DERIVED from the cold C pools via the PFT C:N ratios (Fortran does the same).
    _cn_ok = _pft_ok && !isempty(pftcon.leafcn) && !isempty(pftcon.frootcn) && !isempty(pftcon.deadwdcn)
    isempty(veg.cnveg_nitrogenstate_inst.leafn_patch) ||
        cnveg_nitrogen_state_init_cold!(veg.cnveg_nitrogenstate_inst, bp;
                                        use_matrixcn = use_matrixcn,
                                        use_crop = use_crop,
                                        carbonstate = _cn_ok ? veg.cnveg_carbonstate_inst : nothing,
                                        patch_itype = _cn_ok ? pch.itype : nothing,
                                        leafcn      = _cn_ok ? pftcon.leafcn : nothing,
                                        frootcn     = _cn_ok ? pftcon.frootcn : nothing,
                                        deadwdcn    = _cn_ok ? pftcon.deadwdcn : nothing,
                                        woody       = _cn_ok ? pftcon.woody : nothing)
    isempty(veg.cnveg_carbonflux_inst.gpp_patch) ||
        cnveg_carbon_flux_init_cold!(veg.cnveg_carbonflux_inst, bp, bc;
                                     use_matrixcn = use_matrixcn,
                                     nlevdecomp_full = nlevdecomp_full,
                                     ndecomp_pools = ndecomp_pools)
    isempty(veg.cnveg_nitrogenflux_inst.ndeploy_patch) ||
        cnveg_nitrogen_flux_init_cold!(veg.cnveg_nitrogenflux_inst, bp, bc;
                                       nlevdecomp_full = nlevdecomp_full,
                                       ndecomp_pools = ndecomp_pools)

    return nothing
end

# Small helper: read a config field if the struct has it, else a default. Keeps
# init_cold_biogeochem! decoupled from the exact CNVegetation config schema.
getfield_or(obj, name::Symbol, default) =
    hasproperty(obj, name) ? getproperty(obj, name) : default
