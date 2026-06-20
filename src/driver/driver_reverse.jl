# ==========================================================================
# Productionized clm_drv! reverse-mode phases for the compositional engine.
#
# `compositional_reverse!` (canopy_fluxes_reverse.jl) reverses an ordered list of
# `(phase_fn, const_args)` entries, one Enzyme.autodiff call each, on a checkpointed
# bundle. The canopy energy-balance/photosynthesis solve is exposed there as
# `cf_rev_phases`. This module exposes the subsequent clm_drv! timestep phases that
# reverse-differentiate cleanly — soil_temperature!, soil_water! (Zeng-Decker 2009),
# water_table! — as the SAME kind of entries, operating on a unified bundle
# `(; inst, scratch)` so a multi-module / whole-timestep reverse is one call:
#
#     b   = driver_rev_bundle(inst)
#     ph  = driver_rev_phases(bounds, filt, config; canopy_aux = caux, n_canopy = 12)
#     db  = compositional_reverse!(ph, b, (db,b) -> (db.inst.<field> .= 2 .* b.inst.<field>))
#
# Each wrapper mutates `b.inst` in place and calls the exact production physics
# function; the per-phase aux (filters/bounds/config flags/configs) is passed as a
# Const NamedTuple — NEVER closure-captured (Enzyme analyzes a closure's captured
# environment for activity and trips on nested structs/BitVectors).
#
# Validated (Julia 1.10 LTS + Enzyme, real Bow cold-start inst) by
# scripts/enzyme_driver_reverse_full.jl (canopy+soil_temp+soil_water) and
# scripts/enzyme_driver_reverse_hydro.jl (+ water_table!). The forward wrappers are
# guarded WITHOUT Enzyme in test/test_driver_reverse.jl (runs in the 1.12 suite).
# ==========================================================================

# Unified bundle: the differentiated inst + the canopy scratch arrays (so the canopy
# block and the hydrology phases share one checkpointed object).
driver_rev_bundle(inst) =
    (; inst, scratch = cf_rev_scratch(Float64, length(inst.patch.column)))

# --------------------------------------------------------------------------
# Canopy block — runs the N canopy sub-phases (cf_rev_phases) as ONE chain entry
# on the bundle's cf view. Reversed coarsely (one Enzyme call over the block) in the
# chain; for the tight per-sub-phase canopy gradient use canopy_rev_gradient!.
# --------------------------------------------------------------------------
function canopy_rev_block!(b, aux, n::Int)
    cv = cf_rev_bundle(b.inst.canopystate, b.inst.energyflux, b.inst.frictionvel,
        b.inst.temperature, b.inst.solarabs, b.inst.soilstate, b.inst.water.waterfluxbulk_inst,
        b.inst.water.waterstatebulk_inst, b.inst.water.waterdiagnosticbulk_inst,
        b.inst.photosyns, b.scratch)
    for (f, cargs) in cf_rev_phases(aux, n)
        f(cv, cargs...)
    end
    return nothing
end

# Build the canopy Const aux from a warmed-up inst (downscaled forcing must be set).
# Energy-balance path (use_psn=false) by default — the canopy water/energy gradient. The
# z0v roughness fills match clm_drv!'s "z0v keep core defaults" (clm_driver.jl), so this is
# faithful to the production canopy_fluxes_core! energy balance. Photosynthesis (use_psn=true)
# needs the real LUNA/alb arrays and is a follow-up.
function canopy_rev_aux(inst, bounds, filt; use_psn::Bool = false, dtime = 1800.0)
    FT = Float64; NP = bounds.endp; ev = filt.exposedvegp
    fp = Int[p for p in bounds.begp:bounds.endp if ev[p]]; fn = length(fp)
    a2l = inst.atm2lnd
    forc_q_col = fill!(similar(a2l.forc_pbot_downscaled_col), zero(FT))
    compute_forc_q!(forc_q_col, inst.column.gridcell, a2l.forc_vp_grc, a2l.forc_pbot_downscaled_col)
    forc = (; lwrad = a2l.forc_lwrad_downscaled_col, q = forc_q_col,
              pbot = a2l.forc_pbot_downscaled_col, th = a2l.forc_th_downscaled_col,
              rho = a2l.forc_rho_downscaled_col, t = a2l.forc_t_downscaled_col,
              u_grc = a2l.forc_u_grc, v_grc = a2l.forc_v_grc, pco2 = a2l.forc_pco2_grc, po2 = a2l.forc_po2_grc,
              hgt_t = a2l.forc_hgt_t_grc, hgt_u = a2l.forc_hgt_u_grc, hgt_q = a2l.forc_hgt_q_grc,
              dayl = inst.gridcell.dayl, max_dayl = inst.gridcell.max_dayl,
              downreg = fill(FT(1.0), NP), leafn = fill(FT(1.0), NP))
    MP = MXPFT + 1
    pft = (; dleaf = fill(FT(0.04), MP), z0v_Cr = fill(FT(0.35), MP), z0v_Cs = fill(FT(0.003), MP),
             z0v_c = fill(FT(0.25), MP), z0v_cw = fill(FT(2.0), MP), z0v_LAImax = fill(FT(8.0), MP),
             grnd_ch4 = fill(FT(0.0), NP))
    # Photosynthesis aux. With use_psn=true the psn sub-phases run, so source the REAL
    # canopy radiation / LUNA arrays from the warmed-up inst (surfalb/solarabs/canopystate/
    # pftcon) — but as Const COPIES, never bundle references, or Enzyme sees the same array
    # as both Const (here) and Duplicated (in the differentiated bundle). photosynthesis!'s
    # scalar params come from the module-global params_inst (set during init/warmup). With
    # use_psn=false the psn sub-phases are skipped, so cheap synthetic fills suffice.
    psn = if use_psn
        sa = inst.surfalb; so = inst.solarabs; cs2 = inst.canopystate
        (; c3psn = copy(pftcon.c3psn), leafcn = copy(pftcon.leafcn), flnr = copy(pftcon.flnr),
           fnitr = copy(pftcon.fnitr), slatop = copy(pftcon.slatop), mbbopt = copy(pftcon.mbbopt),
           medlynintercept = fill(FT(100.0), MP), medlynslope = fill(FT(6.0), MP),
           nrad = copy(sa.nrad_patch), tlai_z = copy(sa.tlai_z_patch),
           parsun_z = copy(so.parsun_z_patch), parsha_z = copy(so.parsha_z_patch),
           laisun_z = copy(cs2.laisun_z_patch), laisha_z = copy(cs2.laisha_z_patch),
           vcmaxcint_sun = copy(sa.vcmaxcintsun_patch), vcmaxcint_sha = copy(sa.vcmaxcintsha_patch),
           o3coefv = fill(FT(1.0), NP), o3coefg = fill(FT(1.0), NP),
           t10 = copy(inst.temperature.t_a10_patch))
    else
        (; c3psn = fill(FT(1.0), MP), leafcn = fill(FT(25.0), MP), flnr = fill(FT(0.1), MP),
           fnitr = fill(FT(1.0), MP), slatop = fill(FT(0.01), MP), mbbopt = fill(FT(9.0), MP),
           medlynintercept = fill(FT(100.0), MP), medlynslope = fill(FT(6.0), MP),
           nrad = fill(1, NP), tlai_z = fill(FT(1.0), NP, NLEVCAN),
           parsun_z = fill(FT(0.0), NP, NLEVCAN), parsha_z = fill(FT(0.0), NP, NLEVCAN),
           laisun_z = fill(FT(0.0), NP, NLEVCAN), laisha_z = fill(FT(0.0), NP, NLEVCAN),
           vcmaxcint_sun = fill(FT(1.0), NP), vcmaxcint_sha = fill(FT(0.6), NP),
           o3coefv = fill(FT(1.0), NP), o3coefg = fill(FT(1.0), NP),
           t10 = copy(inst.temperature.t_a10_patch))
    end
    return (; patch = inst.patch, col = inst.column, grid = inst.gridcell, forc = forc, pft = pft,
        psn = psn, filterp = fp, fn = fn, active = Bool[ev[p] for p in 1:NP], mask = ev,
        ivt = inst.patch.itype .+ 1,
        forc_pbot_patch = FT[a2l.forc_pbot_downscaled_col[inst.patch.column[p]] for p in 1:NP],
        soilevap_beta = do_soilevap_beta(), soil_resis_sl14 = do_soil_resistance_sl14(),
        nlevsno = varpar.nlevsno, dtime = FT(dtime), use_psn = use_psn)
end

# --------------------------------------------------------------------------
# soil_temperature!
# --------------------------------------------------------------------------
soiltemp_rev_aux(bounds, filt; dtime = 1800.0) = (;
    urbantv = fill(323.15, bounds.endl), dtime = dtime,
    bc_col = bounds.begc:bounds.endc, bc_lun = bounds.begl:bounds.endl,
    bc_patch = bounds.begp:bounds.endp,
    nolakec = filt.nolakec, nolakep = filt.nolakep, urbanl = filt.urbanl, urbanc = filt.urbanc)

function soiltemp_rev_phase!(b, aux)
    i = b.inst
    soil_temperature!(i.column, i.landunit, i.patch, i.temperature, i.energyflux,
        i.soilstate, i.water.waterstatebulk_inst, i.water.waterdiagnosticbulk_inst,
        i.water.waterfluxbulk_inst, i.solarabs, i.canopystate, i.urbanparams,
        aux.urbantv, i.atm2lnd.forc_lwrad_downscaled_col, aux.nolakec, aux.nolakep,
        aux.urbanl, aux.urbanc, aux.bc_col, aux.bc_lun, aux.bc_patch, aux.dtime)
    return nothing
end

# --------------------------------------------------------------------------
# soil_water! (Zeng-Decker 2009 default; the retention curve is an empty singleton
# and the movement config a small mutable struct — both pass cleanly as Const).
# --------------------------------------------------------------------------
function soilwater_rev_aux(filt, config; dtime = 1800.0)
    cfg = config.use_aquifer_layer ? SoilWaterMovementConfig() :
        SoilWaterMovementConfig(soilwater_movement_method = MOISTURE_FORM,
                                lower_boundary_condition = BC_ZERO_FLUX)
    return (; hydrologyc = filt.hydrologyc, urbanc = filt.urbanc,
              swrc = SoilWaterRetentionCurveClappHornberg1978(), cfg = cfg, dtime = dtime)
end

function soilwater_rev_phase!(b, aux)
    i = b.inst
    soil_water!(i.column, aux.hydrologyc, aux.urbanc,
        i.soilhydrology, i.soilstate, i.water.waterfluxbulk_inst, i.water.waterstatebulk_inst,
        i.temperature, i.canopystate, i.energyflux, aux.swrc, aux.cfg, aux.dtime)
    return nothing
end

# --------------------------------------------------------------------------
# water_table! (aquifer / ThetaBasedWaterTable). The single Bool kwarg is resolved
# inside the wrapper (primal-time), so Enzyme differentiates the body without a
# kwarg-NamedTuple thunk.
# --------------------------------------------------------------------------
watertable_rev_aux(bounds, filt, config; dtime = 1800.0) = (;
    hydrologyc = filt.hydrologyc, bc_col = bounds.begc:bounds.endc,
    nlevsoi = varpar.nlevsoi, dtime = dtime,
    recompute_frost_table = config.use_aquifer_layer)

function watertable_rev_phase!(b, aux)
    i = b.inst
    water_table!(i.soilhydrology, i.soilstate,
        i.temperature.t_soisno_col, i.water.waterstatebulk_inst.ws, i.water.waterfluxbulk_inst,
        i.column.dz, i.column.z, i.column.zi,
        aux.hydrologyc, aux.bc_col, aux.nlevsoi, aux.dtime;
        recompute_frost_table = aux.recompute_frost_table)
    return nothing
end

# --------------------------------------------------------------------------
# hydrology_no_drainage! — the namesake HydrologyNoDrainage end-of-step
# diagnostics (recomputes h2osoi_vol, the 10 cm / total water columns, etc. from
# the updated h2osoi_liq/ice). Runs after water_table! (the intervening snow
# capping/compaction/layer combine-divide ops are discrete and no-ops at snl=0).
# --------------------------------------------------------------------------
hydnodrain_rev_aux(bounds, filt; dtime = 1800.0) = (;
    nolakec = filt.nolakec, hydrologyc = filt.hydrologyc, urbanc = filt.urbanc,
    snowc = filt.snowc, nosnowc = filt.nosnowc, bc_col = bounds.begc:bounds.endc,
    dtime = dtime, nlevsno = varpar.nlevsno, nlevsoi = varpar.nlevsoi,
    nlevgrnd = varpar.nlevgrnd, nlevurb = varpar.nlevurb)

function hydnodrain_rev_phase!(b, aux)
    i = b.inst
    hydrology_no_drainage!(i.temperature, i.soilstate, i.water.waterstatebulk_inst,
        i.water.waterdiagnosticbulk_inst, i.column, i.landunit,
        aux.nolakec, aux.hydrologyc, aux.urbanc, aux.snowc, aux.nosnowc,
        aux.bc_col, aux.dtime, aux.nlevsno, aux.nlevsoi, aux.nlevgrnd, aux.nlevurb)
    return nothing
end

# --------------------------------------------------------------------------
# Pre-soil_water! SURFACE-HYDROLOGY block (HydrologyNoDrainage). The full ordered
# sequence that runs AFTER soil_temperature! and BEFORE soil_water!, partitioning the
# surface water that becomes soil_water!'s top boundary:
#   set_soil_water_fractions! → set_floodc! → saturated_excess_runoff! → set_qflx_inputs!
#   → infiltration_excess_runoff! → route_infiltration_excess! → update_h2osfc!
#   → infiltration! → total_surface_runoff! → compute_effec_rootfrac_and_vert_tran_sink!
# All share one Const aux. saturated_excess_runoff!'s Bool kwargs + update_h2osfc!'s
# dtime/h2osfcflag + compute_effec_rootfrac's use_hydrstress default inside the wrappers;
# qflx_floodg is a zero gridcell vector (no river flood input), captured as Const.
function surfhydro_rev_aux(bounds, filt; dtime = 1800.0)
    return (; hydrologyc = filt.hydrologyc, nolakec = filt.nolakec, urbanc = filt.urbanc,
              bc_col = bounds.begc:bounds.endc, nlevsoi = varpar.nlevsoi, nlevsno = varpar.nlevsno,
              dtime = dtime, qflx_floodg = zeros(Float64, bounds.endg))
end

function setsoilfrac_rev_phase!(b, aux)
    i = b.inst
    set_soil_water_fractions!(i.soilhydrology, i.soilstate, i.water.waterstatebulk_inst.ws,
        i.column.dz, aux.hydrologyc, aux.bc_col, aux.nlevsoi, aux.nlevsno)
    return nothing
end
function setfloodc_rev_phase!(b, aux)
    i = b.inst
    set_floodc!(i.water.waterfluxbulk_inst.wf.qflx_floodc_col, aux.qflx_floodg,
        i.column.gridcell, i.column.itype, aux.nolakec, aux.bc_col)
    return nothing
end
function satexcess_rev_phase!(b, aux)
    i = b.inst
    saturated_excess_runoff!(i.sat_excess_runoff, aux.hydrologyc, aux.bc_col,
        i.column, i.landunit, i.soilhydrology, i.soilstate, i.water.waterfluxbulk_inst)
    return nothing
end
function setqflx_rev_phase!(b, aux)
    i = b.inst
    set_qflx_inputs!(i.water.waterfluxbulk_inst, i.water.waterdiagnosticbulk_inst,
        i.column.snl, aux.hydrologyc, aux.bc_col)
    return nothing
end
function inflexcess_rev_phase!(b, aux)
    i = b.inst
    infiltration_excess_runoff!(i.infilt_excess_runoff, i.soilhydrology, i.soilstate,
        i.sat_excess_runoff.fsat_col, i.water.waterfluxbulk_inst, i.water.waterdiagnosticbulk_inst,
        aux.hydrologyc, aux.bc_col)
    return nothing
end
function routeinfl_rev_phase!(b, aux)
    i = b.inst
    route_infiltration_excess!(i.water.waterfluxbulk_inst, i.soilhydrology,
        i.column.landunit, i.landunit.itype, aux.hydrologyc, aux.bc_col)
    return nothing
end
function updateh2osfc_rev_phase!(b, aux)
    i = b.inst
    update_h2osfc!(i.column, i.soilhydrology, i.energyflux, i.water.waterfluxbulk_inst,
        i.water.waterstatebulk_inst, i.water.waterdiagnosticbulk_inst,
        aux.hydrologyc, aux.bc_col; dtime = aux.dtime)
    return nothing
end
function infil_rev_phase!(b, aux)
    infiltration!(b.inst.water.waterfluxbulk_inst, aux.hydrologyc, aux.bc_col)
    return nothing
end
function totalrunoff_rev_phase!(b, aux)
    i = b.inst
    total_surface_runoff!(i.water.waterfluxbulk_inst, i.soilhydrology, i.water.waterstatebulk_inst.ws,
        i.column.snl, i.column.itype, i.column.landunit, i.landunit.urbpoi,
        aux.hydrologyc, aux.urbanc, aux.bc_col, aux.dtime)
    return nothing
end
function rootsink_rev_phase!(b, aux)
    i = b.inst
    compute_effec_rootfrac_and_vert_tran_sink!(aux.bc_col, aux.nlevsoi, aux.hydrologyc,
        i.soilstate, i.canopystate, i.water.waterfluxbulk_inst, i.energyflux,
        i.column, i.landunit, i.patch)
    return nothing
end

# Ordered pre-soil_water surface-hydrology phase entries.
function surface_hydrology_rev_phases(bounds, filt; dtime = 1800.0)
    a = surfhydro_rev_aux(bounds, filt; dtime)
    return Any[(setsoilfrac_rev_phase!, (a,)), (setfloodc_rev_phase!, (a,)),
               (satexcess_rev_phase!, (a,)),   (setqflx_rev_phase!, (a,)),
               (inflexcess_rev_phase!, (a,)),  (routeinfl_rev_phase!, (a,)),
               (updateh2osfc_rev_phase!, (a,)), (infil_rev_phase!, (a,)),
               (totalrunoff_rev_phase!, (a,)),  (rootsink_rev_phase!, (a,))]
end

# --------------------------------------------------------------------------
# BGC (use_cn) carbon/nitrogen phases. FIRST entry into the biogeochemistry domain:
# cn_gresp! (growth respiration, cpool_*_gr = grperc·grpnow·allocation — smooth/linear)
# reverse-differentiates at machine precision (scripts/enzyme_bgc_reverse.jl), proving the
# CN carbon-flux cascade is reverse-able via the same template. Operates on a use_cn inst's
# bgc_vegetation.cnveg_carbonflux_inst; the PFT growth params come from the global pftcon.
# The full BGC chain (allocation, decomposition, phenology, mortality, fire, N-cycling) is a
# large separate domain — several of its sub-phases (phenology onset/offset, fire, gap
# mortality) are genuinely non-differentiable discrete state machines — so this is exposed as
# a building block, not auto-inserted into driver_rev_phases.
function cngresp_rev_aux(inst, bounds, filt; npcropmin::Int = NPCROPMIN, nrepr::Int = NREPR)
    gr = PftConGrowthResp(woody = Float64.(pftcon.woody),
                          grperc = Float64.(pftcon.grperc), grpnow = Float64.(pftcon.grpnow))
    return (; mask = filt.bgc_vegp, bounds = bounds.begp:bounds.endp, pftcon = gr,
              patch = inst.patch, npcropmin = npcropmin, nrepr = nrepr)
end
function cngresp_rev_phase!(b, aux)
    cn_gresp!(aux.mask, aux.bounds, aux.pftcon, aux.patch,
              b.inst.bgc_vegetation.cnveg_carbonflux_inst;
              npcropmin = aux.npcropmin, nrepr = aux.nrepr)
    return nothing
end

# cn_mresp! — maintenance respiration (the `br`/`br_root` calibration params; e.g.
# livestem_mr = livestemn·br·tc(Q10)). Reverse-validated at machine precision
# (scripts/enzyme_bgc_reverse.jl [P2]). Differentiated state: cnveg_carbonflux_inst (out)
# + cnveg_nitrogenstate_inst (N-pool inputs); the params (MaintRespParams from
# maint_resp_read_params! w/ br_root from cn_shared_params, PftConMaintResp from pftcon) +
# cn_shared_params are Const aux. canopy/soil/temperature/photosyns flow from the bundle.
function cnmresp_rev_aux(inst, bounds, filt; npcropmin::Int = NPCROPMIN, nrepr::Int = NREPR)
    mrp = MaintRespParams(); maint_resp_read_params!(mrp; br_root = inst.cn_shared_params.br_root)
    pftmr = PftConMaintResp{Float64}(woody = Float64.(pftcon.woody))
    return (; mask_c = filt.bgc_soilc, mask_p = filt.bgc_vegp,
              bounds_c = bounds.begc:bounds.endc, bounds_p = bounds.begp:bounds.endp,
              params = mrp, cn_params = inst.cn_shared_params, pftcon = pftmr,
              nlevgrnd = varpar.nlevgrnd, nlevsno = varpar.nlevsno, npcropmin = npcropmin, nrepr = nrepr)
end
function cnmresp_rev_phase!(b, aux)
    i = b.inst
    cn_mresp!(aux.mask_c, aux.mask_p, aux.bounds_c, aux.bounds_p, aux.params, aux.cn_params,
              aux.pftcon, i.patch, i.canopystate, i.soilstate, i.temperature, i.photosyns,
              i.bgc_vegetation.cnveg_carbonflux_inst, i.bgc_vegetation.cnveg_nitrogenstate_inst;
              nlevgrnd = aux.nlevgrnd, nlevsno = aux.nlevsno, npcropmin = aux.npcropmin, nrepr = aux.nrepr)
    return nothing
end

# c_state_update1! — the carbon-cycle pool integration (e.g. leafc += cpool_to_leafc·dt;
# cpool drained by allocation/MR/growth-resp). Reverse-validated at machine precision
# (scripts/enzyme_bgc_reverse.jl [P3]). Differentiated state: cnveg_carbonstate_inst (pools,
# out) ← cnveg_carbonflux_inst (fluxes, in) + soilbiogeochem_carbonflux (litter sink). Every
# index/cascade arg (patch_column, ivt, woody, cascade_donor/receiver_pool, harvdate,
# col_is_fates, i_litr_*, …) is a Const, sourced from inst/config exactly as cn_driver does.
function cstate1_rev_aux(inst, bounds, filt, config)
    return (; mask_c = filt.bgc_soilc, mask_p = filt.bgc_vegp,
              bounds_c = bounds.begc:bounds.endc, bounds_p = bounds.begp:bounds.endp,
              patch_column = inst.patch.column, ivt = inst.patch.itype,
              woody = Float64.(pftcon.woody),
              cascade_donor_pool = inst.decomp_cascade.cascade_donor_pool,
              cascade_receiver_pool = inst.decomp_cascade.cascade_receiver_pool,
              harvdate = inst.crop.harvdate_patch, col_is_fates = inst.column.is_fates,
              nlevdecomp = varpar.nlevdecomp,
              ndecomp_cascade_transitions = config.ndecomp_cascade_transitions,
              i_litr_min = config.i_litr_min, i_litr_max = config.i_litr_max, i_cwd = config.i_cwd,
              npcropmin = config.npcropmin, nrepr = config.nrepr, dt = 1800.0)
end
function cstate1_rev_phase!(b, aux)
    i = b.inst
    c_state_update1!(i.bgc_vegetation.cnveg_carbonstate_inst,
                     i.bgc_vegetation.cnveg_carbonflux_inst, i.soilbiogeochem_carbonflux;
                     mask_soilc = aux.mask_c, mask_soilp = aux.mask_p,
                     bounds_col = aux.bounds_c, bounds_patch = aux.bounds_p,
                     patch_column = aux.patch_column, ivt = aux.ivt, woody = aux.woody,
                     cascade_donor_pool = aux.cascade_donor_pool,
                     cascade_receiver_pool = aux.cascade_receiver_pool,
                     harvdate = aux.harvdate, col_is_fates = aux.col_is_fates,
                     nlevdecomp = aux.nlevdecomp,
                     ndecomp_cascade_transitions = aux.ndecomp_cascade_transitions,
                     i_litr_min = aux.i_litr_min, i_litr_max = aux.i_litr_max, i_cwd = aux.i_cwd,
                     npcropmin = aux.npcropmin, nrepr = aux.nrepr, dt = aux.dt)
    return nothing
end

# n_state_update1! — the nitrogen-cycle pool integration (leafn += npool_to_leafn·dt; the
# N counterpart of c_state_update1!, fewer consts — no cascade-pool/harvdate). Reverse-
# validated at machine precision (scripts/enzyme_bgc_reverse.jl [P4]). Differentiated state:
# cnveg_nitrogenstate_inst (N pools, out) ← cnveg_nitrogenflux_inst (N fluxes, in) +
# soilbiogeochem_nitrogenflux (litter-N sink). Consts sourced from inst/config as cn_driver.
function nstate1_rev_aux(inst, bounds, filt, config)
    return (; mask_c = filt.bgc_soilc, mask_p = filt.bgc_vegp,
              bounds_c = bounds.begc:bounds.endc, bounds_p = bounds.begp:bounds.endp,
              ivt = inst.patch.itype, woody = Float64.(pftcon.woody),
              col_is_fates = inst.column.is_fates, nlevdecomp = varpar.nlevdecomp,
              i_litr_min = config.i_litr_min, i_litr_max = config.i_litr_max, i_cwd = config.i_cwd,
              npcropmin = config.npcropmin, nrepr = config.nrepr, use_fun = config.use_fun, dt = 1800.0)
end
function nstate1_rev_phase!(b, aux)
    i = b.inst
    n_state_update1!(i.bgc_vegetation.cnveg_nitrogenstate_inst,
                     i.bgc_vegetation.cnveg_nitrogenflux_inst, i.soilbiogeochem_nitrogenflux;
                     mask_soilc = aux.mask_c, mask_soilp = aux.mask_p,
                     bounds_col = aux.bounds_c, bounds_patch = aux.bounds_p,
                     ivt = aux.ivt, woody = aux.woody, col_is_fates = aux.col_is_fates,
                     nlevdecomp = aux.nlevdecomp, i_litr_min = aux.i_litr_min,
                     i_litr_max = aux.i_litr_max, i_cwd = aux.i_cwd,
                     npcropmin = aux.npcropmin, nrepr = aux.nrepr, use_fun = aux.use_fun, dt = aux.dt)
    return nothing
end

# decomp_rate_constants_bgc! — soil-carbon decomposition rate constants (the Q10/tau soil-C
# turnover calibration: t_scalar = Q10^((Tsoi−Tref)/10), w_scalar(soilpsi), decomp_k). Reverse-
# validated at machine precision (scripts/enzyme_bgc_reverse.jl [P5]). The DIFFERENTIABLE
# inputs flow LIVE from the bundle: t_soisno (a materialized soil-layer slice of
# inst.temperature.t_soisno_col — a getindex copy Enzyme tracks back to the array) and
# soilpsi (inst.soilstate.soilpsi_col); output is inst.soilbiogeochem_carbonflux (t_scalar/
# w_scalar/decomp_k). params/state/cascade + geometry (zsoi[], col_dz slice) are Const aux,
# sourced from inst/globals exactly as cn_driver_no_leaching! does.
function decomprate_rev_aux(inst, bounds, filt)
    nsno = varpar.nlevsno
    return (; mask = filt.bgc_soilc, bounds = bounds.begc:bounds.endc,
              params = inst.decomp_bgc_params, bgc_state = inst.decomp_bgc_state,
              cn_params = inst.cn_shared_params, cascade_con = inst.decomp_cascade,
              nlevdecomp = varpar.nlevdecomp, nlevsno = nsno, zsoi_vals = zsoi[],
              col_dz = inst.column.dz[:, (nsno + 1):end], days_per_year = 365.0, dt = 1800.0)
end
function decomprate_rev_phase!(b, aux)
    i = b.inst
    t_soisno = i.temperature.t_soisno_col[:, (aux.nlevsno + 1):end]   # differentiable slice copy
    decomp_rate_constants_bgc!(i.soilbiogeochem_carbonflux, aux.bgc_state, aux.params,
                               aux.cn_params, aux.cascade_con;
                               mask_bgc_soilc = aux.mask, bounds = aux.bounds,
                               nlevdecomp = aux.nlevdecomp, t_soisno = t_soisno,
                               soilpsi = i.soilstate.soilpsi_col, days_per_year = aux.days_per_year,
                               dt = aux.dt, zsoi_vals = aux.zsoi_vals, col_dz = aux.col_dz)
    return nothing
end

# --------------------------------------------------------------------------
# Assembler: the ordered (phase_fn, const_args) list for a clm_drv! reverse, in
# forward order: [canopy] → soil_temp → <surface hydrology block> → soil_water →
# water_table → hydrology_no_drainage. `canopy_aux === nothing` skips the canopy block.
# `include_surface=false` reverts to the soil_temp→soil_water direct jump.
# --------------------------------------------------------------------------
function driver_rev_phases(bounds, filt, config; canopy_aux = nothing,
                           n_canopy::Int = 12, dtime = 1800.0, include_surface::Bool = true)
    phases = Any[]
    canopy_aux === nothing || push!(phases, (canopy_rev_block!, (canopy_aux, n_canopy)))
    push!(phases, (soiltemp_rev_phase!, (soiltemp_rev_aux(bounds, filt; dtime),)))
    include_surface && append!(phases, surface_hydrology_rev_phases(bounds, filt; dtime))
    push!(phases, (soilwater_rev_phase!,  (soilwater_rev_aux(filt, config; dtime),)))
    push!(phases, (watertable_rev_phase!, (watertable_rev_aux(bounds, filt, config; dtime),)))
    push!(phases, (hydnodrain_rev_phase!, (hydnodrain_rev_aux(bounds, filt; dtime),)))
    return phases
end
