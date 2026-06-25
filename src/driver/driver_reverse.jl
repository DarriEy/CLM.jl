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

# Robust config-field accessor. CLMDriverConfig forwards most CN settings via a custom
# getproperty (ndecomp_pools/i_litr_*/npcropmin/nrepr/use_c13/…) but NOT use_fun /
# use_nitrif_denitrif; CNDriverConfig has those instead. `_cfgget` returns the field if any
# accessor path exposes it, else the default — so the BGC aux builders work with either config.
_cfgget(config, f::Symbol, default) = try getproperty(config, f) catch; default end

# BGC bundle: the differentiated inst + the shared decomposition SCRATCH arrays
# (cn_decomp_pools / p_decomp_cpool_loss / p_decomp_cn_gain / pmnf_decomp_cascade /
# p_decomp_npool_to_din — the cn_driver `_z3` work arrays). These MUST be part of the
# Duplicated bundle, NOT Const aux: soil_bgc_potential! computes potential_immob THROUGH them
# (decomp_k → p_decomp_cpool_loss → pmnf → potential_immob), and a Const scratch array cuts
# that derivative (Enzyme propagates no gradient through a Const). Sharing them across the
# decomp phases (potential writes, competition + decomp read) also matches cn_driver, where
# the same _z3 locals thread the cascade. Sizes from config; mirrors driver_rev_bundle.
function bgc_rev_bundle(inst, bounds, config)
    nc = length(bounds.begc:bounds.endc); nld = varpar.nlevdecomp
    ndp = config.ndecomp_pools; nct = config.ndecomp_cascade_transitions
    z3(d) = zeros(Float64, nc, nld, d)
    bgc = (; cn_decomp_pools = z3(ndp), p_decomp_cpool_loss = z3(nct), p_decomp_cn_gain = z3(ndp),
             pmnf_decomp_cascade = z3(nct), p_decomp_npool_to_din = z3(nct))
    return (; inst, scratch = cf_rev_scratch(Float64, length(inst.patch.column)), bgc)
end

# --------------------------------------------------------------------------
# Canopy block — runs the N canopy sub-phases (cf_rev_phases) as ONE chain entry
# on the bundle's cf view. Reversed coarsely (one Enzyme call over the block) in the
# chain; for the tight per-sub-phase canopy gradient use canopy_rev_gradient!.
# --------------------------------------------------------------------------
# A canopy `cf_rev_bundle` view onto the driver bundle's inst + shared canopy scratch —
# the object the canopy sub-phases (cf_rev_*) operate on. Aliases the inst's arrays, so
# mutating the view mutates the driver bundle in place.
_canopy_view(b) = cf_rev_bundle(b.inst.canopystate, b.inst.energyflux, b.inst.frictionvel,
    b.inst.temperature, b.inst.solarabs, b.inst.soilstate, b.inst.water.waterfluxbulk_inst,
    b.inst.water.waterstatebulk_inst, b.inst.water.waterdiagnosticbulk_inst,
    b.inst.photosyns, b.scratch)

function canopy_rev_block!(b, aux, n::Int)
    cv = _canopy_view(b)
    for (f, cargs) in cf_rev_phases(aux, n)
        f(cv, cargs...)
    end
    return nothing
end

# Convergence-aware canopy iteration count for a DRIVER bundle: builds the canopy view
# and defers to canopy_rev_converged_n (the per-patch convergence test). Used by
# clm_drv_reverse! to auto-size the canopy block instead of hard-coding n_canopy.
driver_canopy_converged_n(b, canopy_aux) = canopy_rev_converged_n(_canopy_view(b), canopy_aux)

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
              npcropmin = config.npcropmin, nrepr = config.nrepr, use_fun = _cfgget(config, :use_fun, false), dt = 1800.0)
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

# soil_bgc_potential! — potential decomposition + mineral-N flux (the decomposition cascade
# step that consumes decomp_k from decomp_rate_constants_bgc!): p_decomp_cpool_loss[tr] =
# Cpool[donor]·decomp_k[donor]·pathfrac, gross mineralization / potential immobilization,
# floating-vs-fixed C:N. Reverse-validated at machine precision (scripts/enzyme_bgc_reverse.jl
# [P6]). All differentiable C/N structs flow LIVE from the bundle (cf decomp_k in / phr_vr +
# nf potential_immob/gross_nmin out, cs decomp_cpools / ns decomp_npools in); cascade_con +
# the five potential-flux SCRATCH arrays (cn_decomp_pools/p_decomp_*/pmnf — preallocated,
# overwritten each call) are Const aux, sized exactly as cn_driver_no_leaching!'s _z3.
function decomppot_rev_aux(inst, bounds, filt, config)
    return (; cascade_con = inst.decomp_cascade, mask = filt.bgc_soilc, bounds = bounds.begc:bounds.endc,
              nlevdecomp = varpar.nlevdecomp, ndecomp_pools = config.ndecomp_pools,
              ndecomp_cascade_transitions = config.ndecomp_cascade_transitions)
end
function decomppot_rev_phase!(b, aux)
    i = b.inst; s = b.bgc
    soil_bgc_potential!(i.soilbiogeochem_carbonflux, i.soilbiogeochem_carbonstate,
                        i.soilbiogeochem_nitrogenflux, i.soilbiogeochem_nitrogenstate,
                        i.soilbiogeochem_state, aux.cascade_con;
                        mask_bgc_soilc = aux.mask, bounds = aux.bounds, nlevdecomp = aux.nlevdecomp,
                        ndecomp_pools = aux.ndecomp_pools,
                        ndecomp_cascade_transitions = aux.ndecomp_cascade_transitions,
                        cn_decomp_pools = s.cn_decomp_pools,
                        p_decomp_cpool_loss = s.p_decomp_cpool_loss,
                        p_decomp_cn_gain = s.p_decomp_cn_gain,
                        pmnf_decomp_cascade = s.pmnf_decomp_cascade,
                        p_decomp_npool_to_din = s.p_decomp_npool_to_din)
    return nothing
end

# soilbiogeochem_n_state_update1! — the MINERAL-N pool update (applies N deposition+fixation,
# gross mineralization, immobilization, plant uptake, nitrification/denitrification and the
# carbon-only-limitation supplement to smin_nh4_vr / smin_no3_vr; sminn = nh4+no3). The
# downstream sink of the decomposition cascade's N side. Reverse-validated at machine
# precision (scripts/enzyme_bgc_reverse.jl [P7]). All three structs flow LIVE from the bundle
# (ns smin_nh4/no3/sminn out ← nf the mineral-N fluxes in; st ndep/nfixation profiles, read-
# only). use_fun resolved inside the wrapper. Sourced from inst/config as cn_driver does.
function sminnupdate_rev_aux(inst, bounds, filt, config)
    return (; mask = filt.bgc_soilc, bounds = bounds.begc:bounds.endc,
              nlevdecomp = varpar.nlevdecomp, use_fun = _cfgget(config, :use_fun, false), dt = 1800.0)
end
function sminnupdate_rev_phase!(b, aux)
    i = b.inst
    soilbiogeochem_n_state_update1!(i.soilbiogeochem_nitrogenstate, i.soilbiogeochem_nitrogenflux,
                                    i.soilbiogeochem_state;
                                    mask_bgc_soilc = aux.mask, bounds_col = aux.bounds,
                                    nlevdecomp = aux.nlevdecomp, dt = aux.dt, use_fun = aux.use_fun)
    return nothing
end

# cn_gap_mortality! — background gap mortality (m_{leaf,froot,stem,…}c_to_litter = pool · r_mort
# rate). A CONTINUOUS/linear phase: its branches are on the STATIC `woody` PFT flag and `>0`
# true-zeros, so it reverse-differentiates at machine precision with NO smoothing needed
# (scripts/enzyme_gapmort.jl). Operates on bgc_vegetation.cnveg_carbon/nitrogen{state,flux}_inst;
# params/pftcon/dgvs sourced from the global pftcon + zeros exactly as cn_driver (DGVS path off).
function gapmort_rev_aux(inst, bounds, filt, config)
    np = length(bounds.begp:bounds.endp)
    gmp = GapMortalityParams{Float64}(k_mort = 0.3, r_mort = fill(0.02, length(pftcon.woody)))
    pftgm = PftConGapMort{Float64}(woody = Float64.(pftcon.woody), leafcn = Float64.(pftcon.leafcn),
        livewdcn = Float64.(pftcon.livewdcn), lf_f = Float64.(pftcon.lf_f), fr_f = Float64.(pftcon.fr_f))
    dgvs = DgvsGapMortData{Float64}(greffic_patch = zeros(np), heatstress_patch = zeros(np), nind_patch = zeros(np))
    return (; params = gmp, pftcon = pftgm, dgvs = dgvs, mask = filt.bgc_vegp,
              bounds = bounds.begp:bounds.endp, use_matrixcn = config.use_matrixcn,
              npcropmin = config.npcropmin, dt = 1800.0)
end
function gapmort_rev_phase!(b, aux)
    i = b.inst; veg = i.bgc_vegetation
    cn_gap_mortality!(aux.mask, aux.bounds, aux.params, aux.pftcon, aux.dgvs, i.patch, i.canopystate,
                      veg.cnveg_carbonstate_inst, veg.cnveg_carbonflux_inst,
                      veg.cnveg_nitrogenstate_inst, veg.cnveg_nitrogenflux_inst;
                      dt = aux.dt, days_per_year = 365.0, use_cndv = false,
                      use_matrixcn = aux.use_matrixcn, spinup_state = 0, npcropmin = aux.npcropmin)
    return nothing
end

# calc_plant_cn_alloc! — plant C/N ALLOCATION (after nutrient competition): distributes the
# available C/N to leaf/froot/stem/croot/repr pools via the allometry + downregulation logic
# (sminn_to_npool = plant_ndemand·fpg → plant_nalloc → plant_calloc → cpool_to_*; psn down-
# regulated by excess C). The first KA-kernel BGC phase reverse-validated with a MANUAL
# backend launch + struct-grouped (_CnAllocOut/_CnAllocIn) kernel args — reverses at machine
# precision (scripts/enzyme_bgc_reverse.jl [P8]). All C/N structs + fpg_col + cn_shared_params/
# patch/crop flow LIVE from the bundle's inst (the whole inst is Duplicated, so the freshly-
# built PftConNutrientCompetition is the only Const → no Enzyme aliasing). pftcon assembled
# from the global pftcon exactly as cn_driver's _pftnc.
function alloc_rev_aux(inst, bounds, filt, config)
    pftnc = PftConNutrientCompetition{Float64}(
        woody = Float64.(pftcon.woody), froot_leaf = Float64.(pftcon.froot_leaf),
        croot_stem = Float64.(pftcon.croot_stem), stem_leaf = Float64.(pftcon.stem_leaf),
        flivewd = Float64.(pftcon.flivewd), leafcn = Float64.(pftcon.leafcn),
        frootcn = Float64.(pftcon.frootcn), livewdcn = Float64.(pftcon.livewdcn),
        deadwdcn = Float64.(pftcon.deadwdcn), fcur = Float64.(pftcon.fcur),
        graincn = Float64.(pftcon.graincn), grperc = Float64.(pftcon.grperc),
        grpnow = Float64.(pftcon.grpnow), fleafcn = Float64.(pftcon.fleafcn),
        ffrootcn = Float64.(pftcon.ffrootcn), fstemcn = Float64.(pftcon.fstemcn),
        astemf = Float64.(pftcon.astemf), season_decid = Float64.(pftcon.season_decid),
        stress_decid = Float64.(pftcon.stress_decid))
    return (; pftcon = pftnc, mask = filt.bgc_vegp, bounds = bounds.begp:bounds.endp,
              use_c13 = config.use_c13, use_c14 = config.use_c14,
              npcropmin = config.npcropmin, nrepr = config.nrepr, dt = 1800.0)
end
function alloc_rev_phase!(b, aux)
    i = b.inst; veg = i.bgc_vegetation
    calc_plant_cn_alloc!(aux.mask, aux.bounds, aux.pftcon, i.cn_shared_params, i.patch, i.crop,
                         veg.cnveg_state_inst, veg.cnveg_carbonstate_inst,
                         veg.cnveg_carbonflux_inst, veg.cnveg_nitrogenflux_inst;
                         fpg_col = i.soilbiogeochem_state.fpg_col,
                         cnveg_ns = veg.cnveg_nitrogenstate_inst, dt = aux.dt,
                         use_c13 = aux.use_c13, use_c14 = aux.use_c14,
                         npcropmin = aux.npcropmin, nrepr = aux.nrepr)
    return nothing
end

# soil_biogeochem_decomp! — the ACTUAL decomposition flux (the cascade step that consumes the
# competition-resolved potential losses to move C between decomp pools): decomp_cascade_hr_vr =
# rf_decomp_cascade·p_decomp_cpool_loss (heterotrophic respiration) + the inter-pool C transfer
# + N mineralization. Completes the soil-C decomposition cascade ([P5] rate → [P6] potential →
# [P9] actual). Reverse-validated at machine precision (scripts/enzyme_bgc_reverse.jl [P9]).
# cf/nf (written) + cs/ns/st/cascade_con/params (read) all flow LIVE from the bundle's inst;
# the scratch in/out arrays (cn_decomp_pools / p_decomp_cpool_loss / pmnf_decomp_cascade /
# p_decomp_npool_to_din) + dzsoi_decomp are Const aux, sized/sourced as cn_driver does.
function soildecomp_rev_aux(inst, bounds, filt, config)
    return (; mask = filt.bgc_soilc, bounds = bounds.begc:bounds.endc, nlevdecomp = varpar.nlevdecomp,
              ndecomp_pools = config.ndecomp_pools,
              ndecomp_cascade_transitions = config.ndecomp_cascade_transitions, dzsoi_decomp = dzsoi_decomp[],
              use_nitrif_denitrif = _cfgget(config, :use_nitrif_denitrif, false))
end
function soildecomp_rev_phase!(b, aux)
    i = b.inst; s = b.bgc
    soil_biogeochem_decomp!(i.soilbiogeochem_carbonflux, i.soilbiogeochem_carbonstate,
                            i.soilbiogeochem_nitrogenflux, i.soilbiogeochem_nitrogenstate,
                            i.soilbiogeochem_state, i.decomp_cascade, i.decomp_params;
                            mask_bgc_soilc = aux.mask, bounds = aux.bounds, nlevdecomp = aux.nlevdecomp,
                            ndecomp_pools = aux.ndecomp_pools,
                            ndecomp_cascade_transitions = aux.ndecomp_cascade_transitions,
                            cn_decomp_pools = s.cn_decomp_pools,
                            p_decomp_cpool_loss = s.p_decomp_cpool_loss,
                            pmnf_decomp_cascade = s.pmnf_decomp_cascade,
                            p_decomp_npool_to_din = s.p_decomp_npool_to_din,
                            dzsoi_decomp = aux.dzsoi_decomp,
                            use_nitrif_denitrif = aux.use_nitrif_denitrif)
    return nothing
end

# soil_bgc_competition! — the plant-vs-microbe mineral-N COMPETITION split: resolves the
# fraction of immobilization demand met (fpi_vr) and the fraction of plant N demand met
# (fpg_col, consumed by allocation) from the available mineral N vs the total immobilization
# + plant + nitrif/denit demand, then sets the actual_immob / smin_*_to_plant offered fluxes.
# In the N-limited regime fpi = avail_N/total_demand (smooth ratio); reverse-validated at
# machine precision (scripts/enzyme_bgc_reverse.jl [P10]). st/nf/cf (written) + ns (read) flow
# LIVE from the bundle's inst; competition_state/params + cascade_receiver_pool read live too;
# the scratch arrays (pmnf_decomp_cascade / p_decomp_cn_gain) + dzsoi_decomp are Const aux.
# use_fun=false here → fun_hook=nothing (the FUN cost-based-uptake closure is a separate path).
function competition_rev_aux(inst, bounds, filt, config)
    return (; mask = filt.bgc_soilc, bounds = bounds.begc:bounds.endc, nlevdecomp = varpar.nlevdecomp,
              ndecomp_cascade_transitions = config.ndecomp_cascade_transitions, dzsoi_decomp = dzsoi_decomp[],
              use_nitrif_denitrif = _cfgget(config, :use_nitrif_denitrif, false))
end
function competition_rev_phase!(b, aux)
    i = b.inst; s = b.bgc
    soil_bgc_competition!(i.soilbiogeochem_state, i.soilbiogeochem_nitrogenflux,
                          i.soilbiogeochem_carbonflux, i.soilbiogeochem_nitrogenstate,
                          i.competition_state, i.competition_params;
                          mask_bgc_soilc = aux.mask, bounds = aux.bounds, nlevdecomp = aux.nlevdecomp,
                          ndecomp_cascade_transitions = aux.ndecomp_cascade_transitions,
                          dzsoi_decomp = aux.dzsoi_decomp, pmnf_decomp_cascade = s.pmnf_decomp_cascade,
                          p_decomp_cn_gain = s.p_decomp_cn_gain,
                          cascade_receiver_pool = i.decomp_cascade.cascade_receiver_pool,
                          use_nitrif_denitrif = aux.use_nitrif_denitrif, use_fun = false)
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

# --------------------------------------------------------------------------
# clm_drv_reverse! — the TOP-LEVEL clm_drv! reverse entry point. Given a warmed-up
# inst (downscaled forcing set) + a seed, it reverse-differentiates one timestep's
# compositional chain end-to-end and returns the INPUT-gradient bundle `db`
# (db.inst.<field> = dL/d(that field at step entry)). MULTI-PATCH and CONVERGENCE-
# AWARE: the canopy block's Newton iteration count is auto-detected per step from the
# decomposed forward's per-patch convergence (driver_canopy_converged_n) rather than
# hard-coded — so it adapts when the solve converges early or late. The canopy block is
# OPTIONAL (include_canopy): with it the chain is canopy → soil_temp → surface-hydrology
# → soil_water → water_table → hydrology_no_drainage; without it the hydrology-only chain.
# Set use_psn to differentiate through photosynthesis (the stomatal feedback).
#
#     db = clm_drv_reverse!(inst, bounds, filt, config;
#           seed_bang! = (db, b) -> (db.inst.temperature.t_veg_patch .= 2 .* b.inst.temperature.t_veg_patch),
#           include_canopy = true, use_psn = false)
#     # db.inst.temperature.t_grnd_col[c] == dL/d(t_grnd at step entry), L = sum(abs2, t_veg)
#
# `seed_bang!(db, b)` seeds the adjoint from the final forward state (default:
# L = sum(abs2, t_soisno_col[:, nlevsno+1:end]) — the soil thermal profile). The inst
# is left holding the step's final forward state. Validated on Julia 1.10 + Enzyme by
# scripts/enzyme_clm_drv_reverse.jl (FD-checked dL/d(t_grnd) through the whole step,
# auto-detected N); the forward orchestration is guarded (no Enzyme) in
# test/test_driver_reverse.jl.
function clm_drv_reverse!(inst, bounds, filt, config;
        seed_bang! = _default_drv_seed!, include_canopy::Bool = true, use_psn::Bool = false,
        n_canopy::Union{Int,Nothing} = nothing, dtime = 1800.0, include_surface::Bool = true)
    b = driver_rev_bundle(inst)
    canopy_aux = include_canopy ? canopy_rev_aux(inst, bounds, filt; use_psn, dtime) : nothing
    # Convergence-aware N: probe the canopy forward on a throwaway view (no inst mutation).
    Nc = if !include_canopy
        0
    elseif n_canopy === nothing
        driver_canopy_converged_n(b, canopy_aux)
    else
        n_canopy
    end
    phases = driver_rev_phases(bounds, filt, config; canopy_aux = canopy_aux,
                               n_canopy = max(Nc, 1), dtime = dtime, include_surface = include_surface)
    return compositional_reverse!(phases, b, seed_bang!)
end

# Default driver-reverse seed: L = sum(abs2, soil-layer t_soisno) → dL = 2·t_soisno on
# the active soil layers (the soil thermal profile, the natural whole-step energy output).
function _default_drv_seed!(db, b)
    j0 = varpar.nlevsno + 1
    ts  = b.inst.temperature.t_soisno_col
    dts = db.inst.temperature.t_soisno_col
    @views dts[:, j0:end] .= 2 .* ts[:, j0:end]
    return nothing
end

# --------------------------------------------------------------------------
# PER-STEP FORCING injection — for multistep_reverse! over a TRAJECTORY of timesteps with
# TIME-VARYING atmospheric forcing (a real diurnal/seasonal cycle, not the fixed-forcing
# fill(phases,N) demo). `forcingset_rev_phase!` overwrites the live DOWNSCALED forcing arrays
# in b.inst.atm2lnd with this step's snapshot (Const aux). Because the downstream phases read
# forcing LIVE — soiltemp_rev_phase! takes i.atm2lnd.forc_lwrad_downscaled_col directly, and
# canopy_rev_aux's `forc` NamedTuple ALIASES the same a2l arrays by reference (not copies) —
# mutating them in place advances the forcing seen by the whole step WITHOUT rebuilding any aux.
# A pure copy reverses trivially: forcing is exogenous, so its adjoint flows to the Const source
# and is discarded — correctly NOT leaking into the carried state's gradient.
#
# Sets the four radiative/thermodynamic drivers that feed soil_temperature! + the canopy energy
# balance (lwrad, t, th, rho). Keep forc_vp/forc_pbot fixed across the schedule so canopy's
# derived forc_q snapshot stays valid; vary these four for a forcing-driven trajectory.
function forcingset_rev_phase!(b, aux)
    a2l = b.inst.atm2lnd
    copyto!(a2l.forc_lwrad_downscaled_col, aux.lwrad)
    copyto!(a2l.forc_t_downscaled_col,     aux.t)
    copyto!(a2l.forc_th_downscaled_col,    aux.th)
    copyto!(a2l.forc_rho_downscaled_col,   aux.rho)
    return nothing
end

# Snapshot the current downscaled forcing as a forcingset aux (Const copies). Use as the base
# for a schedule (scale/offset per step). The arrays are the COLUMN-indexed downscaled fields.
forcingset_aux(inst) = (; lwrad = copy(inst.atm2lnd.forc_lwrad_downscaled_col),
                          t     = copy(inst.atm2lnd.forc_t_downscaled_col),
                          th    = copy(inst.atm2lnd.forc_th_downscaled_col),
                          rho   = copy(inst.atm2lnd.forc_rho_downscaled_col))

# Build the per-timestep `steps` list for multistep_reverse!/multistep_reverse_binomial! from a
# forcing `schedule` (a Vector of forcingset auxes, one per step). Each step = a forcing-set
# phase followed by the shared driver phase list (canopy+hydrology). The base phases are built
# ONCE and reused (their aux don't capture forcing — soiltemp reads it live, canopy_aux aliases
# the live arrays); only the leading forcingset phase differs per step. With canopy_aux passed,
# the canopy energy/photosynthesis block sees each step's forcing too (live alias).
function forced_driver_steps(bounds, filt, config, schedule; canopy_aux = nothing,
                             n_canopy::Int = 12, dtime = 1800.0)
    base = driver_rev_phases(bounds, filt, config; canopy_aux = canopy_aux,
                             n_canopy = n_canopy, dtime = dtime)
    return Any[vcat(Any[(forcingset_rev_phase!, (fk,))], base) for fk in schedule]
end

# --------------------------------------------------------------------------
# BGC (use_cn) whole-step assembler — the ordered (phase_fn, const_args) list for the
# REVERSE-READY subset of cn_driver_no_leaching!, in the SAME forward order the production
# CN driver runs them (cn_driver.jl line numbers in parens):
#   cn_mresp!(225) → decomp_rate_constants_bgc!(325) → soil_bgc_potential!(343) →
#   soil_bgc_competition!(432) → calc_plant_cn_alloc!(453) → soil_biogeochem_decomp!(464) →
#   cn_gresp!(566) → c_state_update1!(585) → n_state_update1!(612) →
#   soilbiogeochem_n_state_update1!(636)
# i.e. respiration → the soil-C/N decomposition+competition+allocation cascade → the C/N pool
# integrations. Mirrors driver_rev_phases (hydrology) but the BGC aux builders source struct
# fields from `inst`, so this takes inst. The phases couple through inst fields (cf.decomp_k →
# nf.potential_immob → st.fpg_col → cnveg pools → smin_nh4), so a single compositional_reverse!
# over this list flows the gradient across the whole BGC step. Each phase is individually
# FD-validated at machine precision (scripts/enzyme_bgc_reverse.jl [P1..P10]); the coupled
# chain is validated by scripts/enzyme_bgc_wholestep.jl. NOT included: the genuinely discrete
# cn_driver calls (cn_phenology! onset/offset FSMs, cn_gap_mortality!, nitrif_denitrif!
# thresholds) + the demand pre-pass (calc_plant_nutrient_demand!), per the "discontinuities
# deferred to Phase 3" policy.
function bgc_rev_phases(inst, bounds, filt, config)
    return Any[
        (cnmresp_rev_phase!,     (cnmresp_rev_aux(inst, bounds, filt),)),
        (decomprate_rev_phase!,  (decomprate_rev_aux(inst, bounds, filt),)),
        (decomppot_rev_phase!,   (decomppot_rev_aux(inst, bounds, filt, config),)),
        (competition_rev_phase!, (competition_rev_aux(inst, bounds, filt, config),)),
        (alloc_rev_phase!,       (alloc_rev_aux(inst, bounds, filt, config),)),
        (soildecomp_rev_phase!,  (soildecomp_rev_aux(inst, bounds, filt, config),)),
        (cngresp_rev_phase!,     (cngresp_rev_aux(inst, bounds, filt),)),
        (cstate1_rev_phase!,     (cstate1_rev_aux(inst, bounds, filt, config),)),
        (nstate1_rev_phase!,     (nstate1_rev_aux(inst, bounds, filt, config),)),
        (sminnupdate_rev_phase!, (sminnupdate_rev_aux(inst, bounds, filt, config),)),
    ]
end

# --------------------------------------------------------------------------
# FULL clm_drv! whole-step assembler — the COMPLETE default timestep as ONE ordered phase
# list: the canopy+hydrology block (driver_rev_phases) followed by the use_cn BGC block
# (bgc_rev_phases), in the exact clm_drv! order. The production driver runs cn_driver_no_
# leaching! (the BGC step) at clm_driver.jl:1402 — right AFTER the HydrologyNoDrainage block
# (canopy → soil_temp → surface-hydrology → soil_water → water_table → hydrology_no_drainage!)
# and before drainage — so concatenation is faithful. ALL phases share ONE unified bundle
# (bgc_rev_bundle: inst + canopy `scratch` + decomp `bgc` scratch); the canopy/hydrology
# phases ignore b.bgc and the BGC phases ignore b.scratch. A single compositional_reverse!
# over this list flows the gradient across all three domains: e.g. perturb t_grnd → canopy →
# soil_temperature (t_soisno) / soil_water (soilpsi) → BGC decomposition → leaf carbon.
# Validated by scripts/enzyme_fullstep_all.jl (Julia 1.10/Enzyme). Needs a use_cn inst (the
# BGC aux builders read inst.decomp_cascade / competition / bgc_vegetation).
function full_rev_phases(inst, bounds, filt, config; canopy_aux = nothing,
                         n_canopy::Int = 12, dtime = 1800.0)
    return vcat(driver_rev_phases(bounds, filt, config; canopy_aux = canopy_aux,
                                  n_canopy = n_canopy, dtime = dtime),
                bgc_rev_phases(inst, bounds, filt, config))
end
