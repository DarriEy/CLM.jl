# ==========================================================================
# clmdrv_cn_fixture.jl — populate a fixture inst with a consistent finite CN
# cold start, so use_cn=true clm_drv! yields FINITE BGC outputs.
#
# make_driver_data_physical() (clmdrv_make_data.jl) builds an SP-mode fixture:
# the CN-veg facade arrays are left UNALLOCATED (size 0) and the soil-BGC pools
# (decomp_cpools_vr_col) are NaN — so the CN half of the composite runs over empty
# patch arrays (no-ops) and the soil pools stay NaN. This helper does what a real
# CLM cold start does: allocate the CN-veg facade, then run every CN + soil-BGC
# init_cold! routine, giving a self-consistent finite initial state.
#
# NOTE (ceiling): the SP fixture leaves the surface-energy-balance / radiation
# fields NaN (albedo<-hydrology<-energy-balance circular dependency, see the note
# at the end of make_driver_data_physical). Photosynthesis therefore feeds NaN
# psn->cpool fluxes into the CN-VEG update, so leafc can stay NaN even with a cold
# CN state. The SOIL-BGC decomposition (decomp_cpools_vr) does not depend on that
# radiation and becomes finite — that is the parity signal this fixture unlocks.
# ==========================================================================

"""
    populate_cn_cold!(inst; nc, np, ng)

Allocate the CN-veg facade of a fixture `inst` and cold-init every CN + soil-BGC
sub-state to a finite, self-consistent starting point. nc/np/ng are the fixture's
column/patch/gridcell counts. Returns `inst`.
"""
function populate_cn_cold!(inst; nc::Int, np::Int, ng::Int)
    scs = inst.soilbiogeochem_carbonstate
    # decomp geometry is fixed by the already-allocated soil-C pool array
    nlevdecomp   = size(scs.decomp_cpools_vr_col, 2)
    ndecomp_pools = size(scs.decomp_cpools_vr_col, 3)
    # CENTURY cascade transition count (matches the driver's CENTURY decomp config)
    ndecomp_cascade_transitions = 8

    # --- allocate the CN-veg facade (was size 0 in the SP fixture) ---
    # The SP fixture built the facade with config.use_cn=false, so cn_vegetation_init!
    # would skip the carbon/nitrogen allocation branch. Flip it on first.
    veg = inst.bgc_vegetation
    veg.config.use_cn = true
    CLM.cn_vegetation_init!(veg, np, nc, ng;
        nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools,
        ndecomp_cascade_transitions=ndecomp_cascade_transitions)

    bp = 1:np; bc = 1:nc

    # --- CN-veg cold inits (carbon / nitrogen / structure / fluxes) ---
    CLM.cnveg_carbon_state_init_cold!(veg.cnveg_carbonstate_inst, bp)
    CLM.cnveg_nitrogen_state_init_cold!(veg.cnveg_nitrogenstate_inst, bp)
    CLM.cnveg_state_init_cold!(veg.cnveg_state_inst, bc, bp)
    CLM.cnveg_carbon_flux_init_cold!(veg.cnveg_carbonflux_inst, bp, bc;
        nlevdecomp_full=nlevdecomp, ndecomp_pools=ndecomp_pools)
    CLM.cnveg_nitrogen_flux_init_cold!(veg.cnveg_nitrogenflux_inst, bp, bc;
        nlevdecomp_full=nlevdecomp, ndecomp_pools=ndecomp_pools)

    # --- soil-BGC cold inits (carbon / nitrogen state; state; fluxes) ---
    CLM.soil_bgc_carbon_state_init_cold!(scs, bc, 1.0;
        nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools)
    CLM.soil_bgc_nitrogen_state_init_cold!(inst.soilbiogeochem_nitrogenstate, bc;
        nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools,
        decomp_cpools_vr_col=Array(scs.decomp_cpools_vr_col),
        decomp_cpools_col=Array(scs.decomp_cpools_col),
        decomp_cpools_1m_col=Array(scs.decomp_cpools_1m_col))
    CLM.soil_bgc_state_init_cold!(inst.soilbiogeochem_state, bc, nlevdecomp)

    return inst
end
