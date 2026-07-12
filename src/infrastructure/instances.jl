# ============================================================================
# Instance factory — creates and initializes all CLM data instances
#
# Ported from: src/main/clm_instMod.F90
#
# The Fortran module declares global instances and calls Init() on each.
# Here we bundle them into a single CLMInstances container and provide
# clm_instInit!() to allocate/initialise everything in one call.
#
# Skipped from Fortran: file I/O, restart, history, namelist reading,
# FATES, ozone factory, dust factory, snow-cover-fraction factory,
# excess-ice streams, soil-water-retention-curve factory.
# ============================================================================

# Forward-declared supertype for the FATES interface container. The concrete
# `fates_interface_type` (FatesInterfaceMod.jl) is `<: AbstractFatesInterface` and
# is `include`d AFTER this file, so the CLMInstances `fates` field is annotated on
# the abstract supertype to satisfy parse-time ordering. The attached value is
# always a concrete `fates_interface_type`; it is excluded from the AD dual-copy and
# GPU adapt exactly like `surfdata`.
abstract type AbstractFatesInterface end

"""
    CLMInstances

Container holding all CLM data type instances.  Corresponds to the module-level
variables declared in `clm_instMod.F90`.
"""
Base.@kwdef mutable struct CLMInstances
    # --- Hierarchy / grid ---
    gridcell::GridcellData           = GridcellData()
    landunit::LandunitData           = LandunitData()
    column::ColumnData               = ColumnData()
    patch::PatchData                 = PatchData()

    # --- Topography ---
    topo::TopoData                   = TopoData()

    # --- Urban ---
    urbanparams::UrbanParamsData     = UrbanParamsData()

    # --- Physics ---
    temperature::TemperatureData     = TemperatureData()
    energyflux::EnergyFluxData       = EnergyFluxData()
    canopystate::CanopyStateData     = CanopyStateData()
    soilstate::SoilStateData         = SoilStateData()
    soilhydrology::SoilHydrologyData = SoilHydrologyData()
    water::WaterData                 = WaterData()
    frictionvel::FrictionVelocityData = FrictionVelocityData()
    lakestate::LakeStateData         = LakeStateData()
    solarabs::SolarAbsorbedData      = SolarAbsorbedData()
    surfalb::SurfaceAlbedoData       = SurfaceAlbedoData()
    surfalb_con::SurfaceAlbedoConstants = SurfaceAlbedoConstants()

    # --- Radiation diagnostics ---
    surfrad::SurfaceRadiationData    = SurfaceRadiationData()

    # --- Photosynthesis ---
    photosyns::PhotosynthesisData    = PhotosynthesisData()

    # --- Ozone ---
    ozone::OzoneData                 = OzoneData()

    # --- Aerosol deposition ---
    aerosol::AerosolData             = AerosolData()

    # --- Irrigation ---
    irrigation::IrrigationData       = IrrigationData()

    # --- Atmosphere / land coupling ---
    atm2lnd::Atm2LndData             = Atm2LndData()
    lnd2atm::Lnd2AtmData             = Lnd2AtmData()

    # --- Crop ---
    crop::CropData                   = CropData()

    # --- Soil biogeochemistry ---
    soilbiogeochem_state::SoilBiogeochemStateData           = SoilBiogeochemStateData()
    soilbiogeochem_carbonstate::SoilBiogeochemCarbonStateData = SoilBiogeochemCarbonStateData()
    soilbiogeochem_carbonflux::SoilBiogeochemCarbonFluxData   = SoilBiogeochemCarbonFluxData()
    soilbiogeochem_nitrogenstate::SoilBiogeochemNitrogenStateData = SoilBiogeochemNitrogenStateData()
    soilbiogeochem_nitrogenflux::SoilBiogeochemNitrogenFluxData   = SoilBiogeochemNitrogenFluxData()

    # C13/C14 soil-BGC carbon state/flux — only sized when use_c13/use_c14 (empty
    # otherwise); the parallel soil pools the CIsoFlux* cascade routes isotope C
    # into (and that c14_decay! radioactively decays), alongside the C13/C14
    # vegetation state on the CN vegetation facade.
    c13_soilbiogeochem_carbonstate::SoilBiogeochemCarbonStateData = SoilBiogeochemCarbonStateData()
    c13_soilbiogeochem_carbonflux::SoilBiogeochemCarbonFluxData   = SoilBiogeochemCarbonFluxData()
    c14_soilbiogeochem_carbonstate::SoilBiogeochemCarbonStateData = SoilBiogeochemCarbonStateData()
    c14_soilbiogeochem_carbonflux::SoilBiogeochemCarbonFluxData   = SoilBiogeochemCarbonFluxData()

    # --- Methane (CH4) — prognostic state + parameters; only used when use_lch4 ---
    ch4::CH4Data                     = CH4Data{Float64}()
    ch4_params::CH4Params            = CH4Params()
    ch4_varcon::CH4VarCon            = CH4VarCon()

    # --- Emissions ---
    dust_emis::DustEmisBaseData      = DustEmisBaseData()
    vocemis::VOCEmisData             = VOCEmisData()
    # MEGAN static descriptors (factors table + compound/mechanism mappings) live
    # on CLMDriverConfig.megan, NOT here — CLMInstances is dual-copied for AD and
    # cannot traverse those non-numeric struct vectors.

    # --- Dry deposition velocity ---
    drydep::DryDepVelocityData = DryDepVelocityData()

    # --- Satellite phenology ---
    satellite_phenology::SatellitePhenologyData = SatellitePhenologyData()

    # --- Active layer ---
    active_layer::ActiveLayerData = ActiveLayerData()

    # --- Balance check ---
    balcheck::BalanceCheckData = BalanceCheckData()

    # --- CN vegetation ---
    bgc_vegetation::CNVegetationData = CNVegetationData()
    cn_products::CNProductsData      = CNProductsData()

    # --- Dynamic global vegetation (CNDV / DGVM) — only used when use_cndv ---
    # dgvs is prognostic patch-level state (parametric {FT}, device-movable): it
    # rides the normal AD dual-copy. dgv_ecophyscon holds read-only ecophysiological
    # constants derived from pftcon; it is parametric on the *vector* type (not FT)
    # so it is left untouched by make_instances and is excluded from the AD dual-copy
    # (like surfdata/fates) in calibration._make_dual_instances.
    dgvs::DGVSData                   = DGVSData()
    dgv_ecophyscon::DGVEcophysCon    = DGVEcophysCon()

    # --- Decomposition cascade ---
    decomp_cascade::DecompCascadeConData = DecompCascadeConData()

    # --- BGC decomposition state and params ---
    decomp_bgc_state::DecompBGCState = DecompBGCState()
    decomp_bgc_params::DecompBGCParams = DecompBGCParams()
    cn_shared_params::CNSharedParamsData = CNSharedParamsData()
    decomp_params::DecompParams = DecompParams()
    competition_state::SoilBGCCompetitionState = SoilBGCCompetitionState()
    competition_params::SoilBGCCompetitionParams = SoilBGCCompetitionParams()
    litter_params::LitterVertTranspParams = LitterVertTranspParams()

    # --- Saturated / infiltration excess runoff ---
    sat_excess_runoff::SaturatedExcessRunoffData = SaturatedExcessRunoffData()
    infilt_excess_runoff::InfiltrationExcessRunoffData = InfiltrationExcessRunoffData()

    # --- Snow cover fraction ---
    scf_method::SnowCoverFractionBase = SnowCoverFractionSwensonLawrence2012()

    # --- Calibration overrides (NaN = use default) ---
    overrides::CalibrationOverrides = CalibrationOverrides()

    # --- Surface input data (for monthly phenology re-reads) ---
    surfdata::Union{SurfaceInputData, Nothing} = nothing

    # FATES interface state (pointer-linked cohort/patch/site hierarchy) — attached
    # like surfdata: excluded from the AD dual-copy and GPU adapt. FATES columns are
    # CPU-only / non-differentiable by design. nothing unless use_fates.
    fates::Union{AbstractFatesInterface, Nothing} = nothing
end

# Make the whole instance tree device-movable: adapt(backend, inst) recursively
# moves the device-movable sub-structs (state structs registered with
# Adapt.@adapt_structure) to the backend's array type and passes the rest through.
Adapt.@adapt_structure CLMInstances

# --- working-precision factory ----------------------------------------------
# Rebuild a default-constructed field at working precision FT. Parametric state
# types expose a `T{FT}()` convenience constructor (T parametric on `FT<:Real`);
# types without one are returned at their default. The `hasmethod` probe is what
# makes this self-maintaining — newly added parametric state fields are picked up
# with no change here. Guarded so a parametric type whose first type parameter is
# not `<:Real` is simply left alone rather than throwing.
function _instance_at(::Type{FT}, x) where {FT}
    W = Base.typename(typeof(x)).wrapper
    if W isa UnionAll
        try
            hasmethod(W{FT}, Tuple{}) && return (W{FT})()
        catch
        end
    end
    return x
end

"""
    make_instances(::Type{FT} = Float64) -> CLMInstances

Construct a `CLMInstances` whose state types use working precision `FT`. With the
default `Float64` this is equivalent to `CLMInstances()`. With `Float32`, every
parametric state type (temperature, energy flux, soil state/hydrology, water
state, photosynthesis, the biogeochem states/fluxes, …) is allocated in Float32 —
the precision required by the Apple Metal GPU backend, which has no Float64.

Not yet reparametrised (kept Float64): the non-parametric *container* structs
`water` and `bgc_vegetation`, and the decomposition/competition *parameter* sets.
These are physical-parameter or aggregate types, not on the current GPU kernel
path; threading `FT` through them (and a fully-Float32 model run) additionally
requires finishing driver kernelisation — see the PRD.

Pair with [`clm_float_type`](@ref) to choose `FT` from the target backend, e.g.
`make_instances(clm_float_type(backend = :metal))`.
"""
function make_instances(::Type{FT} = Float64) where {FT}
    base = CLMInstances()
    return CLMInstances((_instance_at(FT, getfield(base, f)) for f in fieldnames(CLMInstances))...)
end

"""
    clm_instInit!(inst::CLMInstances;
                  ng, nl, nc, np,
                  nlevdecomp_full, ndecomp_pools, ndecomp_cascade_transitions)

Allocate and cold-initialise every data instance inside `inst`.

Dimension arguments
- `ng` — number of gridcells
- `nl` — number of landunits
- `nc` — number of columns
- `np` — number of patches

Biogeochemistry dimensions (keyword, with defaults)
- `nlevdecomp_full`            — full number of decomposition levels (default 10)
- `ndecomp_pools`              — number of decomposition pools      (default 7)
- `ndecomp_cascade_transitions`— number of cascade transitions      (default 5)
"""
function clm_instInit!(inst::CLMInstances;
                       ng::Int,
                       nl::Int,
                       nc::Int,
                       np::Int,
                       nlevdecomp_full::Int = 10,
                       ndecomp_pools::Int = 7,
                       ndecomp_cascade_transitions::Int = 5,
                       nlevurb::Int = 0,
                       use_luna::Bool = false,
                       use_lch4::Bool = false,
                       use_cndv::Bool = false,
                       use_c13::Bool = false,
                       use_c14::Bool = false)

    # --- Grid hierarchy ---
    gridcell_init!(inst.gridcell, ng)
    landunit_init!(inst.landunit, nl)
    column_init!(inst.column, nc)
    patch_init!(inst.patch, np)

    # --- Urban ---
    urbanparams_init!(inst.urbanparams, nl; nlevurb=nlevurb)

    # --- Topography (allocate only; cold-start requires valid col/lun linkage) ---
    topo_init_allocate!(inst.topo, nc)

    # --- Core physics ---
    temperature_init!(inst.temperature, np, nc, nl, ng)
    energyflux_init!(inst.energyflux, np, nc, nl, ng)
    canopystate_init!(inst.canopystate, np)
    soilstate_init!(inst.soilstate, np, nc)
    soilhydrology_init!(inst.soilhydrology, nc)
    water_init!(inst.water, nc, np, nl, ng)
    frictionvel_init!(inst.frictionvel, np, nc)
    frictionvel_read_nml!(inst.frictionvel)
    frictionvel_read_params!(inst.frictionvel)
    lakestate_init!(inst.lakestate, nc, np)
    solarabs_init!(inst.solarabs, np, nl; use_luna=use_luna)
    surfalb_init!(inst.surfalb, np, nc, ng)

    # --- Radiation diagnostics ---
    surfrad_init!(inst.surfrad, np)

    # --- Photosynthesis ---
    photosynthesis_data_init!(inst.photosyns, np; use_luna=use_luna)

    # --- Ozone ---
    ozone_init!(inst.ozone, np)

    # --- Aerosol deposition ---
    # Fortran `aerosol_inst%Init` = InitAllocate + InitHistory + InitCold, and InitCold
    # zeros every mss_*/mss_cnc_* snow-layer array. aerosol_init! is the InitAllocate
    # analogue: it NaN-fills (matching Fortran's allocate-to-nan), so the InitCold zeroing
    # is what makes those arrays usable. Without it the mass arrays only ever get zeroed
    # inside aerosol_masses!'s "layer j is above snl" branch — which runs at the END of a
    # step. A column that forms its first snow layer during step 1 (e.g. a cold-start
    # glacier under snowfall) therefore reaches aerosol_fluxes! / SNICAR with allocation-
    # time NaN in mss_bcphi[c, top] → NaN mss_cnc → NaN SNICAR albedo + flx_abs → NaN
    # sabg_lyr → NaN soil-temperature solve on the NEXT step. Run it here (not only in
    # cold_start_initialize!) so the restart path is also zeroed before the dump overwrites,
    # exactly as Fortran's Init does.
    aerosol_init!(inst.aerosol, nc)
    aerosol_init_cold!(inst.aerosol, 1:nc)

    # --- Irrigation ---
    irrigation_init_allocate!(inst.irrigation, np, nc, varpar.nlevsoi)

    # --- Atmosphere / land coupling ---
    atm2lnd_init!(inst.atm2lnd, ng, nc, np)
    lnd2atm_init!(inst.lnd2atm, ng, nc)

    # --- Crop ---
    crop_init!(inst.crop, np)

    # --- Soil biogeochemistry ---
    soil_bgc_state_init!(inst.soilbiogeochem_state, nc, np,
                         nlevdecomp_full, ndecomp_cascade_transitions)
    soil_bgc_carbon_state_init!(inst.soilbiogeochem_carbonstate, nc, ng,
                                nlevdecomp_full, ndecomp_pools)
    soil_bgc_carbon_flux_init!(inst.soilbiogeochem_carbonflux, nc,
                               nlevdecomp_full, ndecomp_pools,
                               ndecomp_cascade_transitions)
    soil_bgc_nitrogen_state_init!(inst.soilbiogeochem_nitrogenstate, nc, ng,
                                  nlevdecomp_full, ndecomp_pools)
    soil_bgc_nitrogen_flux_init!(inst.soilbiogeochem_nitrogenflux, nc,
                                 nlevdecomp_full, ndecomp_pools,
                                 ndecomp_cascade_transitions)

    # C13/C14 soil-BGC carbon state/flux — size only when the tracer is active.
    if use_c13
        soil_bgc_carbon_state_init!(inst.c13_soilbiogeochem_carbonstate, nc, ng,
                                    nlevdecomp_full, ndecomp_pools)
        soil_bgc_carbon_flux_init!(inst.c13_soilbiogeochem_carbonflux, nc,
                                   nlevdecomp_full, ndecomp_pools,
                                   ndecomp_cascade_transitions)
    end
    if use_c14
        soil_bgc_carbon_state_init!(inst.c14_soilbiogeochem_carbonstate, nc, ng,
                                    nlevdecomp_full, ndecomp_pools)
        soil_bgc_carbon_flux_init!(inst.c14_soilbiogeochem_carbonflux, nc,
                                   nlevdecomp_full, ndecomp_pools,
                                   ndecomp_cascade_transitions)
    end

    # --- Methane (CH4) — size the prognostic state when active ---
    if use_lch4
        ch4_init_allocate!(inst.ch4, nc, np, ng, varpar.nlevsoi)
    end

    # --- Emissions ---
    dust_emis_init!(inst.dust_emis, np; nc=nc)
    vocemis_init!(inst.vocemis, np, ng, 20, 20)
    # MEGAN factors live on CLMDriverConfig.megan and are initialized there.

    # --- Dry deposition velocity ---
    drydep_init!(inst.drydep, np, 0)  # n_drydep=0 by default; set >0 to activate

    # --- Satellite phenology ---
    satellite_phenology_init!(inst.satellite_phenology, np)

    # --- Active layer ---
    active_layer_init!(inst.active_layer, nc)

    # --- Balance check ---
    balance_check_init!(inst.balcheck, 1800.0)

    # --- CN vegetation ---
    cn_vegetation_init!(inst.bgc_vegetation, np, nc, ng;
                        nlevdecomp=nlevdecomp_full,
                        ndecomp_pools=ndecomp_pools,
                        ndecomp_cascade_transitions=ndecomp_cascade_transitions)
    cn_products_init!(inst.cn_products, ng)

    # --- Dynamic global vegetation (CNDV) — only sized/initialized when active ---
    # Cold-start the DGVS patch state and derive the ecophysiological constants from
    # the global pftcon (must already be allocated). Left at size 0 when use_cndv is
    # false so the default path is byte-identical (the gated driver calls no-op on
    # an empty dgvs via their isempty(bounds) guards).
    if use_cndv
        dgvs_init_cold!(inst.dgvs, np)
        if !isempty(pftcon.pftpar20)
            dgv_ecophyscon_init!(inst.dgv_ecophyscon, pftcon)
        end
    end

    # --- Saturated / infiltration excess runoff ---
    inst.sat_excess_runoff.fsat_col = fill(0.0, nc)
    inst.sat_excess_runoff.fcov_col = fill(0.0, nc)
    inst.infilt_excess_runoff.qinmax_col = fill(0.0, nc)

    # --- Decomposition cascade (allocate with defaults) ---
    decomp_cascade_con_init!(inst.decomp_cascade,
                              ndecomp_pools, ndecomp_cascade_transitions)

    return nothing
end

"""
    decomp_cascade_con_init!(cascade, ndecomp_pools, ndecomp_cascade_transitions)

Allocate decomposition cascade configuration arrays with default sizes.
The actual cascade parameters are set by `init_decomp_cascade_bgc!` or
`init_decomp_cascade_mimics!` during model initialization.
"""
function decomp_cascade_con_init!(cascade::DecompCascadeConData,
                                   ndecomp_pools::Int,
                                   ndecomp_cascade_transitions::Int)
    cascade.cascade_donor_pool             = zeros(Int, ndecomp_cascade_transitions)
    cascade.cascade_receiver_pool          = zeros(Int, ndecomp_cascade_transitions)
    cascade.floating_cn_ratio_decomp_pools = falses(ndecomp_pools)
    cascade.is_litter                      = falses(ndecomp_pools)
    cascade.is_soil                        = falses(ndecomp_pools)
    cascade.is_cwd                         = falses(ndecomp_pools)
    cascade.initial_cn_ratio               = zeros(ndecomp_pools)
    cascade.initial_stock                  = zeros(ndecomp_pools)
    cascade.is_metabolic                   = falses(ndecomp_pools)
    cascade.is_cellulose                   = falses(ndecomp_pools)
    cascade.is_lignin                      = falses(ndecomp_pools)
    cascade.spinup_factor                  = ones(ndecomp_pools)
    cascade.decomp_pool_name_restart       = fill("", ndecomp_pools)
    cascade.decomp_pool_name_history       = fill("", ndecomp_pools)
    cascade.decomp_pool_name_long          = fill("", ndecomp_pools)
    cascade.decomp_pool_name_short         = fill("", ndecomp_pools)
    cascade.cascade_step_name              = fill("", ndecomp_cascade_transitions)
    return nothing
end
