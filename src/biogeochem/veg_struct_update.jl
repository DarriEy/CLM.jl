# ==========================================================================
# Ported from: src/biogeochem/CNVegStructUpdateMod.F90
# Vegetation structure updates (LAI, SAI, htop, hbot) from C state variables.
#
# On the radiation time step, uses C state variables and ecological process
# constants (epc) to diagnose vegetation structure (LAI, SAI, height).
#
# Public functions:
#   cn_veg_struct_update! -- Update vegetation structure from live C pools
#
# GPU kernelization: the per-patch diagnosis is one KernelAbstractions kernel
# (one thread per patch). Every write is to the patch's own index (no patch->
# column scatter), so it is race-free and byte-identical to the sequential host
# loop on the KA CPU backend. The ~30 patch/column arrays + 11 PFT-parameter
# vectors the body touches are grouped into `Adapt.@adapt_structure` device-view
# bundles (mirroring the original variable names) to stay under Metal's ~31-arg
# launch limit; the PFT params (a host-global pftcon) are moved onto the state
# backend at launch via `_to_backend_like`. All literals are carried at the
# working element type `T` so the kernel compiles on Metal (no Float64) while
# the Float64 host path is unchanged.
# ==========================================================================

# ---------------------------------------------------------------------------
# Module-level constants
# ---------------------------------------------------------------------------

const DTSMONTH = 2592000.0            # seconds in a 30-day month (60*60*24*30)
const FRAC_SNO_THRESHOLD = 0.999      # frac_sno values > this treated as 1

# ---------------------------------------------------------------------------
# Device-view bundles (group the many patch/column arrays into few kernel args)
# ---------------------------------------------------------------------------

# Written outputs: canopy structure (V=float patch vec) + the two Int flags (VI).
Base.@kwdef struct _VSUOut{V,VI}
    tlai::V; tsai::V; htop::V; hbot::V; elai::V; esai::V
    stem_biomass::V; leaf_biomass::V; htmx::V
    frac_veg_nosno_alb::VI; peaklai::VI
end
Adapt.@adapt_structure _VSUOut

# Read-only inputs (patch + column vectors).
Base.@kwdef struct _VSUIn{V}
    leafc::V; deadstemc::V; livestemc::V
    forc_hgt_u::V; frac_sno::V; snow_depth::V; farea_burned::V
end
Adapt.@adapt_structure _VSUIn

# PFT-parameter vectors (indexed by ivt+1); all Float in pftcon.
Base.@kwdef struct _VSUPar{V}
    woody::V; slatop::V; dsladlai::V; z0mr::V; displar::V; dwood::V
    ztopmx::V; laimx::V; nstem::V; taper::V; fbw::V
end
Adapt.@adapt_structure _VSUPar

# Integer index vectors.
Base.@kwdef struct _VSUIdx{VI}
    ivt::VI; col::VI; harvdate::VI
end
Adapt.@adapt_structure _VSUIdx

# Scalars (isbits — passes to the kernel directly, no adapt). Floats carried at
# the working precision T; PFT-type indices and mode flags as Int/Bool.
struct _VSUScalars{T}
    dt::T; spinup_factor_deadwood::T; c_to_b::T
    use_cndv::Bool; use_biomass_heat_storage::Bool
    noveg::Int; nc3crop::Int; nc3irrig::Int; nbrdlf_dcd_brl_shrub::Int; npcropmin::Int
    ntmp_corn::Int; nirrig_tmp_corn::Int; ntrp_corn::Int; nirrig_trp_corn::Int
    nsugarcane::Int; nirrig_sugarcane::Int; nmiscanthus::Int; nirrig_miscanthus::Int
    nswitchgrass::Int; nirrig_switchgrass::Int
end

# ---------------------------------------------------------------------------
# Per-patch kernel — one thread per patch (byte-identical to the host loop on
# the KA CPU backend; every write is to the thread's own patch index).
# ---------------------------------------------------------------------------
@kernel function _vsu_kernel!(@Const(mask), idx, inb, out, par, sc)
    p = @index(Global)
    @inbounds if mask[p] && isfinite(inb.leafc[p])
        T = eltype(out.tlai)

        # aliases (match the Fortran/associate names) — compile-time, zero cost
        ivt = idx.ivt; col = idx.col; harvdate = idx.harvdate
        woody = par.woody; slatop = par.slatop; dsladlai = par.dsladlai
        z0mr = par.z0mr; displar = par.displar; dwood = par.dwood
        ztopmx = par.ztopmx; laimx = par.laimx; nstem = par.nstem
        taper = par.taper; fbw = par.fbw
        leafc = inb.leafc; deadstemc = inb.deadstemc; livestemc = inb.livestemc
        forc_hgt_u_patch = inb.forc_hgt_u; frac_sno = inb.frac_sno
        snow_depth = inb.snow_depth; farea_burned = inb.farea_burned
        tlai = out.tlai; tsai = out.tsai; htop = out.htop; hbot = out.hbot
        elai = out.elai; esai = out.esai; stem_biomass = out.stem_biomass
        leaf_biomass = out.leaf_biomass; htmx = out.htmx
        frac_veg_nosno_alb = out.frac_veg_nosno_alb; peaklai = out.peaklai
        dt = sc.dt; spinup_factor_deadwood = sc.spinup_factor_deadwood; c_to_b = sc.c_to_b
        use_cndv = sc.use_cndv; use_biomass_heat_storage = sc.use_biomass_heat_storage
        noveg = sc.noveg; nc3crop = sc.nc3crop; nc3irrig = sc.nc3irrig
        nbrdlf_dcd_brl_shrub = sc.nbrdlf_dcd_brl_shrub; npcropmin = sc.npcropmin
        ntmp_corn = sc.ntmp_corn; nirrig_tmp_corn = sc.nirrig_tmp_corn
        ntrp_corn = sc.ntrp_corn; nirrig_trp_corn = sc.nirrig_trp_corn
        nsugarcane = sc.nsugarcane; nirrig_sugarcane = sc.nirrig_sugarcane
        nmiscanthus = sc.nmiscanthus; nirrig_miscanthus = sc.nirrig_miscanthus
        nswitchgrass = sc.nswitchgrass; nirrig_switchgrass = sc.nirrig_switchgrass

        c = col[p]

        if ivt[p] != noveg

            tlai_old = tlai[p]  # n-1 value
            tsai_old = tsai[p]  # n-1 value

            # Update the leaf area index based on leafC and SLA
            # Eq 3 from Thornton and Zimmerman, 2007, J Clim, 20, 3902-3923.
            if dsladlai[ivt[p] + 1] > zero(T)
                tlai[p] = (slatop[ivt[p] + 1] * (exp(leafc[p] * dsladlai[ivt[p] + 1]) - one(T))) / dsladlai[ivt[p] + 1]
            else
                tlai[p] = slatop[ivt[p] + 1] * leafc[p]
            end
            tlai[p] = max(zero(T), tlai[p])   # HARD: constant-0 branch on an LAI [m2/m2] axis. A leafless/dormant patch is exactly 0, so this ReLU sat on its kink and gave it +0.0139 of permanent phantom LAI.

            # Update the stem area index and height based on LAI, stem mass, and veg type.
            # tsai formula from Zeng et al. 2002, Journal of Climate, p1835
            if ivt[p] == nc3crop || ivt[p] == nc3irrig  # generic crops
                tsai_alpha = one(T) - one(T) * dt / T(DTSMONTH)
                tsai_min = T(0.1)
            else
                tsai_alpha = one(T) - T(0.5) * dt / T(DTSMONTH)
                tsai_min = one(T)
            end
            tsai_min = tsai_min * T(0.5)
            tsai[p] = smooth_max(tsai_alpha * tsai_old + smooth_max(tlai_old - tlai[p], zero(T)), tsai_min)

            # Vegetation physiological parameters used in biomass heat storage
            if use_biomass_heat_storage
                leaf_biomass[p] = smooth_max(T(0.0025), leafc[p]) *
                    c_to_b * T(1.0e-3) / (one(T) - fbw[ivt[p] + 1])
            end

            if woody[ivt[p] + 1] == one(T)

                # Trees and shrubs: simple allometry with hard-wired stem taper
                # and nstem from PFT parameter file. (CNDV and non-CNDV use the
                # same standard allometry here.)
                if use_cndv
                    htop[p] = ((T(3.0) * deadstemc[p] * spinup_factor_deadwood * taper[ivt[p] + 1] * taper[ivt[p] + 1]) /
                        (T(pi) * nstem[ivt[p] + 1] * dwood[ivt[p] + 1]))^(T(1) / T(3))
                else
                    htop[p] = ((T(3.0) * deadstemc[p] * spinup_factor_deadwood * taper[ivt[p] + 1] * taper[ivt[p] + 1]) /
                        (T(pi) * nstem[ivt[p] + 1] * dwood[ivt[p] + 1]))^(T(1) / T(3))
                end

                if use_biomass_heat_storage
                    stem_biomass[p] = (spinup_factor_deadwood * deadstemc[p] + livestemc[p]) *
                        c_to_b * T(1.0e-3) / (one(T) - fbw[ivt[p] + 1])
                end

                # Keep htop from getting too close to forcing height for windspeed
                htop[p] = smooth_min(htop[p], (forc_hgt_u_patch[p] / (displar[ivt[p] + 1] + z0mr[ivt[p] + 1])) - T(3.0))

                # Constraint to keep htop from going to 0.0
                htop[p] = smooth_max(htop[p], T(0.01))

                hbot[p] = smooth_max(zero(T), smooth_min(T(3.0), htop[p] - one(T)))

            elseif ivt[p] >= npcropmin  # prognostic crops

                if tlai[p] >= laimx[ivt[p] + 1]
                    peaklai[p] = 1  # used in CNAllocation
                end

                if ivt[p] == ntmp_corn || ivt[p] == nirrig_tmp_corn ||
                   ivt[p] == ntrp_corn || ivt[p] == nirrig_trp_corn ||
                   ivt[p] == nsugarcane || ivt[p] == nirrig_sugarcane ||
                   ivt[p] == nmiscanthus || ivt[p] == nirrig_miscanthus ||
                   ivt[p] == nswitchgrass || ivt[p] == nirrig_switchgrass
                    tsai[p] = T(0.1) * tlai[p]
                else
                    tsai[p] = T(0.2) * tlai[p]
                end

                # "stubble" after harvest
                if harvdate[p] < 999 && tlai[p] == zero(T)
                    tsai[p] = T(0.25) * (one(T) - farea_burned[c] * T(0.90))
                    htmx[p] = zero(T)
                    peaklai[p] = 0
                end

                # canopy top and bottom heights
                htop[p] = ztopmx[ivt[p] + 1] * (smooth_min(tlai[p] / (laimx[ivt[p] + 1] - one(T)), one(T)))^2
                htmx[p] = smooth_max(htmx[p], htop[p])
                htop[p] = smooth_max(T(0.05), smooth_max(htmx[p], htop[p]))
                hbot[p] = T(0.02)

            else  # generic crops and grasses

                # height for grasses depends only on LAI
                htop[p] = smooth_max(T(0.25), tlai[p] * T(0.25))

                htop[p] = smooth_min(htop[p], (forc_hgt_u_patch[p] / (displar[ivt[p] + 1] + z0mr[ivt[p] + 1])) - T(3.0))

                # Constraint to keep htop from going to 0.0
                htop[p] = smooth_max(htop[p], T(0.01))

                hbot[p] = smooth_max(zero(T), smooth_min(T(0.05), htop[p] - T(0.20)))
            end

        else
            # noveg
            tlai[p] = zero(T)
            tsai[p] = zero(T)
            htop[p] = zero(T)
            hbot[p] = zero(T)
        end

        # Adjust lai and sai for burying by snow (Lombardozzi et al. 2018,
        # GRL 45(18), 9889-9897). NOTE: duplicated in SatellitePhenologyMod.
        if ivt[p] > noveg && ivt[p] <= nbrdlf_dcd_brl_shrub
            ol = smooth_min(smooth_max(snow_depth[c] - hbot[p], zero(T)), htop[p] - hbot[p])
            fb = one(T) - ol / smooth_max(T(1.0e-06), htop[p] - hbot[p])
        else
            fb = one(T) - (smooth_max(smooth_min(snow_depth[c], smooth_max(T(0.05), htop[p] * T(0.8))), zero(T)) /
                        (smooth_max(T(0.05), htop[p] * T(0.8))))
            # depth of snow required for complete burial of grasses
        end

        if frac_sno[c] <= T(FRAC_SNO_THRESHOLD)
            frac_sno_adjusted = frac_sno[c]
        else
            # avoid tiny but non-zero elai and esai that can blow up radiation/photosynthesis
            frac_sno_adjusted = one(T)
        end

        elai[p] = max(tlai[p] * (one(T) - frac_sno_adjusted) + tlai[p] * fb * frac_sno_adjusted, zero(T))
        esai[p] = max(tsai[p] * (one(T) - frac_sno_adjusted) + tsai[p] * fb * frac_sno_adjusted, zero(T))

        # Fraction of vegetation free of snow
        if (elai[p] + esai[p]) > zero(T)
            frac_veg_nosno_alb[p] = 1
        else
            frac_veg_nosno_alb[p] = 0
        end
    end
end

# ---------------------------------------------------------------------------
# cn_veg_struct_update! -- Main vegetation structure update
# ---------------------------------------------------------------------------

"""
    cn_veg_struct_update!(mask_soilp, bounds, patch, canopystate,
        cnveg_carbonstate, waterdiagnosticbulk, frictionvel, cnveg_state,
        crop, pftcon_data;
        dt, use_cndv, use_biomass_heat_storage, spinup_factor_deadwood,
        noveg, nc3crop, nc3irrig, nbrdlf_evr_shrub, nbrdlf_dcd_brl_shrub,
        npcropmin, ntmp_corn, nirrig_tmp_corn, ntrp_corn, nirrig_trp_corn,
        nsugarcane, nirrig_sugarcane, nmiscanthus, nirrig_miscanthus,
        nswitchgrass, nirrig_switchgrass, c_to_b)

On the radiation time step, use C state variables and epc to diagnose
vegetation structure (LAI, SAI, height).

Updates `canopystate` fields: `tlai_patch`, `tsai_patch`, `htop_patch`,
`hbot_patch`, `elai_patch`, `esai_patch`, `frac_veg_nosno_alb_patch`,
`stem_biomass_patch`, `leaf_biomass_patch`.

Updates `cnveg_state` fields: `htmx_patch`, `peaklai_patch`.

Ported from `CNVegStructUpdate` in `CNVegStructUpdateMod.F90`. Runs as one
KernelAbstractions kernel (one thread per patch); the CPU path is byte-identical.
"""
function cn_veg_struct_update!(mask_soilp::AbstractVector{Bool},
                                bounds::UnitRange{Int},
                                patch::PatchData,
                                canopystate::CanopyStateData,
                                cnveg_carbonstate::CNVegCarbonStateData,
                                waterdiagnosticbulk::WaterDiagnosticBulkData,
                                frictionvel::FrictionVelocityData,
                                cnveg_state::CNVegStateData,
                                crop::CropData,
                                pftcon_data::PftconType;
                                dt::Real = 1800.0,
                                use_cndv::Bool = false,
                                use_biomass_heat_storage::Bool = false,
                                spinup_factor_deadwood::Real = SPINUP_FACTOR_DEADWOOD_DEFAULT,
                                noveg::Int = CLM.noveg,
                                nc3crop::Int = CLM.nc3crop,
                                nc3irrig::Int = CLM.nc3irrig,
                                nbrdlf_evr_shrub::Int = CLM.nbrdlf_evr_shrub,
                                nbrdlf_dcd_brl_shrub::Int = CLM.nbrdlf_dcd_brl_shrub,
                                npcropmin::Int = CLM.npcropmin,
                                ntmp_corn::Int = CLM.ntmp_corn,
                                nirrig_tmp_corn::Int = CLM.nirrig_tmp_corn,
                                ntrp_corn::Int = CLM.ntrp_corn,
                                nirrig_trp_corn::Int = CLM.nirrig_trp_corn,
                                nsugarcane::Int = CLM.nsugarcane,
                                nirrig_sugarcane::Int = CLM.nirrig_sugarcane,
                                nmiscanthus::Int = CLM.nmiscanthus,
                                nirrig_miscanthus::Int = CLM.nirrig_miscanthus,
                                nswitchgrass::Int = CLM.nswitchgrass,
                                nirrig_switchgrass::Int = CLM.nirrig_switchgrass,
                                c_to_b::Real = C_TO_B)

    FT = eltype(canopystate.tlai_patch)

    idx = _VSUIdx(; ivt = patch.itype, col = patch.column, harvdate = crop.harvdate_patch)
    inb = _VSUIn(; leafc = cnveg_carbonstate.leafc_patch,
                   deadstemc = cnveg_carbonstate.deadstemc_patch,
                   livestemc = cnveg_carbonstate.livestemc_patch,
                   forc_hgt_u = frictionvel.forc_hgt_u_patch,
                   frac_sno = waterdiagnosticbulk.frac_sno_col,
                   snow_depth = waterdiagnosticbulk.snow_depth_col,
                   farea_burned = cnveg_state.farea_burned_col)
    out = _VSUOut(; tlai = canopystate.tlai_patch, tsai = canopystate.tsai_patch,
                    htop = canopystate.htop_patch, hbot = canopystate.hbot_patch,
                    elai = canopystate.elai_patch, esai = canopystate.esai_patch,
                    stem_biomass = canopystate.stem_biomass_patch,
                    leaf_biomass = canopystate.leaf_biomass_patch,
                    htmx = cnveg_state.htmx_patch,
                    frac_veg_nosno_alb = canopystate.frac_veg_nosno_alb_patch,
                    peaklai = cnveg_state.peaklai_patch)
    # PFT-parameter vectors: host-global pftcon → move onto the state backend +
    # precision (no-op on the host path).
    tl = canopystate.tlai_patch
    par = _VSUPar(; woody = _to_backend_like(tl, FT, pftcon_data.woody),
                    slatop = _to_backend_like(tl, FT, pftcon_data.slatop),
                    dsladlai = _to_backend_like(tl, FT, pftcon_data.dsladlai),
                    z0mr = _to_backend_like(tl, FT, pftcon_data.z0mr),
                    displar = _to_backend_like(tl, FT, pftcon_data.displar),
                    dwood = _to_backend_like(tl, FT, pftcon_data.dwood),
                    ztopmx = _to_backend_like(tl, FT, pftcon_data.ztopmx),
                    laimx = _to_backend_like(tl, FT, pftcon_data.laimx),
                    nstem = _to_backend_like(tl, FT, pftcon_data.nstem),
                    taper = _to_backend_like(tl, FT, pftcon_data.taper),
                    fbw = _to_backend_like(tl, FT, pftcon_data.fbw))
    sc = _VSUScalars{FT}(FT(dt), FT(spinup_factor_deadwood), FT(c_to_b),
                         use_cndv, use_biomass_heat_storage,
                         noveg, nc3crop, nc3irrig, nbrdlf_dcd_brl_shrub, npcropmin,
                         ntmp_corn, nirrig_tmp_corn, ntrp_corn, nirrig_trp_corn,
                         nsugarcane, nirrig_sugarcane, nmiscanthus, nirrig_miscanthus,
                         nswitchgrass, nirrig_switchgrass)

    _launch!(_vsu_kernel!, mask_soilp, idx, inb, out, par, sc)

    return nothing
end
