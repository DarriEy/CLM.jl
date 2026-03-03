# ==========================================================================
# Ported from: src/main/controlMod.F90 (~1350 lines)
# Control module — simulation control settings
#
# Skipped: all namelist reading (control_setNL, NAMELIST I/O),
#           MPI broadcast (control_spmd), file-existence checks
#           (check_missing_initdata_status, apply_use_init_interp).
# Ported:   run-type labels, control_init! (defaults + consistency),
#           control_print.
# ==========================================================================

# Run-type labels (Fortran: runtyp(4), indexed nsrest+1)
const RUNTYP = ["initial", "restart", "branch", "missing"]

"""
    control_init!(ctl::VarCtl)

Apply derived defaults and perform consistency checks on control variables.
Mirrors the non-I/O parts of Fortran `control_init`.
"""
function control_init!(ctl::VarCtl)

    # --- FATES overrides for matrix CN ---
    if ctl.use_fates
        ctl.spinup_matrixcn        = false
        ctl.hist_wrt_matrixcn_diag = false
    end

    # --- FATES BGC vs SP ---
    if ctl.use_fates
        ctl.use_fates_bgc = !ctl.use_fates_sp
    else
        ctl.use_fates_sp  = false
        ctl.use_fates_bgc = false
    end

    # --- Anoxia tied to methane submodel ---
    if ctl.use_lch4
        ctl.anoxia = true
    else
        ctl.anoxia = false
    end

    # --- nfix_timeconst default ---
    if ctl.nfix_timeconst == -1.2345
        ctl.nfix_timeconst = ctl.use_nitrif_denitrif ? 10.0 : 0.0
    end

    # --- Consistency checks ---

    # co2_type
    if !(ctl.co2_type in ("constant", "prognostic", "diagnostic"))
        error("co2_type = '$(ctl.co2_type)' is not supported; choices are constant, prognostic, or diagnostic")
    end

    # nsrest must be set
    if ctl.nsrest == IUNDEF
        error("must set nsrest")
    end

    # Branch requires restart file
    if ctl.nsrest == NSRBRANCH && isempty(strip(ctl.nrevsn))
        error("need to set restart data file name for branch run")
    end

    # co2_ppmv range
    if ctl.co2_ppmv <= 0.0 || ctl.co2_ppmv > 3000.0
        error("co2_ppmv = $(ctl.co2_ppmv) is out of a reasonable range")
    end

    # SNICAR internal mixing exclusivity
    if ctl.snicar_snobc_intmix && ctl.snicar_snodst_intmix
        error("currently dust-snow and BC-snow internal mixing cannot be activated together")
    end

    # nsrest validity
    if !(ctl.nsrest in (NSRSTARTUP, NSRCONTINUE, NSRBRANCH))
        error("nsrest = $(ctl.nsrest) is NOT set to a valid value")
    end

    # nrevsn adjustments
    if ctl.nsrest == NSRSTARTUP
        ctl.nrevsn = ""
    end
    if ctl.nsrest == NSRCONTINUE
        ctl.nrevsn = "set by restart pointer file"
    end

    # Crop requires crop landunit
    if ctl.use_crop && !ctl.create_crop_landunit
        error("prognostic crop patches require create_crop_landunit=true")
    end

    # Single column mode requires lat/lon
    if ctl.single_column && (ctl.scmlat == RUNDEF || ctl.scmlon == RUNDEF)
        error("single column mode on but scmlat and scmlon are NOT set")
    end

    # FATES incompatibilities
    if ctl.use_fates
        ctl.use_cn && error("use_cn and use_fates cannot both be true")
        ctl.use_hydrstress && error("use_hydrstress and use_fates cannot both be true")
        ctl.use_crop && error("use_crop and use_fates cannot both be true")
        ctl.use_luna && error("use_luna is not compatible with FATES")
        ctl.o3_veg_stress_method != "unset" && error("ozone is not compatible with FATES")
        (ctl.use_c13 || ctl.use_c14) && error("C13 and C14 dynamics are not compatible with FATES")
    end

    # LAI streams compatibility
    if ctl.use_lai_streams
        if (ctl.use_fates && !ctl.use_fates_sp) || ctl.use_cn
            error("cannot use LAI streams unless in SP mode (use_cn=false or use_fates_sp=true)")
        end
    end

    # Anoxia requires CH4 submodel
    if !ctl.use_lch4 && ctl.anoxia
        error("anoxia requires the CH4 submodel (use_lch4=true)")
    end

    # Dynamic glacier incompatibilities
    if ctl.glc_do_dynglacier
        if ctl.collapse_urban
            error("glc_do_dynglacier is incompatible with collapse_urban=true")
        end
        if ctl.n_dom_pfts > 0 || ctl.n_dom_landunits > 0 ||
           ctl.toosmall_soil > 0.0 || ctl.toosmall_crop > 0.0 ||
           ctl.toosmall_glacier > 0.0 || ctl.toosmall_lake > 0.0 ||
           ctl.toosmall_wetland > 0.0 || ctl.toosmall_urban > 0.0
            error("glc_do_dynglacier is incompatible with n_dom_pfts/n_dom_landunits/toosmall_* > 0")
        end
    end

    # Spinup state validity
    if !(ctl.spinup_state in (0, 1, 2))
        error("spinup_state can only be 0, 1, or 2; got $(ctl.spinup_state)")
    end

    return nothing
end

"""
    control_print(ctl::VarCtl; io::IO=stdout)

Print run control information.  Mirrors Fortran `control_print`.
"""
function control_print(ctl::VarCtl; io::IO=stdout)
    runlabel = (0 <= ctl.nsrest <= 2) ? RUNTYP[ctl.nsrest + 1] : "unknown"

    println(io, "define run:")
    println(io, "   source                = ", ctl.source)
    println(io, "   model_version         = ", ctl.version)
    println(io, "   run type              = ", runlabel)
    println(io, "   case title            = ", ctl.ctitle)
    println(io, "   username              = ", ctl.username)
    println(io, "   hostname              = ", ctl.hostname)

    println(io, "process control parameters:")
    println(io, "    use_lch4 = ", ctl.use_lch4)
    println(io, "    use_nitrif_denitrif = ", ctl.use_nitrif_denitrif)
    println(io, "    use_extralakelayers = ", ctl.use_extralakelayers)
    println(io, "    use_vichydro = ", ctl.use_vichydro)
    println(io, "    use_excess_ice = ", ctl.use_excess_ice)
    println(io, "    use_cn = ", ctl.use_cn)
    println(io, "    use_cndv = ", ctl.use_cndv)
    println(io, "    use_crop = ", ctl.use_crop)
    println(io, "    flush_gdd20 = ", ctl.flush_gdd20)
    println(io, "    use_fertilizer = ", ctl.use_fertilizer)
    println(io, "    use_grainproduct = ", ctl.use_grainproduct)
    println(io, "    crop_residue_removal_frac = ", ctl.crop_residue_removal_frac)
    println(io, "    o3_veg_stress_method = ", ctl.o3_veg_stress_method)
    println(io, "    use_snicar_frc = ", ctl.use_snicar_frc)
    println(io, "    snicar_use_aerosol = ", ctl.snicar_use_aerosol)
    println(io, "    use_vancouver = ", ctl.use_vancouver)
    println(io, "    use_mexicocity = ", ctl.use_mexicocity)
    println(io, "    use_noio = ", ctl.use_noio)
    println(io, "    use_SSRE = ", ctl.use_SSRE)

    println(io, "input data files:")
    println(io, "   PFT physiology and parameters file = ", ctl.paramfile)
    if isempty(strip(ctl.fsurdat))
        println(io, "   fsurdat, surface dataset not set")
    else
        println(io, "   surface data   = ", ctl.fsurdat)
    end

    println(io, "   Convert ocean to land = ", ctl.convert_ocean_to_land)
    println(io, "   n_dom_pfts = ", ctl.n_dom_pfts)
    println(io, "   n_dom_landunits = ", ctl.n_dom_landunits)
    println(io, "   collapse_urban = ", ctl.collapse_urban)
    println(io, "   toosmall_soil = ", ctl.toosmall_soil)
    println(io, "   toosmall_crop = ", ctl.toosmall_crop)
    println(io, "   toosmall_glacier = ", ctl.toosmall_glacier)
    println(io, "   toosmall_lake = ", ctl.toosmall_lake)
    println(io, "   toosmall_wetland = ", ctl.toosmall_wetland)
    println(io, "   toosmall_urban = ", ctl.toosmall_urban)

    if ctl.use_cn || ctl.use_fates
        println(io, "   spinup_state = ", ctl.spinup_state)
        println(io, "   override_bgc_restart_mismatch_dump = ", ctl.override_bgc_restart_mismatch_dump)
    end

    println(io, "   co2_type = ", ctl.co2_type)
    if ctl.co2_type == "constant"
        println(io, "   CO2 volume mixing ratio (umol/mol) = ", ctl.co2_ppmv)
    end

    println(io, "   use_hillslope = ", ctl.use_hillslope)
    println(io, "   downscale_hillslope_meteorology = ", ctl.downscale_hillslope_meteorology)
    println(io, "   use_hillslope_routing = ", ctl.use_hillslope_routing)
    println(io, "   hillslope_fsat_equals_zero = ", ctl.hillslope_fsat_equals_zero)

    println(io, "   soil_layerstruct_predefined = ", ctl.soil_layerstruct_predefined)
    println(io, "   use_hydrstress = ", ctl.use_hydrstress)
    println(io, "   nsegspc = ", ctl.nsegspc)

    println(io, "   use_flexibleCN = ", ctl.use_flexibleCN)
    println(io, "   use_luna = ", ctl.use_luna)
    println(io, "   use_fates = ", ctl.use_fates)

    if ctl.use_fates
        println(io, "   fates_spitfire_mode = ", ctl.fates_spitfire_mode)
        println(io, "   fates_harvest_mode = ", ctl.fates_harvest_mode)
        println(io, "   fates_paramfile = ", ctl.fates_paramfile)
        println(io, "   fates_parteh_mode = ", ctl.fates_parteh_mode)
        println(io, "   use_fates_planthydro = ", ctl.use_fates_planthydro)
        println(io, "   use_fates_tree_damage = ", ctl.use_fates_tree_damage)
        println(io, "   use_fates_cohort_age_tracking = ", ctl.use_fates_cohort_age_tracking)
        println(io, "   use_fates_ed_st3 = ", ctl.use_fates_ed_st3)
        println(io, "   use_fates_ed_prescribed_phys = ", ctl.use_fates_ed_prescribed_phys)
        println(io, "   use_fates_inventory_init = ", ctl.use_fates_inventory_init)
        println(io, "   use_fates_fixed_biogeog = ", ctl.use_fates_fixed_biogeog)
        println(io, "   use_fates_nocomp = ", ctl.use_fates_nocomp)
        println(io, "   use_fates_sp = ", ctl.use_fates_sp)
        println(io, "   use_fates_luh = ", ctl.use_fates_luh)
        println(io, "   use_fates_lupft = ", ctl.use_fates_lupft)
        println(io, "   use_fates_potentialveg = ", ctl.use_fates_potentialveg)
        println(io, "   fates_seeddisp_cadence = ", ctl.fates_seeddisp_cadence)
    end

    return nothing
end
