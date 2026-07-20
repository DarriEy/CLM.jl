# ==========================================================================
# Drift protection for the driver / control-flag defaults audit.
#
# See docs/DRIVER_DEFAULTS_AUDIT.md for the full sweep and the CTSM evidence
# behind every value pinned here.
#
# WHY THIS FILE EXISTS
# --------------------
# Two merged bugs (#252 use_bedrock, #225 h2osfcflag) were wrong DEFAULTS. Both
# were invisible to the whole test suite, because every test and harness that
# cared passed an explicit value — the wrong default only bit callers who
# trusted it. No conservation check, parity scorecard or physics test can catch
# this class: they all construct their own config.
#
# The only thing that catches a wrong default is an assertion ON THE DEFAULT
# ITSELF. That is what this file is. Each assertion carries the CTSM source
# location that justifies it, so a future change has to argue with the evidence
# rather than silently re-tune.
#
# Resolution conditions for this port: phys="clm5_0", use_fates=.false.,
# configuration="clm", structure="standard", vichydro=0.
# ==========================================================================

using Test
using CLM

# --- Extract a function's keyword-argument defaults straight from the AST. ---
# Parsing the source (rather than hardcoding a mirror list) means this test
# tracks the real signature and cannot drift out of sync with it.
function _kwarg_defaults(file::String, fname::Symbol)
    src = read(file, String)
    defaults = Dict{Symbol,Any}()
    for expr in Meta.parseall(src).args
        found = Ref(false)
        walk(e) = begin
            if e isa Expr
                if e.head === :function && length(e.args) >= 1
                    sig = e.args[1]
                    if sig isa Expr && sig.head === :call && !isempty(sig.args)
                        name = sig.args[1]
                        if name === fname && length(sig.args) >= 2 &&
                           sig.args[2] isa Expr && sig.args[2].head === :parameters
                            found[] = true
                            for kw in sig.args[2].args
                                if kw isa Expr && kw.head === :kw
                                    lhs = kw.args[1]
                                    nm = lhs isa Expr ? lhs.args[1] : lhs
                                    defaults[nm] = kw.args[2]
                                end
                            end
                        end
                    end
                end
                foreach(walk, e.args)
            end
        end
        walk(expr)
        found[] && break
    end
    return defaults
end

const _SRC = normpath(joinpath(@__DIR__, "..", "src"))
const _INIT_KW = _kwarg_defaults(joinpath(_SRC, "driver", "clm_initialize.jl"), :clm_initialize!)
const _RUN_KW  = _kwarg_defaults(joinpath(_SRC, "driver", "clm_run.jl"), :clm_run!)

@testset "driver defaults audit (docs/DRIVER_DEFAULTS_AUDIT.md)" begin

    @testset "AST extraction actually found the signatures" begin
        # Guard against the vacuous-green failure mode: if the parse silently
        # returned nothing, every assertion below would be trivially skipped.
        @test !isempty(_INIT_KW)
        @test !isempty(_RUN_KW)
        @test haskey(_INIT_KW, :h2osfcflag)
        @test haskey(_RUN_KW, :h2osfcflag)
    end

    @testset "h2osfcflag == 1 (fixed in #225)" begin
        # CTSM code fallback: SoilHydrologyType.F90:339,347 `h2osfcflag = 1`
        # (clm_soilhydrology_inparm). 0 silently disables the surface-water
        # store: infiltration excess bypasses ponding/re-infiltration and
        # frac_h2osfc pins at 0.
        @test _INIT_KW[:h2osfcflag] == 1
        @test _RUN_KW[:h2osfcflag] == 1
    end

    @testset "use_aquifer_layer == false (fixed in the defaults campaign)" begin
        # CTSM does not read this from a namelist — it DERIVES it:
        #   namelist_defaults_ctsm.xml:419-428 → soilwater_movement_method=1
        #   (clm5_0) → lower_boundary_condition=2 (bc_zero_flux)
        #   SoilWaterMovementMod.F90:221-236 → use_aquifer_layer() is true only
        #   for bc_aquifer(4)/bc_watertable ⇒ .false. for clm5_0.
        # The port shipped `true`, which is CTSM's CODE FALLBACK
        # (SoilWaterMovementMod.F90:~150, under the comment "Default values for
        # namelist" — i.e. the value the namelist overrides). `true` also pairs
        # with use_bedrock into a config CTSM endrun's on (same file:181).
        # It selects the whole soil-water solver (ZD09+BC_AQUIFER vs the CLM5
        # moisture-form solver + BC_ZERO_FLUX) via clm_driver.jl's swm_cfg switch.
        # 71 of the 75 explicit call sites in scripts/+test/ already passed false.
        @test _INIT_KW[:use_aquifer_layer] === false
        @test _RUN_KW[:use_aquifer_layer] === false
    end

    @testset "use_bedrock is CONDITIONAL, not a bare Bool (fixed in #252)" begin
        # namelist_defaults_ctsm.xml:178-181 — .false. under use_fates /
        # vichydro / clm4_5, .true. otherwise. The `nothing` sentinel is what
        # encodes "resolve the condition"; a literal `true` here was the bug.
        @test _INIT_KW[:use_bedrock] === :nothing
        @test _RUN_KW[:use_bedrock] === :nothing
    end

    @testset "create_crop_landunit derives from use_fates (M5, FIXED)" begin
        # CTSM keys this on use_fates ALONE
        # (namelist_defaults_ctsm.xml:2377-2378), and CLMBuildNamelist.pm:2248-2250
        # makes `.false.` a FATAL error for any non-FATES run. Deriving it from
        # use_crop gave a default non-crop run cft_size=0 and folded crop area
        # into natural veg, where CTSM builds a separate crop landunit.
        #
        # The struct field default stays `false` on purpose: that mirrors CTSM's
        # CODE fallback (clm_varctl.F90:154). It is the DERIVATION that has to
        # key on use_fates, exactly as CLMBuildNamelist does.
        @test CLM.VarCtl().create_crop_landunit == false   # clm_varctl.F90:154

        init_src = read(joinpath(_SRC, "driver", "clm_initialize.jl"), String)
        @test occursin(r"varctl\.create_crop_landunit\s*=\s*!use_fates", init_src)
        @test !occursin(r"varctl\.create_crop_landunit\s*=\s*use_crop\b", init_src)

        # varpar_init! must branch on the FLAG, not on use_crop — this is the
        # port's mirror of clm_varpar.F90:209.
        varpar_src = read(joinpath(_SRC, "constants", "varpar.jl"), String)
        @test occursin("if varctl.create_crop_landunit", varpar_src)

        # ...and behaviourally, not just textually.
        old_clu  = CLM.varctl.create_crop_landunit
        old_crop = CLM.varctl.use_crop
        try
            CLM.varctl.use_crop = false          # a NON-crop run, the default
            CLM.varctl.create_crop_landunit = true
            CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
            @test CLM.varpar.cft_size == 2       # CTSM: cft_size = surf_numcft
            @test CLM.varpar.cft_lb == 15
            @test CLM.varpar.cft_ub == 16

            # The FATES branch (create_crop_landunit=.false.) still collapses.
            CLM.varctl.create_crop_landunit = false
            CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
            @test CLM.varpar.cft_size == 0
        finally
            CLM.varctl.create_crop_landunit = old_clu
            CLM.varctl.use_crop = old_crop
        end
    end

    @testset "use_luna is CONDITIONAL, not a bare Bool (row M4)" begin
        # namelist_defaults_ctsm.xml:578-580 — .true. generally, .false. under
        # use_fates and under phys=clm4_5. The port shipped `false`, which is the
        # CODE FALLBACK (clm_varctl.F90:371) — same root cause as #252/#259.
        # A bare `true` would be wrong too: CTSM endrun's on LUNA+FATES
        # (controlMod.F90:505), mirrored at control.jl:104. Hence `nothing`.
        @test _INIT_KW[:use_luna] === :nothing
        @test _RUN_KW[:use_luna] === :nothing

        # The conditional must actually RESOLVE — pin the resolved values, not
        # just the sentinel, so this cannot go vacuously green.
        @test CLM.CLMDriverConfig().use_luna === true
        @test CLM.CLMDriverConfig(use_fates=true).use_luna === false
        @test CLM.CLMDriverConfig(use_luna=false).use_luna === false
    end

    @testset "use_hydrstress is CONDITIONAL, not a bare Bool (row M3)" begin
        # namelist_defaults_ctsm.xml — .true. for phys=clm5_0, use_fates=.false.,
        # configuration="clm"; .false. otherwise. The port shipped `false`, which is
        # CTSM's CODE FALLBACK (clm_varctl.F90) — the same root cause as #252/#259/#267.
        # A bare `true` would be wrong too: control.jl:102 endruns on PHS+FATES,
        # mirroring CTSM. Hence `nothing`. See docs/PHS_DEFAULT_BLOCKERS.md.
        @test _INIT_KW[:use_hydrstress] === :nothing
        @test _RUN_KW[:use_hydrstress] === :nothing

        # The conditional must actually RESOLVE — pin the resolved values, not just
        # the sentinel, so this cannot go vacuously green (the #267 lesson).
        @test CLM.CLMDriverConfig().use_hydrstress === true
        @test CLM.CLMDriverConfig(use_fates=true).use_hydrstress === false
        @test CLM.CLMDriverConfig(use_hydrstress=false).use_hydrstress === false

        # The varctl STRUCT default stays `false`: that mirrors CTSM's code fallback,
        # the value the namelist overrides. It is the DERIVATION that carries the
        # namelist default — the #265 principle.
        @test CLM.VarCtl().use_hydrstress == false
    end

    @testset "baseflow_scalar == 0.001 — and all FOUR literals agree (row M2)" begin
        # namelist_defaults_ctsm.xml:195-197 keys this on lower_boundary_condition:
        #   lbc=2 -> 0.001 | lbc=1 -> 1.d-2 | phys=clm4_5 + lbc=2 -> 1.d-2
        # For clm5_0 the chain is soilwater_movement_method=1 (defaults:419) +
        # use_bedrock=.true. (defaults:180) => lbc=2 (defaults:426) => 0.001.
        # The port shipped 1.0e-2 — CTSM's CODE FALLBACK
        # (SoilHydrologyMod.F90:71 `= 1.e-2_r8`, the pre-namelist initialiser),
        # the identical trap as #252/#259/#273. Corroborated by two
        # CTSM-emitted lnd_in files (Bow, MerBleue): both `0.001d00`.
        @test _RUN_KW[:baseflow_scalar] == 0.001

        # WHY ALL FOUR, NOT JUST THE ONE PHYSICS READS.
        # Only BASEFLOW_SCALAR[] is read by the drainage kernels, but
        # test_soil_hydrology_mod.jl calls a bare init_soil_hydrology_config()
        # to "reset to defaults". If the struct/kwarg defaults disagreed with the
        # Ref, that call would silently install the OTHER value into the live
        # global — making the physics depend on test execution ORDER. Pin the
        # agreement itself, not just the value.
        @test CLM.BASEFLOW_SCALAR[] == 0.001            # module Ref (live)
        @test CLM.SoilHydrologyConfig().baseflow_scalar == 0.001
        @test CLM.init_soil_hydrology_config().baseflow_scalar == 0.001
        @test CLM.BASEFLOW_SCALAR[] == 0.001            # ...and it stayed put

        # The calibration path already used the CTSM value; it must not drift
        # back apart from the driver (they disagreed for the port's whole life
        # until this row was closed).
        hydro = filter(p -> p.name == "baseflow_scalar", CLM.default_hydrology_params())
        @test length(hydro) == 1
        @test hydro[1].default_value == 0.001
    end

    @testset "int_snow_max == 2000.0" begin
        # namelist_defaults_ctsm.xml:476 (clm4_5 variant is 1.e30).
        @test _INIT_KW[:int_snow_max] == 2000.0
        @test _RUN_KW[:int_snow_max] == 2000.0
    end

    @testset "soil_layerstruct == 20SL_8.5m (clm5_0)" begin
        @test _INIT_KW[:soil_layerstruct] == "20SL_8.5m"
    end

    # ---------------------------------------------------------------------
    # varctl defaults. These are pinned as-declared; the ones that DISAGREE
    # with CTSM are pinned too, with the disagreement named, so that the
    # audit's findings are visible from the test suite rather than only from a
    # markdown file. If someone fixes one of these, this test tells them
    # exactly which audit row they are closing.
    # ---------------------------------------------------------------------
    @testset "varctl — values that MATCH CTSM" begin
        v = CLM.VarCtl()
        @test v.co2_type == "constant"                       # defaults:33
        @test v.use_subgrid_fluxes == true                   # defaults:639
        @test v.snicar_solarspec == "mid_latitude_winter"    # defaults:2051
        @test v.snicar_dust_optics == "sahara"               # defaults:2052
        @test v.snicar_use_aerosol == true                   # defaults:2057
        @test v.snicar_snobc_intmix == false                 # defaults:2054 (clm5_1+ → true)
        @test v.snicar_snodst_intmix == false                # defaults:2053
        @test v.snicar_numrad_snw == 5                       # defaults:2041
        @test v.snicar_snw_shape == "sphere"                 # defaults:2047 (phys=clm5_0)
        @test v.snow_thermal_cond_method == "Jordan1991"     # defaults:503 (clm5_1+ → Sturm1997)
        @test v.do_sno_oc == false                           # defaults:2058
        @test v.use_snicar_frc == false                      # defaults:2050
        @test v.use_biomass_heat_storage == false            # defaults:407 (phys=clm5_0)
        @test v.use_excess_ice == false                      # defaults:2597 (clm6_0 → true)
        @test v.collapse_urban == false                      # defaults:2384 (structure=standard)
        @test v.run_zero_weight_urban == false               # defaults:2380
        @test v.crop_fsat_equals_zero == false               # defaults:137
        @test v.crop_residue_removal_frac == 0.0             # defaults:609 (clm5_1+ → 0.5)
        @test v.spinup_state == 0                            # defaults:56 (accelerated_spinup=off)
        @test v.use_hillslope == false                       # defaults:657
        @test v.use_hillslope_routing == false               # defaults:660
        @test v.nyr_SASU == 1                                # defaults:690
        @test v.iloop_avg == -999                            # defaults:692
        @test v.o3_veg_stress_method == "unset"              # defaults:603
        @test v.use_soil_moisture_streams == false           # defaults:2115
        @test v.use_lai_streams == false                     # defaults:2126
        @test v.use_z0m_snowmelt == false                    # defaults:568
        @test v.irrigate == false                            # defaults:132 (sim_year_range=constant)
    end

    @testset "varctl — not namelist variables in CTSM (code fallback is authority)" begin
        v = CLM.VarCtl()
        @test v.outnc_large_files == true      # clm_varctl.F90:104
        @test v.ndep_from_cpl == false         # clm_varctl.F90:125
        @test v.bound_h2osoi == true           # clm_varctl.F90:131
        @test v.nhillslope == 0                # clm_varctl.F90:157
        @test v.max_columns_hillslope == 1     # clm_varctl.F90:160
        @test v.anoxia == true                 # clm_varctl.F90:210
        @test v.o3_ppbv == 100.0               # clm_varctl.F90:248
        @test v.nfix_timeconst == -1.2345      # clm_varctl.F90:228
    end

    @testset "varctl — KNOWN MISMATCHES vs CTSM (pinned so a fix is deliberate)" begin
        v = CLM.VarCtl()
        # Each of these is a real disagreement with CTSM's NAMELIST default.
        # CLM.jl copied the Fortran CODE FALLBACK in every one — the same root
        # cause as #252. They are pinned at the current (wrong) value so that
        # closing one is an explicit, reviewed act rather than an accident.
        # Severity is recorded in docs/DRIVER_DEFAULTS_AUDIT.md.

        # LIVE (read on a real path) — highest priority to fix:
        @test v.convert_ocean_to_land == false   # CTSM defaults:2382 → .true.  (live: surfdata.jl:306)
        @test v.co2_ppmv == 355.0                # CTSM defaults:25   → 379.0 at sim_year=2000
        @test v.glc_snow_persistence_max_days == 7300  # CTSM defaults:540 → 0; 7300 is the clm4_5 variant

        # INERT tranche — CLOSED to the CTSM values (defaults campaign #2).
        # Inertness was PROVEN, not asserted: with these twelve flipped, the full
        # suite moved by exactly the twelve pin assertions below and nothing else
        # (26150→26140 passed with 12 failing pins, +2 from a new calibration
        # test). Zero physics assertions moved. These now pin the CORRECT value.
        @test v.nsegspc == 35                    # CTSM defaults:2404 (code fallback 20)
        @test v.nyr_forcing == 1                 # CTSM defaults:688  (code fallback 10)
        @test v.h2osno_max == 10000.0            # CTSM defaults:465; controlMod.F90:557 HARD-ERRORS on <=0, so -999.0 was a value CTSM refuses to start on. This varctl copy is a dead duplicate — physics reads varcon.H2OSNO_MAX.
        @test v.n_dom_landunits == 0             # CTSM defaults:2387 (guards are `> 0`, so -1 ≡ 0)
        @test v.n_dom_pfts == 0                  # CTSM defaults:2390
        @test v.toosmall_soil == 0.0             # CTSM defaults:2393
        @test v.toosmall_crop == 0.0             # CTSM defaults:2394
        @test v.toosmall_glacier == 0.0          # CTSM defaults:2395
        @test v.toosmall_lake == 0.0             # CTSM defaults:2396
        @test v.toosmall_wetland == 0.0          # CTSM defaults:2397
        @test v.toosmall_urban == 0.0            # CTSM defaults:2398
        @test v.downscale_hillslope_meteorology == true   # CTSM defaults:665 (inert: hillslope off)
        @test v.z0param_method == "ZengWang2007" # CTSM clm5_0; cosmetic — friction_velocity.jl:322 branches only on "Meier2022"
    end

    @testset "varctl — the use_flexibleCN cascade is OFF (open finding M6)" begin
        # Under phys=clm5_0 WITH use_cn=.true., CTSM turns use_flexibleCN on,
        # which cascades to all of the below. CLM.jl hardcodes the SP-mode
        # values. Pinned so that enabling flexible-CN is a deliberate campaign
        # with its own Fortran validation, not a silent drift.
        v = CLM.VarCtl()
        @test v.use_flexibleCN == false             # CTSM (clm5_0 + use_cn) → .true.
        @test v.MM_Nuptake_opt == false             # → .true.
        @test v.CNratio_floating == false           # → .true.
        @test v.lnc_opt == false                    # → .true.
        @test v.vcmax_opt == 0                      # → 3
        @test v.CN_evergreen_phenology_opt == 0     # → 1
        @test v.carbon_resp_opt == 0                # → 1 (0 if FUN on)
        @test v.use_nguardrail == false             # → .true.
    end
end
