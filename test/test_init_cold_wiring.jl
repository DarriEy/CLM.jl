# ==========================================================================
# test_init_cold_wiring.jl — guard the "ported InitCold, never called" bug CLASS.
#
# Fortran: Init = InitAllocate (NaN-fills) + InitHistory + InitCold (zeroes/seeds).
# CLM.jl faithfully ported ~20 `*_init_cold!` routines — and then never called
# most of them, so the NaN-fill survived into the run. Three separate live
# failures were traced to that single omission before it was recognised as a
# class (waterflux -> qflx_evap_tot NaN on every urban roof/wall column, a live
# lnd2atm output; aerosol -> SNICAR albedo NaN -> the glacier band solve NaN'd;
# activelayer -> altmax never tracked).
#
# The unit tests for those routines all PASSED the whole time — each called its
# own `*_init_cold!` directly. That is exactly why the ports looked tested and
# were dead. So this file asserts the WIRING, not the function:
#
#   1. STATIC  — every `*_init_cold!` defined anywhere in src/ has a real call
#                site on the LIVE INIT PATH. A newly added port that nobody wires
#                fails here, immediately, with no domain-specific NaN needed.
#   2. RUNTIME — after clm_initialize! + one clm_drv! step, the fields those
#                routines OWN are finite on every ACTIVE column / patch, and the
#                total NaN-on-active count stays under a per-domain ceiling.
#
# See src/infrastructure/init_cold.jl for the dispatcher and the full rationale
# (including which fields Fortran DELIBERATELY leaves at NaN/spval as a
# missing-value flag — those are not zeroed, and are excluded here).
# ==========================================================================
using Test, CLM, Dates

const SRC_DIR = normpath(joinpath(@__DIR__, "..", "src"))

# Files that constitute the LIVE INIT PATH. A `*_init_cold!` is "alive" iff it is
# called from one of these. (A call from test/ or scripts/ does NOT count — that
# is precisely the trap this file exists to catch.)
const LIVE_INIT_FILES = [
    joinpath(SRC_DIR, "infrastructure", "init_cold.jl"),      # the dispatcher
    joinpath(SRC_DIR, "infrastructure", "cold_start.jl"),     # cold_start_initialize!
    joinpath(SRC_DIR, "infrastructure", "instances.jl"),      # clm_instInit!
    joinpath(SRC_DIR, "driver", "clm_initialize.jl"),         # clm_initialize!
    joinpath(SRC_DIR, "driver", "clm_driver.jl"),             # irrigation re-init hook
]
# NOTE: src/infrastructure/topo.jl is deliberately NOT in that list. It defines
# `topo_init!` = allocate + history + `topo_init_cold!`, which LOOKS like a live
# call site — but `clm_instInit!` calls `topo_init_allocate!` DIRECTLY and never
# calls `topo_init!`. Including topo.jl would false-green a routine that is
# dead-in-effect behind an alive-looking wrapper. See EXPECTED_DEAD below.

# Routines that are genuinely dead and are STAYING dead for a documented reason.
# Adding to this set requires the reason; it is not an escape hatch.
const EXPECTED_DEAD = Dict(
    "topo_init_cold!" =>
    """
    DEAD-IN-EFFECT and deliberately left so — wiring it would be a REGRESSION.

    `topo_init_allocate!` fills topo_col with NaN; `topo_init_cold!` would set it
    to 0.0 (non-ice) / topo_glc_mec (ice), per TopoMod.F90::InitCold. But CLM.jl's
    `downscale_forcings!` gates on `isnan(topo_col)` as its SKIP-DOWNSCALING
    sentinel, whereas Fortran gates on `needs_downscaling_col`. So seeding
    topo_col = 0.0 would silently switch on elevation downscaling FROM 0 m for
    every column whose caller does not set topo_col per step — a large, wrong
    forcing correction.

    Today no NaN reaches physics: the live run harness sets topo_col from the
    forcing every run (clm_run.jl:277), and the domain smokes set it per step.

    Wiring this properly = first change `downscale_forcings!` to gate on
    `needs_downscaling_col` (the Fortran gate), then call `topo_init_cold!` from
    `init_cold_biogeophys!`. That is its own change with its own parity check, not
    a rider on the InitCold sweep.
    """,
)

"Every `*_init_cold!` FUNCTION DEFINITION under src/, as a Set of names."
function defined_init_colds()
    names = Set{String}()
    for (root, _, files) in walkdir(SRC_DIR), f in files
        endswith(f, ".jl") || continue
        for line in eachline(joinpath(root, f))
            m = match(r"^function\s+([A-Za-z0-9_]+_init_cold!)\s*\(", line)
            m === nothing || push!(names, m.captures[1])
        end
    end
    return names
end

"Concatenated source of the live init path, with comment lines stripped so a
mention inside a doc/comment cannot masquerade as a call site."
function live_init_source()
    io = IOBuffer()
    for f in LIVE_INIT_FILES
        isfile(f) || continue
        for line in eachline(f)
            stripped = lstrip(line)
            startswith(stripped, "#") && continue
            println(io, line)
        end
    end
    return String(take!(io))
end

@testset "InitCold wiring (the dead-port guard)" begin

    @testset "every *_init_cold! in src/ is CALLED on the live init path" begin
        defined = defined_init_colds()
        src = live_init_source()
        @test !isempty(defined)          # the scan itself must find something

        dead = String[]
        for fn in sort(collect(defined))
            haskey(EXPECTED_DEAD, fn) && continue
            # a call site is `name(` appearing outside a comment, in a live file
            occursin(fn * "(", src) || push!(dead, fn)
        end

        if !isempty(dead)
            @error """
            DEAD *_init_cold! PORT(S) — defined in src/ but never called from the live
            init path. Fortran runs InitCold inside Init; if the Julia port skips it, the
            InitAllocate NaN-fill survives into the run. Wire it into
            src/infrastructure/init_cold.jl (init_cold_biogeophys! or
            init_cold_biogeochem!), NOT into a unit test.
            """ dead
        end
        @test isempty(dead)

        # The documented exceptions must still EXIST (a rename should surface here
        # rather than silently emptying the allowlist).
        for fn in keys(EXPECTED_DEAD)
            @test fn in defined
        end
    end

    @testset "photosynthesis_timestep_init! is called from the driver" begin
        # Same dead-port class, per-STEP rather than per-INIT: Fortran calls
        # photosyns_inst%TimeStepInit at the top of CanopyFluxes
        # (CanopyFluxesMod.F90:661) to zero psnsun/psnsha/fpsn on every NON-LAKE
        # patch. The Julia port had the routine and never called it, so a
        # bare-ground patch's psn/fpsn sat at NaN for the whole run — and got
        # papered over by an isfinite() fallback in history_fpsn_patch.
        drv = read(joinpath(SRC_DIR, "driver", "clm_driver.jl"), String)
        drv_code = join(filter(l -> !startswith(lstrip(l), "#"), split(drv, '\n')), '\n')
        @test occursin("photosynthesis_timestep_init!(", drv_code)
    end

    @testset "no *_init_cold! is called ONLY from its own unit test" begin
        # Belt-and-braces on the above: a routine whose only callers live in
        # test/ or scripts/ is dead in production, however green its unit test.
        defined = defined_init_colds()
        src = live_init_source()
        for fn in sort(collect(defined))
            haskey(EXPECTED_DEAD, fn) && continue
            @test occursin(fn * "(", src)
        end
    end
end

# ---------------------------------------------------------------------------
# Runtime half: the fields the wired InitCold routines OWN must be finite on
# every ACTIVE column / patch after clm_initialize! + one driver step.
#
# GATED on the machine-local Bow inputs; absent -> skip (as the other
# domain-driven tests here do). The contrasting-domain sweeps (glacier / lake /
# urban / CN) live in scripts/initcold_nan_probe.jl, which prints the full
# per-domain table.
# ---------------------------------------------------------------------------
@testset "InitCold runtime: no NaN on active columns/patches (gated)" begin
    probe = joinpath(@__DIR__, "..", "scripts", "initcold_nan_probe.jl")
    if !isfile(probe)
        @info "InitCold runtime: probe harness absent, skipping"
        @test_skip isfile(probe)
    else
        mod = Module(:INITCOLD_PROBE)
        Core.eval(mod, :(using Test, CLM, NCDatasets, Dates, Printf))
        Base.include(mod, probe)   # PROGRAM_FILE guard: no auto-run
        BOW_FS = getfield(mod, :BOW_FS)
        BOW_FP = getfield(mod, :BOW_FP)

        if !(isfile(BOW_FS) && isfile(BOW_FP))
            @info "InitCold runtime: Bow inputs absent, skipping" BOW_FS BOW_FP
            @test_skip isfile(BOW_FS)
        else
            r = Base.invokelatest(getfield(mod, :run_domain), "BOW";
                                  fsurdat = BOW_FS, nsteps = 1)
            @test r !== nothing
            nan_fields = Set(String(x.field) for x in r.post_step)

            # Fields OWNED by the routines wired in init_cold.jl. Every one of
            # these was NaN on an active column/patch before the wiring; each is
            # read by physics, by the water/energy balance, or by a live output.
            # (`qflx_evap_tot_col` is the #209 find; `mss_*` the #210 find.)
            must_be_finite = [
                # EnergyFluxType::InitCold — urban anthropogenic terms are added
                # into the ground/canopy energy balance for EVERY patch.
                "eflx_wasteheat_patch", "eflx_traffic_patch", "eflx_ventilation_patch",
                "eflx_heat_from_ac_patch", "eflx_building_lun", "eflx_ventilation_lun",
                "rresis_patch",
                # CanopyStateType::InitCold — vegwp is the PHS plant water potential.
                "vegwp_patch", "leaf_biomass_patch", "stem_biomass_patch",
                "tlai_hist_patch", "tsai_hist_patch", "htop_hist_patch",
                # SoilStateType::InitCold
                "smp_l_col", "hk_l_col",
                # WaterFluxType::InitCold (PR #209's class)
                "qflx_evap_tot_col", "qflx_evap_veg_patch", "qflx_tran_veg_patch",
                "qflx_ev_soil_patch", "qflx_ev_snow_patch", "qflx_ev_h2osfc_patch",
                # WaterDiagnosticBulkType::InitCold
                "qflx_prec_intr_patch", "snow_layer_unity_col",
                # AerosolType::InitCold (PR #210's class)
                "mss_bcphi_col", "mss_bcpho_col", "mss_dst1_col", "mss_ocphi_col",
                "mss_cnc_bcphi_col", "mss_cnc_dst1_col",
                # TemperatureType::InitCold
                "dynbal_baseline_heat_col",
                # PhotosynthesisMod::InitCold — psnsun/psnsha are summed into the
                # patch->column->gridcell GPP aggregation on EVERY patch.
                "psnsun_patch", "psnsha_patch",
                # WaterBalanceType::InitCold
                "wa_reset_nonconservation_gain_col",
            ]
            still_nan = [f for f in must_be_finite if f in nan_fields]
            if !isempty(still_nan)
                @error "InitCold-owned fields still NaN on ACTIVE indices after init+1 step" still_nan
            end
            @test isempty(still_nan)

            # Regression ceiling. The residual NaN-on-active fields are ones
            # FORTRAN ALSO leaves non-finite: spval missing-value flags (the
            # rural/urban split diagnostics, wf/wf2, snw_rds_top), per-branch
            # diagnostics on patches whose branch did not run, arrays Fortran
            # would not even ALLOCATE in this config (the soil-BGC fluxes with
            # use_cn=false), and the accumulator buffers Fortran seeds in a
            # SEPARATE phase (InitAccBuffer/InitAccVars — t_ref2m_max/min, t_mo,
            # annlai), which is a distinct, still-open gap.
            #
            # Measured on Bow: 590 before the wiring, 476 after. The ceiling
            # catches a REGRESSION (a new dead field) without pinning the exact
            # Fortran-faithful residual, which legitimately moves.
            nbad = length(r.post_step)
            @info "InitCold runtime: NaN-on-active fields after init+1 step (Bow)" nbad
            @test nbad <= 480
        end
    end
end
