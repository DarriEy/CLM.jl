# ==========================================================================
# paper_fidelity_collect.jl — collect per-timestep teacher-forced Julia↔Fortran
# parity diffs into a single CSV for the GMD headline fidelity figure.
#
# For a representative set of subsystem CONFIGURATIONS (SP hydrology/energy,
# CN biogeochem in summer + autumn, C13/C14 isotopes), this:
#   1. injects the Fortran before_step dump as the shared IC,
#   2. advances EXACTLY ONE clm_drv! step under matched forcing,
#   3. diffs the live Julia state field-by-field vs the Fortran end-of-step
#      dump (after_hydrologydrainage) via the trusted compare_inst_to_dump,
#   4. appends every (harness/config, var, nabs, nrel, npts) row to
#      $PARITY_CSV_DIR/parity_fields.csv (the env-gated hook in
#      fortran_parity_common.jl).
#
# The 16-field core registry (soil/snow T, soil water liq/ice, water table,
# canopy T, LAI, canopy water) is the integrated prognostic state the model
# carries step-to-step; for CN configs a curated set of column-resolved CN
# pools (litter/soil/CWD C+N, mineral N) is appended via `extra=` so the
# biogeochem cycle is represented too.
#
#   PARITY_CSV_DIR=/path julia +1.12 --project=. scripts/paper_fidelity_collect.jl
# ==========================================================================
include(joinpath(@__DIR__, "fortran_parity_common.jl"))
using CLM, NCDatasets, Printf, Dates

const MAXSTEPS = haskey(ENV, "MAXSTEPS") ? parse(Int, ENV["MAXSTEPS"]) : 1_000_000

# curated column-resolved CN pools (no pft-remap needed; unambiguous) --------
function cn_extra(inst)
    scs = inst.soilbiogeochem_carbonstate
    sns = inst.soilbiogeochem_nitrogenstate
    cp  = scs.decomp_cpools_vr_col      # (col, lev, pool)
    np  = sns.decomp_npools_vr_col      # (col, lev, pool)
    Tuple[
        ("litr1c_vr",   :col2d, () -> @view cp[:, :, 1]),
        ("litr2c_vr",   :col2d, () -> @view cp[:, :, 2]),
        ("litr3c_vr",   :col2d, () -> @view cp[:, :, 3]),
        ("soil1c_vr",   :col2d, () -> @view cp[:, :, 4]),
        ("soil2c_vr",   :col2d, () -> @view cp[:, :, 5]),
        ("soil3c_vr",   :col2d, () -> @view cp[:, :, 6]),
        ("cwdc_vr",     :col2d, () -> @view cp[:, :, 7]),
        ("litr1n_vr",   :col2d, () -> @view np[:, :, 1]),
        ("soil1n_vr",   :col2d, () -> @view np[:, :, 4]),
        ("cwdn_vr",     :col2d, () -> @view np[:, :, 7]),
        ("sminn_vr",    :col2d, () -> sns.sminn_vr_col),
        ("smin_no3_vr", :col2d, () -> sns.smin_no3_vr_col),
        ("smin_nh4_vr", :col2d, () -> sns.smin_nh4_vr_col),
    ]
end

# a group = (label, dumpdir, steps, runner). runner(nstep) -> (inst, edump).
const SUM_BASE = 1753153   # nstep -> DateTime(2202,1,1)+Hour(n-SUM_BASE); forcing yr 2002
const F2002    = replace(FFORCING, "clmforc.2003.nc" => "clmforc.2002.nc")

function run_sp(nstep)
    # The Bow SP reference run is PHS + LUNA (see fortran_parity_validate.jl);
    # a faithful one-step comparison requires use_hydrstress=true (which also
    # seeds vegwp from the dump and turns on use_luna).
    inst, _ = run_one_parity_step!(nstep; use_hydrstress=true)  # 2003 clock + forcing
    edump   = joinpath(DUMPDIR, "pdump_after_hydrologydrainage_n$(nstep).nc")
    (inst, edump, Tuple[])
end

function run_cn(dir)
    (nstep) -> begin
        inst, _ = run_one_parity_step!(nstep; use_cn=true, dumpdir=dir,
            use_hydrstress=true, use_luna=true,
            step_date=DateTime(2002,1,1)+Hour(nstep-SUM_BASE), forcing_file=F2002)
        edump = joinpath(dir, "pdump_after_hydrologydrainage_n$(nstep).nc")
        (inst, edump, cn_extra(inst))
    end
end

const SUMMER = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_summer"
const CLMBGC = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup"
const AUTUMN = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_autumn"

# The BGC summer window (model 2202, forcing 2002) is staged contiguously across
# two dirs: bgc_ref_summer holds n1757845..1757872, the top-level clm_bgc_spinup
# dir continues n1757873..1757917 (same run, same config). Sampled every 2 steps.
const GROUPS = [
    ("CN-summer",       SUMMER,  collect(1757845:2:1757871),  run_cn(SUMMER)),
    ("CN-summer",       CLMBGC,  collect(1757873:2:1757917),  run_cn(CLMBGC)),
    ("CN-autumn",       AUTUMN,  collect(1759900:1759911),    run_cn(AUTUMN)),
    ("SP-hydro-energy", DUMPDIR, collect(8761:8763),          run_sp),
]

function main()
    total_rows = 0; ran = Dict{String,Int}(); failed = Dict{String,Int}()
    for (label, dir, steps, runner) in GROUPS
        isdir(dir) || (println("!! missing dir $dir — skipping $label"); continue)
        nrun = 0
        for nstep in steps
            nrun >= MAXSTEPS && break
            bd = joinpath(dir, "pdump_before_step_n$(nstep).nc")
            ed = joinpath(dir, "pdump_after_hydrologydrainage_n$(nstep).nc")
            (isfile(bd) && isfile(ed)) || continue
            try
                (inst, edump, extra) = runner(nstep)
                res, gmax = compare_inst_to_dump(inst, edump; label=label,
                                                 tol=1e-9, extra=extra)
                total_rows += length(res); nrun += 1
                @printf("  [%s n%d] %d fields, gmax=%.3e\n", label, nstep, length(res), gmax)
            catch e
                failed[label] = get(failed, label, 0) + 1
                @printf("  [%s n%d] FAILED: %s\n", label, nstep,
                        sprint(showerror, e))
            end
        end
        ran[label] = nrun
    end
    println("\n", "="^60)
    println("collected rows: ", total_rows)
    for (label, _, _, _) in GROUPS
        @printf("  %-18s steps ran=%d  failed=%d\n", label,
                get(ran, label, 0), get(failed, label, 0))
    end
    return 0
end

exit(main())
