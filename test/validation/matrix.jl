# CLM.jl — Validation matrix (declarative config list)
#
# Each entry is ONE configuration to validate. The runner (scripts/validation/
# validate.jl) maps each entry → build → run → oracle set → verdict. New coverage
# is added by appending entries here, never by changing the runner.
#
# This is the dependency-free Julia form of the `configs.toml` in
# VALIDATION_HARNESS_DESIGN.md — a Vector of NamedTuples so the pairwise covering
# array (Step 5) can be generated programmatically.
#
# Schema (see VALIDATION_HARNESS_DESIGN.md §2–3):
#   id       :: String          unique slug
#   mode     :: Symbol          :sp | :cn | :fates_sp | :fates_bgc
#   flags    :: NamedTuple      extra use_* toggles, e.g. (use_luna=true, cnfire_method=:li2021)
#   domain   :: Symbol          :bow | :aripuana | :stillwater | :krycklan | :abisko |
#                               :tagus | :massa | :baltimore | :iceland
#   backend  :: Symbol          :cpu | :metal | :cuda | :amdgpu
#   init     :: Symbol          :cold | :warm  (warm = inject Fortran spun-up restart)
#   depth    :: Symbol          :smoke | :day | :season | :year | :multiyear
#   oracles  :: Vector{Symbol}  subset of ORACLE_KINDS (see below)
#   data_dep :: Bool            needs machine-local domain data (auto-skipped in CI)
#   note     :: String
#
# Oracle tiers (VALIDATION_HARNESS_DESIGN.md §1):
#   :conservation  (T2) water/energy/C-N closure + finiteness + bounded multi-year drift
#   :parity        (T1) per-field match to a Fortran reference dump
#   :byte_identity (T3) feature-gated-off run == baseline binary (metamorphic)
#   :matrix_eq     (T3) use_matrixcn == sequential cascade
#   :restart_rt    (T3) write→read→continue == uninterrupted
#   :ad_fd         (T3) AD gradient == finite difference
#   :mpi_serial    (T3) MPI 2-rank gather == serial
#   :determinism   (T3) same config run twice == bit-identical fingerprint
#   :streamflow    (T4) KGE/NSE vs gauge obs
const ORACLE_KINDS = (:conservation, :parity, :byte_identity, :matrix_eq,
                      :determinism,
                      :restart_rt, :ad_fd, :mpi_serial, :streamflow)

const DEPTH_DEFAULT_STEPS = Dict(
    :smoke     => 6,
    :day       => 48,        # 48 × 1800 s = 24 h
    :season    => 48 * 90,
    :year      => 48 * 365,
    :multiyear => 48 * 365 * 2,
)

"""
    vcfg(id; mode=:sp, flags=(;), domain=:bow, backend=:cpu, init=:cold,
         depth=:smoke, oracles=[:conservation], data_dep=true, note="")

Construct one validation-matrix entry with defaults. Validates field domains so a
typo'd oracle/mode/depth fails loudly at matrix-load time, not mid-run.
"""
function vcfg(id::AbstractString; mode::Symbol=:sp, flags::NamedTuple=(;),
              domain::Symbol=:bow, backend::Symbol=:cpu, init::Symbol=:cold,
              depth::Symbol=:smoke, oracles::Vector{Symbol}=[:conservation],
              data_dep::Bool=true, note::AbstractString="")
    mode in (:sp, :cn, :fates_sp, :fates_bgc) || error("bad mode $mode in $id")
    init in (:cold, :warm) || error("bad init $init in $id")
    haskey(DEPTH_DEFAULT_STEPS, depth) || error("bad depth $depth in $id")
    backend in (:cpu, :metal, :cuda, :amdgpu) || error("bad backend $backend in $id")
    for o in oracles
        o in ORACLE_KINDS || error("bad oracle $o in $id")
    end
    (; id=String(id), mode, flags, domain, backend, init, depth,
       oracles, data_dep, note=String(note))
end

"""
    validation_matrix() -> Vector{NamedTuple}

The full config matrix. Step 1 seeds the SP and CN Bow baselines (T2 conservation);
later build steps append single-axis sweeps, metamorphic pairs, domain/landunit
coverage, the pairwise array, and Fortran-parity rows.
"""
function validation_matrix()
    M = NamedTuple[]

    # --- Step 1: reference-domain baselines (T2 conservation, the universal verdict) ---
    push!(M, vcfg("sp-bow-baseline"; mode=:sp, domain=:bow, depth=:smoke,
                  oracles=[:conservation],
                  note="Satellite-phenology baseline on Bow — no BGC, the reference physics path."))
    push!(M, vcfg("cn-bow-baseline"; mode=:cn, domain=:bow, depth=:smoke,
                  oracles=[:conservation],
                  note="CN/BGC baseline on Bow (use_cn, FUN+FlexibleCN as in the Bow lnd_in)."))
    # A deeper SP stability probe (day-length run) — still cheap, exercises a diurnal cycle.
    push!(M, vcfg("sp-bow-day"; mode=:sp, domain=:bow, depth=:day,
                  oracles=[:conservation],
                  note="SP Bow over a full diurnal cycle — finiteness + per-step water closure."))

    # --- Step 2: metamorphic determinism (T3) + deeper-depth stability (T2) ---
    # Determinism: a fresh build run twice must be bit-identical — catches
    # nondeterminism (uninit reads, iteration-order, RNG) that conservation misses.
    push!(M, vcfg("sp-bow-determinism"; mode=:sp, domain=:bow, depth=:smoke,
                  oracles=[:conservation, :determinism],
                  note="SP Bow run twice → bit-identical final-state fingerprint."))
    push!(M, vcfg("cn-bow-determinism"; mode=:cn, domain=:bow, depth=:smoke,
                  oracles=[:conservation, :determinism],
                  note="CN Bow run twice → bit-identical (the CN cascade is deterministic)."))
    # Deeper stability — a CN diurnal cycle exercises the BGC chain over 48 steps.
    push!(M, vcfg("cn-bow-day"; mode=:cn, domain=:bow, depth=:day,
                  oracles=[:conservation],
                  note="CN Bow over a diurnal cycle — BGC + water/energy/C-N closure."))

    # NOTE: A *full* single-axis flag sweep (use_lch4/use_ozone/use_hydrstress/…)
    # needs the instance built WITH each flag's state allocated (CH4/ozone/PHS
    # fields), i.e. the build_for/clm_initialize! generalization done in Step 4
    # (domain/landunit). It is appended there so each swept flag gets a well-posed
    # instance rather than a config toggle over un-allocated state.

    # NOTE: Steps 3–7 append here (remaining metamorphic oracles, domains/landunits,
    # pairwise, Fortran parity, multi-year + streamflow).
    return M
end
