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

# CI cadence tiers (DESIGN §5): which schedule a config runs on.
#   :pr      — fast, every PR (baselines, single-axis smokes, determinism, parity-anchor)
#   :nightly — the full matrix at smoke/day depth (pairwise, bundles, domains)
#   :weekly  — deep tier (multi-year stability, streamflow KGE) — minutes-to-hours
const TIERS = (:pr, :nightly, :weekly)

"""
    vcfg(id; mode=:sp, flags=(;), domain=:bow, backend=:cpu, init=:cold,
         depth=:smoke, oracles=[:conservation], data_dep=true, tier=:nightly, note="")

Construct one validation-matrix entry with defaults. Validates field domains so a
typo'd oracle/mode/depth/tier fails loudly at matrix-load time, not mid-run.
"""
function vcfg(id::AbstractString; mode::Symbol=:sp, flags::NamedTuple=(;),
              domain::Symbol=:bow, backend::Symbol=:cpu, init::Symbol=:cold,
              depth::Symbol=:smoke, oracles::Vector{Symbol}=[:conservation],
              data_dep::Bool=true, tier::Symbol=:nightly, note::AbstractString="")
    mode in (:sp, :cn, :fates_sp, :fates_bgc) || error("bad mode $mode in $id")
    init in (:cold, :warm) || error("bad init $init in $id")
    haskey(DEPTH_DEFAULT_STEPS, depth) || error("bad depth $depth in $id")
    backend in (:cpu, :metal, :cuda, :amdgpu) || error("bad backend $backend in $id")
    tier in TIERS || error("bad tier $tier in $id")
    for o in oracles
        o in ORACLE_KINDS || error("bad oracle $o in $id")
    end
    (; id=String(id), mode, flags, domain, backend, init, depth,
       oracles, data_dep, tier, note=String(note))
end

"""
    pairwise_binary(flags::Vector{Symbol}) -> Vector{NamedTuple}

Deterministic greedy (AETG-style) 2-wise covering array over a set of BINARY flags:
the returned rows guarantee that for every pair of flags, all four value-combinations
(off/off, off/on, on/off, on/on) co-occur in at least one row. Each row is returned as
a NamedTuple of only its TRUE flags (the all-off row is the baseline and is dropped),
so a flag absent from a row means "off". Covers all 2-way flag interactions in O(log)
rows instead of 2^N. Pure/deterministic so the CI schema test can verify it.
"""
function pairwise_binary(flags::Vector{Symbol})
    n = length(flags)
    n >= 2 || return NamedTuple[]
    # the set of pairs still needing coverage: (i, j, vi, vj), i<j, vi/vj ∈ {false,true}
    UC = Set{Tuple{Int,Int,Bool,Bool}}()
    for i in 1:n-1, j in i+1:n, vi in (false, true), vj in (false, true)
        push!(UC, (i, j, vi, vj))
    end
    rows = Vector{Bool}[]
    cap = 4 * binomial(n, 2) + 1          # each row removes ≥1 pair → this many suffices
    while !isempty(UC)
        length(rows) < cap || error("pairwise_binary failed to converge (bug) for $flags")
        # ANCHOR each row to a still-uncovered pair (deterministic scan) so this
        # iteration is guaranteed to remove it — the loop cannot stall.
        anchor = nothing
        for i in 1:n-1, j in i+1:n, vi in (false, true), vj in (false, true)
            if (i, j, vi, vj) in UC; anchor = (i, j, vi, vj); break; end
        end
        (ai, aj, avi, avj) = anchor
        row = fill(false, n); assigned = fill(false, n)
        row[ai] = avi; assigned[ai] = true
        row[aj] = avj; assigned[aj] = true
        # greedily set the remaining params to cover the most additional pairs
        for k in 1:n
            assigned[k] && continue
            bestv = false; bestg = -1
            for v in (false, true)
                g = 0
                for m in 1:n
                    (assigned[m] && m != k) || continue
                    a, b = minmax(k, m); va, vb = k < m ? (v, row[m]) : (row[m], v)
                    (a, b, va, vb) in UC && (g += 1)
                end
                g > bestg && (bestg = g; bestv = v)
            end
            row[k] = bestv; assigned[k] = true
        end
        for i in 1:n-1, j in i+1:n
            delete!(UC, (i, j, row[i], row[j]))
        end
        push!(rows, copy(row))
    end
    # emit each row as a NamedTuple of its TRUE flags; drop the all-off (baseline) row
    out = NamedTuple[]
    for row in rows
        ks = Tuple(flags[k] for k in 1:n if row[k])
        isempty(ks) && continue
        push!(out, NamedTuple{ks}(ntuple(_ -> true, length(ks))))
    end
    out
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
    push!(M, vcfg("sp-bow-baseline"; mode=:sp, domain=:bow, depth=:smoke, tier=:pr,
                  oracles=[:conservation],
                  note="Satellite-phenology baseline on Bow — no BGC, the reference physics path."))
    push!(M, vcfg("cn-bow-baseline"; mode=:cn, domain=:bow, depth=:smoke, tier=:pr,
                  oracles=[:conservation],
                  note="CN/BGC baseline on Bow (use_cn, FUN+FlexibleCN as in the Bow lnd_in)."))
    # A deeper SP stability probe (day-length run) — still cheap, exercises a diurnal cycle.
    push!(M, vcfg("sp-bow-day"; mode=:sp, domain=:bow, depth=:day, tier=:pr,
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

    # --- Step 4: single-axis flag sweep (T2 conservation per flag) ---
    # The keystone payoff: build_for now allocates each flag's state (via the
    # init-gated clm_initialize! flags for cn/luna/lch4/cndv/crop/fates, or the
    # always-allocated config-only flags for hydrstress/voc/ozone/c13/c14/matrixcn/
    # irrigate). Each row flips EXACTLY ONE flag on from the relevant baseline and
    # asserts it still runs finite + balance-clean — catches "this flag alone breaks
    # something." Independent driver flags sweep from the SP baseline; CN sub-flags
    # (gated on use_cn) sweep from the CN baseline.
    #   (flag, init, note) — most flags are well-posed cold; PHS (use_hydrstress)
    #   needs a seeded vegwp (warm restart) or its Newton solve diverges to NaN at
    #   step 1 (the documented coldstart-canopy-nan gap), so it sweeps from a warm IC.
    DRIVER_SWEEP = [
        (:use_luna,          :cold, "LUNA photosynthetic N acclimation (allocates vcmx25_z/jmx25_z)"),
        (:use_hydrstress,    :warm, "plant hydraulic stress (PHS) canopy solve (vegwp seeded from restart)"),
        (:use_voc,           :cold, "MEGAN VOC emissions"),
        (:use_ozone,         :cold, "ozone stress on photosynthesis/conductance"),
        (:irrigate,          :cold, "irrigation demand/application"),
        (:use_aquifer_layer, :cold, "explicit aquifer lower boundary (zwt/recharge path)"),
    ]
    for (flag, init, note) in DRIVER_SWEEP
        push!(M, vcfg("sp-bow-$(replace(String(flag), "use_"=>"", "_"=>""))";
                      mode=:sp, domain=:bow, depth=:smoke, init=init,
                      flags=NamedTuple{(flag,)}((true,)),
                      oracles=[:conservation],
                      note="Single-axis: SP Bow + $note."))
    end
    CN_SWEEP = [
        (:use_lch4,     "methane (CH4) production/oxidation/transport state"),
        (:use_cndv,     "CNDV dynamic vegetation (allocates DGVS state)"),
        (:use_c13,      "C13 isotope tracer through the CN cascade"),
        (:use_c14,      "C14 isotope tracer + radioactive decay"),
        (:use_crop,     "prognostic crop (GDD/phenology/management)"),
        (:use_matrixcn, "matrix-CN solve == sequential cascade path"),
    ]
    for (flag, note) in CN_SWEEP
        push!(M, vcfg("cn-bow-$(replace(String(flag), "use_"=>"", "_"=>""))";
                      mode=:cn, domain=:bow, depth=:smoke,
                      flags=NamedTuple{(flag,)}((true,)),
                      oracles=[:conservation],
                      note="Single-axis: CN Bow + $note."))
    end

    # --- Step 4: domain coverage (T2 finiteness smoke on contrasting climates) ---
    # build_for's domain dispatch reuses the run_clm_streamflow.jl path. data_dep so
    # these auto-skip in CI; locally they run wherever the Symfluence inputs exist.
    for (dom, note) in [(:aripuana,   "tropical Amazon — wet/warm, no snow"),
                        (:stillwater, "semi-arid continental — dry/hot"),
                        (:krycklan,   "boreal Sweden — cold/snowy"),
                        (:abisko,     "arctic Sweden — permafrost-adjacent"),
                        (:tagus,      "Mediterranean Spain — seasonal-dry")]
        push!(M, vcfg("sp-$(dom)-smoke"; mode=:sp, domain=dom, depth=:smoke,
                      oracles=[:conservation],
                      note="Domain coverage: $note (finiteness + water closure)."))
    end

    # --- Step 3: metamorphic oracles (T3), now additive on the generalized build_for ---
    # Equivalent-method, round-trip, and invariance relations — no external truth.
    push!(M, vcfg("cn-bow-matrixeq"; mode=:cn, domain=:bow, depth=:smoke,
                  oracles=[:matrix_eq],
                  note="T3: matrix-CN solve == sequential cascade (soil C/N + veg-C fingerprint)."))
    push!(M, vcfg("sp-bow-restart"; mode=:sp, domain=:bow, depth=:smoke,
                  oracles=[:restart_rt],
                  note="T3: evolved SP prognostic state round-trips restart write→read bit-exact."))
    push!(M, vcfg("cn-bow-restart"; mode=:cn, domain=:bow, depth=:smoke,
                  oracles=[:restart_rt],
                  note="T3: evolved CN state (incl. soil C/N pools) round-trips restart bit-exact."))
    # Invariants validated by dedicated always-green machinery (DESIGN §6) — these
    # oracles record the relation + point at the authoritative check (no in-process
    # re-derivation; MPI needs separate ranks, AD-over-driver is heavy/version-sensitive).
    push!(M, vcfg("bow-invariants"; mode=:sp, domain=:bow, depth=:smoke,
                  oracles=[:mpi_serial, :ad_fd],
                  note="T3 invariants: MPI==serial (CI lane) + AD==FD (AD suite)."))

    # --- Step 5: t-wise (pairwise) covering array + realistic bundles (T2) ---
    # Pairwise covers every 2-way flag interaction in O(log) rows instead of 2^N.
    # Generated over the cold-safe independent flags (each proven solo in Step 4);
    # hydrstress is excluded (needs a warm vegwp seed) and instead appears in the
    # full-biophysics bundle below. Driver flags sweep SP; CN sub-flags sweep CN.
    for (i, fl) in enumerate(pairwise_binary([:use_luna, :use_voc, :use_ozone,
                                              :irrigate, :use_aquifer_layer]))
        push!(M, vcfg("pair-sp-$(lpad(i,2,'0'))"; mode=:sp, domain=:bow, depth=:smoke,
                      flags=fl, oracles=[:conservation],
                      note="Pairwise (SP driver flags): $(keys(fl))."))
    end
    for (i, fl) in enumerate(pairwise_binary([:use_lch4, :use_c13, :use_c14,
                                              :use_crop, :use_cndv, :use_matrixcn]))
        push!(M, vcfg("pair-cn-$(lpad(i,2,'0'))"; mode=:cn, domain=:bow, depth=:smoke,
                      flags=fl, oracles=[:conservation],
                      note="Pairwise (CN sub-flags): $(keys(fl))."))
    end

    # Realistic stacked profiles the science actually runs (DESIGN §2.3).
    push!(M, vcfg("bundle-bgc-production"; mode=:cn, domain=:bow, depth=:day,
                  flags=(use_lch4=true, use_crop=true, use_cndv=true, use_c13=true,
                         use_matrixcn=true),
                  oracles=[:conservation],
                  note="Full BGC production: CN + lch4 + crop + cndv + c13 + matrixcn."))
    push!(M, vcfg("bundle-full-biophysics"; mode=:sp, domain=:bow, depth=:smoke, init=:warm,
                  flags=(use_luna=true, use_hydrstress=true, use_voc=true, use_ozone=true),
                  oracles=[:conservation],
                  note="Full biophysics: SP + hydrstress + luna + voc + ozone (warm IC)."))
    # NOTE: a FATES-bgc + SPITFIRE bundle is deferred — FATES instance bootstrap +
    # spitfire wiring through build_for is its own effort; FATES is validated by the
    # dedicated FATES suite (Tier F: runs live through clm_drv! on 14-PFT params).

    # --- Step 6: Fortran T1 parity (the only tier with external ground truth) ---
    # Inject the Fortran before-step dump, run one clm_drv! step, diff vs the after
    # dump. The :pr-tier physics anchor. Needs DUMPDIR dumps (reports missing if absent).
    push!(M, vcfg("sp-bow-parity"; mode=:sp, domain=:bow, depth=:smoke, tier=:pr,
                  oracles=[:parity],
                  note="T1: Julia state == Fortran dump at the Bow daytime peak-sun step (n13461)."))

    # --- Step 7: deep tier — multi-year stability (T2) + streamflow realism (T4) ---
    # Multi-year drift: a long repeat-forcing run must stay finite + balance-clean every
    # step (the :conservation oracle handles any depth). Weekly tier — season depth is
    # ~thousands of steps. SP keeps it tractable; the verdict is bounded long-run drift.
    push!(M, vcfg("sp-bow-multiyear"; mode=:sp, domain=:bow, depth=:season, tier=:weekly,
                  oracles=[:conservation],
                  note="T2 deep: SP Bow over a season — long-run finiteness + water closure."))
    # Streamflow KGE/NSE vs gauge across the wired eval domains (T4 realism). Scored
    # out-of-band by run_clm_streamflow.jl; the oracle records the check + points there.
    for (dom, note) in [(:bow,       "Bow at Banff (WSC)"),
                        (:krycklan,  "Boreal Krycklan (SITES C16)"),
                        (:abisko,    "Arctic Abisko (SMHI)"),
                        (:tagus,     "Mediterranean Tagus (CEDEX)")]
        push!(M, vcfg("streamflow-$(dom)"; mode=:sp, domain=dom, depth=:smoke, tier=:weekly,
                      oracles=[:streamflow],
                      note="T4 realism: $note — KGE/NSE vs observed daily hydrograph."))
    end

    # NOTE (coverage ledger, DESIGN §6): gated-off byte-identity (T3) is covered by the
    # per-feature suite (the r1==r2 pattern, e.g. test_btran_smoothing.jl /
    # test_cnfire_wiring.jl / test_distributed_driver.jl); run-reproducibility by the
    # :determinism oracle. A FATES-bgc+SPITFIRE bundle is deferred to the FATES suite.
    return M
end
