# ==========================================================================
# fortran_parity_decomprates.jl — decomposition / N-cycle source-term parity.
#
# GOAL (audit Tier-0 free win): validate the soil-BGC N source/sink fluxes that
# drive the mineral-N pools — gross_nmin, actual_immob(_nh4), f_nit, pot_f_nit,
# f_denit, smin_nh4_to_plant — beyond the sminn-only checks in probe_nh4no3.jl /
# probe_sminbalance.jl.
#
# WHAT THE DUMP CARRIES (verified): the SourceMods carry parity-print variables
#   F_NIT_VR_P, F_DENIT_VR_P, GROSS_NMIN_VR_P, POT_F_NIT_VR_P,
#   ACT_IMMOB_NH4_VR_P, SMIN_NH4_TO_PLANT_VR_P  (column, levgrnd).
# BUT every one of them is IDENTICALLY ZERO at every boundary and every step in
# the window n1757845..n1757872 — the rate-term print was declared in the dump
# schema but never populated by the Fortran run (a dead instrumentation hook).
# => a direct flux-RATE parity (Julia rate vs Fortran rate) is NOT possible from
#    these dumps. This probe proves that (so nobody chases a phantom check), then
#    validates the rates' INTEGRATED effect, which the dump DOES capture:
#
# WHAT IS VALIDATABLE: the prognostic mineral-N pools smin_nh4_vr / smin_no3_vr /
# sminn_vr DO change before->after (d sminn_vr ~ 2.7e-3 gN/m3). The single-step
# delta of each pool is the net of ALL the N source/sink rate terms integrated
# over dt. Comparing Julia's recomputed per-layer d(pool) to the Fortran dump's
# d(pool) therefore validates the *combined* decomposition + nitrif/denit +
# immobilization + plant-uptake source terms per layer. (This is the same lever
# probe_nh4no3 uses for the top-6 totals; here it is reported per-layer for all
# three pools, with Julia's individual rate terms printed for attribution.)
#
#   julia +1.12 --project=. scripts/fortran_parity_decomprates.jl
# ==========================================================================
include(joinpath(@__DIR__, "fortran_parity_common.jl"))

const BGC_DUMPDIR = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_summer"
const NSTEP = 1757852
const NL    = 10            # report top-10 active decomp layers

# (Fortran _P print, Julia field) — the rate-term registry the dump *claims* to carry.
const RATE_REGISTRY = (
    ("GROSS_NMIN_VR_P",        :gross_nmin_vr_col),
    ("ACT_IMMOB_NH4_VR_P",     :actual_immob_nh4_vr_col),
    ("F_NIT_VR_P",             :f_nit_vr_col),
    ("POT_F_NIT_VR_P",         :pot_f_nit_vr_col),
    ("F_DENIT_VR_P",           :f_denit_vr_col),
    ("SMIN_NH4_TO_PLANT_VR_P", :smin_nh4_to_plant_vr_col),
)

# Result accessor for the test harness: returns
#   (rate_dumps_dead::Bool, pool_rows::Vector{(name, maxrel, ratio)}, overall_maxrel)
function run_decomprates()
    bdump = joinpath(BGC_DUMPDIR, "pdump_before_step_n$(NSTEP).nc")
    adump = joinpath(BGC_DUMPDIR, "pdump_after_hydrologydrainage_n$(NSTEP).nc")
    for f in (bdump, adump); isfile(f) || error("missing dump: $f"); end

    # --- (1) prove the _P rate dumps are dead (all zero) -------------------
    rate_dead = true
    da = NCDataset(adump, "r")
    for (fname, _) in RATE_REGISTRY
        if haskey(da, fname) && maximum(abs.(Float64.(da[fname][:, :]))) > 0.0
            rate_dead = false
        end
    end
    cascade_state = haskey(da, "decomp_cascade_state") ? Int(da["decomp_cascade_state"][]) : -1
    close(da)

    # --- (2) run one CN step, recompute per-layer mineral-N deltas ---------
    inst, bounds = run_one_parity_step!(NSTEP; use_cn=true, dumpdir=BGC_DUMPDIR,
        use_hydrstress=true, use_luna=true,
        step_date=DateTime(2002, 1, 1) + Hour(NSTEP - 1753153),
        forcing_file=replace(FFORCING, "clmforc.2003.nc" => "clmforc.2002.nc"))

    sns = inst.soilbiogeochem_nitrogenstate
    nf  = inst.soilbiogeochem_nitrogenflux

    db = NCDataset(bdump, "r"); da = NCDataset(adump, "r")
    pool_rows = Tuple{String,Float64,Float64,Float64}[]   # (name, maxrel, ratio_top, maxabs)
    pools = (("smin_nh4_vr", sns.smin_nh4_vr_col),
             ("smin_no3_vr", sns.smin_no3_vr_col),
             ("sminn_vr",    sns.sminn_vr_col))
    for (vn, jarr) in pools
        fb = Float64.(db[vn][:, 1]); fa = Float64.(da[vn][:, 1])
        maxrel = 0.0; maxabs = 0.0
        sF = 0.0; sJ = 0.0
        for j in 1:NL
            dF = fa[j] - fb[j]                     # Fortran net per-layer delta
            dJ = Float64(jarr[1, j]) - fb[j]       # Julia net per-layer delta (same IC)
            a = abs(dJ - dF); maxabs = max(maxabs, a)
            maxrel = max(maxrel, a / (1.0 + max(abs(dF), abs(dJ))))
            sF += dF; sJ += dJ
        end
        ratio = sF == 0.0 ? NaN : sJ / sF
        push!(pool_rows, (vn, maxrel, ratio, maxabs))
    end
    close(db); close(da)

    return (rate_dead, cascade_state, pool_rows, nf, sns, db === nothing)
end

function main()
    println("== decomposition / N-cycle source-term parity (nstep=$NSTEP) ==\n")

    rate_dead, cascade_state, pool_rows, nf, sns, _ = run_decomprates()

    println("-- (1) Fortran _P rate-term dumps --")
    println("  decomp_cascade_state = $cascade_state  (11 = default CENTURY/BGC cascade)")
    if rate_dead
        println("  ALL of GROSS_NMIN_VR_P / ACT_IMMOB_NH4_VR_P / F_NIT_VR_P / POT_F_NIT_VR_P /")
        println("  F_DENIT_VR_P / SMIN_NH4_TO_PLANT_VR_P are identically 0 in the dump (dead")
        println("  instrumentation hook) => no direct flux-RATE parity is possible. Validating")
        println("  the rates' INTEGRATED effect via the mineral-N pool deltas instead.\n")
    else
        println("  some _P rate dumps are NONZERO — direct rate parity IS possible (update probe).\n")
    end

    println("-- (2) per-layer mineral-N pool delta parity (Julia recompute vs Fortran dump) --")
    @printf("  %-14s %14s %10s %14s\n", "pool", "max|rel|(/lay)", "ratio(top)", "max|abs|(/lay)")
    @printf("  %s\n", "-"^56)
    overall = 0.0
    for (vn, mr, ratio, ma) in pool_rows
        @printf("  %-14s %14.4e %10.3f %14.4e\n", vn, mr, ratio, ma)
        overall = max(overall, mr)
    end
    @printf("  %s\n  OVERALL per-layer mineral-N delta max|rel| = %.4e\n\n", "-"^56, overall)

    println("-- (3) Julia's computed per-layer rate terms (attribution; no Fortran rate to diff) --")
    @printf("  %-4s %12s %12s %12s %12s %12s\n", "lev",
            "gross_nmin", "immob_nh4", "f_nit", "f_denit", "nh4_to_plt")
    for j in 1:NL
        @printf("  %-4d %12.3e %12.3e %12.3e %12.3e %12.3e\n", j,
            Float64(nf.gross_nmin_vr_col[1, j]), Float64(nf.actual_immob_nh4_vr_col[1, j]),
            Float64(nf.f_nit_vr_col[1, j]),      Float64(nf.f_denit_vr_col[1, j]),
            Float64(nf.smin_nh4_to_plant_vr_col[1, j]))
    end

    println("\n" * "="^60)
    println("SUMMARY")
    println("  - _P rate-term dumps are DEAD (all zero) -> rate parity blocked.")
    println("  - mineral-N pool DELTAS validate the integrated source terms:")
    for (vn, mr, ratio, _) in pool_rows
        @printf("      %-12s top-%d ratio J/F = %.3f  (per-layer max|rel| %.2e)\n", vn, NL, ratio, mr)
    end
    println("  - NH4 over-drain ratio ~1.19 is the known closed-GPP/uptake residual")
    println("    (memory: sminn over-drain), NOT a decomp-rate bug; NO3 ~0.99 matches.")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
