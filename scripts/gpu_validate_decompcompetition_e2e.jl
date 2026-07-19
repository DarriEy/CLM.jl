# ==========================================================================
# gpu_validate_decompcompetition_e2e.jl — end-to-end GPU parity for the WHOLE
# soil_bgc_competition! BGC driver (plant/heterotroph/nitrifier/denitrifier
# competition for mineral N).
#
# Builds a small multi-column / multi-level instance mirroring
# test/test_decomp_competition.jl (make_competition_data), runs soil_bgc_competition!
# on the CPU, converts every state struct to Float32 + adapts to the GPU, runs the
# SAME call on the device, and compares the mutated outputs field-by-field. The
# scenario deliberately exercises BOTH big config branches in one run:
#   * use_nitrif_denitrif = false  (the non-nitrif pathway: sminn_tot reduction,
#     resolve, residual distribution, fpg/fpi)
#   * use_nitrif_denitrif = true   (the nitrif/denitrif pathway: the big main
#     column/vertical loop + NH4/NO3 second-pass residuals + fpg/fpi)
#
#   julia --project=scripts scripts/gpu_validate_decompcompetition_e2e.jl
#
# The BGC *_init! / Data structs allocate Float64 regardless of struct type param,
# so we build at Float64 and down-convert array fields to Float32 for the device
# snapshot (Metal has no Float64).
# ==========================================================================

using CLM
using Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))

# NaN-aware mixed abs/rel diff; asserts the CPU reference is finite so a both-NaN
# false PASS can't slip through.
function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(A[i] - B[i]) / (1.0 + max(abs(A[i]), abs(B[i]))))
    end
    return m
end
cpu_has_finite(a) = any(isfinite, Array(a))

# Float32-down-converting Metal adaptor. The BGC state structs are concretely
# Float64-typed, so we adapt with a custom storage rule that down-converts float
# arrays to Float32 as it reconstructs the struct (integer/bool arrays move as-is).
struct MetalF32 end
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{<:AbstractFloat}) = device_array_type()(Float32.(x))
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{<:Integer})       = device_array_type()(x)
CLM.Adapt.adapt_storage(::MetalF32, x::AbstractArray{Bool})            = device_array_type()(x)

const NC = 3
const NLEVDECOMP = 2
const NTRANS = 2

# Build the CPU reference state (Float64) — mirrors make_competition_data from the
# oracle test, with a 2-column masked-out column to exercise the in-kernel mask.
function build(; sminn_vr_val, smin_nh4_vr_val, smin_no3_vr_val,
                 plant_ndemand_val, potential_immob_vr_val,
                 pot_f_nit_vr_val, pot_f_denit_vr_val,
                 n2_n2o_ratio_denit_vr_val=1.0, nfixation_prof_val=0.5)
    nc = NC; nlevdecomp = NLEVDECOMP; ntrans = NTRANS
    dzsoi_decomp = fill(1.0, nlevdecomp)

    st = CLM.SoilBiogeochemStateData()
    st.fpg_col            = zeros(nc)
    st.fpi_col            = zeros(nc)
    st.fpi_vr_col         = zeros(nc, nlevdecomp)
    st.nfixation_prof_col = fill(nfixation_prof_val, nc, nlevdecomp)
    st.plant_ndemand_col  = fill(plant_ndemand_val, nc)

    ns = CLM.SoilBiogeochemNitrogenStateData()
    ns.sminn_vr_col    = fill(sminn_vr_val, nc, nlevdecomp)
    ns.smin_nh4_vr_col = fill(smin_nh4_vr_val, nc, nlevdecomp)
    ns.smin_no3_vr_col = fill(smin_no3_vr_val, nc, nlevdecomp)

    nf = CLM.SoilBiogeochemNitrogenFluxData()
    nf.potential_immob_vr_col        = fill(potential_immob_vr_val, nc, nlevdecomp)
    nf.actual_immob_vr_col           = zeros(nc, nlevdecomp)
    nf.sminn_to_plant_vr_col         = zeros(nc, nlevdecomp)
    nf.sminn_to_plant_col            = zeros(nc)
    nf.actual_immob_col              = zeros(nc)
    nf.potential_immob_col           = zeros(nc)
    nf.supplement_to_sminn_vr_col    = zeros(nc, nlevdecomp)
    nf.sminn_to_denit_excess_vr_col  = zeros(nc, nlevdecomp)
    nf.actual_immob_no3_vr_col       = zeros(nc, nlevdecomp)
    nf.actual_immob_nh4_vr_col       = zeros(nc, nlevdecomp)
    nf.smin_no3_to_plant_vr_col      = zeros(nc, nlevdecomp)
    nf.smin_nh4_to_plant_vr_col      = zeros(nc, nlevdecomp)
    nf.pot_f_nit_vr_col              = fill(pot_f_nit_vr_val, nc, nlevdecomp)
    nf.pot_f_denit_vr_col            = fill(pot_f_denit_vr_val, nc, nlevdecomp)
    nf.f_nit_vr_col                  = zeros(nc, nlevdecomp)
    nf.f_denit_vr_col                = zeros(nc, nlevdecomp)
    nf.n2_n2o_ratio_denit_vr_col     = fill(n2_n2o_ratio_denit_vr_val, nc, nlevdecomp)
    nf.f_n2o_denit_vr_col            = zeros(nc, nlevdecomp)
    nf.f_n2o_nit_vr_col              = zeros(nc, nlevdecomp)
    nf.sminn_to_plant_fun_vr_col     = zeros(nc, nlevdecomp)
    nf.sminn_to_plant_fun_no3_vr_col = zeros(nc, nlevdecomp)
    nf.sminn_to_plant_fun_nh4_vr_col = zeros(nc, nlevdecomp)

    cf = CLM.SoilBiogeochemCarbonFluxData()
    cf.c_overflow_vr = zeros(nc, nlevdecomp, ntrans)

    mask_bgc_soilc = trues(nc); mask_bgc_soilc[2] = false   # exercise the mask
    pmnf_decomp_cascade   = zeros(nc, nlevdecomp, ntrans)
    p_decomp_cn_gain      = zeros(nc, nlevdecomp, ntrans)
    cascade_receiver_pool = ones(Int, ntrans)

    return (; st, nf, cf, ns, mask_bgc_soilc, dzsoi_decomp,
              pmnf_decomp_cascade, p_decomp_cn_gain, cascade_receiver_pool)
end

run_comp!(S, mask, dzsoi, pmnf, pcng, crp, state, params, use_nd) =
    CLM.soil_bgc_competition!(S.st, S.nf, S.cf, S.ns, state, params;
        mask_bgc_soilc=mask, bounds=1:NC,
        nlevdecomp=NLEVDECOMP, ndecomp_cascade_transitions=NTRANS,
        dzsoi_decomp=dzsoi, pmnf_decomp_cascade=pmnf,
        p_decomp_cn_gain=pcng, cascade_receiver_pool=crp,
        use_nitrif_denitrif=use_nd)

# Compare one scenario (CPU vs device). Returns number of failed fields.
function compare_scenario(name, kwargs, use_nd, params_kw, to)
    H = build(; kwargs...)   # CPU reference (Float64)
    B = build(; kwargs...)   # source for the device snapshot

    ad(x) = CLM.Adapt.adapt(MetalF32(), x)
    st_d = ad(B.st); nf_d = ad(B.nf); cf_d = ad(B.cf); ns_d = ad(B.ns)
    Sd = (; st=st_d, nf=nf_d, cf=cf_d, ns=ns_d)

    if !(st_d.fpg_col isa device_array_type())
        println("  BLOCKED: a state struct did not move to the device under adapt.")
        return -1
    end

    # masks: BitVector -> device Vector{Bool}
    dmask(m) = to(collect(Bool, m))

    state_cpu = CLM.SoilBGCCompetitionState(dt=1800.0, bdnr=0.01)
    state_dev = CLM.SoilBGCCompetitionState(dt=1800.0, bdnr=0.01)
    params    = CLM.SoilBGCCompetitionParams(; params_kw...)

    # CPU
    run_comp!(H, H.mask_bgc_soilc, H.dzsoi_decomp, H.pmnf_decomp_cascade,
              H.p_decomp_cn_gain, H.cascade_receiver_pool, state_cpu, params, use_nd)
    # Device — float arrays down to Float32, masks/integers move as-is
    run_comp!(Sd, dmask(B.mask_bgc_soilc), to(Float32.(B.dzsoi_decomp)),
              to(Float32.(B.pmnf_decomp_cascade)), to(Float32.(B.p_decomp_cn_gain)),
              to(B.cascade_receiver_pool), state_dev, params, use_nd)

    checks = [
        ("fpg",                  H.st.fpg_col,                       Sd.st.fpg_col),
        ("fpi",                  H.st.fpi_col,                       Sd.st.fpi_col),
        ("fpi_vr",               H.st.fpi_vr_col,                    Sd.st.fpi_vr_col),
        ("actual_immob_vr",      H.nf.actual_immob_vr_col,           Sd.nf.actual_immob_vr_col),
        ("actual_immob",         H.nf.actual_immob_col,              Sd.nf.actual_immob_col),
        ("potential_immob",      H.nf.potential_immob_col,           Sd.nf.potential_immob_col),
        ("sminn_to_plant_vr",    H.nf.sminn_to_plant_vr_col,         Sd.nf.sminn_to_plant_vr_col),
        ("sminn_to_plant",       H.nf.sminn_to_plant_col,            Sd.nf.sminn_to_plant_col),
        ("supplement_to_sminn",  H.nf.supplement_to_sminn_vr_col,    Sd.nf.supplement_to_sminn_vr_col),
        ("denit_excess_vr",      H.nf.sminn_to_denit_excess_vr_col,  Sd.nf.sminn_to_denit_excess_vr_col),
        ("actual_immob_no3_vr",  H.nf.actual_immob_no3_vr_col,       Sd.nf.actual_immob_no3_vr_col),
        ("actual_immob_nh4_vr",  H.nf.actual_immob_nh4_vr_col,       Sd.nf.actual_immob_nh4_vr_col),
        ("smin_no3_to_plant_vr", H.nf.smin_no3_to_plant_vr_col,      Sd.nf.smin_no3_to_plant_vr_col),
        ("smin_nh4_to_plant_vr", H.nf.smin_nh4_to_plant_vr_col,      Sd.nf.smin_nh4_to_plant_vr_col),
        ("f_nit_vr",             H.nf.f_nit_vr_col,                  Sd.nf.f_nit_vr_col),
        ("f_denit_vr",           H.nf.f_denit_vr_col,                Sd.nf.f_denit_vr_col),
        ("f_n2o_nit_vr",         H.nf.f_n2o_nit_vr_col,              Sd.nf.f_n2o_nit_vr_col),
        ("f_n2o_denit_vr",       H.nf.f_n2o_denit_vr_col,            Sd.nf.f_n2o_denit_vr_col),
        ("c_overflow_vr",        H.cf.c_overflow_vr,                 Sd.cf.c_overflow_vr),
    ]
    nfail = 0
    println("  --- scenario: $name  (use_nitrif_denitrif=$use_nd) ---")
    for (nm, a, b) in checks
        if !cpu_has_finite(a)
            @printf("    [WARN] %-22s CPU reference is all-NaN/Inf — skipping\n", nm)
            continue
        end
        d = reldiff(a, b); ok = d < 1f-3
        @printf("    [%s] %-22s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end
    return nfail
end

function main(backend)
    println("=" ^ 70)
    println("END-TO-END GPU parity for soil_bgc_competition! (BGC N competition)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU path exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    nfail = 0

    # Non-nitrif: N-limited competition (exercises resolve, residual, fpg/fpi)
    n = compare_scenario("non-nitrif N-limited",
        (; sminn_vr_val=0.5, smin_nh4_vr_val=0.3, smin_no3_vr_val=0.2,
           plant_ndemand_val=2.0, potential_immob_vr_val=1.0,
           pot_f_nit_vr_val=0.0, pot_f_denit_vr_val=0.0),
        false, (;), to)
    nfail += (n < 0 ? 1 : n)

    # Non-nitrif: excess denitrification path (high mineral N, low demand)
    n = compare_scenario("non-nitrif excess-denit",
        (; sminn_vr_val=1000.0, smin_nh4_vr_val=5.0, smin_no3_vr_val=5.0,
           plant_ndemand_val=0.01, potential_immob_vr_val=0.01,
           pot_f_nit_vr_val=0.0, pot_f_denit_vr_val=0.0),
        false, (; bdnr=0.5), to)
    nfail += (n < 0 ? 1 : n)

    # Nitrif/denitrif: NH4 limited, NO3 picks up the slack (the big main loop)
    n = compare_scenario("nitrif NH4-limited",
        (; sminn_vr_val=0.001, smin_nh4_vr_val=0.0001, smin_no3_vr_val=50.0,
           plant_ndemand_val=5.0, potential_immob_vr_val=2.0,
           pot_f_nit_vr_val=0.5, pot_f_denit_vr_val=0.1),
        true, (;), to)
    nfail += (n < 0 ? 1 : n)

    # Nitrif/denitrif: N not limiting (the not-limited branches + N2O fluxes)
    n = compare_scenario("nitrif N-not-limiting",
        (; sminn_vr_val=10000.0, smin_nh4_vr_val=5000.0, smin_no3_vr_val=5000.0,
           plant_ndemand_val=0.1, potential_immob_vr_val=0.05,
           pot_f_nit_vr_val=0.01, pot_f_denit_vr_val=0.005),
        true, (;), to)
    nfail += (n < 0 ? 1 : n)

    println()
    println(nfail == 0 ? "  WHOLE soil_bgc_competition! MATCHES CPU ON $name ($FT) ✓" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
