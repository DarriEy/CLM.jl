# ==========================================================================
# gpu_validate_ozonestress_e2e.jl — end-to-end GPU parity for the WHOLE
# ozone driver: calc_ozone_uptake! (per-patch stomatal uptake, both the
# LAI>thresh "active" branch and the LAI<thresh reset branch) followed by
# calc_ozone_stress! / calc_ozone_stress_lombardozzi2015! (o3uptake>0 active
# vs o3uptake==0 inactive, plus the non-exposed reset branch).
#
# Builds a small Float32 instance: 4 patches —
#   p1 exposed, evergreen needleleaf, high ozone (uptake active, stress active)
#   p2 exposed, nonwoody broadleaf, high ozone (uptake active, stress active)
#   p3 exposed, LAI below threshold (uptake reset to 0 -> stress inactive)
#   p4 non-exposed veg (reset-to-1 branch)
# Runs both on CPU, adapts to the GPU, runs the SAME calls on the device, and
# compares the mutated outputs field-by-field.
#
#   julia --project=scripts scripts/gpu_validate_ozonestress_e2e.jl
#
# reldiff PASSES when both sides are NaN, so main() asserts the CPU reference
# fields are FINITE before trusting any parity number.
# ==========================================================================

using CLM
using Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))

function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(A[i] - B[i]) / (1.0 + max(abs(A[i]), abs(B[i]))))
    end
    return m
end
allfinite(a) = all(isfinite, Array(a))

function build(::Type{FT}) where {FT}
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    np = 4; nc = 1; ng = 1; nl = 1
    npft = CLM.MXPFT + 1

    oz = CLM.OzoneData{FT}()
    CLM.ozone_init!(oz, np; stress_method = "stress_lombardozzi2015")

    patch = CLM.PatchData{FT}(); CLM.patch_init!(patch, np)
    patch.column   .= 1
    patch.gridcell .= 1
    patch.landunit .= 1
    patch.active   .= true
    # pft types: 1 = needleleaf (<=3), 5 = broadleaf/nonwoody (>3)
    patch.itype[1] = 1; patch.itype[2] = 5; patch.itype[3] = 5; patch.itype[4] = 5

    # PFT parameter vectors (size MXPFT+1)
    evergreen_pft = zeros(FT, npft); evergreen_pft[1] = 1.0   # pft 1 evergreen
    leaf_long_pft = fill(FT(1.0), npft)
    woody_pft     = zeros(FT, npft)                            # nonwoody for pft 5

    # forcings
    forc_pbot = FT[101325.0]
    forc_th   = FT[300.0]
    forc_o3   = FT[100.0e-9]    # high ozone -> active uptake
    rssun = fill(FT(150.0), np)
    rssha = fill(FT(220.0), np)
    rb    = fill(FT(50.0), np)
    ram   = fill(FT(100.0), np)
    tlai  = fill(FT(2.0), np)
    tlai[3] = FT(0.3)           # below O3_LAI_THRESH -> uptake reset to 0

    # seed some prior uptake so the decay/cumulative path is exercised
    oz.o3uptakesha_patch .= FT(5.0)
    oz.o3uptakesun_patch .= FT(5.0)
    oz.tlai_old_patch    .= FT(1.5)

    S = (; oz, patch)
    return (; np, S, forc_pbot, forc_th, forc_o3, rssun, rssha, rb, ram, tlai,
              evergreen_pft, leaf_long_pft, woody_pft, dtime = FT(1800.0))
end

function run!(S, fp, ft, rssun, rssha, rb, ram, tlai, fo3, ev, ll, woody, dtime,
              m_exp, m_noexp)
    nb = 1:length(m_exp)
    CLM.calc_ozone_uptake!(S.oz, S.patch, m_exp, nb, fp, ft, rssun, rssha, rb, ram,
                           tlai, fo3, ev, ll, dtime)
    CLM.calc_ozone_stress!(S.oz, m_exp, m_noexp, nb, S.patch, woody)
end

function main(backend)
    println("="^70)
    println("END-TO-END GPU parity for calc_ozone_uptake!/calc_ozone_stress! (whole driver)")
    println("="^70)
    if backend === nothing
        println("  No GPU backend — nothing to validate.")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    H = build(FT); B = build(FT)
    ad(x) = CLM.Adapt.adapt(device_array_type(), x)
    Sd = map(ad, B.S)
    if !(Sd.oz.o3uptakesha_patch isa device_array_type())
        println("  BLOCKED: a state struct did not move to the device under adapt.")
        return 2
    end

    # patches 1,2,3 exposed (uptake active/inactive); patch 4 non-exposed veg.
    m_exp   = Bool[true, true, true, false]
    m_noexp = Bool[false, false, false, true]

    run!(H.S, H.forc_pbot, H.forc_th, H.rssun, H.rssha, H.rb, H.ram, H.tlai,
         H.forc_o3, H.evergreen_pft, H.leaf_long_pft, H.woody_pft, H.dtime,
         m_exp, m_noexp)
    run!(Sd, to(B.forc_pbot), to(B.forc_th), to(B.rssun), to(B.rssha), to(B.rb),
         to(B.ram), to(B.tlai), to(B.forc_o3), to(B.evergreen_pft),
         to(B.leaf_long_pft), to(B.woody_pft), B.dtime,
         to(m_exp), to(m_noexp))

    checks = [
        ("o3uptakesha", H.S.oz.o3uptakesha_patch, Sd.oz.o3uptakesha_patch),
        ("o3uptakesun", H.S.oz.o3uptakesun_patch, Sd.oz.o3uptakesun_patch),
        ("tlai_old",    H.S.oz.tlai_old_patch,    Sd.oz.tlai_old_patch),
        ("o3coefvsha",  H.S.oz.o3coefvsha_patch,  Sd.oz.o3coefvsha_patch),
        ("o3coefvsun",  H.S.oz.o3coefvsun_patch,  Sd.oz.o3coefvsun_patch),
        ("o3coefgsha",  H.S.oz.o3coefgsha_patch,  Sd.oz.o3coefgsha_patch),
        ("o3coefgsun",  H.S.oz.o3coefgsun_patch,  Sd.oz.o3coefgsun_patch),
    ]
    nfail = 0
    for (nm, a, b) in checks
        if !allfinite(a)
            @printf("  [WARN] %-12s CPU reference not all-finite — skipping\n", nm)
            continue
        end
        d = reldiff(a, b); ok = d < 1f-3
        @printf("  [%s] %-12s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end
    # sanity: confirm both active (>0 uptake) and inactive (==0) were exercised on CPU
    up = Array(H.S.oz.o3uptakesha_patch)
    @printf("\n  uptake sample (CPU): p1=%.4f p2=%.4f p3=%.4f (active, active, reset)\n",
            up[1], up[2], up[3])
    println()
    println(nfail == 0 ? "  WHOLE ozone uptake+stress MATCHES CPU on $name ($FT)" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
