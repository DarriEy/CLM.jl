# ==========================================================================
# gpu_validate_phs_hostloop.jl — byte-identity guard for the PHS (use_hydrstress)
# photosynthesis passes' AD host-loop fallback (mirrors the non-PHS #134 pattern).
#
# Each PHS pass kernel (_psn_phs_pass{1,2,3,4}_kernel!) was factored into a shared
# per-patch `_psn_phs_pass{N}_body!`; the corresponding `psn_phs_pass{N}_update!`
# runs that body in a plain host loop when `_PSN_CI_AD_HOSTLOOP[]` is set (the path
# Enzyme reverse takes on Julia 1.12), else launches the KA kernel. This script runs
# every pass BOTH ways on the SAME input and asserts max|Δ| == 0.0 over all written
# outputs — Step A acceptance. Also checks the full photosynthesis_hydrstress! is
# byte-identical flag-on vs flag-off.
#
#   julia --project=. scripts/gpu_validate_phs_hostloop.jl
# ==========================================================================

using CLM
using Printf

include(joinpath(@__DIR__, "gpu_validate_phs_cpu_oracle.jl"))  # setparams!/build/run_phs!

const FT = Float64

# nan-aware max abs difference between two arrays (both-NaN at a slot ⇒ 0 diff).
function maxdiff(a, b)
    m = 0.0
    @inbounds for i in eachindex(a, b)
        x = a[i]; y = b[i]
        if isnan(x) && isnan(y)
            continue
        end
        d = abs(x - y)
        d > m && (m = d)
    end
    return m
end

# max|Δ| across all Array fields of two PhotosynthesisData structs.
function maxdiff_ps(p1, p2)
    m = 0.0; worst = :none
    for f in fieldnames(typeof(p1))
        a = getfield(p1, f); b = getfield(p2, f)
        a isa AbstractArray || continue
        d = maxdiff(a, b)
        d > m && (m = d; worst = f)
    end
    return m, worst
end

# Allocate the scratch photosynthesis_hydrstress! builds internally.
function make_scratch(np, nlevcan)
    return (; jmax_z_local = zeros(FT, np, 2, nlevcan), kn = zeros(FT, np),
        psn_wc_z_sun = zeros(FT, np, nlevcan), psn_wj_z_sun = zeros(FT, np, nlevcan),
        psn_wp_z_sun = zeros(FT, np, nlevcan), psn_wc_z_sha = zeros(FT, np, nlevcan),
        psn_wj_z_sha = zeros(FT, np, nlevcan), psn_wp_z_sha = zeros(FT, np, nlevcan),
        rh_leaf_sun = zeros(FT, np), rh_leaf_sha = zeros(FT, np),
        bsun_arr = ones(FT, np), bsha_arr = ones(FT, np))
end

# Run all 4 passes on B with scratch S, returning a snapshot (deepcopy ps + scratch +
# B.inp output arrays + vegwp) AFTER each pass.
function run_passes!(B, S; hostloop::Bool)
    prev = CLM._PSN_CI_AD_HOSTLOOP[]
    CLM._PSN_CI_AD_HOSTLOOP[] = hostloop
    snaps = Any[]
    try
        ps = B.ps; i = B.inp; bp = 1:B.np
        overrides = CLM.CalibrationOverrides()
        snap() = push!(snaps, (; ps = deepcopy(ps), S = deepcopy(S),
            k_soil_root = copy(i.k_soil_root), rc = copy(i.root_conductance_out),
            sc = copy(i.soil_conductance_out), vegwp = copy(B.vegwp),
            vegwp_ln = copy(B.vegwp_ln), qflx = copy(i.qflx_tran_veg), btran = copy(i.btran)))

        CLM.psn_phs_pass1_update!(ps, B.mask, i.col_of_patch, i.ivt, i.froot_carbon,
            i.rootfr, i.dz, i.tsai, i.tlai, i.froot_leaf_pft, i.root_radius_pft,
            i.root_density_pft, i.hksat, i.hk_l, i.smp_l, i.z_col, i.k_soil_root,
            i.root_conductance_out, i.soil_conductance_out, 2.0, 0.25, B.nlevsoi, bp)
        snap()

        CLM.psn_phs_pass2_update!(ps, S.kn, S.jmax_z_local, B.mask, i.ivt, i.c3psn_pft,
            i.mbbopt_pft, i.forc_pbot, i.oair, i.slatop_pft, i.leafcn_pft, i.flnr_pft,
            i.fnitr_pft, i.dayl_factor, i.t10, i.t_veg, i.tlai_z, i.par_z_sun_in,
            i.par_z_sha_in, i.vcmaxcint_sun, i.vcmaxcint_sha, i.nrad, false, false,
            0.015, B.nlevcan, ps.stomatalcond_mtd, ps.light_inhibit, ps.leafresp_method,
            overrides, bp; use_luna = false)
        snap()

        CLM.psn_phs_pass3_update!(ps, B.mask, i.col_of_patch, i.ivt, i.nrad,
            i.medlynslope_pft, i.medlynintercept_pft, i.crop_pft,
            i.forc_pbot, i.tgcm, i.rb, i.eair, i.esat_tv, i.cair, i.oair, i.qsatl, i.qaf,
            i.laisun, i.laisha, i.htop, i.tsai, i.elai, i.esai, i.fdry, i.forc_rho,
            i.par_z_sun_in, i.par_z_sha_in, S.jmax_z_local,
            i.o3coefg_sun, i.o3coefg_sha, i.o3coefv_sun, i.o3coefv_sha,
            i.k_soil_root, i.smp_l, i.z_col,
            B.vegwp, B.vegwp_ln, i.qflx_tran_veg, S.bsun_arr, S.bsha_arr,
            S.rh_leaf_sun, S.rh_leaf_sha,
            S.psn_wc_z_sun, S.psn_wj_z_sun, S.psn_wp_z_sun,
            S.psn_wc_z_sha, S.psn_wj_z_sha, S.psn_wp_z_sha,
            B.nlevsoi, ps.modifyphoto_and_lmr_forcrop, ps.stomatalcond_mtd,
            2.0e4, overrides, (p) -> false, bp)
        snap()

        CLM.psn_phs_pass4_update!(ps, B.mask, i.nrad, i.ivt, i.crop_pft,
            i.lai_z_sun_in, i.lai_z_sha_in, i.rb, i.btran,
            S.psn_wc_z_sun, S.psn_wj_z_sun, S.psn_wp_z_sun,
            S.psn_wc_z_sha, S.psn_wj_z_sha, S.psn_wp_z_sha,
            S.bsun_arr, S.bsha_arr, ps.modifyphoto_and_lmr_forcrop, bp)
        snap()
    finally
        CLM._PSN_CI_AD_HOSTLOOP[] = prev
    end
    return snaps
end

function snap_maxdiff(s1, s2)
    m, w = maxdiff_ps(s1.ps, s2.ps)
    for f in (:k_soil_root, :rc, :sc, :vegwp, :vegwp_ln, :qflx, :btran)
        d = maxdiff(getfield(s1, f), getfield(s2, f)); d > m && (m = d; w = f)
    end
    for f in fieldnames(typeof(s1.S))
        d = maxdiff(getfield(s1.S, f), getfield(s2.S, f)); d > m && (m = d; w = f)
    end
    return m, w
end

function main()
    println("="^70)
    println("STEP A byte-identity: PHS pass host-loop (flag ON) vs KA kernel (flag OFF)")
    println("="^70)
    CLM.varpar_init!(CLM.varpar, 5, 16, 0, 5)
    nlevsoi = 5
    setparams!()

    Boff = build(FT, nlevsoi); Soff = make_scratch(Boff.np, Boff.nlevcan)
    Bon  = build(FT, nlevsoi); Son  = make_scratch(Bon.np,  Bon.nlevcan)
    soff = run_passes!(Boff, Soff; hostloop = false)
    son  = run_passes!(Bon,  Son;  hostloop = true)

    worst = 0.0; allzero = true
    for k in 1:4
        m, w = snap_maxdiff(soff[k], son[k])
        @printf("  Pass %d : max|Δ| = %.3e   (worst field: %s)\n", k, m, w)
        worst = max(worst, m)
        m == 0.0 || (allzero = false)
    end

    # ---- full photosynthesis_hydrstress! end-to-end byte-identity ----
    Bf = build(FT, nlevsoi); Bg = build(FT, nlevsoi)
    let prev = CLM._PSN_CI_AD_HOSTLOOP[]
        CLM._PSN_CI_AD_HOSTLOOP[] = false; run_phs!(Bf)
        CLM._PSN_CI_AD_HOSTLOOP[] = true;  run_phs!(Bg)
        CLM._PSN_CI_AD_HOSTLOOP[] = prev
    end
    mfull, wfull = maxdiff_ps(Bf.ps, Bg.ps)
    for (a, b) in ((Bf.vegwp, Bg.vegwp), (Bf.inp.btran, Bg.inp.btran),
                   (Bf.inp.qflx_tran_veg, Bg.inp.qflx_tran_veg),
                   (Bf.inp.k_soil_root, Bg.inp.k_soil_root))
        d = maxdiff(a, b); d > mfull && (mfull = d; wfull = :ext)
    end
    @printf("\n  Full photosynthesis_hydrstress! : max|Δ| = %.3e (worst: %s)\n", mfull, wfull)
    mfull == 0.0 || (allzero = false)

    println("\n", "="^70)
    if allzero
        println("STEP A PASSED ✓ — every PHS pass + the full PHS solve are BYTE-IDENTICAL")
        println("between the KA-kernel path (flag off) and the AD host-loop (flag on).")
        return 0
    else
        @printf("STEP A FAILED ✗ — worst per-pass max|Δ| = %.3e, full = %.3e\n", worst, mfull)
        return 1
    end
end

exit(main())
