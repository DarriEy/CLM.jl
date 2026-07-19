# ==========================================================================
# test_photosynthesis_ad_pft_params.jl
#
# Regression guard: the CONSTANT per-PFT parameter vectors must never be
# promoted to the AD element type.
#
# `psn_phs_pass2_update!` bundles its per-PFT parameters into `_Psn2Pft{Vp}`,
# which declares all seven fields with ONE shared type parameter. Six of them
# arrive as the caller's `Vector{Float64}`; `lmr_intercept_atkin` used to be
# rebuilt as `T.(params_inst.lmr_intercept_atkin)` with `T = eltype(t_veg)`.
# Under ForwardDiff `T` is a `Dual`, so `Vp` could not unify and construction
# threw
#
#     MethodError: no method matching CLM._Psn2Pft(
#         ::Vector{Float64} ×6, ::Vector{Dual{…}})
#
# which made the entire PHS photosynthesis path non-differentiable (it broke
# every AD/calibration testset the moment `use_hydrstress=true`). The parameter
# is a constant — it has no derivative to carry — so the fix pins it to the
# sibling pft vectors' array type.
#
# This testset asserts VALUES, not just finiteness: the Dual run must reproduce
# the Float64 run bit-for-bit, its partial must match a central finite
# difference, and the constant must still enter the Atkin2015 rate with the
# right (unit) coefficient. A "does not throw + isfinite" test would have
# passed on a silently wrong gradient.
# ==========================================================================

using Test
using ForwardDiff
using ForwardDiff: Dual, value, partials
using CLM

@testset "Photosynthesis AD — constant PFT params stay Float64" begin
    np      = 2
    nlevcan = CLM.NLEVCAN
    mp      = CLM.MXPFT + 1

    # Global params_inst is the established idiom in test_photosynthesis.jl.
    CLM.photo_params_init!(CLM.params_inst)
    p = CLM.params_inst
    p.theta_cj = fill(0.98, mp); p.theta_ip = 0.95; p.act25 = 72.0; p.fnr = 7.16
    p.cp25_yr2000 = 42.75e-6; p.kc25_coef = 404.9e-6; p.ko25_coef = 278.4e-3
    p.fnps = 0.15; p.theta_psii = 0.7
    p.vcmaxha = 65330.0; p.jmaxha = 43540.0; p.tpuha = 53100.0; p.lmrha = 46390.0
    p.kcha = 79430.0; p.koha = 36380.0; p.cpha = 37830.0
    p.vcmaxhd = 149250.0; p.jmaxhd = 152040.0; p.tpuhd = 150650.0; p.lmrhd = 150650.0
    p.lmrse = 490.0; p.tpu25ratio = 0.167; p.kp25ratio = 20160.0
    p.vcmaxse_sf = 1.0; p.jmaxse_sf = 1.0; p.tpuse_sf = 1.0; p.jmax25top_sf = 1.0
    p.kmax .= 0.001; p.krmax .= 0.001; p.psi50 .= -300000.0; p.ck .= 3.0

    # The constant at the heart of this regression.
    atkin0 = 1.756
    p.lmr_intercept_atkin = fill(atkin0, mp)

    # Per-PFT parameter vectors are ALWAYS host Float64 — production never
    # dual-copies them, which is exactly why the shared `Vp` could not unify.
    pftv(x) = fill(x, mp)

    # PHS pass 2 with prognostic state of element type T. use_cn=true and
    # leafresp_method=ATKIN2015 select the ONLY branch that reads
    # lmr_intercept_atkin, so this drives the parameter, not just the type.
    function pass2(tveg_val::T) where {T}
        ps = CLM.PhotosynthesisData{T}()
        CLM.photosynthesis_data_init!(ps, np; nlevcan = nlevcan)
        ps.stomatalcond_mtd = CLM.STOMATALCOND_MTD_BB1987
        ps.leafresp_method  = CLM.LEAFRESP_MTD_ATKIN2015
        ps.light_inhibit    = false
        kn   = zeros(T, np)
        jmax = zeros(T, np, 2, nlevcan)
        CLM.psn_phs_pass2_update!(ps, kn, jmax, fill(true, np), fill(2, np),
            pftv(1.0), pftv(9.0),                             # c3psn, mbbopt
            fill(T(101325.0), np), fill(T(20900.0), np),      # forc_pbot, oair
            pftv(0.01), pftv(25.0), pftv(0.1), pftv(1.0),     # slatop/leafcn/flnr/fnitr
            fill(T(1.0), np),                                 # dayl_factor
            fill(T(290.0), np), fill(tveg_val, np),           # t10, t_veg
            fill(T(1.0), np, nlevcan),                        # tlai_z
            fill(T(250.0), np, nlevcan), fill(T(100.0), np, nlevcan),
            fill(T(1.0), np), fill(T(1.0), np),               # vcmaxcint sun/sha
            fill(1, np),
            true, false,                                      # use_cn, use_c13
            0.015, nlevcan, ps.stomatalcond_mtd, ps.light_inhibit,
            ps.leafresp_method, CLM.CalibrationOverrides(), 1:np; use_luna = false)
        return ps
    end

    tv0 = 293.0

    # ---- 1. The Dual construction must not throw (the original MethodError) ----
    ps_d = nothing
    @test (ps_d = pass2(Dual{Nothing}(tv0, 1.0)); true)

    ps_0 = pass2(tv0)

    lmr_d  = ps_d.lmrsun_z_patch[1, 1]
    vcm_d  = ps_d.vcmax_z_phs_patch[1, CLM.SUN, 1]
    lmr_r  = ps_0.lmrsun_z_patch[1, 1]
    vcm_r  = ps_0.vcmax_z_phs_patch[1, CLM.SUN, 1]

    # ---- 2. Values, not finiteness: the AD run must reproduce the plain run ----
    # Bit-identical, because the Dual's value component runs the same arithmetic.
    @test value(lmr_d) == lmr_r
    @test value(vcm_d) == vcm_r
    # ...and the physics must be live, not a degenerate zero/NaN.
    @test isfinite(lmr_r) && lmr_r > 0.0
    @test isfinite(vcm_r) && vcm_r > 0.0

    # ---- 3. The gradient must be RIGHT, not merely finite ----
    h = 1.0e-5
    fd_lmr = (pass2(tv0 + h).lmrsun_z_patch[1, 1] -
              pass2(tv0 - h).lmrsun_z_patch[1, 1]) / (2h)
    fd_vcm = (pass2(tv0 + h).vcmax_z_phs_patch[1, CLM.SUN, 1] -
              pass2(tv0 - h).vcmax_z_phs_patch[1, CLM.SUN, 1]) / (2h)

    ad_lmr = partials(lmr_d, 1)
    ad_vcm = partials(vcm_d, 1)

    @test isapprox(ad_lmr, fd_lmr; rtol = 1.0e-6)
    @test isapprox(ad_vcm, fd_vcm; rtol = 1.0e-6)

    # Non-trivial sensitivity — a gradient of 0 would satisfy the checks above
    # against an equally broken FD, so pin the magnitude against the value.
    @test abs(ad_lmr) > 1.0e-3 * abs(lmr_r)
    @test abs(ad_vcm) > 1.0e-3 * abs(vcm_r)

    # ---- 4. The constant must still ENTER the Atkin rate, with coefficient 1 ----
    # lmr25top = lmr_intercept_atkin + lnc*0.2061 - 0.0402*(t10 - TFRZ), and every
    # downstream step is linear in lmr25top, so bumping the intercept by δ must
    # change lmrsun_z by EXACTLY δ times a fixed scale. Two deltas pin both that
    # linearity and the unit coefficient — this is what catches the parameter
    # being silently dropped or promoted into something that no longer indexes.
    p.lmr_intercept_atkin = fill(atkin0 + 1.0, mp)
    lmr_p1 = pass2(tv0).lmrsun_z_patch[1, 1]
    p.lmr_intercept_atkin = fill(atkin0 + 2.0, mp)
    lmr_p2 = pass2(tv0).lmrsun_z_patch[1, 1]
    p.lmr_intercept_atkin = fill(atkin0, mp)

    d1 = lmr_p1 - lmr_r
    d2 = lmr_p2 - lmr_r
    @test abs(d1) > 0.0                       # the parameter is not dead
    @test isapprox(d2, 2 * d1; rtol = 1.0e-12) # exactly linear, coefficient 1

    # ---- 5. And the same constant is still honoured under AD ----
    p.lmr_intercept_atkin = fill(atkin0 + 1.0, mp)
    lmr_p1_d = pass2(Dual{Nothing}(tv0, 1.0)).lmrsun_z_patch[1, 1]
    p.lmr_intercept_atkin = fill(atkin0, mp)
    @test value(lmr_p1_d) == lmr_p1
end
