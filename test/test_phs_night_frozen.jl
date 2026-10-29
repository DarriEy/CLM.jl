# ==========================================================================
# PHS (use_hydrstress) night path on frozen soil — regression for the
# cn_annual step-1 NaN.
#
# The cn_annual harness (scripts/fortran_parity_cn_annual.jl) starts from the
# 2202-01-01 Bow BGC restart: hour 1 is NIGHT on near-frozen soil, where
# hk_l ~ 3e-14 floors every soil conductance at 1e-16 and k_soil_root sits at
# ~3e-15 per layer (0 in layer 1, by construction). The original defect was in
# the night branch of calcstress: it computed the night soil water flux and
# DISCARDED it, so qflx_tran_veg stayed at its NaN init and the NaN cascaded into
# t_veg / t_grnd / h2osoi at step 1. Fortran sets it in the same place:
#
#   PhotosynthesisMod.F90 (ctsm5.3.012) calcstress, lines 4669-4676:
#     if (night) then
#        gs0sun=bsun*gs_mol_sun ; gs0sha=bsha*gs_mol_sha
#        call getvegwp(p, c, x, gb_mol, gs0sun, gs0sha, qsatl, qaf, soilflux, ...)
#        if (soilflux<0._r8) soilflux = 0._r8
#        qflx_tran_veg(p) = soilflux
#     endif
#
# These tests pin that `_calcstress` returns that night soil flux, finite, and
# equal to getvegwp's per-layer sum (PhotosynthesisMod.F90 5049-5052) at the
# returned root potential, in exactly the regime that produced the NaN.
# ==========================================================================

@testset "PHS night path on frozen soil (cn_annual step-1 regression)" begin
    T = Float64
    nlevsoi = 20
    np, nc = 1, 1
    p, c = 1, 1
    ivt = 1

    # Per-PFT curve params (Bow needleleaf-evergreen-like magnitudes; only one PFT).
    phs = CLM.PsnPhsParams(
        krmax    = fill(1.4e-9, 1),
        theta_cj = fill(0.98, 1),
        kmax     = fill(2.0e-8, 1, CLM.NVEGWCS),
        psi50    = fill(-2.75e5, 1, CLM.NVEGWCS),
        ck       = fill(3.95, 1, CLM.NVEGWCS),
        theta_ip = 0.95)

    # Soil: depths of the CLM5 layer nodes (m), uniform frozen-soil matric potential.
    zsoi = [0.01, 0.04, 0.09, 0.16, 0.26, 0.40, 0.58, 0.80, 1.06, 1.36,
            1.70, 2.08, 2.50, 2.99, 3.58, 4.27, 5.06, 5.95, 6.94, 8.03]
    z_col = reshape(zsoi, nc, nlevsoi)
    smp_l = fill(-3.7e5, nc, nlevsoi)

    # The cn_annual hour-1 state: veg potentials from the restart, the night
    # signal xsun > 0, and a nearly-closed (frozen) soil-root pathway.
    vegwp0 = (1.0, -2.88e5, -2.88e5, -2.80e5)
    gb_mol, gsmin = 2.0e5, 100.0          # boundary-layer and minimum stomatal conductance
    qsatl, qaf = 2.5e-3, 2.0e-3           # dry air => positive potential evaporation
    laisun, laisha, htop, tsai = 0.01, 0.02, 0.5, 0.06
    elai, esai, fdry = 0.03, 0.06, 0.9
    forc_rho, forc_pbot, tgcm = 1.3, 8.7e4, 263.0

    run(ksr) = CLM._calcstress(vegwp0..., gb_mol, gsmin, gsmin, qsatl, qaf,
        ksr, smp_l, z_col, p, c, laisun, laisha, htop, tsai, ivt, nlevsoi,
        elai, esai, fdry, forc_rho, forc_pbot, tgcm, phs, CLM.RGAS)

    # Fortran getvegwp's soil flux at a given root potential (5049-5052).
    fortran_soilflux(ksr, xroot) =
        sum(ksr[p, j] * (smp_l[c, j] - xroot - 1000.0 * z_col[c, j]) for j in 1:nlevsoi)

    for (label, ksr) in (
            # floored conductances: layer 1 is 0 (Fortran sets k=0 for j==1)
            ("floored (cn_annual hour 1)", hcat(0.0, fill(3.0e-15, 1, nlevsoi - 1))),
            # thawing: an open pathway, so the night flux is non-trivial
            ("thawing", hcat(0.0, fill(1.0e-10, 1, nlevsoi - 1))),
            # completely closed pathway: Fortran's sum(k)==0 branch
            ("all zero", zeros(np, nlevsoi)))
        @testset "$label" begin
            xsun, xsha, xxyl, xroot, bsun, bsha, soilflux = run(ksr)
            # Before the fix every one of these fed a NaN qflx_tran_veg downstream.
            @test all(isfinite, (xsun, xsha, xxyl, xroot, bsun, bsha, soilflux))
            # The night branch must hand back getvegwp's soil flux, not discard it.
            @test soilflux ≈ fortran_soilflux(ksr, xroot) atol = 1e-25
            label == "thawing" && @test soilflux > 0.0
            # Night transpiration (qflx_tran_veg = max(soilflux, 0)) is finite and
            # cannot exceed what the near-frozen roots can deliver.
            q = max(soilflux, 0.0)
            @test isfinite(q) && q <= sum(ksr) * 1.0e6
            # Stress factors are bounded plc values.
            @test 0.0 <= bsun <= 1.0 && 0.0 <= bsha <= 1.0
            # Water potentials stay physical (<= 0): the night signal was consumed.
            @test xsun <= 0.0 && xsha <= 0.0 && xxyl <= 0.0
        end
    end
end
