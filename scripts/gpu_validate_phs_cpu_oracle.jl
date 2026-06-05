# ==========================================================================
# gpu_validate_phs_cpu_oracle.jl — STANDALONE CPU oracle harness for the PHS
# (plant-hydraulic-stress) photosynthesis path: CLM.photosynthesis_hydrstress!.
#
# This function is config-gated off (use_hydrstress defaults false) and is never
# called anywhere in src/ or test/, so it had likely never executed. This harness
# constructs a small Float64 case (np=2 patches, single column, nlevsoi=5,
# nlevcan=1), calls photosynthesis_hydrstress! on the CPU, and prints the key
# mutated outputs. It is the CPU baseline / oracle for later GPU kernelization.
#
# Inputs are cribbed from:
#   * test/test_photosynthesis.jl   (soil-hydraulics arrays: k_soil_root, smp, z_c,
#                                     kmax/krmax/psi50/ck params, getvegwp!/spacF! tests)
#   * scripts/gpu_validate_photosynthesis_e2e.jl (PhotosynthesisData build pattern,
#                                                  setparams!, photosynthesis_data_init!)
#
#   julia --project=scripts scripts/gpu_validate_phs_cpu_oracle.jl
# ==========================================================================

using CLM
using Printf

# --------------------------------------------------------------------------
# Populate the module-global photosynthesis params (scalars + PHS arrays).
# Scalars mirror gpu_validate_photosynthesis_e2e.jl; the PHS vulnerability-curve
# / conductance params (kmax, krmax, psi50, ck) mirror the spacF!/getvegwp! tests
# in test_photosynthesis.jl.
# --------------------------------------------------------------------------
function setparams!()
    CLM.photo_params_init!(CLM.params_inst)   # allocates array params (NaN)
    p = CLM.params_inst
    # ---- scalar kinetics (same as the non-PHS e2e harness) ----
    p.theta_cj = fill(0.98, CLM.MXPFT + 1)
    p.theta_ip = 0.95; p.act25 = 72.0; p.fnr = 7.16; p.cp25_yr2000 = 42.75e-6
    p.kc25_coef = 404.9e-6; p.ko25_coef = 278.4e-3; p.fnps = 0.15; p.theta_psii = 0.7
    p.vcmaxha = 65330.0; p.jmaxha = 43540.0; p.tpuha = 53100.0; p.lmrha = 46390.0
    p.kcha = 79430.0; p.koha = 36380.0; p.cpha = 37830.0
    p.vcmaxhd = 149250.0; p.jmaxhd = 152040.0; p.tpuhd = 150650.0; p.lmrhd = 150650.0
    p.lmrse = 490.0; p.tpu25ratio = 0.167; p.kp25ratio = 20160.0
    p.vcmaxse_sf = 1.0; p.jmaxse_sf = 1.0; p.tpuse_sf = 1.0; p.jmax25top_sf = 1.0
    # ---- PHS plant-hydraulics params (per-PFT × per-segment) ----
    p.kmax  .= 0.001        # maximum segment conductance
    p.krmax .= 0.001        # maximum root-segment conductance
    p.psi50 .= -300000.0    # water potential at 50% loss of conductance [mm]
    p.ck    .= 3.0          # vulnerability-curve shape exponent
    return nothing
end

# --------------------------------------------------------------------------
# Build a small Float64 PHS case.
# --------------------------------------------------------------------------
function build(::Type{FT}, nlevsoi::Int) where {FT}
    np      = 2
    nlevcan = CLM.NLEVCAN          # = 1
    nc      = 1                    # single column
    mp      = CLM.MXPFT + 1

    ps = CLM.PhotosynthesisData{FT}()
    CLM.photosynthesis_data_init!(ps, np; nlevcan=nlevcan)
    ps.stomatalcond_mtd            = CLM.STOMATALCOND_MTD_BB1987
    ps.leafresp_method             = CLM.LEAFRESP_MTD_RYAN1991
    ps.light_inhibit               = false
    ps.modifyphoto_and_lmr_forcrop = false

    # vegwp / vegwp_ln must be seeded to a plausible (negative) leaf-water
    # potential — calcstress!/getvegwp! march from this initial guess.
    vegwp    = fill(FT(-50000.0), np, CLM.NVEGWCS)   # (np, NVEGWCS)
    vegwp_ln = fill(FT(0.0),      np, CLM.NVEGWCS)

    inp = (
        # ---- per-patch atmosphere / leaf forcing ----
        esat_tv     = fill(FT(2300.0),   np),
        eair        = fill(FT(1500.0),   np),
        oair        = fill(FT(20900.0),  np),
        cair        = fill(FT(40.0),     np),
        rb          = fill(FT(50.0),     np),
        btran       = fill(FT(0.8),      np),
        dayl_factor = fill(FT(1.0),      np),
        leafn       = fill(FT(1.0),      np),
        qsatl       = fill(FT(0.02),     np),
        qaf         = fill(FT(0.015),    np),
        forc_pbot   = fill(FT(101325.0), np),
        forc_rho    = fill(FT(1.2),      nc),         # column-indexed
        t_veg       = fill(FT(293.0),    np),
        t10         = fill(FT(290.0),    np),
        tgcm        = fill(FT(293.0),    np),
        nrad        = fill(1,            np),         # Int; 1 active canopy layer
        tlai_z      = fill(FT(1.0),      np, nlevcan),
        tlai        = fill(FT(1.0),      np),
        tsai        = fill(FT(0.5),      np),
        par_z_sun_in = fill(FT(250.0),   np, nlevcan), # day time (>0)
        par_z_sha_in = fill(FT(100.0),   np, nlevcan),
        lai_z_sun_in = fill(FT(1.0),     np, nlevcan),
        lai_z_sha_in = fill(FT(1.0),     np, nlevcan),
        vcmaxcint_sun = fill(FT(1.0),    np),
        vcmaxcint_sha = fill(FT(1.0),    np),
        o3coefv_sun  = fill(FT(1.0),     np),
        o3coefg_sun  = fill(FT(1.0),     np),
        o3coefv_sha  = fill(FT(1.0),     np),
        o3coefg_sha  = fill(FT(1.0),     np),
        # ---- per-PFT params (1-based PFT index) ----
        c3psn_pft        = fill(FT(1.0),  mp),        # C3
        leafcn_pft       = fill(FT(25.0), mp),
        flnr_pft         = fill(FT(0.1),  mp),
        fnitr_pft        = fill(FT(1.0),  mp),
        slatop_pft       = fill(FT(0.01), mp),
        mbbopt_pft       = fill(FT(9.0),  mp),
        medlynintercept_pft = fill(FT(100.0), mp),
        medlynslope_pft     = fill(FT(6.0),   mp),
        froot_leaf_pft   = fill(FT(1.0),  mp),
        root_radius_pft  = fill(FT(0.29e-3), mp),     # CLM default fine-root radius [m]
        root_density_pft = fill(FT(0.31e6),  mp),     # CLM default root density [g/m3]
        crop_pft         = fill(FT(0.0),  mp),        # not crop
        # ---- patch / column topology ----
        ivt          = fill(2, np),                   # 1-based PFT index
        col_of_patch = fill(1, np),                   # both patches in column 1
        # ---- soil-hydraulics stack (column-indexed, per soil layer) ----
        froot_carbon = fill(FT(200.0),  np),
        croot_carbon = fill(FT(200.0),  np),
        k_soil_root  = zeros(FT, np, nlevsoi),        # OUTPUT (filled in Pass 1)
        root_conductance_out = zeros(FT, np, nlevsoi),
        soil_conductance_out = zeros(FT, np, nlevsoi),
        rootfr       = fill(FT(1.0 / nlevsoi), np, nlevsoi),
        dz           = fill(FT(0.1),    nc, nlevsoi),
        z_col        = repeat(reshape(collect(FT, range(0.05, 1.0, length=nlevsoi)), 1, nlevsoi), nc, 1),
        hk_l         = fill(FT(1.0e-5), nc, nlevsoi),
        hksat        = fill(FT(1.0e-4), nc, nlevsoi),
        smp_l        = fill(FT(-5000.0), nc, nlevsoi),
        # ---- canopy / leaf geometry ----
        laisun       = fill(FT(2.0),    np),
        laisha       = fill(FT(3.0),    np),
        elai         = fill(FT(4.0),    np),
        esai         = fill(FT(0.5),    np),
        htop         = fill(FT(10.0),   np),
        fdry         = fill(FT(0.8),    np),
        qflx_tran_veg = zeros(FT, np),                # OUTPUT
    )
    mask = fill(true, np)
    return (; np, nc, nlevsoi, nlevcan, ps, vegwp, vegwp_ln, inp, mask)
end

run_phs!(B) = CLM.photosynthesis_hydrstress!(B.ps,
    B.inp.esat_tv, B.inp.eair, B.inp.oair, B.inp.cair, B.inp.rb, B.inp.btran,
    B.inp.dayl_factor, B.inp.leafn, B.inp.qsatl, B.inp.qaf, B.inp.forc_pbot,
    B.inp.forc_rho, B.inp.t_veg, B.inp.t10, B.inp.tgcm, B.inp.nrad, B.inp.tlai_z,
    B.inp.tlai, B.inp.tsai, B.inp.par_z_sun_in, B.inp.par_z_sha_in,
    B.inp.lai_z_sun_in, B.inp.lai_z_sha_in, B.inp.vcmaxcint_sun, B.inp.vcmaxcint_sha,
    B.inp.o3coefv_sun, B.inp.o3coefg_sun, B.inp.o3coefv_sha, B.inp.o3coefg_sha,
    B.inp.c3psn_pft, B.inp.leafcn_pft, B.inp.flnr_pft, B.inp.fnitr_pft,
    B.inp.slatop_pft, B.inp.mbbopt_pft, B.inp.medlynintercept_pft,
    B.inp.medlynslope_pft, B.inp.froot_leaf_pft, B.inp.root_radius_pft,
    B.inp.root_density_pft, B.inp.crop_pft, B.inp.ivt, B.inp.col_of_patch,
    B.mask, 1:B.np, B.inp.froot_carbon, B.inp.croot_carbon, B.inp.k_soil_root,
    B.inp.root_conductance_out, B.inp.soil_conductance_out, B.inp.rootfr,
    B.inp.dz, B.inp.z_col, B.inp.hk_l, B.inp.hksat, B.inp.smp_l,
    B.vegwp, B.vegwp_ln, B.inp.laisun, B.inp.laisha, B.inp.elai, B.inp.esai,
    B.inp.htop, B.inp.fdry, B.inp.qflx_tran_veg;
    nlevcan=B.nlevcan, nlevsoi=B.nlevsoi, use_cn=false)

function main()
    println("=" ^ 70)
    println("CPU ORACLE for photosynthesis_hydrstress! (PHS / plant-hydraulic-stress)")
    println("=" ^ 70)

    # varpar must be initialized before any *_init! that sizes from it.
    CLM.varpar_init!(CLM.varpar, 5, 16, 0, 5)   # UNSET layerstruct -> nlevsoi=20
    nlevsoi = 5                                 # but we drive PHS with a small column

    setparams!()
    FT = Float64
    B  = build(FT, nlevsoi)

    run_phs!(B)

    ps = B.ps
    println("\n--- key mutated outputs (patch 1) ---")
    @printf("  psnsun_patch[1]      = %.8g\n", ps.psnsun_patch[1])
    @printf("  psnsha_patch[1]      = %.8g\n", ps.psnsha_patch[1])
    @printf("  an_sun_patch[1,1]    = %.8g\n", ps.an_sun_patch[1, 1])
    @printf("  an_sha_patch[1,1]    = %.8g\n", ps.an_sha_patch[1, 1])
    @printf("  gs_mol_sun_patch[1,1]= %.8g\n", ps.gs_mol_sun_patch[1, 1])
    @printf("  gs_mol_sha_patch[1,1]= %.8g\n", ps.gs_mol_sha_patch[1, 1])
    @printf("  rssun_patch[1]       = %.8g\n", ps.rssun_patch[1])
    @printf("  rssha_patch[1]       = %.8g\n", ps.rssha_patch[1])
    @printf("  lmrsun_patch[1]      = %.8g\n", ps.lmrsun_patch[1])
    @printf("  vcmax_z_phs[1,SUN,1] = %.8g\n", ps.vcmax_z_phs_patch[1, CLM.SUN, 1])
    @printf("  cisun_z_patch[1,1]   = %.8g\n", ps.cisun_z_patch[1, 1])
    @printf("  vegwp[1,:]           = %s\n", string(B.vegwp[1, :]))
    @printf("  btran[1]             = %.8g\n", B.inp.btran[1])
    @printf("  qflx_tran_veg[1]     = %.8g\n", B.inp.qflx_tran_veg[1])
    @printf("  k_soil_root[1,:]     = %s\n", string(B.inp.k_soil_root[1, :]))

    # finiteness gate over all the headline outputs
    outs = Dict(
        "psnsun"     => ps.psnsun_patch,
        "psnsha"     => ps.psnsha_patch,
        "an_sun"     => ps.an_sun_patch,
        "an_sha"     => ps.an_sha_patch,
        "gs_mol_sun" => ps.gs_mol_sun_patch,
        "gs_mol_sha" => ps.gs_mol_sha_patch,
        "rssun"      => ps.rssun_patch,
        "rssha"      => ps.rssha_patch,
        "lmrsun"     => ps.lmrsun_patch,
        "vegwp"      => B.vegwp,
        "btran"      => B.inp.btran,
        "qflx_tran"  => B.inp.qflx_tran_veg,
    )
    nbad = 0
    for (nm, a) in outs
        if !all(isfinite, a)
            @printf("  NON-FINITE in %s: %s\n", nm, string(Array(a)))
            nbad += 1
        end
    end
    println()
    if nbad == 0
        println("  ALL HEADLINE OUTPUTS FINITE — PHS CPU oracle established ✓")
        return 0
    else
        println("  $nbad output group(s) contained non-finite values — investigate.")
        return 1
    end
end

exit(main())