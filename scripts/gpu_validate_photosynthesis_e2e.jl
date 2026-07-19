# ==========================================================================
# gpu_validate_photosynthesis_e2e.jl — end-to-end GPU parity for the WHOLE
# non-hydrstress photosynthesis! (Pass 1/2/4 + the Ci solve cascade).
#
# There is no full-photosynthesis! test in the suite (test_photosynthesis.jl
# only covers components), so this constructs the full call directly: a small
# Float32 instance (sized ps + all per-patch / per-canopy-layer inputs), runs
# photosynthesis! (sun phase) on the CPU, adapts ps + the inputs to the GPU, runs
# the SAME call on device, and compares the mutated ps outputs.
#
#   julia --project=scripts scripts/gpu_validate_photosynthesis_e2e.jl
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

# Populate the module-global photosynthesis params (host; wrappers resolve to
# scalars before each launch). theta_cj/theta_ip are NaN after init → set real
# values (mirrors test_photosynthesis component tests).
function setparams!()
    CLM.photo_params_init!(CLM.params_inst)   # allocates array params (NaN)
    p = CLM.params_inst
    # photo_params_init! leaves the SCALAR params NaN (read from a param file in
    # production). Set them to the CLM defaults (from read_params.jl _read_scalar
    # defaults) so the full photosynthesis! produces finite values.
    p.theta_cj = fill(0.98, CLM.MXPFT + 1)
    p.theta_ip = 0.95; p.act25 = 72.0; p.fnr = 7.16; p.cp25_yr2000 = 42.75e-6
    p.kc25_coef = 404.9e-6; p.ko25_coef = 278.4e-3; p.fnps = 0.15; p.theta_psii = 0.7
    p.vcmaxha = 65330.0; p.jmaxha = 43540.0; p.tpuha = 53100.0; p.lmrha = 46390.0
    p.kcha = 79430.0; p.koha = 36380.0; p.cpha = 37830.0
    p.vcmaxhd = 149250.0; p.jmaxhd = 152040.0; p.tpuhd = 150650.0; p.lmrhd = 150650.0
    p.lmrse = 490.0; p.tpu25ratio = 0.167; p.kp25ratio = 20160.0
    p.vcmaxse_sf = 1.0; p.jmaxse_sf = 1.0; p.tpuse_sf = 1.0; p.jmax25top_sf = 1.0
end

function build(::Type{FT}) where {FT}
    np = 2
    nlevcan = CLM.NLEVCAN          # = 1
    mp = CLM.MXPFT + 1

    ps = CLM.PhotosynthesisData{FT}()
    CLM.photosynthesis_data_init!(ps, np)
    ps.stomatalcond_mtd = CLM.STOMATALCOND_MTD_BB1987
    ps.leafresp_method  = CLM.LEAFRESP_MTD_RYAN1991
    ps.light_inhibit    = false

    inp = (
        esat_tv     = fill(FT(2300.0), np),
        eair        = fill(FT(1500.0), np),
        oair        = fill(FT(20900.0), np),
        cair        = fill(FT(40.0), np),
        rb          = fill(FT(50.0), np),
        btran       = fill(FT(0.8), np),
        dayl_factor = fill(FT(1.0), np),
        leafn       = fill(FT(1.0), np),
        forc_pbot   = fill(FT(101325.0), np),
        t_veg       = fill(FT(293.0), np),
        t10         = fill(FT(290.0), np),
        tgcm        = fill(FT(293.0), np),
        nrad        = fill(1, np),                       # Int; 1 active layer
        tlai_z      = fill(FT(1.0), np, nlevcan),
        tlai        = fill(FT(1.0), np),
        par_z_in    = fill(FT(250.0), np, nlevcan),      # day time (>0)
        lai_z_in    = fill(FT(1.0), np, nlevcan),
        vcmaxcint   = fill(FT(1.0), np),
        o3coefv     = fill(FT(1.0), np),
        o3coefg     = fill(FT(1.0), np),
        c3psn_pft        = fill(FT(1.0), mp),            # C3
        leafcn_pft       = fill(FT(25.0), mp),
        flnr_pft         = fill(FT(0.1), mp),
        fnitr_pft        = fill(FT(1.0), mp),
        slatop_pft       = fill(FT(0.01), mp),
        mbbopt_pft       = fill(FT(9.0), mp),
        medlynintercept_pft = fill(FT(100.0), mp),
        medlynslope_pft     = fill(FT(6.0), mp),
        ivt         = fill(2, np),                       # 1-based PFT index
        col_of_patch = fill(1, np),
        gridcell_of_patch = fill(1, np),                 # for fractionation!
        forc_pco2   = fill(FT(40.0), 1),                 # per-gridcell CO2 partial pressure
        laisun      = fill(FT(0.5), np),                 # for photosynthesis_total!
        laisha      = fill(FT(0.5), np),
    )
    mask = fill(true, np)                                # Vector{Bool} (device-able)
    return (; np, ps, inp, mask)
end

run_psn!(ps, inp, mask, np, phase) = CLM.photosynthesis!(ps,
    inp.esat_tv, inp.eair, inp.oair, inp.cair, inp.rb, inp.btran, inp.dayl_factor,
    inp.leafn, inp.forc_pbot, inp.t_veg, inp.t10, inp.tgcm, inp.nrad, inp.tlai_z,
    inp.tlai, inp.par_z_in, inp.lai_z_in, inp.vcmaxcint, inp.o3coefv, inp.o3coefg,
    inp.c3psn_pft, inp.leafcn_pft, inp.flnr_pft, inp.fnitr_pft, inp.slatop_pft,
    inp.mbbopt_pft, inp.medlynintercept_pft, inp.medlynslope_pft, inp.ivt,
    inp.col_of_patch, mask, 1:np, phase; use_cn = false)

# Full default-path sequence: timestep-init zero → sun + sha leaf photosynthesis →
# total combine → C13 fractionation. Exercises the 3 newly-kernelized tail fns.
function run_full!(ps, inp, mask, np)
    CLM.photosynthesis_timestep_init!(ps, mask, 1:np; use_c13 = true)
    run_psn!(ps, inp, mask, np, "sun")
    run_psn!(ps, inp, mask, np, "sha")
    CLM.photosynthesis_total!(ps, inp.laisun, inp.laisha, mask, 1:np)
    CLM.fractionation!(ps, inp.forc_pbot, inp.forc_pco2, inp.par_z_in, inp.nrad,
        inp.c3psn_pft, inp.ivt, inp.col_of_patch, inp.gridcell_of_patch, mask, 1:np, "sun")
    return nothing
end

function main(backend)
    println("=" ^ 70)
    println("END-TO-END GPU parity for photosynthesis! + tail (init/total/fractionation)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU path exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)
    sv_nlevcan = CLM.NLEVCAN

    setparams!()
    H = build(FT)            # CPU reference
    setparams!()
    B = build(FT)            # to be adapted to device (fresh initial state)

    ad(x) = CLM.Adapt.adapt(device_array_type(), x)
    ps_d  = ad(B.ps)
    inp_d = (; (k => ad(getfield(B.inp, k)) for k in keys(B.inp))...)
    mask_d = to(collect(Bool, B.mask))

    if !(ps_d.ac_patch isa device_array_type())
        println("  BLOCKED: PhotosynthesisData did not move to the device under adapt.")
        return 2
    end

    run_full!(H.ps, H.inp, H.mask, H.np)       # CPU
    run_full!(ps_d, inp_d, mask_d, B.np)       # device

    checks = [
        ("rssun",     H.ps.rssun_patch,     ps_d.rssun_patch),
        ("psnsun",    H.ps.psnsun_patch,    ps_d.psnsun_patch),
        ("psnsun_wc", H.ps.psnsun_wc_patch, ps_d.psnsun_wc_patch),
        ("psnsun_wj", H.ps.psnsun_wj_patch, ps_d.psnsun_wj_patch),
        ("lmrsun",    H.ps.lmrsun_patch,    ps_d.lmrsun_patch),
        ("psnsha",    H.ps.psnsha_patch,    ps_d.psnsha_patch),    # sha phase ran
        ("ac",        H.ps.ac_patch,        ps_d.ac_patch),
        ("aj",        H.ps.aj_patch,        ps_d.aj_patch),
        ("an",        H.ps.an_patch,        ps_d.an_patch),
        ("gs_mol",    H.ps.gs_mol_patch,    ps_d.gs_mol_patch),
        ("cisun_z",   H.ps.cisun_z_patch,   ps_d.cisun_z_patch),
        ("rssun_z",   H.ps.rssun_z_patch,   ps_d.rssun_z_patch),
        ("vcmax_z",   H.ps.vcmax_z_patch,   ps_d.vcmax_z_patch),
        ("kc",        H.ps.kc_patch,        ps_d.kc_patch),
        ("cp",        H.ps.cp_patch,        ps_d.cp_patch),
        # --- newly-kernelized tail fns ---
        ("fpsn",      H.ps.fpsn_patch,      ps_d.fpsn_patch),      # photosynthesis_total!
        ("fpsn_wc",   H.ps.fpsn_wc_patch,   ps_d.fpsn_wc_patch),
        ("fpsn_wj",   H.ps.fpsn_wj_patch,   ps_d.fpsn_wj_patch),
        ("fpsn_wp",   H.ps.fpsn_wp_patch,   ps_d.fpsn_wp_patch),
        ("alphapsnsun", H.ps.alphapsnsun_patch, ps_d.alphapsnsun_patch),  # fractionation!
    ]
    nfail = 0
    for (nm, a, b) in checks
        d = reldiff(a, b); ok = d < 1f-2
        @printf("  [%s] %-11s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end
    println()
    println(nfail == 0 ? "  WHOLE photosynthesis! MATCHES CPU ON $name ($FT) ✓" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
