# ==========================================================================
# cn_conservation_audit.jl — DOES CARBON CONSERVE?
#
# The C half of the CN balance check could never run: gpp_col/er_col (and
# every other column/gridcell C flux) were dead p2c's. Now that they are
# wired, run the check LIVE over the validated windows and QUANTIFY the
# imbalance rather than merely crash on it.
#
# Thresholds are set to Inf here on purpose: the point is to MEASURE the
# non-closure per step, per term, per season — not to trip on the first one.
# The default (shipping) thresholds are cerror=1e-7 / nerror=1e-3 and the
# check IS fatal; that is exercised by the test suite.
#
#   julia +1.12 --project=. scripts/cn_conservation_audit.jl [summer|autumn|coldstart]
# ==========================================================================
include(joinpath(@__DIR__, "fortran_parity_common.jl"))

const AUD_SUM   = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_summer"
const AUD_AUT   = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_autumn"
const AUD_BASE  = 1753153
const AUD_SUM_N0, AUD_SUM_NSTEPS = 1757845, 28
const AUD_AUT_N0, AUD_AUT_NSTEPS = 1759897, 480

function run_audit(; n0::Int, nsteps::Int, dumpdir::String, label::String)
    start_date = DateTime(2002, 1, 1) + Hour(n0 - AUD_BASE)
    ffile = replace(FFORCING, "clmforc.2003.nc" => "clmforc.2002.nc")
    (inst, bounds, filt, tm) = build_bow_inst(; dtime=3600, start_date=start_date,
                                              use_cn=true, use_luna=true)
    if isempty(inst.photosyns.vcmx25_z_patch)
        inst.photosyns.vcmx25_z_patch = fill(30.0, bounds.endp, CLM.NLEVCAN)
        inst.photosyns.jmx25_z_patch  = fill(60.0, bounds.endp, CLM.NLEVCAN)
    end
    icdump = joinpath(dumpdir, "pdump_before_step_n$(n0).nc")
    inject_dump!(inst, bounds, icdump)
    let ds = NCDataset(icdump, "r")
        if haskey(ds, "vegwp")
            vw = ds["vegwp"][:, :]
            for pd in 1:size(vw, 2), seg in 1:4
                inst.canopystate.vegwp_patch[pd, seg] = Float64(vw[seg, pd])
            end
        end
        close(ds)
    end

    config = CLM.CLMDriverConfig(use_cn=true, use_aquifer_layer=false,
                                 use_hydrstress=true, use_luna=true)
    filt_ia = CLM.clump_filter_inactive_and_active
    ng, nc = bounds.endg, bounds.endc
    fr = CLM.ForcingReader(); CLM.forcing_reader_init!(fr, ffile)

    # Enable the check, but MEASURE instead of throwing.
    CLM.cn_balance_check_enabled!(true)
    bal = inst.bgc_vegetation.cn_balance_inst
    bal.cerror = Inf; bal.nerror = Inf
    bal.cwarning = Inf; bal.nwarning = Inf     # silence the per-step @warn spam

    let calday0 = CLM.get_curr_calday(tm)
        (d0, _)  = CLM.compute_orbital(calday0)
        (dm1, _) = CLM.compute_orbital(calday0 - 3600.0 / CLM.SECSPDAY)
        CLM.init_daylength!(inst.gridcell, d0, dm1, CLM.ORB_OBLIQR_DEFAULT, 1:bounds.endg)
    end

    cvf = inst.bgc_vegetation.cnveg_carbonflux_inst
    scf = inst.soilbiogeochem_carbonflux
    worst_c = 0.0; worst_c_step = 0
    worst_n = 0.0; worst_n_step = 0
    worst_cg = 0.0; worst_ng = 0.0
    sum_gpp = 0.0; sum_er = 0.0; sum_hr = 0.0

    println("\n=== $label: $nsteps steps from n$n0 ($(start_date)) ===")
    println("step |    max|errcb_col|  max|errcb_grc| |    max|errnb_col| |     gpp_col      er_col      hr_col")
    report_every = nsteps > 60 ? 48 : 4

    for i in 0:(nsteps - 1)
        cur = start_date + Hour(i)
        CLM.advance_timestep!(tm)
        calday = CLM.get_curr_calday(tm); (declin, _) = CLM.compute_orbital(calday)
        CLM.read_forcing_step!(fr, inst.atm2lnd, cur, ng, nc)
        CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
        (yr, mon, d, tod) = CLM.get_curr_date(tm)
        CLM.clm_drv!(config, inst, filt, filt_ia, bounds, true, calday, declin, declin,
            CLM.ORB_OBLIQR_DEFAULT, false, false, "", false;
            nstep=tm.nstep, is_first_step=false,
            is_beg_curr_day=CLM.is_beg_curr_day(tm), is_end_curr_day=CLM.is_end_curr_day(tm),
            is_beg_curr_year=CLM.is_beg_curr_year(tm), dtime=3600.0, mon=mon, day=d,
            photosyns=inst.photosyns)

        # Only soil columns carry a CN budget.
        soil = [c for c in 1:nc if filt.bgc_soilc[c]]
        ec = maximum(abs, bal.errcb_col[soil]); en = maximum(abs, bal.errnb_col[soil])
        ecg = maximum(abs, bal.errcb_grc[1:ng]); eng = maximum(abs, bal.errnb_grc[1:ng])
        if ec > worst_c; worst_c = ec; worst_c_step = i + 1; end
        if en > worst_n; worst_n = en; worst_n_step = i + 1; end
        worst_cg = max(worst_cg, ecg); worst_ng = max(worst_ng, eng)
        c1 = soil[1]
        sum_gpp += cvf.gpp_col[c1] * 3600; sum_er += cvf.er_col[c1] * 3600
        sum_hr  += scf.hr_col[c1] * 3600

        if (i % report_every == 0) || i == nsteps - 1
            println("n", lpad(i + 1, 3), " | ", rpad(round(ec, sigdigits=4), 15),
                    rpad(round(ecg, sigdigits=4), 16), " | ", rpad(round(en, sigdigits=4), 17),
                    " | ", rpad(round(cvf.gpp_col[c1], sigdigits=4), 12),
                    rpad(round(cvf.er_col[c1], sigdigits=4), 12),
                    round(scf.hr_col[c1], sigdigits=4))
        end
    end
    CLM.forcing_reader_close!(fr)
    CLM.cn_balance_check_enabled!(false)

    println("\n--- $label VERDICT ---")
    println("  worst |errcb_col| = ", worst_c, "   (step ", worst_c_step, ")   [cerror = 1e-7]")
    println("  worst |errcb_grc| = ", worst_cg)
    println("  worst |errnb_col| = ", worst_n, "   (step ", worst_n_step, ")   [nerror = 1e-3, nwarning = 1e-7]")
    println("  worst |errnb_grc| = ", worst_ng)
    println("  cumulative col-1 gpp = ", round(sum_gpp, sigdigits=6), " gC/m2   er = ",
            round(sum_er, sigdigits=6), "   hr = ", round(sum_hr, sigdigits=6))
    println("  CARBON ", worst_c <= 1e-7 ? "CONSERVES (within cerror)" : "*** DOES NOT CONSERVE ***")
    println("  NITROGEN ", worst_n <= 1e-7 ? "conserves (within nwarning 1e-7)" :
            (worst_n <= 1e-3 ? "passes nerror but EXCEEDS the 1e-7 warning" : "*** DOES NOT CONSERVE ***"))
    return worst_c, worst_n
end

mode = isempty(ARGS) ? "summer" : lowercase(ARGS[1])
if mode == "autumn"
    run_audit(; n0=AUD_AUT_N0, nsteps=AUD_AUT_NSTEPS, dumpdir=AUD_AUT, label="AUTUMN (leaf-offset)")
elseif mode == "both"
    run_audit(; n0=AUD_SUM_N0, nsteps=AUD_SUM_NSTEPS, dumpdir=AUD_SUM, label="SUMMER")
    run_audit(; n0=AUD_AUT_N0, nsteps=AUD_AUT_NSTEPS, dumpdir=AUD_AUT, label="AUTUMN (leaf-offset)")
else
    run_audit(; n0=AUD_SUM_N0, nsteps=AUD_SUM_NSTEPS, dumpdir=AUD_SUM, label="SUMMER")
end
