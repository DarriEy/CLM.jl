# ==========================================================================
# probe_cbalance_decompose.jl — WHERE does the residual carbon go?
#
# The column C balance closes to ~6e-4 gC/m2/step, not 1e-7. Decompose
# totc_col into its constituent pools across ONE step and compare each
# pool's delta against the fluxes that are supposed to move it. Whatever
# moves without a matching flux is the leak.
#
#   julia +1.12 --project=. scripts/probe_cbalance_decompose.jl
# ==========================================================================
include(joinpath(@__DIR__, "fortran_parity_common.jl"))

const SUMDIR = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_summer"
const N0 = 1757845

function poolsnap(inst, bounds)
    ccs = inst.bgc_vegetation.cnveg_carbonstate_inst
    scs = inst.soilbiogeochem_carbonstate
    c = 1
    veg_patch = 0.0
    for p in 1:bounds.endp
        w = inst.patch.wtcol[p]
        inst.patch.active[p] || continue
        (inst.patch.column[p] == c) || continue
        veg_patch += ccs.totc_patch[p] * w
    end
    return (
        totc      = scs.totc_col[c],
        cwdc      = scs.cwdc_col[c],
        totlitc   = scs.totlitc_col[c],
        totsomc   = scs.totsomc_col[c],
        ctrunc    = scs.ctrunc_col[c],
        totc_p2c  = ccs.totc_p2c_col[c],
        veg_p2c   = veg_patch,
        cpool     = ccs.cpool_patch[2],
        xsmrpool  = ccs.xsmrpool_patch[2],
        ctrunc_p  = ccs.ctrunc_patch[2],
        leafc     = ccs.leafc_patch[2],
    )
end

function main()
    start_date = DateTime(2002, 1, 1) + Hour(N0 - 1753153)
    ffile = replace(FFORCING, "clmforc.2003.nc" => "clmforc.2002.nc")
    (inst, bounds, filt, tm) = build_bow_inst(; dtime=3600, start_date=start_date,
                                              use_cn=true, use_luna=true)
    if isempty(inst.photosyns.vcmx25_z_patch)
        inst.photosyns.vcmx25_z_patch = fill(30.0, bounds.endp, CLM.NLEVCAN)
        inst.photosyns.jmx25_z_patch  = fill(60.0, bounds.endp, CLM.NLEVCAN)
    end
    icdump = joinpath(SUMDIR, "pdump_before_step_n$(N0).nc")
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
    CLM.cn_balance_check_enabled!(true)
    bal = inst.bgc_vegetation.cn_balance_inst
    bal.cerror = Inf; bal.nerror = Inf; bal.cwarning = Inf; bal.nwarning = Inf
    let calday0 = CLM.get_curr_calday(tm)
        (d0, _)  = CLM.compute_orbital(calday0)
        (dm1, _) = CLM.compute_orbital(calday0 - 3600.0 / CLM.SECSPDAY)
        CLM.init_daylength!(inst.gridcell, d0, dm1, CLM.ORB_OBLIQR_DEFAULT, 1:bounds.endg)
    end

    # Warm up a few steps so we are not looking at injection transients.
    for i in 0:5
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
    end

    before = poolsnap(inst, bounds)
    i = 6
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
    after = poolsnap(inst, bounds)
    CLM.forcing_reader_close!(fr)

    cvf = inst.bgc_vegetation.cnveg_carbonflux_inst
    scf = inst.soilbiogeochem_carbonflux
    dt = 3600.0; c = 1

    println("\n=== totc_col POOL DECOMPOSITION over one step (col 1) ===")
    println(rpad("pool", 14), rpad("before", 20), rpad("after", 20), "delta")
    for k in (:totc, :cwdc, :totlitc, :totsomc, :ctrunc, :totc_p2c, :veg_p2c,
              :cpool, :xsmrpool, :ctrunc_p, :leafc)
        b = getfield(before, k); a = getfield(after, k)
        println(rpad(String(k), 14), rpad(round(b, sigdigits=10), 20),
                rpad(round(a, sigdigits=10), 20), round(a - b, sigdigits=6))
    end

    dtotc = after.totc - before.totc
    cin  = cvf.gpp_col[c] * dt
    cout = (cvf.er_col[c] + cvf.fire_closs_col[c] + cvf.hrv_xsmrpool_to_atm_col[c] +
            cvf.xsmrpool_to_atm_col[c] + cvf.gru_conv_cflux_col[c] +
            cvf.wood_harvestc_col[c] + cvf.gru_wood_productc_gain_col[c] +
            cvf.crop_harvestc_to_cropprodc_col[c] - scf.som_c_leached_col[c]) * dt

    println("\n=== BUDGET (gC/m2 over the step) ===")
    println("  delta totc_col        = ", round(dtotc, sigdigits=8))
    println("  inputs  (gpp)         = ", round(cin, sigdigits=8))
    println("  outputs (er+fire+...) = ", round(cout, sigdigits=8))
    println("    er                  = ", round(cvf.er_col[c]*dt, sigdigits=8))
    println("    fire_closs          = ", round(cvf.fire_closs_col[c]*dt, sigdigits=8))
    println("    -som_c_leached      = ", round(-scf.som_c_leached_col[c]*dt, sigdigits=8))
    println("  RESIDUAL (in-out-d)   = ", round(cin - cout - dtotc, sigdigits=8))
    println("  bal.errcb_col[1]      = ", round(bal.errcb_col[1], sigdigits=8))

    # Is the residual explained by the veg side or the soil side?
    dveg  = after.totc_p2c - before.totc_p2c
    dsoil = (after.cwdc - before.cwdc) + (after.totlitc - before.totlitc) +
            (after.totsomc - before.totsomc) + (after.ctrunc - before.ctrunc)
    println("\n=== SPLIT ===")
    println("  delta veg (totc_p2c)  = ", round(dveg, sigdigits=8))
    println("  delta soil pools      = ", round(dsoil, sigdigits=8))
    println("  veg + soil            = ", round(dveg + dsoil, sigdigits=8), "  (should equal delta totc_col)")
    println("  ctrunc(col) delta     = ", round(after.ctrunc - before.ctrunc, sigdigits=8))
    println("  ctrunc(patch2) delta  = ", round(after.ctrunc_p - before.ctrunc_p, sigdigits=8))

    # The total residual is identically (C leaving veg) - (C arriving in soil):
    #   veg :  dveg  = gpp - ar - L_out          => L_out = gpp - ar - dveg
    #   soil:  dsoil = L_in - hr + leach         => L_in  = dsoil + hr - leach
    #   resid = gpp - ar - hr + leach - dveg - dsoil = L_out - L_in
    ar = cvf.ar_col[c] * dt
    hr = scf.hr_col[c] * dt
    leach = scf.som_c_leached_col[c] * dt
    L_out = cin - ar - dveg
    L_in  = dsoil + hr - leach
    println("\n=== VEG -> SOIL LITTER TRANSFER (the residual is exactly L_out - L_in) ===")
    println("  ar_col*dt             = ", round(ar, sigdigits=8))
    println("  hr_col*dt             = ", round(hr, sigdigits=8))
    println("  L_out (leaves veg)    = ", round(L_out, sigdigits=8))
    println("  L_in  (enters soil)   = ", round(L_in, sigdigits=8))
    println("  L_out - L_in          = ", round(L_out - L_in, sigdigits=8), "   <-- vs errcb ", round(bal.errcb_col[1], sigdigits=8))
    println("  fractional loss       = ", round((L_out - L_in) / max(abs(L_out), 1e-30) * 100, sigdigits=5), " %")

    # --- Which side leaks? Measure the actual litter FLUX on each side. ---
    nl = CLM.varpar.nlevdecomp
    dz = CLM.dzsoi_decomp[]
    nlit = CLM.varpar.i_litr_max   # litter pools (met/cel/lig)
    icwd = CLM.varpar.i_cwd

    # SOIL SIDE: the vertically-resolved litter/CWD inputs actually handed to the
    # soil pools this step (gC/m3/s -> gC/m2 over the step).
    soil_in = 0.0
    for j in 1:nl
        for i in 1:nlit
            soil_in += (cvf.phenology_c_to_litr_c_col[c, j, i] +
                        cvf.gap_mortality_c_to_litr_c_col[c, j, i]) * dz[j] * dt
        end
        soil_in += (cvf.gap_mortality_c_to_cwdc_col[c, j] +
                    cvf.fire_mortality_c_to_cwdc_col[c, j]) * dz[j] * dt
    end

    # VEG SIDE: the patch-level litterfall fluxes, area-weighted to the column.
    veg_out = 0.0
    for p in 1:bounds.endp
        inst.patch.active[p] || continue
        (inst.patch.column[p] == c) || continue
        w = inst.patch.wtcol[p]
        veg_out += (cvf.litfall_patch[p]) * w * dt
    end

    println("\n=== WHICH SIDE LEAKS? (gC/m2 over the step) ===")
    println("  L_out implied by veg dC  = ", round(L_out, sigdigits=8))
    println("  litfall_patch p2c (flux) = ", round(veg_out, sigdigits=8))
    println("  L_in  implied by soil dC = ", round(L_in, sigdigits=8))
    println("  litter+cwd flux TO soil  = ", round(soil_in, sigdigits=8))
    println("  --> veg-side residual    = ", round(L_out - veg_out, sigdigits=8),
            "   (veg dC vs its own litterfall flux)")
    println("  --> soil-side residual   = ", round(soil_in - L_in, sigdigits=8),
            "   (soil litter flux vs soil dC)")
    println("  --> TRANSFER GAP         = ", round(veg_out - soil_in, sigdigits=8),
            "   (litterfall flux OUT vs litter flux IN)")

    println("\n=== litter-input flux breakdown (gC/m2/step) ===")
    println("  varpar.i_litr_min/max = ", CLM.varpar.i_litr_min, "/", CLM.varpar.i_litr_max,
            "   size(phen,3) = ", size(cvf.phenology_c_to_litr_c_col, 3))
    for (nm, arr3, arr2) in (("phenology", cvf.phenology_c_to_litr_c_col, nothing),
                             ("gap_mort",  cvf.gap_mortality_c_to_litr_c_col, cvf.gap_mortality_c_to_cwdc_col),
                             ("m_c_fire",  cvf.m_c_to_litr_fire_col, cvf.fire_mortality_c_to_cwdc_col),
                             ("harvest",   cvf.harvest_c_to_litr_c_col, cvf.harvest_c_to_cwdc_col),
                             ("dwt_froot", cvf.dwt_frootc_to_litr_c_col, nothing),
                             ("gru",       cvf.gru_c_to_litr_c_col, cvf.gru_c_to_cwdc_col))
        s = 0.0
        for j in 1:nl, i in 1:size(arr3, 3); s += arr3[c, j, i] * dz[j] * dt; end
        scwd = 0.0
        if arr2 !== nothing
            for j in 1:nl; scwd += arr2[c, j] * dz[j] * dt; end
        end
        println("  ", rpad(nm, 12), "litr=", rpad(round(s, sigdigits=6), 16), "cwd=", round(scwd, sigdigits=6))
    end

    # Does hr_col capture ALL the decomposition respiration? Any cascade transition
    # whose donor pool is not flagged is_litter/is_soil/is_cwd/is_microbe respires
    # from the pools but never enters hr_col -> soil loses C the budget cannot see.
    nk = size(scf.decomp_cascade_hr_col, 2)
    hr_all = sum(scf.decomp_cascade_hr_col[c, k] for k in 1:nk) * dt
    println("\n=== HR COMPLETENESS ===")
    println("  hr_col*dt                        = ", round(hr, sigdigits=8))
    println("  sum_k decomp_cascade_hr_col*dt   = ", round(hr_all, sigdigits=8))
    println("  MISSING HR                       = ", round(hr_all - hr, sigdigits=8),
            abs(hr_all - hr) > 1e-12 ? "   <<<<<< hr_col UNDERCOUNTS" : "   ok")
    dc = inst.decomp_cascade
    println("  is_litter=", collect(dc.is_litter), " is_soil=", collect(dc.is_soil),
            " is_cwd=", collect(dc.is_cwd))
    println("  cascade_donor_pool=", dc.cascade_donor_pool)

    # --- VEG-SIDE: enumerate every patch C flux that REMOVES C from the veg pools,
    #     area-weighted to the column, and compare to what the routing delivered. ---
    p2cw(f) = sum(f[p] * inst.patch.wtcol[p] for p in 1:bounds.endp
                  if inst.patch.active[p] && inst.patch.column[p] == c) * dt
    gapnames = (:m_leafc_to_litter_patch, :m_frootc_to_litter_patch,
        :m_livestemc_to_litter_patch, :m_deadstemc_to_litter_patch,
        :m_livecrootc_to_litter_patch, :m_deadcrootc_to_litter_patch,
        :m_leafc_storage_to_litter_patch, :m_frootc_storage_to_litter_patch,
        :m_livestemc_storage_to_litter_patch, :m_deadstemc_storage_to_litter_patch,
        :m_livecrootc_storage_to_litter_patch, :m_deadcrootc_storage_to_litter_patch,
        :m_gresp_storage_to_litter_patch, :m_leafc_xfer_to_litter_patch,
        :m_frootc_xfer_to_litter_patch, :m_livestemc_xfer_to_litter_patch,
        :m_deadstemc_xfer_to_litter_patch, :m_livecrootc_xfer_to_litter_patch,
        :m_deadcrootc_xfer_to_litter_patch, :m_gresp_xfer_to_litter_patch)
    gap_sum = sum(p2cw(getfield(cvf, nm)) for nm in gapnames)
    phen_sum = p2cw(cvf.leafc_to_litter_patch) + p2cw(cvf.frootc_to_litter_patch)

    println("\n=== VEG-SIDE FLUX LEDGER (gC/m2/step, weighted to col 1) ===")
    println("  gap m_*_to_litter (20 terms) = ", round(gap_sum, sigdigits=8),
            "   routed(gap litr+cwd) = ", round(0.0, sigdigits=3))
    println("  phenology leafc+frootc_to_litter = ", round(phen_sum, sigdigits=8))
    println("  livestemc_to_litter          = ", round(p2cw(cvf.livestemc_to_litter_patch), sigdigits=8))
    println("  leafc_to_litter_fun          = ", round(p2cw(cvf.leafc_to_litter_fun_patch), sigdigits=8))
    println("  repr_structurec_to_litter    = ", round(p2cw(cvf.repr_structurec_to_litter_patch), sigdigits=8))
    println("  cpool_to_resp                = ", round(p2cw(cvf.cpool_to_resp_patch), sigdigits=8))
    println("  xsmrpool_recover             = ", round(p2cw(cvf.xsmrpool_recover_patch), sigdigits=8))
    println("  soilc_change (FUN, in ar)    = ", round(p2cw(cvf.soilc_change_patch), sigdigits=8))
    println("  excess_cflux                 = ", round(p2cw(cvf.excess_cflux_patch), sigdigits=8))
    println("  npp_burnedoff                = ", round(p2cw(cvf.npp_burnedoff_patch), sigdigits=8))
    println("  litfall_patch (diagnostic)   = ", round(p2cw(cvf.litfall_patch), sigdigits=8))
    println("  ---")
    println("  gap+phen (should equal L_in) = ", round(gap_sum + phen_sum, sigdigits=8))
    println("  L_in (actually reached soil) = ", round(L_in, sigdigits=8))
    println("  L_out (actually left veg)    = ", round(L_out, sigdigits=8))
    println("  UNROUTED VEG C LOSS          = ", round(L_out - (gap_sum + phen_sum), sigdigits=8))

    # Vertical litter profiles MUST satisfy  sum_j prof(p,j)*dzsoi_decomp(j) == 1,
    # or every gC of litterfall routed through them is silently rescaled.
    sbs = inst.soilbiogeochem_state
    println("\n=== VERTICAL PROFILE INTEGRALS  sum_j prof*dz  (MUST be 1.0) ===")
    for p in 2:3
        for (nm, pr) in (("leaf_prof", sbs.leaf_prof_patch), ("froot_prof", sbs.froot_prof_patch),
                         ("croot_prof", sbs.croot_prof_patch), ("stem_prof", sbs.stem_prof_patch))
            isempty(pr) && continue
            s = sum(pr[p, j] * dz[j] for j in 1:nl)
            println("  p", p, " ", rpad(nm, 12), round(s, sigdigits=12),
                    abs(s - 1) > 1e-8 ? "   <<<<<< NOT NORMALIZED" : "   ok")
        end
    end
    println("  nlevdecomp = ", nl, "   sum(dz[1:nl]) = ", round(sum(dz[1:nl]), sigdigits=8))
    return 0
end

exit(main())
