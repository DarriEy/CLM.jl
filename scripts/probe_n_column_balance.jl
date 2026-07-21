# ==========================================================================
# probe_n_column_balance.jl — WHERE does the column N budget fail to close?
#
# Dump, for the Bow soil column at one summer step, every term the CN N balance
# check reads, plus the crux test:  p2c(sminn_to_npool_patch)  vs  sminn_to_plant_col.
# If those differ, the plant takes up more (or less) N than is debited from the
# soil mineral pool — N created/destroyed in the active-uptake path.
#
#   julia +1.12 --project=. scripts/probe_n_column_balance.jl
# ==========================================================================
include(joinpath(@__DIR__, "fortran_parity_common.jl"))

const SUMDIR = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_summer"
const N0 = 1757845

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

    function step!(i)
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
    for i in 0:6; step!(i); end
    CLM.forcing_reader_close!(fr)

    snf = inst.soilbiogeochem_nitrogenflux
    nvf = inst.bgc_vegetation.cnveg_nitrogenflux_inst
    dt = 3600.0
    soil = [c for c in 1:nc if filt.bgc_soilc[c]]
    c1 = soil[1]

    # p2c of sminn_to_npool
    p2c_smin = 0.0
    for p in bounds.begp:bounds.endp
        (filt.bgc_vegp[p]) || continue
        w = inst.patch.wtcol[p]; isfinite(w) || continue
        (inst.patch.column[p] == c1) || continue
        p2c_smin += nvf.sminn_to_npool_patch[p] * w
    end

    getcol(nm) = (arr = getfield(snf, nm); isempty(arr) ? 0.0 : Float64(arr[c1]))

    # FUN decomposition: veg uptake (sminn_to_plant_fun) vs its components, and the
    # vr soil debit summed over layers.
    dzsoi = CLM.dzsoi_decomp[]
    p2c_fun = 0.0; p2c_nfix = 0.0; p2c_npass = 0.0; p2c_nact = 0.0; p2c_nnonmyc = 0.0
    p2c_retr = 0.0; p2c_freeretr = 0.0
    for p in bounds.begp:bounds.endp
        (filt.bgc_vegp[p]) || continue
        w = inst.patch.wtcol[p]; isfinite(w) || continue
        (inst.patch.column[p] == c1) || continue
        gp(nm) = (a = getfield(nvf, nm); isempty(a) ? 0.0 : Float64(a[p]))
        p2c_fun      += gp(:sminn_to_plant_fun_patch) * w
        p2c_nfix     += gp(:Nfix_patch) * w
        p2c_npass    += gp(:Npassive_patch) * w
        p2c_nact     += gp(:Nactive_patch) * w
        p2c_nnonmyc  += gp(:Nnonmyc_patch) * w
        p2c_retr     += gp(:retransn_to_npool_patch) * w
        p2c_freeretr += gp(:free_retransn_to_npool_patch) * w
    end
    # vr soil debit (what nitrif_denitrif removes from smin pools) summed over layers
    debit_vr = 0.0
    let a3 = snf.sminn_to_plant_fun_no3_vr_col, a4 = snf.sminn_to_plant_fun_nh4_vr_col
        if !isempty(a3)
            for j in 1:length(dzsoi)
                debit_vr += (Float64(a3[c1,j]) + Float64(a4[c1,j])) * dzsoi[j]
            end
        end
    end

    println("\n=== FUN uptake decomposition (p2c, *dt gN/m2) ===")
    println("  sminn_to_plant_fun (veg uptake) = ", round(p2c_fun * dt, sigdigits=8))
    println("    Npassive                      = ", round(p2c_npass * dt, sigdigits=8))
    println("    Nactive(+nonmyc)              = ", round(p2c_nact * dt, sigdigits=8))
    println("    Nnonmyc                       = ", round(p2c_nnonmyc * dt, sigdigits=8))
    println("    Nfix                          = ", round(p2c_nfix * dt, sigdigits=8))
    println("  vr soil debit (passive+act+nm)  = ", round(debit_vr * dt, sigdigits=8))
    println("  veg uptake - soil debit         = ", round((p2c_fun - debit_vr) * dt, sigdigits=8), "  <- should equal Nfix")
    println("  retransn_to_npool (p2c)         = ", round(p2c_retr * dt, sigdigits=8))
    println("  free_retransn_to_npool (p2c)    = ", round(p2c_freeretr * dt, sigdigits=8))

    println("\n=== COLUMN $c1 (soil), one summer step, *dt in gN/m2 ===")
    println("  errnb_col                       = ", round(bal.errnb_col[c1], sigdigits=8))
    println("  --- the crux: veg uptake vs soil debit ---")
    println("  p2c(sminn_to_npool)*dt          = ", round(p2c_smin * dt, sigdigits=8))
    println("  sminn_to_plant_col*dt           = ", round(getcol(:sminn_to_plant_col) * dt, sigdigits=8))
    println("  DIFF (veg credit - soil debit)  = ", round((p2c_smin - getcol(:sminn_to_plant_col)) * dt, sigdigits=8))
    println("  --- N inputs*dt ---")
    for nm in (:ndep_to_sminn_col, :nfix_to_sminn_col, :ffix_to_sminn_col,
               :supplement_to_sminn_col, :fert_to_sminn_col, :soyfixn_to_sminn_col)
        println("    ", rpad(String(nm), 26), round(getcol(nm) * dt, sigdigits=8))
    end
    println("  --- N outputs*dt ---")
    for nm in (:denit_col, :sminn_leached_col, :smin_no3_leached_col,
               :smin_no3_runoff_col, :f_n2o_nit_col, :som_n_leached_col)
        println("    ", rpad(String(nm), 26), round(getcol(nm) * dt, sigdigits=8))
    end
    fnl = inst.bgc_vegetation.cnveg_nitrogenflux_inst.fire_nloss_col
    println("  fire_nloss_col*dt               = ", round(isempty(fnl) ? 0.0 : fnl[c1]*dt, sigdigits=8))
    return 0
end

exit(main())
