# ==========================================================================
# probe_n_split_vegsoil.jl — VEG vs SOIL: which side gains the created N?
#
# Snapshot, across one summer step, the column's total N split into:
#   veg (p2c of totn_patch), soil-mineral (sminn_col), soil-decomp
#   (decomp_npools_col), soil-truncation (ntrunc_col), and totn_col.
# Then attribute d(soil-mineral) and d(soil-decomp) to their fluxes.
#
#   julia +1.12 --project=. scripts/probe_n_split_vegsoil.jl
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

    sns = inst.soilbiogeochem_nitrogenstate
    vns = inst.bgc_vegetation.cnveg_nitrogenstate_inst
    soil = [c for c in 1:nc if filt.bgc_soilc[c]]
    c1 = soil[1]

    vegn(c) = begin
        s = 0.0
        for p in bounds.begp:bounds.endp
            (filt.bgc_vegp[p]) || continue
            w = inst.patch.wtcol[p]; isfinite(w) || continue
            (inst.patch.column[p] == c) || continue
            v = vns.totn_patch[p]
            isfinite(v) && (s += v * w)
        end
        s
    end
    soilmin(c)   = (isempty(sns.sminn_col) ? 0.0 : Float64(sns.sminn_col[c]))
    soildecomp(c)= (isempty(sns.decomp_npools_col) ? 0.0 : sum(Float64, @view sns.decomp_npools_col[c, :]))
    soiltrunc(c) = (isempty(sns.ntrunc_col) ? 0.0 : Float64(sns.ntrunc_col[c]))
    totn(c)      = (isempty(sns.totn_col) ? 0.0 : Float64(sns.totn_col[c]))

    for i in 0:5; step!(i); end
    b = (veg=vegn(c1), min=soilmin(c1), dec=soildecomp(c1), tr=soiltrunc(c1), tot=totn(c1))
    step!(6)
    a = (veg=vegn(c1), min=soilmin(c1), dec=soildecomp(c1), tr=soiltrunc(c1), tot=totn(c1))
    CLM.forcing_reader_close!(fr)

    println("\n=== COLUMN $c1: N split before/after one summer step (gN/m2) ===")
    for (nm, k) in (("veg (p2c totn_patch)", :veg), ("soil mineral (sminn)", :min),
                    ("soil decomp pools", :dec), ("soil ntrunc", :tr), ("TOTN_COL", :tot))
        println("  ", rpad(nm, 24), " d = ", round(getfield(a,k) - getfield(b,k), sigdigits=8),
                "   (", round(getfield(b,k), sigdigits=6), " -> ", round(getfield(a,k), sigdigits=6), ")")
    end
    dsum = (a.veg-b.veg) + (a.min-b.min) + (a.dec-b.dec) + (a.tr-b.tr)
    println("  sum of parts d          = ", round(dsum, sigdigits=8))
    println("  errnb_col               = ", round(bal.errnb_col[c1], sigdigits=8))
    return 0
end

exit(main())
