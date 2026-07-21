# ==========================================================================
# probe_veg_pool_audit_n.jl — WHICH VEG N POOL GAINS NITROGEN FROM NOWHERE?
#
# N analogue of probe_veg_pool_audit.jl. The veg->soil litter transfer
# conserves; the residual (~3.1e-4 gN/m2/step, summer) is a VEG-side gain:
# the veg N pools grow by more than (sminn_to_npool - litterfall) supplies.
#
# Snapshot every patch-level N pool across one step and attribute each pool's
# delta to the fluxes that move it. Then run an NPOOL / RETRANSN LEDGER.
#
#   julia +1.12 --project=. scripts/probe_veg_pool_audit_n.jl
# ==========================================================================
include(joinpath(@__DIR__, "fortran_parity_common.jl"))

const SUMDIR = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_summer"
const N0 = 1757845

const NPOOLS = (:leafn, :leafn_storage, :leafn_xfer,
    :frootn, :frootn_storage, :frootn_xfer,
    :livestemn, :livestemn_storage, :livestemn_xfer,
    :deadstemn, :deadstemn_storage, :deadstemn_xfer,
    :livecrootn, :livecrootn_storage, :livecrootn_xfer,
    :deadcrootn, :deadcrootn_storage, :deadcrootn_xfer,
    :retransn, :npool, :ntrunc)

snap(ns, p) = Dict(k => Float64(getfield(ns, Symbol(String(k) * "_patch"))[p]) for k in NPOOLS)

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
    for i in 0:5; step!(i); end

    ns  = inst.bgc_vegetation.cnveg_nitrogenstate_inst
    nvf = inst.bgc_vegetation.cnveg_nitrogenflux_inst
    dt = 3600.0

    befores = Dict(p => snap(ns, p) for p in 2:3)
    step!(6)
    afters = Dict(p => snap(ns, p) for p in 2:3)
    CLM.forcing_reader_close!(fr)

    fx(p, nm) = (arr = getfield(nvf, Symbol(String(nm) * "_patch")); (isempty(arr) || ndims(arr) != 1) ? 0.0 : Float64(arr[p]) * dt)

    for p in 2:3
        b = befores[p]; a = afters[p]
        println("\n=== PATCH $p (wt=", round(inst.patch.wtcol[p], sigdigits=4),
                ", ivt=", inst.patch.itype[p], ") ===")
        dtot = sum(a[k] - b[k] for k in NPOOLS)
        smin = fx(p, :sminn_to_npool)
        lit  = fx(p, :leafn_to_litter) + fx(p, :frootn_to_litter)
        println("  d(sum of N pools)   = ", round(dtot, sigdigits=8))
        println("  sminn_to_npool      = ", round(smin, sigdigits=8), "   <- only external INPUT to veg N")
        println("  leaf/froot->litter  = ", round(lit, sigdigits=8), "   <- external OUTPUT to soil")
        println("  smin - litter       = ", round(smin - lit, sigdigits=8))
        println("  RESIDUAL (dtot-that)= ", round(dtot - (smin - lit), sigdigits=8))
        println("  --- per-pool deltas (nonzero) ---")
        for k in NPOOLS
            d = a[k] - b[k]
            abs(d) < 1e-14 && continue
            println("    ", rpad(String(k), 22), round(d, sigdigits=8))
        end

        # --- NPOOL LEDGER: every credit/debit n_state_update1! applies ---
        woody = (inst.patch.itype[p] in (1,2,3,4,5,6,7,8))
        credits = fx(p, :sminn_to_npool) + fx(p, :retransn_to_npool) + fx(p, :free_retransn_to_npool)
        debits  = fx(p, :npool_to_leafn) + fx(p, :npool_to_leafn_storage) +
                  fx(p, :npool_to_frootn) + fx(p, :npool_to_frootn_storage)
        if woody
            debits += fx(p, :npool_to_livestemn) + fx(p, :npool_to_livestemn_storage) +
                      fx(p, :npool_to_deadstemn) + fx(p, :npool_to_deadstemn_storage) +
                      fx(p, :npool_to_livecrootn) + fx(p, :npool_to_livecrootn_storage) +
                      fx(p, :npool_to_deadcrootn) + fx(p, :npool_to_deadcrootn_storage)
        end
        dnp = a[:npool] - b[:npool]
        println("  --- NPOOL LEDGER (woody=", woody, ") ---")
        println("    credits (smin+retrans+free) = ", round(credits, sigdigits=8))
        println("    debits (npool_to_*)         = ", round(debits, sigdigits=8))
        println("    expected d(npool)           = ", round(credits - debits, sigdigits=8))
        println("    observed d(npool)           = ", round(dnp, sigdigits=8))
        println("    UNEXPLAINED                 = ", round(dnp - (credits - debits), sigdigits=8))

        # --- RETRANSN LEDGER ---
        rt_cred = fx(p, :leafn_to_retransn) + fx(p, :livestemn_to_retransn) +
                  fx(p, :livecrootn_to_retransn) + fx(p, :frootn_to_retransn)
        rt_deb  = fx(p, :retransn_to_npool) + fx(p, :free_retransn_to_npool)
        drt = a[:retransn] - b[:retransn]
        println("  --- RETRANSN LEDGER ---")
        println("    credits (leaf+wood->retrans)= ", round(rt_cred, sigdigits=8))
        println("    debits (->npool + free)     = ", round(rt_deb, sigdigits=8))
        println("    expected d(retransn)        = ", round(rt_cred - rt_deb, sigdigits=8))
        println("    observed d(retransn)        = ", round(drt, sigdigits=8))
        println("    UNEXPLAINED                 = ", round(drt - (rt_cred - rt_deb), sigdigits=8))

        println("  --- selected N fluxes*dt (nonzero) ---")
        for nm in (:sminn_to_npool, :retransn_to_npool, :free_retransn_to_npool,
                   :leafn_to_retransn, :livestemn_to_retransn, :livecrootn_to_retransn,
                   :frootn_to_retransn, :leafn_to_litter, :frootn_to_litter,
                   :npool_to_leafn, :npool_to_leafn_storage, :npool_to_frootn, :npool_to_frootn_storage,
                   :npool_to_livestemn, :npool_to_livestemn_storage,
                   :npool_to_deadstemn, :npool_to_deadstemn_storage,
                   :npool_to_livecrootn, :npool_to_livecrootn_storage,
                   :npool_to_deadcrootn, :npool_to_deadcrootn_storage,
                   :leafn_xfer_to_leafn, :frootn_xfer_to_frootn,
                   :leafn_storage_to_xfer, :frootn_storage_to_xfer,
                   :livestemn_to_deadstemn, :livecrootn_to_deadcrootn,
                   :supplement_to_plantn)
            v = fx(p, nm)
            abs(v) < 1e-14 && continue
            println("    ", rpad(String(nm), 26), round(v, sigdigits=8))
        end
    end
    return 0
end

exit(main())
