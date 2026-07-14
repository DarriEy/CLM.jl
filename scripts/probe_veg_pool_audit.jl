# ==========================================================================
# probe_veg_pool_audit.jl — WHICH VEG POOL GAINS CARBON FROM NOWHERE?
#
# After the lf_f and phenology-nlevdecomp fixes the veg->soil litter transfer
# conserves exactly, and the soil side closes. The residual is now entirely
# a VEG-side gain: d(totc_patch) exceeds (gpp - ar - litterfall).
#
# Snapshot every patch-level C pool across one step and compare each pool's
# delta with the fluxes that are supposed to move it.
#
#   julia +1.12 --project=. scripts/probe_veg_pool_audit.jl
# ==========================================================================
include(joinpath(@__DIR__, "fortran_parity_common.jl"))

const SUMDIR = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_summer"
const N0 = 1757845

const POOLS = (:cpool, :xsmrpool, :ctrunc, :leafc, :leafc_storage, :leafc_xfer,
    :frootc, :frootc_storage, :frootc_xfer, :livestemc, :livestemc_storage,
    :livestemc_xfer, :deadstemc, :deadstemc_storage, :deadstemc_xfer,
    :livecrootc, :livecrootc_storage, :livecrootc_xfer, :deadcrootc,
    :deadcrootc_storage, :deadcrootc_xfer, :gresp_storage, :gresp_xfer)

snap(cs, p) = Dict(k => Float64(getfield(cs, Symbol(String(k) * "_patch"))[p]) for k in POOLS)

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

    cs  = inst.bgc_vegetation.cnveg_carbonstate_inst
    cvf = inst.bgc_vegetation.cnveg_carbonflux_inst
    dt = 3600.0

    befores = Dict(p => snap(cs, p) for p in 2:3)
    step!(6)
    afters = Dict(p => snap(cs, p) for p in 2:3)
    CLM.forcing_reader_close!(fr)

    for p in 2:3
        b = befores[p]; a = afters[p]
        println("\n=== PATCH $p (wt=", round(inst.patch.wtcol[p], sigdigits=4),
                ", ivt=", inst.patch.itype[p], ") ===")
        dtot = sum(a[k] - b[k] for k in POOLS)
        gpp  = cvf.gpp_patch[p] * dt
        ar   = cvf.ar_patch[p] * dt
        lit  = cvf.litfall_patch[p] * dt
        fire = cvf.fire_closs_patch[p] * dt
        println("  d(sum of pools)  = ", round(dtot, sigdigits=8))
        println("  gpp              = ", round(gpp, sigdigits=8))
        println("  ar               = ", round(ar, sigdigits=8))
        println("  litfall          = ", round(lit, sigdigits=8))
        println("  fire_closs       = ", round(fire, sigdigits=8))
        println("  gpp-ar-lit-fire  = ", round(gpp - ar - lit - fire, sigdigits=8))
        println("  RESIDUAL         = ", round(dtot - (gpp - ar - lit - fire), sigdigits=8))
        println("  --- per-pool deltas (nonzero) ---")
        for k in POOLS
            d = a[k] - b[k]
            abs(d) < 1e-14 && continue
            println("    ", rpad(String(k), 22), round(d, sigdigits=8))
        end
        println("  --- selected fluxes*dt ---")
        for nm in (:psnsun_to_cpool, :psnshade_to_cpool, :cpool_to_resp, :soilc_change,
                   :cpool_to_xsmrpool, :xsmrpool_recover, :leafc_to_litter,
                   :frootc_to_litter, :leafc_to_litter_fun, :excess_cflux,
                   :plant_calloc, :availc, :mr, :gr, :gpp_before_downreg,
                   # MR is split into a "current" part paid from cpool (curmr, via
                   # cpool_to_resp) and an "excess" part paid from xsmrpool (xsmr).
                   # Every gC of MR is counted in ar; it MUST be debited from one
                   # of those two pools or it is created from nothing.
                   :leaf_mr, :froot_mr, :livestem_mr, :livecroot_mr,
                   :leaf_curmr, :froot_curmr, :livestem_curmr, :livecroot_curmr,
                   :leaf_xsmr, :froot_xsmr, :livestem_xsmr, :livecroot_xsmr)
            f = getfield(cvf, Symbol(String(nm) * "_patch"))
            isempty(f) && continue
            v = Float64(f[p]) * dt
            abs(v) < 1e-14 && continue
            println("    ", rpad(String(nm), 22), round(v, sigdigits=8))
        end

        # --- CPOOL LEDGER: every debit c_state_update1! applies, vs the observed delta.
        woody = (inst.patch.itype[p] in (1,2,3,4,5,6,7,8))
        fx(nm) = (arr = getfield(cvf, Symbol(String(nm) * "_patch")); isempty(arr) ? 0.0 : Float64(arr[p]) * dt)
        deb = fx(:cpool_to_xsmrpool) + fx(:leaf_curmr) + fx(:froot_curmr) +
              fx(:cpool_to_resp) + fx(:soilc_change) +
              fx(:cpool_to_leafc) + fx(:cpool_to_leafc_storage) +
              fx(:cpool_to_frootc) + fx(:cpool_to_frootc_storage) +
              fx(:cpool_leaf_gr) + fx(:cpool_froot_gr) +
              fx(:cpool_leaf_storage_gr) + fx(:cpool_froot_storage_gr)
        if woody
            deb += fx(:livestem_curmr) + fx(:livecroot_curmr) +
                   fx(:cpool_to_livestemc) + fx(:cpool_to_livestemc_storage) +
                   fx(:cpool_to_deadstemc) + fx(:cpool_to_deadstemc_storage) +
                   fx(:cpool_to_livecrootc) + fx(:cpool_to_livecrootc_storage) +
                   fx(:cpool_to_deadcrootc) + fx(:cpool_to_deadcrootc_storage) +
                   fx(:cpool_livestem_gr) + fx(:cpool_deadstem_gr) +
                   fx(:cpool_livecroot_gr) + fx(:cpool_deadcroot_gr) +
                   fx(:cpool_livestem_storage_gr) + fx(:cpool_deadstem_storage_gr) +
                   fx(:cpool_livecroot_storage_gr) + fx(:cpool_deadcroot_storage_gr)
        end
        gin = fx(:psnsun_to_cpool) + fx(:psnshade_to_cpool)
        dcp = a[:cpool] - b[:cpool]
        println("  --- CPOOL LEDGER (woody=", woody, ") ---")
        println("    gpp into cpool        = ", round(gin, sigdigits=8))
        println("    debits applied by code= ", round(deb, sigdigits=8))
        println("    expected d(cpool)     = ", round(gin - deb, sigdigits=8))
        println("    observed d(cpool)     = ", round(dcp, sigdigits=8))
        println("    UNEXPLAINED           = ", round(dcp - (gin - deb), sigdigits=8))
        println("    cpool_to_gresp_storage= ", round(fx(:cpool_to_gresp_storage), sigdigits=8),
                "   <-- Fortran debits cpool by this (CNCStateUpdate1Mod.F90:501); Julia does NOT")
        println("    d(gresp_storage)      = ", round(a[:gresp_storage] - b[:gresp_storage], sigdigits=8))
        println("    d(gresp_xfer)         = ", round(a[:gresp_xfer] - b[:gresp_xfer], sigdigits=8))
    end
    return 0
end

exit(main())
