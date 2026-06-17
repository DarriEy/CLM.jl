# ==========================================================================
# test_fortran_parity_cn.jl — gated CN/BGC/PHS/LUNA/FUN parity regression.
#
# Companion to test_fortran_parity.jl (the SP, use_cn=false biogeophysics
# path). This file guards the use_cn=true BGC path that we verified by hand
# against the instrumented Fortran run (PHS + LUNA + FUN, summer Bow window):
#
#   1. SINGLE-STEP parity (nstep 1757852): inject the before_step BGC dump as
#      the shared IC, run ONE use_cn clm_drv! step, and diff the evolved Julia
#      state to the after_hydrologydrainage dump. Verified ground truth:
#        - mineral N (sminn_vr / smin_nh4_vr) per-layer rel err < 0.1%
#          (worst layer was 0.024%),
#        - FUN plant-N uptake per patch within ~2% (tree 1.002×, grass 0.988×),
#        - availc per patch within ~2%.
#
#   2. MULTI-STEP DRIFT (free-running, IC injected ONCE): advance forcing+time
#      and diff the Julia trajectory against the per-step Fortran dumps. The
#      mineral-N rel err grows ~linearly and stays bounded (<1.5% by step 28);
#      long-memory pools (soil1c, leafc) stay <0.2%; nothing runs away.
#
# These tolerances are set to PASS the current code with ~3-5× headroom over
# the verified numbers, yet FAIL if any of the CN parity fixes regress:
#   - LUNA vcmax/jmax injection, CO2/pbot, Jackson rootfr, hksat percolation,
#   - FUN persist/inject leafc_to_litter, br_root / Atkin / light_inhibit /
#     modifyphoto_and_lmr_forcrop leaf-respiration wiring,
#   - the smin N-balance (nitrif/denitrif/leaching/competition) chain.
#
# GATED: the Fortran BGC reference dumps + Bow input files live outside the
# repo (machine-specific). When absent, the test is skipped (not failed) so CI
# and other machines stay green.
# ==========================================================================

using Test
using CLM
using NCDatasets
using Dates

@testset "Fortran CN/BGC parity (gated)" begin
    common  = joinpath(@__DIR__, "..", "scripts", "fortran_parity_common.jl")
    dumpdir = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_summer"

    # The BGC run cycles datm 2002-2009 (year_align=2002): model year 2202 maps
    # to forcing year 2002. nstep -> date base DateTime(2002,1,1)+Hour(n-1753153).
    date_base = 1753153

    # Single-step check uses one daytime summer step; the drift check free-runs a
    # contiguous window of dumps starting at drift_n0.
    nstep_single = 1757852
    drift_n0     = 1757845
    drift_nsteps = 12          # keep runtime reasonable (28 dumps available)

    forcing2002 = nothing      # resolved after `common` is loaded (needs FFORCING)

    have_single = isfile(joinpath(dumpdir, "pdump_before_step_n$(nstep_single).nc")) &&
                  isfile(joinpath(dumpdir, "pdump_after_hydrologydrainage_n$(nstep_single).nc"))
    have_drift = all(
        isfile(joinpath(dumpdir, "pdump_before_step_n$(drift_n0).nc")) &&
        isfile(joinpath(dumpdir, "pdump_after_hydrologydrainage_n$(drift_n0 + i).nc"))
        for i in 0:(drift_nsteps - 1))
    have_dumps = have_single && have_drift

    if !(isfile(common) && have_dumps)
        @info "Fortran CN/BGC parity test SKIPPED (BGC reference dumps not present on this machine)"
        @test_skip have_dumps
    else
        Base.include(@__MODULE__, common)   # build_bow_inst / run_one_parity_step! / FFORCING / _dumpvar
        forcing2002 = replace(FFORCING, "clmforc.2003.nc" => "clmforc.2002.nc")

        # --- per-layer mineral-N relative error vs a dump var (col 1, levdecomp) --
        # Inline the missing->NaN read (instead of calling common's _dumpvar) to
        # keep the comparison self-contained. NOTE: read nlevdecomp at CALL time —
        # it is set by clm_initialize! (run inside the step), so it is still 0 at
        # this point in module load.
        dumpvec(ds, name) = haskey(ds, name) ?
            [ismissing(v) ? NaN : Float64(v) for v in ds[name][:]] : nothing
        minn_relerr(jlcol, ds, name) = begin
            nl = CLM.varpar.nlevdecomp
            d = dumpvec(ds, name); d === nothing && return NaN
            m = 0.0
            for j in 1:min(nl, length(d), size(jlcol, 2))
                f = d[j]; v = Float64(jlcol[1, j])
                (isnan(f) || isnan(v)) && continue
                den = abs(f) > 1e-12 ? abs(f) : 1.0
                m = max(m, abs(v - f) / den)
            end
            m
        end

        # --------------------------------------------------------------------
        # 1. SINGLE-STEP CN parity
        # --------------------------------------------------------------------
        @testset "single-step n$nstep_single" begin
            inst, bounds = run_one_parity_step!(nstep_single; use_cn=true, dumpdir=dumpdir,
                use_hydrstress=true, use_luna=true,
                step_date=DateTime(2002,1,1) + Hour(nstep_single - date_base),
                forcing_file=forcing2002)

            sns = inst.soilbiogeochem_nitrogenstate
            cnf = inst.bgc_vegetation.cnveg_nitrogenflux_inst
            cf  = inst.bgc_vegetation.cnveg_carbonflux_inst

            ds = NCDataset(joinpath(dumpdir, "pdump_after_hydrologydrainage_n$(nstep_single).nc"), "r")

            # (a) mineral N per-layer rel err < 1e-3 (verified worst 0.024%; 4× headroom)
            sminn_err   = minn_relerr(sns.sminn_vr_col,    ds, "sminn_vr")
            smin_nh4_err = minn_relerr(sns.smin_nh4_vr_col, ds, "smin_nh4_vr")
            @info "single-step mineral-N max rel err" sminn=sminn_err smin_nh4=smin_nh4_err
            @test sminn_err   < 1e-3
            @test smin_nh4_err < 1e-3

            # patch remap (Fortran np may differ from Julia np; match by PFT itype)
            f_pfts  = Int.(ds["pfts1d_itypveg"][:])
            jl_pfts = Int.(inst.patch.itype)
            f_of(pj) = (it = jl_pfts[pj]; it == 0 ? nothing : findfirst(==(it), f_pfts))

            # (b) FUN plant-N uptake per patch within 0.03 of 1.0 (verified tree
            #     1.002× / grass 0.988×). Total uptake = NH4 + NO3 active+nonmyc.
            funtot(arrsym, p) = Float64(getfield(cnf, Symbol(arrsym))[p])
            f_dump(name) = haskey(ds, name) ? Float64.(ds[name][:]) : nothing
            f_Nactive_nh4 = f_dump("Nactive_nh4"); f_Nnonmyc_nh4 = f_dump("Nnonmyc_nh4")
            f_Nactive_no3 = f_dump("Nactive_no3"); f_Nnonmyc_no3 = f_dump("Nnonmyc_no3")
            for pj in 1:length(jl_pfts)
                pf = f_of(pj); pf === nothing && continue
                (f_Nactive_nh4 === nothing) && continue
                jl_up = funtot("Nactive_nh4_patch", pj) + funtot("Nnonmyc_nh4_patch", pj) +
                        funtot("Nactive_no3_patch", pj) + funtot("Nnonmyc_no3_patch", pj)
                f_up  = f_Nactive_nh4[pf] + f_Nnonmyc_nh4[pf] +
                        f_Nactive_no3[pf] + f_Nnonmyc_no3[pf]
                (abs(f_up) < 1e-20) && continue
                ratio = jl_up / f_up
                @info "FUN N-uptake ratio" patch=pj pft=jl_pfts[pj] ratio=ratio
                @test abs(ratio - 1.0) < 0.03
            end

            # (c) availc per patch within 0.03 of 1.0 (verified ~2%)
            availc_f = f_dump("availc")
            if availc_f !== nothing
                for pj in 1:length(jl_pfts)
                    pf = f_of(pj); pf === nothing && continue
                    avc_f = availc_f[pf]
                    (abs(avc_f) < 1e-12) && continue
                    ratio = Float64(cf.availc_patch[pj]) / avc_f
                    @info "availc ratio" patch=pj pft=jl_pfts[pj] ratio=ratio
                    @test abs(ratio - 1.0) < 0.03
                end
            end
            close(ds)
        end

        # --------------------------------------------------------------------
        # 2. MULTI-STEP DRIFT (free-running): inject IC once, advance, diff each
        #    step. Mirrors scripts/fortran_parity_drift.jl's loop body.
        # --------------------------------------------------------------------
        @testset "multi-step drift (n$drift_n0 +$(drift_nsteps))" begin
            start_date = DateTime(2002,1,1) + Hour(drift_n0 - date_base)
            (inst, bounds, filt, tm) = build_bow_inst(; dtime=3600,
                                                       start_date=start_date, use_cn=true)
            if isempty(inst.photosyns.vcmx25_z_patch)
                inst.photosyns.vcmx25_z_patch = fill(30.0, bounds.endp, CLM.NLEVCAN)
                inst.photosyns.jmx25_z_patch  = fill(60.0, bounds.endp, CLM.NLEVCAN)
            end
            icdump = joinpath(dumpdir, "pdump_before_step_n$(drift_n0).nc")
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
            ng, nc, np = bounds.endg, bounds.endc, bounds.endp
            fr = CLM.ForcingReader(); CLM.forcing_reader_init!(fr, forcing2002)
            tf = replace(forcing2002, r"clmforc\.[^/]*\.nc$" => "topo_forcing.nc")
            if isfile(tf)
                dt = NCDataset(tf, "r")
                if haskey(dt, "TOPO")
                    ft = Float64(dt["TOPO"][1])
                    for g in 1:ng; inst.atm2lnd.forc_topo_grc[g] = ft; end
                    for c in 1:nc; inst.topo.topo_col[c]        = ft; end
                end
                close(dt)
            end

            sbns = inst.soilbiogeochem_nitrogenstate
            sbcs = inst.soilbiogeochem_carbonstate
            cs   = inst.bgc_vegetation.cnveg_carbonstate_inst

            sminn_drift = Float64[]; nh4_drift = Float64[]
            leafc_drift = Float64[]; soil1c_drift = Float64[]
            for i in 0:(drift_nsteps - 1)
                cur = start_date + Hour(i)
                calday = CLM.get_curr_calday(tm); (declin, _) = CLM.compute_orbital(calday)
                nextsw = calday + 3600.0 / CLM.SECSPDAY
                (declinm1, _) = CLM.compute_orbital(calday - 3600.0 / CLM.SECSPDAY)
                CLM.init_daylength!(inst.gridcell, declin, declinm1, CLM.ORB_OBLIQR_DEFAULT, 1:bounds.endg)
                CLM.advance_timestep!(tm)
                CLM.read_forcing_step!(fr, inst.atm2lnd, cur, ng, nc)
                CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
                (yr, mon, d, tod) = CLM.get_curr_date(tm)
                CLM.clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw, declin, declin,
                    CLM.ORB_OBLIQR_DEFAULT, false, false, "", false;
                    nstep=tm.nstep, is_first_step=false,
                    is_beg_curr_day=CLM.is_beg_curr_day(tm), is_end_curr_day=CLM.is_end_curr_day(tm),
                    is_beg_curr_year=CLM.is_beg_curr_year(tm), dtime=3600.0, mon=mon, day=d,
                    photosyns=inst.photosyns)

                n = drift_n0 + i
                da = NCDataset(joinpath(dumpdir, "pdump_after_hydrologydrainage_n$(n).nc"), "r")
                push!(sminn_drift, minn_relerr(sbns.sminn_vr_col,    da, "sminn_vr"))
                push!(nh4_drift,   minn_relerr(sbns.smin_nh4_vr_col, da, "smin_nh4_vr"))
                # leafc: patches 2,3 (tree,grass); slot 1 is bare
                lc = dumpvec(da, "leafc")
                if lc !== nothing
                    e = 0.0
                    for p in 2:min(3, length(lc), length(cs.leafc_patch))
                        f = lc[p]; v = Float64(cs.leafc_patch[p])
                        (isnan(f) || isnan(v) || abs(f) < 1e-12) && continue
                        e = max(e, abs(v - f) / abs(f))
                    end
                    push!(leafc_drift, e)
                end
                push!(soil1c_drift, minn_relerr(sbcs.decomp_cpools_vr_col[:, :, 4], da, "soil1c_vr"))
                close(da)
            end
            CLM.forcing_reader_close!(fr)

            sminn_final = sminn_drift[end]
            sminn_mid   = sminn_drift[max(1, cld(length(sminn_drift), 2))]
            @info "drift summary" steps=drift_nsteps sminn_final=sminn_final sminn_mid=sminn_mid nh4_final=nh4_drift[end] leafc_final=(isempty(leafc_drift) ? NaN : leafc_drift[end]) soil1c_final=soil1c_drift[end]

            # mineral N stays bounded: < 1.5% by the (short) window end
            @test sminn_final   < 1.5e-2
            @test nh4_drift[end] < 1.5e-2

            # not super-linear: final/midpoint ratio roughly consistent with linear
            # growth (would be ~2 for exact linear over a halved window; allow 2.5×)
            if sminn_mid > 1e-9
                @test (sminn_final / sminn_mid) < 2.5
            end

            # long-memory pools barely move over a short window
            @test soil1c_drift[end] < 2e-3
            isempty(leafc_drift) || @test leafc_drift[end] < 2e-3
        end
    end
end
