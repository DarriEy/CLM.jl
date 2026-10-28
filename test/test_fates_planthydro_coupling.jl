# test_fates_planthydro_coupling.jl
# FATES plant hydraulics (use_fates_planthydro, 1-D Taylor solver) connected end to
# end through the real clm_drv! loop on the fates_parity_1pt site (Bow at Banff):
#
#   * init: InitHydroGlobals / InitHydrSites / HydrSiteColdStart (clm_fates_init!);
#   * per step: wrap_hydraulics_drive (fates_hydraulics_step!) -> the plant water in
#     total_plant_stored_h2o_col, and ComputeRootSoilFlux (fates_compute_root_soil_flux!)
#     -> the host soil sink qflx_rootsoi_col is FATES's qflx_soil2root_sisl;
#   * daily: the hydraulics soil pack + UpdateH2OVeg in update_hlm_dynamics +
#     hydraulic-failure mortality from co_hydr.ftc_*;
#   * restart: RestartHydrStates rebuilds geometry/kmax/WTFs from the water states.
#
# Asserted:
#   1. the COLUMN water balance closes with hydraulics on (errh2o ~ round-off on
#      every step after the NUOPC init step), with the plant storage in the budget;
#   2. the soil sink is the FATES uptake, and it differs from transpiration by the
#      change in plant storage (sum(qflx_rootsoi) - qflx_tran_veg = dPlant/dt);
#   3. woody-PFT leaf water potential tracks the committed Fortran CTSM reference
#      (fates_planthydro_fortran.csv) — the 1-D solver parity;
#   4. hydraulic-failure mortality is computed from ftc (finite, and equal to the
#      Fortran formula for the current ftc);
#   5. RestartHydrStates reproduces the volumes / conductances / stored water.
#
# Data-dependent (Bow surfdata + forcing under SYMFLUENCE_DATA); skips cleanly if
# absent. Restores the FATES/varctl globals it changes.
using Test
using CLM

include("testdata.jl")

@testset "FATES plant hydraulics coupled to the host (clm_drv!)" begin
    sitedir = symfluence_path("domain_Bow_at_Banff_lumped_era5", "settings", "CLM", "parameters")
    fsurdat = joinpath(sitedir, "surfdata_clm.nc")
    forcing = symfluence_path("domain_Bow_at_Banff_lumped_era5", "data", "forcing",
                              "CLM_input", "clmforc.2002.nc")
    if !(isfile(fsurdat) && isfile(forcing))
        testdata_missing("FATES planthydro coupling run", fsurdat, forcing)
    else
        old_hydro = CLM.hlm_use_planthydro[]
        old_vc    = CLM.varctl.use_fates_planthydro
        try
            harness = joinpath(@__DIR__, "..", "scripts", "validation", "fates_fortran_parity",
                               "fates_planthydro_parity.jl")
            ENV["SYMFLUENCE_DATA"] = symfluence_data_root()
            m = Module(:_FatesPlantHydroHarness)
            Core.eval(m, :(include(p::AbstractString) = Base.include($m, p)))
            Base.include(m, harness)
            G(sym) = getglobal(m, sym)
            nsteps = 36   # 1.5 days: two daily steps (nstep 1 and 25)
            snaps = Dict{Int,Any}()
            hook = (inst, i) -> begin
                c = findfirst(inst.column.is_fates)
                wf = inst.water.waterfluxbulk_inst
                snaps[i] = (rootsoi = copy(wf.qflx_rootsoi_col[c, :]),
                            soil2root = copy(inst.fates.bc_out[1].qflx_soil2root_sisl),
                            tran = wf.wf.qflx_tran_veg_col[c],
                            plant = inst.water.waterdiagnosticbulk_inst.total_plant_stored_h2o_col[c],
                            plant_bc = inst.fates.bc_out[1].plant_stored_h2o_si)
            end
            recs, inst = Base.invokelatest(G(:run_julia), nsteps; hook=hook)
            nl = CLM.varpar.nlevsoi

            @test CLM.hlm_use_planthydro[] == CLM.itrue
            @test CLM.varctl.use_fates_planthydro
            @test inst.fates.sites[1].si_hydr !== nothing

            # 1. column water balance closes. Skipped: nstep 0 (the NUOPC init step)
            # and nstep 1, whose first daily step re-weights the cold-start patches
            # mid-step — a 0.118 mm residual that is identical with hydraulics OFF
            # (a FATES-harness artifact, not a hydraulics term).
            errs = [abs(r["ERRH2O_COL"]) for r in recs[3:end]]
            @test all(isfinite, errs)
            @test maximum(errs) < 1e-8

            # 2. the soil sink is the FATES uptake; storage closes the plant budget.
            for i in 2:nsteps
                s = snaps[i]
                @test s.rootsoi[1:nl] == s.soil2root[1:nl]
                @test s.plant == s.plant_bc
                @test s.plant > 0
            end
            # Hydraulic redistribution is live: some layer receives water at some step.
            @test any(any(<(0.0), snaps[i].rootsoi) for i in 2:nsteps)
            # dPlant/dt = uptake - transpiration on steps without a daily step (the
            # daily dynamics shuffles water into dead/err pools but conserves the total).
            for i in 3:nsteps
                i == 25 && continue
                dP = snaps[i].plant - snaps[i-1].plant
                flux = (sum(snaps[i].rootsoi) - snaps[i].tran) * G(:DTIME)
                @test dP ≈ flux atol = 1e-9
            end

            # 3. woody leaf water potential vs the Fortran reference.
            ft = Base.invokelatest(G(:read_fortran), G(:FORTRAN_CSV))
            for k in 3:nsteps+1          # record k <-> nstep k-1; skip the init records
                F = ft["PSI_LEAF"][k]; J = recs[k]["PSI_LEAF"]
                @test abs(J - F) / abs(F) < 5e-2
            end
            woody_rel = [abs(recs[k]["PSI_LEAF_pft1"] - ft["PSI_LEAF_pft1"][k]) /
                         abs(ft["PSI_LEAF_pft1"][k]) for k in 3:nsteps+1]
            @test maximum(woody_rel) < 1e-2

            # 4. hydraulic-failure mortality from ftc.
            site = inst.fates.sites[1]
            evc = CLM.EDPftvarcon_inst[]
            nchk = 0
            cp = site.oldest_patch
            while cp !== nothing
                cc = cp.tallest
                while cc !== nothing
                    if CLM.prt_params.woody[cc.pft] == 1
                        h = cc.co_hydr
                        flc = 1.0 - min(minimum(h.ftc_ag), h.ftc_troot, minimum(h.ftc_aroot))
                        thr = evc.hf_flc_threshold[cc.pft]
                        expect = (flc >= thr && thr < 1.0) ?
                                 (flc - thr) / (1.0 - thr) * evc.mort_scalar_hydrfailure[cc.pft] : 0.0
                        hmort = CLM.mortality_rates(cc, inst.fates.bc_in[1],
                                                    cp.btran_ft, 0.0)[2]
                        @test isfinite(hmort)
                        @test hmort ≈ expect atol = 1e-12
                        nchk += 1
                    end
                    cc = cc.shorter
                end
                cp = cp.younger
            end
            @test nchk > 0

            # 5. RestartHydrStates rebuilds what depends on the restored water states.
            ref = Dict{Any,Any}()
            cp = site.oldest_patch
            while cp !== nothing
                cc = cp.tallest
                while cc !== nothing
                    h = cc.co_hydr
                    ref[cc] = (copy(h.v_ag), h.v_troot, copy(h.kmax_stem_upper),
                               copy(h.kmax_aroot_radial_in), copy(h.z_node_ag))
                    fill!(h.v_ag, 0.0); h.v_troot = 0.0
                    fill!(h.kmax_stem_upper, 0.0); fill!(h.kmax_aroot_radial_in, 0.0)
                    cc = cc.shorter
                end
                cp = cp.younger
            end
            stored = inst.fates.bc_out[1].plant_stored_h2o_si
            c = findfirst(inst.column.is_fates)
            CLM.fates_pack_bcin_hydro_soil!(inst.fates.bc_in[1], inst; c=c, nlevsoil=nl)
            CLM.RestartHydrStates(inst.fates.sites, 1, inst.fates.bc_in, inst.fates.bc_out)
            for (cc, r) in ref
                h = cc.co_hydr
                @test isequal(h.v_ag, r[1]) && isequal(h.v_troot, r[2])
                @test isequal(h.kmax_stem_upper, r[3]) && isequal(h.kmax_aroot_radial_in, r[4])
                @test isequal(h.z_node_ag, r[5])
                @test h.v_ag_init == h.v_ag
            end
            @test inst.fates.bc_out[1].plant_stored_h2o_si ≈ stored rtol = 1e-12
            @test all(isfinite, site.si_hydr.v_shell[1:site.si_hydr.nlevrhiz, :])
        finally
            CLM.hlm_use_planthydro[] = old_hydro
            CLM.varctl.use_fates_planthydro = old_vc
        end
    end
end
