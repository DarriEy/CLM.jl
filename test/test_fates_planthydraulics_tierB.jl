# test_fates_planthydraulics_tierB.jl
# Tests for the FATES plant-hydraulics "Tier B" transpiration/uptake solve added to
# FatesPlantHydraulicsMod.jl (SetMaxCondConnections / Report1DError / RecruitWUptake /
# FillDrainRhizShells / Hydraulics_BC / hydraulics_drive + UpdateH2OVeg!).
#
# These routines are gated behind `hlm_use_planthydro==itrue`. We flip the flag on,
# build a full single-woody-PFT FATES parameter table (allometry + hydraulics), a
# hydro-enabled site whose `si_hydr` is built through InitHydrSites/HydrSiteColdStart,
# the size-dependent rhizosphere geometry (UpdateSizeDepRhizHydProps), and one real-PRT
# cohort with hydro geometry from UpdateSizeDepPlantHydProps! hung off a vegetated patch.
# Then we drive the transpiration solve and assert:
#   * Hydraulics_BC / hydraulics_drive run and POPULATE co_hydr.ftc_ag/ftc_troot/ftc_aroot
#     ∈ [0,1], a finite co_hydr.btran ∈ [0,1], and a finite leaf_psi (psi_ag[1]).
#   * The water mass-balance closes (|errh2o_hyd| small) and the site plant-storage
#     check inside Hydraulics_BC does not throw.
#   * FillDrainRhizShells conserves rhizosphere water.
#   * RecruitWUptake returns finite recruit-water uptake.
# The `hlm_use_planthydro` flag (and the other touched globals) are restored at the end.

using Test
using CLM

const FPHB = CLM   # all FATES symbols live in the CLM module

@testset "FATES Plant Hydraulics — Tier B (transpiration solve)" begin

    # Save + later restore mutable global state we touch.
    saved_planthydro = CLM.hlm_use_planthydro[]
    saved_wrf        = CLM._wrf_plant[]
    saved_wkf        = CLM._wkf_plant[]
    saved_pftvarcon  = CLM.EDPftvarcon_inst[]
    saved_edparams   = CLM.EDParams[]
    saved_numpft     = CLM.numpft[]
    saved_nlevsclass = CLM.nlevsclass[]
    saved_prtglobal  = CLM.prt_global[]
    saved_numel      = CLM.num_elements[]
    saved_ellist     = copy(CLM.element_list)
    saved_parteh     = CLM.hlm_parteh_mode[]
    saved_paramderiv = CLM.ParamDerived[]

    try
        CLM.hlm_use_planthydro[] = CLM.itrue

        numpft = 1
        ft = 1
        CLM.numpft[]     = numpft
        CLM.nlevsclass[] = 1

        # ------------------------------------------------------------------------------
        # EDParams + EDPftvarcon (node thetas/resid + hydraulics) used by the solve.
        # ------------------------------------------------------------------------------
        edp = CLM.ed_params_type()
        edp.hydr_solver       = CLM.hydr_solver_1DTaylor
        edp.hydr_kmax_rsurf1  = 20.0
        edp.hydr_kmax_rsurf2  = 0.0001
        CLM.EDParams[] = edp

        pcon = CLM.EDPftvarcon_type()
        pcon.hydr_thetas_node = fill(0.75, numpft, CLM.n_plant_media)
        pcon.hydr_resid_node  = fill(0.15, numpft, CLM.n_plant_media)
        pcon.hydr_srl         = fill(25.0, numpft)
        pcon.hydr_rs2         = fill(0.0001, numpft)
        pcon.hydr_rfrac_stem  = fill(0.625, numpft)
        pcon.hydr_p_taper     = fill(0.333, numpft)
        pcon.hydr_kmax_node   = fill(1.0e-7, numpft, CLM.n_plant_media)
        CLM.EDPftvarcon_inst[] = pcon

        # ------------------------------------------------------------------------------
        # Plant WRF/WKF globals (TFS), mirroring InitHydroGlobals / the solver test.
        # ------------------------------------------------------------------------------
        function build_tfs_wrf(pm::Int)
            th_sat = 0.75; th_res = 0.15; pinot = -1.5; epsil = 12.0
            rwc_ft = (pm == CLM.leaf_p_media) ? 1.0 : 0.958
            if pm == CLM.leaf_p_media
                cap_slp = 0.0; cap_int = 0.0; cap_corr = 1.0
            else
                hydr_psi0 = 0.0; hydr_psicap = -0.6; rwccap_pm = 0.947
                cap_slp = (hydr_psi0 - hydr_psicap) / (1.0 - rwccap_pm)
                cap_int = -cap_slp + hydr_psi0
                cap_corr = -cap_int / cap_slp
            end
            wrf = CLM.wrf_type_tfs()
            CLM.set_wrf_param!(wrf, [th_sat, th_res, pinot, epsil, rwc_ft,
                                     cap_corr, cap_int, cap_slp, Float64(pm)])
            return wrf
        end
        wrfp = Matrix{Union{CLM.WRFType,Nothing}}(nothing, CLM.n_plant_media + 1, numpft)
        wkfp = Matrix{Union{CLM.WKFType,Nothing}}(nothing, CLM.n_plant_media + 1, numpft)
        for pm in 1:CLM.n_plant_media
            w = build_tfs_wrf(pm)
            wrfp[pm + 1, ft] = w
            wkf = CLM.wkf_type_tfs()
            CLM.set_wkf_param!(wkf, [-1.5, 2.0])
            wkf.wrf = w
            wkfp[pm + 1, ft] = wkf
        end
        sto_wkf = CLM.wkf_type_tfs()
        CLM.set_wkf_param!(sto_wkf, [-1.5, 2.0])
        sto_wkf.wrf = wrfp[CLM.leaf_p_media + 1, ft]
        wkfp[CLM.stomata_p_media + 1, ft] = sto_wkf
        CLM._wrf_plant[] = wrfp
        CLM._wkf_plant[] = wkfp

        # ------------------------------------------------------------------------------
        # Full single woody-PFT allometry param table for the geometry helpers.
        # ------------------------------------------------------------------------------
        p = CLM.prt_params
        CLM.allocate_prt_params!(p, numpft, CLM.num_organ_types, 1)
        p.c2b .= 2.0; p.wood_density .= 0.6; p.slatop .= 0.012; p.slamax .= 0.020
        p.allom_agb_frac .= 0.6; p.allom_dbh_maxheight .= 90.0
        p.allom_la_per_sa_int .= 0.8e3; p.allom_la_per_sa_slp .= 0.0
        p.allom_sai_scaler .= 0.1; p.allom_l2fr .= 1.0; p.cushion .= 1.0
        p.allom_blca_expnt_diff .= 0.0
        p.allom_d2ca_coefficient_min .= 0.3; p.allom_d2ca_coefficient_max .= 0.6
        p.allom_h2cd1 .= 0.5; p.allom_h2cd2 .= 1.0
        p.allom_dmode .= 1; p.woody .= 1; p.allom_cmode .= 1; p.allom_smode .= 1
        p.allom_fmode .= 1; p.allom_stmode .= 1; p.allom_hmode .= 1
        p.allom_amode .= 1; p.allom_lmode .= 1
        p.allom_d2h1 .= 0.64; p.allom_d2h2 .= 0.37; p.allom_d2h3 .= -0.034
        p.allom_agb1 .= 0.06896; p.allom_agb2 .= 0.572; p.allom_agb3 .= 1.94; p.allom_agb4 .= 0.931
        p.allom_d2bl1 .= 0.07; p.allom_d2bl2 .= 1.3; p.allom_d2bl3 .= 0.55
        p.fnrt_prof_mode .= 1; p.fnrt_prof_a .= 0.976; p.fnrt_prof_b .= 0.0
        p.allom_zroot_max_dbh .= 100.0; p.allom_zroot_min_dbh .= 1.0
        p.allom_zroot_max_z .= 0.9;     p.allom_zroot_min_z .= 0.9
        p.allom_zroot_k .= 10.0

        # param_derived (branch_frac used by bsap_allom).
        pd = CLM.param_derived_type()
        pd.branch_frac = fill(0.25, numpft)
        CLM.ParamDerived[] = pd

        # PRT global + carbon-only element registry for the PRT GetState path.
        CLM.InitPRTGlobalAllometricCarbon!()
        CLM.num_elements[] = 1
        empty!(CLM.element_list)
        push!(CLM.element_list, CLM.carbon12_element)
        CLM.hlm_parteh_mode[] = CLM.prt_carbon_allom_hyp

        # ------------------------------------------------------------------------------
        # A minimal soil bc_in: 3 soil layers (combine12 -> 2 rhiz layers).
        # ------------------------------------------------------------------------------
        nlevsoil = 3
        function build_bc_in()
            dz = [0.1, 0.3, 0.6]
            zi = [0.0, 0.1, 0.4, 1.0]   # zero index for the surface -> length nlevsoil+1
            CLM.bc_in_type(
                nlevsoil        = nlevsoil,
                zi_sisl         = zi,
                dz_sisl         = dz,
                eff_porosity_sl = fill(0.45, nlevsoil),
                watsat_sisl     = fill(0.45, nlevsoil),
                sucsat_sisl     = fill(200.0, nlevsoil),   # mm
                bsw_sisl        = fill(5.0, nlevsoil),
                hksat_sisl      = fill(0.01, nlevsoil),    # mm/s
                h2o_liq_sisl    = [0.1 * 0.4 * CLM.dens_fresh_liquid_water,
                                   0.3 * 0.4 * CLM.dens_fresh_liquid_water,
                                   0.6 * 0.4 * CLM.dens_fresh_liquid_water],  # kg/m2 (~theta 0.4)
                smpmin_si       = -2.0e5,                  # mm
                qflx_transp_pa  = [0.0],                   # patch transpiration (mm/s); 1 veg patch
            )
        end

        # ------------------------------------------------------------------------------
        # Build a site whose si_hydr is initialized + cold-started, with rhizosphere
        # geometry sized from a cohort's root length.
        # ------------------------------------------------------------------------------
        sites  = [CLM.ed_site_type()]
        bc_in  = [build_bc_in()]
        bc_out = [CLM.bc_out_type(
                    qflx_soil2root_sisl = zeros(Float64, nlevsoil),
                    qflx_ro_sisl        = zeros(Float64, nlevsoil),
                    btran_pa            = zeros(Float64, 1))]
        CLM.InitHydrSites(sites, bc_in)
        CLM.HydrSiteColdStart(sites, bc_in)
        site = sites[1]
        site.lat = 50.0; site.lon = -115.0
        csh = site.si_hydr

        # Build a real-PRT cohort with hydro geometry from UpdateSizeDepPlantHydProps.
        function seed_prt(dbh0::Float64)
            prt = CLM.callom_prt_vartypes()
            CLM.InitPRTVartype!(prt)
            ct = 1.0; cd = 1; ef = 1.0
            l2fr = p.allom_l2fr[ft]
            tgt_leaf, _    = CLM.bleaf(dbh0, ft, cd, ct, ef)
            tgt_fnrt, _    = CLM.bfineroot(dbh0, ft, ct, l2fr, ef)
            _, tgt_sapw, _ = CLM.bsap_allom(dbh0, ft, cd, ct, ef)
            tgt_store, _   = CLM.bstore_allom(dbh0, ft, cd, ct)
            tgt_agw, _     = CLM.bagw_allom(dbh0, ft, cd, ef)
            tgt_bgw, _     = CLM.bbgw_allom(dbh0, ft, ef)
            tgt_struct, _  = CLM.bdead_allom(tgt_agw, tgt_bgw, tgt_sapw, ft)
            CLM.SetState!(prt, CLM.leaf_organ,   CLM.carbon12_element, tgt_leaf)
            CLM.SetState!(prt, CLM.fnrt_organ,   CLM.carbon12_element, tgt_fnrt)
            CLM.SetState!(prt, CLM.sapw_organ,   CLM.carbon12_element, tgt_sapw)
            CLM.SetState!(prt, CLM.store_organ,  CLM.carbon12_element, tgt_store)
            CLM.SetState!(prt, CLM.struct_organ, CLM.carbon12_element, tgt_struct)
            CLM.SetState!(prt, CLM.repro_organ,  CLM.carbon12_element, 0.0)
            return prt
        end

        function real_cohort(; n::Float64, dbh::Float64, th::Float64,
                              g_sb::Float64, newly::Bool=false)
            c = CLM.fates_cohort_type()
            c.pft = ft; c.n = n; c.dbh = dbh
            c.crowndamage = 1; c.canopy_trim = 1.0
            c.efleaf_coh = 1.0; c.efstem_coh = 1.0
            c.height, _ = CLM.h_allom(dbh, ft)
            c.g_sb_laweight = g_sb
            c.size_class = 1
            c.treelai = 1.0
            c.prt = seed_prt(dbh)
            CLM.InitHydrCohort(site, c)
            CLM.UpdatePlantHydrNodes!(c, ft, c.height, csh)
            CLM.UpdateSizeDepPlantHydProps!(site, c, bc_in[1])
            CLM.SavePreviousCompartmentVolumes!(c.co_hydr)
            ch = c.co_hydr
            ch.th_ag .= th; ch.th_troot = th; ch.th_aroot .= th
            ch.is_newly_recruited = newly
            CLM.InitPlantHydStates!(site, c)
            return c
        end

        # One vegetated patch with a single cohort.
        cohort = real_cohort(n = 0.10, dbh = 10.0, th = 0.55, g_sb = 0.02)
        patch = CLM.fates_patch_type()
        patch.area = CLM.area
        patch.total_canopy_area = CLM.area
        patch.nocomp_pft_label = ft   # not bareground
        patch.tallest = cohort
        patch.shortest = cohort
        site.youngest_patch = patch
        site.oldest_patch   = patch

        # Build the site-level rhizosphere geometry from this cohort's root length.
        CLM.UpdateSizeDepRhizHydProps(site, bc_in[1])

        @test csh.nlevrhiz == max(1, nlevsoil - 1)
        @test all(csh.l_aroot_layer .> 0.0)
        @test all(isfinite, csh.v_shell)

        # --- Test: SetMaxCondConnections returns finite conductances ---
        @testset "SetMaxCondConnections finite" begin
            ch = cohort.co_hydr
            # SetMaxCondConnections spans all rhiz layers; its node index runs over the
            # full plant + per-layer (aroot+shell) continuum.
            ntot_nodes = CLM.n_hypool_ag + 1 +
                         (CLM.n_hypool_aroot + CLM.nshell) * csh.nlevrhiz
            ncnx = CLM.n_hypool_leaf + CLM.n_hypool_stem + CLM.n_hypool_troot - 1 +
                   (CLM.n_hypool_aroot + CLM.nshell) * csh.nlevrhiz
            h_node = zeros(Float64, ntot_nodes)
            kdn, kup = CLM.SetMaxCondConnections(csh, ch, h_node)
            @test length(kdn) == ncnx
            @test length(kup) == ncnx
            @test all(isfinite, kdn)
            @test all(isfinite, kup)
            @test all(kdn .> 0.0)
            @test all(kup .> 0.0)
        end

        # --- Test: FillDrainRhizShells conserves rhizosphere water ---
        @testset "FillDrainRhizShells conserves water" begin
            # Total rhizosphere water before (kg/m2): the host h2o_liq for the mapped
            # soil layers IS the conservation target (FillDrain redistributes shells to
            # match it). Run + assert per-layer closure (the routine itself throws if
            # |errh2o| > 1e-9; here we additionally check the post-fill shell water
            # matches the host soil water).
            CLM.FillDrainRhizShells(1, sites, bc_in, bc_out)
            for j in 1:csh.nlevrhiz
                j_t = csh.map_r2s[j, 1]; j_b = csh.map_r2s[j, 2]
                shellw = sum(csh.h2osoi_liqvol_shell[j, :] .* csh.v_shell[j, :]) *
                         CLM.dens_fresh_liquid_water * CLM.AREA_INV
                hostw = sum(bc_in[1].h2o_liq_sisl[j_t:j_b])
                @test isapprox(shellw, hostw; atol = 1.0e-8)
            end
        end

        # --- Test: RecruitWUptake finite ---
        @testset "RecruitWUptake finite" begin
            rsites = [CLM.ed_site_type()]
            rbc    = [build_bc_in()]
            CLM.InitHydrSites(rsites, rbc)
            CLM.HydrSiteColdStart(rsites, rbc)
            rsite = rsites[1]
            rcsh = rsite.si_hydr
            rsite.lat = 50.0; rsite.lon = -115.0

            # A newly-recruited cohort on a vegetated patch.
            saved_site = (; site)  # keep `site` global intact while we reuse helpers
            # Temporarily point module-level `site`/`csh` at the recruit site via closure
            # is awkward; instead build the recruit cohort inline with rsite geometry.
            rc = CLM.fates_cohort_type()
            rc.pft = ft; rc.n = 0.05; rc.dbh = 9.0
            rc.crowndamage = 1; rc.canopy_trim = 1.0
            rc.efleaf_coh = 1.0; rc.efstem_coh = 1.0
            rc.height, _ = CLM.h_allom(9.0, ft)
            rc.size_class = 1
            rc.prt = seed_prt(9.0)
            CLM.InitHydrCohort(rsite, rc)
            CLM.UpdatePlantHydrNodes!(rc, ft, rc.height, rcsh)
            CLM.UpdateSizeDepPlantHydProps!(rsite, rc, rbc[1])
            rch = rc.co_hydr
            rch.th_ag .= 0.5; rch.th_troot = 0.5; rch.th_aroot .= 0.5
            rch.is_newly_recruited = true

            rpatch = CLM.fates_patch_type()
            rpatch.area = CLM.area
            rpatch.total_canopy_area = CLM.area
            rpatch.nocomp_pft_label = ft
            rpatch.tallest = rc; rpatch.shortest = rc
            rsite.oldest_patch = rpatch; rsite.youngest_patch = rpatch
            CLM.UpdateSizeDepRhizHydProps(rsite, rbc[1])

            recruitflag = CLM.RecruitWUptake(1, rsites, rbc, 1800.0)
            @test recruitflag == true
            @test all(isfinite, rcsh.recruit_w_uptake[1:rcsh.nlevrhiz])
            @test sum(rcsh.recruit_w_uptake[1:rcsh.nlevrhiz]) > 0.0
            # Flag cleared after uptake.
            @test rch.is_newly_recruited == false
        end

        # --- Test: hydraulics_drive populates ftc_*/btran/leaf_psi + closes balance ---
        @testset "hydraulics_drive populates ftc_*/btran/leaf_psi" begin
            ch = cohort.co_hydr

            # Reset the site stored-water bookkeeping for a clean balance check, and set
            # a moderate transpiration demand on the patch.
            CLM.UpdateH2OVeg!(site, bc_out[1])
            bc_in[1].qflx_transp_pa[1] = 1.0e-6   # mm H2O/s into root

            dtime = 1800.0
            # This is the full dispatcher (FillDrainRhizShells + Hydraulics_BC). It must
            # run without throwing (the internal plant-water balance check is active).
            CLM.hydraulics_drive(1, sites, bc_in, bc_out, dtime)

            # ftc_* populated and in [0,1]
            @test all(isfinite, ch.ftc_ag)
            @test all(0.0 .<= ch.ftc_ag .<= 1.0)
            @test isfinite(ch.ftc_troot)
            @test 0.0 <= ch.ftc_troot <= 1.0
            @test all(isfinite, ch.ftc_aroot)
            @test all(0.0 .<= ch.ftc_aroot .<= 1.0)

            # btran populated and in [0,1]
            @test isfinite(ch.btran)
            @test 0.0 <= ch.btran <= 1.0

            # leaf_psi (psi_ag[1]) finite + non-positive (suction)
            @test isfinite(ch.psi_ag[1])
            @test ch.psi_ag[1] <= 0.0

            # Site mass balance bookkeeping is finite + small (Hydraulics_BC would have
            # thrown if the plant water balance did not close to error_thresh).
            @test isfinite(csh.errh2o_hyd)
            @test abs(csh.errh2o_hyd) < 1.0e-4
            @test isfinite(csh.h2oveg)
            @test isfinite(bc_out[1].plant_stored_h2o_si)

            # Root uptake disaggregated onto the host soil layers (finite).
            @test all(isfinite, bc_out[1].qflx_soil2root_sisl)
        end

    finally
        # Restore global state.
        CLM.hlm_use_planthydro[] = saved_planthydro
        CLM._wrf_plant[]         = saved_wrf
        CLM._wkf_plant[]         = saved_wkf
        CLM.EDPftvarcon_inst[]   = saved_pftvarcon
        CLM.EDParams[]           = saved_edparams
        CLM.numpft[]             = saved_numpft
        CLM.nlevsclass[]         = saved_nlevsclass
        CLM.prt_global[]         = saved_prtglobal
        CLM.num_elements[]       = saved_numel
        empty!(CLM.element_list); append!(CLM.element_list, saved_ellist)
        CLM.hlm_parteh_mode[]    = saved_parteh
        CLM.ParamDerived[]       = saved_paramderiv
    end
end
