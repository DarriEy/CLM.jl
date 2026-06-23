# test_fates_planthydraulics_tierA.jl
# Tests for the FATES plant-hydraulics "Tier A" init/lifecycle routines added to
# FatesPlantHydraulicsMod.jl (InitHydrCohort / DeallocateHydrCohort /
# SavePreviousRhizVolumes / constrain_water_contents / AccumulateMortalityWaterStorage /
# RecruitWaterStorage / UpdateSizeDepRhizVolLenCon / UpdateSizeDepRhizHydProps /
# UpdateSizeDepPlantHydStates / FuseCohortHydraulics / InitHydrSites / HydrSiteColdStart).
#
# These routines are gated behind `hlm_use_planthydro==itrue`. We flip the flag on,
# construct hydro-enabled site/cohort objects directly (mirroring
# test_fates_planthydraulics.jl's "build geometry directly" style so the test stays
# self-contained, without a fully populated FATES parameter file), and assert:
#   * InitHydrCohort allocates the per-layer cohort hydro arrays.
#   * AccumulateMortalityWaterStorage accumulates a positive dead-vegetation water mass.
#   * UpdateSizeDepRhizVolLenCon / UpdateSizeDepRhizHydProps build finite, positive
#     rhizosphere volumes.
#   * FuseCohortHydraulics conserves total plant water mass across a fuse.
#   * HydrSiteColdStart builds a finite site hydro state.
# The `hlm_use_planthydro` flag is restored at the end.

using Test
using CLM

const FPHA = CLM   # all FATES symbols live in the CLM module

@testset "FATES Plant Hydraulics — Tier A (init/lifecycle)" begin

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
        # Minimal EDParams (hydr_solver) + EDPftvarcon (node thetas/resid) so the
        # state-clamp + factory routines have the params they read.
        # ------------------------------------------------------------------------------
        edp = CLM.ed_params_type()
        edp.hydr_solver = CLM.hydr_solver_1DTaylor
        CLM.EDParams[] = edp

        pcon = CLM.EDPftvarcon_type()
        # Indexed [ft, pm] with pm = leaf(1)/stem(2)/troot(3)/aroot(4).
        pcon.hydr_thetas_node = fill(0.75, numpft, CLM.n_plant_media)
        pcon.hydr_resid_node  = fill(0.15, numpft, CLM.n_plant_media)
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
        # A minimal soil bc_in: 3 soil layers (combine12 -> 2 rhiz layers).
        # ------------------------------------------------------------------------------
        nlevsoil = 3
        function build_bc_in()
            dz = [0.1, 0.3, 0.6]
            # zi carries a zero index for the surface -> length nlevsoil+1.
            zi = [0.0, 0.1, 0.4, 1.0]
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
            )
        end

        # ------------------------------------------------------------------------------
        # Build a real ed_site_type whose si_hydr is initialized through InitHydrSites.
        # ------------------------------------------------------------------------------
        sites  = [CLM.ed_site_type()]
        bc_in  = [build_bc_in()]
        CLM.InitHydrSites(sites, bc_in)
        site = sites[1]
        csh  = site.si_hydr

        # Build a cohort (on `site`) with directly-set hydraulic geometry/volumes.
        function build_cohort(; n::Float64, th::Float64)
            nlevrhiz = csh.nlevrhiz
            cohort = CLM.fates_cohort_type()
            cohort.pft    = ft
            cohort.n      = n
            cohort.dbh    = 10.0
            cohort.height = 8.0
            CLM.InitHydrCohort(site, cohort)   # allocate co_hydr
            ch = cohort.co_hydr
            ch.z_node_ag[1] = 7.96; ch.z_lower_ag[1] = 7.92; ch.z_upper_ag[1] = 8.0
            ch.z_node_ag[2] = 3.95; ch.z_lower_ag[2] = 0.0;  ch.z_upper_ag[2] = 7.9
            ch.z_node_troot = -0.5
            ch.v_ag[1] = 2.0e-4; ch.v_ag[2] = 2.0e-3
            ch.v_troot = 1.0e-3
            ch.v_ag_init    .= ch.v_ag
            ch.v_troot_init  = ch.v_troot
            for j in 1:nlevrhiz
                ch.v_aroot_layer[j]      = 1.0e-3
                ch.v_aroot_layer_init[j] = 1.0e-3
                ch.l_aroot_layer[j]      = 25.0
                ch.th_aroot[j]           = th
            end
            ch.th_ag .= th
            ch.th_troot = th
            return cohort
        end

        @test csh.nlevrhiz == max(1, nlevsoil - 1)   # combine12
        @test length(csh.l_aroot_layer) == csh.nlevrhiz
        @test all(isfinite, csh.dz_rhiz)
        @test all(isfinite, csh.zi_rhiz)

        # --- Test: InitHydrCohort allocates cohort hydro arrays ---
        @testset "InitHydrCohort allocates arrays" begin
            cohort = build_cohort(n = 0.1, th = 0.6)
            ch = cohort.co_hydr
            @test ch !== nothing
            @test length(ch.th_aroot)       == csh.nlevrhiz
            @test length(ch.v_aroot_layer)  == csh.nlevrhiz
            @test length(ch.kmax_aroot_upper) == csh.nlevrhiz
            @test ch.is_newly_recruited == false

            # DeallocateHydrCohort releases them.
            CLM.DeallocateHydrCohort(cohort)
            @test cohort.co_hydr === nothing
        end

        # --- Test: constrain_water_contents clamps into (thr+δ, ths-δ) ---
        @testset "constrain_water_contents clamps" begin
            δ = 1.0e-7
            @test CLM.constrain_water_contents(10.0, δ, ft, CLM.leaf_p_media) ≈ 0.75 - δ
            @test CLM.constrain_water_contents(-10.0, δ, ft, CLM.leaf_p_media) ≈ 0.15 + δ
            @test CLM.constrain_water_contents(0.5, δ, ft, CLM.leaf_p_media) ≈ 0.5
        end

        # --- Build a couple cohorts on a patch hung off the site. ---
        function attach_patch_with_cohorts(c1, c2)
            patch = CLM.fates_patch_type()
            patch.area = CLM.area
            patch.tallest = c1
            c1.shorter = c2
            patch.shortest = c2
            c2.taller = c1
            site.youngest_patch = patch
            site.oldest_patch   = patch
            return patch
        end

        # --- Test: AccumulateMortalityWaterStorage accumulates positive water ---
        @testset "AccumulateMortalityWaterStorage" begin
            csh.h2oveg      = 100.0
            csh.h2oveg_dead = 0.0
            cohort = build_cohort(n = 0.2, th = 0.6)
            CLM.AccumulateMortalityWaterStorage(site, cohort, cohort.n)
            @test csh.h2oveg_dead > 0.0
            @test isfinite(csh.h2oveg_dead)
            # Live pool drops by the same amount.
            @test csh.h2oveg ≈ 100.0 - csh.h2oveg_dead
        end

        # --- Test: RecruitWaterStorage diagnoses positive recruit water ---
        @testset "RecruitWaterStorage" begin
            c1 = build_cohort(n = 0.1, th = 0.6)
            c2 = build_cohort(n = 0.1, th = 0.6)
            c1.co_hydr.is_newly_recruited = true
            c2.co_hydr.is_newly_recruited = false
            attach_patch_with_cohorts(c1, c2)
            CLM.RecruitWaterStorage(1, sites, [nothing])
            @test csh.h2oveg_recruit > 0.0
            @test isfinite(csh.h2oveg_recruit)
        end

        # --- Test: UpdateSizeDepRhizVolLenCon / UpdateSizeDepRhizHydProps ---
        @testset "UpdateSizeDepRhiz* finite positive volumes" begin
            c1 = build_cohort(n = 0.1, th = 0.6)
            c2 = build_cohort(n = 0.1, th = 0.6)
            attach_patch_with_cohorts(c1, c2)

            CLM.UpdateSizeDepRhizHydProps(site, bc_in[1])

            @test all(csh.l_aroot_layer .> 0.0)
            @test all(isfinite, csh.l_aroot_layer)
            @test all(csh.v_shell .> 0.0)
            @test all(isfinite, csh.v_shell)
            @test all(csh.r_node_shell .> 0.0)
            @test all(csh.r_out_shell  .> 0.0)
            @test all(isfinite, csh.kmax_upper_shell)
            @test all(isfinite, csh.kmax_lower_shell)
            @test all(csh.kmax_upper_shell .> 0.0)
            @test all(csh.kmax_lower_shell .> 0.0)

            # SavePreviousRhizVolumes snapshots the current geometry.
            CLM.SavePreviousRhizVolumes(site)
            @test csh.v_shell_init      == csh.v_shell
            @test csh.r_node_shell_init == csh.r_node_shell
            @test csh.l_aroot_layer_init == csh.l_aroot_layer
        end

        # --- Test: UpdateSizeDepPlantHydStates rescales + clamps water contents ---
        @testset "UpdateSizeDepPlantHydStates" begin
            cohort = build_cohort(n = 0.1, th = 0.6)
            ch = cohort.co_hydr
            # Grow the leaf volume so init/new differ; conservation rescales theta.
            ch.v_ag[1]  = 2.0 * ch.v_ag_init[1]
            csh.h2oveg_growturn_err = 0.0
            CLM.UpdateSizeDepPlantHydStates(site, cohort)
            @test all(isfinite, ch.th_ag)
            @test all(0.15 .<= ch.th_ag .<= 0.75)
            @test isfinite(ch.th_troot)
            @test all(isfinite, ch.th_aroot)
            @test isfinite(csh.h2oveg_growturn_err)
        end

        # --- Test: FuseCohortHydraulics conserves total plant water mass ---
        # FuseCohortHydraulics re-derives the fused cohort's compartment geometry from
        # the PFT allometry/hydraulics params (it calls UpdatePlantHydr{Nodes,LenVol}),
        # so this block installs a full single-PFT param table + real PRT cohorts.
        @testset "FuseCohortHydraulics conserves water" begin
            # Full single woody-PFT param table for the allometry/hydraulics helpers.
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
            # Dynamic rooting-depth (MaximumRootingDepth).
            p.allom_zroot_max_dbh .= 100.0; p.allom_zroot_min_dbh .= 1.0
            p.allom_zroot_max_z .= 0.9;     p.allom_zroot_min_z .= 0.9
            p.allom_zroot_k .= 10.0

            # Hydraulics PFT params (UpdatePlantHydrLenVol / UpdatePlantKmax).
            pc = CLM.EDPftvarcon_inst[]
            pc.hydr_srl        = fill(25.0, numpft)
            pc.hydr_rs2        = fill(0.0001, numpft)
            pc.hydr_rfrac_stem = fill(0.625, numpft)
            pc.hydr_p_taper    = fill(0.333, numpft)
            pc.hydr_kmax_node  = fill(1.0e-7, numpft, CLM.n_plant_media)
            CLM.EDPftvarcon_inst[] = pc
            edp2 = CLM.EDParams[]
            edp2.hydr_kmax_rsurf1 = 20.0
            edp2.hydr_kmax_rsurf2 = 0.0001
            CLM.EDParams[] = edp2

            # param_derived: branch_frac (used by bsap_allom).
            pd = CLM.param_derived_type()
            pd.branch_frac = fill(0.25, numpft)
            CLM.ParamDerived[] = pd

            # PRT global + carbon-only element registry for the PRT GetState path.
            CLM.InitPRTGlobalAllometricCarbon!()
            CLM.num_elements[] = 1
            empty!(CLM.element_list)
            push!(CLM.element_list, CLM.carbon12_element)
            CLM.hlm_parteh_mode[] = CLM.prt_carbon_allom_hyp

            # Helper: a real carbon PRT object for one cohort (woody PFT1, dbh0).
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

            # Build a real-PRT cohort with hydro geometry from UpdateSizeDepPlantHydProps.
            function real_cohort(; n::Float64, dbh::Float64, th::Float64)
                c = CLM.fates_cohort_type()
                c.pft = ft; c.n = n; c.dbh = dbh
                c.crowndamage = 1; c.canopy_trim = 1.0
                c.efleaf_coh = 1.0; c.efstem_coh = 1.0
                c.height, _ = CLM.h_allom(dbh, ft)
                c.prt = seed_prt(dbh)
                CLM.InitHydrCohort(site, c)
                CLM.UpdatePlantHydrNodes!(c, ft, c.height, csh)
                CLM.UpdateSizeDepPlantHydProps!(site, c, bc_in[1])
                CLM.SavePreviousCompartmentVolumes!(c.co_hydr)
                ch = c.co_hydr
                ch.th_ag .= th; ch.th_troot = th; ch.th_aroot .= th
                return c
            end

            cur = real_cohort(n = 0.10, dbh = 10.0, th = 0.62)
            nxt = real_cohort(n = 0.05, dbh = 9.5,  th = 0.55)
            cur_h = cur.co_hydr; nxt_h = nxt.co_hydr

            # Pre-fusion total plant water [kg]. The current cohort contributes via its
            # *_init volumes (its geometry is re-derived inside the fuse), the donor via v.
            w_before =
                (sum(cur_h.th_ag .* cur_h.v_ag_init) +
                 cur_h.th_troot * cur_h.v_troot_init +
                 sum(cur_h.th_aroot .* cur_h.v_aroot_layer_init)) *
                 CLM.dens_fresh_liquid_water * cur.n +
                (sum(nxt_h.th_ag .* nxt_h.v_ag) +
                 nxt_h.th_troot * nxt_h.v_troot +
                 sum(nxt_h.th_aroot .* nxt_h.v_aroot_layer)) *
                 CLM.dens_fresh_liquid_water * nxt.n

            newn = cur.n + nxt.n
            CLM.FuseCohortHydraulics(site, cur, nxt, bc_in[1], newn)
            cur.n = newn

            # Post-fusion total water in the fused cohort (uses the re-derived volumes).
            w_after = (sum(cur_h.th_ag .* cur_h.v_ag) +
                       cur_h.th_troot * cur_h.v_troot +
                       sum(cur_h.th_aroot .* cur_h.v_aroot_layer)) *
                      CLM.dens_fresh_liquid_water * cur.n

            @test w_after > 0.0
            @test isapprox(w_after, w_before; rtol = 1.0e-9)
            @test all(isfinite, cur_h.th_ag)
            @test 0.0 <= cur_h.btran <= 1.0
            @test all(isfinite, cur_h.psi_ag)
            @test all(isfinite, cur_h.ftc_aroot)
        end

        # --- Test: HydrSiteColdStart builds a finite site hydro state ---
        @testset "HydrSiteColdStart finite state" begin
            cs_sites = [CLM.ed_site_type()]
            cs_bc    = [build_bc_in()]
            CLM.InitHydrSites(cs_sites, cs_bc)
            CLM.HydrSiteColdStart(cs_sites, cs_bc)
            csh2 = cs_sites[1].si_hydr

            @test all(isfinite, csh2.h2osoi_liqvol_shell)
            @test all(csh2.h2osoi_liqvol_shell .> 0.0)
            @test all(csh2.l_aroot_layer .== 0.0)
            @test length(csh2.wrf_soil) == csh2.nlevrhiz
            @test length(csh2.wkf_soil) == csh2.nlevrhiz
            for j in 1:csh2.nlevrhiz
                @test csh2.wrf_soil[j] isa CLM.WRFType
                @test csh2.wkf_soil[j] isa CLM.WKFType
                @test csh2.wkf_soil[j].wrf === csh2.wrf_soil[j]
                # psi from the cold-start shell water should be finite.
                psi = CLM.psi_from_th(csh2.wrf_soil[j], csh2.h2osoi_liqvol_shell[j, 1])
                @test isfinite(psi)
            end
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
