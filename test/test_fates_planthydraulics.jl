# test_fates_planthydraulics.jl
# Tests for the FATES plant-hydraulics solver (FatesPlantHydraulicsMod.jl).
#
# We build a synthetic single-cohort / single-patch site with a small soil column and one
# rhizosphere shell, set up TFS water-transfer functions, initialize the hydraulic states,
# and run one implicit 1D-Taylor hydraulics solve with a prescribed transpiration demand.
#
# Asserts:
#   1. Water mass conservation across the network (plant storage change == soil uptake -
#      transpiration, to a tight tolerance).
#   2. A monotonic total-potential gradient soil -> leaf under transpiration.
#   3. btran in [0,1], and btran decreasing as the soil dries.
#
# We construct the cohort hydraulic geometry / kmax's directly rather than running the full
# carbon-pool-driven UpdatePlantHydr* path, so the test isolates the solver math without
# needing a fully-populated FATES parameter file. The solver, mass-balance bookkeeping,
# psi/ftc updates, and btran logic exercised here are the ported module's core.

using Test
using CLM

const FPH = CLM   # all FATES symbols live in the CLM module

@testset "FATES Plant Hydraulics" begin

    # ----------------------------------------------------------------------------------
    # Helpers to build TFS WRF/WKF functions for one PFT across the plant media.
    # ----------------------------------------------------------------------------------
    # TFS WRF params: [th_sat, th_res, pinot, epsil, rwc_ft, cap_corr, cap_int, cap_slp, pmedia]
    # TFS WKF params: [p50, avuln]
    function build_tfs_wrf(pm::Int)
        th_sat = 0.75
        th_res = 0.15
        pinot  = -1.5      # osmotic potential at full turgor [MPa]
        epsil  = 12.0      # bulk elastic modulus [MPa]
        rwc_ft = (pm == FPH.leaf_p_media) ? 1.0 : 0.958
        if pm == FPH.leaf_p_media
            cap_slp = 0.0; cap_int = 0.0; cap_corr = 1.0
        else
            # mirror InitHydroGlobals: cap from hydr_psi0 / hydr_psicap
            hydr_psi0   = 0.0
            hydr_psicap = -0.6
            rwccap_pm   = 0.947
            cap_slp = (hydr_psi0 - hydr_psicap) / (1.0 - rwccap_pm)
            cap_int = -cap_slp + hydr_psi0
            cap_corr = -cap_int / cap_slp
        end
        wrf = CLM.wrf_type_tfs()
        CLM.set_wrf_param!(wrf, [th_sat, th_res, pinot, epsil, rwc_ft,
                                 cap_corr, cap_int, cap_slp, Float64(pm)])
        return wrf
    end
    function build_tfs_wkf(wrf)
        wkf = CLM.wkf_type_tfs()
        CLM.set_wkf_param!(wkf, [-1.5, 2.0])   # p50 = -1.5 MPa, avuln = 2
        wkf.wrf = wrf
        return wkf
    end

    # Inject global plant WRF/WKF for 1 PFT (rows = pm+1, pm = 0..n_plant_media)
    numpft = 1
    ft = 1
    wrfp = Matrix{Union{CLM.WRFType,Nothing}}(nothing, FPH.n_plant_media + 1, numpft)
    wkfp = Matrix{Union{CLM.WKFType,Nothing}}(nothing, FPH.n_plant_media + 1, numpft)
    for pm in 1:FPH.n_plant_media
        w = build_tfs_wrf(pm)
        wrfp[pm + 1, ft] = w
        wkfp[pm + 1, ft] = build_tfs_wkf(w)
    end
    # stomata media (pm=0): conductance points at leaf retention function
    leaf_wrf = wrfp[FPH.leaf_p_media + 1, ft]
    sto_wkf = CLM.wkf_type_tfs()
    CLM.set_wkf_param!(sto_wkf, [-1.5, 2.0])
    sto_wkf.wrf = leaf_wrf
    wkfp[FPH.stomata_p_media + 1, ft] = sto_wkf
    CLM._wrf_plant[] = wrfp
    CLM._wkf_plant[] = wkfp

    # Soil WRF/WKF (Campbell), used by the rhizosphere shell.
    function build_soil_cch()
        wrf = CLM.wrf_type_cch()
        CLM.set_wrf_param!(wrf, [0.45, -0.0001, 4.0])   # [th_sat, psi_sat (MPa), beta]
        wkf = CLM.wkf_type_cch()
        CLM.set_wkf_param!(wkf, [0.45, -0.0001, 4.0])
        wkf.wrf = wrf
        return wrf, wkf
    end

    # ----------------------------------------------------------------------------------
    # Build a single-layer site hydraulics object + a cohort, with geometry & kmax set
    # directly. nshell == 1, nlevrhiz == 1, n_hypool_tot == 5.
    # ----------------------------------------------------------------------------------
    function build_site_and_cohort(; soil_theta::Float64)
        nlevrhiz = 1
        csite = CLM.ed_site_hydr_type()
        csite.nlevrhiz = nlevrhiz
        CLM.InitHydrSite!(csite, numpft, 1, CLM.hydr_solver_1DTaylor, nlevrhiz)

        # Rhizosphere geometry for the single layer
        dz = 1.0                 # layer width [m]
        zi = 1.0                 # bottom-edge depth [m] (positive down)
        csite.dz_rhiz[1] = dz
        csite.zi_rhiz[1] = zi
        csite.map_r2s[1, 1] = 1
        csite.map_r2s[1, 2] = 1
        csite.map_s2r[1] = 1

        # Soil functions
        swrf, swkf = build_soil_cch()
        csite.wrf_soil[1] = swrf
        csite.wkf_soil[1] = swkf

        # Absorbing-root length on the site (total over all plants) [m]
        l_aroot_site = 50.0
        csite.l_aroot_layer[1] = l_aroot_site
        csite.l_aroot_layer_init[1] = l_aroot_site

        # Rhizosphere shell geometry for this layer (single shell)
        rs1 = FPH.fine_root_radius_const
        CLM.shellGeom!(l_aroot_site, rs1, CLM.area, dz,
                       view(csite.r_out_shell, 1, :),
                       view(csite.r_node_shell, 1, :),
                       view(csite.v_shell, 1, :))
        csite.v_shell_init[1, 1] = csite.v_shell[1, 1]

        # Shell water content [m3/m3]
        csite.h2osoi_liqvol_shell[1, 1] = soil_theta

        # Max soil conductance to upper / lower shell boundary [kg s-1 MPa-1].
        # (Site-level path conductance the solver scales by aroot_frac_plant.)
        csite.kmax_upper_shell[1, 1] = 1.0e-3
        csite.kmax_lower_shell[1, 1] = 1.0e-3

        # ---- Cohort ----
        cohort = CLM.fates_cohort_type()
        cohort.pft = ft
        cohort.n = 0.1            # individuals / m2
        cohort.dbh = 10.0
        cohort.height = 8.0
        ch = CLM.ed_cohort_hydr_type()
        CLM.AllocateHydrCohortArrays!(ch, nlevrhiz)
        cohort.co_hydr = ch

        # Node heights (leaf high, stem mid, troot below ground)
        ch.z_node_ag[1] = 7.96     # leaf
        ch.z_lower_ag[1] = 7.92
        ch.z_upper_ag[1] = 8.0
        ch.z_node_ag[2] = 3.95     # stem
        ch.z_lower_ag[2] = 0.0
        ch.z_upper_ag[2] = 7.9
        ch.z_node_troot = -0.5     # transporting root (below surface)

        # Compartment volumes [m3 / indiv]
        ch.v_ag[1] = 2.0e-4        # leaf
        ch.v_ag[2] = 2.0e-3        # stem
        ch.v_troot = 1.0e-3
        ch.v_aroot_layer[1] = 1.0e-3

        # Absorbing-root length for this cohort in this layer [m / plant]
        ch.l_aroot_layer[1] = l_aroot_site   # one cohort holds all the site root

        # Maximum conductances [kg s-1 MPa-1] (directly set, plausible magnitudes)
        ch.kmax_petiole_to_leaf = 1.0e8
        ch.kmax_stem_upper[1] = 1.0e-3
        ch.kmax_stem_lower[1] = 1.0e-3
        ch.kmax_troot_upper = 1.0e-3
        ch.kmax_troot_lower[1] = 1.0e-3
        ch.kmax_aroot_upper[1] = 1.0e-3
        ch.kmax_aroot_lower[1] = 1.0e-3
        ch.kmax_aroot_radial_in[1]  = 1.0e-3
        ch.kmax_aroot_radial_out[1] = 1.0e-3

        return csite, cohort
    end

    # ----------------------------------------------------------------------------------
    # Run one solve and return diagnostics.
    # ----------------------------------------------------------------------------------
    function run_one_solve(; soil_theta::Float64, q_top::Float64, dtime::Float64=1800.0)
        csite, cohort = build_site_and_cohort(soil_theta=soil_theta)
        ch = cohort.co_hydr

        # Initialize plant hydraulic states from soil
        site = (; si_hydr = csite)   # InitPlantHydStates! only uses .si_hydr
        CLM.InitPlantHydStates!(site, cohort)

        # Total plant water before [kg/indiv]
        w_plant_beg = (sum(ch.th_ag .* ch.v_ag) + ch.th_troot * ch.v_troot +
                       ch.th_aroot[1] * ch.v_aroot_layer[1]) * CLM.dens_fresh_liquid_water

        # Solve
        ordered = collect(1:csite.nlevrhiz)
        kbg = zeros(Float64, csite.nlevrhiz)
        CLM.OrderLayersForSolve1D(csite, cohort, ch, ordered, kbg)

        dth_layershell = zeros(Float64, CLM.nlevsoi_hyd_max, CLM.nshell)
        sapflow, rootuptake, wb_err, dwat_plant =
            CLM.ImTaylorSolve1D(0.0, 0.0, false, csite, cohort, ch, dtime, q_top,
                                ordered, kbg, dth_layershell)

        # Update psi/ftc and btran from new theta
        CLM.UpdatePlantPsiFTCFromTheta!(cohort, csite)

        w_plant_end = (sum(ch.th_ag .* ch.v_ag) + ch.th_troot * ch.v_troot +
                       ch.th_aroot[1] * ch.v_aroot_layer[1]) * CLM.dens_fresh_liquid_water

        # Soil water removed from the rhizosphere over the step [kg/indiv]
        # dth_layershell stores change scaled by (l_aroot_layer * n / l_aroot_site);
        # here that scale == n, so undo it to per-plant.
        soil_uptake = -dth_layershell[1, 1] * csite.v_shell[1, 1] *
                      CLM.dens_fresh_liquid_water / cohort.n
        transp = q_top * dtime   # kg/indiv

        return (; ch, csite, sapflow, rootuptake, wb_err, dwat_plant,
                w_plant_beg, w_plant_end, soil_uptake, transp)
    end

    # ----------------------------------------------------------------------------------
    # Test 1: water mass conservation under a moderate transpiration demand.
    # ----------------------------------------------------------------------------------
    @testset "mass conservation" begin
        q_top = 1.0e-7   # kg/indiv/s
        r = run_one_solve(soil_theta=0.40, q_top=q_top)

        # The solver tracks its own per-step balance error; it must be tiny.
        @test abs(r.wb_err) < 1.0e-4

        # Change in plant storage == (soil uptake - transpiration), to tolerance.
        delta_storage = r.w_plant_end - r.w_plant_beg
        @test isapprox(delta_storage, r.soil_uptake - r.transp;
                       atol = 1.0e-4, rtol = 1.0e-3)

        # dwat_plant (returned by the solver) must match the measured storage change.
        @test isapprox(r.dwat_plant, delta_storage; atol = 1.0e-6)

        # Roots took up water (positive uptake under transpiration from moist soil).
        @test r.soil_uptake > 0.0
    end

    # ----------------------------------------------------------------------------------
    # Test 2: monotonic potential gradient soil -> leaf under transpiration.
    # ----------------------------------------------------------------------------------
    @testset "potential gradient" begin
        q_top = 2.0e-7
        r = run_one_solve(soil_theta=0.40, q_top=q_top)
        ch = r.ch
        csite = r.csite

        geo = CLM.mpa_per_pa * CLM.dens_fresh_liquid_water * CLM.grav_earth

        # Total potentials (matric + gravitational) along the path leaf<-stem<-troot<-aroot<-soil
        h_leaf  = ch.psi_ag[1]    + geo * ch.z_node_ag[1]
        h_stem  = ch.psi_ag[2]    + geo * ch.z_node_ag[2]
        h_troot = ch.psi_troot    + geo * ch.z_node_troot
        h_aroot = ch.psi_aroot[1] + geo * (-csite.zi_rhiz[1] + 0.5 * csite.dz_rhiz[1])
        psi_soil = CLM.psi_from_th(csite.wrf_soil[1], csite.h2osoi_liqvol_shell[1, 1])
        h_soil  = psi_soil + geo * (-csite.zi_rhiz[1] + 0.5 * csite.dz_rhiz[1])

        # Under transpiration, water flows soil -> leaf, so total potential must
        # decrease monotonically from soil to leaf.
        @test h_soil  >= h_aroot
        @test h_aroot >= h_troot
        @test h_troot >= h_stem
        @test h_stem  >= h_leaf
    end

    # ----------------------------------------------------------------------------------
    # Test 3: btran in [0,1] and decreasing under drier soil.
    # ----------------------------------------------------------------------------------
    @testset "btran range and drought response" begin
        q_top = 2.0e-7
        r_wet = run_one_solve(soil_theta=0.42, q_top=q_top)
        r_dry = run_one_solve(soil_theta=0.20, q_top=q_top)

        @test 0.0 <= r_wet.ch.btran <= 1.0
        @test 0.0 <= r_dry.ch.btran <= 1.0

        # Drier soil -> lower leaf potential -> lower btran (more stomatal stress).
        @test r_dry.ch.btran < r_wet.ch.btran
    end

end
