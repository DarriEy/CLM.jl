# test_fates_planthydro_2d_solvers.jl
# The FATES plant-hydraulics 2-D solvers (FatesPlantHydraulicsMod.F90 MatSolve2D /
# PicardSolve2D, fates_hydro_solver = 3 / 2) and their helpers:
#
#   1. `_dgesv!` (the Fortran DGESV): LU with partial pivoting against Julia's `\`,
#      a row-swap-only system, the LAPACK `info` for a singular matrix, and a type-
#      generic (non-Float64-storage) call.
#   2. The 2-D connection tables (SetConnections!) and SetMaxCondConnections! (the
#      in-place form equals the allocating one).
#   3. MatSolve2D on a small synthetic site (1 cohort, 2 rhizosphere layers): the
#      converged state satisfies the implicit mass balance it solves (residual
#      recomputed independently from the node arrays), the water balance closes, and
#      the shell change is the cohort's share.
#   4. PicardSolve2D on the same site: conservation, finite states, CFL substeps.
#   5. Both through hydraulics_drive (Hydraulics_BC dispatch on EDParams.hydr_solver),
#      plus the clm_fates_init! / clm_initialize! solver-option validation.
#   6. Data-gated parity vs the Fortran fates_parity_1pt runs (see README §"PLANT
#      HYDRAULICS: 2-D solvers"): Picard (solver 2) per-step vs the committed extract,
#      and Newton (solver 3) reproducing the Fortran's fatal non-convergence on the
#      first call, iteration by iteration against the committed Fortran trace.
using Test
using CLM
using LinearAlgebra

include("testdata.jl")

# Per-field max relative error tolerances of the Picard (solver 2) run vs the Fortran
# extract over nsteps 2-24 (day 1, before the first daily FATES step). Measured maxima
# (README §"PLANT HYDRAULICS: 2-D solvers"): PSI_LEAF 2.1e-5, PSI_LEAF_pft1 2.2e-5,
# PSI_LEAF_pft12 7.2e-5, TH_LEAF 1.3e-5, BTRAN_CO 5.0e-5, FATES_VEGH2O 1.1e-5; ~10x headroom.
const PICARD_TOL = Dict("PSI_LEAF" => 3e-4, "PSI_LEAF_pft1" => 3e-4, "PSI_LEAF_pft12" => 1e-3,
                        "TH_LEAF" => 2e-4, "BTRAN_CO" => 5e-4, "FATES_VEGH2O" => 2e-4)

# ---------------------------------------------------------------------------------
# Synthetic single-cohort site (same construction as test_fates_planthydraulics_tierB)
# ---------------------------------------------------------------------------------
function _ph2d_setup!(solver::Int)
    numpft = 1; ft = 1
    CLM.hlm_use_planthydro[] = CLM.itrue
    CLM.numpft[] = numpft
    CLM.nlevsclass[] = 1

    edp = CLM.ed_params_type()
    edp.hydr_solver      = solver
    edp.hydr_kmax_rsurf1 = 20.0
    edp.hydr_kmax_rsurf2 = 0.0001
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

    function build_tfs_wrf(pm::Int)
        th_sat = 0.75; th_res = 0.15; pinot = -1.5; epsil = 12.0
        rwc_ft = (pm == CLM.leaf_p_media) ? 1.0 : 0.958
        if pm == CLM.leaf_p_media
            cap_slp = 0.0; cap_int = 0.0; cap_corr = 1.0
        else
            cap_slp = (0.0 - -0.6) / (1.0 - 0.947)
            cap_int = -cap_slp
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
    sto = CLM.wkf_type_tfs()
    CLM.set_wkf_param!(sto, [-1.5, 2.0])
    sto.wrf = wrfp[CLM.leaf_p_media + 1, ft]
    wkfp[CLM.stomata_p_media + 1, ft] = sto
    CLM._wrf_plant[] = wrfp
    CLM._wkf_plant[] = wkfp

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
    p.allom_zroot_max_z .= 0.9; p.allom_zroot_min_z .= 0.9
    p.allom_zroot_k .= 10.0
    pd = CLM.param_derived_type()
    pd.branch_frac = fill(0.25, numpft)
    CLM.ParamDerived[] = pd
    CLM.InitPRTGlobalAllometricCarbon!()
    CLM.num_elements[] = 1
    empty!(CLM.element_list); push!(CLM.element_list, CLM.carbon12_element)
    CLM.hlm_parteh_mode[] = CLM.prt_carbon_allom_hyp

    nlevsoil = 3
    rho = CLM.dens_fresh_liquid_water
    bc_in = [CLM.bc_in_type(nlevsoil = nlevsoil, zi_sisl = [0.0, 0.1, 0.4, 1.0],
        dz_sisl = [0.1, 0.3, 0.6], eff_porosity_sl = fill(0.45, nlevsoil),
        watsat_sisl = fill(0.45, nlevsoil), sucsat_sisl = fill(200.0, nlevsoil),
        bsw_sisl = fill(5.0, nlevsoil), hksat_sisl = fill(0.01, nlevsoil),
        h2o_liq_sisl = [0.1 * 0.4 * rho, 0.3 * 0.35 * rho, 0.6 * 0.3 * rho],
        smpmin_si = -2.0e5, qflx_transp_pa = [0.0])]
    bc_out = [CLM.bc_out_type(qflx_soil2root_sisl = zeros(nlevsoil),
                              qflx_ro_sisl = zeros(nlevsoil), btran_pa = zeros(1))]
    sites = [CLM.ed_site_type()]
    CLM.InitHydrSites(sites, bc_in)
    CLM.HydrSiteColdStart(sites, bc_in)
    site = sites[1]; site.lat = 50.0; site.lon = -115.0
    csh = site.si_hydr

    dbh0 = 10.0
    prt = CLM.callom_prt_vartypes(); CLM.InitPRTVartype!(prt)
    tl, _ = CLM.bleaf(dbh0, ft, 1, 1.0, 1.0)
    tf, _ = CLM.bfineroot(dbh0, ft, 1.0, p.allom_l2fr[ft], 1.0)
    _, tsw, _ = CLM.bsap_allom(dbh0, ft, 1, 1.0, 1.0)
    tst, _ = CLM.bstore_allom(dbh0, ft, 1, 1.0)
    tag, _ = CLM.bagw_allom(dbh0, ft, 1, 1.0)
    tbg, _ = CLM.bbgw_allom(dbh0, ft, 1.0)
    tsd, _ = CLM.bdead_allom(tag, tbg, tsw, ft)
    for (o, v) in ((CLM.leaf_organ, tl), (CLM.fnrt_organ, tf), (CLM.sapw_organ, tsw),
                   (CLM.store_organ, tst), (CLM.struct_organ, tsd), (CLM.repro_organ, 0.0))
        CLM.SetState!(prt, o, CLM.carbon12_element, v)
    end
    c = CLM.fates_cohort_type()
    c.pft = ft; c.n = 0.10; c.dbh = dbh0; c.crowndamage = 1; c.canopy_trim = 1.0
    c.efleaf_coh = 1.0; c.efstem_coh = 1.0
    c.height, _ = CLM.h_allom(dbh0, ft)
    c.g_sb_laweight = 0.02; c.size_class = 1; c.treelai = 1.0; c.prt = prt
    CLM.InitHydrCohort(site, c)
    CLM.UpdatePlantHydrNodes!(c, ft, c.height, csh)
    CLM.UpdateSizeDepPlantHydProps!(site, c, bc_in[1])
    CLM.SavePreviousCompartmentVolumes!(c.co_hydr)
    ch = c.co_hydr
    ch.th_ag .= 0.55; ch.th_troot = 0.55; ch.th_aroot .= 0.55
    CLM.InitPlantHydStates!(site, c)

    patch = CLM.fates_patch_type()
    patch.area = CLM.area; patch.total_canopy_area = CLM.area
    patch.nocomp_pft_label = ft
    patch.tallest = c; patch.shortest = c
    site.youngest_patch = patch; site.oldest_patch = patch
    CLM.UpdateSizeDepRhizHydProps(site, bc_in[1])
    CLM.FillDrainRhizShells(1, sites, bc_in, bc_out)
    return (; sites, site, csh, cohort = c, ch, bc_in, bc_out)
end

"Mass-balance residual of the 2-D network at the solver's final node state (one substep)."
function _ph2d_residual(csh, ch, ft, dtime, qtop)
    nn = csh.num_nodes
    rho = CLM.dens_fresh_liquid_water
    res = zeros(nn); h = zeros(nn); ftc = zeros(nn); dftc = zeros(nn)
    for k in 1:nn
        wrf, wkf = CLM._node_wtf(csh, k, ft)
        res[k] = rho * csh.v_node[k] * (csh.th_node[k] - csh.th_node_init[k]) / dtime
        psi = CLM.psi_from_th(wrf, csh.th_node[k])
        h[k] = CLM.mpa_per_pa * rho * CLM.grav_earth * csh.z_node[k] + psi
        ftc[k] = CLM.ftc_from_psi(wkf, psi)
        dftc[k] = CLM.dftcdpsi_from_psi(wkf, psi)
    end
    kdn, kup = CLM.SetMaxCondConnections(csh, ch, h)
    for icnx in 1:csh.num_connections
        dn = csh.conn_dn[icnx]; up = csh.conn_up[icnx]
        _, _, keff = CLM.GetKAndDKDPsi(kdn[icnx], kup[icnx], h[dn], h[up], ftc[dn], ftc[up],
                                       dftc[dn], dftc[up])
        q = keff * (h[up] - h[dn])
        res[dn] -= q; res[up] += q
    end
    res[1] += qtop
    return res
end

@testset "FATES plant hydraulics: 2-D solvers (MatSolve2D / PicardSolve2D)" begin

    @testset "_dgesv! (DGESV)" begin
        for n in (1, 2, 5, 17, 40)
            A = randn(n, n) + n * I
            b = randn(n)
            a = copy(A); x = copy(b); ip = zeros(Int, n)
            @test CLM._dgesv!(a, x, ip) == 0
            @test x ≈ A \ b rtol = 1e-10
            @test all(1 .<= ip .<= n)
        end
        # Zero leading pivot: only solvable with a row interchange.
        A = [0.0 2.0 1.0; 1.0 1.0 0.0; 3.0 0.0 4.0]
        b = [3.0, 2.0, 7.0]
        a = copy(A); x = copy(b); ip = zeros(Int, 3)
        @test CLM._dgesv!(a, x, ip) == 0
        @test ip[1] == 3                     # IDAMAX picks |3| in column 1
        @test A * x ≈ b
        # Exactly singular: LAPACK info = first zero pivot of U; b untouched.
        S = [1.0 2.0 3.0; 2.0 4.0 6.0; 1.0 0.0 1.0]
        a = copy(S); x = [1.0, 2.0, 3.0]; ip = zeros(Int, 3)
        info = CLM._dgesv!(a, x, ip)
        @test info == 3
        @test x == [1.0, 2.0, 3.0]
        # Tridiagonal (the 1-layer network shape) agrees with the 1-D Thomas solver.
        n = 6
        lo = rand(n) .+ 0.1; di = rand(n) .+ 3.0; up = rand(n) .+ 0.1; r = randn(n)
        T = zeros(n, n)
        for i in 1:n
            T[i, i] = di[i]
            i > 1 && (T[i, i-1] = lo[i])
            i < n && (T[i, i+1] = up[i])
        end
        u, _ = CLM.Hydraulics_Tridiagonal(lo, di, up, r)
        a = copy(T); x = copy(r)
        @test CLM._dgesv!(a, x, zeros(Int, n)) == 0
        @test x ≈ u rtol = 1e-12
    end

    saved = (CLM.hlm_use_planthydro[], CLM._wrf_plant[], CLM._wkf_plant[],
             CLM.EDPftvarcon_inst[], CLM.EDParams[], CLM.numpft[], CLM.nlevsclass[],
             CLM.prt_global[], CLM.num_elements[], copy(CLM.element_list),
             CLM.hlm_parteh_mode[], CLM.ParamDerived[])
    try
        @testset "connection tables" begin
            s = _ph2d_setup!(CLM.hydr_solver_2DNewton)
            csh = s.csh
            nl = csh.nlevrhiz
            @test nl == 2
            @test csh.num_nodes == CLM.n_hypool_ag + CLM.n_hypool_troot +
                                   (CLM.n_hypool_aroot + CLM.nshell) * nl
            @test csh.num_connections == csh.num_nodes - 1      # a tree
            @test size(csh.ajac) == (csh.num_nodes, csh.num_nodes)
            # Every node is reached by exactly one connection from the leaf side, and
            # the troot (node n_hypool_ag+1) is the junction to every layer's aroot.
            itr = CLM.n_hypool_ag + 1
            @test count(==(itr), csh.conn_dn) == nl
            @test sort(csh.conn_up) == collect(2:csh.num_nodes)
            @test csh.pm_node[itr] == CLM.troot_p_media
            for j in 1:nl
                ia = itr + (j - 1) * (CLM.n_hypool_aroot + CLM.nshell) + 1
                @test csh.pm_node[ia] == CLM.aroot_p_media
                @test csh.pm_node[ia + 1] == CLM.rhiz_p_media
                @test csh.node_layer[ia + 1] == j
                # the aroot <-> 1st shell connection index used for rootuptake
                icnx_ar = CLM.n_hypool_ag + (j - 1) * (CLM.nshell + 1) + 2
                @test csh.conn_dn[icnx_ar] == ia && csh.conn_up[icnx_ar] == ia + 1
            end
            h = collect(range(-2.0, -0.1; length = csh.num_nodes))
            kdn, kup = CLM.SetMaxCondConnections(csh, s.ch, h)
            kdn2 = similar(kdn); kup2 = similar(kup)
            CLM.SetMaxCondConnections!(csh, s.ch, h, kdn2, kup2)
            @test kdn2 == kdn && kup2 == kup
            @test all(isfinite, kdn) && all(>(0), kup)
        end

        @testset "MatSolve2D (Newton) on a 2-layer network" begin
            s = _ph2d_setup!(CLM.hydr_solver_2DNewton)
            csh, ch, c = s.csh, s.ch, s.cohort
            dtime = 1800.0
            qtop = 2.0e-6                        # kg/s/plant
            dth = zeros(CLM.nlevsoi_hyd_max, CLM.nshell)
            th0 = (copy(ch.th_ag), ch.th_troot, copy(ch.th_aroot))
            shell0 = copy(csh.h2osoi_liqvol_shell)
            sapflow, ru, wb_err, dwat = CLM.MatSolve2D(csh, c, ch, dtime, qtop, dth)
            @test all(isfinite, (sapflow, wb_err, dwat)) && all(isfinite, ru)
            @test ch.iterh1 == 0                  # one substep, no failed round
            @test ch.iterh2 >= 1
            # The converged state satisfies the implicit balance it solved.
            res = _ph2d_residual(csh, ch, c.pft, dtime, qtop)
            @test maximum(abs, res) < 1.0e-7
            # Plant + rhizosphere water balance (kg/plant): flux out = change in storage.
            # (Newton stops at a signed-sum residual < 1e-8 kg/s, so <= 1e-8 kg/s x dt)
            @test abs(wb_err) < 1.0e-8 * dtime
            # Sapflow / uptake are the converged connection fluxes x dt; the troot
            # balance ties them together: sapflow = sum(uptake) - d(troot+aroot storage).
            rho = CLM.dens_fresh_liquid_water
            dstore_bg = (ch.th_troot - th0[2]) * ch.v_troot +
                        sum((ch.th_aroot .- th0[3]) .* ch.v_aroot_layer)
            @test sapflow ≈ sum(ru) - dstore_bg * rho atol = 2.0e-8 * dtime
            # Cohort water moved; the shell change is this cohort's (n-weighted) share.
            @test ch.th_ag[1] < th0[1][1]
            nl = csh.nlevrhiz
            for j in 1:nl
                inode = CLM.n_hypool_ag + CLM.n_hypool_troot + (j - 1) * (1 + CLM.nshell) + 2
                @test dth[j, 1] ≈ csh.dth_node[inode] * ch.l_aroot_layer[j] * c.n /
                                  csh.l_aroot_layer[j]
            end
            @test csh.h2osoi_liqvol_shell == shell0      # applied later by Hydraulics_BC
        end

        @testset "PicardSolve2D on a 2-layer network" begin
            s = _ph2d_setup!(CLM.hydr_solver_2DPicard)
            csh, ch, c = s.csh, s.ch, s.cohort
            dtime = 1800.0
            qtop = 2.0e-6
            dth = zeros(CLM.nlevsoi_hyd_max, CLM.nshell)
            th_leaf0 = ch.th_ag[1]
            sapflow, ru, wb_err, dwat = CLM.PicardSolve2D(csh, c, ch, dtime, qtop, dth,
                                                         csh.num_nodes)
            @test all(isfinite, (sapflow, wb_err, dwat)) && all(isfinite, ru)
            # Flux-form linear system: the column sums of the Jacobian cancel, so mass is
            # conserved up to the LU round-off (no theta clip at this demand).
            @test abs(wb_err) < 1.0e-6 * qtop * dtime
            @test ch.th_ag[1] < th_leaf0
            @test all(ch.th_ag .> 0)
            @test all(isfinite, csh.th_node)
            # A demand large enough for CFL > 1 splits the step (and still conserves).
            s2 = _ph2d_setup!(CLM.hydr_solver_2DPicard)
            ntrace = Ref(0); dts = Float64[]
            CLM._HYDRO2D_TRACE[] = (kind, args...) -> (ntrace[] += 1; push!(dts, args[5]))
            try
                qbig = 5.0 * s2.ch.v_ag[1] * CLM.dens_fresh_liquid_water / dtime
                _, _, wb2, _ = CLM.PicardSolve2D(s2.csh, s2.cohort, s2.ch, dtime, qbig,
                                                 zeros(CLM.nlevsoi_hyd_max, CLM.nshell),
                                                 s2.csh.num_nodes)
                @test isfinite(wb2)
            catch e
                # The leaf can empty under this demand; the Fortran negative-water trap.
                @test occursin("negative water content", sprint(showerror, e))
            finally
                CLM._HYDRO2D_TRACE[] = nothing
            end
            @test ntrace[] >= 2
            @test minimum(dts) < dtime
        end

        @testset "hydraulics_drive dispatch (solver $solver)" for solver in
                (CLM.hydr_solver_1DTaylor, CLM.hydr_solver_2DPicard, CLM.hydr_solver_2DNewton)
            s = _ph2d_setup!(solver)
            CLM.UpdateH2OVeg!(s.site, s.bc_out[1])
            # 2e-11 mm/s x canopy area / n = 2e-6 kg/s per plant (1e-6 mm/s would be
            # 0.1 kg/s per plant: it empties the leaf, Newton then fails to converge
            # exactly as the Fortran does, and Picard clips at theta_min).
            s.bc_in[1].qflx_transp_pa[1] = 2.0e-11
            CLM.hydraulics_drive(1, s.sites, s.bc_in, s.bc_out, 1800.0)
            ch = s.ch
            @test all(0.0 .<= ch.ftc_ag .<= 1.0) && 0.0 <= ch.btran <= 1.0
            @test isfinite(ch.psi_ag[1]) && ch.psi_ag[1] <= 0.0
            @test abs(s.csh.errh2o_hyd) < 1.0e-6
            @test all(isfinite, s.bc_out[1].qflx_soil2root_sisl)
            if solver != CLM.hydr_solver_1DTaylor
                @test all(isfinite, s.csh.q_flux)
            end
        end
    finally
        (CLM.hlm_use_planthydro[], CLM._wrf_plant[], CLM._wkf_plant[],
         CLM.EDPftvarcon_inst[], CLM.EDParams[], CLM.numpft[], CLM.nlevsclass[],
         CLM.prt_global[], CLM.num_elements[], _el, CLM.hlm_parteh_mode[],
         CLM.ParamDerived[]) = saved
        empty!(CLM.element_list); append!(CLM.element_list, _el)
    end

    @testset "solver option validation" begin
        old = (CLM.varctl.use_fates, CLM.varctl.use_fates_planthydro)
        try
            @test_throws "requires use_fates_planthydro" CLM.clm_initialize!(;
                fsurdat = "none.nc",
                paramfile = "none.nc", use_fates = true, use_fates_planthydro = false,
                fates_hydro_solver = 2)
            @test_throws "must be 1 (1-D Taylor)" CLM.clm_initialize!(; fsurdat = "none.nc",
                paramfile = "none.nc", use_fates = true, use_fates_planthydro = true,
                fates_hydro_solver = 4)
        finally
            CLM.varctl.use_fates, CLM.varctl.use_fates_planthydro = old
        end
    end

    # -------------------------------------------------------------------------------
    # Data-gated parity vs the Fortran fates_parity_1pt runs with each 2-D solver.
    # -------------------------------------------------------------------------------
    @testset "parity vs Fortran (data-gated)" begin
        sitedir = symfluence_path("domain_Bow_at_Banff_lumped_era5", "settings", "CLM", "parameters")
        fsurdat = joinpath(sitedir, "surfdata_clm.nc")
        forcing = symfluence_path("domain_Bow_at_Banff_lumped_era5", "data", "forcing",
                                  "CLM_input", "clmforc.2002.nc")
        if !(isfile(fsurdat) && isfile(forcing))
            testdata_missing("FATES planthydro 2-D solver parity", fsurdat, forcing)
        else
            pdir = joinpath(@__DIR__, "..", "scripts", "validation", "fates_fortran_parity")
            harness = joinpath(pdir, "fates_planthydro_parity.jl")
            ENV["SYMFLUENCE_DATA"] = symfluence_data_root()
            old = (CLM.hlm_use_planthydro[], CLM.varctl.use_fates_planthydro)
            function load_harness(solver)
                ENV["FATES_HYDRO_SOLVER"] = string(solver)
                m = Module(Symbol("_PH2DHarness", solver))
                Core.eval(m, :(include(p::AbstractString) = Base.include($m, p)))
                try
                    Base.include(m, harness)
                finally
                    delete!(ENV, "FATES_HYDRO_SOLVER")
                end
                return sym -> getglobal(m, sym)
            end
            try
                # --- Newton (3): the Fortran endruns on the first hydraulics call. ---
                # Every round's third Newton iteration has a NaN residual, so the substep
                # is halved 11 times (3600 s -> 3.5 s) and the solve aborts for the first
                # cohort (pft 1, dbh 0.764 cm). The port must fail the same way, along the
                # same iterations (fates_planthydro_newton_fortran_trace.txt).
                G3 = load_harness(CLM.hydr_solver_2DNewton)
                tr = Any[]
                CLM._HYDRO2D_TRACE[] = (k, a...) -> (k === :newton && push!(tr, a); nothing)
                err = try
                    Base.invokelatest(G3(:run_julia), 2); nothing
                catch e
                    e
                finally
                    CLM._HYDRO2D_TRACE[] = nothing
                end
                @test err !== nothing
                @test occursin("Newton hydraulics solve could not converge", sprint(showerror, err))
                ftr = [split(l) for l in eachline(joinpath(pdir, "fates_planthydro_newton_fortran_trace.txt"))]
                @test length(tr) == length(ftr) == 33
                for (a, f) in zip(tr, ftr)
                    ft, dbh, nsteps, it, amax, ssum, dt, tm = a
                    F = parse.(Float64, f[2:end])   # ft dbh nsteps iter nsd amax sum dtime tm
                    @test ft == Int(F[1]) && nsteps == Int(F[3]) && it == Int(F[4])
                    @test dbh ≈ F[2] rtol = 1e-12
                    @test amax ≈ F[6] rtol = 1e-6
                    @test isnan(ssum) == isnan(F[7])
                    @test dt == F[8] && tm == F[9]
                end

                # --- Picard (2): stepped on the host vs the committed hourly extract. ---
                G2 = load_harness(CLM.hydr_solver_2DPicard)
                nsteps = 24
                recs, _ = Base.invokelatest(G2(:run_julia), nsteps)
                ftp = Base.invokelatest(G2(:read_fortran),
                                        joinpath(pdir, "fates_planthydro_picard_fortran.csv"))
                rel(f, k) = abs(recs[k][f] - ftp[f][k]) / abs(ftp[f][k])
                ks = 3:nsteps+1                  # record k <-> nstep k-1; skip init records
                for (f, tol) in PICARD_TOL
                    @test maximum(rel(f, k) for k in ks) < tol
                end
                errs = [abs(r["ERRH2O_COL"]) for r in recs[3:end]]
                @test maximum(errs) < 1e-8
            finally
                CLM.hlm_use_planthydro[], CLM.varctl.use_fates_planthydro = old
            end
        end
    end
end
