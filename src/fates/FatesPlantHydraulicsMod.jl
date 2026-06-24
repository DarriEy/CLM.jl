# FatesPlantHydraulicsMod.jl
# Julia port of FATES src/fates/biogeophys/FatesPlantHydraulicsMod.F90
#
# The FATES plant-hydraulics engine (Christoffersen et al., GMD 2016; TFS v.1-Hydro).
# This module initializes the per-cohort hydraulic node network (soil -> absorbing
# root -> transporting root -> stem -> leaf, plus rhizosphere shells), computes water
# potentials / fractional conductances via the WRF/WKF functions (FatesHydroWTFMod),
# solves the transient water-flow (Richards-style) ODE down each soil-layer path with
# an implicit first-order Taylor (1D) solver, applies the resulting `btran` water
# stress, and performs the water mass-balance bookkeeping.
#
# --------------------------------------------------------------------------------------
# Translation notes (per project conventions)
#   * fates_r8 -> Float64. Fortran subroutine names are preserved.
#   * Operates on the merged ed_site_type / fates_patch_type / fates_cohort_type and the
#     ed_site_hydr_type / ed_cohort_hydr_type from FatesHydraulicsMemMod, and calls the
#     merged WRF/WKF methods (th_from_psi / psi_from_th / ftc_from_psi / derivatives) from
#     FatesHydroWTFMod, allometry from FatesAllometryMod, and PRT GetState.
#   * Fortran's polymorphic global pointer arrays `wrf_plant(pm,ft)` / `wkf_plant(pm,ft)`
#     (pm runs from stomata_p_media=0 .. n_plant_media=4) become module-global
#     `Matrix{Union{WRFType,Nothing}}` indexed [pm+1, ft] (pm 0 -> row 1), wrapped in a
#     Ref so InitHydroGlobals! can (re)allocate. Helper getters `wrf_plant(pm,ft)` /
#     `wkf_plant(pm,ft)` preserve the Fortran call shape and apply the +1 offset.
#   * The Fortran dense / tri-diagonal linear solves are replaced by a Thomas tri-diagonal
#     solve (Hydraulics_Tridiagonal) whose ASSEMBLY math is identical to the Fortran. The
#     Newton/Picard convergence loop and the water mass-balance + error handling of the
#     default 1D solver are preserved verbatim.
#
# --------------------------------------------------------------------------------------
# Solver coverage
#   The FATES default `hydr_solver` is hydr_solver_1DTaylor. This module ports that path
#   FULLY end-to-end (OrderLayersForSolve1D + ImTaylorSolve1D + GetImTaylorKAB +
#   Hydraulics_Tridiagonal + all geometry/kmax setup + init + psi/ftc update + btran +
#   mass balance). The two NON-DEFAULT 2D solvers (hydr_solver_2DNewton -> MatSolve2D,
#   hydr_solver_2DPicard -> PicardSolve2D) are STUBBED: calling hydraulics_bc with those
#   solver settings raises a clear error. GetKAndDKDPsi / SetMaxCondConnections (used only
#   by the 2D solvers) are ported as faithful helpers for future 2D work.
#
# Deps: FatesConstantsMod, FatesGlobals, EDParamsMod (EDParams), EDTypesMod (ed_site_type,
#       area), FatesPatchMod, FatesCohortMod, FatesHydraulicsMemMod, FatesHydroWTFMod,
#       FatesAllometryMod (bleaf/bsap_allom), PRTGenericMod (GetState + organ/element ids),
#       PRTParametersMod (prt_params), EDPftvarcon (EDPftvarcon_inst).

# ---------------------------------------------------------------------------
# Module-level parameters (Fortran `real(r8), parameter` / scalars)
# ---------------------------------------------------------------------------

# 0 => use vanilla btran; 1 => use BC hydraulics; 2 => use CX hydraulics
const use_ed_planthydraulics = 1

# Pressure-volume / pressure-conductivity hypothesis index flags
const van_genuchten_type    = 2
const campbell_type         = 3
const smooth1_campbell_type = 31
const smooth2_campbell_type = 32
const tfs_type              = 1

# Soil hypothesis (hard-coded to Campbell, matching Fortran ELM/ALM default)
const soil_wrf_type = campbell_type
const soil_wkf_type = campbell_type

# Behavior switches (Fortran module logicals/parameters)
const purge_supersaturation = false
const do_growthrecruiteffects = true
const do_upstream_k    = true     # parameter, true
const do_parallel_stem = true     # treat the conduit as closed root-layer -> stomata
const trap_neg_wc      = false
const trap_supersat_psi = false
const hydr_debug       = false

const error_thresh   = 1.0e-5     # site conservation error threshold [mm = kg/m2]
const thsat_buff     = 0.001      # buffer left between soil moisture and saturation [m3/m3]
const max_wb_step_err = 2.0e-6    # max allowable water-balance error over a continuum step [kg]

# Testing parameters for Van Genuchten soil WRTs (unused unless VG soil is selected)
const alpha_vg      = 0.001
const th_sat_vg     = 0.65
const th_res_vg     = 0.15
const psd_vg        = 2.7
const m_vg          = 0.62963
const soil_tort_vg  = 0.5
const plant_tort_vg = 0.0

# AREA / AREA_INV: notional forest area [m2]. EDTypesMod defines `area`.
const AREA_INV = 1.0 / area

# ---------------------------------------------------------------------------
# Global plant Water Transfer Functions, indexed [pm+1, ft]
# (Fortran wrf_plant(stomata_p_media:n_plant_media, numpft))
# ---------------------------------------------------------------------------
const _wrf_plant = Ref{Matrix{Union{WRFType,Nothing}}}(Matrix{Union{WRFType,Nothing}}(undef, 0, 0))
const _wkf_plant = Ref{Matrix{Union{WKFType,Nothing}}}(Matrix{Union{WKFType,Nothing}}(undef, 0, 0))

"Access the plant water-retention function for porous-media `pm` (0..4) and PFT `ft`."
wrf_plant(pm::Integer, ft::Integer) = _wrf_plant[][pm + 1, ft]
"Access the plant water-conductance function for porous-media `pm` (0..4) and PFT `ft`."
wkf_plant(pm::Integer, ft::Integer) = _wkf_plant[][pm + 1, ft]

# ===========================================================================
# InitHydroGlobals : allocate plant WRF/WKF (pm x ft) and set their parameters
# ===========================================================================

"""
    InitHydroGlobals!(numpft::Integer)

Allocate and parameterize the plant Water Transfer Functions for every porous-media
compartment (leaf/stem/troot/aroot) and PFT, plus the stomata-media conductance function.
Mirrors the Fortran `InitHydroGlobals`. Parameters come from `EDPftvarcon_inst[]` and the
`hydr_htftype_node` hypothesis switch in `EDParams[]`.
"""
function InitHydroGlobals!(numpft::Integer)
    pcon = EDPftvarcon_inst[]
    eparm = EDParams[]
    hydr_htftype_node = eparm.hydr_htftype_node

    # rows 1..(n_plant_media+1) correspond to pm = 0 (stomata) .. n_plant_media (=4)
    wrfp = Matrix{Union{WRFType,Nothing}}(nothing, n_plant_media + 1, numpft)
    wkfp = Matrix{Union{WKFType,Nothing}}(nothing, n_plant_media + 1, numpft)

    # --- Water Retention Functions (plant media 1..n_plant_media) ---
    for pm in 1:n_plant_media
        htype = hydr_htftype_node[pm]
        if htype == van_genuchten_type
            for ft in 1:numpft
                wrf = wrf_type_vg()
                set_wrf_param!(wrf, [pcon.hydr_vg_alpha_node[ft, pm],
                                     pcon.hydr_vg_n_node[ft, pm],
                                     pcon.hydr_vg_m_node[ft, pm],
                                     pcon.hydr_thetas_node[ft, pm],
                                     pcon.hydr_resid_node[ft, pm]])
                wrfp[pm + 1, ft] = wrf
            end
        elseif htype == tfs_type
            for ft in 1:numpft
                wrf = wrf_type_tfs()
                if pm == leaf_p_media
                    cap_slp  = 0.0
                    cap_int  = 0.0
                    cap_corr = 1.0
                else
                    cap_slp  = (eparm.hydr_psi0 - eparm.hydr_psicap) / (1.0 - rwccap[pm])
                    cap_int  = -cap_slp + eparm.hydr_psi0
                    cap_corr = -cap_int / cap_slp
                end
                set_wrf_param!(wrf, [pcon.hydr_thetas_node[ft, pm],
                                     pcon.hydr_resid_node[ft, pm],
                                     pcon.hydr_pinot_node[ft, pm],
                                     pcon.hydr_epsil_node[ft, pm],
                                     rwcft[pm],
                                     cap_corr,
                                     cap_int,
                                     cap_slp, Float64(pm)])
                wrfp[pm + 1, ft] = wrf
            end
        else
            fates_endrun("undefined water retention type for plants, pm=$pm type=$htype")
        end
    end

    # --- Water Conductance Functions (plant media 1..n_plant_media) ---
    for pm in 1:n_plant_media
        htype = hydr_htftype_node[pm]
        if htype == van_genuchten_type
            for ft in 1:numpft
                wkf = wkf_type_vg()
                set_wkf_param!(wkf, [pcon.hydr_vg_alpha_node[ft, pm],
                                     pcon.hydr_vg_n_node[ft, pm],
                                     pcon.hydr_vg_m_node[ft, pm],
                                     pcon.hydr_thetas_node[ft, pm],
                                     pcon.hydr_resid_node[ft, pm],
                                     plant_tort_vg])
                wkfp[pm + 1, ft] = wkf
            end
        elseif htype == tfs_type
            for ft in 1:numpft
                wkf = wkf_type_tfs()
                set_wkf_param!(wkf, [pcon.hydr_p50_node[ft, pm],
                                     pcon.hydr_avuln_node[ft, pm]])
                wkfp[pm + 1, ft] = wkf
            end
        else
            fates_endrun("undefined water conductance type for plants, pm=$pm type=$htype")
        end
        # The FTC functions need psi_min -> point each WKF at its matching WRF.
        for ft in 1:numpft
            wkfp[pm + 1, ft].wrf = wrfp[pm + 1, ft]
        end
    end

    # --- Stomatal conductance (single hypothesis: TFS p50/avuln_gs) ---
    for ft in 1:numpft
        wkf = wkf_type_tfs()
        set_wkf_param!(wkf, [pcon.hydr_p50_gs[ft], pcon.hydr_avuln_gs[ft]])
        wkfp[stomata_p_media + 1, ft] = wkf
    end
    # Stomata FTC points at the internal leaf retention structure (pm=1)
    for ft in 1:numpft
        wkfp[stomata_p_media + 1, ft].wrf = wrfp[leaf_p_media + 1, ft]
    end

    _wrf_plant[] = wrfp
    _wkf_plant[] = wkfp
    return nothing
end

# ===========================================================================
# Geometry helpers : rooting depth, root fraction, rhizosphere shells, xylem taper
# ===========================================================================

"""
    MaximumRootingDepth(dbh, ft, z_max_soil) -> z_fr

Maximum rooting depth of the plant (positive-down [m]), an exponential constrained by the
maximum soil depth. Mirrors the Fortran `MaximumRootingDepth` (Ding 2021 dynamic roots).
"""
function MaximumRootingDepth(dbh::Real, ft::Integer, z_max_soil::Real)
    p = prt_params
    dbh_max  = p.allom_zroot_max_dbh[ft]
    dbh_0    = p.allom_zroot_min_dbh[ft]
    z_fr_max = p.allom_zroot_max_z[ft]
    z_fr_0   = p.allom_zroot_min_z[ft]
    frk      = p.allom_zroot_k[ft]

    dbh_rel = min(1.0, (max(dbh, dbh_0) - dbh_0) / (dbh_max - dbh_0))
    z_fr = min(z_max_soil,
               z_fr_max / (1.0 + ((z_fr_max - z_fr_0) / z_fr_0) * exp(-frk * dbh_rel)))
    return z_fr
end

"""
    zeng2001_crootfr(a, b, z; z_max=nothing) -> crootfr

Cumulative root fraction at depth `z` (Zeng 2001). If `z_max` is given, normalize so the
profile sums to unity at the maximum rooting depth. Mirrors the Fortran `zeng2001_crootfr`.
"""
function zeng2001_crootfr(a::Real, b::Real, z::Real, z_max::Union{Real,Nothing}=nothing)
    crootfr = 1.0 - 0.5 * (exp(-a * z) + exp(-b * z))
    if z_max !== nothing
        zc = min(z, z_max)
        crootfr = 1.0 - 0.5 * (exp(-a * zc) + exp(-b * zc))
        crootfr_max = 1.0 - 0.5 * (exp(-a * z_max) + exp(-b * z_max))
        crootfr = crootfr / crootfr_max
    end
    return crootfr
end

"""
    bisect_rootfr(a, b, z_max, lower_init, upper_init, xtol, ytol, crootfr) -> x_new

Bisection inverse of the cumulative root distribution: depth `x_new` at which the
cumulative root fraction equals `crootfr`. Mirrors the Fortran `bisect_rootfr`.
"""
function bisect_rootfr(a::Real, b::Real, z_max::Real, lower_init::Real, upper_init::Real,
                       xtol::Real, ytol::Real, crootfr::Real)
    lower = lower_init
    upper = upper_init
    f_lo  = zeng2001_crootfr(a, b, lower, z_max) - crootfr
    f_hi  = zeng2001_crootfr(a, b, upper, z_max) - crootfr
    chg   = upper - lower
    nitr  = 0
    x_new = 0.5 * (lower + upper)
    while abs(chg) > xtol
        x_new = 0.5 * (lower + upper)
        f_new = zeng2001_crootfr(a, b, x_new, z_max) - crootfr
        if abs(f_new) <= ytol
            break
        end
        if (f_lo * f_new) < 0.0
            upper = x_new
        end
        if (f_hi * f_new) < 0.0
            lower = x_new
        end
        chg = upper - lower
        nitr += 1
        nitr > 100 && break
    end
    if nitr > 100
        @warn "Warning: number of iterations exceeds 100 for bisect_rootfr"
    end
    return x_new
end

"""
    shellGeom!(l_aroot_in, rs1_in, area_site, dz, r_out_shell, r_node_shell, v_shell)

Update the 'representative' rhizosphere geometry (outer radii, node radii, shell volumes)
for the given total absorbing-root length in a soil layer. Writes into the passed shell
arrays (length `nshell`). Mirrors the Fortran `shellGeom` (Sperry 1998 eqns 7,8).
"""
function shellGeom!(l_aroot_in::Real, rs1_in::Real, area_site::Real, dz::Real,
                    r_out_shell::AbstractVector, r_node_shell::AbstractVector,
                    v_shell::AbstractVector)
    nshells = length(r_out_shell)

    if l_aroot_in <= nearzero
        r_out_shell  .= 0.0
        r_node_shell .= 0.0
        v_shell      .= 0.0
        return nothing
    end
    rs1 = rs1_in
    l_aroot = l_aroot_in

    # Outer radii (eqn 8 then eqn 7, Sperry 1998)
    r_out_shell[nshells] = (pi_const * l_aroot / (area_site * dz))^(-0.5)
    if nshells > 1
        for k in 1:(nshells - 1)
            r_out_shell[k] = rs1 * (r_out_shell[nshells] / rs1)^(Float64(k) / Float64(nshells))
        end
    end

    # Nodal (midpoint) radii
    r_node_shell[1] = 0.5 * (rs1 + r_out_shell[1])
    for k in 2:nshells
        r_node_shell[k] = 0.5 * (r_out_shell[k-1] + r_out_shell[k])
    end

    # Volumes
    for k in 1:nshells
        if k == 1
            v_shell[k] = pi_const * l_aroot * (r_out_shell[k]^2 - rs1^2)
        else
            v_shell[k] = pi_const * l_aroot * (r_out_shell[k]^2 - r_out_shell[k-1]^2)
        end
    end
    return nothing
end

"""
    xylemtaper(pexp, dz) -> chi_tapnotap

Ratio of total tree conductance accounting for xylem taper to that without, over the
hydraulic distance `dz` (Savage et al. 2010). Mirrors the Fortran `xylemtaper`.
"""
function xylemtaper(pexp::Real, dz::Real)
    lN    = 0.005    # petiole length [m]
    n_ext = 2.0      # daughter branches per parent branch
    a5 = -3.555547; a4 = 9.760275; a3 = -8.468005
    a2 = 1.096488;  a1 = 1.844792; a0 = 1.320732   # a0 unused (matches Fortran)

    qexp = a5 * pexp^5 + a4 * pexp^4 + a3 * pexp^3 + a2 * pexp^2 + a1 * pexp
    num  = 3.0 * log(1.0 - dz / lN * (1.0 - n_ext^(1.0 / 3.0)))
    den  = log(n_ext)
    big_N = num / den - 1.0
    r0rN  = n_ext^(big_N / 2.0)
    return r0rN^qexp
end

# ===========================================================================
# Plant node geometry, volumes/lengths, and maximum conductances
# ===========================================================================

"""
    UpdatePlantHydrNodes!(ccohort, ft, plant_height, csite_hydr)

Compute the nodal heights of the leaf, stem, and transporting-root compartments for a
cohort. Writes `z_node_ag`, `z_lower_ag`, `z_upper_ag`, `z_node_troot`. Mirrors the
Fortran `UpdatePlantHydrNodes`.
"""
function UpdatePlantHydrNodes!(ccohort, ft::Integer, plant_height::Real, csite_hydr)
    ch = ccohort.co_hydr
    roota = prt_params.fnrt_prof_a[ft]
    rootb = prt_params.fnrt_prof_b[ft]
    nlevrhiz = csite_hydr.nlevrhiz

    crown_depth = min(plant_height, 0.1)

    # Crown (leaf) nodes
    dz_canopy = crown_depth / Float64(n_hypool_leaf)
    for k in 1:n_hypool_leaf
        ch.z_lower_ag[k] = plant_height - dz_canopy * Float64(k)
        ch.z_node_ag[k]  = ch.z_lower_ag[k] + 0.5 * dz_canopy
        ch.z_upper_ag[k] = ch.z_lower_ag[k] + dz_canopy
    end

    # Stem nodes
    z_stem  = plant_height - crown_depth
    dz_stem = z_stem / Float64(n_hypool_stem)
    for k in (n_hypool_leaf + 1):n_hypool_ag
        ch.z_upper_ag[k] = Float64(n_hypool_stem - (k - 1 - n_hypool_leaf)) * dz_stem
        ch.z_node_ag[k]  = ch.z_upper_ag[k] - 0.5 * dz_stem
        ch.z_lower_ag[k] = ch.z_upper_ag[k] - dz_stem
    end

    # Transporting-root node depth (negative from surface)
    z_fr = MaximumRootingDepth(ccohort.dbh, ft, csite_hydr.zi_rhiz[nlevrhiz])
    z_cumul_rf = bisect_rootfr(roota, rootb, z_fr, 0.0, 1.0e10, 0.001, 0.001, 0.5)

    if z_cumul_rf > csite_hydr.zi_rhiz[nlevrhiz]
        fates_endrun("z_cumul_rf > zi_rhiz(nlevrhiz)? $z_cumul_rf $(csite_hydr.zi_rhiz[nlevrhiz])")
    end
    z_cumul_rf = min(z_cumul_rf, abs(csite_hydr.zi_rhiz[nlevrhiz]))
    ch.z_node_troot = -z_cumul_rf
    return nothing
end

"""
    SavePreviousCompartmentVolumes!(ccohort_hydr)

Save the current plant compartment volumes into the `*_init` save-space, so size changes
can be detected. Mirrors the Fortran `SavePreviousCompartmentVolumes`.
"""
function SavePreviousCompartmentVolumes!(ch::ed_cohort_hydr_type)
    ch.v_ag_init          = copy(ch.v_ag)
    ch.v_troot_init       = ch.v_troot
    ch.v_aroot_layer_init = copy(ch.v_aroot_layer)
    return nothing
end

"""
    UpdatePlantHydrLenVol!(ccohort, csite_hydr)

Compute the plant compartment storage volumes [m3] and absorbing-root lengths [m] (per
soil layer) from the cohort's carbon pools and geometry. Mirrors the Fortran
`UpdatePlantHydrLenVol`.
"""
function UpdatePlantHydrLenVol!(ccohort, csite_hydr)
    ch = ccohort.co_hydr
    ft = ccohort.pft
    p  = prt_params
    nlevrhiz = csite_hydr.nlevrhiz

    leaf_c   = GetState(ccohort.prt, leaf_organ, carbon12_element)
    sapw_c   = GetState(ccohort.prt, sapw_organ, carbon12_element)
    fnrt_c   = GetState(ccohort.prt, fnrt_organ, carbon12_element)
    struct_c = GetState(ccohort.prt, struct_organ, carbon12_element)
    roota = p.fnrt_prof_a[ft]
    rootb = p.fnrt_prof_b[ft]

    l2sap_vol_donate_frac = 0.5
    min_leaf_frac = 0.1
    min_trim = 0.1

    # --- Leaf volumes ---
    sla     = p.slatop[ft] * cm2_per_m2                 # cm2/gC
    denleaf = -2.3231 * sla / p.c2b[ft] + 781.899       # kg/m3

    # Target leaf carbon (efleaf_coh hard-coded to 1 to avoid zero leaf volume)
    leaf_c_target, _ = bleaf(ccohort.dbh, ccohort.pft, ccohort.crowndamage,
                             max(ccohort.canopy_trim, min_trim), 1.0)

    for k in 1:n_hypool_leaf
        ch.v_ag[k] = max(leaf_c, min_leaf_frac * leaf_c_target) *
                     p.c2b[ft] / denleaf / Float64(n_hypool_leaf)
    end

    # --- Sapwood / stem volume ---
    a_sapwood_target, _ = bsap_allom(ccohort.dbh, ccohort.pft, ccohort.crowndamage,
                                     ccohort.canopy_trim, ccohort.efstem_coh)
    a_sapwood = a_sapwood_target

    crown_depth = min(ccohort.height, 0.1)
    z_stem      = ccohort.height - crown_depth
    v_sapwood   = a_sapwood * z_stem

    # Above-ground node volume: woody plants keep stem; grass donates half leaf to xylem
    if p.woody[ft] == 1
        for k in (n_hypool_leaf + 1):n_hypool_ag
            ch.v_ag[k] = v_sapwood / n_hypool_stem
        end
    else
        v_leaf_donate = zeros(Float64, n_hypool_leaf)
        for k in 1:n_hypool_leaf
            v_leaf_donate[k] = ch.v_ag[k] * l2sap_vol_donate_frac
            ch.v_ag[k] = ch.v_ag[k] - v_leaf_donate[k]
        end
        for k in (n_hypool_leaf + 1):n_hypool_ag
            ch.v_ag[k] = (v_sapwood + sum(v_leaf_donate)) / n_hypool_stem
        end
    end

    # --- Transporting (coarse) root volume ---
    woody_bg_c = (1.0 - p.allom_agb_frac[ft]) * (sapw_c + struct_c)
    v_troot = woody_bg_c * p.c2b[ft] / (p.wood_density[ft] * kg_per_g * cm3_per_m3)

    # --- Absorbing root total length and volume ---
    pcon = EDPftvarcon_inst[]
    l_aroot_tot = fnrt_c * g_per_kg * p.c2b[ft] * pcon.hydr_srl[ft]
    v_aroot_tot = pi_const * (pcon.hydr_rs2[ft]^2) * l_aroot_tot

    # New method: troot & aroot each get 50% of the sum of both
    ch.v_troot = 0.5 * (v_troot + v_aroot_tot)

    # Partition root length/volume into the active soil layers (Zeng + max depth)
    z_fr = MaximumRootingDepth(ccohort.dbh, ft, csite_hydr.zi_rhiz[nlevrhiz])
    for j in 1:nlevrhiz
        rootfr = zeng2001_crootfr(roota, rootb, csite_hydr.zi_rhiz[j], z_fr) -
                 zeng2001_crootfr(roota, rootb, csite_hydr.zi_rhiz[j] - csite_hydr.dz_rhiz[j], z_fr)
        ch.l_aroot_layer[j] = rootfr * l_aroot_tot
        ch.v_aroot_layer[j] = rootfr * 0.5 * (v_aroot_tot + v_troot)
    end
    return nothing
end

"""
    UpdatePlantKmax!(ccohort_hydr, ccohort, csite_hydr)

Set the maximum hydraulic conductances [kg H2O s-1 MPa-1] of every plant compartment from
leaf -> stem -> transporting root -> absorbing roots, given the compartment geometry.
Mirrors the Fortran `UpdatePlantKmax`.
"""
function UpdatePlantKmax!(ch::ed_cohort_hydr_type, ccohort, csite_hydr)
    pft  = ccohort.pft
    pcon = EDPftvarcon_inst[]
    min_pet_stem_dz = 0.00001

    a_sapwood, _ = bsap_allom(ccohort.dbh, pft, ccohort.crowndamage,
                              ccohort.canopy_trim, ccohort.efstem_coh)

    # Leaf maximum conductance: assume negligible internal resistance (nominally very high)
    ch.kmax_petiole_to_leaf = 1.0e8

    kmax_node_stem = pcon.hydr_kmax_node[pft, 2]   # stem media (pm index 2)
    p_taper        = pcon.hydr_p_taper[pft]

    # Stem maximum conductances
    for k in 1:n_hypool_stem
        k_ag = k + n_hypool_leaf
        z_lower = ch.z_node_ag[n_hypool_leaf] - ch.z_lower_ag[k_ag]
        z_node  = ch.z_node_ag[n_hypool_leaf] - ch.z_node_ag[k_ag]
        z_upper = max(min_pet_stem_dz, ch.z_node_ag[n_hypool_leaf] - ch.z_upper_ag[k_ag])

        kmax_upper = kmax_node_stem * xylemtaper(p_taper, z_upper) * a_sapwood / z_upper
        kmax_node  = kmax_node_stem * xylemtaper(p_taper, z_node)  * a_sapwood / z_node
        kmax_lower = kmax_node_stem * xylemtaper(p_taper, z_lower) * a_sapwood / z_lower

        ch.kmax_stem_upper[k] = (1.0 / kmax_node - 1.0 / kmax_upper)^(-1.0)
        ch.kmax_stem_lower[k] = (1.0 / kmax_lower - 1.0 / kmax_node)^(-1.0)
    end

    # Upper transporting-root compartment (connects to the lowest stem)
    z_upper = ch.z_lower_ag[n_hypool_leaf]
    z_node  = ch.z_lower_ag[n_hypool_leaf] - ch.z_node_troot
    kmax_node  = kmax_node_stem * xylemtaper(p_taper, z_node)  * a_sapwood / z_node
    kmax_upper = kmax_node_stem * xylemtaper(p_taper, z_upper) * a_sapwood / z_upper
    ch.kmax_troot_upper = (1.0 / kmax_node - 1.0 / kmax_upper)^(-1.0)

    # Below-ground kmax as a residual of total resistance, split per layer by root length
    rmin_ag = 1.0 / ch.kmax_petiole_to_leaf +
              sum(1.0 ./ ch.kmax_stem_upper[1:n_hypool_stem]) +
              sum(1.0 ./ ch.kmax_stem_lower[1:n_hypool_stem]) +
              1.0 / ch.kmax_troot_upper

    kmax_bg = 1.0 / (rmin_ag * (1.0 / pcon.hydr_rfrac_stem[pft] - 1.0))

    sum_l_aroot = sum(ch.l_aroot_layer)
    for j in 1:csite_hydr.nlevrhiz
        kmax_layer = kmax_bg * ch.l_aroot_layer[j] / sum_l_aroot
        ch.kmax_troot_lower[j] = 3.0 * kmax_layer
        ch.kmax_aroot_upper[j] = 3.0 * kmax_layer
        ch.kmax_aroot_lower[j] = 3.0 * kmax_layer
    end

    # Radial conductance (root surface -> center node), both gradient directions
    for j in 1:csite_hydr.nlevrhiz
        surfarea_aroot_layer = 2.0 * pi_const * pcon.hydr_rs2[ccohort.pft] * ch.l_aroot_layer[j]
        ch.kmax_aroot_radial_in[j]  = EDParams[].hydr_kmax_rsurf1 * surfarea_aroot_layer
        ch.kmax_aroot_radial_out[j] = EDParams[].hydr_kmax_rsurf2 * surfarea_aroot_layer
    end
    return nothing
end

"""
    UpdateSizeDepPlantHydProps!(currentSite, ccohort, bc_in=nothing)

Update absorbing-root length, plant node positions, compartment volumes/lengths, and
maximum conductances after the plant grows. Mirrors the Fortran
`UpdateSizeDepPlantHydProps`.
"""
function UpdateSizeDepPlantHydProps!(currentSite, ccohort, bc_in=nothing)
    csite_hydr = currentSite.si_hydr
    ch = ccohort.co_hydr
    ft = ccohort.pft
    SavePreviousCompartmentVolumes!(ch)
    UpdatePlantHydrNodes!(ccohort, ft, ccohort.height, csite_hydr)
    UpdatePlantHydrLenVol!(ccohort, csite_hydr)
    UpdatePlantKmax!(ch, ccohort, csite_hydr)
    return nothing
end

# ===========================================================================
# State initialization & psi/ftc updates
# ===========================================================================

"""
    InitPlantHydStates!(site, cohort)

Initialize the cohort's compartment water contents (`th_*`), water potentials (`psi_*`),
and fractional conductances (`ftc_*`) in (near) hydrostatic equilibrium with the soil, and
set the initial `btran`. Uses init_mode 2 (match soil totals, then offset by a small
gradient). Mirrors the Fortran `InitPlantHydStates`.
"""
function InitPlantHydStates!(site, cohort)
    csite_hydr = site.si_hydr
    ch = cohort.co_hydr
    ft = cohort.pft
    wrfa = wrf_plant(aroot_p_media, ft)
    wkfa = wkf_plant(aroot_p_media, ft)
    wrft = wrf_plant(troot_p_media, ft)
    wkft = wkf_plant(troot_p_media, ft)

    psi_aroot_init = -0.2     # MPa
    dh_dz = 0.02              # MPa/m offset added downstream
    init_mode = 2

    # Absorbing roots
    if init_mode == 2
        for j in 1:csite_hydr.nlevrhiz
            if ch.l_aroot_layer[j] > nearzero
                ch.psi_aroot[j] = psi_from_th(csite_hydr.wrf_soil[j],
                                              csite_hydr.h2osoi_liqvol_shell[j, 1])
                ch.th_aroot[j]  = max(th_from_psi(wrfa, ch.psi_aroot[j]), get_thmin(wrfa))
                ch.ftc_aroot[j] = ftc_from_psi(wkfa, ch.psi_aroot[j])
            else
                ch.psi_aroot[j] = psi_aroot_init
                ch.th_aroot[j]  = 0.0
            end
        end
    else
        for j in 1:csite_hydr.nlevrhiz
            ch.psi_aroot[j] = psi_aroot_init
            ch.th_aroot[j]  = max(th_from_psi(wrfa, ch.psi_aroot[j]), get_thmin(wrfa))
            ch.ftc_aroot[j] = ftc_from_psi(wkfa, ch.psi_aroot[j])
        end
    end

    # Mean total potential of absorbing roots (using layer centers)
    h_aroot_mean = minimum(ch.psi_aroot .+ mpa_per_pa * dens_fresh_liquid_water * grav_earth .*
                           (-csite_hydr.zi_rhiz .+ 0.5 .* csite_hydr.dz_rhiz))

    # Transporting root in equilibrium with mean aroot potential, minus offset
    ch.psi_troot = h_aroot_mean -
                   mpa_per_pa * dens_fresh_liquid_water * grav_earth * ch.z_node_troot - dh_dz
    ch.th_troot  = max(th_from_psi(wrft, ch.psi_troot), get_thmin(wrft))
    ch.ftc_troot = ftc_from_psi(wkft, ch.psi_troot)

    # Top above-ground (stem) compartment
    dz = ch.z_node_ag[n_hypool_ag] - ch.z_node_troot
    ch.psi_ag[n_hypool_ag] = ch.psi_troot -
                             mpa_per_pa * dens_fresh_liquid_water * grav_earth * dz - dh_dz
    wrf_stem = wrf_plant(stem_p_media, ft)
    wkf_stem = wkf_plant(stem_p_media, ft)
    ch.th_ag[n_hypool_ag]  = max(get_thmin(wrf_stem), th_from_psi(wrf_stem, ch.psi_ag[n_hypool_ag]))
    ch.ftc_ag[n_hypool_ag] = ftc_from_psi(wkf_stem, ch.psi_ag[n_hypool_ag])

    # Work up the plant (stem -> leaf), hydrostatic equilibrium minus dh_dz each step
    for k in (n_hypool_ag - 1):-1:1
        dz = ch.z_node_ag[k] - ch.z_node_ag[k+1]
        ch.psi_ag[k] = ch.psi_ag[k+1] -
                       mpa_per_pa * dens_fresh_liquid_water * grav_earth * dz - dh_dz
        pm = csite_hydr.pm_node[k]
        wrfk = wrf_plant(pm, ft)
        wkfk = wkf_plant(pm, ft)
        ch.th_ag[k]  = max(th_from_psi(wrfk, ch.psi_ag[k]), get_thmin(wrfk))
        ch.ftc_ag[k] = ftc_from_psi(wkfk, ch.psi_ag[k])
    end

    # Cohort-level btran (stomata FTC from leaf psi)
    ch.btran = ftc_from_psi(wkf_plant(stomata_p_media, ft), ch.psi_ag[1])

    if ch.psi_troot > 0.0 || any(ch.psi_ag .> 0.0) || any(ch.psi_aroot .> 0.0)
        fates_endrun("Initialized plant compartments with positive pressure?")
    end
    return nothing
end

"""
    UpdatePlantPsiFTCFromTheta!(ccohort, csite_hydr)

Update the water potential `psi_*` and fractional conductance `ftc_*` of every plant
compartment from the current water content `th_*`, and recompute `btran`. Mirrors the
Fortran `UpdatePlantPsiFTCFromTheta`.
"""
function UpdatePlantPsiFTCFromTheta!(ccohort, csite_hydr)
    ch = ccohort.co_hydr
    ft = ccohort.pft

    for k in 1:n_hypool_leaf
        ch.psi_ag[k] = psi_from_th(wrf_plant(leaf_p_media, ft), ch.th_ag[k])
        ch.ftc_ag[k] = ftc_from_psi(wkf_plant(leaf_p_media, ft), ch.psi_ag[k])
    end

    ch.btran = ftc_from_psi(wkf_plant(stomata_p_media, ft), ch.psi_ag[1])

    for k in (n_hypool_leaf + 1):n_hypool_ag
        ch.psi_ag[k] = psi_from_th(wrf_plant(stem_p_media, ft), ch.th_ag[k])
        ch.ftc_ag[k] = ftc_from_psi(wkf_plant(stem_p_media, ft), ch.psi_ag[k])
    end

    ch.psi_troot = psi_from_th(wrf_plant(troot_p_media, ft), ch.th_troot)
    ch.ftc_troot = ftc_from_psi(wkf_plant(troot_p_media, ft), ch.psi_troot)

    for j in 1:csite_hydr.nlevrhiz
        ch.psi_aroot[j] = psi_from_th(wrf_plant(aroot_p_media, ft), ch.th_aroot[j])
        ch.ftc_aroot[j] = ftc_from_psi(wkf_plant(aroot_p_media, ft), ch.psi_aroot[j])
    end
    return nothing
end

# ===========================================================================
# Tri-diagonal solver (Thomas algorithm)
# ===========================================================================

"""
    Hydraulics_Tridiagonal(a, b, c, r) -> (u, ierr)

Solve `a(i)*u(i-1) + b(i)*u(i) + c(i)*u(i+1) = r(i)` for `u` (Thomas algorithm), with
`a(1)` and `c(N)` non-existent. Returns the solution and an error flag (0=ok). Mirrors the
Fortran `Hydraulics_Tridiagonal` (assembly identical; debug forward-error check omitted as
in the Fortran non-debug path).
"""
function Hydraulics_Tridiagonal(a::AbstractVector, b::AbstractVector,
                                c::AbstractVector, r::AbstractVector)
    N = length(r)
    u = zeros(Float64, N)
    gam = zeros(Float64, N)
    bet = b[1]
    for k in 1:N
        if k == 1
            u[k] = r[k] / bet
        else
            gam[k] = c[k-1] / bet
            bet    = b[k] - a[k] * gam[k]
            u[k]   = (r[k] - a[k] * u[k-1]) / bet
        end
    end
    for k in (N-1):-1:1
        u[k] = u[k] - gam[k+1] * u[k+1]
    end
    return u, 0
end

# ===========================================================================
# Conductance-and-Taylor-term builders for the implicit (1D) solve
# ===========================================================================

"""
    GetImTaylorKAB(kmax_up, kmax_dn, ftc_up, ftc_dn, h_up, h_dn,
                   dftc_dtheta_up, dftc_dtheta_dn, dpsi_dtheta_up, dpsi_dtheta_dn)
        -> (k_eff, a_term, b_term)

Return the effective conductance over a node-to-node path and the two first-order Taylor
("A","B") terms used by the implicit solve. Applies the upstream-K convention. Mirrors the
Fortran `GetImTaylorKAB` (returns values; does not mutate caller ftc).
"""
function GetImTaylorKAB(kmax_up::Real, kmax_dn::Real,
                        ftc_up::Real, ftc_dn::Real,
                        h_up::Real, h_dn::Real,
                        dftc_dtheta_up::Real, dftc_dtheta_dn::Real,
                        dpsi_dtheta_up::Real, dpsi_dtheta_dn::Real)
    h_diff = h_up - h_dn

    # Upstream K: the side governing FTC; zero the derivative on the other side
    if do_upstream_k
        if h_diff > 0.0
            ftc_dn = ftc_up
            dftc_dtheta_dn = 0.0
        else
            ftc_up = ftc_dn
            dftc_dtheta_up = 0.0
        end
    end

    if (ftc_up * kmax_up) > nearzero && (ftc_dn * kmax_dn) > nearzero
        k_eff = 1.0 / (1.0 / (ftc_up * kmax_up) + 1.0 / (ftc_dn * kmax_dn))
        a_term = k_eff^2 * h_diff * kmax_dn^(-1.0) * ftc_dn^(-2.0) * dftc_dtheta_dn -
                 k_eff * dpsi_dtheta_dn
        b_term = k_eff^2 * h_diff * kmax_up^(-1.0) * ftc_up^(-2.0) * dftc_dtheta_up +
                 k_eff * dpsi_dtheta_up
    else
        k_eff = 0.0; a_term = 0.0; b_term = 0.0
    end
    return k_eff, a_term, b_term
end

"""
    GetKAndDKDPsi(kmax_dn, kmax_up, h_dn, h_up, ftc_dn, ftc_up, dftc_dpsi_dn, dftc_dpsi_up)
        -> (dk_dpsi_dn, dk_dpsi_up, k_eff)

Effective conductance and its sensitivities to downstream/upstream psi over a path. Used by
the (non-default) 2D solvers. Mirrors the Fortran `GetKAndDKDPsi`.
"""
function GetKAndDKDPsi(kmax_dn::Real, kmax_up::Real, h_dn::Real, h_up::Real,
                       ftc_dn::Real, ftc_up::Real, dftc_dpsi_dn::Real, dftc_dpsi_up::Real)
    ftc_dnx = ftc_dn; ftc_upx = ftc_up
    dftc_dpsi_dnx = dftc_dpsi_dn; dftc_dpsi_upx = dftc_dpsi_up
    h_diff = h_up - h_dn
    if do_upstream_k
        if h_diff > 0.0
            ftc_dnx = ftc_up; dftc_dpsi_dnx = 0.0
        else
            ftc_upx = ftc_dn; dftc_dpsi_upx = 0.0
        end
    end
    k_eff = 1.0 / (1.0 / (ftc_upx * kmax_up) + 1.0 / (ftc_dnx * kmax_dn))
    dk_dpsi_dn = k_eff^2 * kmax_dn^(-1.0) * ftc_dnx^(-2.0) * dftc_dpsi_dnx
    dk_dpsi_up = k_eff^2 * kmax_up^(-1.0) * ftc_upx^(-2.0) * dftc_dpsi_upx
    return dk_dpsi_dn, dk_dpsi_up, k_eff
end

# ===========================================================================
# Layer ordering & the default 1D implicit (Taylor) solver
# ===========================================================================

"""
    OrderLayersForSolve1D(csite_hydr, cohort, cohort_hydr, ordered, kbg_layer)

Compute the relative below-ground conductance `kbg_layer` per soil layer and reorder
`ordered` so layers are solved in order of decreasing total root-soil conductance. Both
arrays are mutated in place. Mirrors the Fortran `OrderLayersForSolve1D`.
"""
function OrderLayersForSolve1D(csite_hydr, cohort, cohort_hydr,
                               ordered::AbstractVector{<:Integer}, kbg_layer::AbstractVector)
    neglibible_cond = 1.0e-10
    ft = cohort.pft
    kbg_tot = 0.0
    fill!(kbg_layer, 0.0)

    for j in 1:csite_hydr.nlevrhiz
        if cohort_hydr.l_aroot_layer[j] > nearzero
            psi_inner_shell = psi_from_th(csite_hydr.wrf_soil[j], csite_hydr.h2osoi_liqvol_shell[j, 1])

            if cohort_hydr.psi_aroot[j] < psi_inner_shell
                kmax_aroot = cohort_hydr.kmax_aroot_radial_in[j]
            else
                kmax_aroot = cohort_hydr.kmax_aroot_radial_out[j]
            end

            ftc_aroot = ftc_from_psi(wkf_plant(aroot_p_media, ft), cohort_hydr.psi_aroot[j])
            r_bg = 1.0 / (kmax_aroot * max(ftc_aroot, neglibible_cond))

            aroot_frac_plant = cohort_hydr.l_aroot_layer[j] / csite_hydr.l_aroot_layer[j]
            for k in 1:nshell
                kmax_up = csite_hydr.kmax_upper_shell[j, k] * aroot_frac_plant
                kmax_lo = csite_hydr.kmax_lower_shell[j, k] * aroot_frac_plant
                psi_shell = psi_from_th(csite_hydr.wrf_soil[j], csite_hydr.h2osoi_liqvol_shell[j, k])
                ftc_shell = ftc_from_psi(csite_hydr.wkf_soil[j], psi_shell)
                r_bg += 1.0 / (kmax_up * max(ftc_shell, neglibible_cond))
                if k < nshell
                    r_bg += 1.0 / (kmax_lo * max(ftc_shell, neglibible_cond))
                end
            end
            kbg_layer[j] = 1.0 / r_bg
        else
            kbg_layer[j] = 0.0
        end
        kbg_tot += kbg_layer[j]
    end

    kbg_layer ./= kbg_tot

    # Order layers by decreasing kbg_layer (bubble sort over indices, matching Fortran)
    for j in (csite_hydr.nlevrhiz - 1):-1:1
        for jj in 1:j
            if kbg_layer[ordered[jj]] <= kbg_layer[ordered[jj+1]]
                ordered[jj], ordered[jj+1] = ordered[jj+1], ordered[jj]
            end
        end
    end
    return nothing
end

"""
    ImTaylorSolve1D(slat, slon, recruitflag, csite_hydr, cohort, cohort_hydr, dtime, q_top,
                    ordered, kbg_layer, dth_layershell_col)
        -> (sapflow, rootuptake, wb_err_ps, dwat_plant)

Default FATES hydraulics solver: an implicit first-order-Taylor solution of the transient
flux equations along each soil-layer path (leaf -> stem -> troot -> aroot -> rhizosphere
shells), solved sequentially in `ordered` layer order. Sub-steps are halved/refined until
the per-step water-balance error is acceptable. Updates the cohort `th_*`, accumulates
`dth_layershell_col` (rhizosphere change shared across cohorts), and returns the integrated
sapflow / root-uptake fluxes, the plant-soil balance error, and the plant water change [kg].
Mirrors the Fortran `ImTaylorSolve1D`.
"""
function ImTaylorSolve1D(slat::Real, slon::Real, recruitflag::Bool, csite_hydr, cohort,
                         cohort_hydr, dtime::Real, q_top::Real,
                         ordered::AbstractVector{<:Integer}, kbg_layer::AbstractVector,
                         dth_layershell_col::AbstractMatrix)
    imult = 3
    max_iter = 20
    max_wb_err = 2.0e-5
    no_ftc_radialk = false
    weight_serial_dt = true

    pm_node = csite_hydr.pm_node
    denh2o = dens_fresh_liquid_water
    geo = mpa_per_pa * denh2o * grav_earth

    cohort_hydr.iterh1 = 0.0
    cohort_hydr.iterh2 = 0.0
    wb_err_ps  = 0.0
    dwat_plant = 0.0
    sapflow = 0.0
    rootuptake = zeros(Float64, csite_hydr.nlevrhiz)
    ft = cohort.pft
    sum_l_aroot = sum(cohort_hydr.l_aroot_layer)

    # Node-based scratch arrays (length n_hypool_tot)
    th_node_init = zeros(Float64, n_hypool_tot)
    th_node      = zeros(Float64, n_hypool_tot)
    z_node       = zeros(Float64, n_hypool_tot)
    v_node       = zeros(Float64, n_hypool_tot)
    psi_node     = zeros(Float64, n_hypool_tot)
    ftc_node     = zeros(Float64, n_hypool_tot)
    h_node       = zeros(Float64, n_hypool_tot)
    dftc_dtheta_node = zeros(Float64, n_hypool_tot)
    dpsi_dtheta_node = zeros(Float64, n_hypool_tot)
    k_eff  = zeros(Float64, n_hypool_tot - 1)
    a_term = zeros(Float64, n_hypool_tot - 1)
    b_term = zeros(Float64, n_hypool_tot - 1)
    tris_a = zeros(Float64, n_hypool_tot)
    tris_b = zeros(Float64, n_hypool_tot)
    tris_c = zeros(Float64, n_hypool_tot)
    tris_r = zeros(Float64, n_hypool_tot)

    for jj in 1:csite_hydr.nlevrhiz
        ilayer = ordered[jj]
        cohort_hydr.l_aroot_layer[ilayer] <= nearzero && continue

        if do_parallel_stem
            dt_step = dtime
        else
            dt_step = weight_serial_dt ? dtime * kbg_layer[ilayer] :
                                         dtime / Float64(csite_hydr.nlevrhiz)
        end

        aroot_frac_plant = cohort_hydr.l_aroot_layer[ilayer] / csite_hydr.l_aroot_layer[ilayer]
        wb_err_layer = 0.0
        rootfr_scaler = do_parallel_stem ? cohort_hydr.l_aroot_layer[ilayer] / sum_l_aroot : 1.0
        q_top_eff = q_top * rootfr_scaler

        # Fill node arrays for this layer's path
        for i in 1:n_hypool_tot
            if i <= n_hypool_ag
                z_node[i] = cohort_hydr.z_node_ag[i]
                v_node[i] = cohort_hydr.v_ag[i]
                th_node_init[i] = cohort_hydr.th_ag[i]
            elseif i == n_hypool_ag + 1
                z_node[i] = cohort_hydr.z_node_troot
                v_node[i] = cohort_hydr.v_troot
                th_node_init[i] = cohort_hydr.th_troot
            elseif i == n_hypool_ag + 2
                z_node[i] = -csite_hydr.zi_rhiz[ilayer] + 0.5 * csite_hydr.dz_rhiz[ilayer]
                v_node[i] = cohort_hydr.v_aroot_layer[ilayer]
                th_node_init[i] = cohort_hydr.th_aroot[ilayer]
            else
                ishell = i - (n_hypool_ag + 2)
                z_node[i] = -csite_hydr.zi_rhiz[ilayer] + 0.5 * csite_hydr.dz_rhiz[ilayer]
                v_node[i] = csite_hydr.v_shell[ilayer, ishell] * aroot_frac_plant
                th_node_init[i] = csite_hydr.h2osoi_liqvol_shell[ilayer, ishell]
                if th_node_init[i] < -nearzero
                    fates_endrun("ImTaylorSolve1D: negative shell theta layer=$ilayer shell=$ishell")
                end
            end
        end

        solution_found = false
        iter = 0
        nsteps = 1
        w_tot_beg = 0.0; w_tot_end = 0.0
        sapflow_lyr = 0.0; rootuptake_lyr = 0.0
        dth_node = zeros(Float64, n_hypool_tot)

        while !solution_found
            if iter > max_iter
                fates_endrun("ImTaylorSolve1D: exceeded max_iter; water balance not converging" *
                             " (layer=$ilayer, lat=$slat, lon=$slon)")
            end

            sapflow_lyr = 0.0
            rootuptake_lyr = 0.0
            th_node .= th_node_init

            nsteps = max(imult * iter, 1)
            dt_substep = dt_step / Float64(nsteps)

            for _istep in 1:nsteps
                w_tot_beg = sum(th_node .* v_node) * denh2o

                # On-node quantities for plant nodes
                for i in 1:n_hypool_plant
                    wrfi = wrf_plant(pm_node[i], ft)
                    wkfi = wkf_plant(pm_node[i], ft)
                    psi_node[i] = psi_from_th(wrfi, th_node[i])
                    h_node[i]   = geo * z_node[i] + psi_node[i]
                    ftc_node[i] = ftc_from_psi(wkfi, psi_node[i])
                    dpsi_dtheta_node[i] = dpsidth_from_th(wrfi, th_node[i])
                    dftc_dpsi = dftcdpsi_from_psi(wkfi, psi_node[i])
                    dftc_dtheta_node[i] = dftc_dpsi * dpsi_dtheta_node[i]
                    if i == n_hypool_ag + 2 && no_ftc_radialk
                        ftc_node[i] = 1.0
                        dftc_dtheta_node[i] = 0.0
                    end
                end
                # Rhizosphere shell nodes
                for i in (n_hypool_plant + 1):n_hypool_tot
                    psi_node[i] = psi_from_th(csite_hydr.wrf_soil[ilayer], th_node[i])
                    h_node[i]   = geo * z_node[i] + psi_node[i]
                    ftc_node[i] = ftc_from_psi(csite_hydr.wkf_soil[ilayer], psi_node[i])
                    dpsi_dtheta_node[i] = dpsidth_from_th(csite_hydr.wrf_soil[ilayer], th_node[i])
                    dftc_dpsi = dftcdpsi_from_psi(csite_hydr.wkf_soil[ilayer], psi_node[i])
                    dftc_dtheta_node[i] = dftc_dpsi * dpsi_dtheta_node[i]
                end

                # --- Effective conductances and Taylor terms over each path ---
                # Path 1: leaf -> first stem
                j = 1; i_up = 2; i_dn = 1
                kmax_dn = rootfr_scaler * cohort_hydr.kmax_petiole_to_leaf
                kmax_up = rootfr_scaler * cohort_hydr.kmax_stem_upper[1]
                k_eff[j], a_term[j], b_term[j] = GetImTaylorKAB(kmax_up, kmax_dn,
                    ftc_node[i_up], ftc_node[i_dn], h_node[i_up], h_node[i_dn],
                    dftc_dtheta_node[i_up], dftc_dtheta_node[i_dn],
                    dpsi_dtheta_node[i_up], dpsi_dtheta_node[i_dn])

                # Paths between stem nodes
                for j in 2:(n_hypool_ag - 1)
                    i_up = j + 1; i_dn = j
                    kmax_dn = rootfr_scaler * cohort_hydr.kmax_stem_lower[i_dn - n_hypool_leaf]
                    kmax_up = rootfr_scaler * cohort_hydr.kmax_stem_upper[i_up - n_hypool_leaf]
                    k_eff[j], a_term[j], b_term[j] = GetImTaylorKAB(kmax_up, kmax_dn,
                        ftc_node[i_up], ftc_node[i_dn], h_node[i_up], h_node[i_dn],
                        dftc_dtheta_node[i_up], dftc_dtheta_node[i_dn],
                        dpsi_dtheta_node[i_up], dpsi_dtheta_node[i_dn])
                end

                # Path: lowest stem -> transporting root
                j = n_hypool_ag; i_up = j + 1; i_dn = j
                kmax_dn = rootfr_scaler * cohort_hydr.kmax_stem_lower[n_hypool_stem]
                kmax_up = rootfr_scaler * cohort_hydr.kmax_troot_upper
                k_eff[j], a_term[j], b_term[j] = GetImTaylorKAB(kmax_up, kmax_dn,
                    ftc_node[i_up], ftc_node[i_dn], h_node[i_up], h_node[i_dn],
                    dftc_dtheta_node[i_up], dftc_dtheta_node[i_dn],
                    dpsi_dtheta_node[i_up], dpsi_dtheta_node[i_dn])

                # Path: transporting root -> absorbing root (already parallel, no scaling)
                j = n_hypool_ag + 1; i_up = j + 1; i_dn = j
                kmax_dn = cohort_hydr.kmax_troot_lower[ilayer]
                kmax_up = cohort_hydr.kmax_aroot_upper[ilayer]
                k_eff[j], a_term[j], b_term[j] = GetImTaylorKAB(kmax_up, kmax_dn,
                    ftc_node[i_up], ftc_node[i_dn], h_node[i_up], h_node[i_dn],
                    dftc_dtheta_node[i_up], dftc_dtheta_node[i_dn],
                    dpsi_dtheta_node[i_up], dpsi_dtheta_node[i_dn])

                # Path: absorbing root -> first rhizosphere shell (kmax depends on gradient)
                j = n_hypool_ag + 2; i_up = j + 1; i_dn = j
                if h_node[i_up] > h_node[i_dn]
                    kmax_dn = 1.0 / (1.0 / cohort_hydr.kmax_aroot_lower[ilayer] +
                                     1.0 / cohort_hydr.kmax_aroot_radial_in[ilayer])
                else
                    kmax_dn = 1.0 / (1.0 / cohort_hydr.kmax_aroot_lower[ilayer] +
                                     1.0 / cohort_hydr.kmax_aroot_radial_out[ilayer])
                end
                kmax_up = csite_hydr.kmax_upper_shell[ilayer, 1] * aroot_frac_plant
                k_eff[j], a_term[j], b_term[j] = GetImTaylorKAB(kmax_up, kmax_dn,
                    ftc_node[i_up], ftc_node[i_dn], h_node[i_up], h_node[i_dn],
                    dftc_dtheta_node[i_up], dftc_dtheta_node[i_dn],
                    dpsi_dtheta_node[i_up], dpsi_dtheta_node[i_dn])

                # Paths between rhizosphere shells
                for j in (n_hypool_ag + 3):(n_hypool_tot - 1)
                    i_up = j + 1; i_dn = j
                    ishell_up = i_up - (n_hypool_tot - nshell)
                    ishell_dn = i_dn - (n_hypool_tot - nshell)
                    kmax_dn = csite_hydr.kmax_lower_shell[ilayer, ishell_dn] * aroot_frac_plant
                    kmax_up = csite_hydr.kmax_upper_shell[ilayer, ishell_up] * aroot_frac_plant
                    k_eff[j], a_term[j], b_term[j] = GetImTaylorKAB(kmax_up, kmax_dn,
                        ftc_node[i_up], ftc_node[i_dn], h_node[i_up], h_node[i_dn],
                        dftc_dtheta_node[i_up], dftc_dtheta_node[i_dn],
                        dpsi_dtheta_node[i_up], dpsi_dtheta_node[i_dn])
                end

                # --- Build the tri-diagonal matrix ---
                tris_a[1] = 0.0
                tris_b[1] = a_term[1] - denh2o * v_node[1] / dt_substep
                tris_c[1] = b_term[1]
                tris_r[1] = q_top_eff - k_eff[1] * (h_node[2] - h_node[1])

                for i in 2:(n_hypool_tot - 1)
                    j = i
                    tris_a[i] = -a_term[j-1]
                    tris_b[i] = a_term[j] - b_term[j-1] - denh2o * v_node[i] / dt_substep
                    tris_c[i] = b_term[j]
                    tris_r[i] = -k_eff[j] * (h_node[i+1] - h_node[i]) +
                                 k_eff[j-1] * (h_node[i] - h_node[i-1])
                end

                i = n_hypool_tot; j = n_hypool_tot
                tris_a[i] = -a_term[j-1]
                tris_b[i] = -b_term[j-1] - denh2o * v_node[i] / dt_substep
                tris_c[i] = 0.0
                tris_r[i] = k_eff[j-1] * (h_node[i] - h_node[i-1])

                dth_node, tri_ierr = Hydraulics_Tridiagonal(tris_a, tris_b, tris_c, tris_r)
                if tri_ierr == 1
                    solution_found = false
                    break
                end

                th_node .= th_node .+ dth_node
                w_tot_end = sum(th_node .* v_node) * denh2o
                wb_step_err = (q_top_eff * dt_substep) - (w_tot_beg - w_tot_end)

                if abs(wb_step_err) > max_wb_step_err || any(isnan, dth_node)
                    solution_found = false
                    break
                else
                    solution_found = true
                end

                if trap_neg_wc && any(<(0.0), th_node)
                    solution_found = false
                    break
                end

                wb_err_layer += wb_step_err

                # Diagnostics: sapflow (troot<->lowest stem) and root uptake (shell<->aroot)
                i = n_hypool_ag
                sapflow_lyr += dt_substep * (k_eff[i] * (h_node[i+1] - h_node[i]) +
                                             a_term[i] * dth_node[i] +
                                             b_term[i] * dth_node[i+1])
                i = n_hypool_ag + 2
                rootuptake_lyr += dt_substep * (k_eff[i] * (h_node[i+1] - h_node[i]) +
                                                a_term[i] * dth_node[i] +
                                                b_term[i] * dth_node[i+1])
            end  # substep loop

            iter += 1
        end  # solution iteration

        if abs(wb_err_layer) > max_wb_err
            fates_endrun("EDPlantHydraulics water balance error exceeds threshold " *
                         "$max_wb_err (layer=$ilayer)")
        end

        dth_node = th_node .- th_node_init

        sapflow += sapflow_lyr
        rootuptake[ilayer] = rootuptake_lyr

        if Float64(iter) > cohort_hydr.iterh1 && iter > 1
            cohort_hydr.iterlayer = Float64(ilayer)
        end
        cohort_hydr.iterh1 = max(cohort_hydr.iterh1, Float64(iter))
        cohort_hydr.iterh2 = max(cohort_hydr.iterh2, Float64(nsteps))

        # Update plant compartment water contents
        for k in 1:n_hypool_ag
            cohort_hydr.th_ag[k] += dth_node[k]
        end
        cohort_hydr.th_troot += dth_node[n_hypool_ag + 1]
        cohort_hydr.th_aroot[ilayer] += dth_node[n_hypool_ag + 2]

        dwat_plant += (sum(dth_node[1:n_hypool_ag] .* cohort_hydr.v_ag[1:n_hypool_ag]) +
                       dth_node[n_hypool_ag + 1] * cohort_hydr.v_troot +
                       dth_node[n_hypool_ag + 2] * cohort_hydr.v_aroot_layer[ilayer]) * denh2o

        wb_err_ps += wb_err_layer

        # Accumulate rhizosphere change (shared across cohorts), applied later
        for ishell in 1:nshell
            inode = (n_hypool_tot - nshell) + ishell
            dth_layershell_col[ilayer, ishell] += dth_node[inode] *
                cohort_hydr.l_aroot_layer[ilayer] * cohort.n / csite_hydr.l_aroot_layer[ilayer]
        end
    end  # loop over root layers

    return sapflow, rootuptake, wb_err_ps, dwat_plant
end

# ===========================================================================
# Diagnostics
# ===========================================================================

"""
    SumBetweenDepths(csite_hydr, depth_t, depth_b, array_in) -> depth_sum

Depth-normalized sum of a per-rhizosphere-layer quantity between `depth_t` and `depth_b`
(positive-down). Mirrors the Fortran `SumBetweenDepths`.
"""
function SumBetweenDepths(csite_hydr, depth_t::Real, depth_b::Real, array_in::AbstractVector)
    zi = csite_hydr.zi_rhiz
    dz = csite_hydr.dz_rhiz
    nlevrhiz = csite_hydr.nlevrhiz
    i_rhiz_t = count(<(depth_t), zi .- dz) + 1
    i_rhiz_b = count(<(depth_b), zi)
    depth_sum = 0.0
    if i_rhiz_t > nlevrhiz
        return depth_sum
    end
    if i_rhiz_b >= i_rhiz_t
        depth_sum += sum(@view array_in[i_rhiz_t:i_rhiz_b])
    end
    if i_rhiz_t > 1
        frac = (zi[i_rhiz_t-1] - depth_t) / dz[i_rhiz_t-1]
        depth_sum += frac * array_in[i_rhiz_t-1]
    end
    if i_rhiz_b < nlevrhiz
        # i_rhiz_b == 0 means depth_b falls within the first rhiz layer: the previous
        # interface is the surface (depth 0). Fortran reads zi_rhiz(0) (1-based, out of
        # bounds) here; the physically intended interface depth is 0.0.
        zi_prev = i_rhiz_b >= 1 ? zi[i_rhiz_b] : 0.0
        frac = (depth_b - zi_prev) / dz[i_rhiz_b+1]
        depth_sum += frac * array_in[i_rhiz_b+1]
    end
    depth_sum /= (min(depth_b, zi[nlevrhiz]) - depth_t)
    return depth_sum
end

"""
    BTranForHLMDiagnosticsFromCohortHydr(sites, bc_out)

Aggregate cohort-level `btran` to a balive-weighted patch `btran` for HLM diagnostics.
`bc_out[s].btran_pa` is filled per patch index. Mirrors the Fortran routine.
"""
function BTranForHLMDiagnosticsFromCohortHydr(sites::AbstractVector, bc_out::AbstractVector)
    for s in eachindex(sites)
        ifp = 0
        cpatch = sites[s].oldest_patch
        while cpatch !== nothing
            ifp += 1
            balive_patch = 0.0
            ccohort = cpatch.tallest
            while ccohort !== nothing
                balive_patch += (GetState(ccohort.prt, fnrt_organ, carbon12_element) +
                                 GetState(ccohort.prt, sapw_organ, carbon12_element) +
                                 GetState(ccohort.prt, leaf_organ, carbon12_element)) * ccohort.n
                ccohort = ccohort.shorter
            end
            bc_out[s].btran_pa[ifp] = 0.0
            ccohort = cpatch.tallest
            while ccohort !== nothing
                bc_out[s].btran_pa[ifp] += ccohort.co_hydr.btran *
                    (GetState(ccohort.prt, fnrt_organ, carbon12_element) +
                     GetState(ccohort.prt, sapw_organ, carbon12_element) +
                     GetState(ccohort.prt, leaf_organ, carbon12_element)) * ccohort.n / balive_patch
                ccohort = ccohort.shorter
            end
            cpatch = cpatch.younger
        end
    end
    return nothing
end

# ===========================================================================
# Tier A — init / lifecycle routines
# ===========================================================================

"""
    InitHydrCohort(currentSite, currentCohort)

Allocate the per-cohort plant-hydraulics object (`co_hydr`) and its per-
rhizosphere-layer arrays, sized from the site's `nlevrhiz`. Returns immediately
(no-op) when plant hydraulics is disabled. Mirrors the Fortran `InitHydrCohort`.
"""
function InitHydrCohort(currentSite::ed_site_type, currentCohort::fates_cohort_type)
    hlm_use_planthydro[] == ifalse && return nothing
    ccohort_hydr = ed_cohort_hydr_type()
    currentCohort.co_hydr = ccohort_hydr
    AllocateHydrCohortArrays!(ccohort_hydr, currentSite.si_hydr.nlevrhiz)
    ccohort_hydr.is_newly_recruited = false
    return nothing
end

"""
    DeallocateHydrCohort(currentCohort)

Release the per-cohort plant-hydraulics object (`co_hydr`) and its per-layer
arrays. No-op when plant hydraulics is disabled. Mirrors the Fortran
`DeallocateHydrCohort`.
"""
function DeallocateHydrCohort(currentCohort::fates_cohort_type)
    hlm_use_planthydro[] == ifalse && return nothing
    ccohort_hydr = currentCohort.co_hydr
    if ccohort_hydr !== nothing
        DeallocateHydrCohortArrays!(ccohort_hydr)
    end
    currentCohort.co_hydr = nothing
    return nothing
end

"""
    SavePreviousRhizVolumes(currentSite)

Snapshot the current site-level rhizosphere geometry (absorbing-root length,
nodal shell radius, shell volume) into their `*_init` companions, used to
compute the effect of size change on plant water states. Mirrors the Fortran
`SavePreviousRhizVolumes`.
"""
function SavePreviousRhizVolumes(currentSite::ed_site_type)
    csite_hydr = currentSite.si_hydr
    csite_hydr.l_aroot_layer_init  .= csite_hydr.l_aroot_layer
    csite_hydr.r_node_shell_init   .= csite_hydr.r_node_shell
    csite_hydr.v_shell_init        .= csite_hydr.v_shell
    return nothing
end

"""
    constrain_water_contents(th_uncorr, delta, ft, pm_type) -> th_corr

Clamp an uncorrected volumetric water content into the open interval
`(thr+delta, ths-delta)` defined by the PFT's residual / saturated node water
contents for the given porous-medium type. Mirrors the Fortran (private)
`constrain_water_contents`.
"""
function constrain_water_contents(th_uncorr::Real, delta::Real, ft::Integer, pm_type::Integer)
    pcon = EDPftvarcon_inst[]
    ths = pcon.hydr_thetas_node[ft, pm_type]
    thr = pcon.hydr_resid_node[ft, pm_type]
    return max((thr + delta), min((ths - delta), th_uncorr))
end

"""
    AccumulateMortalityWaterStorage(csite, ccohort, delta_n)

Accumulate the water bound in plants that have just died (number-density loss
`delta_n`, per area) into the site-level dead-vegetation water pool, removing it
from the live stored-vegetation pool. Mirrors the Fortran
`AccumulateMortalityWaterStorage`.
"""
function AccumulateMortalityWaterStorage(csite::ed_site_type, ccohort::fates_cohort_type,
                                         delta_n::Real)
    ccohort_hydr = ccohort.co_hydr
    csite_hydr   = csite.si_hydr
    delta_w = (sum(ccohort_hydr.th_ag .* ccohort_hydr.v_ag) +
               ccohort_hydr.th_troot * ccohort_hydr.v_troot +
               sum(ccohort_hydr.th_aroot .* ccohort_hydr.v_aroot_layer)) *
              dens_fresh_liquid_water * delta_n * AREA_INV

    csite_hydr.h2oveg_dead += delta_w
    csite_hydr.h2oveg      -= delta_w
    return nothing
end

"""
    RecruitWaterStorage(nsites, sites, bc_out)

Diagnose the site-level water bound in newly recruited plants (the
`h2oveg_recruit` pool). No mass is moved here; the actual soil->plant uptake
happens later in `RecruitWUptake`. No-op when plant hydraulics is disabled.
Mirrors the Fortran `RecruitWaterStorage`.
"""
function RecruitWaterStorage(nsites::Integer, sites::AbstractVector, bc_out::AbstractVector)
    hlm_use_planthydro[] == ifalse && return nothing
    for s in 1:nsites
        csite_hydr = sites[s].si_hydr
        csite_hydr.h2oveg_recruit = 0.0
        currentPatch = sites[s].oldest_patch
        while currentPatch !== nothing
            currentCohort = currentPatch.tallest
            while currentCohort !== nothing
                ccohort_hydr = currentCohort.co_hydr
                if ccohort_hydr.is_newly_recruited
                    csite_hydr.h2oveg_recruit += (sum(ccohort_hydr.th_ag .* ccohort_hydr.v_ag) +
                        ccohort_hydr.th_troot * ccohort_hydr.v_troot +
                        sum(ccohort_hydr.th_aroot .* ccohort_hydr.v_aroot_layer)) *
                        dens_fresh_liquid_water * currentCohort.n
                end
                currentCohort = currentCohort.shorter
            end
            currentPatch = currentPatch.younger
        end
        csite_hydr.h2oveg_recruit *= AREA_INV
    end
    return nothing
end

"""
    UpdateSizeDepRhizVolLenCon(currentSite, bc_in)

Update the site's representative-rhizosphere geometry (outer/nodal shell radii
and shell volumes via [`shellGeom!`](@ref)) and the max soil->shell hydraulic
conductances, as the site-aggregated absorbing-root length (summed over all
cohorts and patches) changes. Mirrors the Fortran `UpdateSizeDepRhizVolLenCon`.
"""
function UpdateSizeDepRhizVolLenCon(currentSite::ed_site_type, bc_in)
    csite_hydr = currentSite.si_hydr
    nlevrhiz = csite_hydr.nlevrhiz

    large_kmax_bound = 1.0e4
    k_inner = 1

    # Accumulate cohort-level root length to the site level across patches/cohorts.
    fill!(csite_hydr.l_aroot_layer, 0.0)
    cPatch = currentSite.youngest_patch
    while cPatch !== nothing
        cCohort = cPatch.tallest
        while cCohort !== nothing
            ccohort_hydr = cCohort.co_hydr
            csite_hydr.l_aroot_layer .+= ccohort_hydr.l_aroot_layer .* cCohort.n
            cCohort = cCohort.shorter
        end
        cPatch = cPatch.older
    end

    # Update outer/nodal radii and shell volumes (handles no-root layers inside).
    for j in 1:nlevrhiz
        shellGeom!(csite_hydr.l_aroot_layer[j], csite_hydr.rs1[j], area, csite_hydr.dz_rhiz[j],
                   view(csite_hydr.r_out_shell, j, :), view(csite_hydr.r_node_shell, j, :),
                   view(csite_hydr.v_shell, j, :))
    end

    for j in 1:nlevrhiz
        # hksat [mm s-1] -> [kg s-1 m-1 MPa-1]
        hksat_s = AggBCToRhiz(csite_hydr, bc_in.hksat_sisl, j, bc_in.dz_sisl) *
                  m_per_mm * (1.0 / grav_earth) * pa_per_mpa

        # Only recompute conductances where absorbing-root length changed.
        if (csite_hydr.l_aroot_layer[j] != csite_hydr.l_aroot_layer_init[j]) &&
           (csite_hydr.l_aroot_layer[j] > nearzero)

            if csite_hydr.r_node_shell[j, k_inner] <= csite_hydr.rs1[j]
                csite_hydr.kmax_upper_shell[j, k_inner] = large_kmax_bound
            else
                csite_hydr.kmax_upper_shell[j, k_inner] = 2.0 * pi_const * csite_hydr.l_aroot_layer[j] /
                    log(csite_hydr.r_node_shell[j, k_inner] / csite_hydr.rs1[j]) * hksat_s
            end

            csite_hydr.kmax_lower_shell[j, k_inner] = 2.0 * pi_const * csite_hydr.l_aroot_layer[j] /
                log(csite_hydr.r_out_shell[j, k_inner] / csite_hydr.r_node_shell[j, k_inner]) * hksat_s

            for k in 2:nshell
                csite_hydr.kmax_upper_shell[j, k] = 2.0 * pi_const * csite_hydr.l_aroot_layer[j] /
                    log(csite_hydr.r_node_shell[j, k] / csite_hydr.r_out_shell[j, k - 1]) * hksat_s
                csite_hydr.kmax_lower_shell[j, k] = 2.0 * pi_const * csite_hydr.l_aroot_layer[j] /
                    log(csite_hydr.r_out_shell[j, k] / csite_hydr.r_node_shell[j, k]) * hksat_s
            end
        end
    end
    return nothing
end

"""
    UpdateSizeDepRhizHydProps(currentSite, bc_in)

Thin wrapper: snapshot the previous rhizosphere volumes
([`SavePreviousRhizVolumes`](@ref)) then update the size-dependent rhizosphere
geometry/conductances ([`UpdateSizeDepRhizVolLenCon`](@ref)). Mirrors the
Fortran `UpdateSizeDepRhizHydProps`.
"""
function UpdateSizeDepRhizHydProps(currentSite::ed_site_type, bc_in)
    SavePreviousRhizVolumes(currentSite)
    UpdateSizeDepRhizVolLenCon(currentSite, bc_in)
    return nothing
end

"""
    UpdateSizeDepPlantHydStates(currentSite, ccohort)

Apply water mass conservation to a cohort whose compartment volumes have just
changed (growth/turnover): rescale each compartment's water content by the
volume ratio, clamp via [`constrain_water_contents`](@ref), and book the
resulting tissue-volume-change water error into the site's `h2oveg_growturn_err`.
Mirrors the Fortran `UpdateSizeDepPlantHydStates`.
"""
function UpdateSizeDepPlantHydStates(currentSite::ed_site_type, ccohort::fates_cohort_type)
    ccohort_hydr = ccohort.co_hydr
    ft = ccohort.pft
    csite_hydr = currentSite.si_hydr
    small_theta_num = 1.0e-7

    for k in 1:n_hypool_leaf
        if ccohort_hydr.v_ag[k] > nearzero
            th_uncorr = ccohort_hydr.th_ag[k] * ccohort_hydr.v_ag_init[k] / ccohort_hydr.v_ag[k]
            ccohort_hydr.th_ag[k] = constrain_water_contents(th_uncorr, small_theta_num, ft, leaf_p_media)
        else
            th_uncorr = ccohort_hydr.th_ag[k]
        end
        csite_hydr.h2oveg_growturn_err += dens_fresh_liquid_water * ccohort.n * AREA_INV *
            (ccohort_hydr.th_ag[k] - th_uncorr) * ccohort_hydr.v_ag[k]
    end

    for k in (n_hypool_leaf + 1):n_hypool_ag
        th_uncorr = ccohort_hydr.th_ag[k] * ccohort_hydr.v_ag_init[k] / ccohort_hydr.v_ag[k]
        ccohort_hydr.th_ag[k] = constrain_water_contents(th_uncorr, small_theta_num, ft, stem_p_media)
        csite_hydr.h2oveg_growturn_err += dens_fresh_liquid_water * ccohort.n * AREA_INV *
            (ccohort_hydr.th_ag[k] - th_uncorr) * ccohort_hydr.v_ag[k]
    end

    th_uncorr = ccohort_hydr.th_troot * ccohort_hydr.v_troot_init / ccohort_hydr.v_troot
    ccohort_hydr.th_troot = constrain_water_contents(th_uncorr, small_theta_num, ft, troot_p_media)
    csite_hydr.h2oveg_growturn_err += dens_fresh_liquid_water * ccohort.n * AREA_INV *
        (ccohort_hydr.th_troot - th_uncorr) * ccohort_hydr.v_troot

    for j in 1:csite_hydr.nlevrhiz
        if ccohort_hydr.v_aroot_layer[j] > nearzero
            th_uncorr = ccohort_hydr.th_aroot[j] *
                ccohort_hydr.v_aroot_layer_init[j] / ccohort_hydr.v_aroot_layer[j]
            ccohort_hydr.th_aroot[j] = constrain_water_contents(th_uncorr, small_theta_num, ft, aroot_p_media)
            csite_hydr.h2oveg_growturn_err += dens_fresh_liquid_water * ccohort.n * AREA_INV *
                (ccohort_hydr.th_aroot[j] - th_uncorr) * ccohort_hydr.v_aroot_layer[j]
        end
    end
    return nothing
end

"""
    FuseCohortHydraulics(currentSite, currentCohort, nextCohort, bc_in, newn)

Fuse the donor cohort `nextCohort`'s hydraulics into `currentCohort` (the
post-fusion cohort), conserving total water volume per compartment (the fused
`th` is the number-and-volume-weighted blend), then recompute psi/ftc/btran from
the blended water contents. `newn` is the fused number density. Mirrors the
Fortran `FuseCohortHydraulics`. Assumes `currentCohort`'s size/biomass were
already updated (i.e. called after parteh).
"""
function FuseCohortHydraulics(currentSite::ed_site_type, currentCohort::fates_cohort_type,
                              nextCohort::fates_cohort_type, bc_in, newn::Real)
    csite_hydr   = currentSite.si_hydr
    ccohort_hydr = currentCohort.co_hydr
    ncohort_hydr = nextCohort.co_hydr
    ft = currentCohort.pft

    # Save old volumes (needed for the pre-fusion water volume of each cohort),
    # then update node positions, volumes/lengths, and max conductances.
    SavePreviousCompartmentVolumes!(ccohort_hydr)
    UpdatePlantHydrNodes!(currentCohort, ft, currentCohort.height, csite_hydr)
    UpdatePlantHydrLenVol!(currentCohort, csite_hydr)

    # Conserve total water volume.
    for k in 1:n_hypool_ag
        vol_c1 = currentCohort.n * ccohort_hydr.th_ag[k] * ccohort_hydr.v_ag_init[k]
        vol_c2 = nextCohort.n    * ncohort_hydr.th_ag[k] * ncohort_hydr.v_ag[k]
        ccohort_hydr.th_ag[k] = (vol_c1 + vol_c2) / (ccohort_hydr.v_ag[k] * newn)
    end

    vol_c1 = currentCohort.n * ccohort_hydr.th_troot * ccohort_hydr.v_troot_init
    vol_c2 = nextCohort.n    * ncohort_hydr.th_troot * ncohort_hydr.v_troot
    ccohort_hydr.th_troot = (vol_c1 + vol_c2) / (ccohort_hydr.v_troot * newn)

    for j in 1:csite_hydr.nlevrhiz
        vol_c1 = currentCohort.n * ccohort_hydr.th_aroot[j] * ccohort_hydr.v_aroot_layer_init[j]
        vol_c2 = nextCohort.n    * ncohort_hydr.th_aroot[j] * ncohort_hydr.v_aroot_layer[j]
        ccohort_hydr.th_aroot[j] = (vol_c1 + vol_c2) / (ccohort_hydr.v_aroot_layer[j] * newn)
    end

    ccohort_hydr.supsub_flag = 0

    # Keep the iteration counters of the worse of the two cohorts.
    if ncohort_hydr.iterh1 > ccohort_hydr.iterh1
        ccohort_hydr.iterh1    = ncohort_hydr.iterh1
        ccohort_hydr.iterh2    = ncohort_hydr.iterh2
        ccohort_hydr.iterlayer = ncohort_hydr.iterlayer
    end

    for k in 1:n_hypool_leaf
        ccohort_hydr.psi_ag[k] = psi_from_th(wrf_plant(leaf_p_media, ft), ccohort_hydr.th_ag[k])
        ccohort_hydr.ftc_ag[k] = ftc_from_psi(wkf_plant(leaf_p_media, ft), ccohort_hydr.psi_ag[k])
    end
    for k in (n_hypool_leaf + 1):n_hypool_ag
        ccohort_hydr.psi_ag[k] = psi_from_th(wrf_plant(stem_p_media, ft), ccohort_hydr.th_ag[k])
        ccohort_hydr.ftc_ag[k] = ftc_from_psi(wkf_plant(stem_p_media, ft), ccohort_hydr.psi_ag[k])
    end

    ccohort_hydr.psi_troot = psi_from_th(wrf_plant(troot_p_media, ft), ccohort_hydr.th_troot)
    ccohort_hydr.ftc_troot = ftc_from_psi(wkf_plant(troot_p_media, ft), ccohort_hydr.psi_troot)

    for j in 1:csite_hydr.nlevrhiz
        ccohort_hydr.psi_aroot[j] = psi_from_th(wrf_plant(aroot_p_media, ft), ccohort_hydr.th_aroot[j])
        ccohort_hydr.ftc_aroot[j] = ftc_from_psi(wkf_plant(aroot_p_media, ft), ccohort_hydr.psi_aroot[j])
    end

    ccohort_hydr.btran = ftc_from_psi(wkf_plant(stomata_p_media, ft), ccohort_hydr.psi_ag[1])

    ccohort_hydr.qtop = (currentCohort.n * ccohort_hydr.qtop +
                         nextCohort.n * ncohort_hydr.qtop) / newn
    ccohort_hydr.errh2o = (currentCohort.n * ccohort_hydr.errh2o +
                           nextCohort.n * ncohort_hydr.errh2o) / newn
    return nothing
end

"""
    ConstrainRecruitNumber(csite, ccohort, cpatch, bc_in, mean_temp)

Constrain a new recruit cohort's number density so that the rhizosphere can
supply the water that would be subsumed in the new plant tissues (50% of the
above-residual water per layer is deemed available). When freezing
(`mean_temp <= 273.15`) recruitment is fully suppressed. The carbon/nutrient
mass of any reduced individuals is returned to the patch germination seed pool.
Mirrors the Fortran `ConstrainRecruitNumber`.
"""
function ConstrainRecruitNumber(csite::ed_site_type, ccohort::fates_cohort_type,
                                cpatch::fates_patch_type, bc_in, mean_temp::Real)
    csite_hydr   = csite.si_hydr
    ccohort_hydr = ccohort.co_hydr

    recruitw = (sum(ccohort_hydr.th_ag .* ccohort_hydr.v_ag) +
                ccohort_hydr.th_troot * ccohort_hydr.v_troot +
                sum(ccohort_hydr.th_aroot .* ccohort_hydr.v_aroot_layer)) *
               dens_fresh_liquid_water

    sum_l_aroot = sum(ccohort_hydr.l_aroot_layer)
    for j in 1:csite_hydr.nlevrhiz
        csite_hydr.cohort_recruit_water_layer[j] = recruitw * ccohort_hydr.l_aroot_layer[j] / sum_l_aroot
    end

    for j in 1:csite_hydr.nlevrhiz
        watres_local = th_from_psi(csite_hydr.wrf_soil[j],
            bc_in.smpmin_si * dens_fresh_liquid_water * grav_earth * m_per_mm * mpa_per_pa)
        total_water     = sum(view(csite_hydr.v_shell, j, :) .* view(csite_hydr.h2osoi_liqvol_shell, j, :))
        total_water_min = sum(view(csite_hydr.v_shell, j, :) .* watres_local)
        # Only 50% is deemed available for recruit water.
        csite_hydr.recruit_water_avail_layer[j] = 0.5 * max(0.0, total_water - total_water_min)
    end

    nmin = 1.0e36
    for j in 1:csite_hydr.nlevrhiz
        if csite_hydr.cohort_recruit_water_layer[j] > nearzero
            n = csite_hydr.recruit_water_avail_layer[j] / csite_hydr.cohort_recruit_water_layer[j]
            nmin = min(n, nmin)
        end
    end

    # Prevent recruitment when temperatures are freezing or below.
    if mean_temp <= 273.15
        nmin = 0.0
    end

    # Reduce number density (and return carbon/nutrient mass to germination pool)
    # if water-limited recruitment is below the carbon-allowed number density.
    if nmin < ccohort.n
        for el in 1:num_elements[]
            element_id = element_list[el]
            leaf_m   = GetState(ccohort.prt, leaf_organ, element_id)
            store_m  = GetState(ccohort.prt, store_organ, element_id)
            sapw_m   = GetState(ccohort.prt, sapw_organ, element_id)
            fnrt_m   = GetState(ccohort.prt, fnrt_organ, element_id)
            struct_m = GetState(ccohort.prt, struct_organ, element_id)
            repro_m  = GetState(ccohort.prt, repro_organ, element_id)

            cpatch.litter[el].seed_germ[ccohort.pft] += (ccohort.n - nmin) / cpatch.area *
                (leaf_m + store_m + sapw_m + fnrt_m + struct_m + repro_m)
        end
        ccohort.n = nmin
    end
    return nothing
end

"""
    InitHydrSites(sites, bc_in)

Allocate each site's `si_hydr` object, choose the rhizosphere-layer aggregation
(the FATES default `combine12`: aggregate the top two soil layers, the rest 1:1),
and size/initialize the site hydraulics arrays via [`InitHydrSite!`](@ref) with
the soil->rhizosphere layer map and depths. No-op when plant hydraulics is
disabled. Mirrors the Fortran `InitHydrSites`.
"""
function InitHydrSites(sites::AbstractVector, bc_in::AbstractVector)
    hlm_use_planthydro[] == ifalse && return nothing

    rhizlayer_aggmeth_none      = 1
    rhizlayer_aggmeth_combine12 = 2
    rhizlayer_aggmeth_balN      = 3

    solver = EDParams[].hydr_solver

    nsites = length(sites)
    for s in 1:nsites
        csite_hydr = ed_site_hydr_type()
        sites[s].si_hydr = csite_hydr
        if bc_in[s].nlevsoil > nlevsoi_hyd_max
            fates_endrun("FATES-hydro: host soil has $(bc_in[s].nlevsoil) layers, " *
                         "exceeding nlevsoi_hyd_max = $nlevsoi_hyd_max")
        end

        aggmeth = rhizlayer_aggmeth_combine12
        aggN    = 10

        # NOTE: bc_in.zi_sisl carries a zero index for the surface, so Fortran
        # zi_sisl(j) (bottom of soil layer j) is zi_sisl[j+1] here.
        if aggmeth == rhizlayer_aggmeth_none
            csite_hydr.nlevrhiz = bc_in[s].nlevsoil
            InitHydrSite!(csite_hydr, numpft[], nlevsclass[], solver, bc_in[s].nlevsoil)
            for j in 1:csite_hydr.nlevrhiz
                csite_hydr.map_r2s[j, 1] = j
                csite_hydr.map_r2s[j, 2] = j
                csite_hydr.zi_rhiz[j]    = bc_in[s].zi_sisl[j + 1]
                csite_hydr.dz_rhiz[j]    = bc_in[s].dz_sisl[j]
            end

        elseif aggmeth == rhizlayer_aggmeth_combine12
            csite_hydr.nlevrhiz = max(1, bc_in[s].nlevsoil - 1)
            InitHydrSite!(csite_hydr, numpft[], nlevsclass[], solver, bc_in[s].nlevsoil)

            csite_hydr.map_r2s[1, 1] = 1
            j_bc = min(2, bc_in[s].nlevsoil)   # protects the single-soil-layer case
            csite_hydr.map_r2s[1, 2] = j_bc
            csite_hydr.zi_rhiz[1]    = bc_in[s].zi_sisl[j_bc + 1]
            csite_hydr.dz_rhiz[1]    = sum(@view bc_in[s].dz_sisl[1:j_bc])

            for j in 2:csite_hydr.nlevrhiz
                csite_hydr.map_r2s[j, 1] = j + 1
                csite_hydr.map_r2s[j, 2] = j + 1
                csite_hydr.zi_rhiz[j]    = bc_in[s].zi_sisl[(j + 1) + 1]
                csite_hydr.dz_rhiz[j]    = bc_in[s].dz_sisl[j + 1]
            end

        elseif aggmeth == rhizlayer_aggmeth_balN
            csite_hydr.nlevrhiz = min(aggN, bc_in[s].nlevsoil)
            InitHydrSite!(csite_hydr, numpft[], nlevsclass[], solver, bc_in[s].nlevsoil)

            ntoagg = Int(ceil(bc_in[s].nlevsoil / csite_hydr.nlevrhiz - nearzero))
            if ntoagg < 1
                fates_endrun("rhizlayer_aggmeth_balN: bad soil-layers-per-rhiz estimate $ntoagg")
            end
            ns_per_rhiz = fill(ntoagg, csite_hydr.nlevrhiz)
            while sum(ns_per_rhiz) > bc_in[s].nlevsoil
                for j in csite_hydr.nlevrhiz:-1:1
                    ns_per_rhiz[j] -= 1
                    sum(ns_per_rhiz) <= bc_in[s].nlevsoil && break
                    if ns_per_rhiz[j] == 0
                        fates_endrun("rhizlayer_aggmeth_balN produced a 0-soil-layer rhiz layer")
                    end
                end
            end

            csite_hydr.map_r2s[1, 1] = 1
            for j in 1:(csite_hydr.nlevrhiz - 1)
                j_t = csite_hydr.map_r2s[j, 1]
                j_b = j_t + ns_per_rhiz[j] - 1
                csite_hydr.map_r2s[j, 2]     = j_b
                csite_hydr.map_r2s[j + 1, 1] = j_b + 1
                csite_hydr.zi_rhiz[j]        = bc_in[s].zi_sisl[j_b + 1]
                csite_hydr.dz_rhiz[j]        = sum(@view bc_in[s].dz_sisl[j_t:j_b])
            end
            j_t = csite_hydr.map_r2s[csite_hydr.nlevrhiz, 1]
            j_b = j_t + ns_per_rhiz[csite_hydr.nlevrhiz] - 1
            csite_hydr.map_r2s[csite_hydr.nlevrhiz, 2] = j_b
            csite_hydr.zi_rhiz[csite_hydr.nlevrhiz]    = bc_in[s].zi_sisl[j_b + 1]
            csite_hydr.dz_rhiz[csite_hydr.nlevrhiz]    = sum(@view bc_in[s].dz_sisl[j_t:j_b])
        else
            fates_endrun("undefined rhizosphere layer aggregation method: $aggmeth")
        end
    end
    return nothing
end

"""
    HydrSiteColdStart(sites, bc_in)

Cold-start each site's hydraulics: set the rhizosphere shell liquid water from
the soil effective porosity / liquid water, zero the site absorbing-root length,
and construct the per-soil-layer soil water-retention (WRF) and conductance (WKF)
functions (Campbell/Clapp-Hornberger by default) with parameters aggregated to
each rhizosphere layer. Mirrors the Fortran `HydrSiteColdStart`.
"""
function HydrSiteColdStart(sites::AbstractVector, bc_in::AbstractVector)
    nsites = length(sites)
    for s in 1:nsites
        csite_hydr = sites[s].si_hydr
        nlevrhiz = csite_hydr.nlevrhiz

        for j in 1:nlevrhiz
            j_t = csite_hydr.map_r2s[j, 1]
            j_b = csite_hydr.map_r2s[j, 2]
            eff_por = AggBCToRhiz(csite_hydr, bc_in[s].eff_porosity_sl, j, bc_in[s].dz_sisl)
            # [kg/m2] / ([m] * [kg/m3]) = [m3/m3]
            h2osoi_liqvol = min(eff_por,
                sum(@view bc_in[s].h2o_liq_sisl[j_t:j_b]) / (csite_hydr.dz_rhiz[j] * dens_fresh_liquid_water))
            csite_hydr.h2osoi_liqvol_shell[j, 1:nshell] .= h2osoi_liqvol
        end

        csite_hydr.l_aroot_layer[1:nlevrhiz] .= 0.0

        # --- Soil Water Retention Functions (WRFs) ---
        if soil_wrf_type == van_genuchten_type
            for j in 1:nlevrhiz
                wrf = wrf_type_vg()
                set_wrf_param!(wrf, [alpha_vg, psd_vg, m_vg, th_sat_vg, th_res_vg])
                csite_hydr.wrf_soil[j] = wrf
            end
        elseif soil_wrf_type == campbell_type
            for j in 1:nlevrhiz
                wrf = wrf_type_cch()
                watsat = AggBCToRhiz(csite_hydr, bc_in[s].watsat_sisl, j, bc_in[s].dz_sisl)
                sucsat = AggBCToRhiz(csite_hydr, bc_in[s].sucsat_sisl, j, bc_in[s].dz_sisl)
                bsw    = AggBCToRhiz(csite_hydr, bc_in[s].bsw_sisl, j, bc_in[s].dz_sisl)
                set_wrf_param!(wrf, [watsat,
                    (-1.0) * sucsat * dens_fresh_liquid_water * grav_earth * mpa_per_pa * m_per_mm, bsw])
                csite_hydr.wrf_soil[j] = wrf
            end
        elseif soil_wrf_type == smooth1_campbell_type || soil_wrf_type == smooth2_campbell_type
            smooth_order = soil_wrf_type == smooth1_campbell_type ? 1.0 : 2.0
            for j in 1:nlevrhiz
                wrf = wrf_type_smooth_cch()
                watsat = AggBCToRhiz(csite_hydr, bc_in[s].watsat_sisl, j, bc_in[s].dz_sisl)
                sucsat = AggBCToRhiz(csite_hydr, bc_in[s].sucsat_sisl, j, bc_in[s].dz_sisl)
                bsw    = AggBCToRhiz(csite_hydr, bc_in[s].bsw_sisl, j, bc_in[s].dz_sisl)
                set_wrf_param!(wrf, [watsat,
                    (-1.0) * sucsat * dens_fresh_liquid_water * grav_earth * mpa_per_pa * m_per_mm, bsw, smooth_order])
                csite_hydr.wrf_soil[j] = wrf
            end
        elseif soil_wrf_type == tfs_type
            fates_endrun("TFS water retention curves not available for soil")
        end

        # --- Soil Water Conductance (K) Functions (WKFs) ---
        if soil_wkf_type == van_genuchten_type
            for j in 1:nlevrhiz
                wkf = wkf_type_vg()
                set_wkf_param!(wkf, [alpha_vg, psd_vg, m_vg, th_sat_vg, th_res_vg, soil_tort_vg])
                csite_hydr.wkf_soil[j] = wkf
            end
        elseif soil_wkf_type == campbell_type
            for j in 1:nlevrhiz
                wkf = wkf_type_cch()
                watsat = AggBCToRhiz(csite_hydr, bc_in[s].watsat_sisl, j, bc_in[s].dz_sisl)
                sucsat = AggBCToRhiz(csite_hydr, bc_in[s].sucsat_sisl, j, bc_in[s].dz_sisl)
                bsw    = AggBCToRhiz(csite_hydr, bc_in[s].bsw_sisl, j, bc_in[s].dz_sisl)
                set_wkf_param!(wkf, [watsat,
                    (-1.0) * sucsat * dens_fresh_liquid_water * grav_earth * mpa_per_pa * m_per_mm, bsw])
                csite_hydr.wkf_soil[j] = wkf
            end
        elseif soil_wkf_type == smooth1_campbell_type || soil_wkf_type == smooth2_campbell_type
            smooth_order = soil_wkf_type == smooth1_campbell_type ? 1.0 : 2.0
            for j in 1:nlevrhiz
                wkf = wkf_type_smooth_cch()
                watsat = AggBCToRhiz(csite_hydr, bc_in[s].watsat_sisl, j, bc_in[s].dz_sisl)
                sucsat = AggBCToRhiz(csite_hydr, bc_in[s].sucsat_sisl, j, bc_in[s].dz_sisl)
                bsw    = AggBCToRhiz(csite_hydr, bc_in[s].bsw_sisl, j, bc_in[s].dz_sisl)
                set_wkf_param!(wkf, [watsat,
                    (-1.0) * sucsat * dens_fresh_liquid_water * grav_earth * mpa_per_pa * m_per_mm, bsw, smooth_order])
                csite_hydr.wkf_soil[j] = wkf
            end
        elseif soil_wkf_type == tfs_type
            fates_endrun("TFS conductance not used in soil")
        end

        # The FTC functions need psi_min -> point each WKF at its matching WRF.
        for j in 1:nlevrhiz
            csite_hydr.wkf_soil[j].wrf = csite_hydr.wrf_soil[j]
        end
    end
    return nothing
end

# ===========================================================================
# Site stored-vegetation water bookkeeping
# ===========================================================================

"""
    UpdateH2OVeg!(csite, bc_out; prev_site_h2o=nothing, icall=-1)

Re-assess the liquid water bound in the plants of a site (summing each cohort's
compartment water content × volume × number density), store it in
`csite.si_hydr.h2oveg` [kg/m2], and report the total stored live+dead vegetation
water (less the growth/turnover and hydrodynamics error pools) into
`bc_out.plant_stored_h2o_si`. Newly recruited cohorts are excluded (their water
is accounted via the recruitment-uptake path). When `prev_site_h2o` is given, a
conservation check (threshold `error_thresh`) is performed. No-op when plant
hydraulics is disabled (still zeroes `plant_stored_h2o_si`). Mirrors the Fortran
`UpdateH2OVeg`.
"""
function UpdateH2OVeg!(csite::ed_site_type, bc_out;
                       prev_site_h2o::Union{Real,Nothing}=nothing, icall::Integer=-1)
    bc_out.plant_stored_h2o_si = 0.0
    hlm_use_planthydro[] == ifalse && return nothing

    csite_hydr = csite.si_hydr
    csite_hydr.h2oveg = 0.0
    currentPatch = csite.oldest_patch
    while currentPatch !== nothing
        currentCohort = currentPatch.tallest
        while currentCohort !== nothing
            ccohort_hydr = currentCohort.co_hydr
            if !ccohort_hydr.is_newly_recruited
                csite_hydr.h2oveg += (sum(ccohort_hydr.th_ag .* ccohort_hydr.v_ag) +
                    ccohort_hydr.th_troot * ccohort_hydr.v_troot +
                    sum(ccohort_hydr.th_aroot .* ccohort_hydr.v_aroot_layer)) *
                    dens_fresh_liquid_water * currentCohort.n
            end
            currentCohort = currentCohort.shorter
        end
        currentPatch = currentPatch.younger
    end

    # Convert from kg/site to kg/m2.
    csite_hydr.h2oveg *= AREA_INV

    bc_out.plant_stored_h2o_si = csite_hydr.h2oveg + csite_hydr.h2oveg_dead -
                                 csite_hydr.h2oveg_growturn_err -
                                 csite_hydr.h2oveg_hydro_err

    if prev_site_h2o !== nothing
        if abs(bc_out.plant_stored_h2o_si - prev_site_h2o) > error_thresh
            fates_endrun("Total FATES site level water was not conserved during a check " *
                         "where it was supposed to be conserved (call index $icall): " *
                         "old=$(prev_site_h2o) new=$(bc_out.plant_stored_h2o_si)")
        end
    end
    return nothing
end

# ===========================================================================
# Tier B — transpiration / uptake solve
# ===========================================================================

"""
    SetMaxCondConnections(csite_hydr, cohort_hydr, h_node) -> (kmax_dn, kmax_up)

Set the maximum conductance on the downstream (towards-atmosphere) and upstream
(towards-soil) side of every node-to-node connection of the leaf->stem->troot->
aroot->rhizosphere continuum. The aroot<->soil connection's downstream kmax
depends on the sign of the (matric+height) potential gradient `h_node` (radial-in
vs radial-out root surface conductance). Used by the 2D solvers; ported as a
faithful helper. Returns the two conductance vectors (length `num_connections`).
Mirrors the Fortran `SetMaxCondConnections` (preserving `do_parallel_stem`'s
single-conduit connection layout via the rhizosphere-shell loop).
"""
function SetMaxCondConnections(csite_hydr, cohort_hydr, h_node::AbstractVector)
    # The full multi-layer (2D) connection count this routine fills: leaf->stem,
    # stem->stem, lowest-stem->troot, then (aroot + shells) per rhizosphere layer.
    # (The 1D solver stores a per-layer `num_connections`; SetMaxCondConnections is a
    # 2D-only helper that spans all layers, so size from the explicit formula.)
    num_cnxs = n_hypool_leaf + n_hypool_stem + n_hypool_troot - 1 +
               (n_hypool_aroot + nshell) * csite_hydr.nlevrhiz
    kmax_dn = fill(fates_unset_r8, num_cnxs)
    kmax_up = fill(fates_unset_r8, num_cnxs)

    # Leaf -> stem (single leaf layer).
    icnx = 1
    kmax_dn[icnx] = cohort_hydr.kmax_petiole_to_leaf
    kmax_up[icnx] = cohort_hydr.kmax_stem_upper[1]

    # Stem -> stem.
    for istem in 1:(n_hypool_stem - 1)
        icnx += 1
        kmax_dn[icnx] = cohort_hydr.kmax_stem_lower[istem]
        kmax_up[icnx] = cohort_hydr.kmax_stem_upper[istem + 1]
    end

    # Lowest stem -> transporting root.
    icnx += 1
    kmax_dn[icnx] = cohort_hydr.kmax_stem_lower[n_hypool_stem]
    kmax_up[icnx] = cohort_hydr.kmax_troot_upper

    # Transporting root -> absorbing roots -> rhizosphere shells, per soil layer.
    inode = n_hypool_ag
    for j in 1:csite_hydr.nlevrhiz
        aroot_frac_plant = cohort_hydr.l_aroot_layer[j] / csite_hydr.l_aroot_layer[j]
        for k in 1:(n_hypool_aroot + nshell)
            icnx += 1
            inode += 1
            if k == 1            # troot -> aroot
                kmax_dn[icnx] = cohort_hydr.kmax_troot_lower[j]
                kmax_up[icnx] = cohort_hydr.kmax_aroot_upper[j]
            elseif k == 2        # aroot -> soil (gradient-dependent)
                if h_node[inode] < h_node[inode + 1]
                    kmax_dn[icnx] = 1.0 / (1.0 / cohort_hydr.kmax_aroot_lower[j] +
                                           1.0 / cohort_hydr.kmax_aroot_radial_in[j])
                else
                    kmax_dn[icnx] = 1.0 / (1.0 / cohort_hydr.kmax_aroot_lower[j] +
                                           1.0 / cohort_hydr.kmax_aroot_radial_out[j])
                end
                kmax_up[icnx] = csite_hydr.kmax_upper_shell[j, 1] * aroot_frac_plant
            else                 # soil -> soil
                kmax_dn[icnx] = csite_hydr.kmax_lower_shell[j, k - 2] * aroot_frac_plant
                kmax_up[icnx] = csite_hydr.kmax_upper_shell[j, k - 1] * aroot_frac_plant
            end
        end
    end
    return kmax_dn, kmax_up
end

"""
    Report1DError(cohort, csite_hydr, ilayer, z_node, v_node, th_node,
                  q_top_eff, dt_step, w_tot_beg, w_tot_end, rootfr_scaler,
                  aroot_frac_plant, err_code, err_arr, slat, slon, recruitflag)
        -> String

Build a faithful diagnostic record describing the failed initial condition of a
1D solve (node water contents, potentials, total-potentials, per-compartment
water, kmaxes) for a cohort/soil-layer that could not find a stable solution.
Returns the assembled report as a String (also emitted via `fates_log`). The
Fortran routine writes to the model log and the caller decides whether to
`endrun`; here the report is returned so callers (e.g. tests) can inspect it.
Mirrors the Fortran `Report1DError`.
"""
function Report1DError(cohort, csite_hydr, ilayer::Integer,
                       z_node::AbstractVector, v_node::AbstractVector,
                       th_node::AbstractVector, q_top_eff::Real, dt_step::Real,
                       w_tot_beg::Real, w_tot_end::Real, rootfr_scaler::Real,
                       aroot_frac_plant::Real, err_code::Integer,
                       err_arr::AbstractVector, slat::Real, slon::Real,
                       recruitflag::Bool)
    cohort_hydr = cohort.co_hydr
    ft = cohort.pft
    geo = mpa_per_pa * dens_fresh_liquid_water * grav_earth

    nn = length(z_node)
    psi_node = zeros(Float64, nn)
    h_node   = zeros(Float64, nn)
    for i in 1:n_hypool_plant
        psi_node[i] = psi_from_th(wrf_plant(csite_hydr.pm_node[i], ft), th_node[i])
        h_node[i]   = geo * z_node[i] + psi_node[i]
    end
    for i in (n_hypool_plant + 1):n_hypool_tot
        psi_node[i] = psi_from_th(csite_hydr.wrf_soil[ilayer], th_node[i])
        h_node[i]   = geo * z_node[i] + psi_node[i]
    end

    leaf_water  = sum(cohort_hydr.th_ag[1:n_hypool_leaf] .* cohort_hydr.v_ag[1:n_hypool_leaf]) *
                  dens_fresh_liquid_water
    stem_water  = sum(cohort_hydr.th_ag[(n_hypool_leaf + 1):n_hypool_ag] .*
                      cohort_hydr.v_ag[(n_hypool_leaf + 1):n_hypool_ag]) * dens_fresh_liquid_water
    troot_water = cohort_hydr.th_troot * cohort_hydr.v_troot * dens_fresh_liquid_water
    aroot_water = sum(cohort_hydr.th_aroot .* cohort_hydr.v_aroot_layer) * dens_fresh_liquid_water

    io = IOBuffer()
    println(io, "Could not find a stable solution for hydro 1D solve")
    println(io, "error code: ", err_code)
    println(io, "error diag: ", err_arr)
    println(io, "lat: ", slat, " lon: ", slon)
    println(io, "is recruitment: ", recruitflag)
    println(io, "layer: ", ilayer)
    println(io, "wb_step_err = ", (q_top_eff * dt_step) - (w_tot_beg - w_tot_end))
    println(io, "q_top_eff*dt_step = ", q_top_eff * dt_step)
    println(io, "w_tot_beg = ", w_tot_beg)
    println(io, "w_tot_end = ", w_tot_end)
    println(io, "leaf water: ", leaf_water, " kg/plant")
    println(io, "stem_water: ", stem_water, " kg/plant")
    println(io, "troot_water: ", troot_water)
    println(io, "aroot_water: ", aroot_water)
    println(io, "LWP: ", cohort_hydr.psi_ag[1])
    println(io, "dbh: ", cohort.dbh)
    println(io, "pft: ", cohort.pft)
    println(io, "z nodes: ", z_node)
    println(io, "psi_z: ", h_node .- psi_node)
    for i in 1:n_hypool_leaf
        println(io, "leaf node ", i, " ", v_node[i], " ", th_node[i], " ", h_node[i], " ", psi_node[i])
    end
    for i in 1:n_hypool_stem
        k = i + n_hypool_leaf
        println(io, "stem node ", k, " ", v_node[k], " ", th_node[k], " ", h_node[k], " ", psi_node[k])
    end
    k = n_hypool_leaf + n_hypool_stem + 1
    println(io, "troot node: ", k, " ", v_node[k], " ", th_node[k], " ", h_node[k])
    k = n_hypool_leaf + n_hypool_stem + 2
    println(io, "aroot node: ", k, " ", v_node[k], " ", th_node[k], " ", h_node[k])
    for i in 1:nshell
        k = n_hypool_leaf + n_hypool_stem + 2 + i
        println(io, "rhizo shell k: ", k, " ", v_node[k], " ", th_node[k], " ", h_node[k])
    end
    println(io, "kmax_aroot_radial_out: ", cohort_hydr.kmax_aroot_radial_out[ilayer])
    println(io, "aroot_frac_plant: ", aroot_frac_plant, " ",
            cohort_hydr.l_aroot_layer[ilayer], " ", csite_hydr.l_aroot_layer[ilayer])
    println(io, "tree lai: ", cohort.treelai, " m2/m2 crown")
    msg = String(take!(io))
    @warn msg
    return msg
end

"""
    RecruitWUptake(nsites, sites, bc_in, dtime) -> recruitflag

For all newly recruited cohorts across every site, compute the water demand that
must be drawn from the rhizosphere as the plants pop into existence and tally it
into `csite_hydr.recruit_w_uptake` [kg/m2/s] per rhizosphere layer (apportioned by
absorbing-root length), zeroing the `is_newly_recruited` flag afterwards. A
per-site mass-balance check (threshold `1e-10`) is performed. Returns `true` if
any cohort was newly recruited. Mirrors the Fortran `RecruitWUptake`.
"""
function RecruitWUptake(nsites::Integer, sites::AbstractVector, bc_in::AbstractVector,
                        dtime::Real)
    recruitflag = false
    for s in 1:nsites
        csite_hydr = sites[s].si_hydr
        fill!(csite_hydr.recruit_w_uptake, 0.0)
        currentPatch = sites[s].oldest_patch
        recruitw_total = 0.0
        while currentPatch !== nothing
            currentCohort = currentPatch.tallest
            while currentCohort !== nothing
                ccohort_hydr = currentCohort.co_hydr
                if ccohort_hydr.is_newly_recruited
                    recruitflag = true
                    recruitw = (sum(ccohort_hydr.th_ag .* ccohort_hydr.v_ag) +
                        ccohort_hydr.th_troot * ccohort_hydr.v_troot +
                        sum(ccohort_hydr.th_aroot .* ccohort_hydr.v_aroot_layer)) *
                        dens_fresh_liquid_water * currentCohort.n * AREA_INV / dtime
                    recruitw_total += recruitw
                    sum_l_aroot = sum(ccohort_hydr.l_aroot_layer)
                    for j in 1:csite_hydr.nlevrhiz
                        rootfr = ccohort_hydr.l_aroot_layer[j] / sum_l_aroot
                        csite_hydr.recruit_w_uptake[j] += recruitw * rootfr
                    end
                    ccohort_hydr.is_newly_recruited = false
                end
                currentCohort = currentCohort.shorter
            end
            currentPatch = currentPatch.younger
        end

        # Mass balance check.
        sumrw_uptake = sum(@view csite_hydr.recruit_w_uptake[1:csite_hydr.nlevrhiz])
        err = recruitw_total - sumrw_uptake
        if abs(err) > 1.0e-10
            for j in 1:csite_hydr.nlevrhiz
                csite_hydr.recruit_w_uptake[j] += err * csite_hydr.recruit_w_uptake[j] / sumrw_uptake
            end
            fates_endrun("math check on recruit water failed with err= $err " *
                         "$sumrw_uptake $recruitw_total")
        end
    end
    return recruitflag
end

"""
    FillDrainRhizShells(nsites, sites, bc_in, bc_out)

Heuristically redistribute the change in each rhizosphere layer's mean liquid
water (from infiltration/drainage/vertical movement on the host soil grid since
the last hydraulics call) across that layer's radial shells: water draining out
comes preferentially from the wettest shells, water filling in goes preferentially
to the driest shells. The shells are ordered by water content and filled/drained
up/down to the next shell's content until the layer's total change is accounted.
Conserves the per-layer water mass (a balance check, threshold `1e-9`, runs per
layer). Mirrors the Fortran `FillDrainRhizShells`.
"""
function FillDrainRhizShells(nsites::Integer, sites::AbstractVector,
                             bc_in::AbstractVector, bc_out::AbstractVector)
    for s in 1:nsites
        csite_hydr = sites[s].si_hydr

        # If there are just no plants in this site, don't shuffle water.
        sum(@view csite_hydr.l_aroot_layer[1:csite_hydr.nlevrhiz]) <= nearzero && continue

        for j in 1:csite_hydr.nlevrhiz
            j_t = csite_hydr.map_r2s[j, 1]   # top soil layer matching rhiz layer
            j_b = csite_hydr.map_r2s[j, 2]   # bottom soil layer matching rhiz layer

            csite_hydr.l_aroot_layer[j] <= nearzero && continue

            cumShellH2O = sum(view(csite_hydr.h2osoi_liqvol_shell, j, :) .*
                              view(csite_hydr.v_shell, j, :)) * dens_fresh_liquid_water * AREA_INV

            dwat_kgm2 = sum(@view bc_in[s].h2o_liq_sisl[j_t:j_b]) - cumShellH2O   # [kg/m2]
            dwat_kg   = dwat_kgm2 * area

            # Order shells by increasing (fill) or decreasing (drain) water content.
            ordered = collect(1:nshell)
            if nshell > 1
                for k in (nshell - 1):-1:1
                    for kk in 1:k
                        if csite_hydr.h2osoi_liqvol_shell[j, ordered[kk]] >
                           csite_hydr.h2osoi_liqvol_shell[j, ordered[kk + 1]]
                            if dwat_kg > 0.0   # order increasing
                                ordered[kk], ordered[kk + 1] = ordered[kk + 1], ordered[kk]
                            end
                        else
                            if dwat_kg < 0.0   # order decreasing
                                ordered[kk], ordered[kk + 1] = ordered[kk + 1], ordered[kk]
                            end
                        end
                    end
                end
            end

            # Fill (dwat>0) or drain (dwat<0) shells up/down to the next shell's content.
            k = 1
            while (dwat_kg != 0.0) && (k < nshell)
                thdiff = csite_hydr.h2osoi_liqvol_shell[j, ordered[k + 1]] -
                         csite_hydr.h2osoi_liqvol_shell[j, ordered[k]]
                v_cum  = sum(csite_hydr.v_shell[j, ordered[kk]] for kk in 1:k)
                wdiff  = thdiff * v_cum * dens_fresh_liquid_water   # [kg] over shells ordered[1:k]
                if abs(dwat_kg) >= abs(wdiff)
                    for kk in 1:k
                        csite_hydr.h2osoi_liqvol_shell[j, ordered[kk]] =
                            csite_hydr.h2osoi_liqvol_shell[j, ordered[k + 1]]
                    end
                    dwat_kg -= wdiff
                else
                    for kk in 1:k
                        csite_hydr.h2osoi_liqvol_shell[j, ordered[kk]] += dwat_kg / dens_fresh_liquid_water / v_cum
                    end
                    dwat_kg = 0.0
                end
                k += 1
            end

            if dwat_kg != 0.0
                v_cum  = sum(csite_hydr.v_shell[j, ordered[kk]] for kk in 1:nshell)
                thdiff = dwat_kg / v_cum / dens_fresh_liquid_water
                for k in nshell:-1:1
                    csite_hydr.h2osoi_liqvol_shell[j, k] += thdiff
                end
            end

            # Per-layer water balance check [kg/m2].
            h2osoi_liq_shell = sum(view(csite_hydr.h2osoi_liqvol_shell, j, :) .*
                                   view(csite_hydr.v_shell, j, :)) * dens_fresh_liquid_water
            errh2o = h2osoi_liq_shell * AREA_INV - sum(@view bc_in[s].h2o_liq_sisl[j_t:j_b])
            if abs(errh2o) > 1.0e-9
                fates_endrun("water balance error in FillDrainRhizShells: errh2o= $errh2o [kg/m2]")
            end
        end
    end
    return nothing
end

"""
    Hydraulics_BC(nsites, sites, bc_in, bc_out, dtime)

The core per-cohort transpiration/uptake orchestrator (default `1DTaylor` solver
path). For every cohort of every vegetated patch it: partitions the HLM-supplied
patch transpiration (`bc_in.qflx_transp_pa`) by leaf-area-weighted stomatal
conductance, runs `OrderLayersForSolve1D` + `ImTaylorSolve1D` to advance the
plant-soil water-flow network, accumulates the site-level stored-plant-water
change / error, the sapflow + per-depth root-uptake diagnostics, and the
rhizosphere shell-water change shared across cohorts. After the cohort loop it
applies the rhizosphere change, disaggregates the per-rhiz-layer root uptake onto
the host soil grid (conductance-weighted) into `bc_out.qflx_soil2root_sisl`, and
runs the site plant-water + total-water balance checks. Crucially it populates
`co_hydr.ftc_ag/ftc_troot/ftc_aroot` (via `UpdatePlantPsiFTCFromTheta!`),
`co_hydr.btran`, and `leaf_psi` (`psi_ag[1]`) — the fields the mortality and
photosynthesis branches consume. The 2D solver branches raise (not ported).
Mirrors the Fortran `hydraulics_bc`.
"""
function Hydraulics_BC(nsites::Integer, sites::AbstractVector,
                       bc_in::AbstractVector, bc_out::AbstractVector, dtime::Real)
    soilz_disagg = 0
    soilk_disagg = 1
    rootflux_disagg = soilk_disagg

    denh2o = dens_fresh_liquid_water
    solver = EDParams[].hydr_solver

    ordered_init = collect(1:nlevsoi_hyd_max)
    kbg_layer = zeros(Float64, nlevsoi_hyd_max)

    # For newly recruited cohorts, add their water demand to recruit_w_uptake.
    recruitflag = RecruitWUptake(nsites, sites, bc_in, dtime)

    # Update stored veg water after incorporating newly recruited cohorts.
    if recruitflag
        for s in 1:nsites
            UpdateH2OVeg!(sites[s], bc_out[s])
        end
    end

    for s in 1:nsites
        csite_hydr = sites[s].si_hydr

        fill!(bc_out[s].qflx_soil2root_sisl, 0.0)
        fill!(csite_hydr.sapflow_scpf, 0.0)
        fill!(csite_hydr.rootuptake_sl, 0.0)
        fill!(csite_hydr.rootuptake0_scpf, 0.0)
        fill!(csite_hydr.rootuptake10_scpf, 0.0)
        fill!(csite_hydr.rootuptake50_scpf, 0.0)
        fill!(csite_hydr.rootuptake100_scpf, 0.0)

        sum(@view csite_hydr.l_aroot_layer[1:csite_hydr.nlevrhiz]) == 0.0 && continue

        lat = sites[s].lat
        lon = sites[s].lon
        nlevrhiz = csite_hydr.nlevrhiz

        # Average root water uptake (by rhizosphere shell) across all cohorts.
        dth_layershell_col = zeros(Float64, nlevsoi_hyd_max, nshell)
        csite_hydr.dwat_veg   = 0.0
        csite_hydr.errh2o_hyd = 0.0
        prev_h2oveg  = csite_hydr.h2oveg
        prev_h2osoil = sum(view(csite_hydr.h2osoi_liqvol_shell, 1:nlevrhiz, :) .*
                           view(csite_hydr.v_shell, 1:nlevrhiz, :)) * denh2o * AREA_INV

        fill!(bc_out[s].qflx_ro_sisl, 0.0)

        transp_flux = 0.0
        root_flux   = 0.0

        ifp = 0
        cpatch = sites[s].oldest_patch
        while cpatch !== nothing
            if cpatch.nocomp_pft_label != nocomp_bareground
                ifp += 1

                # Partition the patch transpiration to cohorts by g_sb_laweight.
                gscan_patch = 0.0
                ccohort = cpatch.tallest
                while ccohort !== nothing
                    ccohort_hydr = ccohort.co_hydr
                    ccohort_hydr.psi_ag[1] = psi_from_th(wrf_plant(leaf_p_media, ccohort.pft),
                                                         ccohort_hydr.th_ag[1])
                    gscan_patch += ccohort.g_sb_laweight
                    ccohort = ccohort.shorter
                end

                if bc_in[s].qflx_transp_pa[ifp] > 1.0e-10 && gscan_patch < nearzero
                    fates_endrun("ERROR in plant hydraulics: the HLM predicted a non-zero " *
                                 "total transpiration flux for this patch, yet there is no " *
                                 "leaf-area-weighted conductance. transp=$(bc_in[s].qflx_transp_pa[ifp]) " *
                                 "gscan_patch=$gscan_patch")
                end

                ccohort = cpatch.tallest
                while ccohort !== nothing
                    ccohort_hydr = ccohort.co_hydr
                    ft = ccohort.pft

                    # Relative transpiration of this cohort from the whole patch.
                    if ccohort.g_sb_laweight > nearzero
                        qflx_tran_veg_indiv = bc_in[s].qflx_transp_pa[ifp] * cpatch.total_canopy_area *
                            (ccohort.g_sb_laweight / gscan_patch) / ccohort.n
                    else
                        qflx_tran_veg_indiv = 0.0
                    end

                    ccohort_hydr.qtop = qflx_tran_veg_indiv * dtime

                    transp_flux += (qflx_tran_veg_indiv * dtime) * ccohort.n * AREA_INV

                    # Advance the cohort's flow path (default 1D Taylor solver).
                    if solver == hydr_solver_2DNewton
                        MatSolve2D()
                        sapflow = 0.0; rootuptake = zeros(Float64, nlevrhiz)
                        wb_err_plant = 0.0; dwat_plant = 0.0
                    elseif solver == hydr_solver_2DPicard
                        PicardSolve2D()
                        sapflow = 0.0; rootuptake = zeros(Float64, nlevrhiz)
                        wb_err_plant = 0.0; dwat_plant = 0.0
                    else  # hydr_solver_1DTaylor
                        ordered = copy(ordered_init)
                        OrderLayersForSolve1D(csite_hydr, ccohort, ccohort_hydr, ordered, kbg_layer)
                        sapflow, rootuptake, wb_err_plant, dwat_plant =
                            ImTaylorSolve1D(lat, lon, recruitflag, csite_hydr, ccohort,
                                            ccohort_hydr, dtime, qflx_tran_veg_indiv,
                                            ordered, kbg_layer, dth_layershell_col)
                    end

                    # Remember the cohort error, accumulate site-level diagnostics.
                    ccohort_hydr.errh2o += wb_err_plant
                    csite_hydr.errh2o_hyd += wb_err_plant * ccohort.n * AREA_INV
                    csite_hydr.dwat_veg   += dwat_plant   * ccohort.n * AREA_INV
                    csite_hydr.h2oveg     += dwat_plant   * ccohort.n * AREA_INV

                    sc = ccohort.size_class

                    csite_hydr.sapflow_scpf[sc, ft] += sapflow * ccohort.n / dtime

                    ru = @view rootuptake[1:nlevrhiz]
                    csite_hydr.rootuptake0_scpf[sc, ft]   += SumBetweenDepths(csite_hydr, 0.0,  0.1,   ru) * ccohort.n / dtime
                    csite_hydr.rootuptake10_scpf[sc, ft]  += SumBetweenDepths(csite_hydr, 0.1,  0.5,   ru) * ccohort.n / dtime
                    csite_hydr.rootuptake50_scpf[sc, ft]  += SumBetweenDepths(csite_hydr, 0.5,  1.0,   ru) * ccohort.n / dtime
                    csite_hydr.rootuptake100_scpf[sc, ft] += SumBetweenDepths(csite_hydr, 1.0,  1.0e10, ru) * ccohort.n / dtime

                    # Update water potential / frac total conductivity + btran.
                    UpdatePlantPsiFTCFromTheta!(ccohort, csite_hydr)
                    ccohort_hydr.btran = ftc_from_psi(wkf_plant(stomata_p_media, ft),
                                                      ccohort_hydr.psi_ag[1])

                    ccohort = ccohort.shorter
                end
            end
            cpatch = cpatch.younger
        end

        # ----- Site-level mass-balance / disaggregation -----

        root_flux = -sum(view(dth_layershell_col, 1:nlevrhiz, :) .*
                         view(csite_hydr.v_shell, 1:nlevrhiz, :)) * denh2o * AREA_INV

        fill!(bc_out[s].qflx_soil2root_sisl, 0.0)
        fill!(bc_out[s].qflx_ro_sisl, 0.0)

        # Root density (length) on the soil layer grid, for disaggregation.
        fill!(csite_hydr.rootl_sl, 0.0)
        cpatch = sites[s].oldest_patch
        while cpatch !== nothing
            ccohort = cpatch.tallest
            while ccohort !== nothing
                ccohort_hydr = ccohort.co_hydr
                sum_l_aroot = sum(ccohort_hydr.l_aroot_layer)
                ft = ccohort.pft
                z_fr = MaximumRootingDepth(ccohort.dbh, ft, bc_in[s].zi_sisl[bc_in[s].nlevsoil + 1])
                for j_bc in 1:bc_in[s].nlevsoil
                    rootfr = zeng2001_crootfr(prt_params.fnrt_prof_a[ft], prt_params.fnrt_prof_b[ft],
                                              bc_in[s].zi_sisl[j_bc + 1], z_fr) -
                             zeng2001_crootfr(prt_params.fnrt_prof_a[ft], prt_params.fnrt_prof_b[ft],
                                              bc_in[s].zi_sisl[j_bc + 1] - bc_in[s].dz_sisl[j_bc], z_fr)
                    csite_hydr.rootl_sl[j_bc] += sum_l_aroot * rootfr * ccohort.n *
                                                 prt_params.c2b[ft] * EDPftvarcon_inst[].hydr_srl[ft]
                end
                ccohort = ccohort.shorter
            end
            cpatch = cpatch.younger
        end

        weight_sl = zeros(Float64, nlevsoi_hyd_max)
        for j in 1:nlevrhiz
            if csite_hydr.l_aroot_layer[j] > nearzero
                # Update the rhizosphere shell water content.
                for ish in 1:nshell
                    csite_hydr.h2osoi_liqvol_shell[j, ish] += dth_layershell_col[j, ish]
                end

                # Total root uptake at this rhiz layer [kg/m2/s].
                qflx_soil2root_rhiz =
                    -(sum(view(dth_layershell_col, j, :) .* view(csite_hydr.v_shell, j, :)) *
                      denh2o * AREA_INV / dtime) + csite_hydr.recruit_w_uptake[j]

                # Disaggregate onto the host soil layers.
                j_t = csite_hydr.map_r2s[j, 1]
                j_b = csite_hydr.map_r2s[j, 2]

                sumweight = 0.0
                for j_bc in j_t:j_b
                    if rootflux_disagg == soilk_disagg
                        if qflx_soil2root_rhiz > 0.0
                            eff_por = bc_in[s].eff_porosity_sl[j_bc]
                            h2osoi_liqvol = min(eff_por,
                                bc_in[s].h2o_liq_sisl[j_bc] / (bc_in[s].dz_sisl[j_bc] * denh2o))
                            psi_layer = psi_from_th(csite_hydr.wrf_soil[j], h2osoi_liqvol)
                            ftc_layer = ftc_from_psi(csite_hydr.wkf_soil[j], psi_layer)
                            weight_sl[j_bc] = bc_in[s].hksat_sisl[j_bc] * ftc_layer * csite_hydr.rootl_sl[j_bc]
                        else
                            weight_sl[j_bc] = csite_hydr.rootl_sl[j_bc]
                        end
                    elseif rootflux_disagg == soilz_disagg
                        weight_sl[j_bc] = csite_hydr.rootl_sl[j_bc]
                    else
                        fates_endrun("Unknown rhiz->soil disaggregation method $rootflux_disagg")
                    end
                    sumweight += weight_sl[j_bc]
                end

                # Fall back to root-length weighting if conductance weighting is ~0.
                if sumweight < nearzero
                    sumweight = 0.0
                    for j_bc in j_t:j_b
                        weight_sl[j_bc] = csite_hydr.rootl_sl[j_bc]
                        sumweight += weight_sl[j_bc]
                    end
                end

                for j_bc in j_t:j_b
                    bc_out[s].qflx_soil2root_sisl[j_bc] = qflx_soil2root_rhiz * weight_sl[j_bc] / sumweight
                    csite_hydr.rootuptake_sl[j_bc]      = qflx_soil2root_rhiz * weight_sl[j_bc] / sumweight
                end
            end
        end

        fill!(bc_out[s].qflx_ro_sisl, 0.0)
        site_runoff = sum(bc_out[s].qflx_ro_sisl) * dtime

        delta_plant_storage = csite_hydr.h2oveg - prev_h2oveg
        delta_soil_storage  = sum(view(csite_hydr.h2osoi_liqvol_shell, 1:nlevrhiz, :) .*
                                  view(csite_hydr.v_shell, 1:nlevrhiz, :)) * denh2o * AREA_INV - prev_h2osoil

        # Plant water balance closure (the known solver error is removed).
        if abs(delta_plant_storage - (root_flux + csite_hydr.errh2o_hyd - transp_flux)) > error_thresh
            fates_endrun("Site plant water balance does not close. " *
                         "delta_plant_storage=$delta_plant_storage root_flux=$root_flux " *
                         "transp_flux=$transp_flux errh2o_hyd=$(csite_hydr.errh2o_hyd) " *
                         "dwat_veg=$(csite_hydr.dwat_veg)")
        end

        csite_hydr.h2oveg_hydro_err += csite_hydr.errh2o_hyd

        UpdateH2OVeg!(sites[s], bc_out[s])
    end
    return nothing
end

"""
    hydraulics_drive(nsites, sites, bc_in, bc_out, dtime)

Dispatcher for the FATES plant-hydraulics timestep. For the default
`use_ed_planthydraulics == 1` (BC hydraulics) it first reconciles the
rhizosphere shell water with the host soil column (`FillDrainRhizShells`) then
runs the transpiration/uptake solve (`Hydraulics_BC`). Mirrors the Fortran
`hydraulics_drive`.
"""
function hydraulics_drive(nsites::Integer, sites::AbstractVector,
                          bc_in::AbstractVector, bc_out::AbstractVector, dtime::Real)
    if use_ed_planthydraulics == 1
        FillDrainRhizShells(nsites, sites, bc_in, bc_out)
        Hydraulics_BC(nsites, sites, bc_in, bc_out, dtime)
    elseif use_ed_planthydraulics == 2
        # Hydraulics_CX() — not implemented (CX hydraulics).
    end
    return nothing
end

# ===========================================================================
# Non-default 2D solvers (STUBBED — see module header)
# ===========================================================================

"""
    MatSolve2D(args...)

NON-DEFAULT solver (hydr_solver_2DNewton). Not ported; the default FATES solver is
`hydr_solver_1DTaylor` (`ImTaylorSolve1D`). Calling this raises an error.
"""
MatSolve2D(args...; kwargs...) =
    fates_endrun("MatSolve2D (hydr_solver_2DNewton) is not ported in CLM.jl; " *
                 "use the default hydr_solver_1DTaylor path (ImTaylorSolve1D).")

"""
    PicardSolve2D(args...)

NON-DEFAULT solver (hydr_solver_2DPicard). Not ported; see [`MatSolve2D`](@ref).
"""
PicardSolve2D(args...; kwargs...) =
    fates_endrun("PicardSolve2D (hydr_solver_2DPicard) is not ported in CLM.jl; " *
                 "use the default hydr_solver_1DTaylor path (ImTaylorSolve1D).")
