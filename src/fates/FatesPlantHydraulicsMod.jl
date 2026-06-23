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
        frac = (depth_b - zi[i_rhiz_b]) / dz[i_rhiz_b+1]
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
