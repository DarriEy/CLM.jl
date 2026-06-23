# PRTParamsFATESMod.jl
# Julia port of FATES src/fates/parteh/PRTParamsFATESMod.F90
# (Fortran module name `PRTInitParamsFatesMod`).
#
# This is the PARTEH parameter INITIALIZATION / REGISTRATION layer: it registers
# and receives all the PARTEH / PFT allocation / stoichiometry / turnover
# parameters into the merged `prt_params` (`prt_param_type`, PRTParametersMod.jl)
# via the FATES parameter reader (`fates_parameters_type`,
# FatesParametersInterface.jl), plus the PARTEH consistency checks
# (`PRTCheckParams`) and the derived-parameter setup (`PRTDerivedParams`).
#
# This is the param-wiring counterpart to the already-merged PRTParametersMod
# storage struct. The parameter-file READ stays stubbed (host supplied), exactly
# as in the foundation FatesParametersInterface; here we register names and
# retrieve already-stored data.
#
# Translation notes:
#   * Fortran subroutine + param-name-string names preserved EXACTLY. All
#     `"fates_*"` strings are verbatim.
#   * `fates_r8`  -> Float64.
#   * Fortran `class(fates_parameters_type)` arg -> `fates_parameters_type`.
#   * `ArrayNint` (real -> nearest-int conversion) ported as a helper; the
#     integer-typed `prt_params` fields receive `round`-to-Int.
#   * Fortran `RetrieveParameterAllocate` (allocate-on-receive) -> the
#     `RetrieveParameter{1D,2D}Allocate` helpers, which return freshly allocated
#     arrays; we assign them straight into the `prt_params` fields.
#   * `endrun`/`call endrun` consistency aborts -> Julia `error`. To mirror the
#     Fortran "collect all errors, then abort once" behaviour, PRTCheckParams
#     accumulates messages and throws a single `error` at the end if nerror>0.
#   * `hlm_parteh_mode` is the module Ref{Int} from FatesInterfaceTypesMod.
#
# Deps: FatesConstantsMod (fates_r8, itrue/ifalse, nearzero, years_per_day,
# i*_stress_decid), FatesInterfaceTypesMod (hlm_parteh_mode), PRTParametersMod
# (prt_params), PRTGenericMod (organ/element index consts + StorageNutrientTarget
# + hypothesis consts), FatesParametersInterface (register/retrieve API),
# FatesGlobals (fates_log, fates_endrun), EDPftvarcon (EDPftvarcon_inst),
# FatesAllometryMod (allometry chain), EDTypesMod (init_recruit_trim).

# Lower bounds (Fortran: lower_bound_pft, lower_bound_general; both 1).
const lower_bound_pft     = 1
const lower_bound_general = 1

# =====================================================================================
# ArrayNint — round a real array to nearest integers (Fortran nint).
# =====================================================================================
"""
    ArrayNint(realarr) -> Vector{Int}

Round each element of `realarr` to the nearest integer (mirrors Fortran `nint`).
"""
function ArrayNint(realarr::AbstractVector)
    intarr = Vector{Int}(undef, length(realarr))
    for i in eachindex(realarr)
        intarr[i] = round(Int, realarr[i])
    end
    return intarr
end

# =====================================================================================
# Top-level register / receive orchestrators
# =====================================================================================

"""
    PRTRegisterParams!(fates_params::fates_parameters_type)

Register all PARTEH / PFT / organ / leaf-age parameters with the FATES
parameter reader. Mirrors the Fortran `PRTRegisterParams`.
"""
function PRTRegisterParams!(fates_params::fates_parameters_type)
    PRTRegisterPFT!(fates_params)
    PRTRegisterPFTOrgans!(fates_params)
    PRTRegisterPFTLeafAge!(fates_params)
    Register_PFT_nvariants!(fates_params)
    PRTRegisterOrgan!(fates_params)
    return nothing
end

"""
    PRTReceiveParams!(fates_params::fates_parameters_type)

Receive all PARTEH / PFT / organ / leaf-age parameters from the FATES parameter
reader into the merged `prt_params` struct. Mirrors `PRTReceiveParams`.
"""
function PRTReceiveParams!(fates_params::fates_parameters_type)
    PRTReceivePFT!(fates_params)
    PRTReceivePFTOrgans!(fates_params)
    PRTReceivePFTLeafAge!(fates_params)
    Receive_PFT_nvariants!(fates_params)
    PRTReceiveOrgan!(fates_params)
    return nothing
end

# =====================================================================================
# Organ-dimensioned parameters
# =====================================================================================

"""
    PRTRegisterOrgan!(fates_params)

Register the organ-dimensioned parameter (`fates_alloc_organ_id`).
"""
function PRTRegisterOrgan!(fates_params::fates_parameters_type)
    dim_names = [dimension_name_prt_organs]
    dim_lower_bound = [lower_bound_general]
    RegisterParameter!(fates_params, "fates_alloc_organ_id", dimension_shape_1d,
                       dim_names; lower_bounds=dim_lower_bound)
    return nothing
end

"""
    PRTReceiveOrgan!(fates_params)

Receive the organ-id mapping. Must be called after PRTReceivePFTOrgans (so the
organ dimension is known). Converts the real param to integer organ ids.
"""
function PRTReceiveOrgan!(fates_params::fates_parameters_type)
    tmpreal = RetrieveParameter1DAllocate(fates_params, "fates_alloc_organ_id")
    prt_params.organ_id = ArrayNint(tmpreal)
    return nothing
end

# =====================================================================================
# PFT-dimensioned parameters
# =====================================================================================

# The PFT-indexed parameters that map directly to a real Vector field. Each entry
# is (param-name-string, prt_params field symbol). Names verbatim from Fortran.
const _prt_pft_real_params = [
    ("fates_phen_stem_drop_fraction",          :phen_stem_drop_fraction),
    ("fates_phen_fnrt_drop_fraction",          :phen_fnrt_drop_fraction),
    ("fates_phen_mindaysoff",                  :phen_doff_time),
    ("fates_phen_drought_threshold",           :phen_drought_threshold),
    ("fates_phen_moist_threshold",             :phen_moist_threshold),
    ("fates_leaf_slamax",                      :slamax),
    ("fates_leaf_slatop",                      :slatop),
    ("fates_allom_sai_scaler",                 :allom_sai_scaler),
    ("fates_allom_fnrt_prof_a",                :fnrt_prof_a),
    ("fates_allom_fnrt_prof_b",                :fnrt_prof_b),
    ("fates_allom_fnrt_prof_mode",             :fnrt_prof_mode),
    ("fates_wood_density",                     :wood_density),
    ("fates_recruit_seed_dbh_repro_threshold", :dbh_repro_threshold),
    ("fates_alloc_storage_cushion",            :cushion),
    ("fates_alloc_store_priority_frac",        :leaf_stor_priority),
    ("fates_turnover_senleaf_fdrought",        :senleaf_long_fdrought),
    ("fates_turnover_fnrt",                    :root_long),
    ("fates_leafn_vert_scaler_coeff1",         :leafn_vert_scaler_coeff1),
    ("fates_leafn_vert_scaler_coeff2",         :leafn_vert_scaler_coeff2),
    ("fates_recruit_seed_alloc_mature",        :seed_alloc_mature),
    ("fates_recruit_seed_alloc",               :seed_alloc),
    ("fates_trs_repro_alloc_a",                :repro_alloc_a),
    ("fates_trs_repro_alloc_b",                :repro_alloc_b),
    ("fates_c2b",                              :c2b),
    ("fates_allom_l2fr",                       :allom_l2fr),
    ("fates_cnp_pid_kp",                       :pid_kp),
    ("fates_cnp_pid_ki",                       :pid_ki),
    ("fates_cnp_pid_kd",                       :pid_kd),
    ("fates_cnp_store_ovrflw_frac",            :store_ovrflw_frac),
    ("fates_cnp_nfix1",                        :nfix_mresp_scfrac),
    ("fates_grperc",                           :grperc),
    ("fates_allom_dbh_maxheight",              :allom_dbh_maxheight),
    ("fates_allom_la_per_sa_int",              :allom_la_per_sa_int),
    ("fates_allom_la_per_sa_slp",              :allom_la_per_sa_slp),
    ("fates_allom_agb_frac",                   :allom_agb_frac),
    ("fates_allom_d2h1",                       :allom_d2h1),
    ("fates_allom_d2h2",                       :allom_d2h2),
    ("fates_allom_d2h3",                       :allom_d2h3),
    ("fates_allom_d2bl1",                      :allom_d2bl1),
    ("fates_allom_d2bl2",                      :allom_d2bl2),
    ("fates_allom_d2bl3",                      :allom_d2bl3),
    ("fates_allom_blca_expnt_diff",            :allom_blca_expnt_diff),
    ("fates_allom_d2ca_coefficient_max",       :allom_d2ca_coefficient_max),
    ("fates_allom_d2ca_coefficient_min",       :allom_d2ca_coefficient_min),
    ("fates_allom_agb1",                       :allom_agb1),
    ("fates_allom_agb2",                       :allom_agb2),
    ("fates_allom_agb3",                       :allom_agb3),
    ("fates_allom_agb4",                       :allom_agb4),
    ("fates_allom_h2cd1",                      :allom_h2cd1),
    ("fates_allom_h2cd2",                      :allom_h2cd2),
    ("fates_allom_zroot_max_dbh",              :allom_zroot_max_dbh),
    ("fates_allom_zroot_max_z",                :allom_zroot_max_z),
    ("fates_allom_zroot_min_dbh",              :allom_zroot_min_dbh),
    ("fates_allom_zroot_min_z",                :allom_zroot_min_z),
    ("fates_allom_zroot_k",                    :allom_zroot_k),
    ("fates_turnover_branch",                  :branch_long),
    ("fates_cnp_nitr_store_ratio",             :nitr_store_ratio),
    ("fates_cnp_phos_store_ratio",             :phos_store_ratio),
]

# The PFT-indexed parameters that are stored as INTEGER fields (real on file ->
# nint -> Int). Each entry is (param-name-string, prt_params field symbol).
const _prt_pft_int_params = [
    ("fates_phen_stress_decid", :stress_decid),
    ("fates_phen_season_decid", :season_decid),
    ("fates_phen_evergreen",    :evergreen),
    ("fates_woody",             :woody),
    ("fates_allom_hmode",       :allom_hmode),
    ("fates_allom_lmode",       :allom_lmode),
    ("fates_allom_fmode",       :allom_fmode),
    ("fates_allom_amode",       :allom_amode),
    ("fates_allom_stmode",      :allom_stmode),
    ("fates_allom_cmode",       :allom_cmode),
    ("fates_allom_smode",       :allom_smode),
    ("fates_allom_dmode",       :allom_dmode),
]

"""
    PRTRegisterPFT!(fates_params)

Register all PFT-dimensioned (1-D) PARTEH/allometry parameters. Mirrors the
Fortran `PRTRegisterPFT`.
"""
function PRTRegisterPFT!(fates_params::fates_parameters_type)
    dim_names = [dimension_name_pft]
    dim_lower_bound = [lower_bound_pft]
    reg(nm) = RegisterParameter!(fates_params, nm, dimension_shape_1d, dim_names;
                                 lower_bounds=dim_lower_bound)
    for (nm, _) in _prt_pft_int_params
        reg(nm)
    end
    for (nm, _) in _prt_pft_real_params
        reg(nm)
    end
    return nothing
end

"""
    PRTReceivePFT!(fates_params)

Receive all PFT-dimensioned parameters into `prt_params`. Mirrors
`PRTReceivePFT`.
"""
function PRTReceivePFT!(fates_params::fates_parameters_type)
    for (nm, fld) in _prt_pft_int_params
        tmpreal = RetrieveParameter1DAllocate(fates_params, nm)
        setfield!(prt_params, fld, ArrayNint(tmpreal))
    end
    for (nm, fld) in _prt_pft_real_params
        setfield!(prt_params, fld, RetrieveParameter1DAllocate(fates_params, nm))
    end
    return nothing
end

# =====================================================================================
# PFT x leaf-age parameters
# =====================================================================================

"""
    PRTRegisterPFTLeafAge!(fates_params)

Register the leaf-longevity (PFT x leaf-age) parameters. Mirrors
`PRTRegisterPFTLeafAge`.
"""
function PRTRegisterPFTLeafAge!(fates_params::fates_parameters_type)
    dim_names = [dimension_name_pft, dimension_name_leaf_age]
    dim_lower_bound = [lower_bound_pft, lower_bound_general]
    RegisterParameter!(fates_params, "fates_turnover_leaf_canopy", dimension_shape_2d,
                       dim_names; lower_bounds=dim_lower_bound)
    RegisterParameter!(fates_params, "fates_turnover_leaf_ustory", dimension_shape_2d,
                       dim_names; lower_bounds=dim_lower_bound)
    return nothing
end

"""
    PRTReceivePFTLeafAge!(fates_params)

Receive the leaf-longevity parameters into `prt_params`. Mirrors
`PRTReceivePFTLeafAge`.
"""
function PRTReceivePFTLeafAge!(fates_params::fates_parameters_type)
    prt_params.leaf_long        = RetrieveParameter2DAllocate(fates_params, "fates_turnover_leaf_canopy")
    prt_params.leaf_long_ustory = RetrieveParameter2DAllocate(fates_params, "fates_turnover_leaf_ustory")
    return nothing
end

# =====================================================================================
# PFT x variants (currently unused in FATES; kept as no-op for fidelity)
# =====================================================================================

"""
    Register_PFT_nvariants!(fates_params)

No-op (the Fortran body registers no parameters; the variants dimension is
unused). Kept for one-to-one fidelity with `Register_PFT_nvariants`.
"""
function Register_PFT_nvariants!(fates_params::fates_parameters_type)
    return nothing
end

"""
    Receive_PFT_nvariants!(fates_params)

No-op counterpart of [`Register_PFT_nvariants!`](@ref).
"""
function Receive_PFT_nvariants!(fates_params::fates_parameters_type)
    return nothing
end

# =====================================================================================
# PFT x organ parameters
# =====================================================================================

# (param-name-string, prt_params field symbol). All 2-D (PFT x organ).
const _prt_pft_organ_params = [
    ("fates_stoich_nitr",               :nitr_stoich_p1),
    ("fates_stoich_phos",               :phos_stoich_p1),
    ("fates_alloc_organ_priority",      :alloc_priority),
    ("fates_cnp_turnover_nitr_retrans", :turnover_nitr_retrans),
    ("fates_cnp_turnover_phos_retrans", :turnover_phos_retrans),
]

"""
    PRTRegisterPFTOrgans!(fates_params)

Register the 2-D (PFT x organ) PARTEH stoichiometry / priority / retranslocation
parameters. Mirrors `PRTRegisterPFTOrgans`.
"""
function PRTRegisterPFTOrgans!(fates_params::fates_parameters_type)
    dim_names = [dimension_name_pft, dimension_name_prt_organs]
    dim_lower_bound = [lower_bound_pft, lower_bound_general]
    for (nm, _) in _prt_pft_organ_params
        RegisterParameter!(fates_params, nm, dimension_shape_2d, dim_names;
                           lower_bounds=dim_lower_bound)
    end
    return nothing
end

"""
    PRTReceivePFTOrgans!(fates_params)

Receive the 2-D (PFT x organ) parameters into `prt_params`. Mirrors
`PRTReceivePFTOrgans`.
"""
function PRTReceivePFTOrgans!(fates_params::fates_parameters_type)
    for (nm, fld) in _prt_pft_organ_params
        setfield!(prt_params, fld, RetrieveParameter2DAllocate(fates_params, nm))
    end
    return nothing
end

# =====================================================================================
# Derived parameters
# =====================================================================================

"""
    PRTDerivedParams!()

Set the reverse-lookup map from a global PRTGenericMod organ index to the
parameter-file organ index (`prt_params.organ_param_id`). Entries with no
parameter-file organ stay at -1 (invalid). Mirrors `PRTDerivedParams`.
"""
function PRTDerivedParams!()
    norgans = length(prt_params.organ_id)

    # Allocate the reverse-lookup map and initialise as invalid (-1).
    prt_params.organ_param_id = fill(-1, num_organ_types)

    for i in 1:norgans
        prt_params.organ_param_id[prt_params.organ_id[i]] = i
    end
    return nothing
end

# =====================================================================================
# Consistency checks
# =====================================================================================

"""
    PRTCheckParams(is_master::Bool=true)

Perform logical cross-checks on the user-supplied PARTEH parameters. Mirrors the
Fortran `PRTCheckParams`: all checks are run, errors accumulated, and a single
`error` raised at the end if any check failed (so the caller sees every problem
at once). When `is_master` is false this is a no-op (matching the Fortran early
return).
"""
function PRTCheckParams(is_master::Bool=true)
    is_master || return nothing

    npft    = length(prt_params.evergreen)
    norgans = length(prt_params.organ_id)
    parteh_mode = hlm_parteh_mode[]

    nerror = 0
    msgs = String[]
    note(m) = (nerror += 1; push!(msgs, m))

    # organ ids must be valid global organ indices.
    if any(prt_params.organ_id .< 1) || any(prt_params.organ_id .> num_organ_types)
        note("prt_organ_ids should match the global ids of organ types found in " *
             "PRTGenericMod. organ_ids: $(prt_params.organ_id)")
    end

    # With carbon-only or flexible-CNP allocation, repro & storage are special
    # cases and must NOT appear in the parameter-file organ list.
    if parteh_mode == prt_carbon_allom_hyp || parteh_mode == prt_cnp_flex_allom_hyp
        for io in 1:norgans
            if prt_params.organ_id[io] == repro_organ
                note("with flexible cnp or c-only alloc hypotheses, reproductive " *
                     "tissues are a special case and should not be in the parameter " *
                     "file organ list. organ_id: $(prt_params.organ_id)")
            end
            if prt_params.organ_id[io] == store_organ
                note("with flexible cnp or c-only alloc hypotheses, storage is a " *
                     "special case and should not be in the parameter file organ " *
                     "list. organ_id: $(prt_params.organ_id)")
            end
        end
    end

    # N-fixation respiration surcharge fraction must be in [0,1] (CNP flex).
    if parteh_mode == prt_cnp_flex_allom_hyp
        if any(prt_params.nfix_mresp_scfrac .< 0.0) || any(prt_params.nfix_mresp_scfrac .> 1.0)
            note("The N fixation surcharge nfix_mresp_scfrac (fates_cnp_nfix1) must " *
                 "be between 0-1. values: $(prt_params.nfix_mresp_scfrac)")
        end
    end

    for ipft in 1:npft
        is_evergreen    = prt_params.evergreen[ipft]    == itrue
        is_season_decid = prt_params.season_decid[ipft] == itrue
        is_stress_decid = prt_params.stress_decid[ipft] in (ihard_stress_decid, isemi_stress_decid)
        is_semi_decid   = prt_params.stress_decid[ipft] == isemi_stress_decid

        # Phenology habit must be mutually exclusive.
        if (is_evergreen && is_season_decid) ||
           (is_evergreen && is_stress_decid) ||
           (is_season_decid && is_stress_decid)
            note("PFT #$ipft must have exactly one phenology habit. " *
                 "stress_decid=$(prt_params.stress_decid[ipft]) " *
                 "season_decid=$(prt_params.season_decid[ipft]) " *
                 "evergreen=$(prt_params.evergreen[ipft])")
        end

        # Drought semi-deciduous: moist & dry thresholds must be consistent.
        if is_semi_decid
            if prt_params.phen_drought_threshold[ipft] * prt_params.phen_moist_threshold[ipft] < 0.0
                note("PFT #$ipft: with drought semi-deciduous phenology, the moist " *
                     "threshold must have the same sign as the dry threshold " *
                     "(positive=soil water content, negative=matric potential). " *
                     "drought_threshold=$(prt_params.phen_drought_threshold[ipft]) " *
                     "moist_threshold=$(prt_params.phen_moist_threshold[ipft])")
            elseif prt_params.phen_drought_threshold[ipft] >= prt_params.phen_moist_threshold[ipft]
                note("PFT #$ipft: with drought semi-deciduous phenology, the moist " *
                     "threshold must be greater than (more positive / less negative " *
                     "than, and not equal to) the dry threshold. " *
                     "drought_threshold=$(prt_params.phen_drought_threshold[ipft]) " *
                     "moist_threshold=$(prt_params.phen_moist_threshold[ipft])")
            end
        end

        # Deciduous PFTs: abscission fractions must be bounded.
        if prt_params.evergreen[ipft] == ifalse
            if prt_params.phen_fnrt_drop_fraction[ipft] < 0.0 ||
               prt_params.phen_fnrt_drop_fraction[ipft] > 1.0
                note("PFT #$ipft: fine-root abscission rate must be between 0 and 1 " *
                     "for deciduous PFTs. phen_fnrt_drop_fraction=" *
                     "$(prt_params.phen_fnrt_drop_fraction[ipft])")
            end

            if prt_params.woody[ipft] == itrue &&
               (prt_params.phen_stem_drop_fraction[ipft] < 0.0 ||
                prt_params.phen_stem_drop_fraction[ipft] > nearzero)
                note("PFT #$ipft: non-zero stem-drop fractions are not allowed for " *
                     "woody plants. woody=$(prt_params.woody[ipft]) " *
                     "phen_stem_drop_fraction=$(prt_params.phen_stem_drop_fraction[ipft])")
            elseif prt_params.phen_stem_drop_fraction[ipft] < 0.0 ||
                   prt_params.phen_stem_drop_fraction[ipft] > 1.0
                note("PFT #$ipft: deciduous non-wood plants must keep 0-100% of their " *
                     "stems during the deciduous period. phen_stem_drop_fraction=" *
                     "$(prt_params.phen_stem_drop_fraction[ipft])")
            end
        end

        # Sum of base + mature seed allocation cannot exceed 1.
        if (prt_params.seed_alloc[ipft] + prt_params.seed_alloc_mature[ipft]) > 1.0
            note("PFT #$ipft: the sum of seed allocation from base and mature trees " *
                 "may not exceed 1. seed_alloc=$(prt_params.seed_alloc[ipft]) " *
                 "seed_alloc_mature=$(prt_params.seed_alloc_mature[ipft])")
        end

        # Woody plants must have a non-zero AGB intercept.
        if prt_params.allom_agb1[ipft] <= floatmin(Float64) && prt_params.woody[ipft] == 1
            note("PFT #$ipft: woody plants must have a non-zero intercept in the " *
                 "diameter to AGB allometry. allom_agb1=$(prt_params.allom_agb1[ipft]) " *
                 "woody=$(prt_params.woody[ipft])")
        end

        # allom_d2h2 must be negative when allom_hmode==2 (Weibull / Poorter 2006).
        if prt_params.allom_hmode[ipft] == 2 && prt_params.allom_d2h2[ipft] > 0.0
            note("PFT #$ipft: incorrect height-allometry settings. allom_d2h2 must be " *
                 "negative when allom_hmode==2. allom_hmode=$(prt_params.allom_hmode[ipft]) " *
                 "allom_d2h2=$(prt_params.allom_d2h2[ipft])")
        end

        # Leaf storage priority must be in [0,1].
        if prt_params.leaf_stor_priority[ipft] < 0.0 ||
           prt_params.leaf_stor_priority[ipft] > 1.0
            note("PFT #$ipft: prioritization of carbon allocation to leaf and root " *
                 "turnover replacement must be between 0 and 1. leaf_stor_priority=" *
                 "$(prt_params.leaf_stor_priority[ipft])")
        end

        # CNP-flex hypothesis: nutrient storage ratios and retranslocation checks.
        if parteh_mode == prt_cnp_flex_allom_hyp
            if prt_params.nitr_store_ratio[ipft] < 0.0
                note("PFT #$ipft: with PARTEH CNP hypothesis, nitr_store_ratio must " *
                     "be > 0. nitr_store_ratio=$(prt_params.nitr_store_ratio[ipft])")
            end
            if prt_params.phos_store_ratio[ipft] < 0.0
                note("PFT #$ipft: with PARTEH CNP hypothesis, phos_store_ratio must " *
                     "be > 0. phos_store_ratio=$(prt_params.phos_store_ratio[ipft])")
            end

            for i in 1:norgans
                io = prt_params.organ_id[i]

                if io == sapw_organ
                    if prt_params.turnover_nitr_retrans[ipft, i] > nearzero
                        note("PFT #$ipft: retranslocation of sapwood tissues should " *
                             "be zero. nitrogen retrans=" *
                             "$(prt_params.turnover_nitr_retrans[ipft, i])")
                    end
                    if prt_params.turnover_phos_retrans[ipft, i] > nearzero
                        note("PFT #$ipft: retranslocation of sapwood tissues should " *
                             "be zero. phosphorus retrans=" *
                             "$(prt_params.turnover_phos_retrans[ipft, i])")
                    end
                elseif io == struct_organ
                    if prt_params.turnover_nitr_retrans[ipft, i] > nearzero
                        note("PFT #$ipft: retranslocation of structural tissues " *
                             "should be zero. nitrogen retrans=" *
                             "$(prt_params.turnover_nitr_retrans[ipft, i])")
                    end
                    if prt_params.turnover_phos_retrans[ipft, i] > nearzero
                        note("PFT #$ipft: retranslocation of structural tissues " *
                             "should be zero. phosphorus retrans=" *
                             "$(prt_params.turnover_phos_retrans[ipft, i])")
                    end
                end

                # All other retranslocations must be in [0,1].
                if prt_params.turnover_nitr_retrans[ipft, i] > 1.0 ||
                   prt_params.turnover_phos_retrans[ipft, i] > 1.0 ||
                   prt_params.turnover_nitr_retrans[ipft, i] < 0.0 ||
                   prt_params.turnover_phos_retrans[ipft, i] < 0.0
                    note("PFT #$ipft: retranslocation should range from 0 to 1 " *
                         "(param organ index $i, global index $io). nitr=" *
                         "$(prt_params.turnover_nitr_retrans[ipft, i]) phos=" *
                         "$(prt_params.turnover_phos_retrans[ipft, i])")
                end
            end
        end

        # Growth respiration must be in [0,1].
        if prt_params.grperc[ipft] < 0.0 || prt_params.grperc[ipft] > 1.0
            note("PFT #$ipft: growth respiration must be between 0 and 1: " *
                 "$(prt_params.grperc[ipft])")
        end

        # First nitrogen stoichiometry must be in [0,1) (carbon-only & CNP flex).
        if parteh_mode == prt_carbon_allom_hyp || parteh_mode == prt_cnp_flex_allom_hyp
            if any(@view(prt_params.nitr_stoich_p1[ipft, :]) .< 0.0) ||
               any(@view(prt_params.nitr_stoich_p1[ipft, :]) .>= 1.0)
                note("PFT #$ipft: N per C stoichiometry must be between 0-1. " *
                     "$(prt_params.nitr_stoich_p1[ipft, :])")
            end
        end

        # CNP-flex: per-organ stoichiometry and allocation priorities.
        if parteh_mode == prt_cnp_flex_allom_hyp
            for i in 1:norgans
                if prt_params.nitr_stoich_p1[ipft, i] < 0.0 ||
                   prt_params.phos_stoich_p1[ipft, i] < 0.0 ||
                   prt_params.nitr_stoich_p1[ipft, i] > 1.0 ||
                   prt_params.phos_stoich_p1[ipft, i] > 1.0
                    note("PFT #$ipft: with the CNP flexible-stoichiometry hypothesis, " *
                         "all stoichiometries must be in [0,1] (organ index $i). " *
                         "nitr_stoich=$(prt_params.nitr_stoich_p1[ipft, i]) " *
                         "phos_stoich=$(prt_params.phos_stoich_p1[ipft, i])")
                end
            end

            if any(@view(prt_params.alloc_priority[ipft, :]) .< 0) ||
               any(@view(prt_params.alloc_priority[ipft, :]) .> 6)
                note("PFT #$ipft: allocation priorities should be 0-6 for the CNP " *
                     "flex hypothesis. $(prt_params.alloc_priority[ipft, :])")
            end
        end

        # Leaf turnover time-scales.
        nleafage = size(prt_params.leaf_long, 2)
        for iage in 1:nleafage
            if prt_params.leaf_long[ipft, iage] > nearzero
                if (years_per_day / prt_params.leaf_long[ipft, iage]) > 1.0
                    note("PFT #$ipft, age $iage: leaf turnover time-scale is greater " *
                         "than 1 day! leaf_long=$(prt_params.leaf_long[ipft, iage]) [years]")
                end
                if any(@view(prt_params.leaf_long[ipft, :]) .<= nearzero)
                    note("PFT #$ipft: a leaf_long age-class is zero/invalid yet other " *
                         "age classes are non-zero. leaf_long=$(prt_params.leaf_long[ipft, :])")
                end
            else
                if prt_params.evergreen[ipft] == itrue
                    note("PFT #$ipft, age $iage: zero leaf turnover specified yet this " *
                         "is an evergreen PFT. leaf_long=$(prt_params.leaf_long[ipft, iage])")
                end
            end
        end

        # Senescing-leaf-pool turnover.
        if prt_params.leaf_long[ipft, nleafage] > nearzero
            if (years_per_day /
                (prt_params.leaf_long[ipft, nleafage] *
                 prt_params.senleaf_long_fdrought[ipft])) > 1.0
                note("PFT #$ipft: drought-senescent turnover time-scale is greater " *
                     "than 1 day! leaf_long*senleaf_long_fdrought=" *
                     "$(prt_params.leaf_long[ipft, nleafage] * prt_params.senleaf_long_fdrought[ipft]) [years]")
            end
        end

        if prt_params.senleaf_long_fdrought[ipft] < nearzero ||
           prt_params.senleaf_long_fdrought[ipft] > 1.0
            note("PFT #$ipft: senleaf_long_fdrought must be > 0 and <= 1 " *
                 "(set to 1 for no accelerated senescence). " *
                 "senleaf_long_fdrought=$(prt_params.senleaf_long_fdrought[ipft])")
        end

        # Root turnover time-scale.
        if prt_params.root_long[ipft] > nearzero
            if (years_per_day / prt_params.root_long[ipft]) > 1.0
                note("PFT #$ipft: root turnover time-scale is greater than 1 day! " *
                     "root_long=$(prt_params.root_long[ipft]) [years]")
            end
        else
            if prt_params.evergreen[ipft] == itrue
                note("PFT #$ipft: zero root turnover specified yet this is an " *
                     "evergreen PFT. root_long=$(prt_params.root_long[ipft])")
            end
        end

        # Branch turnover time-scale.
        if prt_params.branch_long[ipft] > nearzero
            if (years_per_day / prt_params.branch_long[ipft]) > 1.0
                note("PFT #$ipft: branch turnover time-scale is greater than 1 day! " *
                     "branch_long=$(prt_params.branch_long[ipft]) [years]")
            end
        end
    end

    if nerror > 0
        fates_endrun("One or more PARTEH parameter errors found. Aborting.\n" *
                     join(msgs, "\n"))
    end
    return nothing
end

# =====================================================================================
# NewRecruitTotalStoichiometry
# =====================================================================================

"""
    NewRecruitTotalStoichiometry(ft, l2fr, element_id) -> recruit_stoich

Total N:C or P:C ratio for a newly recruited plant of PFT `ft`. Identifies the
recruit dbh (from `hgt_min`), uses allometry to compute the starting target
carbon in each organ, then applies the stoichiometry parameters to get the total
nutrient. Returns nutrient-to-carbon ratio. Mirrors the Fortran
`NewRecruitTotalStoichiometry`.
"""
function NewRecruitTotalStoichiometry(ft::Integer, l2fr::Real, element_id::Integer)
    not_damaged = 1  # also in DamageMainMod; here for dependency purposes.

    edpft = EDPftvarcon_inst[]

    # For all tissues, assume tissues of new recruits will be fully flushed.
    dbh, _ = h2d_allom(edpft.hgt_min[ft], ft)
    c_leaf, _        = bleaf(dbh, ft, not_damaged, init_recruit_trim, 1.0)
    c_fnrt, _        = bfineroot(dbh, ft, init_recruit_trim, l2fr, 1.0)
    a_sapw, c_sapw, _ = bsap_allom(dbh, ft, not_damaged, init_recruit_trim, 1.0)
    c_agw, dbagwdd   = bagw_allom(dbh, ft, not_damaged, 1.0)
    c_bgw, _         = bbgw_allom(dbh, ft, 1.0)
    c_struct, _      = bdead_allom(c_agw, c_bgw, c_sapw, ft)
    c_store, _       = bstore_allom(dbh, ft, not_damaged, init_recruit_trim)

    # Total carbon in a newly recruited plant.
    c_total = c_leaf + c_fnrt + c_sapw + c_struct + c_store

    n_leaf   = prt_params.organ_param_id[leaf_organ]
    n_fnrt   = prt_params.organ_param_id[fnrt_organ]
    n_sapw   = prt_params.organ_param_id[sapw_organ]
    n_struct = prt_params.organ_param_id[struct_organ]

    if element_id == nitrogen_element
        stoich = prt_params.nitr_stoich_p1
    elseif element_id == phosphorus_element
        stoich = prt_params.phos_stoich_p1
    else
        fates_endrun("NewRecruitTotalStoichiometry: unsupported element_id $element_id")
    end

    nutr_total =
        c_struct * stoich[ft, n_struct] +
        c_leaf   * stoich[ft, n_leaf] +
        c_fnrt   * stoich[ft, n_fnrt] +
        c_sapw   * stoich[ft, n_sapw] +
        StorageNutrientTarget(ft, element_id,
            c_leaf   * stoich[ft, n_leaf],
            c_fnrt   * stoich[ft, n_fnrt],
            c_sapw   * stoich[ft, n_sapw],
            c_struct * stoich[ft, n_struct])

    return nutr_total / c_total
end
