# ==========================================================================
# Ported from: src/soilbiogeochem/SoilBiogeochemDecompCascadeMIMICSMod.F90
# Sets the coefficients used in the decomposition cascade submodel.
# This uses the MIMICS parameters.
#
# Public functions:
#   decomp_mimics_read_params!       — Read MIMICS decomposition parameters
#   init_decompcascade_mimics!       — Initialize rate constants and pathways
#   decomp_rates_mimics!             — Calculate decomposition rates
# ==========================================================================

# ---------------------------------------------------------------------------
# DecompMIMICSParams — MIMICS decomposition parameters
# Ported from params_type in SoilBiogeochemDecompCascadeMIMICSMod.F90
# ---------------------------------------------------------------------------

"""
    DecompMIMICSParams

MIMICS decomposition cascade parameters. Holds microbial growth efficiencies,
Vmax/Km regression parameters, turnover parameters, and initial C stocks.

Ported from `params_type` in `SoilBiogeochemDecompCascadeMIMICSMod.F90`.
"""
Base.@kwdef mutable struct DecompMIMICSParams
    mimics_nue_into_mic           ::Float64          = 0.0
    mimics_desorpQ10              ::Float64          = 0.0
    mimics_densdep                ::Float64          = 1.0
    mimics_tau_mod_factor         ::Float64          = 0.0
    mimics_tau_mod_min            ::Float64          = 0.0
    mimics_tau_mod_max            ::Float64          = 0.0
    mimics_ko_r                   ::Float64          = 0.0
    mimics_ko_k                   ::Float64          = 0.0
    mimics_cn_r                   ::Float64          = 0.0
    mimics_cn_k                   ::Float64          = 0.0
    mimics_cn_mod_num             ::Float64          = 0.0
    mimics_t_soi_ref              ::Float64          = 25.0
    mimics_initial_Cstocks_depth  ::Float64          = 0.3
    mimics_initial_Cstocks        ::Vector{Float64}  = Float64[]
    mimics_mge                    ::Vector{Float64}  = Float64[]
    mimics_vmod                   ::Vector{Float64}  = Float64[]
    mimics_vint                   ::Vector{Float64}  = Float64[]
    mimics_vslope                 ::Vector{Float64}  = Float64[]
    mimics_kmod                   ::Vector{Float64}  = Float64[]
    mimics_kint                   ::Vector{Float64}  = Float64[]
    mimics_kslope                 ::Vector{Float64}  = Float64[]
    mimics_fmet                   ::Vector{Float64}  = Float64[]
    mimics_p_scalar               ::Vector{Float64}  = Float64[]
    mimics_fphys_r                ::Vector{Float64}  = Float64[]
    mimics_fphys_k                ::Vector{Float64}  = Float64[]
    mimics_fchem_r                ::Vector{Float64}  = Float64[]
    mimics_fchem_k                ::Vector{Float64}  = Float64[]
    mimics_desorp                 ::Vector{Float64}  = Float64[]
    mimics_tau_r                  ::Vector{Float64}  = Float64[]
    mimics_tau_k                  ::Vector{Float64}  = Float64[]
end

# ---------------------------------------------------------------------------
# DecompMIMICSState — module-level state for MIMICS decomposition
# Ported from module-level private variables in
# SoilBiogeochemDecompCascadeMIMICSMod.F90
# ---------------------------------------------------------------------------

"""
    DecompMIMICSState

Module-level persistent state for the MIMICS decomposition cascade.
Pool indices, transition indices, respiration fractions, Vmax/Km regression
coefficients, and spatially-varying arrays are set once during
`init_decompcascade_mimics!` and used by `decomp_rates_mimics!`.

Ported from module-level private variables in
`SoilBiogeochemDecompCascadeMIMICSMod.F90`.
"""
Base.@kwdef mutable struct DecompMIMICSState
    # Spatially-varying arrays (col × nlevdecomp)
    desorp    ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    fphys_m1  ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    fphys_m2  ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    p_scalar  ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)

    # Pool indices
    i_phys_som ::Int = 0
    i_chem_som ::Int = 0
    i_avl_som  ::Int = 0
    i_str_lit  ::Int = 0
    i_met_lit  ::Int = 0
    i_cop_mic  ::Int = 0
    i_oli_mic  ::Int = 0

    # Transition indices
    i_l1m1 ::Int = 0
    i_l1m2 ::Int = 0
    i_l2m1 ::Int = 0
    i_l2m2 ::Int = 0
    i_s1m1 ::Int = 0
    i_s1m2 ::Int = 0
    i_m1s1 ::Int = 0
    i_m1s2 ::Int = 0
    i_m2s1 ::Int = 0
    i_m2s2 ::Int = 0
    i_s2s1 ::Int = 0
    i_s3s1 ::Int = 0
    i_m1s3 ::Int = 0
    i_m2s3 ::Int = 0

    # Respiration fractions
    rf_l1m1 ::Float64 = 0.0
    rf_l1m2 ::Float64 = 0.0
    rf_l2m1 ::Float64 = 0.0
    rf_l2m2 ::Float64 = 0.0
    rf_s1m1 ::Float64 = 0.0
    rf_s1m2 ::Float64 = 0.0

    # Vmax regression parameters
    vint_l1_m1   ::Float64 = 0.0
    vint_l2_m1   ::Float64 = 0.0
    vint_s1_m1   ::Float64 = 0.0
    vint_l1_m2   ::Float64 = 0.0
    vint_l2_m2   ::Float64 = 0.0
    vint_s1_m2   ::Float64 = 0.0
    kint_l1_m1   ::Float64 = 0.0
    kint_l2_m1   ::Float64 = 0.0
    kint_s1_m1   ::Float64 = 0.0
    kint_l1_m2   ::Float64 = 0.0
    kint_l2_m2   ::Float64 = 0.0
    kint_s1_m2   ::Float64 = 0.0
    vmod_l1_m1   ::Float64 = 0.0
    vmod_l2_m1   ::Float64 = 0.0
    vmod_s1_m1   ::Float64 = 0.0
    vmod_l1_m2   ::Float64 = 0.0
    vmod_l2_m2   ::Float64 = 0.0
    vmod_s1_m2   ::Float64 = 0.0
    kmod_l1_m1   ::Float64 = 0.0
    kmod_l2_m1   ::Float64 = 0.0
    kmod_s1_m1   ::Float64 = 0.0
    kmod_l1_m2   ::Float64 = 0.0
    kmod_l2_m2   ::Float64 = 0.0
    kmod_s1_m2   ::Float64 = 0.0
    vslope_l1_m1 ::Float64 = 0.0
    vslope_l2_m1 ::Float64 = 0.0
    vslope_s1_m1 ::Float64 = 0.0
    vslope_l1_m2 ::Float64 = 0.0
    vslope_l2_m2 ::Float64 = 0.0
    vslope_s1_m2 ::Float64 = 0.0
    kslope_l1_m1 ::Float64 = 0.0
    kslope_l2_m1 ::Float64 = 0.0
    kslope_s1_m1 ::Float64 = 0.0
    kslope_l1_m2 ::Float64 = 0.0
    kslope_l2_m2 ::Float64 = 0.0
    kslope_s1_m2 ::Float64 = 0.0
end

# ---------------------------------------------------------------------------
# decomp_mimics_read_params! — Read MIMICS decomposition parameters
# Ported from readParams in SoilBiogeochemDecompCascadeMIMICSMod.F90
# ---------------------------------------------------------------------------

"""
    decomp_mimics_read_params!(params; kwargs...)

Populate `DecompMIMICSParams` from keyword arguments (replaces NetCDF file
reading). Each keyword maps to a parameter variable name from the Fortran
`readParams` subroutine.

Ported from `readParams` in `SoilBiogeochemDecompCascadeMIMICSMod.F90`.
"""
function decomp_mimics_read_params!(params::DecompMIMICSParams;
                                     mimics_initial_Cstocks_depth::Real,
                                     mimics_initial_Cstocks::Vector{<:Real},
                                     mimics_mge::Vector{<:Real},
                                     mimics_vmod::Vector{<:Real},
                                     mimics_vslope::Vector{<:Real},
                                     mimics_vint::Vector{<:Real},
                                     mimics_kmod::Vector{<:Real},
                                     mimics_kslope::Vector{<:Real},
                                     mimics_kint::Vector{<:Real},
                                     mimics_p_scalar::Vector{<:Real},
                                     mimics_desorp::Vector{<:Real},
                                     mimics_fphys_r::Vector{<:Real},
                                     mimics_fphys_k::Vector{<:Real},
                                     mimics_fmet::Vector{<:Real},
                                     mimics_fchem_r::Vector{<:Real},
                                     mimics_fchem_k::Vector{<:Real},
                                     mimics_tau_r::Vector{<:Real},
                                     mimics_tau_k::Vector{<:Real},
                                     mimics_nue_into_mic::Real,
                                     mimics_tau_mod_factor::Real,
                                     mimics_tau_mod_min::Real,
                                     mimics_tau_mod_max::Real,
                                     mimics_ko_r::Real,
                                     mimics_ko_k::Real,
                                     mimics_densdep::Real,
                                     mimics_desorpQ10::Real,
                                     mimics_t_soi_ref::Real,
                                     mimics_cn_mod_num::Real,
                                     mimics_cn_r::Real,
                                     mimics_cn_k::Real)
    params.mimics_initial_Cstocks_depth = mimics_initial_Cstocks_depth
    params.mimics_initial_Cstocks       = copy(mimics_initial_Cstocks)
    params.mimics_mge                   = copy(mimics_mge)
    params.mimics_vmod                  = copy(mimics_vmod)
    params.mimics_vslope                = copy(mimics_vslope)
    params.mimics_vint                  = copy(mimics_vint)
    params.mimics_kmod                  = copy(mimics_kmod)
    params.mimics_kslope                = copy(mimics_kslope)
    params.mimics_kint                  = copy(mimics_kint)
    params.mimics_p_scalar              = copy(mimics_p_scalar)
    params.mimics_desorp                = copy(mimics_desorp)
    params.mimics_fphys_r               = copy(mimics_fphys_r)
    params.mimics_fphys_k               = copy(mimics_fphys_k)
    params.mimics_fmet                  = copy(mimics_fmet)
    params.mimics_fchem_r               = copy(mimics_fchem_r)
    params.mimics_fchem_k               = copy(mimics_fchem_k)
    params.mimics_tau_r                 = copy(mimics_tau_r)
    params.mimics_tau_k                 = copy(mimics_tau_k)
    params.mimics_nue_into_mic          = mimics_nue_into_mic
    params.mimics_tau_mod_factor        = mimics_tau_mod_factor
    params.mimics_tau_mod_min           = mimics_tau_mod_min
    params.mimics_tau_mod_max           = mimics_tau_mod_max
    params.mimics_ko_r                  = mimics_ko_r
    params.mimics_ko_k                  = mimics_ko_k
    params.mimics_densdep               = mimics_densdep
    params.mimics_desorpQ10             = mimics_desorpQ10
    params.mimics_t_soi_ref             = mimics_t_soi_ref
    params.mimics_cn_mod_num            = mimics_cn_mod_num
    params.mimics_cn_r                  = mimics_cn_r
    params.mimics_cn_k                  = mimics_cn_k
    return nothing
end

# ---------------------------------------------------------------------------
# init_decompcascade_mimics! — Initialize cascade pathways and coefficients
# Ported from init_decompcascade_mimics in
# SoilBiogeochemDecompCascadeMIMICSMod.F90
# ---------------------------------------------------------------------------

"""
    init_decompcascade_mimics!(mimics_state, cascade_con, params, cn_params;
                                cellclay, bounds, nlevdecomp,
                                ndecomp_pools_max,
                                ndecomp_cascade_transitions_max,
                                spinup_state, use_fates)

Initialize rate constants and decomposition pathways following the
decomposition cascade of the MIMICS model.

Ported from `init_decompcascade_mimics` in
`SoilBiogeochemDecompCascadeMIMICSMod.F90`.

# Arguments
- `mimics_state::DecompMIMICSState`: module-level state to populate
- `cascade_con::DecompCascadeConData`: cascade configuration to populate
- `params::DecompMIMICSParams`: MIMICS parameters
- `cn_params::CNSharedParamsData`: shared CN parameters

# Keyword Arguments
- `cellclay::Matrix{<:Real}`: column clay fraction (%), (col × nlevdecomp)
- `bounds::UnitRange{Int}`: column bounds
- `nlevdecomp::Int`: number of decomposition levels
- `ndecomp_pools_max::Int`: maximum number of decomposition pools
- `ndecomp_cascade_transitions_max::Int`: max number of cascade transitions
- `spinup_state::Int`: spinup state flag (0=normal)
- `use_fates::Bool`: whether FATES is enabled
"""
function init_decompcascade_mimics!(mimics_state::DecompMIMICSState,
                                     cascade_con::DecompCascadeConData,
                                     params::DecompMIMICSParams,
                                     cn_params::CNSharedParamsData;
                                     cellclay::Matrix{<:Real},
                                     bounds::UnitRange{Int},
                                     nlevdecomp::Int,
                                     ndecomp_pools_max::Int,
                                     ndecomp_cascade_transitions_max::Int=15,
                                     spinup_state::Int=0,
                                     use_fates::Bool=false)

    nc = length(bounds)

    # Allocate cascade_con arrays
    FT = typeof(cascade_con.initial_stock_soildepth)
    cascade_con.cascade_donor_pool             = zeros(Int, ndecomp_cascade_transitions_max)
    cascade_con.cascade_receiver_pool          = zeros(Int, ndecomp_cascade_transitions_max)
    cascade_con.floating_cn_ratio_decomp_pools = falses(ndecomp_pools_max)
    cascade_con.is_litter                      = falses(ndecomp_pools_max)
    cascade_con.is_soil                        = falses(ndecomp_pools_max)
    cascade_con.is_cwd                         = falses(ndecomp_pools_max)
    cascade_con.initial_cn_ratio               = zeros(FT, ndecomp_pools_max)
    cascade_con.initial_stock                  = zeros(FT, ndecomp_pools_max)
    cascade_con.is_metabolic                   = falses(ndecomp_pools_max)
    cascade_con.is_cellulose                   = falses(ndecomp_pools_max)
    cascade_con.is_lignin                      = falses(ndecomp_pools_max)
    cascade_con.spinup_factor                  = ones(FT, ndecomp_pools_max)
    cascade_con.decomp_pool_name_restart       = fill("", ndecomp_pools_max)
    cascade_con.decomp_pool_name_history       = fill("", ndecomp_pools_max)
    cascade_con.decomp_pool_name_long          = fill("", ndecomp_pools_max)
    cascade_con.decomp_pool_name_short         = fill("", ndecomp_pools_max)
    cascade_con.cascade_step_name              = fill("", ndecomp_cascade_transitions_max)

    # Allocate spatially-varying arrays
    mimics_state.desorp   = zeros(nc, nlevdecomp)
    mimics_state.fphys_m1 = zeros(nc, nlevdecomp)
    mimics_state.fphys_m2 = zeros(nc, nlevdecomp)
    mimics_state.p_scalar = zeros(nc, nlevdecomp)

    # --- Time-constant coefficients ---
    mimics_nue_into_mic = params.mimics_nue_into_mic
    mimics_p_scalar_p1  = params.mimics_p_scalar[1]
    mimics_p_scalar_p2  = params.mimics_p_scalar[2]
    mimics_fphys_r_p1   = params.mimics_fphys_r[1]
    mimics_fphys_r_p2   = params.mimics_fphys_r[2]
    mimics_fphys_k_p1   = params.mimics_fphys_k[1]
    mimics_fphys_k_p2   = params.mimics_fphys_k[2]
    mimics_desorp_p1    = params.mimics_desorp[1]
    mimics_desorp_p2    = params.mimics_desorp[2]

    # Set respiration fractions for fluxes between compartments
    mimics_state.rf_l1m1 = 1.0 - params.mimics_mge[1]
    mimics_state.rf_l2m1 = 1.0 - params.mimics_mge[2]
    mimics_state.rf_s1m1 = 1.0 - params.mimics_mge[3]
    mimics_state.rf_l1m2 = 1.0 - params.mimics_mge[4]
    mimics_state.rf_l2m2 = 1.0 - params.mimics_mge[5]
    mimics_state.rf_s1m2 = 1.0 - params.mimics_mge[6]

    # vmod = "old" vmod * av  AND  kmod = ak / "old" kmod
    mimics_state.vmod_l1_m1 = params.mimics_vmod[1]
    mimics_state.vmod_l2_m1 = params.mimics_vmod[2]
    mimics_state.vmod_s1_m1 = params.mimics_vmod[3]
    mimics_state.vmod_l1_m2 = params.mimics_vmod[4]
    mimics_state.vmod_l2_m2 = params.mimics_vmod[5]
    mimics_state.vmod_s1_m2 = params.mimics_vmod[6]
    mimics_state.kmod_l1_m1 = params.mimics_kmod[1]
    mimics_state.kmod_l2_m1 = params.mimics_kmod[2]
    mimics_state.kmod_s1_m1 = params.mimics_kmod[3]
    mimics_state.kmod_l1_m2 = params.mimics_kmod[4]
    mimics_state.kmod_l2_m2 = params.mimics_kmod[5]
    mimics_state.kmod_s1_m2 = params.mimics_kmod[6]
    mimics_state.vslope_l1_m1 = params.mimics_vslope[1]
    mimics_state.vslope_l2_m1 = params.mimics_vslope[2]
    mimics_state.vslope_s1_m1 = params.mimics_vslope[3]
    mimics_state.vslope_l1_m2 = params.mimics_vslope[4]
    mimics_state.vslope_l2_m2 = params.mimics_vslope[5]
    mimics_state.vslope_s1_m2 = params.mimics_vslope[6]
    mimics_state.kslope_l1_m1 = params.mimics_kslope[1]
    mimics_state.kslope_l2_m1 = params.mimics_kslope[2]
    mimics_state.kslope_s1_m1 = params.mimics_kslope[3]
    mimics_state.kslope_l1_m2 = params.mimics_kslope[4]
    mimics_state.kslope_l2_m2 = params.mimics_kslope[5]
    mimics_state.kslope_s1_m2 = params.mimics_kslope[6]
    mimics_state.vint_l1_m1 = params.mimics_vint[1]
    mimics_state.vint_l2_m1 = params.mimics_vint[2]
    mimics_state.vint_s1_m1 = params.mimics_vint[3]
    mimics_state.vint_l1_m2 = params.mimics_vint[4]
    mimics_state.vint_l2_m2 = params.mimics_vint[5]
    mimics_state.vint_s1_m2 = params.mimics_vint[6]
    mimics_state.kint_l1_m1 = params.mimics_kint[1]
    mimics_state.kint_l2_m1 = params.mimics_kint[2]
    mimics_state.kint_s1_m1 = params.mimics_kint[3]
    mimics_state.kint_l1_m2 = params.mimics_kint[4]
    mimics_state.kint_l2_m2 = params.mimics_kint[5]
    mimics_state.kint_s1_m2 = params.mimics_kint[6]

    # Compute soil texture-dependent parameters
    for c in 1:nc
        for j in 1:nlevdecomp
            clay_frac = PCT_TO_FRAC * smooth_min(100.0, cellclay[c, j])
            mimics_state.desorp[c, j] = mimics_desorp_p1 * exp(mimics_desorp_p2 * clay_frac)
            mimics_state.fphys_m1[c, j] = smooth_min(1.0, mimics_fphys_r_p1 *
                                               exp(mimics_fphys_r_p2 * clay_frac))
            mimics_state.fphys_m2[c, j] = smooth_min(1.0, mimics_fphys_k_p1 *
                                               exp(mimics_fphys_k_p2 * clay_frac))
            mimics_state.p_scalar[c, j] = 1.0 / (mimics_p_scalar_p1 *
                                                   exp(mimics_p_scalar_p2 * sqrt(clay_frac)))
        end
    end
    cascade_con.initial_stock_soildepth = params.mimics_initial_Cstocks_depth

    # --- List of pools and their attributes ---

    # i_met_lit (metabolic litter)
    i_met_lit = 1
    mimics_state.i_met_lit = i_met_lit
    cascade_con.floating_cn_ratio_decomp_pools[i_met_lit] = true
    cascade_con.decomp_pool_name_restart[i_met_lit]       = "litr1"
    cascade_con.decomp_pool_name_history[i_met_lit]       = "LIT_MET"
    cascade_con.decomp_pool_name_long[i_met_lit]          = "metabolic litter"
    cascade_con.decomp_pool_name_short[i_met_lit]         = "L1"
    cascade_con.is_litter[i_met_lit]    = true
    cascade_con.is_soil[i_met_lit]      = false
    cascade_con.is_cwd[i_met_lit]       = false
    cascade_con.initial_cn_ratio[i_met_lit] = 10.0
    cascade_con.initial_stock[i_met_lit]    = params.mimics_initial_Cstocks[i_met_lit]
    cascade_con.is_metabolic[i_met_lit] = true
    cascade_con.is_cellulose[i_met_lit] = false
    cascade_con.is_lignin[i_met_lit]    = false

    # i_str_lit (structural litter)
    i_str_lit = i_met_lit + 1
    mimics_state.i_str_lit = i_str_lit
    cascade_con.floating_cn_ratio_decomp_pools[i_str_lit] = true
    cascade_con.decomp_pool_name_restart[i_str_lit]       = "litr2"
    cascade_con.decomp_pool_name_history[i_str_lit]       = "LIT_STR"
    cascade_con.decomp_pool_name_long[i_str_lit]          = "structural litter"
    cascade_con.decomp_pool_name_short[i_str_lit]         = "L2"
    cascade_con.is_litter[i_str_lit]    = true
    cascade_con.is_soil[i_str_lit]      = false
    cascade_con.is_cwd[i_str_lit]       = false
    cascade_con.initial_cn_ratio[i_str_lit] = 10.0
    cascade_con.initial_stock[i_str_lit]    = params.mimics_initial_Cstocks[i_str_lit]
    cascade_con.is_metabolic[i_str_lit] = false
    cascade_con.is_cellulose[i_str_lit] = true
    cascade_con.is_lignin[i_str_lit]    = true

    # i_avl_som (available SOM)
    i_avl_som = i_str_lit + 1
    mimics_state.i_avl_som = i_avl_som
    cascade_con.floating_cn_ratio_decomp_pools[i_avl_som] = true
    cascade_con.decomp_pool_name_restart[i_avl_som]       = "soil1"
    cascade_con.decomp_pool_name_history[i_avl_som]       = "SOM_AVL"
    cascade_con.decomp_pool_name_long[i_avl_som]          = "available soil organic matter"
    cascade_con.decomp_pool_name_short[i_avl_som]         = "S1"
    cascade_con.is_litter[i_avl_som]    = false
    cascade_con.is_soil[i_avl_som]      = true
    cascade_con.is_cwd[i_avl_som]       = false
    cascade_con.initial_cn_ratio[i_avl_som] = 10.0
    cascade_con.initial_stock[i_avl_som]    = params.mimics_initial_Cstocks[i_avl_som]
    cascade_con.is_metabolic[i_avl_som] = false
    cascade_con.is_cellulose[i_avl_som] = false
    cascade_con.is_lignin[i_avl_som]    = false

    # i_chem_som (chemically protected SOM)
    i_chem_som = i_avl_som + 1
    mimics_state.i_chem_som = i_chem_som
    cascade_con.floating_cn_ratio_decomp_pools[i_chem_som] = true
    cascade_con.decomp_pool_name_restart[i_chem_som]       = "soil2"
    cascade_con.decomp_pool_name_history[i_chem_som]       = "SOM_CHEM"
    cascade_con.decomp_pool_name_long[i_chem_som]          = "chemically protected soil organic matter"
    cascade_con.decomp_pool_name_short[i_chem_som]         = "S2"
    cascade_con.is_litter[i_chem_som]    = false
    cascade_con.is_soil[i_chem_som]      = true
    cascade_con.is_cwd[i_chem_som]       = false
    cascade_con.initial_cn_ratio[i_chem_som] = 10.0
    cascade_con.initial_stock[i_chem_som]    = params.mimics_initial_Cstocks[i_chem_som]
    cascade_con.is_metabolic[i_chem_som] = false
    cascade_con.is_cellulose[i_chem_som] = false
    cascade_con.is_lignin[i_chem_som]    = false

    # i_phys_som (physically protected SOM)
    i_phys_som = i_chem_som + 1
    mimics_state.i_phys_som = i_phys_som
    cascade_con.floating_cn_ratio_decomp_pools[i_phys_som] = true
    cascade_con.decomp_pool_name_restart[i_phys_som]       = "soil3"
    cascade_con.decomp_pool_name_history[i_phys_som]       = "SOM_PHYS"
    cascade_con.decomp_pool_name_long[i_phys_som]          = "physically protected soil organic matter"
    cascade_con.decomp_pool_name_short[i_phys_som]         = "S3"
    cascade_con.is_litter[i_phys_som]    = false
    cascade_con.is_soil[i_phys_som]      = true
    cascade_con.is_cwd[i_phys_som]       = false
    cascade_con.initial_cn_ratio[i_phys_som] = 10.0
    cascade_con.initial_stock[i_phys_som]    = params.mimics_initial_Cstocks[i_phys_som]
    cascade_con.is_metabolic[i_phys_som] = false
    cascade_con.is_cellulose[i_phys_som] = false
    cascade_con.is_lignin[i_phys_som]    = false

    # i_cop_mic (copiotrophic microbes, M1)
    i_cop_mic = i_phys_som + 1
    mimics_state.i_cop_mic = i_cop_mic
    cascade_con.floating_cn_ratio_decomp_pools[i_cop_mic] = true
    cascade_con.decomp_pool_name_restart[i_cop_mic]       = "micr1"
    cascade_con.decomp_pool_name_history[i_cop_mic]       = "MIC_COP"
    cascade_con.decomp_pool_name_long[i_cop_mic]          = "copiotrophic microbes"
    cascade_con.decomp_pool_name_short[i_cop_mic]         = "M1"
    cascade_con.is_litter[i_cop_mic]    = false
    cascade_con.is_soil[i_cop_mic]      = false
    cascade_con.is_cwd[i_cop_mic]       = false
    cascade_con.initial_cn_ratio[i_cop_mic] = 10.0
    cascade_con.initial_stock[i_cop_mic]    = params.mimics_initial_Cstocks[i_cop_mic]
    cascade_con.is_metabolic[i_cop_mic] = false
    cascade_con.is_cellulose[i_cop_mic] = false
    cascade_con.is_lignin[i_cop_mic]    = false

    # i_oli_mic (oligotrophic microbes, M2)
    i_oli_mic = i_cop_mic + 1
    mimics_state.i_oli_mic = i_oli_mic
    cascade_con.floating_cn_ratio_decomp_pools[i_oli_mic] = true
    cascade_con.decomp_pool_name_restart[i_oli_mic]       = "micr2"
    cascade_con.decomp_pool_name_history[i_oli_mic]       = "MIC_OLI"
    cascade_con.decomp_pool_name_long[i_oli_mic]          = "oligotrophic microbes"
    cascade_con.decomp_pool_name_short[i_oli_mic]         = "M2"
    cascade_con.is_litter[i_oli_mic]    = false
    cascade_con.is_soil[i_oli_mic]      = false
    cascade_con.is_cwd[i_oli_mic]       = false
    cascade_con.initial_cn_ratio[i_oli_mic] = 10.0
    cascade_con.initial_stock[i_oli_mic]    = params.mimics_initial_Cstocks[i_oli_mic]
    cascade_con.is_metabolic[i_oli_mic] = false
    cascade_con.is_cellulose[i_oli_mic] = false
    cascade_con.is_lignin[i_oli_mic]    = false

    # CWD (only if FATES not enabled)
    i_cwd_local = 0
    if !use_fates
        i_cwd_local = i_oli_mic + 1
        cascade_con.floating_cn_ratio_decomp_pools[i_cwd_local] = true
        cascade_con.decomp_pool_name_restart[i_cwd_local]       = "cwd"
        cascade_con.decomp_pool_name_history[i_cwd_local]       = "CWD"
        cascade_con.decomp_pool_name_long[i_cwd_local]          = "coarse woody debris"
        cascade_con.decomp_pool_name_short[i_cwd_local]         = "CWD"
        cascade_con.is_litter[i_cwd_local]    = false
        cascade_con.is_soil[i_cwd_local]      = false
        cascade_con.is_cwd[i_cwd_local]       = true
        cascade_con.initial_cn_ratio[i_cwd_local] = 10.0
        cascade_con.initial_stock[i_cwd_local]    = params.mimics_initial_Cstocks[i_cwd_local]
        cascade_con.is_metabolic[i_cwd_local] = false
        cascade_con.is_cellulose[i_cwd_local] = false
        cascade_con.is_lignin[i_cwd_local]    = false
    end

    speedup_fac = 1.0

    # Spinup factors
    cascade_con.spinup_factor[i_met_lit] = 1.0
    cascade_con.spinup_factor[i_str_lit] = 1.0
    if !use_fates
        cascade_con.spinup_factor[i_cwd_local] = smooth_max(1.0, speedup_fac * cn_params.tau_cwd * 0.5)
    end
    cascade_con.spinup_factor[i_avl_som]  = 1.0
    cascade_con.spinup_factor[i_chem_som] = 1.0
    cascade_con.spinup_factor[i_phys_som] = 1.0
    cascade_con.spinup_factor[i_cop_mic]  = 1.0
    cascade_con.spinup_factor[i_oli_mic]  = 1.0

    # --- List of transitions and their time-independent coefficients ---

    i_l1m1 = 1
    mimics_state.i_l1m1 = i_l1m1
    cascade_con.cascade_step_name[i_l1m1] = "L1M1"
    cascade_con.cascade_donor_pool[i_l1m1]    = i_met_lit
    cascade_con.cascade_receiver_pool[i_l1m1] = i_cop_mic

    i_l1m2 = 2
    mimics_state.i_l1m2 = i_l1m2
    cascade_con.cascade_step_name[i_l1m2] = "L1M2"
    cascade_con.cascade_donor_pool[i_l1m2]    = i_met_lit
    cascade_con.cascade_receiver_pool[i_l1m2] = i_oli_mic

    i_l2m1 = 3
    mimics_state.i_l2m1 = i_l2m1
    cascade_con.cascade_step_name[i_l2m1] = "L2M1"
    cascade_con.cascade_donor_pool[i_l2m1]    = i_str_lit
    cascade_con.cascade_receiver_pool[i_l2m1] = i_cop_mic

    i_l2m2 = 4
    mimics_state.i_l2m2 = i_l2m2
    cascade_con.cascade_step_name[i_l2m2] = "L2M2"
    cascade_con.cascade_donor_pool[i_l2m2]    = i_str_lit
    cascade_con.cascade_receiver_pool[i_l2m2] = i_oli_mic

    i_s1m1 = 5
    mimics_state.i_s1m1 = i_s1m1
    cascade_con.cascade_step_name[i_s1m1] = "S1M1"
    cascade_con.cascade_donor_pool[i_s1m1]    = i_avl_som
    cascade_con.cascade_receiver_pool[i_s1m1] = i_cop_mic

    i_s1m2 = 6
    mimics_state.i_s1m2 = i_s1m2
    cascade_con.cascade_step_name[i_s1m2] = "S1M2"
    cascade_con.cascade_donor_pool[i_s1m2]    = i_avl_som
    cascade_con.cascade_receiver_pool[i_s1m2] = i_oli_mic

    i_s2s1 = 7
    mimics_state.i_s2s1 = i_s2s1
    cascade_con.cascade_step_name[i_s2s1] = "S2S1"
    cascade_con.cascade_donor_pool[i_s2s1]    = i_chem_som
    cascade_con.cascade_receiver_pool[i_s2s1] = i_avl_som

    i_s3s1 = 8
    mimics_state.i_s3s1 = i_s3s1
    cascade_con.cascade_step_name[i_s3s1] = "S3S1"
    cascade_con.cascade_donor_pool[i_s3s1]    = i_phys_som
    cascade_con.cascade_receiver_pool[i_s3s1] = i_avl_som

    i_m1s1 = 9
    mimics_state.i_m1s1 = i_m1s1
    cascade_con.cascade_step_name[i_m1s1] = "M1S1"
    cascade_con.cascade_donor_pool[i_m1s1]    = i_cop_mic
    cascade_con.cascade_receiver_pool[i_m1s1] = i_avl_som

    i_m1s2 = 10
    mimics_state.i_m1s2 = i_m1s2
    cascade_con.cascade_step_name[i_m1s2] = "M1S2"
    cascade_con.cascade_donor_pool[i_m1s2]    = i_cop_mic
    cascade_con.cascade_receiver_pool[i_m1s2] = i_chem_som

    i_m1s3 = 11
    mimics_state.i_m1s3 = i_m1s3
    cascade_con.cascade_step_name[i_m1s3] = "M1S3"
    cascade_con.cascade_donor_pool[i_m1s3]    = i_cop_mic
    cascade_con.cascade_receiver_pool[i_m1s3] = i_phys_som

    i_m2s1 = 12
    mimics_state.i_m2s1 = i_m2s1
    cascade_con.cascade_step_name[i_m2s1] = "M2S1"
    cascade_con.cascade_donor_pool[i_m2s1]    = i_oli_mic
    cascade_con.cascade_receiver_pool[i_m2s1] = i_avl_som

    i_m2s2 = 13
    mimics_state.i_m2s2 = i_m2s2
    cascade_con.cascade_step_name[i_m2s2] = "M2S2"
    cascade_con.cascade_donor_pool[i_m2s2]    = i_oli_mic
    cascade_con.cascade_receiver_pool[i_m2s2] = i_chem_som

    i_m2s3 = 14
    mimics_state.i_m2s3 = i_m2s3
    cascade_con.cascade_step_name[i_m2s3] = "M2S3"
    cascade_con.cascade_donor_pool[i_m2s3]    = i_oli_mic
    cascade_con.cascade_receiver_pool[i_m2s3] = i_phys_som

    # Set NUE for all transitions
    nue_decomp_cascade = zeros(FT, ndecomp_cascade_transitions_max)
    nue_decomp_cascade[i_l1m1] = mimics_nue_into_mic
    nue_decomp_cascade[i_l1m2] = mimics_nue_into_mic
    nue_decomp_cascade[i_l2m1] = mimics_nue_into_mic
    nue_decomp_cascade[i_l2m2] = mimics_nue_into_mic
    nue_decomp_cascade[i_s1m1] = mimics_nue_into_mic
    nue_decomp_cascade[i_s1m2] = mimics_nue_into_mic
    nue_decomp_cascade[i_s2s1] = 1.0
    nue_decomp_cascade[i_s3s1] = 1.0
    nue_decomp_cascade[i_m1s1] = 1.0
    nue_decomp_cascade[i_m1s2] = 1.0
    nue_decomp_cascade[i_m1s3] = 1.0
    nue_decomp_cascade[i_m2s1] = 1.0
    nue_decomp_cascade[i_m2s2] = 1.0
    nue_decomp_cascade[i_m2s3] = 1.0

    i_cwdl2_local = 0
    if !use_fates
        i_cwdl2_local = 15
        cascade_con.cascade_step_name[i_cwdl2_local] = "CWDL2"
        cascade_con.cascade_donor_pool[i_cwdl2_local]    = i_cwd_local
        cascade_con.cascade_receiver_pool[i_cwdl2_local] = i_str_lit
        nue_decomp_cascade[i_cwdl2_local] = 1.0
    end

    return nue_decomp_cascade
end

# ---------------------------------------------------------------------------
# decomp_rates_mimics! — Calculate decomposition rates
# Ported from decomp_rates_mimics in
# SoilBiogeochemDecompCascadeMIMICSMod.F90
# ---------------------------------------------------------------------------

"""
    decomp_rates_mimics!(cf, mimics_state, params, cn_params, cascade_con;
        mask_bgc_soilc, bounds, nlevdecomp, t_soisno, soilpsi,
        decomp_cpools_vr, col_dz, ligninNratioAvg, annsum_npp_col,
        days_per_year, dt, ...)

Calculate rates and decomposition pathways for the MIMICS decomposition
cascade model.

Ported from `decomp_rates_mimics` in
`SoilBiogeochemDecompCascadeMIMICSMod.F90`.

# Arguments
- `cf::SoilBiogeochemCarbonFluxData`: carbon flux data (modified in place)
- `mimics_state::DecompMIMICSState`: MIMICS module state
- `params::DecompMIMICSParams`: MIMICS parameters
- `cn_params::CNSharedParamsData`: shared CN parameters
- `cascade_con::DecompCascadeConData`: cascade configuration

# Keyword Arguments
- `mask_bgc_soilc::BitVector`: mask for BGC soil columns
- `bounds::UnitRange{Int}`: column bounds
- `nlevdecomp::Int`: number of decomposition levels
- `t_soisno::Matrix{<:Real}`: soil temperature (K), (col × nlev)
- `soilpsi::Matrix{<:Real}`: soil water potential (MPa), (col × nlev)
- `decomp_cpools_vr::Array{<:Real,3}`: decomposing C pools (gC/m3), (col × nlev × npool)
- `col_dz::Matrix{<:Real}`: column layer thicknesses (m), (col × nlev)
- `ligninNratioAvg::Vector{<:Real}`: average lignin C:N ratio per column
- `annsum_npp_col::Vector{<:Real}`: annual sum of NPP at column level (gC/m2/s)
- `days_per_year::Float64`: days per year
- `dt::Float64`: timestep (seconds)
- `spinup_state::Int`: spinup state (0=normal)
- `use_lch4::Bool`: whether LCH4 model is active
- `anoxia::Bool`: whether anoxia is enabled
- `use_fates::Bool`: whether FATES is enabled
- `o2stress_unsat::Matrix{<:Real}`: O2 stress ratio (col × nlev)
- `col_gridcell::Vector{Int}`: column-to-gridcell mapping
- `latdeg::Vector{<:Real}`: gridcell latitudes (degrees)
"""
function decomp_rates_mimics!(cf::SoilBiogeochemCarbonFluxData,
                               mimics_state::DecompMIMICSState,
                               params::DecompMIMICSParams,
                               cn_params::CNSharedParamsData,
                               cascade_con::DecompCascadeConData;
                               mask_bgc_soilc::BitVector,
                               bounds::UnitRange{Int},
                               nlevdecomp::Int,
                               t_soisno::Matrix{<:Real},
                               soilpsi::Matrix{<:Real},
                               decomp_cpools_vr::Array{<:Real,3},
                               col_dz::Matrix{<:Real},
                               ligninNratioAvg::Vector{<:Real},
                               annsum_npp_col::Vector{<:Real},
                               days_per_year::Real,
                               dt::Real,
                               spinup_state::Int=0,
                               use_lch4::Bool=false,
                               anoxia::Bool=false,
                               use_fates::Bool=false,
                               o2stress_unsat::Matrix{<:Real}=Matrix{Float64}(undef, 0, 0),
                               col_gridcell::Vector{Int}=Int[],
                               latdeg::Vector{<:Real}=Float64[])

    eps_val = 1.0e-6
    nc = length(bounds)

    # Unpack shared CN parameters
    rf_cwdl2 = cn_params.rf_cwdl2
    minpsi   = cn_params.minpsi
    maxpsi   = cn_params.maxpsi
    mino2lim = cn_params.mino2lim
    nlev_soildecomp_standard = cn_params.nlev_soildecomp_standard

    # Translate CWD fragmentation to per-second rate constant
    k_frag = 1.0 / (SECSPDAY * days_per_year * cn_params.tau_cwd)

    # Unpack pool/transition indices
    i_met_lit  = mimics_state.i_met_lit
    i_str_lit  = mimics_state.i_str_lit
    i_avl_som  = mimics_state.i_avl_som
    i_chem_som = mimics_state.i_chem_som
    i_phys_som = mimics_state.i_phys_som
    i_cop_mic  = mimics_state.i_cop_mic
    i_oli_mic  = mimics_state.i_oli_mic
    i_l1m1     = mimics_state.i_l1m1
    i_l1m2     = mimics_state.i_l1m2
    i_l2m1     = mimics_state.i_l2m1
    i_l2m2     = mimics_state.i_l2m2
    i_s1m1     = mimics_state.i_s1m1
    i_s1m2     = mimics_state.i_s1m2
    i_m1s1     = mimics_state.i_m1s1
    i_m1s2     = mimics_state.i_m1s2
    i_m2s1     = mimics_state.i_m2s1
    i_m2s2     = mimics_state.i_m2s2
    i_s2s1     = mimics_state.i_s2s1
    i_s3s1     = mimics_state.i_s3s1
    i_m1s3     = mimics_state.i_m1s3
    i_m2s3     = mimics_state.i_m2s3

    spinup_factor = cascade_con.spinup_factor

    # --- Compute spinup geographic terms ---
    FT = eltype(t_soisno)
    spinup_geogterm_l1  = ones(FT, nc)
    spinup_geogterm_l2  = ones(FT, nc)
    spinup_geogterm_cwd = ones(FT, nc)
    spinup_geogterm_s1  = ones(FT, nc)
    spinup_geogterm_s2  = ones(FT, nc)
    spinup_geogterm_s3  = ones(FT, nc)
    spinup_geogterm_m1  = ones(FT, nc)
    spinup_geogterm_m2  = ones(FT, nc)

    if spinup_state >= 1
        for c in 1:nc
            mask_bgc_soilc[c] || continue
            lat = latdeg[col_gridcell[c]]

            if abs(spinup_factor[i_met_lit] - 1.0) > eps_val
                spinup_geogterm_l1[c] = spinup_factor[i_met_lit] * get_spinup_latitude_term(lat)
            end
            if abs(spinup_factor[i_str_lit] - 1.0) > eps_val
                spinup_geogterm_l2[c] = spinup_factor[i_str_lit] * get_spinup_latitude_term(lat)
            end
            if !use_fates
                i_cwd_local = i_oli_mic + 1
                if abs(spinup_factor[i_cwd_local] - 1.0) > eps_val
                    spinup_geogterm_cwd[c] = spinup_factor[i_cwd_local] * get_spinup_latitude_term(lat)
                end
            end
            if abs(spinup_factor[i_avl_som] - 1.0) > eps_val
                spinup_geogterm_s1[c] = spinup_factor[i_avl_som] * get_spinup_latitude_term(lat)
            end
            if abs(spinup_factor[i_chem_som] - 1.0) > eps_val
                spinup_geogterm_s2[c] = spinup_factor[i_chem_som] * get_spinup_latitude_term(lat)
            end
            if abs(spinup_factor[i_phys_som] - 1.0) > eps_val
                spinup_geogterm_s3[c] = spinup_factor[i_phys_som] * get_spinup_latitude_term(lat)
            end
            if abs(spinup_factor[i_cop_mic] - 1.0) > eps_val
                spinup_geogterm_m1[c] = spinup_factor[i_cop_mic] * get_spinup_latitude_term(lat)
            end
            if abs(spinup_factor[i_oli_mic] - 1.0) > eps_val
                spinup_geogterm_m2[c] = spinup_factor[i_oli_mic] * get_spinup_latitude_term(lat)
            end
        end
    end

    # --- Time-dependent coefficients ---
    w_scalar                = cf.w_scalar_col
    o_scalar                = cf.o_scalar_col
    decomp_k                = cf.decomp_k_col
    pathfrac_decomp_cascade = cf.pathfrac_decomp_cascade_col
    rf_decomp_cascade       = cf.rf_decomp_cascade_col
    cn_col                  = cf.cn_col

    depth_scalar = ones(FT, nc, nlevdecomp)

    if nlevdecomp == 1
        # --- Single-level decomposition ---
        frw = zeros(FT, nc)
        fr  = zeros(FT, nc, nlev_soildecomp_standard)

        for j in 1:nlev_soildecomp_standard
            for c in 1:nc
                mask_bgc_soilc[c] || continue
                frw[c] += col_dz[c, j]
            end
        end
        for j in 1:nlev_soildecomp_standard
            for c in 1:nc
                mask_bgc_soilc[c] || continue
                if frw[c] != 0.0
                    fr[c, j] = col_dz[c, j] / frw[c]
                else
                    fr[c, j] = 0.0
                end
            end
        end

        # Water scalar
        for j in 1:nlev_soildecomp_standard
            for c in 1:nc
                mask_bgc_soilc[c] || continue
                if j == 1
                    w_scalar[c, :] .= 0.0
                end
                psi = smooth_min(soilpsi[c, j], maxpsi)
                if psi > minpsi
                    w_scalar[c, 1] += (log(minpsi / psi) / log(minpsi / maxpsi)) * fr[c, j]
                end
            end
        end

        # O2 scalar
        if anoxia
            for j in 1:nlev_soildecomp_standard
                for c in 1:nc
                    mask_bgc_soilc[c] || continue
                    if j == 1
                        o_scalar[c, :] .= 0.0
                    end
                    o_scalar[c, 1] += fr[c, j] * smooth_max(o2stress_unsat[c, j], mino2lim)
                end
            end
        else
            for c in 1:nc
                for j in 1:nlevdecomp
                    o_scalar[c, j] = 1.0
                end
            end
        end

    else
        # --- Multi-level decomposition ---

        # Water scalar
        for j in 1:nlevdecomp
            for c in 1:nc
                mask_bgc_soilc[c] || continue
                psi = smooth_min(soilpsi[c, j], maxpsi)
                if psi > minpsi
                    w_scalar[c, j] = log(minpsi / psi) / log(minpsi / maxpsi)
                else
                    w_scalar[c, j] = 0.0
                end
            end
        end

        # O2 scalar
        if anoxia
            for j in 1:nlevdecomp
                for c in 1:nc
                    mask_bgc_soilc[c] || continue
                    o_scalar[c, j] = smooth_max(o2stress_unsat[c, j], mino2lim)
                end
            end
        else
            for c in 1:nc
                for j in 1:nlevdecomp
                    o_scalar[c, j] = 1.0
                end
            end
        end
    end

    # Depth scalar — currently set to 1.0 (placeholder, as in Fortran)
    for j in 1:nlevdecomp
        for c in 1:nc
            mask_bgc_soilc[c] || continue
            depth_scalar[c, j] = 1.0
        end
    end

    # Unpack MIMICS parameters for decomp rate calculation
    mimics_fmet_p1       = params.mimics_fmet[1]
    mimics_fmet_p2       = params.mimics_fmet[2]
    mimics_fmet_p3       = params.mimics_fmet[3]
    mimics_fmet_p4       = params.mimics_fmet[4]
    mimics_fchem_r_p1    = params.mimics_fchem_r[1]
    mimics_fchem_r_p2    = params.mimics_fchem_r[2]
    mimics_fchem_k_p1    = params.mimics_fchem_k[1]
    mimics_fchem_k_p2    = params.mimics_fchem_k[2]
    mimics_tau_mod_min   = params.mimics_tau_mod_min
    mimics_tau_mod_max   = params.mimics_tau_mod_max
    mimics_tau_mod_factor = params.mimics_tau_mod_factor
    mimics_tau_r_p1      = params.mimics_tau_r[1]
    mimics_tau_r_p2      = params.mimics_tau_r[2]
    mimics_tau_k_p1      = params.mimics_tau_k[1]
    mimics_tau_k_p2      = params.mimics_tau_k[2]
    mimics_ko_r          = params.mimics_ko_r
    mimics_ko_k          = params.mimics_ko_k
    mimics_densdep       = params.mimics_densdep
    mimics_desorpQ10     = params.mimics_desorpQ10
    mimics_t_soi_ref     = params.mimics_t_soi_ref
    mimics_cn_mod_num    = params.mimics_cn_mod_num
    mimics_cn_r          = params.mimics_cn_r
    mimics_cn_k          = params.mimics_cn_k

    # Unpack regression coefficients from state
    vslope_l1_m1 = mimics_state.vslope_l1_m1
    vslope_l2_m1 = mimics_state.vslope_l2_m1
    vslope_s1_m1 = mimics_state.vslope_s1_m1
    vslope_l1_m2 = mimics_state.vslope_l1_m2
    vslope_l2_m2 = mimics_state.vslope_l2_m2
    vslope_s1_m2 = mimics_state.vslope_s1_m2
    kslope_l1_m1 = mimics_state.kslope_l1_m1
    kslope_l2_m1 = mimics_state.kslope_l2_m1
    kslope_s1_m1 = mimics_state.kslope_s1_m1
    kslope_l1_m2 = mimics_state.kslope_l1_m2
    kslope_l2_m2 = mimics_state.kslope_l2_m2
    kslope_s1_m2 = mimics_state.kslope_s1_m2
    vint_l1_m1   = mimics_state.vint_l1_m1
    vint_l2_m1   = mimics_state.vint_l2_m1
    vint_s1_m1   = mimics_state.vint_s1_m1
    vint_l1_m2   = mimics_state.vint_l1_m2
    vint_l2_m2   = mimics_state.vint_l2_m2
    vint_s1_m2   = mimics_state.vint_s1_m2
    kint_l1_m1   = mimics_state.kint_l1_m1
    kint_l2_m1   = mimics_state.kint_l2_m1
    kint_s1_m1   = mimics_state.kint_s1_m1
    kint_l1_m2   = mimics_state.kint_l1_m2
    kint_l2_m2   = mimics_state.kint_l2_m2
    kint_s1_m2   = mimics_state.kint_s1_m2
    vmod_l1_m1   = mimics_state.vmod_l1_m1
    vmod_l2_m1   = mimics_state.vmod_l2_m1
    vmod_s1_m1   = mimics_state.vmod_s1_m1
    vmod_l1_m2   = mimics_state.vmod_l1_m2
    vmod_l2_m2   = mimics_state.vmod_l2_m2
    vmod_s1_m2   = mimics_state.vmod_s1_m2
    kmod_l1_m1   = mimics_state.kmod_l1_m1
    kmod_l2_m1   = mimics_state.kmod_l2_m1
    kmod_s1_m1   = mimics_state.kmod_s1_m1
    kmod_l1_m2   = mimics_state.kmod_l1_m2
    kmod_l2_m2   = mimics_state.kmod_l2_m2
    kmod_s1_m2   = mimics_state.kmod_s1_m2

    fphys_m1 = mimics_state.fphys_m1
    fphys_m2 = mimics_state.fphys_m2
    desorp_arr = mimics_state.desorp
    p_scalar_arr = mimics_state.p_scalar

    # --- Calculate rates for all litter and SOM pools ---
    for c in 1:nc
        mask_bgc_soilc[c] || continue

        annsum_npp_col_scalar = smooth_max(0.0, annsum_npp_col[c])

        # fmet: fraction metabolic litter
        fmet = mimics_fmet_p1 * (mimics_fmet_p2 - mimics_fmet_p3 *
            smooth_min(mimics_fmet_p4, ligninNratioAvg[c]))
        tau_mod = smooth_min(mimics_tau_mod_max, smooth_max(mimics_tau_mod_min,
            sqrt(mimics_tau_mod_factor * annsum_npp_col_scalar)))

        # tau_m1 is tauR and tau_m2 is tauK in Wieder et al. 2015
        # Convert from per hour to per second
        tau_m1 = mimics_tau_r_p1 * exp(mimics_tau_r_p2 * fmet) * tau_mod / SECSPHR
        tau_m2 = mimics_tau_k_p1 * exp(mimics_tau_k_p2 * fmet) * tau_mod / SECSPHR

        # Microbial C:N ratios
        cn_col[c, i_cop_mic] = mimics_cn_r * sqrt(mimics_cn_mod_num / fmet)
        cn_col[c, i_oli_mic] = mimics_cn_k * sqrt(mimics_cn_mod_num / fmet)

        # Chemical fractionation
        fchem_m1 = smooth_min(1.0, smooth_max(0.0, mimics_fchem_r_p1 *
                    exp(mimics_fchem_r_p2 * fmet)))
        fchem_m2 = smooth_min(1.0, smooth_max(0.0, mimics_fchem_k_p1 *
                    exp(mimics_fchem_k_p2 * fmet)))

        for j in 1:nlevdecomp
            # Temperature-dependent Vmax and Km
            t_soi_degC = t_soisno[c, j] - TFRZ

            # Vmax values (convert from per hour to per second)
            vmax_l1_m1 = exp(vslope_l1_m1 * t_soi_degC + vint_l1_m1) * vmod_l1_m1 / SECSPHR
            vmax_l1_m2 = exp(vslope_l1_m2 * t_soi_degC + vint_l1_m2) * vmod_l1_m2 / SECSPHR
            vmax_l2_m1 = exp(vslope_l2_m1 * t_soi_degC + vint_l2_m1) * vmod_l2_m1 / SECSPHR
            vmax_l2_m2 = exp(vslope_l2_m2 * t_soi_degC + vint_l2_m2) * vmod_l2_m2 / SECSPHR
            vmax_s1_m1 = exp(vslope_s1_m1 * t_soi_degC + vint_s1_m1) * vmod_s1_m1 / SECSPHR
            vmax_s1_m2 = exp(vslope_s1_m2 * t_soi_degC + vint_s1_m2) * vmod_s1_m2 / SECSPHR

            # Km values
            km_l1_m1 = exp(kslope_l1_m1 * t_soi_degC + kint_l1_m1) * kmod_l1_m1
            km_l1_m2 = exp(kslope_l1_m2 * t_soi_degC + kint_l1_m2) * kmod_l1_m2
            km_l2_m1 = exp(kslope_l2_m1 * t_soi_degC + kint_l2_m1) * kmod_l2_m1
            km_l2_m2 = exp(kslope_l2_m2 * t_soi_degC + kint_l2_m2) * kmod_l2_m2
            km_s1_m1 = exp(kslope_s1_m1 * t_soi_degC + kint_s1_m1) * kmod_s1_m1 * p_scalar_arr[c, j]
            km_s1_m2 = exp(kslope_s1_m2 * t_soi_degC + kint_s1_m2) * kmod_s1_m2 * p_scalar_arr[c, j]

            # Desorption: function of soil temperature and Q10
            desorption = (desorp_arr[c, j] / SECSPHR) * mimics_desorpQ10 *
                         exp((t_soi_degC - mimics_t_soi_ref) / 10.0)

            # Microbial concentration (mgC/cm3)
            m1_conc = (decomp_cpools_vr[c, j, i_cop_mic] / col_dz[c, j]) * G_TO_MG * CM3_TO_M3
            m2_conc = (decomp_cpools_vr[c, j, i_oli_mic] / col_dz[c, j]) * G_TO_MG * CM3_TO_M3

            # Product of w_scalar * depth_scalar * o_scalar
            w_d_o_scalars = w_scalar[c, j] * depth_scalar[c, j] * o_scalar[c, j]

            # Metabolic litter decomposition
            term_1 = vmax_l1_m1 * m1_conc / (km_l1_m1 + m1_conc)
            term_2 = vmax_l1_m2 * m2_conc / (km_l1_m2 + m2_conc)
            decomp_k[c, j, i_met_lit] = (term_1 + term_2) * w_d_o_scalars
            if (term_1 + term_2) != 0.0
                pathfrac_decomp_cascade[c, j, i_l1m1] = term_1 / (term_1 + term_2)
                pathfrac_decomp_cascade[c, j, i_l1m2] = term_2 / (term_1 + term_2)
            else
                pathfrac_decomp_cascade[c, j, i_l1m1] = 0.0
                pathfrac_decomp_cascade[c, j, i_l1m2] = 0.0
            end

            # Structural litter decomposition
            term_1 = vmax_l2_m1 * m1_conc / (km_l2_m1 + m1_conc)
            term_2 = vmax_l2_m2 * m2_conc / (km_l2_m2 + m2_conc)
            decomp_k[c, j, i_str_lit] = (term_1 + term_2) * w_d_o_scalars
            if (term_1 + term_2) != 0.0
                pathfrac_decomp_cascade[c, j, i_l2m1] = term_1 / (term_1 + term_2)
                pathfrac_decomp_cascade[c, j, i_l2m2] = term_2 / (term_1 + term_2)
            else
                pathfrac_decomp_cascade[c, j, i_l2m1] = 0.0
                pathfrac_decomp_cascade[c, j, i_l2m2] = 0.0
            end

            # Available SOM decomposition
            term_1 = vmax_s1_m1 * m1_conc / (km_s1_m1 + m1_conc)
            term_2 = vmax_s1_m2 * m2_conc / (km_s1_m2 + m2_conc)
            decomp_k[c, j, i_avl_som] = (term_1 + term_2) * w_d_o_scalars
            if (term_1 + term_2) != 0.0
                pathfrac_decomp_cascade[c, j, i_s1m1] = term_1 / (term_1 + term_2)
                pathfrac_decomp_cascade[c, j, i_s1m2] = term_2 / (term_1 + term_2)
            else
                pathfrac_decomp_cascade[c, j, i_s1m1] = 0.0
                pathfrac_decomp_cascade[c, j, i_s1m2] = 0.0
            end

            # Physically protected SOM: desorption only
            decomp_k[c, j, i_phys_som] = desorption * depth_scalar[c, j]

            # Chemically protected SOM: oxidation
            term_1 = vmax_l2_m1 * m1_conc / (mimics_ko_r * km_l2_m1 + m1_conc)
            term_2 = vmax_l2_m2 * m2_conc / (mimics_ko_k * km_l2_m2 + m2_conc)
            decomp_k[c, j, i_chem_som] = (term_1 + term_2) * w_d_o_scalars

            # Copiotrophic microbe turnover (M1)
            decomp_k[c, j, i_cop_mic] = tau_m1 *
                   m1_conc^(mimics_densdep - 1.0) * w_d_o_scalars
            favl = smooth_min(1.0, smooth_max(0.0, 1.0 - fphys_m1[c, j] - fchem_m1))
            pathfrac_decomp_cascade[c, j, i_m1s1] = favl
            pathfrac_decomp_cascade[c, j, i_m1s2] = fchem_m1

            # Oligotrophic microbe turnover (M2)
            decomp_k[c, j, i_oli_mic] = tau_m2 *
                   m2_conc^(mimics_densdep - 1.0) * w_d_o_scalars
            favl = smooth_min(1.0, smooth_max(0.0, 1.0 - fphys_m2[c, j] - fchem_m2))
            pathfrac_decomp_cascade[c, j, i_m2s1] = favl
            pathfrac_decomp_cascade[c, j, i_m2s2] = fchem_m2

            # CWD fragmentation (only if FATES not enabled)
            if !use_fates
                i_cwd_local = i_oli_mic + 1
                decomp_k[c, j, i_cwd_local] = k_frag * w_d_o_scalars
            end
        end
    end

    # pathfrac terms not calculated in the previous loop
    for j in 1:nlevdecomp
        for c in 1:nc
            pathfrac_decomp_cascade[c, j, i_s2s1] = 1.0
            pathfrac_decomp_cascade[c, j, i_s3s1] = 1.0
            pathfrac_decomp_cascade[c, j, i_m1s3] = fphys_m1[c, j]
            pathfrac_decomp_cascade[c, j, i_m2s3] = fphys_m2[c, j]
        end
    end
    if !use_fates
        i_cwdl2_local = 15
        for j in 1:nlevdecomp
            for c in 1:nc
                pathfrac_decomp_cascade[c, j, i_cwdl2_local] = 1.0
                rf_decomp_cascade[c, j, i_cwdl2_local] = rf_cwdl2
            end
        end
    end

    # Respiration fractions
    for j in 1:nlevdecomp
        for c in 1:nc
            rf_decomp_cascade[c, j, i_l1m1] = mimics_state.rf_l1m1
            rf_decomp_cascade[c, j, i_l1m2] = mimics_state.rf_l1m2
            rf_decomp_cascade[c, j, i_l2m1] = mimics_state.rf_l2m1
            rf_decomp_cascade[c, j, i_l2m2] = mimics_state.rf_l2m2
            rf_decomp_cascade[c, j, i_s1m1] = mimics_state.rf_s1m1
            rf_decomp_cascade[c, j, i_s1m2] = mimics_state.rf_s1m2
            rf_decomp_cascade[c, j, i_s2s1] = 0.0
            rf_decomp_cascade[c, j, i_s3s1] = 0.0
            rf_decomp_cascade[c, j, i_m1s1] = 0.0
            rf_decomp_cascade[c, j, i_m1s2] = 0.0
            rf_decomp_cascade[c, j, i_m1s3] = 0.0
            rf_decomp_cascade[c, j, i_m2s1] = 0.0
            rf_decomp_cascade[c, j, i_m2s2] = 0.0
            rf_decomp_cascade[c, j, i_m2s3] = 0.0
        end
    end

    return nothing
end
