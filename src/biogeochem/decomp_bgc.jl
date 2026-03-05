# ==========================================================================
# Ported from: src/soilbiogeochem/SoilBiogeochemDecompCascadeBGCMod.F90
# Sets the coefficients used in the decomposition cascade submodel.
# This uses the CENTURY/BGC parameters.
#
# Public functions:
#   decomp_bgc_read_params!       — Read BGC decomposition parameters
#   init_decomp_cascade_bgc!      — Initialize rate constants and pathways
#   decomp_rate_constants_bgc!    — Calculate decomposition rate constants
# ==========================================================================

# ---------------------------------------------------------------------------
# DecompCascadeConData — cascade configuration (pools, transitions, names)
# Ported from decomp_cascade_con in SoilBiogeochemDecompCascadeConType.F90
# ---------------------------------------------------------------------------

"""
    DecompCascadeConData

Decomposition cascade configuration. Holds pool properties, transition
definitions, and spinup factors for the BGC decomposition cascade.

Ported from `decomp_cascade_con` in `SoilBiogeochemDecompCascadeConType.F90`.
"""
Base.@kwdef mutable struct DecompCascadeConData{FT<:AbstractFloat}
    cascade_donor_pool             ::Vector{Int}     = Int[]
    cascade_receiver_pool          ::Vector{Int}     = Int[]
    floating_cn_ratio_decomp_pools ::BitVector        = BitVector()
    is_litter                      ::BitVector        = BitVector()
    is_soil                        ::BitVector        = BitVector()
    is_cwd                         ::BitVector        = BitVector()
    initial_cn_ratio               ::Vector{FT} = Float64[]
    initial_stock                  ::Vector{FT} = Float64[]
    initial_stock_soildepth        ::FT           = 0.3
    is_metabolic                   ::BitVector        = BitVector()
    is_cellulose                   ::BitVector        = BitVector()
    is_lignin                      ::BitVector        = BitVector()
    spinup_factor                  ::Vector{FT} = Float64[]
    decomp_pool_name_restart       ::Vector{String}   = String[]
    decomp_pool_name_history       ::Vector{String}   = String[]
    decomp_pool_name_long          ::Vector{String}   = String[]
    decomp_pool_name_short         ::Vector{String}   = String[]
    cascade_step_name              ::Vector{String}   = String[]
end

# ---------------------------------------------------------------------------
# DecompBGCParams — BGC decomposition parameters
# Ported from params_type in SoilBiogeochemDecompCascadeBGCMod.F90
# ---------------------------------------------------------------------------

"""
    DecompBGCParams

BGC decomposition cascade parameters. Holds C:N ratios, respiration
fractions, turnover times, cellulose fraction, and initial C stocks.

Ported from `params_type` in `SoilBiogeochemDecompCascadeBGCMod.F90`.
"""
Base.@kwdef mutable struct DecompBGCParams{FT<:AbstractFloat}
    cn_s1_bgc               ::FT          = 12.0
    cn_s2_bgc               ::FT          = 12.0
    cn_s3_bgc               ::FT          = 10.0
    rf_l1s1_bgc             ::FT          = 0.39
    rf_l2s1_bgc             ::FT          = 0.55
    rf_l3s2_bgc             ::FT          = 0.29
    rf_s2s1_bgc             ::FT          = 0.55
    rf_s2s3_bgc             ::FT          = 0.55
    rf_s3s1_bgc             ::FT          = 0.55
    rf_cwdl3_bgc            ::FT          = 0.0
    tau_l1_bgc              ::FT          = 1.0 / 18.5
    tau_l2_l3_bgc           ::FT          = 1.0 / 4.9
    tau_s1_bgc              ::FT          = 1.0 / 7.3
    tau_s2_bgc              ::FT          = 1.0 / 0.2
    tau_s3_bgc              ::FT          = 1.0 / 0.0045
    cwd_fcel_bgc            ::FT          = 0.0
    bgc_initial_Cstocks     ::Vector{FT} = Float64[]
    bgc_initial_Cstocks_depth ::FT        = 0.3
end

# ---------------------------------------------------------------------------
# DecompBGCState — module-level state for BGC decomposition
# Holds pool indices, respiration fractions, path fractions, transition
# indices, and configuration flags set during initialization.
# ---------------------------------------------------------------------------

"""
    DecompBGCState

Module-level persistent state for the BGC decomposition cascade.
Pool indices, respiration fractions, path fractions, and transition indices
are set once during `init_decomp_cascade_bgc!` and used by
`decomp_rate_constants_bgc!`.

Ported from module-level private variables in
`SoilBiogeochemDecompCascadeBGCMod.F90`.
"""
Base.@kwdef mutable struct DecompBGCState{FT<:AbstractFloat}
    # Pool indices
    i_pas_som ::Int = 0
    i_slo_som ::Int = 0
    i_act_som ::Int = 0
    i_cel_lit ::Int = 0
    i_lig_lit ::Int = 0

    # Scalar respiration fractions
    cwd_fcel  ::FT = 0.0
    rf_l1s1   ::FT = 0.0
    rf_l2s1   ::FT = 0.0
    rf_l3s2   ::FT = 0.0
    rf_s2s1   ::FT = 0.0
    rf_s2s3   ::FT = 0.0
    rf_s3s1   ::FT = 0.0
    rf_cwdl3  ::FT = 0.0

    # Spatially-varying respiration fractions (col × nlevdecomp)
    rf_s1s2   ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    rf_s1s3   ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)

    # Path fractions
    f_s1s2    ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    f_s1s3    ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    f_s2s1    ::FT = 0.0
    f_s2s3    ::FT = 0.0

    # Transition indices
    i_l1s1    ::Int = 0
    i_l2s1    ::Int = 0
    i_l3s2    ::Int = 0
    i_s1s2    ::Int = 0
    i_s1s3    ::Int = 0
    i_s2s1    ::Int = 0
    i_s2s3    ::Int = 0
    i_s3s1    ::Int = 0
    i_cwdl3   ::Int = 0

    # Public configuration flags
    normalize_q10_to_century_tfunc ::Bool    = true
    use_century_tfunc              ::Bool    = false
    normalization_tref             ::FT = 15.0
end

# ---------------------------------------------------------------------------
# decomp_bgc_read_params! — Read BGC decomposition parameters
# Ported from readParams in SoilBiogeochemDecompCascadeBGCMod.F90
# ---------------------------------------------------------------------------

"""
    decomp_bgc_read_params!(params; kwargs...)

Populate `DecompBGCParams` from keyword arguments (replaces NetCDF file
reading). Each keyword maps to a parameter variable name from the Fortran
`readParams` subroutine.

Ported from `readParams` in `SoilBiogeochemDecompCascadeBGCMod.F90`.
"""
function decomp_bgc_read_params!(params::DecompBGCParams;
                                  tau_l1::Float64,
                                  tau_l2_l3::Float64,
                                  tau_s1::Float64,
                                  tau_s2::Float64,
                                  tau_s3::Float64,
                                  cn_s1::Float64,
                                  cn_s2::Float64,
                                  cn_s3::Float64,
                                  rf_l1s1::Float64,
                                  rf_l2s1::Float64,
                                  rf_l3s2::Float64,
                                  rf_s2s1::Float64,
                                  rf_s2s3::Float64,
                                  rf_s3s1::Float64,
                                  rf_cwdl3::Float64,
                                  cwd_fcel::Float64,
                                  bgc_initial_Cstocks::Vector{Float64},
                                  bgc_initial_Cstocks_depth::Float64)
    params.tau_l1_bgc    = tau_l1
    params.tau_l2_l3_bgc = tau_l2_l3
    params.tau_s1_bgc    = tau_s1
    params.tau_s2_bgc    = tau_s2
    params.tau_s3_bgc    = tau_s3
    params.cn_s1_bgc     = cn_s1
    params.cn_s2_bgc     = cn_s2
    params.cn_s3_bgc     = cn_s3
    params.rf_l1s1_bgc   = rf_l1s1
    params.rf_l2s1_bgc   = rf_l2s1
    params.rf_l3s2_bgc   = rf_l3s2
    params.rf_s2s1_bgc   = rf_s2s1
    params.rf_s2s3_bgc   = rf_s2s3
    params.rf_s3s1_bgc   = rf_s3s1
    params.rf_cwdl3_bgc  = rf_cwdl3
    params.cwd_fcel_bgc  = cwd_fcel
    params.bgc_initial_Cstocks       = copy(bgc_initial_Cstocks)
    params.bgc_initial_Cstocks_depth = bgc_initial_Cstocks_depth
    return nothing
end

# ---------------------------------------------------------------------------
# init_decomp_cascade_bgc! — Initialize cascade pathways and coefficients
# Ported from init_decompcascade_bgc in SoilBiogeochemDecompCascadeBGCMod.F90
# ---------------------------------------------------------------------------

"""
    init_decomp_cascade_bgc!(bgc_state, cascade_con, params, cn_params;
                              cellsand, bounds, nlevdecomp,
                              ndecomp_pools_max, ndecomp_cascade_transitions_max,
                              spinup_state, use_fates)

Initialize rate constants and decomposition pathways following the
decomposition cascade of the BGC model.

Ported from `init_decompcascade_bgc` in `SoilBiogeochemDecompCascadeBGCMod.F90`.

# Arguments
- `bgc_state::DecompBGCState`: module-level state to populate
- `cascade_con::DecompCascadeConData`: cascade configuration to populate
- `params::DecompBGCParams`: BGC parameters
- `cn_params::CNSharedParamsData`: shared CN parameters

# Keyword Arguments
- `cellsand::Matrix{Float64}`: column sand fraction (col × nlevdecomp)
- `bounds::UnitRange{Int}`: column bounds
- `nlevdecomp::Int`: number of decomposition levels
- `ndecomp_pools_max::Int`: maximum number of decomposition pools
- `ndecomp_cascade_transitions_max::Int`: max number of cascade transitions
- `spinup_state::Int`: spinup state flag (0=normal)
- `use_fates::Bool`: whether FATES is enabled
"""
function init_decomp_cascade_bgc!(bgc_state::DecompBGCState,
                                   cascade_con::DecompCascadeConData,
                                   params::DecompBGCParams,
                                   cn_params::CNSharedParamsData;
                                   cellsand::Matrix{Float64},
                                   bounds::UnitRange{Int},
                                   nlevdecomp::Int,
                                   ndecomp_pools_max::Int,
                                   ndecomp_cascade_transitions_max::Int=10,
                                   spinup_state::Int=0,
                                   use_fates::Bool=false)

    nc = length(bounds)

    # Allocate cascade_con arrays
    cascade_con.cascade_donor_pool             = zeros(Int, ndecomp_cascade_transitions_max)
    cascade_con.cascade_receiver_pool          = zeros(Int, ndecomp_cascade_transitions_max)
    cascade_con.floating_cn_ratio_decomp_pools = falses(ndecomp_pools_max)
    cascade_con.is_litter                      = falses(ndecomp_pools_max)
    cascade_con.is_soil                        = falses(ndecomp_pools_max)
    cascade_con.is_cwd                         = falses(ndecomp_pools_max)
    cascade_con.initial_cn_ratio               = zeros(ndecomp_pools_max)
    cascade_con.initial_stock                  = zeros(ndecomp_pools_max)
    cascade_con.is_metabolic                   = falses(ndecomp_pools_max)
    cascade_con.is_cellulose                   = falses(ndecomp_pools_max)
    cascade_con.is_lignin                      = falses(ndecomp_pools_max)
    cascade_con.spinup_factor                  = ones(ndecomp_pools_max)
    cascade_con.decomp_pool_name_restart       = fill("", ndecomp_pools_max)
    cascade_con.decomp_pool_name_history       = fill("", ndecomp_pools_max)
    cascade_con.decomp_pool_name_long          = fill("", ndecomp_pools_max)
    cascade_con.decomp_pool_name_short         = fill("", ndecomp_pools_max)
    cascade_con.cascade_step_name              = fill("", ndecomp_cascade_transitions_max)

    # Allocate spatially-varying arrays
    bgc_state.rf_s1s2 = zeros(nc, nlevdecomp)
    bgc_state.rf_s1s3 = zeros(nc, nlevdecomp)
    bgc_state.f_s1s2  = zeros(nc, nlevdecomp)
    bgc_state.f_s1s3  = zeros(nc, nlevdecomp)

    # --- Time-constant coefficients ---
    cn_s1 = params.cn_s1_bgc
    cn_s2 = params.cn_s2_bgc
    cn_s3 = params.cn_s3_bgc

    # Set respiration fractions
    bgc_state.rf_l1s1  = params.rf_l1s1_bgc
    bgc_state.rf_l2s1  = params.rf_l2s1_bgc
    bgc_state.rf_l3s2  = params.rf_l3s2_bgc
    bgc_state.rf_s2s1  = params.rf_s2s1_bgc
    bgc_state.rf_s2s3  = params.rf_s2s3_bgc
    bgc_state.rf_s3s1  = params.rf_s3s1_bgc
    bgc_state.rf_cwdl3 = params.rf_cwdl3_bgc

    # Set cellulose fraction for CWD
    bgc_state.cwd_fcel = params.cwd_fcel_bgc

    # Set path fractions
    bgc_state.f_s2s1 = 0.42 / 0.45
    bgc_state.f_s2s3 = 0.03 / 0.45

    # Sand-dependent fractions
    for c in 1:nc
        for j in 1:nlevdecomp
            t = 0.85 - 0.68 * 0.01 * (100.0 - cellsand[c, j])
            bgc_state.f_s1s2[c, j]  = 1.0 - 0.004 / (1.0 - t)
            bgc_state.f_s1s3[c, j]  = 0.004 / (1.0 - t)
            bgc_state.rf_s1s2[c, j] = t
            bgc_state.rf_s1s3[c, j] = t
        end
    end
    cascade_con.initial_stock_soildepth = params.bgc_initial_Cstocks_depth

    # --- List of pools and their attributes ---

    # i_met_lit (metabolic litter)
    i_met_lit = 1
    cascade_con.floating_cn_ratio_decomp_pools[i_met_lit] = true
    cascade_con.decomp_pool_name_restart[i_met_lit]       = "litr1"
    cascade_con.decomp_pool_name_history[i_met_lit]       = "LIT_MET"
    cascade_con.decomp_pool_name_long[i_met_lit]          = "metabolic litter"
    cascade_con.decomp_pool_name_short[i_met_lit]         = "L1"
    cascade_con.is_litter[i_met_lit]    = true
    cascade_con.is_soil[i_met_lit]      = false
    cascade_con.is_cwd[i_met_lit]       = false
    cascade_con.initial_cn_ratio[i_met_lit] = 90.0
    cascade_con.initial_stock[i_met_lit]    = params.bgc_initial_Cstocks[i_met_lit]
    cascade_con.is_metabolic[i_met_lit] = true
    cascade_con.is_cellulose[i_met_lit] = false
    cascade_con.is_lignin[i_met_lit]    = false

    # i_cel_lit (cellulosic litter)
    i_cel_lit = i_met_lit + 1
    bgc_state.i_cel_lit = i_cel_lit
    cascade_con.floating_cn_ratio_decomp_pools[i_cel_lit] = true
    cascade_con.decomp_pool_name_restart[i_cel_lit]       = "litr2"
    cascade_con.decomp_pool_name_history[i_cel_lit]       = "LIT_CEL"
    cascade_con.decomp_pool_name_long[i_cel_lit]          = "cellulosic litter"
    cascade_con.decomp_pool_name_short[i_cel_lit]         = "L2"
    cascade_con.is_litter[i_cel_lit]    = true
    cascade_con.is_soil[i_cel_lit]      = false
    cascade_con.is_cwd[i_cel_lit]       = false
    cascade_con.initial_cn_ratio[i_cel_lit] = 90.0
    cascade_con.initial_stock[i_cel_lit]    = params.bgc_initial_Cstocks[i_cel_lit]
    cascade_con.is_metabolic[i_cel_lit] = false
    cascade_con.is_cellulose[i_cel_lit] = true
    cascade_con.is_lignin[i_cel_lit]    = false

    # i_lig_lit (lignin litter)
    i_lig_lit = i_cel_lit + 1
    bgc_state.i_lig_lit = i_lig_lit
    cascade_con.floating_cn_ratio_decomp_pools[i_lig_lit] = true
    cascade_con.decomp_pool_name_restart[i_lig_lit]       = "litr3"
    cascade_con.decomp_pool_name_history[i_lig_lit]       = "LIT_LIG"
    cascade_con.decomp_pool_name_long[i_lig_lit]          = "lignin litter"
    cascade_con.decomp_pool_name_short[i_lig_lit]         = "L3"
    cascade_con.is_litter[i_lig_lit]    = true
    cascade_con.is_soil[i_lig_lit]      = false
    cascade_con.is_cwd[i_lig_lit]       = false
    cascade_con.initial_cn_ratio[i_lig_lit] = 90.0
    cascade_con.initial_stock[i_lig_lit]    = params.bgc_initial_Cstocks[i_lig_lit]
    cascade_con.is_metabolic[i_lig_lit] = false
    cascade_con.is_cellulose[i_lig_lit] = false
    cascade_con.is_lignin[i_lig_lit]    = true

    # i_act_som (active SOM)
    i_act_som = i_lig_lit + 1
    bgc_state.i_act_som = i_act_som
    cascade_con.floating_cn_ratio_decomp_pools[i_act_som] = false
    cascade_con.decomp_pool_name_restart[i_act_som]       = "soil1"
    cascade_con.decomp_pool_name_history[i_act_som]       = "SOM_ACT"
    cascade_con.decomp_pool_name_long[i_act_som]          = "active soil organic matter"
    cascade_con.decomp_pool_name_short[i_act_som]         = "S1"
    cascade_con.is_litter[i_act_som]    = false
    cascade_con.is_soil[i_act_som]      = true
    cascade_con.is_cwd[i_act_som]       = false
    cascade_con.initial_cn_ratio[i_act_som] = cn_s1
    cascade_con.initial_stock[i_act_som]    = params.bgc_initial_Cstocks[i_act_som]
    cascade_con.is_metabolic[i_act_som] = false
    cascade_con.is_cellulose[i_act_som] = false
    cascade_con.is_lignin[i_act_som]    = false

    # i_slo_som (slow SOM)
    i_slo_som = i_act_som + 1
    bgc_state.i_slo_som = i_slo_som
    cascade_con.floating_cn_ratio_decomp_pools[i_slo_som] = false
    cascade_con.decomp_pool_name_restart[i_slo_som]       = "soil2"
    cascade_con.decomp_pool_name_history[i_slo_som]       = "SOM_SLO"
    cascade_con.decomp_pool_name_long[i_slo_som]          = "slow soil organic matter"
    cascade_con.decomp_pool_name_short[i_slo_som]         = "S2"
    cascade_con.is_litter[i_slo_som]    = false
    cascade_con.is_soil[i_slo_som]      = true
    cascade_con.is_cwd[i_slo_som]       = false
    cascade_con.initial_cn_ratio[i_slo_som] = cn_s2
    cascade_con.initial_stock[i_slo_som]    = params.bgc_initial_Cstocks[i_slo_som]
    cascade_con.is_metabolic[i_slo_som] = false
    cascade_con.is_cellulose[i_slo_som] = false
    cascade_con.is_lignin[i_slo_som]    = false

    # i_pas_som (passive SOM)
    i_pas_som = i_slo_som + 1
    bgc_state.i_pas_som = i_pas_som
    cascade_con.floating_cn_ratio_decomp_pools[i_pas_som] = false
    cascade_con.decomp_pool_name_restart[i_pas_som]       = "soil3"
    cascade_con.decomp_pool_name_history[i_pas_som]       = "SOM_PAS"
    cascade_con.decomp_pool_name_long[i_pas_som]          = "passive soil organic matter"
    cascade_con.decomp_pool_name_short[i_pas_som]         = "S3"
    cascade_con.is_litter[i_pas_som]    = false
    cascade_con.is_soil[i_pas_som]      = true
    cascade_con.is_cwd[i_pas_som]       = false
    cascade_con.initial_cn_ratio[i_pas_som] = cn_s3
    cascade_con.initial_stock[i_pas_som]    = params.bgc_initial_Cstocks[i_pas_som]
    cascade_con.is_metabolic[i_pas_som] = false
    cascade_con.is_cellulose[i_pas_som] = false
    cascade_con.is_lignin[i_pas_som]    = false

    # CWD (only if FATES not enabled)
    i_cwd_local = 0
    if !use_fates
        i_cwd_local = i_pas_som + 1
        cascade_con.floating_cn_ratio_decomp_pools[i_cwd_local] = true
        cascade_con.decomp_pool_name_restart[i_cwd_local]       = "cwd"
        cascade_con.decomp_pool_name_history[i_cwd_local]       = "CWD"
        cascade_con.decomp_pool_name_long[i_cwd_local]          = "coarse woody debris"
        cascade_con.decomp_pool_name_short[i_cwd_local]         = "CWD"
        cascade_con.is_litter[i_cwd_local]    = false
        cascade_con.is_soil[i_cwd_local]      = false
        cascade_con.is_cwd[i_cwd_local]       = true
        cascade_con.initial_cn_ratio[i_cwd_local] = 90.0
        cascade_con.initial_stock[i_cwd_local]    = params.bgc_initial_Cstocks[i_cwd_local]
        cascade_con.is_metabolic[i_cwd_local] = false
        cascade_con.is_cellulose[i_cwd_local] = false
        cascade_con.is_lignin[i_cwd_local]    = false
    end

    speedup_fac = 1.0

    # Spinup factors
    cascade_con.spinup_factor[i_met_lit] = 1.0
    cascade_con.spinup_factor[i_cel_lit] = 1.0
    cascade_con.spinup_factor[i_lig_lit] = 1.0
    if !use_fates
        cascade_con.spinup_factor[i_cwd_local] = max(1.0, speedup_fac * cn_params.tau_cwd / 2.0)
    end
    cascade_con.spinup_factor[i_act_som] = 1.0
    cascade_con.spinup_factor[i_slo_som] = max(1.0, speedup_fac * params.tau_s2_bgc)
    cascade_con.spinup_factor[i_pas_som] = max(1.0, speedup_fac * params.tau_s3_bgc)

    # --- List of transitions and their time-independent coefficients ---
    i_l1s1 = 1
    bgc_state.i_l1s1 = i_l1s1
    cascade_con.cascade_step_name[i_l1s1] = "L1S1"
    cascade_con.cascade_donor_pool[i_l1s1]    = i_met_lit
    cascade_con.cascade_receiver_pool[i_l1s1] = i_act_som

    i_l2s1 = 2
    bgc_state.i_l2s1 = i_l2s1
    cascade_con.cascade_step_name[i_l2s1] = "L2S1"
    cascade_con.cascade_donor_pool[i_l2s1]    = i_cel_lit
    cascade_con.cascade_receiver_pool[i_l2s1] = i_act_som

    i_l3s2 = 3
    bgc_state.i_l3s2 = i_l3s2
    cascade_con.cascade_step_name[i_l3s2] = "L3S2"
    cascade_con.cascade_donor_pool[i_l3s2]    = i_lig_lit
    cascade_con.cascade_receiver_pool[i_l3s2] = i_slo_som

    i_s1s2 = 4
    bgc_state.i_s1s2 = i_s1s2
    cascade_con.cascade_step_name[i_s1s2] = "S1S2"
    cascade_con.cascade_donor_pool[i_s1s2]    = i_act_som
    cascade_con.cascade_receiver_pool[i_s1s2] = i_slo_som

    i_s1s3 = 5
    bgc_state.i_s1s3 = i_s1s3
    cascade_con.cascade_step_name[i_s1s3] = "S1S3"
    cascade_con.cascade_donor_pool[i_s1s3]    = i_act_som
    cascade_con.cascade_receiver_pool[i_s1s3] = i_pas_som

    i_s2s1 = 6
    bgc_state.i_s2s1 = i_s2s1
    cascade_con.cascade_step_name[i_s2s1] = "S2S1"
    cascade_con.cascade_donor_pool[i_s2s1]    = i_slo_som
    cascade_con.cascade_receiver_pool[i_s2s1] = i_act_som

    i_s2s3 = 7
    bgc_state.i_s2s3 = i_s2s3
    cascade_con.cascade_step_name[i_s2s3] = "S2S3"
    cascade_con.cascade_donor_pool[i_s2s3]    = i_slo_som
    cascade_con.cascade_receiver_pool[i_s2s3] = i_pas_som

    i_s3s1 = 8
    bgc_state.i_s3s1 = i_s3s1
    cascade_con.cascade_step_name[i_s3s1] = "S3S1"
    cascade_con.cascade_donor_pool[i_s3s1]    = i_pas_som
    cascade_con.cascade_receiver_pool[i_s3s1] = i_act_som

    i_cwdl2_local = 0
    if !use_fates
        i_cwdl2_local = 9
        cascade_con.cascade_step_name[i_cwdl2_local] = "CWDL2"
        cascade_con.cascade_donor_pool[i_cwdl2_local]    = i_cwd_local
        cascade_con.cascade_receiver_pool[i_cwdl2_local] = i_cel_lit

        i_cwdl3_local = 10
        bgc_state.i_cwdl3 = i_cwdl3_local
        cascade_con.cascade_step_name[i_cwdl3_local] = "CWDL3"
        cascade_con.cascade_donor_pool[i_cwdl3_local]    = i_cwd_local
        cascade_con.cascade_receiver_pool[i_cwdl3_local] = i_lig_lit
    end

    return nothing
end

# ---------------------------------------------------------------------------
# decomp_rate_constants_bgc! — Calculate decomposition rate constants
# Ported from decomp_rate_constants_bgc in SoilBiogeochemDecompCascadeBGCMod.F90
# ---------------------------------------------------------------------------

"""
    decomp_rate_constants_bgc!(cf, bgc_state, params, cn_params, cascade_con;
        mask_bgc_soilc, bounds, nlevdecomp, t_soisno, soilpsi,
        days_per_year, dt, zsoi_vals, ...)

Calculate rate constants and decomposition pathways for the CENTURY
decomposition cascade model.

Ported from `decomp_rate_constants_bgc` in
`SoilBiogeochemDecompCascadeBGCMod.F90`.

# Arguments
- `cf::SoilBiogeochemCarbonFluxData`: carbon flux data (modified in place)
- `bgc_state::DecompBGCState`: BGC module state (pool/transition indices, etc.)
- `params::DecompBGCParams`: BGC parameters
- `cn_params::CNSharedParamsData`: shared CN parameters
- `cascade_con::DecompCascadeConData`: cascade configuration

# Keyword Arguments
- `mask_bgc_soilc::BitVector`: mask for BGC soil columns
- `bounds::UnitRange{Int}`: column bounds
- `nlevdecomp::Int`: number of decomposition levels
- `t_soisno::Matrix{Float64}`: soil temperature (K), (col × nlev)
- `soilpsi::Matrix{Float64}`: soil water potential (MPa), (col × nlev)
- `days_per_year::Float64`: days per year
- `dt::Float64`: timestep (seconds)
- `zsoi_vals::Vector{Float64}`: soil layer depths (m)
- `spinup_state::Int`: spinup state (0=normal, >=1=accelerated)
- `use_lch4::Bool`: whether LCH4 model is active
- `anoxia::Bool`: whether anoxia is enabled
- `use_fates::Bool`: whether FATES is enabled
- `o2stress_unsat::Matrix{Float64}`: O2 stress ratio (col × nlev), for anoxia
- `col_dz::Matrix{Float64}`: column layer thicknesses (col × nlev), for nlevdecomp==1
- `col_gridcell::Vector{Int}`: column-to-gridcell mapping, for spinup
- `latdeg::Vector{Float64}`: gridcell latitudes (degrees), for spinup
"""
function decomp_rate_constants_bgc!(cf::SoilBiogeochemCarbonFluxData,
                                     bgc_state::DecompBGCState,
                                     params::DecompBGCParams,
                                     cn_params::CNSharedParamsData,
                                     cascade_con::DecompCascadeConData;
                                     mask_bgc_soilc::BitVector,
                                     bounds::UnitRange{Int},
                                     nlevdecomp::Int,
                                     t_soisno::Matrix{Float64},
                                     soilpsi::Matrix{Float64},
                                     days_per_year::Float64,
                                     dt::Float64,
                                     zsoi_vals::Vector{Float64},
                                     spinup_state::Int=0,
                                     use_lch4::Bool=false,
                                     anoxia::Bool=false,
                                     use_fates::Bool=false,
                                     o2stress_unsat::Matrix{Float64}=Matrix{Float64}(undef, 0, 0),
                                     col_dz::Matrix{Float64}=Matrix{Float64}(undef, 0, 0),
                                     col_gridcell::Vector{Int}=Int[],
                                     latdeg::Vector{Float64}=Float64[])

    eps_val = 1.0e-6
    nc = length(bounds)

    # CENTURY temperature response function
    catanf(t1) = 11.75 + (29.7 / RPI) * atan(RPI * 0.031 * (t1 - 15.4))

    # Unpack shared CN parameters
    cwd_flig = cn_params.cwd_flig
    rf_cwdl2 = cn_params.rf_cwdl2
    minpsi   = cn_params.minpsi
    maxpsi   = cn_params.maxpsi
    Q10      = cn_params.Q10
    froz_q10 = cn_params.froz_q10
    decomp_depth_efolding = cn_params.decomp_depth_efolding
    mino2lim = cn_params.mino2lim

    # Validate config
    if bgc_state.use_century_tfunc && bgc_state.normalize_q10_to_century_tfunc
        error("Cannot have both use_century_tfunc and normalize_q10_to_century_tfunc set as true")
    end

    # Translate turnover times to per-second rate constants
    k_l1    = 1.0 / (SECSPDAY * days_per_year * params.tau_l1_bgc)
    k_l2_l3 = 1.0 / (SECSPDAY * days_per_year * params.tau_l2_l3_bgc)
    k_s1    = 1.0 / (SECSPDAY * days_per_year * params.tau_s1_bgc)
    k_s2    = 1.0 / (SECSPDAY * days_per_year * params.tau_s2_bgc)
    k_s3    = 1.0 / (SECSPDAY * days_per_year * params.tau_s3_bgc)
    k_frag  = 1.0 / (SECSPDAY * days_per_year * cn_params.tau_cwd)

    # Reference rate at 30°C
    catanf_30 = catanf(30.0)

    # Unpack pool/transition indices
    i_met_lit = 1  # always 1
    i_cel_lit = bgc_state.i_cel_lit
    i_lig_lit = bgc_state.i_lig_lit
    i_act_som = bgc_state.i_act_som
    i_slo_som = bgc_state.i_slo_som
    i_pas_som = bgc_state.i_pas_som
    i_l1s1    = bgc_state.i_l1s1
    i_l2s1    = bgc_state.i_l2s1
    i_l3s2    = bgc_state.i_l3s2
    i_s1s2    = bgc_state.i_s1s2
    i_s1s3    = bgc_state.i_s1s3
    i_s2s1    = bgc_state.i_s2s1
    i_s2s3    = bgc_state.i_s2s3
    i_s3s1    = bgc_state.i_s3s1
    i_cwdl3   = bgc_state.i_cwdl3

    spinup_factor = cascade_con.spinup_factor

    # --- Compute spinup geographic terms ---
    spinup_geogterm_l1  = ones(nc)
    spinup_geogterm_l23 = ones(nc)
    spinup_geogterm_cwd = ones(nc)
    spinup_geogterm_s1  = ones(nc)
    spinup_geogterm_s2  = ones(nc)
    spinup_geogterm_s3  = ones(nc)

    if spinup_state >= 1
        for c in 1:nc
            if !mask_bgc_soilc[c]
                continue
            end
            lat = latdeg[col_gridcell[c]]

            if abs(spinup_factor[i_met_lit] - 1.0) > eps_val
                spinup_geogterm_l1[c] = spinup_factor[i_met_lit] * get_spinup_latitude_term(lat)
            end

            if abs(spinup_factor[i_cel_lit] - 1.0) > eps_val
                spinup_geogterm_l23[c] = spinup_factor[i_cel_lit] * get_spinup_latitude_term(lat)
            end

            if !use_fates
                i_cwd_local = i_pas_som + 1
                if abs(spinup_factor[i_cwd_local] - 1.0) > eps_val
                    spinup_geogterm_cwd[c] = spinup_factor[i_cwd_local] * get_spinup_latitude_term(lat)
                end
            end

            if abs(spinup_factor[i_act_som] - 1.0) > eps_val
                spinup_geogterm_s1[c] = spinup_factor[i_act_som] * get_spinup_latitude_term(lat)
            end

            if abs(spinup_factor[i_slo_som] - 1.0) > eps_val
                spinup_geogterm_s2[c] = spinup_factor[i_slo_som] * get_spinup_latitude_term(lat)
            end

            if abs(spinup_factor[i_pas_som] - 1.0) > eps_val
                spinup_geogterm_s3[c] = spinup_factor[i_pas_som] * get_spinup_latitude_term(lat)
            end
        end
    end

    # --- Time-dependent coefficients ---
    t_scalar = cf.t_scalar_col
    w_scalar = cf.w_scalar_col
    o_scalar = cf.o_scalar_col
    decomp_k = cf.decomp_k_col

    depth_scalar = zeros(nc, nlevdecomp)

    if nlevdecomp == 1
        # ---- Single-level decomposition ----
        # Weight temperature and water potential scalars by rooting fraction
        nlev_soildecomp_standard = 5
        frw = zeros(nc)
        fr  = zeros(nc, nlev_soildecomp_standard)

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

        if !bgc_state.use_century_tfunc
            # Q10 temperature scalar
            for j in 1:nlev_soildecomp_standard
                for c in 1:nc
                    mask_bgc_soilc[c] || continue
                    if j == 1
                        t_scalar[c, :] .= 0.0
                    end
                    if t_soisno[c, j] >= TFRZ
                        t_scalar[c, 1] += (Q10^((t_soisno[c, j] - (TFRZ + 25.0)) / 10.0)) * fr[c, j]
                    else
                        t_scalar[c, 1] += (Q10^(-25.0 / 10.0)) *
                            (froz_q10^((t_soisno[c, j] - TFRZ) / 10.0)) * fr[c, j]
                    end
                end
            end
        else
            # CENTURY arctangent temperature function
            for j in 1:nlev_soildecomp_standard
                for c in 1:nc
                    mask_bgc_soilc[c] || continue
                    if j == 1
                        t_scalar[c, :] .= 0.0
                    end
                    t_scalar[c, 1] += max(catanf(t_soisno[c, j] - TFRZ) / catanf_30 * fr[c, j], 0.01)
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
                psi = min(soilpsi[c, j], maxpsi)
                if psi > minpsi
                    w_scalar[c, 1] += (log(minpsi / psi) / log(minpsi / maxpsi)) * fr[c, j]
                end
            end
        end

        # O2 scalar
        if use_lch4 && anoxia
            for j in 1:nlev_soildecomp_standard
                for c in 1:nc
                    mask_bgc_soilc[c] || continue
                    if j == 1
                        o_scalar[c, :] .= 0.0
                    end
                    o_scalar[c, 1] += fr[c, j] * max(o2stress_unsat[c, j], mino2lim)
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
        # ---- Multi-level decomposition ----

        if !bgc_state.use_century_tfunc
            # Q10 temperature scalar
            for j in 1:nlevdecomp
                for c in 1:nc
                    mask_bgc_soilc[c] || continue
                    if t_soisno[c, j] >= TFRZ
                        t_scalar[c, j] = Q10^((t_soisno[c, j] - (TFRZ + 25.0)) / 10.0)
                    else
                        t_scalar[c, j] = (Q10^(-25.0 / 10.0)) *
                            (froz_q10^((t_soisno[c, j] - TFRZ) / 10.0))
                    end
                end
            end
        else
            # CENTURY arctangent temperature function
            for j in 1:nlevdecomp
                for c in 1:nc
                    mask_bgc_soilc[c] || continue
                    t_scalar[c, j] = max(catanf(t_soisno[c, j] - TFRZ) / catanf_30, 0.01)
                end
            end
        end

        # Water scalar
        for j in 1:nlevdecomp
            for c in 1:nc
                mask_bgc_soilc[c] || continue
                psi = min(soilpsi[c, j], maxpsi)
                if psi > minpsi
                    w_scalar[c, j] = log(minpsi / psi) / log(minpsi / maxpsi)
                else
                    w_scalar[c, j] = 0.0
                end
            end
        end

        # O2 scalar
        if use_lch4 && anoxia
            for j in 1:nlevdecomp
                for c in 1:nc
                    mask_bgc_soilc[c] || continue
                    o_scalar[c, j] = max(o2stress_unsat[c, j], mino2lim)
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

    # Normalize Q10 to CENTURY temperature function if requested
    if bgc_state.normalize_q10_to_century_tfunc
        normalization_factor = (catanf(bgc_state.normalization_tref) / catanf_30) /
            (Q10^((bgc_state.normalization_tref - 25.0) / 10.0))
        for j in 1:nlevdecomp
            for c in 1:nc
                mask_bgc_soilc[c] || continue
                t_scalar[c, j] *= normalization_factor
            end
        end
    end

    # Depth scalar — fixed e-folding depth
    for j in 1:nlevdecomp
        for c in 1:nc
            mask_bgc_soilc[c] || continue
            depth_scalar[c, j] = exp(-zsoi_vals[j] / decomp_depth_efolding)
        end
    end

    # --- Calculate rate constants for all pools ---
    for j in 1:nlevdecomp
        for c in 1:nc
            mask_bgc_soilc[c] || continue

            rate = t_scalar[c, j] * w_scalar[c, j] * depth_scalar[c, j] * o_scalar[c, j]

            decomp_k[c, j, i_met_lit] = k_l1    * rate * spinup_geogterm_l1[c]
            decomp_k[c, j, i_cel_lit] = k_l2_l3 * rate * spinup_geogterm_l23[c]
            decomp_k[c, j, i_lig_lit] = k_l2_l3 * rate * spinup_geogterm_l23[c]
            decomp_k[c, j, i_act_som] = k_s1    * rate * spinup_geogterm_s1[c]
            decomp_k[c, j, i_slo_som] = k_s2    * rate * spinup_geogterm_s2[c]
            decomp_k[c, j, i_pas_som] = k_s3    * rate * spinup_geogterm_s3[c]

            if !use_fates
                i_cwd_local = i_pas_som + 1
                decomp_k[c, j, i_cwd_local] = k_frag * rate * spinup_geogterm_cwd[c]
            end
        end
    end

    # --- Set pathfrac and rf for the cascade ---
    rf_decomp_cascade       = cf.rf_decomp_cascade_col
    pathfrac_decomp_cascade = cf.pathfrac_decomp_cascade_col

    for j in 1:nlevdecomp
        for c in 1:nc
            pathfrac_decomp_cascade[c, j, i_l1s1] = 1.0
            pathfrac_decomp_cascade[c, j, i_l2s1] = 1.0
            pathfrac_decomp_cascade[c, j, i_l3s2] = 1.0
            pathfrac_decomp_cascade[c, j, i_s1s2] = bgc_state.f_s1s2[c, j]
            pathfrac_decomp_cascade[c, j, i_s1s3] = bgc_state.f_s1s3[c, j]
            pathfrac_decomp_cascade[c, j, i_s2s1] = bgc_state.f_s2s1
            pathfrac_decomp_cascade[c, j, i_s2s3] = bgc_state.f_s2s3
            pathfrac_decomp_cascade[c, j, i_s3s1] = 1.0

            rf_decomp_cascade[c, j, i_l1s1] = bgc_state.rf_l1s1
            rf_decomp_cascade[c, j, i_l2s1] = bgc_state.rf_l2s1
            rf_decomp_cascade[c, j, i_l3s2] = bgc_state.rf_l3s2
            rf_decomp_cascade[c, j, i_s1s2] = bgc_state.rf_s1s2[c, j]
            rf_decomp_cascade[c, j, i_s1s3] = bgc_state.rf_s1s3[c, j]
            rf_decomp_cascade[c, j, i_s2s1] = bgc_state.rf_s2s1
            rf_decomp_cascade[c, j, i_s2s3] = bgc_state.rf_s2s3
            rf_decomp_cascade[c, j, i_s3s1] = bgc_state.rf_s3s1
        end
    end

    if !use_fates
        i_cwdl2_local = i_cwdl3 - 1  # i_cwdl2 = i_cwdl3 - 1 = 9
        for j in 1:nlevdecomp
            for c in 1:nc
                pathfrac_decomp_cascade[c, j, i_cwdl2_local] = bgc_state.cwd_fcel
                pathfrac_decomp_cascade[c, j, i_cwdl3]       = cwd_flig
                rf_decomp_cascade[c, j, i_cwdl2_local]       = rf_cwdl2
                rf_decomp_cascade[c, j, i_cwdl3]             = bgc_state.rf_cwdl3
            end
        end
    end

    return nothing
end
