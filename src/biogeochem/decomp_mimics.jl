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

# --------------------------------------------------------------------------
# GPU kernels for per-(column, level) MIMICS loops with fully independent
# iterations. @kernel/@index/@Const and _launch! come from
# src/infrastructure/kernels.jl. Each kernel is backend-agnostic.
# --------------------------------------------------------------------------

# Soil texture-dependent MIMICS coefficients (per column AND decomp level).
# Each (c,j) writes desorp/fphys_m1/fphys_m2/p_scalar from cellclay[c,j] only —
# fully independent. desorp is the launch `out` (sets backend/ndrange).
@kernel function _mimics_texture_kernel!(desorp, fphys_m1, fphys_m2, p_scalar,
                                         @Const(cellclay),
                                         desorp_p1, desorp_p2,
                                         fphys_r_p1, fphys_r_p2,
                                         fphys_k_p1, fphys_k_p2,
                                         p_scalar_p1, p_scalar_p2, pct_to_frac)
    c, j = @index(Global, NTuple)
    @inbounds begin
        T = eltype(desorp)
        clay_frac = pct_to_frac * smooth_min(T(100.0), cellclay[c, j])
        desorp[c, j]   = desorp_p1 * exp(desorp_p2 * clay_frac)
        fphys_m1[c, j] = smooth_min(one(T), fphys_r_p1 * exp(fphys_r_p2 * clay_frac))
        fphys_m2[c, j] = smooth_min(one(T), fphys_k_p1 * exp(fphys_k_p2 * clay_frac))
        p_scalar[c, j] = one(T) / (p_scalar_p1 * exp(p_scalar_p2 * sqrt(clay_frac)))
    end
end

"""
    mimics_texture_params!(desorp, fphys_m1, fphys_m2, p_scalar, cellclay,
        desorp_p1, desorp_p2, fphys_r_p1, fphys_r_p2, fphys_k_p1, fphys_k_p2,
        p_scalar_p1, p_scalar_p2, pct_to_frac; nc, nlevdecomp)

Soil texture-dependent MIMICS coefficients over all (column, decomp level)
pairs — a 2D backend-agnostic kernel. One thread per (c,j).
"""
function mimics_texture_params!(desorp, fphys_m1, fphys_m2, p_scalar, cellclay,
                                desorp_p1, desorp_p2, fphys_r_p1, fphys_r_p2,
                                fphys_k_p1, fphys_k_p2, p_scalar_p1, p_scalar_p2,
                                pct_to_frac; nc::Int, nlevdecomp::Int)
    _T = eltype(desorp)   # convert the texture scalar params to working precision
    _launch!(_mimics_texture_kernel!, desorp, fphys_m1, fphys_m2, p_scalar, cellclay,
             _T(desorp_p1), _T(desorp_p2), _T(fphys_r_p1), _T(fphys_r_p2),
             _T(fphys_k_p1), _T(fphys_k_p2), _T(p_scalar_p1), _T(p_scalar_p2),
             _T(pct_to_frac); ndrange = (nc, nlevdecomp))
end

# Multi-level MIMICS water scalar (per column AND decomp level), masked.
# w_scalar[c,j] depends only on soilpsi[c,j] — fully independent.
@kernel function _mimics_wscalar_kernel!(w_scalar, @Const(mask), @Const(soilpsi),
                                         minpsi, maxpsi)
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        psi = smooth_min(soilpsi[c, j], maxpsi)
        if psi > minpsi
            w_scalar[c, j] = log(minpsi / psi) / log(minpsi / maxpsi)
        else
            w_scalar[c, j] = zero(eltype(w_scalar))
        end
    end
end

"""
    mimics_water_scalar!(w_scalar, mask, soilpsi, minpsi, maxpsi; nc, nlevdecomp)

Multi-level MIMICS water scalar over all (column, decomp level) pairs (masked).
Backend-agnostic; one thread per (c,j).
"""
function mimics_water_scalar!(w_scalar, mask, soilpsi, minpsi, maxpsi;
                              nc::Int, nlevdecomp::Int)
    _T = eltype(w_scalar)
    _launch!(_mimics_wscalar_kernel!, w_scalar, mask, soilpsi, _T(minpsi), _T(maxpsi);
             ndrange = (nc, nlevdecomp))
end

# Multi-level MIMICS O2 scalar under anoxia (per column AND level), masked.
# o_scalar[c,j] depends only on o2stress_unsat[c,j] — fully independent.
@kernel function _mimics_oscalar_anoxia_kernel!(o_scalar, @Const(mask),
                                                @Const(o2stress_unsat), mino2lim)
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        o_scalar[c, j] = smooth_max(o2stress_unsat[c, j], mino2lim)
    end
end

"""
    mimics_oscalar_anoxia!(o_scalar, mask, o2stress_unsat, mino2lim; nc, nlevdecomp)

Multi-level MIMICS O2 scalar under anoxia over all (column, decomp level) pairs
(masked). Backend-agnostic; one thread per (c,j).
"""
function mimics_oscalar_anoxia!(o_scalar, mask, o2stress_unsat, mino2lim;
                                nc::Int, nlevdecomp::Int)
    _launch!(_mimics_oscalar_anoxia_kernel!, o_scalar, mask, o2stress_unsat,
             eltype(o_scalar)(mino2lim); ndrange = (nc, nlevdecomp))
end

# Multi-level MIMICS depth scalar (per column AND level), masked: set to 1.0.
@kernel function _mimics_depth_scalar_kernel!(depth_scalar, @Const(mask))
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        depth_scalar[c, j] = one(eltype(depth_scalar))
    end
end

"""
    mimics_depth_scalar!(depth_scalar, mask; nc, nlevdecomp)

Set the masked columns' MIMICS depth scalar to 1.0 over all (column, level)
pairs (placeholder, as in Fortran). Backend-agnostic; one thread per (c,j).
"""
function mimics_depth_scalar!(depth_scalar, mask; nc::Int, nlevdecomp::Int)
    _launch!(_mimics_depth_scalar_kernel!, depth_scalar, mask;
             ndrange = (nc, nlevdecomp))
end

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
Base.@kwdef mutable struct DecompMIMICSState{M<:AbstractMatrix}
    # Spatially-varying arrays (col × nlevdecomp). Array-type param M stays LOOSE so
    # adapt(MetalF32/MtlArray, ·) can reconstruct the struct with device Float32 arrays
    # (pinning ::Matrix{Float64} forces a host-Float64 convert on assignment). The scalar
    # coeff fields stay Float64 — they never reach a kernel (bundled into _DMRCoef{FT} at FT).
    desorp    ::M = Matrix{Float64}(undef, 0, 0)
    fphys_m1  ::M = Matrix{Float64}(undef, 0, 0)
    fphys_m2  ::M = Matrix{Float64}(undef, 0, 0)
    p_scalar  ::M = Matrix{Float64}(undef, 0, 0)

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

# Device-movable: lets the static texture arrays (desorp/fphys_m1/fphys_m2/p_scalar)
# ride to the GPU for decomp_rates_mimics! (scalar/Int fields pass through unchanged).
Adapt.@adapt_structure DecompMIMICSState

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

    # Compute soil texture-dependent parameters (2D kernel; each (c,j) independent)
    mimics_texture_params!(mimics_state.desorp, mimics_state.fphys_m1,
                           mimics_state.fphys_m2, mimics_state.p_scalar, cellclay,
                           mimics_desorp_p1, mimics_desorp_p2,
                           mimics_fphys_r_p1, mimics_fphys_r_p2,
                           mimics_fphys_k_p1, mimics_fphys_k_p2,
                           mimics_p_scalar_p1, mimics_p_scalar_p2, PCT_TO_FRAC;
                           nc = nc, nlevdecomp = nlevdecomp)
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
# Device-view bundles + kernels for decomp_rates_mimics!
#
# decomp_rates_mimics! touches ~16 column/level arrays and ~60 scalar
# coefficients/indices. To stay under the Metal ~31 kernel-arg limit the arrays
# are bundled into one Adapt.@adapt_structure'd struct, the (Int) pool/transition
# indices into one isbits struct, and the (working-precision) scalar coefficients
# into a second isbits struct. Field names mirror the local variable paths so the
# kernel bodies read verbatim. Every per-column compute loop becomes a per-column
# kernel with internal (sequential, own-[c,j]) j/k loops — race-free and (on the
# Float64 CPU path) byte-identical to the original host loops.
# ---------------------------------------------------------------------------

# Column/level array bundle (M = 2D col×lev, A3 = 3D col×lev×pool/transition,
# V = per-column vectors). All writes flow back into the live SoA structs.
Base.@kwdef struct _DMRArr{V, M, A3}
    # outputs (written)
    cn_col::M
    w_scalar::M
    o_scalar::M
    decomp_k::A3
    pathfrac_decomp_cascade::A3
    rf_decomp_cascade::A3
    depth_scalar::M
    # inputs (read)
    t_soisno::M
    decomp_cpools_vr::A3
    col_dz::M
    fphys_m1::M
    fphys_m2::M
    desorp_arr::M
    p_scalar_arr::M
    ligninNratioAvg::V
    annsum_npp_col::V
end
Adapt.@adapt_structure _DMRArr

_dmr_arr(cf, mimics_state, depth_scalar, t_soisno, decomp_cpools_vr, col_dz,
         ligninNratioAvg, annsum_npp_col) = _DMRArr(;
    cn_col                  = cf.cn_col,
    w_scalar                = cf.w_scalar_col,
    o_scalar                = cf.o_scalar_col,
    decomp_k                = cf.decomp_k_col,
    pathfrac_decomp_cascade = cf.pathfrac_decomp_cascade_col,
    rf_decomp_cascade       = cf.rf_decomp_cascade_col,
    depth_scalar            = depth_scalar,
    t_soisno                = t_soisno,
    decomp_cpools_vr        = decomp_cpools_vr,
    col_dz                  = col_dz,
    fphys_m1                = mimics_state.fphys_m1,
    fphys_m2                = mimics_state.fphys_m2,
    desorp_arr              = mimics_state.desorp,
    p_scalar_arr            = mimics_state.p_scalar,
    ligninNratioAvg         = ligninNratioAvg,
    annsum_npp_col          = annsum_npp_col)

# Pool/transition indices — isbits (all Int); passed by value to kernels.
Base.@kwdef struct _DMRIdx
    i_met_lit::Int; i_str_lit::Int; i_avl_som::Int; i_chem_som::Int
    i_phys_som::Int; i_cop_mic::Int; i_oli_mic::Int
    i_l1m1::Int; i_l1m2::Int; i_l2m1::Int; i_l2m2::Int
    i_s1m1::Int; i_s1m2::Int; i_m1s1::Int; i_m1s2::Int
    i_m2s1::Int; i_m2s2::Int; i_s2s1::Int; i_s3s1::Int
    i_m1s3::Int; i_m2s3::Int
end

# Scalar MIMICS regression coefficients — isbits, parametric on working precision T.
# A Float64 scalar kernel arg is invalid Metal IR, so every coefficient is stored at
# the state eltype; on Float64 T(x) === x (byte-identical).
Base.@kwdef struct _DMRCoef{T}
    eps_val::T
    mimics_fmet_p1::T; mimics_fmet_p2::T; mimics_fmet_p3::T; mimics_fmet_p4::T
    mimics_fchem_r_p1::T; mimics_fchem_r_p2::T; mimics_fchem_k_p1::T; mimics_fchem_k_p2::T
    mimics_tau_mod_min::T; mimics_tau_mod_max::T; mimics_tau_mod_factor::T
    mimics_tau_r_p1::T; mimics_tau_r_p2::T; mimics_tau_k_p1::T; mimics_tau_k_p2::T
    mimics_ko_r::T; mimics_ko_k::T; mimics_densdep::T; mimics_desorpQ10::T
    mimics_t_soi_ref::T; mimics_cn_mod_num::T; mimics_cn_r::T; mimics_cn_k::T
    vslope_l1_m1::T; vslope_l2_m1::T; vslope_s1_m1::T; vslope_l1_m2::T; vslope_l2_m2::T; vslope_s1_m2::T
    kslope_l1_m1::T; kslope_l2_m1::T; kslope_s1_m1::T; kslope_l1_m2::T; kslope_l2_m2::T; kslope_s1_m2::T
    vint_l1_m1::T; vint_l2_m1::T; vint_s1_m1::T; vint_l1_m2::T; vint_l2_m2::T; vint_s1_m2::T
    kint_l1_m1::T; kint_l2_m1::T; kint_s1_m1::T; kint_l1_m2::T; kint_l2_m2::T; kint_s1_m2::T
    vmod_l1_m1::T; vmod_l2_m1::T; vmod_s1_m1::T; vmod_l1_m2::T; vmod_l2_m2::T; vmod_s1_m2::T
    kmod_l1_m1::T; kmod_l2_m1::T; kmod_s1_m1::T; kmod_l1_m2::T; kmod_l2_m2::T; kmod_s1_m2::T
    k_frag::T; secsphr::T; tfrz::T; g_to_mg::T; cm3_to_m3::T
end

# --- nlevdecomp==1 helper kernels (per-column reductions over standard soil levels) ---

# frw[c] = sum_j col_dz[c,j]  (own-column reduction over j=1..nls; ascending order)
@kernel function _dmr_frw_kernel!(frw, @Const(mask), @Const(col_dz), nls::Int)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(frw)
        s = zero(T)
        for j in 1:nls
            s += col_dz[c, j]
        end
        frw[c] = s
    end
end

# w_scalar[c,1] = sum_j (log(minpsi/psi)/log(minpsi/maxpsi)) * fr[c,j]  (own-column).
# fr[c,j] = col_dz[c,j]/frw[c] (or 0). Mirrors the host loop exactly: w_scalar[c,:]
# is zeroed (here over the single level), then accumulated in ascending j order.
@kernel function _dmr_wscalar1_kernel!(w_scalar, @Const(mask), @Const(soilpsi),
                                       @Const(col_dz), @Const(frw),
                                       nls::Int, minpsi, maxpsi)
    c = @index(Global)
    @inbounds if mask[c]
        T = typeof(minpsi)
        # Zero the FULL allocated row (host did `w_scalar[c, :] .= 0.0`); the array
        # is sized (nc, nlevdecomp_full=5) even when nls==1, and leaving levels 2..5
        # stale (NaN) would diverge from the host.
        for jj in 1:size(w_scalar, 2)
            w_scalar[c, jj] = zero(T)
        end
        for j in 1:nls
            fr_cj = frw[c] != zero(T) ? col_dz[c, j] / frw[c] : zero(T)
            psi = smooth_min(soilpsi[c, j], maxpsi)
            if psi > minpsi
                w_scalar[c, 1] += (log(minpsi / psi) / log(minpsi / maxpsi)) * fr_cj
            end
        end
    end
end

# o_scalar[c,1] under anoxia = sum_j fr[c,j]*smooth_max(o2stress_unsat[c,j], mino2lim)
@kernel function _dmr_oscalar1_anoxia_kernel!(o_scalar, @Const(mask),
                                              @Const(o2stress_unsat), @Const(col_dz),
                                              @Const(frw), nls::Int, mino2lim)
    c = @index(Global)
    @inbounds if mask[c]
        T = typeof(mino2lim)
        # Zero the FULL allocated row (host did `o_scalar[c, :] .= 0.0`).
        for jj in 1:size(o_scalar, 2)
            o_scalar[c, jj] = zero(T)
        end
        for j in 1:nls
            fr_cj = frw[c] != zero(T) ? col_dz[c, j] / frw[c] : zero(T)
            o_scalar[c, 1] += fr_cj * smooth_max(o2stress_unsat[c, j], mino2lim)
        end
    end
end

# o_scalar[c,j] = 1 over all levels (no anoxia path); unmasked, as in the host loop.
@kernel function _dmr_oscalar_one_kernel!(o_scalar, nlevdecomp::Int)
    c = @index(Global)
    @inbounds begin
        T = eltype(o_scalar)
        for j in 1:nlevdecomp
            o_scalar[c, j] = one(T)
        end
    end
end

# --- Spinup geographic terms (per-column; only run when spinup_state >= 1) ---
@kernel function _dmr_spinup_geogterm_kernel!(geo_l1, geo_l2, geo_cwd, geo_s1, geo_s2,
        geo_s3, geo_m1, geo_m2, @Const(mask), @Const(latdeg), @Const(col_gridcell),
        @Const(spinup_factor), idx::_DMRIdx, use_fates::Bool, eps_val)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(geo_l1)
        one_T = one(T)
        lat = latdeg[col_gridcell[c]]
        slt = T(1.0) + T(50.0) / (T(1.0) + exp(T(-0.15) * (abs(lat) - T(60.0))))
        if abs(spinup_factor[idx.i_met_lit] - one_T) > eps_val
            geo_l1[c] = spinup_factor[idx.i_met_lit] * slt
        end
        if abs(spinup_factor[idx.i_str_lit] - one_T) > eps_val
            geo_l2[c] = spinup_factor[idx.i_str_lit] * slt
        end
        if !use_fates
            i_cwd_local = idx.i_oli_mic + 1
            if abs(spinup_factor[i_cwd_local] - one_T) > eps_val
                geo_cwd[c] = spinup_factor[i_cwd_local] * slt
            end
        end
        if abs(spinup_factor[idx.i_avl_som] - one_T) > eps_val
            geo_s1[c] = spinup_factor[idx.i_avl_som] * slt
        end
        if abs(spinup_factor[idx.i_chem_som] - one_T) > eps_val
            geo_s2[c] = spinup_factor[idx.i_chem_som] * slt
        end
        if abs(spinup_factor[idx.i_phys_som] - one_T) > eps_val
            geo_s3[c] = spinup_factor[idx.i_phys_som] * slt
        end
        if abs(spinup_factor[idx.i_cop_mic] - one_T) > eps_val
            geo_m1[c] = spinup_factor[idx.i_cop_mic] * slt
        end
        if abs(spinup_factor[idx.i_oli_mic] - one_T) > eps_val
            geo_m2[c] = spinup_factor[idx.i_oli_mic] * slt
        end
    end
end

# --- Main per-column rate kernel: internal j loop, own-[c,j] writes ---
# `_out` (a per-column array) only sets the launch ndrange/backend; all data lives
# in the `arr` bundle (same shared refs), mirroring the soil_fluxes! kernel pattern.
@kernel function _dmr_rates_kernel!(@Const(_out), arr::_DMRArr, @Const(mask), idx::_DMRIdx,
                                    coef::_DMRCoef, nlevdecomp::Int, use_fates::Bool)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(arr.cn_col)
        zero_T = zero(T); one_T = one(T)

        annsum_npp_col_scalar = smooth_max(zero_T, arr.annsum_npp_col[c])

        # fmet: fraction metabolic litter
        fmet = coef.mimics_fmet_p1 * (coef.mimics_fmet_p2 - coef.mimics_fmet_p3 *
            smooth_min(coef.mimics_fmet_p4, arr.ligninNratioAvg[c]))
        tau_mod = smooth_min(coef.mimics_tau_mod_max, smooth_max(coef.mimics_tau_mod_min,
            sqrt(coef.mimics_tau_mod_factor * annsum_npp_col_scalar)))

        tau_m1 = coef.mimics_tau_r_p1 * exp(coef.mimics_tau_r_p2 * fmet) * tau_mod / coef.secsphr
        tau_m2 = coef.mimics_tau_k_p1 * exp(coef.mimics_tau_k_p2 * fmet) * tau_mod / coef.secsphr

        # Microbial C:N ratios
        arr.cn_col[c, idx.i_cop_mic] = coef.mimics_cn_r * sqrt(coef.mimics_cn_mod_num / fmet)
        arr.cn_col[c, idx.i_oli_mic] = coef.mimics_cn_k * sqrt(coef.mimics_cn_mod_num / fmet)

        # Chemical fractionation
        fchem_m1 = smooth_min(one_T, smooth_max(zero_T, coef.mimics_fchem_r_p1 *
                    exp(coef.mimics_fchem_r_p2 * fmet)))
        fchem_m2 = smooth_min(one_T, smooth_max(zero_T, coef.mimics_fchem_k_p1 *
                    exp(coef.mimics_fchem_k_p2 * fmet)))

        for j in 1:nlevdecomp
            t_soi_degC = arr.t_soisno[c, j] - coef.tfrz

            vmax_l1_m1 = exp(coef.vslope_l1_m1 * t_soi_degC + coef.vint_l1_m1) * coef.vmod_l1_m1 / coef.secsphr
            vmax_l1_m2 = exp(coef.vslope_l1_m2 * t_soi_degC + coef.vint_l1_m2) * coef.vmod_l1_m2 / coef.secsphr
            vmax_l2_m1 = exp(coef.vslope_l2_m1 * t_soi_degC + coef.vint_l2_m1) * coef.vmod_l2_m1 / coef.secsphr
            vmax_l2_m2 = exp(coef.vslope_l2_m2 * t_soi_degC + coef.vint_l2_m2) * coef.vmod_l2_m2 / coef.secsphr
            vmax_s1_m1 = exp(coef.vslope_s1_m1 * t_soi_degC + coef.vint_s1_m1) * coef.vmod_s1_m1 / coef.secsphr
            vmax_s1_m2 = exp(coef.vslope_s1_m2 * t_soi_degC + coef.vint_s1_m2) * coef.vmod_s1_m2 / coef.secsphr

            km_l1_m1 = exp(coef.kslope_l1_m1 * t_soi_degC + coef.kint_l1_m1) * coef.kmod_l1_m1
            km_l1_m2 = exp(coef.kslope_l1_m2 * t_soi_degC + coef.kint_l1_m2) * coef.kmod_l1_m2
            km_l2_m1 = exp(coef.kslope_l2_m1 * t_soi_degC + coef.kint_l2_m1) * coef.kmod_l2_m1
            km_l2_m2 = exp(coef.kslope_l2_m2 * t_soi_degC + coef.kint_l2_m2) * coef.kmod_l2_m2
            km_s1_m1 = exp(coef.kslope_s1_m1 * t_soi_degC + coef.kint_s1_m1) * coef.kmod_s1_m1 * arr.p_scalar_arr[c, j]
            km_s1_m2 = exp(coef.kslope_s1_m2 * t_soi_degC + coef.kint_s1_m2) * coef.kmod_s1_m2 * arr.p_scalar_arr[c, j]

            desorption = (arr.desorp_arr[c, j] / coef.secsphr) * coef.mimics_desorpQ10 *
                         exp((t_soi_degC - coef.mimics_t_soi_ref) / T(10.0))

            m1_conc = (arr.decomp_cpools_vr[c, j, idx.i_cop_mic] / arr.col_dz[c, j]) * coef.g_to_mg * coef.cm3_to_m3
            m2_conc = (arr.decomp_cpools_vr[c, j, idx.i_oli_mic] / arr.col_dz[c, j]) * coef.g_to_mg * coef.cm3_to_m3

            w_d_o_scalars = arr.w_scalar[c, j] * arr.depth_scalar[c, j] * arr.o_scalar[c, j]

            # Metabolic litter decomposition
            term_1 = vmax_l1_m1 * m1_conc / (km_l1_m1 + m1_conc)
            term_2 = vmax_l1_m2 * m2_conc / (km_l1_m2 + m2_conc)
            arr.decomp_k[c, j, idx.i_met_lit] = (term_1 + term_2) * w_d_o_scalars
            if (term_1 + term_2) != zero_T
                arr.pathfrac_decomp_cascade[c, j, idx.i_l1m1] = term_1 / (term_1 + term_2)
                arr.pathfrac_decomp_cascade[c, j, idx.i_l1m2] = term_2 / (term_1 + term_2)
            else
                arr.pathfrac_decomp_cascade[c, j, idx.i_l1m1] = zero_T
                arr.pathfrac_decomp_cascade[c, j, idx.i_l1m2] = zero_T
            end

            # Structural litter decomposition
            term_1 = vmax_l2_m1 * m1_conc / (km_l2_m1 + m1_conc)
            term_2 = vmax_l2_m2 * m2_conc / (km_l2_m2 + m2_conc)
            arr.decomp_k[c, j, idx.i_str_lit] = (term_1 + term_2) * w_d_o_scalars
            if (term_1 + term_2) != zero_T
                arr.pathfrac_decomp_cascade[c, j, idx.i_l2m1] = term_1 / (term_1 + term_2)
                arr.pathfrac_decomp_cascade[c, j, idx.i_l2m2] = term_2 / (term_1 + term_2)
            else
                arr.pathfrac_decomp_cascade[c, j, idx.i_l2m1] = zero_T
                arr.pathfrac_decomp_cascade[c, j, idx.i_l2m2] = zero_T
            end

            # Available SOM decomposition
            term_1 = vmax_s1_m1 * m1_conc / (km_s1_m1 + m1_conc)
            term_2 = vmax_s1_m2 * m2_conc / (km_s1_m2 + m2_conc)
            arr.decomp_k[c, j, idx.i_avl_som] = (term_1 + term_2) * w_d_o_scalars
            if (term_1 + term_2) != zero_T
                arr.pathfrac_decomp_cascade[c, j, idx.i_s1m1] = term_1 / (term_1 + term_2)
                arr.pathfrac_decomp_cascade[c, j, idx.i_s1m2] = term_2 / (term_1 + term_2)
            else
                arr.pathfrac_decomp_cascade[c, j, idx.i_s1m1] = zero_T
                arr.pathfrac_decomp_cascade[c, j, idx.i_s1m2] = zero_T
            end

            # Physically protected SOM: desorption only
            arr.decomp_k[c, j, idx.i_phys_som] = desorption * arr.depth_scalar[c, j]

            # Chemically protected SOM: oxidation
            term_1 = vmax_l2_m1 * m1_conc / (coef.mimics_ko_r * km_l2_m1 + m1_conc)
            term_2 = vmax_l2_m2 * m2_conc / (coef.mimics_ko_k * km_l2_m2 + m2_conc)
            arr.decomp_k[c, j, idx.i_chem_som] = (term_1 + term_2) * w_d_o_scalars

            # Copiotrophic microbe turnover (M1)
            arr.decomp_k[c, j, idx.i_cop_mic] = tau_m1 *
                   m1_conc^(coef.mimics_densdep - one_T) * w_d_o_scalars
            favl = smooth_min(one_T, smooth_max(zero_T, one_T - arr.fphys_m1[c, j] - fchem_m1))
            arr.pathfrac_decomp_cascade[c, j, idx.i_m1s1] = favl
            arr.pathfrac_decomp_cascade[c, j, idx.i_m1s2] = fchem_m1

            # Oligotrophic microbe turnover (M2)
            arr.decomp_k[c, j, idx.i_oli_mic] = tau_m2 *
                   m2_conc^(coef.mimics_densdep - one_T) * w_d_o_scalars
            favl = smooth_min(one_T, smooth_max(zero_T, one_T - arr.fphys_m2[c, j] - fchem_m2))
            arr.pathfrac_decomp_cascade[c, j, idx.i_m2s1] = favl
            arr.pathfrac_decomp_cascade[c, j, idx.i_m2s2] = fchem_m2

            # CWD fragmentation (only if FATES not enabled)
            if !use_fates
                i_cwd_local = idx.i_oli_mic + 1
                arr.decomp_k[c, j, i_cwd_local] = coef.k_frag * w_d_o_scalars
            end
        end
    end
end

# Post-assignment pathfrac terms (own-[c,j], all levels; unmasked as in host loop).
@kernel function _dmr_pathfrac_post_kernel!(@Const(_out), arr::_DMRArr, idx::_DMRIdx, nlevdecomp::Int)
    c = @index(Global)
    @inbounds begin
        T = eltype(arr.pathfrac_decomp_cascade)
        for j in 1:nlevdecomp
            arr.pathfrac_decomp_cascade[c, j, idx.i_s2s1] = one(T)
            arr.pathfrac_decomp_cascade[c, j, idx.i_s3s1] = one(T)
            arr.pathfrac_decomp_cascade[c, j, idx.i_m1s3] = arr.fphys_m1[c, j]
            arr.pathfrac_decomp_cascade[c, j, idx.i_m2s3] = arr.fphys_m2[c, j]
        end
    end
end

# CWD->L2 pathfrac + respiration fraction (only if FATES not enabled).
@kernel function _dmr_cwdl2_kernel!(@Const(_out), arr::_DMRArr, i_cwdl2_local::Int, nlevdecomp::Int, rf_cwdl2)
    c = @index(Global)
    @inbounds begin
        T = eltype(arr.pathfrac_decomp_cascade)
        for j in 1:nlevdecomp
            arr.pathfrac_decomp_cascade[c, j, i_cwdl2_local] = one(T)
            arr.rf_decomp_cascade[c, j, i_cwdl2_local] = rf_cwdl2
        end
    end
end

# Respiration fractions for all transitions (own-[c,j], all levels).
@kernel function _dmr_rf_kernel!(@Const(_out), arr::_DMRArr, idx::_DMRIdx, nlevdecomp::Int,
        rf_l1m1, rf_l1m2, rf_l2m1, rf_l2m2, rf_s1m1, rf_s1m2)
    c = @index(Global)
    @inbounds begin
        T = eltype(arr.rf_decomp_cascade)
        z = zero(T)
        for j in 1:nlevdecomp
            arr.rf_decomp_cascade[c, j, idx.i_l1m1] = rf_l1m1
            arr.rf_decomp_cascade[c, j, idx.i_l1m2] = rf_l1m2
            arr.rf_decomp_cascade[c, j, idx.i_l2m1] = rf_l2m1
            arr.rf_decomp_cascade[c, j, idx.i_l2m2] = rf_l2m2
            arr.rf_decomp_cascade[c, j, idx.i_s1m1] = rf_s1m1
            arr.rf_decomp_cascade[c, j, idx.i_s1m2] = rf_s1m2
            arr.rf_decomp_cascade[c, j, idx.i_s2s1] = z
            arr.rf_decomp_cascade[c, j, idx.i_s3s1] = z
            arr.rf_decomp_cascade[c, j, idx.i_m1s1] = z
            arr.rf_decomp_cascade[c, j, idx.i_m1s2] = z
            arr.rf_decomp_cascade[c, j, idx.i_m1s3] = z
            arr.rf_decomp_cascade[c, j, idx.i_m2s1] = z
            arr.rf_decomp_cascade[c, j, idx.i_m2s2] = z
            arr.rf_decomp_cascade[c, j, idx.i_m2s3] = z
        end
    end
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
- `mask_bgc_soilc::AbstractVector{Bool}`: mask for BGC soil columns
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
                               mask_bgc_soilc::AbstractVector{Bool},
                               bounds::UnitRange{Int},
                               nlevdecomp::Int,
                               t_soisno::AbstractMatrix{<:Real},
                               soilpsi::AbstractMatrix{<:Real},
                               decomp_cpools_vr::AbstractArray{<:Real,3},
                               col_dz::AbstractMatrix{<:Real},
                               ligninNratioAvg::AbstractVector{<:Real},
                               annsum_npp_col::AbstractVector{<:Real},
                               days_per_year::Real,
                               dt::Real,
                               spinup_state::Int=0,
                               use_lch4::Bool=false,
                               anoxia::Bool=false,
                               use_fates::Bool=false,
                               o2stress_unsat::AbstractMatrix{<:Real}=Matrix{Float64}(undef, 0, 0),
                               col_gridcell::AbstractVector{<:Integer}=Int[],
                               latdeg::AbstractVector{<:Real}=Float64[])

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

    # Pool/transition index bundle (isbits; passed by value to every kernel).
    idx = _DMRIdx(; i_met_lit, i_str_lit, i_avl_som, i_chem_som, i_phys_som,
                  i_cop_mic, i_oli_mic, i_l1m1, i_l1m2, i_l2m1, i_l2m2,
                  i_s1m1, i_s1m2, i_m1s1, i_m1s2, i_m2s1, i_m2s2,
                  i_s2s1, i_s3s1, i_m1s3, i_m2s3)

    # Working precision: take the eltype of a mutated output array so device
    # (Float32/Metal) and CPU (Float64) both pick the right scalar precision.
    FT = eltype(cf.decomp_k_col)

    # --- Compute spinup geographic terms ---
    # (Allocated to mirror the Fortran; not referenced downstream in this port.)
    spinup_geogterm_l1  = similar(cf.somc_fire_col, FT, nc); fill!(spinup_geogterm_l1, one(FT))
    spinup_geogterm_l2  = similar(spinup_geogterm_l1); fill!(spinup_geogterm_l2, one(FT))
    spinup_geogterm_cwd = similar(spinup_geogterm_l1); fill!(spinup_geogterm_cwd, one(FT))
    spinup_geogterm_s1  = similar(spinup_geogterm_l1); fill!(spinup_geogterm_s1, one(FT))
    spinup_geogterm_s2  = similar(spinup_geogterm_l1); fill!(spinup_geogterm_s2, one(FT))
    spinup_geogterm_s3  = similar(spinup_geogterm_l1); fill!(spinup_geogterm_s3, one(FT))
    spinup_geogterm_m1  = similar(spinup_geogterm_l1); fill!(spinup_geogterm_m1, one(FT))
    spinup_geogterm_m2  = similar(spinup_geogterm_l1); fill!(spinup_geogterm_m2, one(FT))

    if spinup_state >= 1
        # Move host-built index/factor vectors to the device, preserving Int eltype
        # for col_gridcell (a Float-coerced index would break device indexing).
        latdeg_d        = similar(cf.decomp_k_col, FT, length(latdeg)); copyto!(latdeg_d, latdeg)
        col_gridcell_d  = similar(cf.somc_fire_col, eltype(col_gridcell), length(col_gridcell))
        copyto!(col_gridcell_d, col_gridcell)
        spinup_factor_d = similar(cf.decomp_k_col, FT, length(spinup_factor))
        copyto!(spinup_factor_d, spinup_factor)
        _launch!(_dmr_spinup_geogterm_kernel!, spinup_geogterm_l1,
                 spinup_geogterm_l2, spinup_geogterm_cwd, spinup_geogterm_s1,
                 spinup_geogterm_s2, spinup_geogterm_s3, spinup_geogterm_m1,
                 spinup_geogterm_m2, mask_bgc_soilc, latdeg_d, col_gridcell_d,
                 spinup_factor_d, idx, use_fates, FT(eps_val))
    end

    # --- Time-dependent coefficients ---
    w_scalar                = cf.w_scalar_col
    o_scalar                = cf.o_scalar_col
    decomp_k                = cf.decomp_k_col
    pathfrac_decomp_cascade = cf.pathfrac_decomp_cascade_col
    rf_decomp_cascade       = cf.rf_decomp_cascade_col
    cn_col                  = cf.cn_col

    depth_scalar = similar(cf.decomp_k_col, FT, nc, nlevdecomp)
    fill!(depth_scalar, one(FT))

    # Working-precision scalar thresholds at the launch site (a Float64 scalar
    # kernel arg is invalid Metal IR).
    minpsi_ft   = FT(minpsi)
    maxpsi_ft   = FT(maxpsi)
    mino2lim_ft = FT(mino2lim)

    if nlevdecomp == 1
        # --- Single-level decomposition ---
        # frw[c] = sum_j col_dz[c,j] over standard soil levels (per-column reduction).
        frw = similar(cf.decomp_k_col, FT, nc); fill!(frw, zero(FT))
        _launch!(_dmr_frw_kernel!, frw, mask_bgc_soilc, col_dz, nlev_soildecomp_standard)

        # Water scalar — per-column accumulation over standard soil levels (fr[c,j]
        # = col_dz/frw computed in-thread; ascending-j order preserved, byte-identical).
        _launch!(_dmr_wscalar1_kernel!, w_scalar, mask_bgc_soilc, soilpsi, col_dz,
                 frw, nlev_soildecomp_standard, minpsi_ft, maxpsi_ft; ndrange = nc)

        # O2 scalar
        if anoxia
            _launch!(_dmr_oscalar1_anoxia_kernel!, o_scalar, mask_bgc_soilc,
                     o2stress_unsat, col_dz, frw, nlev_soildecomp_standard, mino2lim_ft;
                     ndrange = nc)
        else
            _launch!(_dmr_oscalar_one_kernel!, o_scalar, nlevdecomp; ndrange = nc)
        end

    else
        # --- Multi-level decomposition ---

        # Water scalar (2D kernel; each (c,j) independent)
        mimics_water_scalar!(w_scalar, mask_bgc_soilc, soilpsi, minpsi, maxpsi;
                             nc = nc, nlevdecomp = nlevdecomp)

        # O2 scalar
        if anoxia
            mimics_oscalar_anoxia!(o_scalar, mask_bgc_soilc, o2stress_unsat, mino2lim;
                                   nc = nc, nlevdecomp = nlevdecomp)
        else
            _launch!(_dmr_oscalar_one_kernel!, o_scalar, nlevdecomp; ndrange = nc)
        end
    end

    # Depth scalar — currently set to 1.0 (placeholder, as in Fortran)
    # (2D kernel; each (c,j) independent)
    mimics_depth_scalar!(depth_scalar, mask_bgc_soilc;
                         nc = nc, nlevdecomp = nlevdecomp)

    # Bundle the column/level arrays the rate + post kernels touch (one Adapt-moved
    # struct keeps us under the Metal ~31-arg limit; writes flow back into cf).
    arr = _dmr_arr(cf, mimics_state, depth_scalar, t_soisno, decomp_cpools_vr,
                   col_dz, ligninNratioAvg, annsum_npp_col)

    # Bundle the scalar coefficients at the working precision (a Float64 scalar
    # kernel arg is invalid Metal IR; on Float64 FT(x) === x byte-identically).
    coef = _DMRCoef{FT}(;
        eps_val               = FT(eps_val),
        mimics_fmet_p1        = FT(params.mimics_fmet[1]),
        mimics_fmet_p2        = FT(params.mimics_fmet[2]),
        mimics_fmet_p3        = FT(params.mimics_fmet[3]),
        mimics_fmet_p4        = FT(params.mimics_fmet[4]),
        mimics_fchem_r_p1     = FT(params.mimics_fchem_r[1]),
        mimics_fchem_r_p2     = FT(params.mimics_fchem_r[2]),
        mimics_fchem_k_p1     = FT(params.mimics_fchem_k[1]),
        mimics_fchem_k_p2     = FT(params.mimics_fchem_k[2]),
        mimics_tau_mod_min    = FT(params.mimics_tau_mod_min),
        mimics_tau_mod_max    = FT(params.mimics_tau_mod_max),
        mimics_tau_mod_factor = FT(params.mimics_tau_mod_factor),
        mimics_tau_r_p1       = FT(params.mimics_tau_r[1]),
        mimics_tau_r_p2       = FT(params.mimics_tau_r[2]),
        mimics_tau_k_p1       = FT(params.mimics_tau_k[1]),
        mimics_tau_k_p2       = FT(params.mimics_tau_k[2]),
        mimics_ko_r           = FT(params.mimics_ko_r),
        mimics_ko_k           = FT(params.mimics_ko_k),
        mimics_densdep        = FT(params.mimics_densdep),
        mimics_desorpQ10      = FT(params.mimics_desorpQ10),
        mimics_t_soi_ref      = FT(params.mimics_t_soi_ref),
        mimics_cn_mod_num     = FT(params.mimics_cn_mod_num),
        mimics_cn_r           = FT(params.mimics_cn_r),
        mimics_cn_k           = FT(params.mimics_cn_k),
        vslope_l1_m1 = FT(mimics_state.vslope_l1_m1), vslope_l2_m1 = FT(mimics_state.vslope_l2_m1),
        vslope_s1_m1 = FT(mimics_state.vslope_s1_m1), vslope_l1_m2 = FT(mimics_state.vslope_l1_m2),
        vslope_l2_m2 = FT(mimics_state.vslope_l2_m2), vslope_s1_m2 = FT(mimics_state.vslope_s1_m2),
        kslope_l1_m1 = FT(mimics_state.kslope_l1_m1), kslope_l2_m1 = FT(mimics_state.kslope_l2_m1),
        kslope_s1_m1 = FT(mimics_state.kslope_s1_m1), kslope_l1_m2 = FT(mimics_state.kslope_l1_m2),
        kslope_l2_m2 = FT(mimics_state.kslope_l2_m2), kslope_s1_m2 = FT(mimics_state.kslope_s1_m2),
        vint_l1_m1 = FT(mimics_state.vint_l1_m1), vint_l2_m1 = FT(mimics_state.vint_l2_m1),
        vint_s1_m1 = FT(mimics_state.vint_s1_m1), vint_l1_m2 = FT(mimics_state.vint_l1_m2),
        vint_l2_m2 = FT(mimics_state.vint_l2_m2), vint_s1_m2 = FT(mimics_state.vint_s1_m2),
        kint_l1_m1 = FT(mimics_state.kint_l1_m1), kint_l2_m1 = FT(mimics_state.kint_l2_m1),
        kint_s1_m1 = FT(mimics_state.kint_s1_m1), kint_l1_m2 = FT(mimics_state.kint_l1_m2),
        kint_l2_m2 = FT(mimics_state.kint_l2_m2), kint_s1_m2 = FT(mimics_state.kint_s1_m2),
        vmod_l1_m1 = FT(mimics_state.vmod_l1_m1), vmod_l2_m1 = FT(mimics_state.vmod_l2_m1),
        vmod_s1_m1 = FT(mimics_state.vmod_s1_m1), vmod_l1_m2 = FT(mimics_state.vmod_l1_m2),
        vmod_l2_m2 = FT(mimics_state.vmod_l2_m2), vmod_s1_m2 = FT(mimics_state.vmod_s1_m2),
        kmod_l1_m1 = FT(mimics_state.kmod_l1_m1), kmod_l2_m1 = FT(mimics_state.kmod_l2_m1),
        kmod_s1_m1 = FT(mimics_state.kmod_s1_m1), kmod_l1_m2 = FT(mimics_state.kmod_l1_m2),
        kmod_l2_m2 = FT(mimics_state.kmod_l2_m2), kmod_s1_m2 = FT(mimics_state.kmod_s1_m2),
        k_frag    = FT(k_frag),
        secsphr   = FT(SECSPHR),
        tfrz      = FT(TFRZ),
        g_to_mg   = FT(G_TO_MG),
        cm3_to_m3 = FT(CM3_TO_M3))

    # --- Calculate rates for all litter and SOM pools (per-column kernel, internal j loop) ---
    _launch!(_dmr_rates_kernel!, decomp_k, arr, mask_bgc_soilc, idx, coef, nlevdecomp,
             use_fates; ndrange = nc)

    # pathfrac terms not calculated in the rate kernel (own-[c,j], all levels).
    _launch!(_dmr_pathfrac_post_kernel!, pathfrac_decomp_cascade, arr, idx,
             nlevdecomp; ndrange = nc)
    if !use_fates
        i_cwdl2_local = 15
        _launch!(_dmr_cwdl2_kernel!, pathfrac_decomp_cascade, arr, i_cwdl2_local,
                 nlevdecomp, FT(rf_cwdl2); ndrange = nc)
    end

    # Respiration fractions for all transitions (own-[c,j], all levels).
    _launch!(_dmr_rf_kernel!, rf_decomp_cascade, arr, idx, nlevdecomp,
             FT(mimics_state.rf_l1m1), FT(mimics_state.rf_l1m2),
             FT(mimics_state.rf_l2m1), FT(mimics_state.rf_l2m2),
             FT(mimics_state.rf_s1m1), FT(mimics_state.rf_s1m2); ndrange = nc)

    return nothing
end
