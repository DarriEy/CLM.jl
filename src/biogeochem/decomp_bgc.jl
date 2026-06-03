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
Base.@kwdef mutable struct DecompCascadeConData{FT<:Real}
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
Base.@kwdef mutable struct DecompBGCParams{FT<:Real}
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
Base.@kwdef mutable struct DecompBGCState{FT<:Real, M<:AbstractMatrix{FT}}
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

    # Spatially-varying respiration fractions (col × nlevdecomp). Array-type param M
    # stays LOOSE so adapt(MtlArray/CuArray, ·) can swap device storage in (and AD can
    # swap Dual arrays); pinning ::Matrix{FT} breaks device-movability of CLMInstances.
    rf_s1s2   ::M = Matrix{Float64}(undef, 0, 0)
    rf_s1s3   ::M = Matrix{Float64}(undef, 0, 0)

    # Path fractions
    f_s1s2    ::M = Matrix{Float64}(undef, 0, 0)
    f_s1s3    ::M = Matrix{Float64}(undef, 0, 0)
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

# The `{FT}` convenience ctor defaults M to Matrix{FT} at the working precision; @kwdef
# still gives the no-FT `DecompBGCState()` (infers Float64 from the defaults).
DecompBGCState{FT}(; kwargs...) where {FT<:Real} =
    DecompBGCState{FT, Matrix{FT}}(; kwargs...)

# Device-movable: lets the spatially-varying f_s1s2/rf_s1s2/… arrays ride to the GPU
# for decomp_rate_constants_bgc! (scalar/Int/Bool fields pass through unchanged).
Adapt.@adapt_structure DecompBGCState

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
                                  tau_l1::Real,
                                  tau_l2_l3::Real,
                                  tau_s1::Real,
                                  tau_s2::Real,
                                  tau_s3::Real,
                                  cn_s1::Real,
                                  cn_s2::Real,
                                  cn_s3::Real,
                                  rf_l1s1::Real,
                                  rf_l2s1::Real,
                                  rf_l3s2::Real,
                                  rf_s2s1::Real,
                                  rf_s2s3::Real,
                                  rf_s3s1::Real,
                                  rf_cwdl3::Real,
                                  cwd_fcel::Real,
                                  bgc_initial_Cstocks::Vector{<:Real},
                                  bgc_initial_Cstocks_depth::Real)
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
# Kernels for fully-independent per-element loops in this file.
# ---------------------------------------------------------------------------

# Sand-dependent path/respiration fractions (per column AND decomp level — 2D).
# Each (c,j) is independent. Four output arrays share the intermediate `t`;
# the backend is taken from the first (f_s1s2). All written arrays passed plain.
@kernel function _decompb_sandfrac_kernel!(f_s1s2, f_s1s3, rf_s1s2, rf_s1s3,
                                           @Const(cellsand))
    c, j = @index(Global, NTuple)
    @inbounds begin
        T = eltype(f_s1s2)
        t = T(0.85) - T(0.68) * T(0.01) * (T(100.0) - cellsand[c, j])
        f_s1s2[c, j]  = one(T) - T(0.004) / (one(T) - t)
        f_s1s3[c, j]  = T(0.004) / (one(T) - t)
        rf_s1s2[c, j] = t
        rf_s1s3[c, j] = t
    end
end

"""
    decompb_sandfrac!(f_s1s2, f_s1s3, rf_s1s2, rf_s1s3, cellsand, nc, nlevdecomp)

Sand-dependent S1->S2/S3 path and respiration fractions over all
(column, decomp level) pairs. Backend-agnostic 2D kernel.
"""
decompb_sandfrac!(f_s1s2, f_s1s3, rf_s1s2, rf_s1s3, cellsand, nc::Int, nlevdecomp::Int) =
    _launch!(_decompb_sandfrac_kernel!, f_s1s2, f_s1s3, rf_s1s2, rf_s1s3, cellsand;
             ndrange = (nc, nlevdecomp))

# Multi-level Q10 temperature scalar (per column AND decomp level — 2D, masked).
# Each (c,j) independent; freezing branch uses froz_q10.
@kernel function _decompb_tscalar_q10_kernel!(t_scalar, @Const(mask), @Const(t_soisno),
                                              Q10, froz_q10, tfrz)
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        T = eltype(t_scalar)
        ts = t_soisno[c, j]
        if ts >= tfrz
            t_scalar[c, j] = Q10^((ts - (tfrz + T(25.0))) / T(10.0))
        else
            t_scalar[c, j] = (Q10^(T(-25.0) / T(10.0))) * (froz_q10^((ts - tfrz) / T(10.0)))
        end
    end
end

"""
    decompb_tscalar_q10!(t_scalar, mask, t_soisno, Q10, froz_q10, tfrz, nc, nlevdecomp)

Q10 temperature scalar for multi-level decomposition. Backend-agnostic 2D kernel.
"""
function decompb_tscalar_q10!(t_scalar, mask, t_soisno, Q10, froz_q10, tfrz, nc::Int, nlevdecomp::Int)
    _T = eltype(t_scalar)   # convert scalar args to working precision (no Float64 on Metal)
    _launch!(_decompb_tscalar_q10_kernel!, t_scalar, mask, t_soisno,
             _T(Q10), _T(froz_q10), _T(tfrz); ndrange = (nc, nlevdecomp))
end

# Multi-level water scalar (per column AND decomp level — 2D, masked).
@kernel function _decompb_wscalar_kernel!(w_scalar, @Const(mask), @Const(soilpsi),
                                          minpsi, maxpsi)
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        psi = min(soilpsi[c, j], maxpsi)
        if psi > minpsi
            w_scalar[c, j] = log(minpsi / psi) / log(minpsi / maxpsi)
        else
            w_scalar[c, j] = zero(eltype(w_scalar))
        end
    end
end

"""
    decompb_wscalar!(w_scalar, mask, soilpsi, minpsi, maxpsi, nc, nlevdecomp)

Soil-water-potential scalar for multi-level decomposition. Backend-agnostic 2D kernel.
"""
function decompb_wscalar!(w_scalar, mask, soilpsi, minpsi, maxpsi, nc::Int, nlevdecomp::Int)
    _T = eltype(w_scalar)
    _launch!(_decompb_wscalar_kernel!, w_scalar, mask, soilpsi, _T(minpsi), _T(maxpsi);
             ndrange = (nc, nlevdecomp))
end

# Multi-level anoxia O2 scalar (per column AND decomp level — 2D, masked).
@kernel function _decompb_oscalar_anox_kernel!(o_scalar, @Const(mask),
                                               @Const(o2stress_unsat), mino2lim)
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        o_scalar[c, j] = max(o2stress_unsat[c, j], mino2lim)
    end
end

"""
    decompb_oscalar_anox!(o_scalar, mask, o2stress_unsat, mino2lim, nc, nlevdecomp)

O2 stress scalar (anoxia path) for multi-level decomposition. Backend-agnostic.
"""
function decompb_oscalar_anox!(o_scalar, mask, o2stress_unsat, mino2lim, nc::Int, nlevdecomp::Int)
    _launch!(_decompb_oscalar_anox_kernel!, o_scalar, mask, o2stress_unsat,
             eltype(o_scalar)(mino2lim); ndrange = (nc, nlevdecomp))
end

# Unconditional unit O2 scalar (per column AND decomp level — 2D, no mask).
@kernel function _decompb_oscalar_one_kernel!(o_scalar)
    c, j = @index(Global, NTuple)
    @inbounds o_scalar[c, j] = one(eltype(o_scalar))
end

"""
    decompb_oscalar_one!(o_scalar, nc, nlevdecomp)

Set the O2 scalar to 1.0 for all (column, decomp level). Backend-agnostic.
"""
decompb_oscalar_one!(o_scalar, nc::Int, nlevdecomp::Int) =
    _launch!(_decompb_oscalar_one_kernel!, o_scalar; ndrange = (nc, nlevdecomp))

# Normalize t_scalar by a precomputed scalar factor (per column AND level — 2D, masked).
@kernel function _decompb_tscalar_norm_kernel!(t_scalar, @Const(mask), factor)
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        t_scalar[c, j] *= factor
    end
end

"""
    decompb_tscalar_norm!(t_scalar, mask, factor, nc, nlevdecomp)

Multiply t_scalar by the Q10->CENTURY normalization factor. Backend-agnostic.
"""
decompb_tscalar_norm!(t_scalar, mask, factor, nc::Int, nlevdecomp::Int) =
    _launch!(_decompb_tscalar_norm_kernel!, t_scalar, mask, eltype(t_scalar)(factor);
             ndrange = (nc, nlevdecomp))

# Depth scalar — fixed e-folding depth (per column AND decomp level — 2D, masked).
@kernel function _decompb_depthscalar_kernel!(depth_scalar, @Const(mask),
                                              @Const(zsoi_vals), efolding)
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        depth_scalar[c, j] = exp(-zsoi_vals[j] / efolding)
    end
end

"""
    decompb_depthscalar!(depth_scalar, mask, zsoi_vals, efolding, nc, nlevdecomp)

Fixed e-folding depth scalar for multi-level decomposition. Backend-agnostic.
"""
decompb_depthscalar!(depth_scalar, mask, zsoi_vals, efolding, nc::Int, nlevdecomp::Int) =
    _launch!(_decompb_depthscalar_kernel!, depth_scalar, mask, zsoi_vals,
             eltype(depth_scalar)(efolding); ndrange = (nc, nlevdecomp))

# ---------------------------------------------------------------------------
# Kernels for the remaining per-column / per-(column,level) loops of
# decomp_rate_constants_bgc!.  These keep the function whole-fn on the GPU.
# ---------------------------------------------------------------------------

# Device-safe spinup latitude term (mirrors get_spinup_latitude_term, with all
# Float64 literals converted to the working eltype T so it compiles under Metal).
@inline _spinup_lat_term(lat, ::Type{T}) where {T} =
    one(T) + T(50.0) / (one(T) + exp(T(-0.15) * (abs(T(lat)) - T(60.0))))

# Spinup geographic terms (per column).  Each thread owns its column c and writes
# its own [c] entries — race-free.  Only the columns whose spinup_factor differs
# from 1 by more than eps_val get a non-unit term (others keep the preset 1.0).
@kernel function _decompb_spinup_geogterm_kernel!(
        spinup_geogterm_l1, spinup_geogterm_l23, spinup_geogterm_cwd,
        spinup_geogterm_s1, spinup_geogterm_s2, spinup_geogterm_s3,
        @Const(mask), @Const(col_gridcell), @Const(latdeg),
        sf_l1, sf_l23, sf_cwd, sf_s1, sf_s2, sf_s3, eps_val, use_fates)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(spinup_geogterm_l1)
        lat = latdeg[col_gridcell[c]]
        slt = _spinup_lat_term(lat, T)
        one_T = one(T)
        if abs(sf_l1 - one_T) > eps_val
            spinup_geogterm_l1[c] = sf_l1 * slt
        end
        if abs(sf_l23 - one_T) > eps_val
            spinup_geogterm_l23[c] = sf_l23 * slt
        end
        if !use_fates
            if abs(sf_cwd - one_T) > eps_val
                spinup_geogterm_cwd[c] = sf_cwd * slt
            end
        end
        if abs(sf_s1 - one_T) > eps_val
            spinup_geogterm_s1[c] = sf_s1 * slt
        end
        if abs(sf_s2 - one_T) > eps_val
            spinup_geogterm_s2[c] = sf_s2 * slt
        end
        if abs(sf_s3 - one_T) > eps_val
            spinup_geogterm_s3[c] = sf_s3 * slt
        end
    end
end

decompb_spinup_geogterm!(gl1, gl23, gcwd, gs1, gs2, gs3, mask, col_gridcell, latdeg,
                         sf_l1, sf_l23, sf_cwd, sf_s1, sf_s2, sf_s3, eps_val, use_fates) =
    _launch!(_decompb_spinup_geogterm_kernel!, gl1, gl23, gcwd, gs1, gs2, gs3,
             mask, col_gridcell, latdeg,
             sf_l1, sf_l23, sf_cwd, sf_s1, sf_s2, sf_s3, eps_val, use_fates;
             ndrange = length(gl1))

# Single-level rooting-fraction weights (per column).  Thread c sums col_dz over
# the 5 standard levels into its own frw[c], then fills fr[c,1..5] — race-free.
@kernel function _decompb_singlelev_fr_kernel!(frw, fr, @Const(mask), @Const(col_dz),
                                               nlev_std)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(fr)
        s = zero(T)
        for j in 1:nlev_std
            s += col_dz[c, j]
        end
        frw[c] = s
        for j in 1:nlev_std
            if s != zero(T)
                fr[c, j] = col_dz[c, j] / s
            else
                fr[c, j] = zero(T)
            end
        end
    end
end

decompb_singlelev_fr!(frw, fr, mask, col_dz, nlev_std) =
    _launch!(_decompb_singlelev_fr_kernel!, frw, fr, mask, col_dz, nlev_std;
             ndrange = length(frw))

# Single-level Q10 temperature scalar (per column; sum over 5 standard levels into
# own t_scalar[c,1]).  Ascending-j order preserved — race-free.
@kernel function _decompb_singlelev_tscalar_q10_kernel!(
        t_scalar, @Const(mask), @Const(t_soisno), @Const(fr),
        Q10, froz_q10, tfrz, nlev_std)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(t_scalar)
        acc = zero(T)
        for j in 1:nlev_std
            ts = t_soisno[c, j]
            if ts >= tfrz
                acc += (Q10^((ts - (tfrz + T(25.0))) / T(10.0))) * fr[c, j]
            else
                acc += (Q10^(T(-25.0) / T(10.0))) *
                       (froz_q10^((ts - tfrz) / T(10.0))) * fr[c, j]
            end
        end
        t_scalar[c, 1] = acc
    end
end

decompb_singlelev_tscalar_q10!(t_scalar, mask, t_soisno, fr, Q10, froz_q10, tfrz, nlev_std) =
    _launch!(_decompb_singlelev_tscalar_q10_kernel!, t_scalar, mask, t_soisno, fr,
             Q10, froz_q10, tfrz, nlev_std; ndrange = length(mask))

# Single-level CENTURY arctangent temperature scalar (per column; sum over levels).
@kernel function _decompb_singlelev_tscalar_century_kernel!(
        t_scalar, @Const(mask), @Const(t_soisno), @Const(fr),
        catanf_30, tfrz, rpi, nlev_std)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(t_scalar)
        acc = zero(T)
        for j in 1:nlev_std
            t1 = t_soisno[c, j] - tfrz
            cat = T(11.75) + (T(29.7) / rpi) * atan(rpi * T(0.031) * (t1 - T(15.4)))
            acc += max(cat / catanf_30 * fr[c, j], T(0.01))
        end
        t_scalar[c, 1] = acc
    end
end

decompb_singlelev_tscalar_century!(t_scalar, mask, t_soisno, fr, catanf_30, tfrz, rpi, nlev_std) =
    _launch!(_decompb_singlelev_tscalar_century_kernel!, t_scalar, mask, t_soisno, fr,
             catanf_30, tfrz, rpi, nlev_std; ndrange = length(mask))

# Single-level water scalar (per column; sum over levels into own w_scalar[c,1]).
@kernel function _decompb_singlelev_wscalar_kernel!(
        w_scalar, @Const(mask), @Const(soilpsi), @Const(fr), minpsi, maxpsi, nlev_std)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(w_scalar)
        acc = zero(T)
        for j in 1:nlev_std
            psi = min(soilpsi[c, j], maxpsi)
            if psi > minpsi
                acc += (log(minpsi / psi) / log(minpsi / maxpsi)) * fr[c, j]
            end
        end
        w_scalar[c, 1] = acc
    end
end

decompb_singlelev_wscalar!(w_scalar, mask, soilpsi, fr, minpsi, maxpsi, nlev_std) =
    _launch!(_decompb_singlelev_wscalar_kernel!, w_scalar, mask, soilpsi, fr,
             minpsi, maxpsi, nlev_std; ndrange = length(mask))

# Single-level anoxia O2 scalar (per column; sum over levels into own o_scalar[c,1]).
@kernel function _decompb_singlelev_oscalar_anox_kernel!(
        o_scalar, @Const(mask), @Const(o2stress_unsat), @Const(fr), mino2lim, nlev_std)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(o_scalar)
        acc = zero(T)
        for j in 1:nlev_std
            acc += fr[c, j] * max(o2stress_unsat[c, j], mino2lim)
        end
        o_scalar[c, 1] = acc
    end
end

decompb_singlelev_oscalar_anox!(o_scalar, mask, o2stress_unsat, fr, mino2lim, nlev_std) =
    _launch!(_decompb_singlelev_oscalar_anox_kernel!, o_scalar, mask, o2stress_unsat, fr,
             mino2lim, nlev_std; ndrange = length(mask))

# Multi-level CENTURY arctangent temperature scalar (per (column, level), masked).
@kernel function _decompb_tscalar_century_kernel!(
        t_scalar, @Const(mask), @Const(t_soisno), catanf_30, tfrz, rpi)
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        T = eltype(t_scalar)
        t1 = t_soisno[c, j] - tfrz
        cat = T(11.75) + (T(29.7) / rpi) * atan(rpi * T(0.031) * (t1 - T(15.4)))
        t_scalar[c, j] = max(cat / catanf_30, T(0.01))
    end
end

decompb_tscalar_century!(t_scalar, mask, t_soisno, catanf_30, tfrz, rpi, nc::Int, nlevdecomp::Int) =
    _launch!(_decompb_tscalar_century_kernel!, t_scalar, mask, t_soisno, catanf_30, tfrz, rpi;
             ndrange = (nc, nlevdecomp))

# Final rate-constant assembly (per (column, level), masked).  Each (c,j) writes
# its own decomp_k entries — independent.
@kernel function _decompb_decompk_kernel!(
        decomp_k, @Const(mask), @Const(t_scalar), @Const(w_scalar),
        @Const(depth_scalar), @Const(o_scalar),
        @Const(geo_l1), @Const(geo_l23), @Const(geo_cwd),
        @Const(geo_s1), @Const(geo_s2), @Const(geo_s3),
        k_l1, k_l2_l3, k_s1, k_s2, k_s3, k_frag,
        i_met_lit, i_cel_lit, i_lig_lit, i_act_som, i_slo_som, i_pas_som,
        i_cwd_local, use_fates)
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        rate = t_scalar[c, j] * w_scalar[c, j] * depth_scalar[c, j] * o_scalar[c, j]
        decomp_k[c, j, i_met_lit] = k_l1    * rate * geo_l1[c]
        decomp_k[c, j, i_cel_lit] = k_l2_l3 * rate * geo_l23[c]
        decomp_k[c, j, i_lig_lit] = k_l2_l3 * rate * geo_l23[c]
        decomp_k[c, j, i_act_som] = k_s1    * rate * geo_s1[c]
        decomp_k[c, j, i_slo_som] = k_s2    * rate * geo_s2[c]
        decomp_k[c, j, i_pas_som] = k_s3    * rate * geo_s3[c]
        if !use_fates
            decomp_k[c, j, i_cwd_local] = k_frag * rate * geo_cwd[c]
        end
    end
end

decompb_decompk!(decomp_k, mask, t_scalar, w_scalar, depth_scalar, o_scalar,
                 geo_l1, geo_l23, geo_cwd, geo_s1, geo_s2, geo_s3,
                 k_l1, k_l2_l3, k_s1, k_s2, k_s3, k_frag,
                 i_met_lit, i_cel_lit, i_lig_lit, i_act_som, i_slo_som, i_pas_som,
                 i_cwd_local, use_fates, nc::Int, nlevdecomp::Int) =
    _launch!(_decompb_decompk_kernel!, decomp_k, mask, t_scalar, w_scalar, depth_scalar, o_scalar,
             geo_l1, geo_l23, geo_cwd, geo_s1, geo_s2, geo_s3,
             k_l1, k_l2_l3, k_s1, k_s2, k_s3, k_frag,
             i_met_lit, i_cel_lit, i_lig_lit, i_act_som, i_slo_som, i_pas_som,
             i_cwd_local, use_fates; ndrange = (nc, nlevdecomp))

# Pathfrac / respiration-fraction cascade assembly (per (column, level)).  No mask:
# the original host loop writes all columns (mirrors the unmasked host loop).
@kernel function _decompb_pathfrac_rf_kernel!(
        pathfrac, rf, @Const(f_s1s2), @Const(f_s1s3), @Const(rf_s1s2), @Const(rf_s1s3),
        f_s2s1, f_s2s3, rf_l1s1, rf_l2s1, rf_l3s2, rf_s2s1, rf_s2s3, rf_s3s1,
        i_l1s1, i_l2s1, i_l3s2, i_s1s2, i_s1s3, i_s2s1, i_s2s3, i_s3s1)
    c, j = @index(Global, NTuple)
    @inbounds begin
        T = eltype(pathfrac)
        one_T = one(T)
        pathfrac[c, j, i_l1s1] = one_T
        pathfrac[c, j, i_l2s1] = one_T
        pathfrac[c, j, i_l3s2] = one_T
        pathfrac[c, j, i_s1s2] = f_s1s2[c, j]
        pathfrac[c, j, i_s1s3] = f_s1s3[c, j]
        pathfrac[c, j, i_s2s1] = f_s2s1
        pathfrac[c, j, i_s2s3] = f_s2s3
        pathfrac[c, j, i_s3s1] = one_T

        rf[c, j, i_l1s1] = rf_l1s1
        rf[c, j, i_l2s1] = rf_l2s1
        rf[c, j, i_l3s2] = rf_l3s2
        rf[c, j, i_s1s2] = rf_s1s2[c, j]
        rf[c, j, i_s1s3] = rf_s1s3[c, j]
        rf[c, j, i_s2s1] = rf_s2s1
        rf[c, j, i_s2s3] = rf_s2s3
        rf[c, j, i_s3s1] = rf_s3s1
    end
end

decompb_pathfrac_rf!(pathfrac, rf, f_s1s2, f_s1s3, rf_s1s2, rf_s1s3,
                     f_s2s1, f_s2s3, rf_l1s1, rf_l2s1, rf_l3s2, rf_s2s1, rf_s2s3, rf_s3s1,
                     i_l1s1, i_l2s1, i_l3s2, i_s1s2, i_s1s3, i_s2s1, i_s2s3, i_s3s1,
                     nc::Int, nlevdecomp::Int) =
    _launch!(_decompb_pathfrac_rf_kernel!, pathfrac, rf, f_s1s2, f_s1s3, rf_s1s2, rf_s1s3,
             f_s2s1, f_s2s3, rf_l1s1, rf_l2s1, rf_l3s2, rf_s2s1, rf_s2s3, rf_s3s1,
             i_l1s1, i_l2s1, i_l3s2, i_s1s2, i_s1s3, i_s2s1, i_s2s3, i_s3s1;
             ndrange = (nc, nlevdecomp))

# CWD->litter pathfrac / rf assembly (per (column, level), no mask).
@kernel function _decompb_pathfrac_rf_cwd_kernel!(
        pathfrac, rf, cwd_fcel, cwd_flig, rf_cwdl2, rf_cwdl3, i_cwdl2_local, i_cwdl3)
    c, j = @index(Global, NTuple)
    @inbounds begin
        pathfrac[c, j, i_cwdl2_local] = cwd_fcel
        pathfrac[c, j, i_cwdl3]       = cwd_flig
        rf[c, j, i_cwdl2_local]       = rf_cwdl2
        rf[c, j, i_cwdl3]             = rf_cwdl3
    end
end

decompb_pathfrac_rf_cwd!(pathfrac, rf, cwd_fcel, cwd_flig, rf_cwdl2, rf_cwdl3,
                         i_cwdl2_local, i_cwdl3, nc::Int, nlevdecomp::Int) =
    _launch!(_decompb_pathfrac_rf_cwd_kernel!, pathfrac, rf, cwd_fcel, cwd_flig,
             rf_cwdl2, rf_cwdl3, i_cwdl2_local, i_cwdl3; ndrange = (nc, nlevdecomp))

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
- `cellsand::Matrix{<:Real}`: column sand fraction (col × nlevdecomp)
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
                                   cellsand::Matrix{<:Real},
                                   bounds::UnitRange{Int},
                                   nlevdecomp::Int,
                                   ndecomp_pools_max::Int,
                                   ndecomp_cascade_transitions_max::Int=10,
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
    bgc_state.rf_s1s2 = zeros(FT, nc, nlevdecomp)
    bgc_state.rf_s1s3 = zeros(FT, nc, nlevdecomp)
    bgc_state.f_s1s2  = zeros(FT, nc, nlevdecomp)
    bgc_state.f_s1s3  = zeros(FT, nc, nlevdecomp)

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
    decompb_sandfrac!(bgc_state.f_s1s2, bgc_state.f_s1s3,
                      bgc_state.rf_s1s2, bgc_state.rf_s1s3,
                      cellsand, nc, nlevdecomp)
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
- `t_soisno::Matrix{<:Real}`: soil temperature (K), (col × nlev)
- `soilpsi::Matrix{<:Real}`: soil water potential (MPa), (col × nlev)
- `days_per_year::Float64`: days per year
- `dt::Float64`: timestep (seconds)
- `zsoi_vals::Vector{<:Real}`: soil layer depths (m)
- `spinup_state::Int`: spinup state (0=normal, >=1=accelerated)
- `use_lch4::Bool`: whether LCH4 model is active
- `anoxia::Bool`: whether anoxia is enabled
- `use_fates::Bool`: whether FATES is enabled
- `o2stress_unsat::Matrix{<:Real}`: O2 stress ratio (col × nlev), for anoxia
- `col_dz::Matrix{<:Real}`: column layer thicknesses (col × nlev), for nlevdecomp==1
- `col_gridcell::Vector{Int}`: column-to-gridcell mapping, for spinup
- `latdeg::Vector{<:Real}`: gridcell latitudes (degrees), for spinup
"""
function decomp_rate_constants_bgc!(cf::SoilBiogeochemCarbonFluxData,
                                     bgc_state::DecompBGCState,
                                     params::DecompBGCParams,
                                     cn_params::CNSharedParamsData,
                                     cascade_con::DecompCascadeConData;
                                     mask_bgc_soilc::AbstractVector{Bool},
                                     bounds::UnitRange{Int},
                                     nlevdecomp::Int,
                                     t_soisno::AbstractMatrix{<:Real},
                                     soilpsi::AbstractMatrix{<:Real},
                                     days_per_year::Real,
                                     dt::Real,
                                     zsoi_vals::AbstractVector{<:Real},
                                     spinup_state::Int=0,
                                     use_lch4::Bool=false,
                                     anoxia::Bool=false,
                                     use_fates::Bool=false,
                                     o2stress_unsat::AbstractMatrix{<:Real}=Matrix{Float64}(undef, 0, 0),
                                     col_dz::AbstractMatrix{<:Real}=Matrix{Float64}(undef, 0, 0),
                                     col_gridcell::AbstractVector{<:Integer}=Int[],
                                     latdeg::AbstractVector{<:Real}=Float64[])

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
    # Allocated like a state array (t_soisno) so they live on the state's backend
    # (device on GPU): they are passed to _decompb_decompk_kernel!, and a host ones()
    # vector would be a non-bitstype kernel argument on Metal.
    FT = eltype(t_soisno)
    _geo() = fill!(similar(t_soisno, FT, nc), one(FT))
    spinup_geogterm_l1  = _geo()
    spinup_geogterm_l23 = _geo()
    spinup_geogterm_cwd = _geo()
    spinup_geogterm_s1  = _geo()
    spinup_geogterm_s2  = _geo()
    spinup_geogterm_s3  = _geo()

    if spinup_state >= 1
        i_cwd_local = i_pas_som + 1
        # spinup_factor entries are host scalars; pull them out and convert to the
        # working eltype so they are valid kernel scalar args (Metal rejects raw
        # Float64). use_fates guards the CWD pool index (i_cwd_local valid only when
        # CWD exists); index it conditionally to stay in bounds.
        sf_l1  = FT(spinup_factor[i_met_lit])
        sf_l23 = FT(spinup_factor[i_cel_lit])
        sf_cwd = use_fates ? one(FT) : FT(spinup_factor[i_cwd_local])
        sf_s1  = FT(spinup_factor[i_act_som])
        sf_s2  = FT(spinup_factor[i_slo_som])
        sf_s3  = FT(spinup_factor[i_pas_som])
        decompb_spinup_geogterm!(spinup_geogterm_l1, spinup_geogterm_l23,
                                 spinup_geogterm_cwd, spinup_geogterm_s1,
                                 spinup_geogterm_s2, spinup_geogterm_s3,
                                 mask_bgc_soilc, col_gridcell, latdeg,
                                 sf_l1, sf_l23, sf_cwd, sf_s1, sf_s2, sf_s3,
                                 FT(eps_val), use_fates)
    end

    # --- Time-dependent coefficients ---
    t_scalar = cf.t_scalar_col
    w_scalar = cf.w_scalar_col
    o_scalar = cf.o_scalar_col
    decomp_k = cf.decomp_k_col

    FT_rate = eltype(t_soisno)
    # Scratch allocated like a state array (t_soisno) so it lands on the state's backend
    # (device on GPU) — zeros() would put it on the host and a device kernel writing it
    # trips scalar indexing.
    depth_scalar = fill!(similar(t_soisno, FT_rate, nc, nlevdecomp), zero(FT_rate))

    if nlevdecomp == 1
        # ---- Single-level decomposition ----
        # Weight temperature and water potential scalars by rooting fraction
        nlev_soildecomp_standard = 5
        frw = fill!(similar(t_soisno, FT_rate, nc), zero(FT_rate))
        fr  = fill!(similar(t_soisno, FT_rate, nc, nlev_soildecomp_standard), zero(FT_rate))

        # Rooting-fraction weights (per column; in-thread sum over the 5 levels).
        decompb_singlelev_fr!(frw, fr, mask_bgc_soilc, col_dz, nlev_soildecomp_standard)

        if !bgc_state.use_century_tfunc
            # Q10 temperature scalar
            decompb_singlelev_tscalar_q10!(t_scalar, mask_bgc_soilc, t_soisno, fr,
                                           FT_rate(Q10), FT_rate(froz_q10), FT_rate(TFRZ),
                                           nlev_soildecomp_standard)
        else
            # CENTURY arctangent temperature function
            decompb_singlelev_tscalar_century!(t_scalar, mask_bgc_soilc, t_soisno, fr,
                                               FT_rate(catanf_30), FT_rate(TFRZ), FT_rate(RPI),
                                               nlev_soildecomp_standard)
        end

        # Water scalar
        decompb_singlelev_wscalar!(w_scalar, mask_bgc_soilc, soilpsi, fr,
                                   FT_rate(minpsi), FT_rate(maxpsi), nlev_soildecomp_standard)

        # O2 scalar
        if use_lch4 && anoxia
            decompb_singlelev_oscalar_anox!(o_scalar, mask_bgc_soilc, o2stress_unsat, fr,
                                            FT_rate(mino2lim), nlev_soildecomp_standard)
        else
            decompb_oscalar_one!(o_scalar, nc, nlevdecomp)
        end

    else
        # ---- Multi-level decomposition ----

        if !bgc_state.use_century_tfunc
            # Q10 temperature scalar
            decompb_tscalar_q10!(t_scalar, mask_bgc_soilc, t_soisno,
                                 Q10, froz_q10, TFRZ, nc, nlevdecomp)
        else
            # CENTURY arctangent temperature function
            decompb_tscalar_century!(t_scalar, mask_bgc_soilc, t_soisno,
                                     FT_rate(catanf_30), FT_rate(TFRZ), FT_rate(RPI),
                                     nc, nlevdecomp)
        end

        # Water scalar
        decompb_wscalar!(w_scalar, mask_bgc_soilc, soilpsi, minpsi, maxpsi, nc, nlevdecomp)

        # O2 scalar
        if use_lch4 && anoxia
            decompb_oscalar_anox!(o_scalar, mask_bgc_soilc, o2stress_unsat, mino2lim,
                                  nc, nlevdecomp)
        else
            decompb_oscalar_one!(o_scalar, nc, nlevdecomp)
        end
    end

    # Normalize Q10 to CENTURY temperature function if requested
    if bgc_state.normalize_q10_to_century_tfunc
        normalization_factor = (catanf(bgc_state.normalization_tref) / catanf_30) /
            (Q10^((bgc_state.normalization_tref - 25.0) / 10.0))
        decompb_tscalar_norm!(t_scalar, mask_bgc_soilc, normalization_factor, nc, nlevdecomp)
    end

    # Depth scalar — fixed e-folding depth
    decompb_depthscalar!(depth_scalar, mask_bgc_soilc, zsoi_vals,
                         decomp_depth_efolding, nc, nlevdecomp)

    # --- Calculate rate constants for all pools ---
    # i_cwd_local is only valid when CWD exists (!use_fates); pass a safe in-bounds
    # index (1) under use_fates since the kernel guards the CWD write on use_fates.
    i_cwd_local = use_fates ? 1 : (i_pas_som + 1)
    decompb_decompk!(decomp_k, mask_bgc_soilc, t_scalar, w_scalar, depth_scalar, o_scalar,
                     spinup_geogterm_l1, spinup_geogterm_l23, spinup_geogterm_cwd,
                     spinup_geogterm_s1, spinup_geogterm_s2, spinup_geogterm_s3,
                     FT_rate(k_l1), FT_rate(k_l2_l3), FT_rate(k_s1), FT_rate(k_s2),
                     FT_rate(k_s3), FT_rate(k_frag),
                     i_met_lit, i_cel_lit, i_lig_lit, i_act_som, i_slo_som, i_pas_som,
                     i_cwd_local, use_fates, nc, nlevdecomp)

    # --- Set pathfrac and rf for the cascade ---
    rf_decomp_cascade       = cf.rf_decomp_cascade_col
    pathfrac_decomp_cascade = cf.pathfrac_decomp_cascade_col

    decompb_pathfrac_rf!(pathfrac_decomp_cascade, rf_decomp_cascade,
                         bgc_state.f_s1s2, bgc_state.f_s1s3,
                         bgc_state.rf_s1s2, bgc_state.rf_s1s3,
                         FT_rate(bgc_state.f_s2s1), FT_rate(bgc_state.f_s2s3),
                         FT_rate(bgc_state.rf_l1s1), FT_rate(bgc_state.rf_l2s1),
                         FT_rate(bgc_state.rf_l3s2), FT_rate(bgc_state.rf_s2s1),
                         FT_rate(bgc_state.rf_s2s3), FT_rate(bgc_state.rf_s3s1),
                         i_l1s1, i_l2s1, i_l3s2, i_s1s2, i_s1s3, i_s2s1, i_s2s3, i_s3s1,
                         nc, nlevdecomp)

    if !use_fates
        i_cwdl2_local = i_cwdl3 - 1  # i_cwdl2 = i_cwdl3 - 1 = 9
        decompb_pathfrac_rf_cwd!(pathfrac_decomp_cascade, rf_decomp_cascade,
                                 FT_rate(bgc_state.cwd_fcel), FT_rate(cwd_flig),
                                 FT_rate(rf_cwdl2), FT_rate(bgc_state.rf_cwdl3),
                                 i_cwdl2_local, i_cwdl3, nc, nlevdecomp)
    end

    return nothing
end
