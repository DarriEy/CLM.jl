# ==========================================================================
# Ported from: src/soilbiogeochem/SoilBiogeochemPotentialMod.F90
# Calculate potential decomp rates and total immobilization demand.
#
# Public functions:
#   soil_bgc_potential!  -- Calculate potential decomposition rates and
#                           total immobilization demand
# ==========================================================================

# ---------------------------------------------------------------------------
# KernelAbstractions kernels for per-element loops in soil_bgc_potential!.
# @kernel/@index/@Const and _launch! are module-wide (infrastructure/kernels.jl).
# Each kernel below has fully independent iterations (no reduction / loop-carried
# dependency); semantics are identical to the scalar loops they replace.
# ---------------------------------------------------------------------------

# Zero a 3D [col x lev x trans] array (one thread per element).
@kernel function _decpot_zero3d_kernel!(out)
    c, j, k = @index(Global, NTuple)
    @inbounds out[c, j, k] = zero(eltype(out))
end
decpot_zero3d!(out) = _launch!(_decpot_zero3d_kernel!, out; ndrange = size(out))

# Zero a 2D [col x lev] array on active (masked) columns over an explicit
# (nc, nlevdecomp) ndrange (the output may have more levels than nlevdecomp).
@kernel function _decpot_zero2d_kernel!(out, @Const(mask))
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        out[c, j] = zero(eltype(out))
    end
end
decpot_zero2d!(out, mask, nc::Int, nlev::Int) =
    _launch!(_decpot_zero2d_kernel!, out, mask; ndrange = (nc, nlev))

# Copy a 2D [col x lev] array on active (masked) columns over an explicit ndrange.
@kernel function _decpot_copy2d_kernel!(out, @Const(mask), @Const(src))
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        out[c, j] = src[c, j]
    end
end
decpot_copy2d!(out, mask, src, nc::Int, nlev::Int) =
    _launch!(_decpot_copy2d_kernel!, out, mask, src; ndrange = (nc, nlev))

# C:N ratios of decomposing pools, per (col, lev, pool); masked per column.
# Floating-C:N pools take decomp_cpools/decomp_npools where N>0 (else unchanged);
# fixed-C:N pools take initial_cn_ratio[l]. The per-pool floating flag and the
# initial-ratio vector are passed as @Const arrays so the branch is in-kernel.
@kernel function _decpot_cn_pools_kernel!(cn_decomp_pools, @Const(mask),
        @Const(floating), @Const(initial_cn_ratio),
        @Const(decomp_cpools_vr), @Const(decomp_npools_vr))
    c, j, l = @index(Global, NTuple)
    @inbounds if mask[c]
        if floating[l]
            if decomp_npools_vr[c, j, l] > zero(eltype(decomp_npools_vr))
                cn_decomp_pools[c, j, l] = decomp_cpools_vr[c, j, l] / decomp_npools_vr[c, j, l]
            end
        else
            cn_decomp_pools[c, j, l] = initial_cn_ratio[l]
        end
    end
end
decpot_cn_pools!(cn_decomp_pools, mask, floating, initial_cn_ratio,
                 decomp_cpools_vr, decomp_npools_vr) =
    _launch!(_decpot_cn_pools_kernel!, cn_decomp_pools, mask, floating,
             initial_cn_ratio, decomp_cpools_vr, decomp_npools_vr;
             ndrange = size(cn_decomp_pools))

# Non-nitrogen-limited potential C and mineral-N fluxes, per (col, lev, trans).
# Each thread writes only its own k-index of the output arrays; reads donor/
# receiver pool slices read-only. No loop-carried dependency across k. Mirrors
# the scalar loop exactly, including the use_mimics messy-eating branch.
@kernel function _decpot_fluxes_kernel!(p_decomp_cpool_loss, pmnf_decomp_cascade,
        p_decomp_cpool_gain, p_decomp_npool_gain, p_decomp_npool_to_din,
        @Const(mask), @Const(cascade_donor_pool), @Const(cascade_receiver_pool),
        @Const(floating), @Const(decomp_cpools_vr), @Const(decomp_k),
        @Const(pathfrac_decomp_cascade), @Const(decomp_npools_vr),
        @Const(cn_decomp_pools), @Const(rf_decomp_cascade),
        @Const(nue_decomp_cascade), use_mimics::Bool)
    c, j, k = @index(Global, NTuple)
    @inbounds if mask[c]
        T = eltype(p_decomp_cpool_loss)
        dpool = cascade_donor_pool[k]
        if decomp_cpools_vr[c, j, dpool] > zero(T) && decomp_k[c, j, dpool] > zero(T)
            p_decomp_cpool_loss[c, j, k] = decomp_cpools_vr[c, j, dpool] *
                decomp_k[c, j, dpool] * pathfrac_decomp_cascade[c, j, k]

            recv_pool = cascade_receiver_pool[k]
            recv_floating = (recv_pool != I_ATM) && floating[recv_pool]

            if !recv_floating
                if recv_pool != I_ATM
                    ratio = zero(T)
                    if decomp_npools_vr[c, j, dpool] > zero(T)
                        ratio = cn_decomp_pools[c, j, recv_pool] /
                                cn_decomp_pools[c, j, dpool]
                    end
                    pmnf_decomp_cascade[c, j, k] = (p_decomp_cpool_loss[c, j, k] *
                        (one(T) - rf_decomp_cascade[c, j, k] - ratio) /
                        cn_decomp_pools[c, j, recv_pool])
                else
                    pmnf_decomp_cascade[c, j, k] = -p_decomp_cpool_loss[c, j, k] /
                        cn_decomp_pools[c, j, dpool]
                end
            else
                pmnf_decomp_cascade[c, j, k] = zero(T)
                if use_mimics
                    decomp_nc_loss_donor =
                        decomp_npools_vr[c, j, dpool] / decomp_cpools_vr[c, j, dpool]
                    p_decomp_npool_loss = p_decomp_cpool_loss[c, j, k] * decomp_nc_loss_donor
                    p_decomp_cpool_gain[c, j, k] =
                        p_decomp_cpool_loss[c, j, k] * (one(T) - rf_decomp_cascade[c, j, k])
                    p_decomp_npool_gain[c, j, k] = p_decomp_npool_loss * nue_decomp_cascade[k]
                    p_decomp_npool_to_din[c, j, k] =
                        p_decomp_npool_loss - p_decomp_npool_gain[c, j, k]
                end
            end
        else
            p_decomp_cpool_gain[c, j, k] = zero(T)
            p_decomp_npool_gain[c, j, k] = zero(T)
            p_decomp_npool_to_din[c, j, k] = zero(T)
        end
    end
end
function decpot_fluxes!(p_decomp_cpool_loss, pmnf_decomp_cascade,
        p_decomp_cpool_gain, p_decomp_npool_gain, p_decomp_npool_to_din,
        mask, cascade_donor_pool, cascade_receiver_pool, floating,
        decomp_cpools_vr, decomp_k, pathfrac_decomp_cascade, decomp_npools_vr,
        cn_decomp_pools, rf_decomp_cascade, nue_decomp_cascade, use_mimics::Bool)
    _launch!(_decpot_fluxes_kernel!, p_decomp_cpool_loss, pmnf_decomp_cascade,
             p_decomp_cpool_gain, p_decomp_npool_gain, p_decomp_npool_to_din,
             mask, cascade_donor_pool, cascade_receiver_pool, floating,
             decomp_cpools_vr, decomp_k, pathfrac_decomp_cascade, decomp_npools_vr,
             cn_decomp_pools, rf_decomp_cascade, nue_decomp_cascade, use_mimics;
             ndrange = size(p_decomp_cpool_loss))
end

# MIMICS receiver-pool gain + pmnf: ONE THREAD PER COLUMN, with internal j and k
# loops accumulating sequentially in-thread (the c_state_update1! cascade pattern —
# each thread owns its column c, so reading/writing [c,j,rp] is race-free and
# byte-identical to the scalar j->c->k loop). Because the receiver pool `rp` is only
# ever `i_cop_mic` or `i_oli_mic`, the two per-pool running sums are held in two
# scalar accumulators (cop/oli), reset per soil level j — avoiding a per-thread
# pool-sized scratch array (not allocatable on Metal). This loop touches < 20 arrays
# so they are passed loose (no device-view bundle). cn_gain[c,j,rp] is written each
# matching k so after the first k-loop the array holds the final running ratio.
@kernel function _decpot_mimics_kernel!(p_decomp_cn_gain, pmnf_decomp_cascade,
        @Const(mask), @Const(cascade_donor_pool), @Const(cascade_receiver_pool),
        @Const(decomp_cpools_vr), @Const(decomp_k),
        @Const(p_decomp_cpool_gain), @Const(p_decomp_npool_gain), @Const(cn_col),
        nlevdecomp::Int, ndecomp_cascade_transitions::Int,
        i_cop_mic::Int, i_oli_mic::Int)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(p_decomp_cn_gain)
        for j in 1:nlevdecomp
            # Sum C & N fluxes from all transitions into the cop & oli microbe pools.
            cop_c = zero(T); cop_n = zero(T)
            oli_c = zero(T); oli_n = zero(T)
            for k in 1:ndecomp_cascade_transitions
                rp = cascade_receiver_pool[k]
                if rp == i_cop_mic || rp == i_oli_mic
                    dp = cascade_donor_pool[k]
                    if decomp_cpools_vr[c, j, dp] > zero(T) && decomp_k[c, j, dp] > zero(T)
                        if rp == i_cop_mic
                            cop_c += p_decomp_cpool_gain[c, j, k]
                            cop_n += p_decomp_npool_gain[c, j, k]
                            if cop_n > zero(T)
                                p_decomp_cn_gain[c, j, rp] = cop_c / cop_n
                            else
                                p_decomp_cn_gain[c, j, rp] = zero(T)
                            end
                        else
                            oli_c += p_decomp_cpool_gain[c, j, k]
                            oli_n += p_decomp_npool_gain[c, j, k]
                            if oli_n > zero(T)
                                p_decomp_cn_gain[c, j, rp] = oli_c / oli_n
                            else
                                p_decomp_cn_gain[c, j, rp] = zero(T)
                            end
                        end
                    end
                end
            end  # first transitions loop

            for k in 1:ndecomp_cascade_transitions
                rp = cascade_receiver_pool[k]
                if rp == i_cop_mic || rp == i_oli_mic
                    dp = cascade_donor_pool[k]
                    if decomp_cpools_vr[c, j, dp] > zero(T) && decomp_k[c, j, dp] > zero(T)
                        # if p_decomp_cn_diff < 0 => N mineralization
                        #                     > 0 => immobilization
                        # "min" turns off immobilization flux
                        p_decomp_cn_diff_ratio = min(zero(T),
                            (p_decomp_cn_gain[c, j, rp] - cn_col[c, rp]) / cn_col[c, rp])
                        pmnf_decomp_cascade[c, j, k] =
                            p_decomp_cn_diff_ratio * p_decomp_npool_gain[c, j, k]
                    end
                end
            end  # second transitions loop
        end
    end
end
function decpot_mimics!(p_decomp_cn_gain, pmnf_decomp_cascade, mask,
        cascade_donor_pool, cascade_receiver_pool, decomp_cpools_vr, decomp_k,
        p_decomp_cpool_gain, p_decomp_npool_gain, cn_col,
        nlevdecomp::Int, ndecomp_cascade_transitions::Int,
        i_cop_mic::Int, i_oli_mic::Int)
    # First positional (p_decomp_cn_gain) sets the launch backend; ndrange is the
    # column count (one thread per column), passed explicitly.
    _launch!(_decpot_mimics_kernel!, p_decomp_cn_gain, pmnf_decomp_cascade, mask,
             cascade_donor_pool, cascade_receiver_pool, decomp_cpools_vr, decomp_k,
             p_decomp_cpool_gain, p_decomp_npool_gain, cn_col,
             nlevdecomp, ndecomp_cascade_transitions, i_cop_mic, i_oli_mic;
             ndrange = length(mask))
    return nothing
end

# Per-column immobilization / gross-mineralization accumulation: ONE THREAD PER
# COLUMN, internal j and k loops summing into the thread's own [c,j] (own-column
# reduction — race-free, ascending-k/j order preserved -> byte-identical). The MIMICS
# DIN term is added per (c,j,k) when use_mimics. immob is the local accumulator the
# scalar code copies into potential_immob_vr; here we write potential_immob_vr
# directly (it is zeroed on active columns just before).
@kernel function _decpot_immob_kernel!(potential_immob_vr, gross_nmin_vr,
        @Const(mask), @Const(pmnf_decomp_cascade), @Const(p_decomp_npool_to_din),
        nlevdecomp::Int, ndecomp_cascade_transitions::Int, use_mimics::Bool)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(potential_immob_vr)
        for j in 1:nlevdecomp
            immob = zero(T)
            gmin  = gross_nmin_vr[c, j]
            for k in 1:ndecomp_cascade_transitions
                if pmnf_decomp_cascade[c, j, k] > zero(T)
                    immob += pmnf_decomp_cascade[c, j, k]
                else
                    gmin -= pmnf_decomp_cascade[c, j, k]
                end
                if use_mimics
                    gmin += p_decomp_npool_to_din[c, j, k]
                end
            end
            potential_immob_vr[c, j] = immob
            gross_nmin_vr[c, j] = gmin
        end
    end
end
function decpot_immob!(potential_immob_vr, gross_nmin_vr, mask,
        pmnf_decomp_cascade, p_decomp_npool_to_din,
        nlevdecomp::Int, ndecomp_cascade_transitions::Int, use_mimics::Bool)
    _launch!(_decpot_immob_kernel!, potential_immob_vr, gross_nmin_vr, mask,
             pmnf_decomp_cascade, p_decomp_npool_to_din,
             nlevdecomp, ndecomp_cascade_transitions, use_mimics;
             ndrange = length(potential_immob_vr) == 0 ? 0 : size(potential_immob_vr, 1))
end

# Per-column potential-HR accumulation for methane: ONE THREAD PER COLUMN, internal
# j and k loops summing rf*p_decomp_cpool_loss into the thread's own phr_vr[c,j]
# (own-column reduction, ascending order -> byte-identical).
@kernel function _decpot_phr_kernel!(phr_vr, @Const(mask),
        @Const(rf_decomp_cascade), @Const(p_decomp_cpool_loss),
        nlevdecomp::Int, ndecomp_cascade_transitions::Int)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(phr_vr)
        for j in 1:nlevdecomp
            s = zero(T)
            for k in 1:ndecomp_cascade_transitions
                s += rf_decomp_cascade[c, j, k] * p_decomp_cpool_loss[c, j, k]
            end
            phr_vr[c, j] = s
        end
    end
end
function decpot_phr!(phr_vr, mask, rf_decomp_cascade, p_decomp_cpool_loss,
        nlevdecomp::Int, ndecomp_cascade_transitions::Int)
    _launch!(_decpot_phr_kernel!, phr_vr, mask, rf_decomp_cascade, p_decomp_cpool_loss,
             nlevdecomp, ndecomp_cascade_transitions;
             ndrange = length(phr_vr) == 0 ? 0 : size(phr_vr, 1))
end

# ---------------------------------------------------------------------------
# DecompPotentialParams -- parameters for potential decomposition
# Ported from params_type in SoilBiogeochemPotentialMod.F90
# ---------------------------------------------------------------------------

"""
    DecompPotentialParams

Parameters for potential decomposition calculations.
Holds denitrification proportion parameter.

Ported from `params_type` in `SoilBiogeochemPotentialMod.F90`.
"""
Base.@kwdef mutable struct DecompPotentialParams
    dnp::Float64 = 0.01  # denitrification proportion
end

"""
    decomp_potential_read_params!(params; dnp)

Read potential decomposition parameters.
Corresponds to `readParams` in `SoilBiogeochemPotentialMod.F90`.
"""
function decomp_potential_read_params!(params::DecompPotentialParams; dnp::Real)
    params.dnp = dnp
    return nothing
end

# ---------------------------------------------------------------------------
# soil_bgc_potential! -- Main potential decomposition subroutine
# Ported from SoilBiogeochemPotential in SoilBiogeochemPotentialMod.F90
#
# Calculates potential decomp rates, C:N ratios of pools, and total
# immobilization demand prior to competition for mineral N.
# ---------------------------------------------------------------------------

"""
    soil_bgc_potential!(cf, cs, nf, ns, st, cascade_con;
        mask_bgc_soilc, bounds, nlevdecomp, ndecomp_pools,
        ndecomp_cascade_transitions,
        cn_decomp_pools, p_decomp_cpool_loss, p_decomp_cn_gain,
        pmnf_decomp_cascade, p_decomp_npool_to_din,
        use_mimics)

Calculate potential decomposition rates and total immobilization demand.

Given decomposition rate coefficients (`decomp_k_col`) and respiration fractions
(`rf_decomp_cascade_col`) previously computed, this routine calculates:
- C:N ratios of decomposing pools
- Potential C loss from each pool for each transition
- Potential mineral N flux (positive = immobilization, negative = mineralization)
- For MIMICS: potential N flux to dissolved inorganic N from messy eating
- Total potential immobilization and gross mineralization
- Potential heterotrophic respiration (for methane calculations)

Corresponds to `SoilBiogeochemPotential` in `SoilBiogeochemPotentialMod.F90`.

# Arguments
- `cf::SoilBiogeochemCarbonFluxData`: carbon flux data (in/out)
- `cs::SoilBiogeochemCarbonStateData`: carbon state data (in)
- `nf::SoilBiogeochemNitrogenFluxData`: nitrogen flux data (out: potential_immob_vr_col, gross_nmin_vr_col)
- `ns::SoilBiogeochemNitrogenStateData`: nitrogen state data (in)
- `st::SoilBiogeochemStateData`: biogeochem state data (in: nue_decomp_cascade_col)
- `cascade_con::DecompCascadeConData`: cascade configuration
- `mask_bgc_soilc::BitVector`: mask for BGC soil columns
- `bounds::UnitRange{Int}`: column bounds
- `nlevdecomp::Int`: number of decomposition levels
- `ndecomp_pools::Int`: number of decomposition pools
- `ndecomp_cascade_transitions::Int`: number of cascade transitions
- `cn_decomp_pools::Array{<:Real,3}`: (out) C:N ratios of pools [col x lev x pool]
- `p_decomp_cpool_loss::Array{<:Real,3}`: (out) potential C loss [col x lev x trans]
- `p_decomp_cn_gain::Array{<:Real,3}`: (out) C:N of flux gained by receiver [col x lev x pool]
- `pmnf_decomp_cascade::Array{<:Real,3}`: (out) potential mineral N flux [col x lev x trans]
- `p_decomp_npool_to_din::Array{<:Real,3}`: (out) potential flux to DIN [col x lev x trans]
- `use_mimics::Bool`: whether MIMICS decomposition is active
- `i_cop_mic::Int`: index of copiotrophic microbe pool (only used when use_mimics=true)
- `i_oli_mic::Int`: index of oligotrophic microbe pool (only used when use_mimics=true)
"""
function soil_bgc_potential!(
        cf::SoilBiogeochemCarbonFluxData,
        cs::SoilBiogeochemCarbonStateData,
        nf::SoilBiogeochemNitrogenFluxData,
        ns::SoilBiogeochemNitrogenStateData,
        st::SoilBiogeochemStateData,
        cascade_con::DecompCascadeConData;
        mask_bgc_soilc::AbstractVector{Bool},
        bounds::UnitRange{Int},
        nlevdecomp::Int,
        ndecomp_pools::Int,
        ndecomp_cascade_transitions::Int,
        cn_decomp_pools::AbstractArray{<:Real,3},
        p_decomp_cpool_loss::AbstractArray{<:Real,3},
        p_decomp_cn_gain::AbstractArray{<:Real,3},
        pmnf_decomp_cascade::AbstractArray{<:Real,3},
        p_decomp_npool_to_din::AbstractArray{<:Real,3},
        use_mimics::Bool=false,
        i_cop_mic::Int=0,
        i_oli_mic::Int=0)

    nc = length(bounds)

    # Unpack cascade configuration
    cascade_donor_pool             = cascade_con.cascade_donor_pool
    cascade_receiver_pool          = cascade_con.cascade_receiver_pool
    floating_cn_ratio_decomp_pools = cascade_con.floating_cn_ratio_decomp_pools
    initial_cn_ratio               = cascade_con.initial_cn_ratio

    # Unpack state
    nue_decomp_cascade = st.nue_decomp_cascade_col

    # Unpack carbon flux/state
    rf_decomp_cascade       = cf.rf_decomp_cascade_col
    pathfrac_decomp_cascade = cf.pathfrac_decomp_cascade_col
    decomp_k                = cf.decomp_k_col
    cn_col                  = cf.cn_col
    phr_vr                  = cf.phr_vr_col

    # Unpack carbon state
    decomp_cpools_vr = cs.decomp_cpools_vr_col

    # Move the per-pool floating-CN flag (a BitVector in cascade_con) onto the state's
    # backend as a plain Bool vector: a BitVector is non-bitstype and cannot be passed
    # to a device kernel. similar(state, Bool, n)+copyto! lands it on the right backend
    # (Array{Bool} on CPU — byte-identical Bool indexing; MtlArray{Bool} on Metal).
    floating_cn_ratio_decomp_pools = let fh = collect(Bool, floating_cn_ratio_decomp_pools)
        d = similar(decomp_cpools_vr, Bool, length(fh))
        copyto!(d, fh)
        d
    end

    # Unpack nitrogen state/flux
    decomp_npools_vr = ns.decomp_npools_vr_col
    potential_immob_vr = nf.potential_immob_vr_col
    gross_nmin_vr      = nf.gross_nmin_vr_col

    # Local arrays for MIMICS pathway (device-resident scratch; allocated like the
    # state arrays so they live on whatever backend the state does).
    FT = eltype(decomp_cpools_vr)

    # Move the per-pool initial C:N ratios (host Vector{FT} in cascade_con) onto the
    # state backend at working precision — a host Float64 vector is non-bitstype for a
    # device kernel; on CPU this is byte-identical.
    initial_cn_ratio = let r = collect(FT, initial_cn_ratio)
        d = similar(decomp_cpools_vr, FT, length(r))
        copyto!(d, r)
        d
    end

    # Cascade donor/receiver pool indices (host Vector{Int}) also feed the flux kernel —
    # move to the state backend preserving Int eltype.
    cascade_donor_pool = let p = collect(Int, cascade_donor_pool)
        d = similar(decomp_cpools_vr, Int, length(p)); copyto!(d, p); d
    end
    cascade_receiver_pool = let p = collect(Int, cascade_receiver_pool)
        d = similar(decomp_cpools_vr, Int, length(p)); copyto!(d, p); d
    end

    p_decomp_cpool_gain = similar(decomp_cpools_vr, FT, nc, nlevdecomp, ndecomp_cascade_transitions)
    p_decomp_npool_gain = similar(decomp_cpools_vr, FT, nc, nlevdecomp, ndecomp_cascade_transitions)
    fill!(p_decomp_cpool_gain, zero(FT))
    fill!(p_decomp_npool_gain, zero(FT))

    # -------------------------------------------------------------------
    # Set initial values for potential C and N fluxes
    # -------------------------------------------------------------------
    decpot_zero3d!(p_decomp_cpool_loss)
    decpot_zero3d!(pmnf_decomp_cascade)

    # -------------------------------------------------------------------
    # Calculate C:N ratios of applicable pools
    # -------------------------------------------------------------------
    decpot_cn_pools!(cn_decomp_pools, mask_bgc_soilc, floating_cn_ratio_decomp_pools,
                     initial_cn_ratio, decomp_cpools_vr, decomp_npools_vr)

    # -------------------------------------------------------------------
    # Calculate the non-nitrogen-limited fluxes
    # These fluxes include the "/ dt" term to put them on a per second
    # basis, since the rate constants have been calculated on a per
    # timestep basis.
    # -------------------------------------------------------------------
    decpot_fluxes!(p_decomp_cpool_loss, pmnf_decomp_cascade,
                   p_decomp_cpool_gain, p_decomp_npool_gain, p_decomp_npool_to_din,
                   mask_bgc_soilc, cascade_donor_pool, cascade_receiver_pool,
                   floating_cn_ratio_decomp_pools, decomp_cpools_vr, decomp_k,
                   pathfrac_decomp_cascade, decomp_npools_vr, cn_decomp_pools,
                   rf_decomp_cascade, nue_decomp_cascade, use_mimics)

    # -------------------------------------------------------------------
    # MIMICS: Calculate cn_gain into microbial biomass and determine
    # immobilization vs. mineralization
    # -------------------------------------------------------------------
    if use_mimics
        decpot_mimics!(p_decomp_cn_gain, pmnf_decomp_cascade, mask_bgc_soilc,
                       cascade_donor_pool, cascade_receiver_pool, decomp_cpools_vr,
                       decomp_k, p_decomp_cpool_gain, p_decomp_npool_gain, cn_col,
                       nlevdecomp, ndecomp_cascade_transitions, i_cop_mic, i_oli_mic)
    end  # use_mimics

    # -------------------------------------------------------------------
    # Sum up all potential immobilization fluxes (positive pmnf)
    # and all mineralization fluxes (negative pmnf)
    # -------------------------------------------------------------------
    decpot_immob!(potential_immob_vr, gross_nmin_vr, mask_bgc_soilc,
                  pmnf_decomp_cascade, p_decomp_npool_to_din,
                  nlevdecomp, ndecomp_cascade_transitions, use_mimics)

    # -------------------------------------------------------------------
    # Add up potential HR for methane calculations
    # -------------------------------------------------------------------
    decpot_phr!(phr_vr, mask_bgc_soilc, rf_decomp_cascade, p_decomp_cpool_loss,
                nlevdecomp, ndecomp_cascade_transitions)

    return nothing
end
