# ==========================================================================
# Ported from: src/soilbiogeochem/SoilBiogeochemCompetitionMod.F90
# Resolve plant/heterotroph competition for mineral N
#
# Public functions:
#   soil_bgc_competition_params_read!  -- Read competition parameters
#   soil_bgc_competition_init!         -- Initialize module-level state
#   soil_bgc_competition!             -- Main competition routine
# ==========================================================================

# ---------------------------------------------------------------------------
# SoilBGCCompetitionParams -- competition parameters read from file
# Ported from params_type in SoilBiogeochemCompetitionMod.F90
# ---------------------------------------------------------------------------

"""
    SoilBGCCompetitionParams

Parameters for soil N competition. Holds bulk denitrification rate and
relative competitiveness factors for plants, decomposers, nitrifiers,
and denitrifiers.

Ported from `params_type` in `SoilBiogeochemCompetitionMod.F90`.
"""
Base.@kwdef mutable struct SoilBGCCompetitionParams
    bdnr              ::Float64 = 0.5   # bulk denitrification rate (1/s)
    compet_plant_no3  ::Float64 = 1.0   # (unitless) relative competitiveness of plants for NO3
    compet_plant_nh4  ::Float64 = 1.0   # (unitless) relative competitiveness of plants for NH4
    compet_decomp_no3 ::Float64 = 1.0   # (unitless) relative competitiveness of immobilizers for NO3
    compet_decomp_nh4 ::Float64 = 1.0   # (unitless) relative competitiveness of immobilizers for NH4
    compet_denit      ::Float64 = 1.0   # (unitless) relative competitiveness of denitrifiers for NO3
    compet_nit        ::Float64 = 1.0   # (unitless) relative competitiveness of nitrifiers for NH4
end

# ---------------------------------------------------------------------------
# SoilBGCCompetitionState -- module-level state
# Holds timestep-scaled bdnr, carbon_only flag, and supplemental N mode
# ---------------------------------------------------------------------------

"""
    SoilBGCCompetitionState

Module-level persistent state for soil N competition.
Holds the timestep-scaled bulk denitrification rate and carbon-only flag.

Ported from module-level variables in `SoilBiogeochemCompetitionMod.F90`.
"""
Base.@kwdef mutable struct SoilBGCCompetitionState
    dt            ::Float64 = 1800.0   # decomp timestep (seconds)
    bdnr          ::Float64 = 0.0      # bulk denitrification rate scaled to timestep
    carbon_only   ::Bool    = false    # if true, supplement N to eliminate N limitation
end

# --- Supplemental nitrogen mode constants ---
const SUPLN_ALL  = "ALL"
const SUPLN_NONE = "NONE"

# ---------------------------------------------------------------------------
# soil_bgc_competition_params_read! -- Read competition parameters
# Ported from readParams in SoilBiogeochemCompetitionMod.F90
# ---------------------------------------------------------------------------

"""
    soil_bgc_competition_params_read!(params; bdnr, compet_plant_no3,
        compet_plant_nh4, compet_decomp_no3, compet_decomp_nh4,
        compet_denit, compet_nit)

Read competition parameters from keyword arguments (replaces NetCDF file reading).
Corresponds to `readParams` in the Fortran source.
"""
function soil_bgc_competition_params_read!(params::SoilBGCCompetitionParams;
                                            bdnr::Real,
                                            compet_plant_no3::Real,
                                            compet_plant_nh4::Real,
                                            compet_decomp_no3::Real,
                                            compet_decomp_nh4::Real,
                                            compet_denit::Real,
                                            compet_nit::Real)
    params.bdnr              = bdnr
    params.compet_plant_no3  = compet_plant_no3
    params.compet_plant_nh4  = compet_plant_nh4
    params.compet_decomp_no3 = compet_decomp_no3
    params.compet_decomp_nh4 = compet_decomp_nh4
    params.compet_denit      = compet_denit
    params.compet_nit        = compet_nit
    return nothing
end

# ---------------------------------------------------------------------------
# soil_bgc_competition_init! -- Initialization
# Ported from SoilBiogeochemCompetitionInit in SoilBiogeochemCompetitionMod.F90
# ---------------------------------------------------------------------------

"""
    soil_bgc_competition_init!(state, params; dt, secspday, suplnitro)

Initialize the competition module: scale bdnr to the decomp timestep and
set the carbon_only flag based on the supplemental nitrogen mode.

Corresponds to `SoilBiogeochemCompetitionInit` in the Fortran source.
"""
function soil_bgc_competition_init!(state::SoilBGCCompetitionState,
                                     params::SoilBGCCompetitionParams;
                                     dt::Real,
                                     secspday::Real=SECSPDAY,
                                     suplnitro::String=SUPLN_NONE)
    state.dt   = dt
    state.bdnr = params.bdnr * (dt / secspday)

    if suplnitro == SUPLN_NONE
        state.carbon_only = false
    elseif suplnitro == SUPLN_ALL
        state.carbon_only = true
    else
        error("Supplemental Nitrogen flag (suplnitro) can only be: " *
              SUPLN_NONE * " or " * SUPLN_ALL *
              ", got: " * suplnitro)
    end

    return nothing
end

# ---------------------------------------------------------------------------
# KernelAbstractions kernels for the non-nitrif/denitrif competition pathway.
# Each kernel covers a (column, decomp-level) loop nest whose iterations are
# fully independent (no accumulation / loop-carried dependency); reductions
# remain as scalar loops. Mask is applied in-kernel (one thread per (c,j)).
# ---------------------------------------------------------------------------

# N uptake profile: nuptake_prof[c,j] = sminn_vr/sminn_tot, else nfixation_prof.
@kernel function _decompc_nuptake_prof_kernel!(nuptake_prof, @Const(mask),
        @Const(sminn_vr), @Const(sminn_tot), @Const(nfixation_prof))
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        T = eltype(nuptake_prof)
        if sminn_tot[c] > zero(T)
            nuptake_prof[c, j] = sminn_vr[c, j] / sminn_tot[c]
        else
            nuptake_prof[c, j] = nfixation_prof[c, j]
        end
    end
end

decompc_nuptake_prof!(nuptake_prof, mask, sminn_vr, sminn_tot, nfixation_prof, nlevdecomp::Int) =
    _launch!(_decompc_nuptake_prof_kernel!, nuptake_prof, mask, sminn_vr, sminn_tot,
             nfixation_prof; ndrange = (length(mask), nlevdecomp))

# Total N demand per level: sum_ndemand_vr = plant_ndemand*nuptake_prof + potential_immob_vr.
@kernel function _decompc_sum_ndemand_kernel!(sum_ndemand_vr, @Const(mask),
        @Const(plant_ndemand), @Const(nuptake_prof), @Const(potential_immob_vr))
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        sum_ndemand_vr[c, j] = plant_ndemand[c] * nuptake_prof[c, j] + potential_immob_vr[c, j]
    end
end

decompc_sum_ndemand!(sum_ndemand_vr, mask, plant_ndemand, nuptake_prof, potential_immob_vr, nlevdecomp::Int) =
    _launch!(_decompc_sum_ndemand_kernel!, sum_ndemand_vr, mask, plant_ndemand,
             nuptake_prof, potential_immob_vr; ndrange = (length(mask), nlevdecomp))

# Resolve plant/decomposer competition at each (column, level). Writes nlimit,
# fpi_vr, actual_immob_vr, sminn_to_plant_vr, supplement_to_sminn_vr; all inputs
# are per-(c,j) or per-column. Fully independent across (c,j).
@kernel function _decompc_resolve_kernel!(actual_immob_vr, @Const(mask), nlimit,
        fpi_vr, sminn_to_plant_vr, supplement_to_sminn_vr, @Const(sum_ndemand_vr),
        @Const(sminn_vr), @Const(potential_immob_vr), @Const(plant_ndemand),
        @Const(nuptake_prof), dt, carbon_only::Bool)
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        T = eltype(actual_immob_vr)
        sd   = sum_ndemand_vr[c, j]
        smn  = sminn_vr[c, j]
        pim  = potential_immob_vr[c, j]
        pdem = plant_ndemand[c] * nuptake_prof[c, j]
        nlimit[c, j] = (smn > sd * dt) ? 0 : 1         # diagnostic flag (hard integer)
        # The supply-vs-N-limited transition (at supply == demand, i.e. frac == 1) carried a
        # d(fpi)/d(sminn) KINK — the hard-branch SUBGRADIENT behind the reverse-AD cross-domain
        # residual. Reformulated via smooth_min on the availability ratio. smooth_min is TYPE-BASED
        # (exact min for Float32/Float64-:auto → GPU-safe + BYTE-IDENTICAL forward; smooth only for
        # Float64-:always and ForwardDiff.Dual), so this is reverse-AD-smooth without changing the
        # default/GPU physics. (The lone non-identical case under :auto is fpi at the N-limited +
        # potential_immob==0 corner, where fpi is a don't-care: actual_immob = fpi·0 = 0 regardless.)
        if carbon_only
            fpi_vr[c, j] = one(T)
            actual_immob_vr[c, j] = pim
            sminn_to_plant_vr[c, j] = pdem
            supplement_to_sminn_vr[c, j] = smooth_max(zero(T), sd - smn / dt)
        else
            frac = sd > zero(T) ? smn / (sd * dt) : T(2)     # availability ratio (≥1 ⇒ unlimited)
            fpi_vr[c, j] = smooth_min(one(T), frac)
            actual_immob_vr[c, j] = fpi_vr[c, j] * pim
            sminn_to_plant_vr[c, j] = smooth_min(pdem, smn / dt - actual_immob_vr[c, j])
        end
    end
end

function decompc_resolve!(actual_immob_vr, mask, nlimit, fpi_vr, sminn_to_plant_vr,
        supplement_to_sminn_vr, sum_ndemand_vr, sminn_vr, potential_immob_vr,
        plant_ndemand, nuptake_prof, dt, carbon_only::Bool, nlevdecomp::Int)
    _launch!(_decompc_resolve_kernel!, actual_immob_vr, mask, nlimit, fpi_vr,
             sminn_to_plant_vr, supplement_to_sminn_vr, sum_ndemand_vr, sminn_vr,
             potential_immob_vr, plant_ndemand, nuptake_prof, dt, carbon_only;
             ndrange = (length(mask), nlevdecomp))
end

# Distribute residual N to plants: per-(c,j) increment using per-column residual
# totals (computed by the preceding scalar reduction). Independent across (c,j).
@kernel function _decompc_distribute_residual_kernel!(sminn_to_plant_vr, @Const(mask),
        @Const(residual_plant_ndemand), @Const(residual_sminn), @Const(residual_sminn_vr),
        @Const(nlimit), dt)
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        T = eltype(sminn_to_plant_vr)
        if residual_plant_ndemand[c] > zero(T) && residual_sminn[c] > zero(T) && nlimit[c, j] == 0
            sminn_to_plant_vr[c, j] += residual_sminn_vr[c, j] *
                min((residual_plant_ndemand[c] * dt) / residual_sminn[c], one(T)) / dt
        end
    end
end

decompc_distribute_residual!(sminn_to_plant_vr, mask, residual_plant_ndemand,
        residual_sminn, residual_sminn_vr, nlimit, dt, nlevdecomp::Int) =
    _launch!(_decompc_distribute_residual_kernel!, sminn_to_plant_vr, mask,
             residual_plant_ndemand, residual_sminn, residual_sminn_vr, nlimit, dt;
             ndrange = (length(mask), nlevdecomp))

# Excess denitrification per (column, level). Independent across (c,j).
@kernel function _decompc_denit_excess_kernel!(sminn_to_denit_excess_vr, @Const(mask),
        @Const(sminn_to_plant_vr), @Const(sminn_to_plant_fun_vr), @Const(actual_immob_vr),
        @Const(sminn_vr), @Const(sum_ndemand_vr), dt, bdnr, local_use_fun::Bool)
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        T = eltype(sminn_to_denit_excess_vr)
        if !local_use_fun
            if (sminn_to_plant_vr[c, j] + actual_immob_vr[c, j]) * dt < sminn_vr[c, j]
                sminn_to_denit_excess_vr[c, j] = max(bdnr * ((sminn_vr[c, j] / dt) - sum_ndemand_vr[c, j]), zero(T))
            else
                sminn_to_denit_excess_vr[c, j] = zero(T)
            end
        else
            if (sminn_to_plant_fun_vr[c, j] + actual_immob_vr[c, j]) * dt < sminn_vr[c, j]
                sminn_to_denit_excess_vr[c, j] = max(bdnr * ((sminn_vr[c, j] / dt) - sum_ndemand_vr[c, j]), zero(T))
            else
                sminn_to_denit_excess_vr[c, j] = zero(T)
            end
        end
    end
end

function decompc_denit_excess!(sminn_to_denit_excess_vr, mask, sminn_to_plant_vr,
        sminn_to_plant_fun_vr, actual_immob_vr, sminn_vr, sum_ndemand_vr, dt, bdnr,
        local_use_fun::Bool, nlevdecomp::Int)
    _launch!(_decompc_denit_excess_kernel!, sminn_to_denit_excess_vr, mask,
             sminn_to_plant_vr, sminn_to_plant_fun_vr, actual_immob_vr, sminn_vr,
             sum_ndemand_vr, dt, bdnr, local_use_fun; ndrange = (length(mask), nlevdecomp))
end

# ---------------------------------------------------------------------------
# Additional per-COLUMN kernels (one thread per column with internal j-loops).
# These cover the remaining HOST reductions/per-column scalars so the whole
# function runs on the device. Each thread owns its column c, accumulating into
# its own [c] (ascending-j order preserved) — race-free + byte-identical.
# ---------------------------------------------------------------------------

# --- Non-nitrif: sum total mineral N over levels (per-column reduction over j) ---
@kernel function _decompc_sminn_tot_kernel!(@Const(mask), sminn_tot,
        @Const(sminn_vr), @Const(dzsoi_decomp), nlevdecomp::Int)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(sminn_tot)
        s = zero(T)
        for j in 1:nlevdecomp
            s += sminn_vr[c, j] * dzsoi_decomp[j]
        end
        sminn_tot[c] = s
    end
end

# --- Non-nitrif: first pass — sum N fluxes to plant + use_fun clamp ---
@kernel function _decompc_sum_to_plant1_kernel!(@Const(mask), sminn_to_plant,
        sminn_to_plant_fun_vr, @Const(sminn_to_plant_vr), @Const(dzsoi_decomp),
        nlevdecomp::Int, local_use_fun::Bool)
    c = @index(Global)
    @inbounds if mask[c]
        # Accumulate DIRECTLY into the own-column field (thread c owns [c]) — matches
        # the host loop's association exactly (byte-identical even if [c] is nonzero
        # on entry; a local-sum + single add would reassociate).
        for j in 1:nlevdecomp
            sminn_to_plant[c] += sminn_to_plant_vr[c, j] * dzsoi_decomp[j]
            if local_use_fun
                if sminn_to_plant_fun_vr[c, j] > sminn_to_plant_vr[c, j]
                    sminn_to_plant_fun_vr[c, j] = sminn_to_plant_vr[c, j]
                end
            end
        end
    end
end

# --- Non-nitrif: residual init + residual_sminn reduction over levels ---
@kernel function _decompc_residual_kernel!(@Const(mask), residual_sminn,
        residual_plant_ndemand, residual_sminn_vr, @Const(plant_ndemand),
        @Const(sminn_to_plant), @Const(sminn_vr), @Const(actual_immob_vr),
        @Const(sminn_to_plant_vr), @Const(nlimit), @Const(dzsoi_decomp),
        nlevdecomp::Int, dt)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(residual_sminn)
        residual_sminn[c] = zero(T)
        residual_plant_ndemand[c] = plant_ndemand[c] - sminn_to_plant[c]
        for j in 1:nlevdecomp
            if residual_plant_ndemand[c] > zero(T)
                if nlimit[c, j] == 0
                    residual_sminn_vr[c, j] = max(sminn_vr[c, j] -
                        (actual_immob_vr[c, j] + sminn_to_plant_vr[c, j]) * dt, zero(T))
                    residual_sminn[c] += residual_sminn_vr[c, j] * dzsoi_decomp[j]
                else
                    residual_sminn_vr[c, j] = zero(T)
                end
            end
        end
    end
end

# --- Non-nitrif: re-sum N fluxes to plant + recompute sum_ndemand_vr ---
@kernel function _decompc_resum_to_plant_kernel!(@Const(mask), sminn_to_plant,
        sminn_to_plant_new, sum_ndemand_vr, @Const(sminn_to_plant_vr),
        @Const(sminn_to_plant_fun_vr), @Const(potential_immob_vr),
        @Const(dzsoi_decomp), nlevdecomp::Int, local_use_fun::Bool)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(sminn_to_plant)
        sminn_to_plant[c] = zero(T)
        for j in 1:nlevdecomp
            sminn_to_plant[c] += sminn_to_plant_vr[c, j] * dzsoi_decomp[j]
            if !local_use_fun
                sum_ndemand_vr[c, j] = potential_immob_vr[c, j] + sminn_to_plant_vr[c, j]
            else
                sminn_to_plant_new[c] += sminn_to_plant_fun_vr[c, j] * dzsoi_decomp[j]
                sum_ndemand_vr[c, j] = potential_immob_vr[c, j] + sminn_to_plant_fun_vr[c, j]
            end
        end
    end
end

# --- Non-nitrif: sum N to immobilization + final fpg/fpi per-column scalars ---
@kernel function _decompc_immob_fpgfpi_kernel!(@Const(mask), actual_immob,
        potential_immob, fpg, fpi, @Const(actual_immob_vr), @Const(potential_immob_vr),
        @Const(plant_ndemand), @Const(sminn_to_plant), @Const(sminn_to_plant_new),
        @Const(dzsoi_decomp), nlevdecomp::Int, local_use_fun::Bool)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(actual_immob)
        for j in 1:nlevdecomp
            actual_immob[c] += actual_immob_vr[c, j] * dzsoi_decomp[j]
            potential_immob[c] += potential_immob_vr[c, j] * dzsoi_decomp[j]
        end
        # fraction of potential growth achieved with available N
        if plant_ndemand[c] > zero(T)
            if !local_use_fun
                fpg[c] = sminn_to_plant[c] / plant_ndemand[c]
            else
                fpg[c] = sminn_to_plant_new[c] / plant_ndemand[c]
            end
        else
            fpg[c] = one(T)
        end
        # fraction of immobilization realized
        if potential_immob[c] > zero(T)
            fpi[c] = actual_immob[c] / potential_immob[c]
        else
            fpi[c] = one(T)
        end
    end
end

# ---------------------------------------------------------------------------
# NITRIF_DENITRIF per-column kernels (one thread per column, internal j/k loops).
# ---------------------------------------------------------------------------

# --- Nitrif: sum total mineral N (no3 + nh4) over levels (per-column reduction) ---
@kernel function _decompc_nd_sminn_tot_kernel!(@Const(mask), sminn_tot,
        @Const(smin_no3_vr), @Const(smin_nh4_vr), @Const(dzsoi_decomp), nlevdecomp::Int)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(sminn_tot)
        s = zero(T)
        for j in 1:nlevdecomp
            s += (smin_no3_vr[c, j] + smin_nh4_vr[c, j]) * dzsoi_decomp[j]
        end
        sminn_tot[c] = s
    end
end

# --- Nitrif: N uptake profile (per-(c,j) own write) ---
@kernel function _decompc_nd_nuptake_prof_kernel!(nuptake_prof, @Const(mask),
        @Const(sminn_vr), @Const(sminn_tot), @Const(nfixation_prof))
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        T = eltype(nuptake_prof)
        if sminn_tot[c] > zero(T)
            nuptake_prof[c, j] = sminn_vr[c, j] / sminn_tot[c]
        else
            nuptake_prof[c, j] = nfixation_prof[c, j]
        end
    end
end

# --- Nitrif main column/vertical loop. Device-view bundles group the ~30 N-flux
#     arrays + the loose state arrays to stay under Metal's arg limit. The internal
#     j-loop runs sequentially in-thread; every write is to own [c,j]. ---
# `M` = float column x level matrices; `I` = integer nlimit matrices (separate
# type param so the Int eltype is preserved when adapting to the device).
Base.@kwdef struct _DECOMPND{M,I}
    # N flux arrays the main loop reads/writes (all 2D column x level)
    potential_immob_vr::M
    actual_immob_vr::M; actual_immob_no3_vr::M; actual_immob_nh4_vr::M
    sminn_to_plant_vr::M; smin_no3_to_plant_vr::M; smin_nh4_to_plant_vr::M
    pot_f_nit_vr::M; pot_f_denit_vr::M; f_nit_vr::M; f_denit_vr::M
    n2_n2o_ratio_denit_vr::M; f_n2o_denit_vr::M; f_n2o_nit_vr::M
    supplement_to_sminn_vr::M; fpi_vr::M
    # N state arrays
    sminn_vr::M; smin_nh4_vr::M; smin_no3_vr::M
    # working / local arrays
    fpi_no3_vr_local::M; fpi_nh4_vr_local::M
    sum_nh4_demand::M; sum_nh4_demand_scaled::M
    sum_no3_demand::M; sum_no3_demand_scaled::M
    nlimit_no3::I; nlimit_nh4::I
end
Adapt.@adapt_structure _DECOMPND

@kernel function _decompc_nd_main_kernel!(@Const(mask), nd, @Const(plant_ndemand),
        @Const(nuptake_prof), dt,
        compet_plant_no3, compet_plant_nh4, compet_decomp_no3, compet_decomp_nh4,
        compet_denit, compet_nit, nitrif_n2o_loss_frac,
        nlevdecomp::Int, local_use_fun::Bool, mimics_decomp::Bool, carbon_only::Bool)
    c = @index(Global)
    @inbounds if mask[c]
        T = typeof(dt)
        for j in 1:nlevdecomp
            # --- First compete for NH4 ---
            nd.sum_nh4_demand[c, j] = plant_ndemand[c] * nuptake_prof[c, j] +
                nd.potential_immob_vr[c, j] + nd.pot_f_nit_vr[c, j]
            nd.sum_nh4_demand_scaled[c, j] = plant_ndemand[c] * nuptake_prof[c, j] * compet_plant_nh4 +
                nd.potential_immob_vr[c, j] * compet_decomp_nh4 + nd.pot_f_nit_vr[c, j] * compet_nit

            if nd.sum_nh4_demand[c, j] * dt < nd.smin_nh4_vr[c, j]
                # NH4 not limiting
                nd.nlimit_nh4[c, j] = 0
                nd.fpi_nh4_vr_local[c, j] = one(T)
                nd.actual_immob_nh4_vr[c, j] = nd.potential_immob_vr[c, j]
                nd.f_nit_vr[c, j] = nd.pot_f_nit_vr[c, j]

                if !local_use_fun
                    nd.smin_nh4_to_plant_vr[c, j] = plant_ndemand[c] * nuptake_prof[c, j]
                else
                    nd.smin_nh4_to_plant_vr[c, j] = nd.smin_nh4_vr[c, j] / dt -
                        nd.actual_immob_nh4_vr[c, j] - nd.f_nit_vr[c, j]
                end
            else
                # NH4 limited: competition
                nd.nlimit_nh4[c, j] = 1
                if nd.sum_nh4_demand[c, j] > zero(T)
                    nd.actual_immob_nh4_vr[c, j] = min(
                        (nd.smin_nh4_vr[c, j] / dt) * (nd.potential_immob_vr[c, j] *
                            compet_decomp_nh4 / nd.sum_nh4_demand_scaled[c, j]),
                        nd.potential_immob_vr[c, j])

                    nd.f_nit_vr[c, j] = min(
                        (nd.smin_nh4_vr[c, j] / dt) * (nd.pot_f_nit_vr[c, j] * compet_nit /
                            nd.sum_nh4_demand_scaled[c, j]),
                        nd.pot_f_nit_vr[c, j])

                    if !local_use_fun
                        nd.smin_nh4_to_plant_vr[c, j] = min(
                            (nd.smin_nh4_vr[c, j] / dt) * (plant_ndemand[c] *
                                nuptake_prof[c, j] * compet_plant_nh4 / nd.sum_nh4_demand_scaled[c, j]),
                            plant_ndemand[c] * nuptake_prof[c, j])
                    else
                        nd.smin_nh4_to_plant_vr[c, j] = nd.smin_nh4_vr[c, j] / dt -
                            nd.actual_immob_nh4_vr[c, j] - nd.f_nit_vr[c, j]
                    end
                else
                    nd.actual_immob_nh4_vr[c, j] = zero(T)
                    nd.smin_nh4_to_plant_vr[c, j] = zero(T)
                    nd.f_nit_vr[c, j] = zero(T)
                end

                if nd.potential_immob_vr[c, j] > zero(T)
                    nd.fpi_nh4_vr_local[c, j] = nd.actual_immob_nh4_vr[c, j] / nd.potential_immob_vr[c, j]
                else
                    nd.fpi_nh4_vr_local[c, j] = zero(T)
                end
            end

            if mimics_decomp
                nd.fpi_nh4_vr_local[c, j] = one(T)
                nd.actual_immob_nh4_vr[c, j] = nd.potential_immob_vr[c, j]
            end

            # --- Then compete for NO3 ---
            if !local_use_fun
                nd.sum_no3_demand[c, j] = (plant_ndemand[c] * nuptake_prof[c, j] -
                    nd.smin_nh4_to_plant_vr[c, j]) +
                    (nd.potential_immob_vr[c, j] - nd.actual_immob_nh4_vr[c, j]) +
                    nd.pot_f_denit_vr[c, j]
                nd.sum_no3_demand_scaled[c, j] = (plant_ndemand[c] * nuptake_prof[c, j] -
                    nd.smin_nh4_to_plant_vr[c, j]) * compet_plant_no3 +
                    (nd.potential_immob_vr[c, j] - nd.actual_immob_nh4_vr[c, j]) * compet_decomp_no3 +
                    nd.pot_f_denit_vr[c, j] * compet_denit
            else
                nd.sum_no3_demand[c, j] = plant_ndemand[c] * nuptake_prof[c, j] +
                    (nd.potential_immob_vr[c, j] - nd.actual_immob_nh4_vr[c, j]) +
                    nd.pot_f_denit_vr[c, j]
                nd.sum_no3_demand_scaled[c, j] = (plant_ndemand[c] * nuptake_prof[c, j]) * compet_plant_no3 +
                    (nd.potential_immob_vr[c, j] - nd.actual_immob_nh4_vr[c, j]) * compet_decomp_no3 +
                    nd.pot_f_denit_vr[c, j] * compet_denit
            end

            if nd.sum_no3_demand[c, j] * dt < nd.smin_no3_vr[c, j]
                # NO3 not limiting
                nd.nlimit_no3[c, j] = 0
                nd.fpi_no3_vr_local[c, j] = one(T) - nd.fpi_nh4_vr_local[c, j]
                nd.actual_immob_no3_vr[c, j] = nd.potential_immob_vr[c, j] - nd.actual_immob_nh4_vr[c, j]
                nd.f_denit_vr[c, j] = nd.pot_f_denit_vr[c, j]

                if !local_use_fun
                    nd.smin_no3_to_plant_vr[c, j] = plant_ndemand[c] * nuptake_prof[c, j] -
                        nd.smin_nh4_to_plant_vr[c, j]
                else
                    nd.smin_no3_to_plant_vr[c, j] = nd.smin_no3_vr[c, j] / dt -
                        nd.actual_immob_no3_vr[c, j] - nd.f_denit_vr[c, j]
                end
            else
                # NO3 limited: competition
                nd.nlimit_no3[c, j] = 1

                if nd.sum_no3_demand[c, j] > zero(T)
                    if !local_use_fun
                        nd.actual_immob_no3_vr[c, j] = min(
                            (nd.smin_no3_vr[c, j] / dt) * ((nd.potential_immob_vr[c, j] -
                                nd.actual_immob_nh4_vr[c, j]) * compet_decomp_no3 /
                                nd.sum_no3_demand_scaled[c, j]),
                            nd.potential_immob_vr[c, j] - nd.actual_immob_nh4_vr[c, j])

                        nd.smin_no3_to_plant_vr[c, j] = min(
                            (nd.smin_no3_vr[c, j] / dt) * ((plant_ndemand[c] *
                                nuptake_prof[c, j] - nd.smin_nh4_to_plant_vr[c, j]) * compet_plant_no3 /
                                nd.sum_no3_demand_scaled[c, j]),
                            plant_ndemand[c] * nuptake_prof[c, j] - nd.smin_nh4_to_plant_vr[c, j])

                        nd.f_denit_vr[c, j] = min(
                            (nd.smin_no3_vr[c, j] / dt) * (nd.pot_f_denit_vr[c, j] * compet_denit /
                                nd.sum_no3_demand_scaled[c, j]),
                            nd.pot_f_denit_vr[c, j])
                    else
                        nd.actual_immob_no3_vr[c, j] = min(
                            (nd.smin_no3_vr[c, j] / dt) * ((nd.potential_immob_vr[c, j] -
                                nd.actual_immob_nh4_vr[c, j]) * compet_decomp_no3 /
                                nd.sum_no3_demand_scaled[c, j]),
                            nd.potential_immob_vr[c, j] - nd.actual_immob_nh4_vr[c, j])

                        nd.f_denit_vr[c, j] = min(
                            (nd.smin_no3_vr[c, j] / dt) * (nd.pot_f_denit_vr[c, j] * compet_denit /
                                nd.sum_no3_demand_scaled[c, j]),
                            nd.pot_f_denit_vr[c, j])

                        nd.smin_no3_to_plant_vr[c, j] = (nd.smin_no3_vr[c, j] / dt) -
                            nd.actual_immob_no3_vr[c, j] - nd.f_denit_vr[c, j]
                    end
                else
                    nd.actual_immob_no3_vr[c, j] = zero(T)
                    nd.smin_no3_to_plant_vr[c, j] = zero(T)
                    nd.f_denit_vr[c, j] = zero(T)
                end

                if nd.potential_immob_vr[c, j] > zero(T)
                    nd.fpi_no3_vr_local[c, j] = nd.actual_immob_no3_vr[c, j] / nd.potential_immob_vr[c, j]
                else
                    nd.fpi_no3_vr_local[c, j] = zero(T)
                end
            end

            if mimics_decomp
                nd.fpi_no3_vr_local[c, j] = one(T) - nd.fpi_nh4_vr_local[c, j]
                nd.actual_immob_no3_vr[c, j] = nd.potential_immob_vr[c, j] - nd.actual_immob_nh4_vr[c, j]
            end

            # N2O emissions
            nd.f_n2o_nit_vr[c, j] = nd.f_nit_vr[c, j] * nitrif_n2o_loss_frac
            nd.f_n2o_denit_vr[c, j] = nd.f_denit_vr[c, j] / (one(T) + nd.n2_n2o_ratio_denit_vr[c, j])

            # Carbon-only supplement
            if carbon_only
                if nd.fpi_no3_vr_local[c, j] + nd.fpi_nh4_vr_local[c, j] < one(T)
                    nd.fpi_nh4_vr_local[c, j] = one(T) - nd.fpi_no3_vr_local[c, j]
                    nd.supplement_to_sminn_vr[c, j] = (nd.potential_immob_vr[c, j] -
                        nd.actual_immob_no3_vr[c, j]) - nd.actual_immob_nh4_vr[c, j]
                    nd.actual_immob_nh4_vr[c, j] = nd.potential_immob_vr[c, j] - nd.actual_immob_no3_vr[c, j]
                end
                if nd.smin_no3_to_plant_vr[c, j] + nd.smin_nh4_to_plant_vr[c, j] < plant_ndemand[c] * nuptake_prof[c, j]
                    nd.supplement_to_sminn_vr[c, j] += (plant_ndemand[c] * nuptake_prof[c, j] -
                        nd.smin_no3_to_plant_vr[c, j]) - nd.smin_nh4_to_plant_vr[c, j]
                    nd.smin_nh4_to_plant_vr[c, j] = plant_ndemand[c] * nuptake_prof[c, j] -
                        nd.smin_no3_to_plant_vr[c, j]
                end
                nd.sminn_to_plant_vr[c, j] = nd.smin_no3_to_plant_vr[c, j] + nd.smin_nh4_to_plant_vr[c, j]
            end

            # sum up no3 and nh4 fluxes
            nd.fpi_vr[c, j] = nd.fpi_no3_vr_local[c, j] + nd.fpi_nh4_vr_local[c, j]
            nd.sminn_to_plant_vr[c, j] = nd.smin_no3_to_plant_vr[c, j] + nd.smin_nh4_to_plant_vr[c, j]
            nd.actual_immob_vr[c, j] = nd.actual_immob_no3_vr[c, j] + nd.actual_immob_nh4_vr[c, j]
        end
    end
end

# --- Nitrif: zero sminn_to_plant + (non-fun) sum it over levels ---
@kernel function _decompc_nd_postsum_kernel!(@Const(mask), sminn_to_plant,
        @Const(sminn_to_plant_vr), @Const(dzsoi_decomp), nlevdecomp::Int, local_use_fun::Bool)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(sminn_to_plant)
        sminn_to_plant[c] = zero(T)
        if !local_use_fun
            for j in 1:nlevdecomp
                sminn_to_plant[c] += sminn_to_plant_vr[c, j] * dzsoi_decomp[j]
            end
        end
    end
end

# --- Nitrif: MIMICS c_overflow_vr (per-column, internal j/k loops) ---
@kernel function _decompc_nd_mimics_overflow_kernel!(@Const(mask), c_overflow_vr,
        sum_ndemand_vr, @Const(sum_no3_demand_scaled), @Const(sum_nh4_demand_scaled),
        @Const(pmnf_decomp_cascade), @Const(p_decomp_cn_gain), @Const(sminn_vr),
        @Const(cascade_receiver_pool), nlevdecomp::Int, ndecomp_cascade_transitions::Int,
        i_cop_mic::Int, i_oli_mic::Int, dt)
    c = @index(Global)
    @inbounds if mask[c]
        T = typeof(dt)
        for j in 1:nlevdecomp
            for k in 1:ndecomp_cascade_transitions
                if cascade_receiver_pool[k] == i_cop_mic || cascade_receiver_pool[k] == i_oli_mic
                    sum_ndemand_vr[c, j] = sum_no3_demand_scaled[c, j] + sum_nh4_demand_scaled[c, j]
                    if pmnf_decomp_cascade[c, j, k] > zero(T) && sum_ndemand_vr[c, j] > zero(T)
                        amnf_immob_vr = (sminn_vr[c, j] / dt) *
                            (pmnf_decomp_cascade[c, j, k] / sum_ndemand_vr[c, j])
                        n_deficit_vr = pmnf_decomp_cascade[c, j, k] - amnf_immob_vr
                        c_overflow_vr[c, j, k] = n_deficit_vr *
                            p_decomp_cn_gain[c, j, cascade_receiver_pool[k]]
                    else
                        c_overflow_vr[c, j, k] = zero(T)
                    end
                else
                    c_overflow_vr[c, j, k] = zero(T)
                end
            end
        end
    end
end

# --- Nitrif (non-fun): second-pass NH4 residual distribution (per-column, internal j) ---
@kernel function _decompc_nd_residual_nh4_kernel!(@Const(mask), residual_plant_ndemand,
        residual_smin_nh4, residual_smin_nh4_vr, smin_nh4_to_plant_vr,
        @Const(plant_ndemand), @Const(sminn_to_plant), @Const(smin_nh4_vr),
        @Const(actual_immob_nh4_vr), @Const(f_nit_vr), @Const(nlimit_nh4),
        @Const(dzsoi_decomp), nlevdecomp::Int, dt)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(residual_smin_nh4)
        residual_plant_ndemand[c] = plant_ndemand[c] - sminn_to_plant[c]
        residual_smin_nh4[c] = zero(T)
        for j in 1:nlevdecomp
            if residual_plant_ndemand[c] > zero(T)
                if nlimit_nh4[c, j] == 0
                    residual_smin_nh4_vr[c, j] = max(smin_nh4_vr[c, j] -
                        (actual_immob_nh4_vr[c, j] + smin_nh4_to_plant_vr[c, j] +
                         f_nit_vr[c, j]) * dt, zero(T))
                    residual_smin_nh4[c] += residual_smin_nh4_vr[c, j] * dzsoi_decomp[j]
                else
                    residual_smin_nh4_vr[c, j] = zero(T)
                end

                if residual_smin_nh4[c] > zero(T) && nlimit_nh4[c, j] == 0
                    smin_nh4_to_plant_vr[c, j] += residual_smin_nh4_vr[c, j] *
                        min((residual_plant_ndemand[c] * dt) / residual_smin_nh4[c], one(T)) / dt
                end
            end
        end
    end
end

# --- Nitrif (non-fun): re-sum sminn_to_plant_vr + sminn_to_plant after NH4 pass ---
@kernel function _decompc_nd_resum_kernel!(@Const(mask), sminn_to_plant,
        sminn_to_plant_vr, @Const(smin_nh4_to_plant_vr), @Const(smin_no3_to_plant_vr),
        @Const(dzsoi_decomp), nlevdecomp::Int)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(sminn_to_plant)
        sminn_to_plant[c] = zero(T)
        for j in 1:nlevdecomp
            sminn_to_plant_vr[c, j] = smin_nh4_to_plant_vr[c, j] + smin_no3_to_plant_vr[c, j]
            sminn_to_plant[c] += sminn_to_plant_vr[c, j] * dzsoi_decomp[j]
        end
    end
end

# --- Nitrif (non-fun): second-pass NO3 residual distribution (per-column, internal j) ---
@kernel function _decompc_nd_residual_no3_kernel!(@Const(mask), residual_plant_ndemand,
        residual_smin_no3, residual_smin_no3_vr, smin_no3_to_plant_vr,
        @Const(plant_ndemand), @Const(sminn_to_plant), @Const(smin_no3_vr),
        @Const(actual_immob_no3_vr), @Const(f_denit_vr), @Const(nlimit_no3),
        @Const(dzsoi_decomp), nlevdecomp::Int, dt)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(residual_smin_no3)
        residual_plant_ndemand[c] = plant_ndemand[c] - sminn_to_plant[c]
        residual_smin_no3[c] = zero(T)
        for j in 1:nlevdecomp
            if residual_plant_ndemand[c] > zero(T)
                if nlimit_no3[c, j] == 0
                    residual_smin_no3_vr[c, j] = max(smin_no3_vr[c, j] -
                        (actual_immob_no3_vr[c, j] + smin_no3_to_plant_vr[c, j] +
                         f_denit_vr[c, j]) * dt, zero(T))
                    residual_smin_no3[c] += residual_smin_no3_vr[c, j] * dzsoi_decomp[j]
                else
                    residual_smin_no3_vr[c, j] = zero(T)
                end

                if residual_smin_no3[c] > zero(T) && nlimit_no3[c, j] == 0
                    smin_no3_to_plant_vr[c, j] += residual_smin_no3_vr[c, j] *
                        min((residual_plant_ndemand[c] * dt) / residual_smin_no3[c], one(T)) / dt
                end
            end
        end
    end
end

# --- Nitrif (use_fun): re-sum sminn_to_plant_vr + add FUN fluxes ---
@kernel function _decompc_nd_fun_resum_kernel!(@Const(mask), sminn_to_plant,
        sminn_to_plant_new, sminn_to_plant_vr, @Const(smin_nh4_to_plant_vr),
        @Const(smin_no3_to_plant_vr), @Const(sminn_to_plant_fun_no3_vr),
        @Const(sminn_to_plant_fun_nh4_vr), @Const(dzsoi_decomp), nlevdecomp::Int)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(sminn_to_plant)
        sminn_to_plant[c] = zero(T)
        for j in 1:nlevdecomp
            sminn_to_plant_vr[c, j] = smin_nh4_to_plant_vr[c, j] + smin_no3_to_plant_vr[c, j]
            sminn_to_plant[c] += sminn_to_plant_vr[c, j] * dzsoi_decomp[j]
        end
        for j in 1:nlevdecomp
            sminn_to_plant_new[c] += (sminn_to_plant_fun_no3_vr[c, j] +
                sminn_to_plant_fun_nh4_vr[c, j]) * dzsoi_decomp[j]
        end
    end
end

# --- Nitrif: sum immobilization + final fpg/fpi per-column scalars ---
@kernel function _decompc_nd_immob_fpgfpi_kernel!(@Const(mask), actual_immob,
        potential_immob, fpg, fpi, @Const(actual_immob_vr), @Const(potential_immob_vr),
        @Const(plant_ndemand), @Const(sminn_to_plant), @Const(dzsoi_decomp),
        nlevdecomp::Int, local_use_fun::Bool)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(actual_immob)
        actual_immob[c] = zero(T)
        potential_immob[c] = zero(T)
        for j in 1:nlevdecomp
            actual_immob[c] += actual_immob_vr[c, j] * dzsoi_decomp[j]
            potential_immob[c] += potential_immob_vr[c, j] * dzsoi_decomp[j]
        end
        if !local_use_fun
            if plant_ndemand[c] > zero(T)
                fpg[c] = sminn_to_plant[c] / plant_ndemand[c]
            else
                fpg[c] = one(T)
            end
        end
        if potential_immob[c] > zero(T)
            fpi[c] = actual_immob[c] / potential_immob[c]
        else
            fpi[c] = one(T)
        end
    end
end

# ---------------------------------------------------------------------------
# soil_bgc_competition! -- Main competition routine
# Ported from SoilBiogeochemCompetition in SoilBiogeochemCompetitionMod.F90
#
# Resolves plant/heterotroph (and optionally nitrifier/denitrifier)
# competition for mineral N.
#
# GPU kernelization: the whole per-timestep compute now runs as
# KernelAbstractions kernels (whole-fn on Metal). Every per-column reduction over
# levels (sminn_tot, sminn_to_plant, residual_*, actual/potential_immob, fpg/fpi)
# is a per-column kernel with an internal sequential j-loop — each thread owns its
# column c, so it is race-free + byte-identical to the host loop. The big nitrif
# main column/vertical loop is a per-column kernel with an internal j-loop, its
# ~27 N-flux/state/working arrays grouped into the `_DECOMPND` device-view bundle
# to stay under Metal's ~31-arg limit. Float64 literals are eltype-converted
# (`zero(T)`/`one(T)`) so the kernels carry no Float64 on a Float32-only backend;
# on Float64 this is byte-identical. The five `_decompc_*` (c,j) kernels are
# preserved.
# ---------------------------------------------------------------------------

"""
    soil_bgc_competition!(st, nf, cf, ns, state, params; ...)

Resolve plant/heterotroph/nitrifier/denitrifier competition for mineral N.

# Arguments
- `st::SoilBiogeochemStateData` -- soil BGC state (fpg, fpi, fpi_vr, nfixation_prof, plant_ndemand)
- `nf::SoilBiogeochemNitrogenFluxData` -- N fluxes (immob, plant uptake, nitrif/denitrif)
- `cf::SoilBiogeochemCarbonFluxData` -- C fluxes (c_overflow_vr for MIMICS)
- `ns::SoilBiogeochemNitrogenStateData` -- N state (sminn_vr, smin_nh4_vr, smin_no3_vr)
- `state::SoilBGCCompetitionState` -- module-level state (dt, bdnr)
- `params::SoilBGCCompetitionParams` -- competition parameters (competitiveness factors)
- `mask_bgc_soilc::AbstractVector{Bool}` -- mask for BGC soil columns
- `bounds::UnitRange{Int}` -- column bounds
- `nlevdecomp::Int` -- number of decomposition levels
- `ndecomp_cascade_transitions::Int` -- number of decomp cascade transitions
- `dzsoi_decomp::AbstractVector{<:Real}` -- decomposition level thicknesses
- `pmnf_decomp_cascade::AbstractArray{<:Real,3}` -- potential mineral N flux from decomp
- `p_decomp_cn_gain::AbstractArray{<:Real,3}` -- C:N ratio of flux gained by receiver pool
- `cascade_receiver_pool::AbstractVector{<:Integer}` -- receiver pool index for each transition
- `use_nitrif_denitrif::Bool` -- use nitrification/denitrification model
- `use_fun::Bool` -- use FUN model
- `carbon_only::Bool` -- carbon-only mode (supplement N)
- `mimics_decomp::Bool` -- using MIMICS decomposition
- `i_cop_mic::Int` -- copiotrophic microbe pool index (MIMICS)
- `i_oli_mic::Int` -- oligotrophic microbe pool index (MIMICS)
- `nitrif_n2o_loss_frac::Real` -- fraction of nitrification lost as N2O

Corresponds to `SoilBiogeochemCompetition` in the Fortran source.
"""
function soil_bgc_competition!(
        st::SoilBiogeochemStateData,
        nf::SoilBiogeochemNitrogenFluxData,
        cf::SoilBiogeochemCarbonFluxData,
        ns::SoilBiogeochemNitrogenStateData,
        state::SoilBGCCompetitionState,
        params::SoilBGCCompetitionParams;
        mask_bgc_soilc::AbstractVector{Bool},
        bounds::UnitRange{Int},
        nlevdecomp::Int,
        ndecomp_cascade_transitions::Int,
        dzsoi_decomp::AbstractVector{<:Real},
        pmnf_decomp_cascade::AbstractArray{<:Real,3},
        p_decomp_cn_gain::AbstractArray{<:Real,3},
        cascade_receiver_pool::AbstractVector{<:Integer},
        use_nitrif_denitrif::Bool=false,
        use_fun::Bool=false,
        carbon_only::Bool=false,
        mimics_decomp::Bool=false,
        i_cop_mic::Int=0,
        i_oli_mic::Int=0,
        nitrif_n2o_loss_frac::Real=NITRIF_N2O_LOSS_FRAC,
        # FUN hook: when use_fun, this zero-arg closure runs after the per-layer
        # competition has set the *offered* N (smin_{nh4,no3}_to_plant_vr) and
        # before the FUN re-sum, to compute the cost-based *actual* uptake
        # (sminn_to_plant_fun_{nh4,no3}_vr) — i.e. it calls cnfun! + p2c. This
        # mirrors Fortran calling CNFUN inside SoilBiogeochemCompetition.
        fun_hook::Union{Function, Nothing}=nothing)

    dt   = state.dt
    bdnr = state.bdnr

    # Unpack state variables
    fpg              = st.fpg_col
    fpi              = st.fpi_col
    fpi_vr           = st.fpi_vr_col
    nfixation_prof   = st.nfixation_prof_col
    plant_ndemand    = st.plant_ndemand_col

    # Unpack nitrogen state
    sminn_vr     = ns.sminn_vr_col
    smin_nh4_vr  = ns.smin_nh4_vr_col
    smin_no3_vr  = ns.smin_no3_vr_col

    # Unpack nitrogen flux
    potential_immob_vr           = nf.potential_immob_vr_col
    actual_immob_vr              = nf.actual_immob_vr_col
    sminn_to_plant_vr            = nf.sminn_to_plant_vr_col
    sminn_to_plant               = nf.sminn_to_plant_col
    actual_immob                 = nf.actual_immob_col
    potential_immob              = nf.potential_immob_col
    supplement_to_sminn_vr       = nf.supplement_to_sminn_vr_col
    sminn_to_denit_excess_vr     = nf.sminn_to_denit_excess_vr_col
    actual_immob_no3_vr          = nf.actual_immob_no3_vr_col
    actual_immob_nh4_vr          = nf.actual_immob_nh4_vr_col
    smin_no3_to_plant_vr         = nf.smin_no3_to_plant_vr_col
    smin_nh4_to_plant_vr         = nf.smin_nh4_to_plant_vr_col
    pot_f_nit_vr                 = nf.pot_f_nit_vr_col
    pot_f_denit_vr               = nf.pot_f_denit_vr_col
    f_nit_vr                     = nf.f_nit_vr_col
    f_denit_vr                   = nf.f_denit_vr_col
    n2_n2o_ratio_denit_vr        = nf.n2_n2o_ratio_denit_vr_col
    f_n2o_denit_vr               = nf.f_n2o_denit_vr_col
    f_n2o_nit_vr                 = nf.f_n2o_nit_vr_col
    sminn_to_plant_fun_vr        = nf.sminn_to_plant_fun_vr_col
    sminn_to_plant_fun_no3_vr    = nf.sminn_to_plant_fun_no3_vr_col
    sminn_to_plant_fun_nh4_vr    = nf.sminn_to_plant_fun_nh4_vr_col

    # Unpack carbon flux (for MIMICS)
    c_overflow_vr = cf.c_overflow_vr

    local_use_fun = use_fun

    # Local arrays
    FT = eltype(sminn_vr)
    n_arr = length(mask_bgc_soilc)
    sminn_tot              = similar(sminn_vr, FT, n_arr); fill!(sminn_tot, zero(FT))
    sminn_to_plant_new     = similar(sminn_vr, FT, n_arr); fill!(sminn_to_plant_new, zero(FT))
    nuptake_prof           = similar(sminn_vr, FT, n_arr, nlevdecomp); fill!(nuptake_prof, zero(FT))
    sum_ndemand_vr         = similar(sminn_vr, FT, n_arr, nlevdecomp); fill!(sum_ndemand_vr, zero(FT))

    # Scalar params/thresholds converted to the state eltype at the launch site:
    # a Float64 scalar kernel arg is invalid Metal IR (rule 5).
    _T   = FT
    dt_k   = _T(dt)
    bdnr_k = _T(bdnr)

    if !use_nitrif_denitrif
        # ====================================================================
        # NON-NITRIF_DENITRIF PATHWAY
        # ====================================================================

        nlimit             = similar(sminn_vr, Int, n_arr, nlevdecomp); fill!(nlimit, 0)
        residual_sminn_vr  = similar(sminn_vr, FT, n_arr, nlevdecomp); fill!(residual_sminn_vr, zero(FT))
        residual_sminn     = similar(sminn_vr, FT, n_arr); fill!(residual_sminn, zero(FT))
        residual_plant_ndemand = similar(sminn_vr, FT, n_arr); fill!(residual_plant_ndemand, zero(FT))

        # sum total mineral N (per-column reduction over levels; sminn_tot starts 0)
        _launch!(_decompc_sminn_tot_kernel!, mask_bgc_soilc, sminn_tot,
                 sminn_vr, dzsoi_decomp, nlevdecomp)

        # define N uptake profile
        decompc_nuptake_prof!(nuptake_prof, mask_bgc_soilc, sminn_vr, sminn_tot,
                              nfixation_prof, nlevdecomp)

        # total N demand at each level
        decompc_sum_ndemand!(sum_ndemand_vr, mask_bgc_soilc, plant_ndemand,
                             nuptake_prof, potential_immob_vr, nlevdecomp)

        # resolve competition at each level
        decompc_resolve!(actual_immob_vr, mask_bgc_soilc, nlimit, fpi_vr,
                         sminn_to_plant_vr, supplement_to_sminn_vr, sum_ndemand_vr,
                         sminn_vr, potential_immob_vr, plant_ndemand, nuptake_prof,
                         dt_k, carbon_only, nlevdecomp)

        # sum up N fluxes to plant (per-column reduction + use_fun clamp)
        _launch!(_decompc_sum_to_plant1_kernel!, mask_bgc_soilc, sminn_to_plant,
                 sminn_to_plant_fun_vr, sminn_to_plant_vr, dzsoi_decomp,
                 nlevdecomp, local_use_fun)

        # --- Second pass: distribute residual N to plants ---
        _launch!(_decompc_residual_kernel!, mask_bgc_soilc, residual_sminn,
                 residual_plant_ndemand, residual_sminn_vr, plant_ndemand,
                 sminn_to_plant, sminn_vr, actual_immob_vr, sminn_to_plant_vr,
                 nlimit, dzsoi_decomp, nlevdecomp, dt_k)

        # distribute residual N to plants
        decompc_distribute_residual!(sminn_to_plant_vr, mask_bgc_soilc,
                                     residual_plant_ndemand, residual_sminn,
                                     residual_sminn_vr, nlimit, dt_k, nlevdecomp)

        # re-sum up N fluxes to plant + recompute sum_ndemand_vr
        _launch!(_decompc_resum_to_plant_kernel!, mask_bgc_soilc, sminn_to_plant,
                 sminn_to_plant_new, sum_ndemand_vr, sminn_to_plant_vr,
                 sminn_to_plant_fun_vr, potential_immob_vr, dzsoi_decomp,
                 nlevdecomp, local_use_fun)

        # excess denitrification
        decompc_denit_excess!(sminn_to_denit_excess_vr, mask_bgc_soilc,
                              sminn_to_plant_vr, sminn_to_plant_fun_vr, actual_immob_vr,
                              sminn_vr, sum_ndemand_vr, dt_k, bdnr_k, local_use_fun, nlevdecomp)

        # sum up N fluxes to immobilization + fraction of potential growth/immob
        _launch!(_decompc_immob_fpgfpi_kernel!, mask_bgc_soilc, actual_immob,
                 potential_immob, fpg, fpi, actual_immob_vr, potential_immob_vr,
                 plant_ndemand, sminn_to_plant, sminn_to_plant_new, dzsoi_decomp,
                 nlevdecomp, local_use_fun)

    else
        # ====================================================================
        # NITRIF_DENITRIF PATHWAY
        # ====================================================================

        # Competition parameters (converted to state eltype at launch site)
        compet_plant_no3  = _T(params.compet_plant_no3)
        compet_plant_nh4  = _T(params.compet_plant_nh4)
        compet_decomp_no3 = _T(params.compet_decomp_no3)
        compet_decomp_nh4 = _T(params.compet_decomp_nh4)
        compet_denit      = _T(params.compet_denit)
        compet_nit        = _T(params.compet_nit)
        nitrif_n2o_loss_frac_k = _T(nitrif_n2o_loss_frac)

        fpi_no3_vr_local       = similar(sminn_vr, FT, n_arr, nlevdecomp); fill!(fpi_no3_vr_local, zero(FT))
        fpi_nh4_vr_local       = similar(sminn_vr, FT, n_arr, nlevdecomp); fill!(fpi_nh4_vr_local, zero(FT))
        sum_nh4_demand         = similar(sminn_vr, FT, n_arr, nlevdecomp); fill!(sum_nh4_demand, zero(FT))
        sum_nh4_demand_scaled  = similar(sminn_vr, FT, n_arr, nlevdecomp); fill!(sum_nh4_demand_scaled, zero(FT))
        sum_no3_demand         = similar(sminn_vr, FT, n_arr, nlevdecomp); fill!(sum_no3_demand, zero(FT))
        sum_no3_demand_scaled  = similar(sminn_vr, FT, n_arr, nlevdecomp); fill!(sum_no3_demand_scaled, zero(FT))
        nlimit_no3             = similar(sminn_vr, Int, n_arr, nlevdecomp); fill!(nlimit_no3, 0)
        nlimit_nh4             = similar(sminn_vr, Int, n_arr, nlevdecomp); fill!(nlimit_nh4, 0)
        residual_plant_ndemand = similar(sminn_vr, FT, n_arr); fill!(residual_plant_ndemand, zero(FT))
        residual_smin_nh4_vr   = similar(sminn_vr, FT, n_arr, nlevdecomp); fill!(residual_smin_nh4_vr, zero(FT))
        residual_smin_no3_vr   = similar(sminn_vr, FT, n_arr, nlevdecomp); fill!(residual_smin_no3_vr, zero(FT))
        residual_smin_nh4      = similar(sminn_vr, FT, n_arr); fill!(residual_smin_nh4, zero(FT))
        residual_smin_no3      = similar(sminn_vr, FT, n_arr); fill!(residual_smin_no3, zero(FT))

        # sum up total mineral N pools (no3 + nh4) — per-column reduction
        _launch!(_decompc_nd_sminn_tot_kernel!, mask_bgc_soilc, sminn_tot,
                 smin_no3_vr, smin_nh4_vr, dzsoi_decomp, nlevdecomp)

        # define N uptake profile (per-(c,j) own write)
        _launch!(_decompc_nd_nuptake_prof_kernel!, nuptake_prof, mask_bgc_soilc,
                 sminn_vr, sminn_tot, nfixation_prof; ndrange = (length(mask_bgc_soilc), nlevdecomp))

        # main column/vertical loop (per-column kernel, internal j-loop)
        nd = _DECOMPND(;
            potential_immob_vr = potential_immob_vr,
            actual_immob_vr = actual_immob_vr, actual_immob_no3_vr = actual_immob_no3_vr,
            actual_immob_nh4_vr = actual_immob_nh4_vr,
            sminn_to_plant_vr = sminn_to_plant_vr, smin_no3_to_plant_vr = smin_no3_to_plant_vr,
            smin_nh4_to_plant_vr = smin_nh4_to_plant_vr,
            pot_f_nit_vr = pot_f_nit_vr, pot_f_denit_vr = pot_f_denit_vr,
            f_nit_vr = f_nit_vr, f_denit_vr = f_denit_vr,
            n2_n2o_ratio_denit_vr = n2_n2o_ratio_denit_vr,
            f_n2o_denit_vr = f_n2o_denit_vr, f_n2o_nit_vr = f_n2o_nit_vr,
            supplement_to_sminn_vr = supplement_to_sminn_vr, fpi_vr = fpi_vr,
            sminn_vr = sminn_vr, smin_nh4_vr = smin_nh4_vr, smin_no3_vr = smin_no3_vr,
            fpi_no3_vr_local = fpi_no3_vr_local, fpi_nh4_vr_local = fpi_nh4_vr_local,
            sum_nh4_demand = sum_nh4_demand, sum_nh4_demand_scaled = sum_nh4_demand_scaled,
            sum_no3_demand = sum_no3_demand, sum_no3_demand_scaled = sum_no3_demand_scaled,
            nlimit_no3 = nlimit_no3, nlimit_nh4 = nlimit_nh4)
        _launch!(_decompc_nd_main_kernel!, mask_bgc_soilc, nd, plant_ndemand,
                 nuptake_prof, dt_k,
                 compet_plant_no3, compet_plant_nh4, compet_decomp_no3, compet_decomp_nh4,
                 compet_denit, compet_nit, nitrif_n2o_loss_frac_k,
                 nlevdecomp, local_use_fun, mimics_decomp, carbon_only)

        # --- Post-competition sums (zero + non-fun sum) ---
        _launch!(_decompc_nd_postsum_kernel!, mask_bgc_soilc, sminn_to_plant,
                 sminn_to_plant_vr, dzsoi_decomp, nlevdecomp, local_use_fun)

        # --- MIMICS c_overflow_vr calculation ---
        if mimics_decomp
            _launch!(_decompc_nd_mimics_overflow_kernel!, mask_bgc_soilc, c_overflow_vr,
                     sum_ndemand_vr, sum_no3_demand_scaled, sum_nh4_demand_scaled,
                     pmnf_decomp_cascade, p_decomp_cn_gain, sminn_vr,
                     cascade_receiver_pool, nlevdecomp, ndecomp_cascade_transitions,
                     i_cop_mic, i_oli_mic, dt_k)
        else
            c_overflow_vr .= zero(FT)
        end

        # --- Second pass residual distribution ---
        if !local_use_fun
            # Second pass for NH4
            _launch!(_decompc_nd_residual_nh4_kernel!, mask_bgc_soilc, residual_plant_ndemand,
                     residual_smin_nh4, residual_smin_nh4_vr, smin_nh4_to_plant_vr,
                     plant_ndemand, sminn_to_plant, smin_nh4_vr, actual_immob_nh4_vr,
                     f_nit_vr, nlimit_nh4, dzsoi_decomp, nlevdecomp, dt_k)

            # re-sum N fluxes after second pass for nh4
            _launch!(_decompc_nd_resum_kernel!, mask_bgc_soilc, sminn_to_plant,
                     sminn_to_plant_vr, smin_nh4_to_plant_vr, smin_no3_to_plant_vr,
                     dzsoi_decomp, nlevdecomp)

            # Second pass for NO3
            _launch!(_decompc_nd_residual_no3_kernel!, mask_bgc_soilc, residual_plant_ndemand,
                     residual_smin_no3, residual_smin_no3_vr, smin_no3_to_plant_vr,
                     plant_ndemand, sminn_to_plant, smin_no3_vr, actual_immob_no3_vr,
                     f_denit_vr, nlimit_no3, dzsoi_decomp, nlevdecomp, dt_k)

            # re-sum N fluxes after second passes
            _launch!(_decompc_nd_resum_kernel!, mask_bgc_soilc, sminn_to_plant,
                     sminn_to_plant_vr, smin_nh4_to_plant_vr, smin_no3_to_plant_vr,
                     dzsoi_decomp, nlevdecomp)
        else
            # use_fun: the per-layer loop above set the *offered* N
            # (smin_{nh4,no3}_to_plant_vr = smin_*_vr/dt - immob - nit/denit).
            # Run FUN now to turn offered → cost-based actual uptake
            # (sminn_to_plant_fun_{nh4,no3}_vr), then sum the FUN fluxes.
            fun_hook !== nothing && fun_hook()
            _launch!(_decompc_nd_fun_resum_kernel!, mask_bgc_soilc, sminn_to_plant,
                     sminn_to_plant_new, sminn_to_plant_vr, smin_nh4_to_plant_vr,
                     smin_no3_to_plant_vr, sminn_to_plant_fun_no3_vr,
                     sminn_to_plant_fun_nh4_vr, dzsoi_decomp, nlevdecomp)
        end

        # --- Sum up N fluxes to immobilization + fraction of potential growth/immob ---
        _launch!(_decompc_nd_immob_fpgfpi_kernel!, mask_bgc_soilc, actual_immob,
                 potential_immob, fpg, fpi, actual_immob_vr, potential_immob_vr,
                 plant_ndemand, sminn_to_plant, dzsoi_decomp, nlevdecomp, local_use_fun)

    end  # end if_nitrif

    return nothing
end
