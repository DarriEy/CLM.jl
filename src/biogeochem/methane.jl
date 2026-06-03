# ==========================================================================
# Ported from: src/biogeochem/ch4Mod.F90 (4450 lines)
# Methane biogeochemistry module
#
# Calculates methane fluxes including production, oxidation, aerenchyma
# transport, ebullition, and diffusive transport through soil layers.
#
# Public types:
#   CH4Params       — Methane model parameters
#   CH4VarCon       — Methane control/configuration flags
#   CH4Data         — Methane state and flux data
#
# Public functions:
#   ch4_init_column_balance_check!   — Column balance check init
#   ch4_init_gridcell_balance_check! — Gridcell balance check init
#   ch4!                             — Main methane driver
#   ch4_prod!                        — Methane production
#   ch4_oxid!                        — Methane oxidation
#   ch4_aere!                        — Aerenchyma transport
#   site_ox_aere!                    — Site-level aerenchyma helper
#   ch4_ebul!                        — Ebullition
#   ch4_tran!                        — Reaction-diffusion transport solver
#   get_jwt!                         — Water table layer identification
#   ch4_annualupdate!                — Annual average update
#   ch4_totcolch4!                   — Total column CH4 calculation
# ==========================================================================

# ---------------------------------------------------------------------------
# Non-tunable constants
# ---------------------------------------------------------------------------
const rgasLatm = 0.0821  # L.atm/mol.K

# ---------------------------------------------------------------------------
# KernelAbstractions kernels for fully-independent per-(column, layer) loops.
# Module-level constants (TFRZ, C_H_INV, ...) are passed in as scalar/array
# args because GPU kernels cannot close over globals. Semantics are identical
# to the inline loops they replace (validated vs Fortran parity and AD).
# ---------------------------------------------------------------------------

# Ebullition rate per (column, soil layer). One thread per (c, j); each element
# writes only ch4_ebul_depth[c, j] from inputs at that index — fully independent.
@kernel function _meth_ebul_kernel!(ch4_ebul_depth, @Const(mask), @Const(jwt),
                                    @Const(t_soisno), @Const(conc_ch4), @Const(watsat),
                                    @Const(forc_pbot), @Const(h2osfc), @Const(frac_h2osfc),
                                    @Const(lake_icefrac), @Const(lakedepth),
                                    @Const(z), @Const(zi),
                                    sat::Int, lake::Bool, vgc_max, vgc_min, bubble_f,
                                    ebul_timescale, rgasm, rgaslatm,
                                    c_h_inv1, kh_theta1, kh_tbase, tfrz, denh2o, grav)
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        T = eltype(ch4_ebul_depth)
        if j > jwt[c] && t_soisno[c, j] > tfrz
            k_h_inv = exp(-c_h_inv1 * (one(T) / t_soisno[c, j] - one(T) / kh_tbase) + log(kh_theta1))
            k_h = one(T) / k_h_inv
            k_h_cc = t_soisno[c, j] * k_h * rgaslatm

            if !lake
                zi_jwt = jwt[c] > 0 ? zi[c, jwt[c]] : zero(T)
                pressure = forc_pbot[c] + denh2o * grav * (z[c, j] - zi_jwt)
                if sat == 1 && frac_h2osfc[c] > zero(T)
                    pressure += denh2o * grav * h2osfc[c] / T(1000.0) / frac_h2osfc[c]
                end
            else
                pressure = forc_pbot[c] + denh2o * grav * (z[c, j] + lakedepth[c])
            end

            vgc = conc_ch4[c, j] / watsat[c, j] / k_h_cc * rgasm * t_soisno[c, j] / pressure

            if vgc > vgc_max * bubble_f
                ch4_ebul_depth[c, j] = (vgc - vgc_min * bubble_f) * conc_ch4[c, j] / ebul_timescale
            else
                ch4_ebul_depth[c, j] = zero(T)
            end
        else
            ch4_ebul_depth[c, j] = zero(T)
        end

        if lake && lake_icefrac[c, 1] > T(0.1)
            ch4_ebul_depth[c, j] = zero(T)
        end
    end
end

function meth_ebul!(ch4_ebul_depth, mask, jwt, t_soisno, conc_ch4, watsat,
                    forc_pbot, h2osfc, frac_h2osfc, lake_icefrac, lakedepth, z, zi,
                    sat::Int, lake::Bool, vgc_max, vgc_min, bubble_f,
                    ebul_timescale, rgasm, rgaslatm, nc::Int, nlevsoi::Int)
    _T = eltype(ch4_ebul_depth)   # convert all scalar args to working precision (no Float64 on Metal)
    _launch!(_meth_ebul_kernel!, ch4_ebul_depth, mask, jwt, t_soisno, conc_ch4, watsat,
             forc_pbot, h2osfc, frac_h2osfc, lake_icefrac, lakedepth, z, zi,
             sat, lake, _T(vgc_max), _T(vgc_min), _T(bubble_f), _T(ebul_timescale),
             _T(rgasm), _T(rgasLatm), _T(C_H_INV[1]), _T(KH_THETA[1]), _T(KH_TBASE),
             _T(TFRZ), _T(DENH2O), _T(GRAV); ndrange = (nc, nlevsoi))
end

# Henry's-law solubility coefficients k_h_cc_arr[c, jj, s] for the transport
# solver. One thread per (c, jj) where jj == j+1 (Fortran j runs 0:nlevsoi, with
# j==0 using t_grnd). Each element is independent; the species loop (s=1:2) is a
# fixed 2-iteration inner loop with no cross-element coupling.
@kernel function _meth_khcc_kernel!(k_h_cc_arr, @Const(mask), @Const(t_grnd),
                                    @Const(t_soisno), @Const(c_h_inv), @Const(kh_theta),
                                    kh_tbase, rgaslatm)
    c, jj = @index(Global, NTuple)
    @inbounds if mask[c]
        T = eltype(k_h_cc_arr)
        for s in 1:2
            if jj == 1  # Fortran j == 0 → ground temperature
                k_h_inv = exp(-c_h_inv[s] * (one(T) / t_grnd[c] - one(T) / kh_tbase) + log(kh_theta[s]))
                k_h_cc_arr[c, jj, s] = t_grnd[c] / k_h_inv * rgaslatm
            else
                k_h_inv = exp(-c_h_inv[s] * (one(T) / t_soisno[c, jj-1] - one(T) / kh_tbase) + log(kh_theta[s]))
                k_h_cc_arr[c, jj, s] = t_soisno[c, jj-1] / k_h_inv * rgaslatm
            end
        end
    end
end

function meth_khcc!(k_h_cc_arr, mask, t_grnd, t_soisno, nc::Int, nlevsoi::Int)
    # Move the const Henry-law coefficient vectors onto the state backend at working
    # precision (host const arrays are non-bitstype on a device kernel).
    _T = eltype(k_h_cc_arr)
    chinv = similar(k_h_cc_arr, _T, length(C_H_INV)); copyto!(chinv, _T.(C_H_INV))
    khth  = similar(k_h_cc_arr, _T, length(KH_THETA)); copyto!(khth, _T.(KH_THETA))
    _launch!(_meth_khcc_kernel!, k_h_cc_arr, mask, t_grnd, t_soisno,
             chinv, khth, _T(KH_TBASE), _T(rgasLatm);
             ndrange = (nc, nlevsoi + 1))
end

# ---------------------------------------------------------------------------
# ch4varcon — Methane control flags
# ---------------------------------------------------------------------------

"""
    CH4VarCon

Methane model control flags and configuration.
Ported from `ch4varcon` module in CLM Fortran.
"""
Base.@kwdef mutable struct CH4VarCon
    allowlakeprod       ::Bool    = false   # allow CH4 production in lakes
    replenishlakec      ::Bool    = true    # replenish lake soil C (keep it constant)
    anoxicmicrosites    ::Bool    = true    # allow CH4 production above WT in anoxic microsites
    usephfact           ::Bool    = false   # use pH factor for CH4 production
    ch4rmcnlim          ::Bool    = true    # remove CN N limitation for methanogenesis
    ch4offline          ::Bool    = true    # use prescribed atmospheric CH4
    transpirationloss   ::Bool    = true    # include transpiration CH4 loss
    use_aereoxid_prog   ::Bool    = false   # prognostic aerenchyma oxidation via O2 diffusion
    ch4frzout           ::Bool    = false   # freeze-out effect on CH4 diffusion
    finundation_mtd_h2osfc ::Int  = 1       # finundation method h2osfc
end

# ---------------------------------------------------------------------------
# CH4Params — Methane model parameters (from params_type)
# ---------------------------------------------------------------------------

"""
    CH4Params

Methane model parameters read from parameter file.
Ported from `params_type` in `ch4Mod.F90`.
"""
Base.@kwdef mutable struct CH4Params
    # CH4 production constants
    q10ch4              ::Float64 = 1.5      # additional Q10 for CH4 production
    q10ch4base          ::Float64 = 295.0    # temperature at which effective f_ch4 equals constant f_ch4 (K)
    f_ch4               ::Float64 = 0.2      # ratio of CH4 production to total C mineralization
    rootlitfrac         ::Float64 = 0.5      # fraction of SOM associated with roots
    cnscalefactor       ::Float64 = 1.0      # scale factor on CN decomp for CH4 flux
    redoxlag            ::Float64 = 30.0     # days to lag in finundated_lag calculation
    lake_decomp_fact    ::Float64 = 2.0e-8   # base decomposition rate (1/s) at 25C
    redoxlag_vertical   ::Float64 = 30.0     # time lag (days) to inhibit production for newly unsaturated layers
    pHmax               ::Float64 = 9.0      # maximum pH for CH4 production
    pHmin               ::Float64 = 2.2      # minimum pH for CH4 production
    oxinhib             ::Float64 = 10.0     # inhibition of CH4 production by oxygen (m^3/mol)

    # CH4 oxidation constants
    vmax_ch4_oxid       ::Float64 = 0.0125   # oxidation rate constant [mol/m3-w/s]
    k_m                 ::Float64 = 5.0e-3   # Michaelis-Menten rate constant for CH4 conc
    q10_ch4oxid         ::Float64 = 1.9      # Q10 oxidation constant
    smp_crit            ::Float64 = -2.4e5   # critical soil moisture potential (mm)
    k_m_o2              ::Float64 = 2.0e-2   # Michaelis-Menten rate constant for O2 conc
    k_m_unsat           ::Float64 = 5.0e-3   # Michaelis-Menten rate constant for CH4 conc (unsat)
    vmax_oxid_unsat     ::Float64 = 1.25e-3  # oxidation vmax for unsaturated conditions [mol/m3-w/s]

    # CH4 aerenchyma constants
    aereoxid            ::Float64 = 0.5      # fraction of CH4 in aerenchyma oxidized
    scale_factor_aere   ::Float64 = 1.0      # scale factor on aerenchyma area
    nongrassporosratio  ::Float64 = 0.33     # ratio of root porosity non-grass/grass
    unsat_aere_ratio    ::Float64 = 0.167    # ratio for upland aerenchyma porosity (0.05/0.3)
    porosmin            ::Float64 = 0.05     # minimum aerenchyma porosity

    # CH4 ebullition constants
    vgc_max             ::Float64 = 0.15     # ratio of saturation pressure triggering ebullition

    # CH4 transport constants
    satpow              ::Float64 = 2.0      # exponent on watsat for saturated soil solute diffusion
    scale_factor_gasdiff::Float64 = 1.0      # scale factor for gas diffusion
    scale_factor_liqdiff::Float64 = 1.0      # scale factor for liquid diffusion
    capthick            ::Float64 = 100.0    # min thickness before h2osfc is impermeable (mm)

    # Additional constants
    f_sat               ::Float64 = 0.95     # volumetric soil water defining top of WT
    qflxlagd            ::Float64 = 30.0     # days to lag qflx_surf_lag in tropics
    highlatfact         ::Float64 = 2.0      # multiple of qflxlagd for high latitudes
    q10lakebase         ::Float64 = 298.0    # base temperature for lake CH4 production (K)
    atmch4              ::Float64 = 1.7e-6   # atmospheric CH4 mixing ratio (mol/mol)
    rob                 ::Float64 = 3.0      # ratio of root length to vertical depth
    om_frac_sf          ::Float64 = 1.0      # scale factor for organic matter fraction
end

# ---------------------------------------------------------------------------
# CH4Data — Methane state and flux data (from ch4_type)
# ---------------------------------------------------------------------------

"""
    CH4Data

Methane state and flux data arrays.
Ported from `ch4_type` in `ch4Mod.F90`.
"""
Base.@kwdef mutable struct CH4Data{FT<:Real, V<:AbstractVector{FT}, M<:AbstractMatrix{FT}, Vb<:AbstractVector{Bool}}
    # Per-layer production/consumption/transport rates (nc × nlevsoi)
    ch4_prod_depth_sat_col      ::M          = Matrix{FT}(undef, 0, 0)
    ch4_prod_depth_unsat_col    ::M          = Matrix{FT}(undef, 0, 0)
    ch4_prod_depth_lake_col     ::M          = Matrix{FT}(undef, 0, 0)
    ch4_oxid_depth_sat_col      ::M          = Matrix{FT}(undef, 0, 0)
    ch4_oxid_depth_unsat_col    ::M          = Matrix{FT}(undef, 0, 0)
    ch4_oxid_depth_lake_col     ::M          = Matrix{FT}(undef, 0, 0)
    ch4_aere_depth_sat_col      ::M          = Matrix{FT}(undef, 0, 0)
    ch4_aere_depth_unsat_col    ::M          = Matrix{FT}(undef, 0, 0)
    ch4_tran_depth_sat_col      ::M          = Matrix{FT}(undef, 0, 0)
    ch4_tran_depth_unsat_col    ::M          = Matrix{FT}(undef, 0, 0)
    ch4_ebul_depth_sat_col      ::M          = Matrix{FT}(undef, 0, 0)
    ch4_ebul_depth_unsat_col    ::M          = Matrix{FT}(undef, 0, 0)
    o2_oxid_depth_sat_col       ::M          = Matrix{FT}(undef, 0, 0)
    o2_oxid_depth_unsat_col     ::M          = Matrix{FT}(undef, 0, 0)
    o2_aere_depth_sat_col       ::M          = Matrix{FT}(undef, 0, 0)
    o2_aere_depth_unsat_col     ::M          = Matrix{FT}(undef, 0, 0)
    co2_decomp_depth_sat_col    ::M          = Matrix{FT}(undef, 0, 0)
    co2_decomp_depth_unsat_col  ::M          = Matrix{FT}(undef, 0, 0)
    co2_oxid_depth_sat_col      ::M          = Matrix{FT}(undef, 0, 0)
    co2_oxid_depth_unsat_col    ::M          = Matrix{FT}(undef, 0, 0)
    co2_aere_depth_sat_col      ::M          = Matrix{FT}(undef, 0, 0)
    co2_aere_depth_unsat_col    ::M          = Matrix{FT}(undef, 0, 0)

    # Column-level surface fluxes (nc)
    ch4_ebul_total_sat_col      ::V          = FT[]
    ch4_ebul_total_unsat_col    ::V          = FT[]
    ch4_surf_aere_sat_col       ::V          = FT[]
    ch4_surf_aere_unsat_col     ::V          = FT[]
    ch4_surf_ebul_sat_col       ::V          = FT[]
    ch4_surf_ebul_unsat_col     ::V          = FT[]
    ch4_surf_ebul_lake_col      ::V          = FT[]
    ch4_surf_diff_sat_col       ::V          = FT[]
    ch4_surf_diff_unsat_col     ::V          = FT[]
    ch4_surf_diff_lake_col      ::V          = FT[]
    ch4_dfsat_flux_col          ::V          = FT[]
    ch4_surf_flux_tot_col       ::V          = FT[]

    # Concentrations (nc × nlevsoi)
    conc_ch4_sat_col            ::M          = Matrix{FT}(undef, 0, 0)
    conc_ch4_unsat_col          ::M          = Matrix{FT}(undef, 0, 0)
    conc_ch4_lake_col           ::M          = Matrix{FT}(undef, 0, 0)
    conc_o2_sat_col             ::M          = Matrix{FT}(undef, 0, 0)
    conc_o2_unsat_col           ::M          = Matrix{FT}(undef, 0, 0)
    conc_o2_lake_col            ::M          = Matrix{FT}(undef, 0, 0)
    o2_decomp_depth_sat_col     ::M          = Matrix{FT}(undef, 0, 0)
    o2_decomp_depth_unsat_col   ::M          = Matrix{FT}(undef, 0, 0)

    # Stress factors (nc × nlevsoi)
    o2stress_sat_col            ::M          = Matrix{FT}(undef, 0, 0)
    o2stress_unsat_col          ::M          = Matrix{FT}(undef, 0, 0)
    ch4stress_sat_col           ::M          = Matrix{FT}(undef, 0, 0)
    ch4stress_unsat_col         ::M          = Matrix{FT}(undef, 0, 0)

    # Column-level state (nc)
    zwt_ch4_unsat_col           ::V          = FT[]
    lake_soilc_col              ::M          = Matrix{FT}(undef, 0, 0)
    totcolch4_col               ::V          = FT[]
    totcolch4_bef_col           ::V          = FT[]
    annsum_counter_col          ::V          = FT[]
    tempavg_somhr_col           ::V          = FT[]
    annavg_somhr_col            ::V          = FT[]
    tempavg_finrw_col           ::V          = FT[]
    annavg_finrw_col            ::V          = FT[]
    sif_col                     ::V          = FT[]
    qflx_surf_lag_col           ::V          = FT[]
    finundated_col              ::V          = FT[]
    finundated_pre_snow_col     ::V          = FT[]
    finundated_lag_col          ::V          = FT[]
    layer_sat_lag_col           ::M          = Matrix{FT}(undef, 0, 0)
    pH_col                      ::V          = FT[]

    # Gridcell-level (ng)
    c_atm_grc                   ::M          = Matrix{FT}(undef, 0, 0)
    ch4co2f_grc                 ::V          = FT[]
    ch4prodg_grc                ::V          = FT[]
    totcolch4_grc               ::V          = FT[]
    totcolch4_bef_grc           ::V          = FT[]

    # Patch-level aerenchyma (np)
    annavg_agnpp_patch          ::V          = FT[]
    annavg_bgnpp_patch          ::V          = FT[]
    tempavg_agnpp_patch         ::V          = FT[]
    tempavg_bgnpp_patch         ::V          = FT[]

    # Boundary layer conductance
    grnd_ch4_cond_patch         ::V          = FT[]
    grnd_ch4_cond_col           ::V          = FT[]

    # First-time flag (ng)
    ch4_first_time_grc          ::Vb               = Bool[]
end

# Array-type params V/M/Vb stay loose so adapt(MtlArray/CuArray, ch4) can swap the
# storage to device arrays (and AD can swap Dual arrays in). The `{FT}` convenience
# ctor defaults them to Vector{FT}/Matrix{FT}/Vector{Bool} at the working precision.
CH4Data{FT}(; kwargs...) where {FT<:Real} =
    CH4Data{FT, Vector{FT}, Matrix{FT}, Vector{Bool}}(; kwargs...)
Adapt.@adapt_structure CH4Data

# ---------------------------------------------------------------------------
# get_jwt! — Water table layer identification
# ---------------------------------------------------------------------------

# get_jwt!: one thread per column; the water-table search (frozen-saturated
# perched scan + ascending unsaturated boundary scan) runs as in-thread
# sequential j-loops writing only the thread's own jwt[c]. Byte-identical to the
# masked scalar loop. f_sat/tfrz are passed at the working precision.
@kernel function _ch4diag_jwt_kernel!(jwt, @Const(mask), @Const(watsat),
                                      @Const(h2osoi_vol), @Const(t_soisno),
                                      nlevsoi::Int, f_sat, tfrz)
    c = @index(Global)
    @inbounds if mask[c]
        # Check for frozen saturated layers → perched water table
        perch = nlevsoi
        for j in nlevsoi:-1:1
            if t_soisno[c, j] < tfrz && h2osoi_vol[c, j] > f_sat * watsat[c, j]
                perch = j - 1
            end
        end
        jwt[c] = perch

        for j in perch:-1:2
            if h2osoi_vol[c, j] > f_sat * watsat[c, j] && h2osoi_vol[c, j-1] < f_sat * watsat[c, j-1]
                jwt[c] = j - 1
                break
            end
        end
        if jwt[c] == perch && h2osoi_vol[c, 1] > f_sat * watsat[c, 1]
            jwt[c] = 0
        end
    end
end

"""
    get_jwt!(jwt, mask_soil, watsat, h2osoi_vol, t_soisno, nlevsoi, params)

Find the first unsaturated layer going up (layer right above water table).
Also allows perched water table over ice.
Ported from `get_jwt` in `ch4Mod.F90`.
"""
function get_jwt!(jwt::AbstractVector{<:Integer},
                  mask_soil::AbstractVector{Bool},
                  watsat::AbstractMatrix{<:Real},
                  h2osoi_vol::AbstractMatrix{<:Real},
                  t_soisno::AbstractMatrix{<:Real},
                  nlevsoi::Int,
                  params::CH4Params)

    FT = eltype(watsat)
    _launch!(_ch4diag_jwt_kernel!, jwt, mask_soil, watsat, h2osoi_vol, t_soisno,
             nlevsoi, FT(params.f_sat), FT(TFRZ); ndrange = length(mask_soil))
    nothing
end

# ---------------------------------------------------------------------------
# ch4_annualupdate! — Annual average update
# ---------------------------------------------------------------------------

# ch4_annualupdate! kernels. All own-index (no cross-thread coupling):
#   1) per-column: bump the annual-sum counter, then the somhr/finrw column update
#      (the counter bump must complete before the >=secsperyear branch reads it,
#       so it is a separate kernel from the update — matches the scalar two-loop
#       order exactly).
#   2) per-patch: agnpp/bgnpp annual update (column flag/counter read-only).
#   3) per-column: reset the counter when a year has elapsed.
# Literals 0.0 become zero(T); dt/secsperyear passed at the working precision.
@kernel function _ch4annual_counter_kernel!(annsum_counter, @Const(mask), dt)
    c = @index(Global)
    @inbounds if mask[c]
        annsum_counter[c] += dt
    end
end

@kernel function _ch4annual_col_kernel!(annsum_counter, annavg_somhr, tempavg_somhr,
                                        annavg_finrw, tempavg_finrw, @Const(finundated),
                                        @Const(mask), @Const(somhr),
                                        dt, secsperyear)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(annavg_somhr)
        if annsum_counter[c] >= secsperyear
            annavg_somhr[c] = tempavg_somhr[c]
            tempavg_somhr[c] = zero(T)
            if annavg_somhr[c] > zero(T)
                annavg_finrw[c] = tempavg_finrw[c] / annavg_somhr[c]
            else
                annavg_finrw[c] = zero(T)
            end
            tempavg_finrw[c] = zero(T)
        else
            tempavg_somhr[c] += dt / secsperyear * somhr[c]
            tempavg_finrw[c] += dt / secsperyear * finundated[c] * somhr[c]
        end
    end
end

@kernel function _ch4annual_patch_kernel!(annavg_agnpp, tempavg_agnpp, annavg_bgnpp,
                                          tempavg_bgnpp, @Const(annsum_counter),
                                          @Const(mask), @Const(patch_column),
                                          @Const(is_fates), @Const(agnpp), @Const(bgnpp),
                                          dt, secsperyear)
    p = @index(Global)
    @inbounds if mask[p]
        c = patch_column[p]
        if !is_fates[c]
            T = eltype(annavg_agnpp)
            if annsum_counter[c] >= secsperyear
                annavg_agnpp[p] = tempavg_agnpp[p]
                tempavg_agnpp[p] = zero(T)
                annavg_bgnpp[p] = tempavg_bgnpp[p]
                tempavg_bgnpp[p] = zero(T)
            else
                tempavg_agnpp[p] += dt / secsperyear * agnpp[p]
                tempavg_bgnpp[p] += dt / secsperyear * bgnpp[p]
            end
        end
    end
end

@kernel function _ch4annual_reset_kernel!(annsum_counter, @Const(mask), secsperyear)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(annsum_counter)
        if annsum_counter[c] >= secsperyear
            annsum_counter[c] = zero(T)
        end
    end
end

"""
    ch4_annualupdate!(ch4, mask_soil, mask_soilp, patch_column, is_fates,
                      somhr, finundated, agnpp, bgnpp,
                      dt, secsperyear, nlevsoi)

Update annual mean fields for methane-related variables.
Ported from `ch4_annualupdate` in `ch4Mod.F90`.
"""
function ch4_annualupdate!(ch4::CH4Data,
                           mask_soil::AbstractVector{Bool},
                           mask_soilp::AbstractVector{Bool},
                           patch_column::AbstractVector{<:Integer},
                           is_fates::AbstractVector{Bool},
                           somhr::AbstractVector{<:Real},
                           agnpp::AbstractVector{<:Real},
                           bgnpp::AbstractVector{<:Real},
                           dt::Real,
                           secsperyear::Real)

    FT  = eltype(ch4.annsum_counter_col)
    dtF = FT(dt)
    syF = FT(secsperyear)

    _launch!(_ch4annual_counter_kernel!, ch4.annsum_counter_col, mask_soil, dtF;
             ndrange = length(mask_soil))

    _launch!(_ch4annual_col_kernel!, ch4.annsum_counter_col, ch4.annavg_somhr_col,
             ch4.tempavg_somhr_col, ch4.annavg_finrw_col, ch4.tempavg_finrw_col,
             ch4.finundated_col, mask_soil, somhr, dtF, syF;
             ndrange = length(mask_soil))

    _launch!(_ch4annual_patch_kernel!, ch4.annavg_agnpp_patch, ch4.tempavg_agnpp_patch,
             ch4.annavg_bgnpp_patch, ch4.tempavg_bgnpp_patch, ch4.annsum_counter_col,
             mask_soilp, patch_column, is_fates, agnpp, bgnpp, dtF, syF;
             ndrange = length(mask_soilp))

    _launch!(_ch4annual_reset_kernel!, ch4.annsum_counter_col, mask_soil, syF;
             ndrange = length(mask_soil))
    nothing
end

# ---------------------------------------------------------------------------
# ch4_totcolch4! — Total column CH4 calculation
# ---------------------------------------------------------------------------

# ch4_totcolch4!: per-column REDUCTION over soil layers. Each thread accumulates
# DIRECTLY into its own totcolch4[c] via an ascending in-thread j-loop (race-free,
# byte-identical to the j-outer/c-inner scalar order — addition for a single c is
# performed in the same ascending-j sequence). Nolake and lake are disjoint masks
# handled by two kernels; lake only contributes when allowlakeprod. CATOMW/1.0 are
# converted to the working precision.
@kernel function _ch4tot_nolake_kernel!(totcolch4, @Const(mask), @Const(finundated),
                                        @Const(conc_ch4_sat), @Const(conc_ch4_unsat),
                                        @Const(dz), nlevsoi::Int, catomw)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(totcolch4)
        totcolch4[c] = zero(T)
        for j in 1:nlevsoi
            totcolch4[c] += (finundated[c] * conc_ch4_sat[c, j] +
                             (one(T) - finundated[c]) * conc_ch4_unsat[c, j]) *
                            dz[c, j] * catomw
        end
    end
end

@kernel function _ch4tot_lake_kernel!(totcolch4, @Const(mask), @Const(conc_ch4_sat),
                                      @Const(dz), nlevsoi::Int, allowlakeprod::Bool, catomw)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(totcolch4)
        totcolch4[c] = zero(T)
        if allowlakeprod
            for j in 1:nlevsoi
                totcolch4[c] += conc_ch4_sat[c, j] * dz[c, j] * catomw
            end
        end
    end
end

"""
    ch4_totcolch4!(totcolch4, ch4, mask_nolake, mask_lake,
                   dz, col_landunit, lun_itype, nlevsoi, allowlakeprod)

Compute total column CH4 by integrating concentrations across soil layers.
Ported from `ch4_totcolch4` in `ch4Mod.F90`.
"""
function ch4_totcolch4!(totcolch4::AbstractVector{<:Real},
                        ch4::CH4Data,
                        mask_nolake::AbstractVector{Bool},
                        mask_lake::AbstractVector{Bool},
                        dz::AbstractMatrix{<:Real},
                        nlevsoi::Int,
                        allowlakeprod::Bool)

    FT = eltype(totcolch4)
    catomwF = FT(CATOMW)

    _launch!(_ch4tot_nolake_kernel!, totcolch4, mask_nolake, ch4.finundated_col,
             ch4.conc_ch4_sat_col, ch4.conc_ch4_unsat_col, dz, nlevsoi, catomwF;
             ndrange = length(mask_nolake))

    _launch!(_ch4tot_lake_kernel!, totcolch4, mask_lake, ch4.conc_ch4_sat_col,
             dz, nlevsoi, allowlakeprod, catomwF; ndrange = length(mask_lake))
    nothing
end

# ---------------------------------------------------------------------------
# ch4_init_column_balance_check!
# ---------------------------------------------------------------------------

"""
    ch4_init_column_balance_check!(ch4, mask_nolake, mask_lake, dz, nlevsoi, allowlakeprod)

Calculate beginning column-level CH4 balance.
Ported from `ch4_init_column_balance_check` in `ch4Mod.F90`.
"""
function ch4_init_column_balance_check!(ch4::CH4Data,
                                        mask_nolake::BitVector,
                                        mask_lake::BitVector,
                                        dz::Matrix{<:Real},
                                        nlevsoi::Int,
                                        allowlakeprod::Bool)
    ch4_totcolch4!(ch4.totcolch4_bef_col, ch4, mask_nolake, mask_lake,
                   dz, nlevsoi, allowlakeprod)
    nothing
end

# ---------------------------------------------------------------------------
# ch4_prod! — Methane production : KernelAbstractions kernels + device-view bundles
#
# Two kernels mirror the original two loop nests verbatim:
#   1. _ch4prod_rrvr_kernel!  — per-PATCH scatter of root respiration into the
#      column-resolved rr_vr[c,j] (many patches -> one column) via _scatter_add!.
#      rr_vr is zeroed by _ch4prod_zero_rrvr_kernel! first (masked-zero).
#   2. _ch4prod_main_kernel!  — per-(column, level) production, reading rr_vr.
#
# The main loop touches ~14 CH4Data array fields + several loose arrays, so the
# CH4Data views are bundled into immutable _Ch4ProdState{V,M} (field names mirror
# CH4Data so the body reads verbatim); scalar params/flags are bundled into
# _Ch4ProdScal{T} and _Ch4ProdFlags. Every Float64 literal/const is converted to
# eltype(out) so no Float64 reaches a Float32-only backend (Metal).
# ---------------------------------------------------------------------------

# Bundle of CH4Data array views read/written by the main production kernel.
# Vectors (V) and matrices (M) kept as distinct type params.
Base.@kwdef struct _Ch4ProdState{V,M}
    sif_col              ::V
    annavg_finrw_col     ::V
    finundated_col       ::V
    finundated_lag_col   ::V
    pH_col               ::V
    lake_soilc_col       ::M
    layer_sat_lag_col    ::M
    conc_o2              ::M
    ch4_prod_depth       ::M
    o2_decomp_depth      ::M
end
Adapt.@adapt_structure _Ch4ProdState

# Bundle of scalar Float64 params (converted to working precision at launch).
Base.@kwdef struct _Ch4ProdScal{T}
    q10ch4base    ::T
    f_ch4         ::T
    rootlitfrac   ::T
    cnscalefactor ::T
    lake_decomp_fact ::T
    pHmax         ::T
    pHmin         ::T
    oxinhib       ::T
    q10ch4        ::T
    q10lake       ::T
    q10lakebase   ::T
    mino2lim      ::T
    catomw        ::T
    tfrz          ::T
    spval         ::T
end
Adapt.@adapt_structure _Ch4ProdScal

# Bundle of Bool/Int config flags (stay scalar args — small enough to pass loose,
# but bundled to respect the Metal ~31-arg cap together with the array bundles).
Base.@kwdef struct _Ch4ProdFlags
    sat                ::Int
    lake               ::Bool
    use_cn             ::Bool
    use_nitrif_denitrif::Bool
    anoxia             ::Bool
    usephfact          ::Bool
    ch4rmcnlim         ::Bool
    anoxicmicrosites   ::Bool
    noveg              ::Int
    nlevdecomp         ::Int
    nlevdecomp_full    ::Int
    nlev_soildecomp_standard ::Int
end

# Masked-zero of rr_vr[c, j] (per-(column, level)).
@kernel function _ch4prod_zero_rrvr_kernel!(rr_vr, @Const(mask))
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        rr_vr[c, j] = zero(eltype(rr_vr))
    end
end

# Per-PATCH scatter: rr_vr[col, j] += rr[p]*crootfr[p,j]*patch_wtcol[p].
# col = patch_column[p]; MANY patches -> ONE column -> atomic via _scatter_add!.
@kernel function _ch4prod_rrvr_kernel!(rr_vr, @Const(mask_soilp),
                                       @Const(patch_column), @Const(patch_itype),
                                       @Const(patch_wtcol), @Const(is_fates),
                                       @Const(rr), @Const(crootfr),
                                       noveg::Int, nlevsoi::Int)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        c = patch_column[p]
        if !is_fates[c]
            if patch_wtcol[p] > zero(eltype(patch_wtcol)) && patch_itype[p] != noveg
                for j in 1:nlevsoi
                    _scatter_add!(rr_vr, c, j, rr[p] * crootfr[p, j] * patch_wtcol[p])
                end
            end
        end
    end
end

# Per-(column, level) production. Body copied verbatim from the scalar loop, with
# every Float64 literal/const -> T() and CH4Data fields read through the bundle `S`.
@kernel function _ch4prod_main_kernel!(@Const(mask), S::_Ch4ProdState, P::_Ch4ProdScal,
                                       F::_Ch4ProdFlags,
                                       @Const(somhr), @Const(lithr), @Const(hr_vr),
                                       @Const(o_scalar), @Const(fphr),
                                       @Const(pot_f_nit_vr), @Const(rr_vr),
                                       @Const(rootfr_col), @Const(t_soisno),
                                       @Const(dz), @Const(zi), @Const(jwt))
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        T = eltype(S.ch4_prod_depth)

        base_decomp = zero(T)
        partition_z = one(T)

        if !F.lake
            if F.use_cn
                base_decomp = (somhr[c] + lithr[c]) / P.catomw
                if F.sat == 1
                    S.sif_col[c] = one(T)
                    if !F.anoxia
                        if S.annavg_finrw_col[c] != P.spval
                            seasonalfin = smooth_max(S.finundated_col[c] - S.annavg_finrw_col[c], zero(T))
                            if seasonalfin > zero(T)
                                S.sif_col[c] = (S.annavg_finrw_col[c] + P.mino2lim * seasonalfin) / S.finundated_col[c]
                                base_decomp *= S.sif_col[c]
                            end
                        end
                    end
                end
            end
            base_decomp *= P.cnscalefactor
        else
            base_decomp = P.lake_decomp_fact * S.lake_soilc_col[c, j] * dz[c, j] *
                          P.q10lake^((t_soisno[c, j] - P.q10lakebase) / T(10)) / P.catomw
        end

        if t_soisno[c, j] <= P.tfrz && (F.nlevdecomp == 1 || F.lake)
            base_decomp = zero(T)
        end

        # Depth dependence
        if !F.lake
            if F.nlevdecomp == 1
                if j <= F.nlev_soildecomp_standard
                    partition_z = rootfr_col[c, j] * P.rootlitfrac +
                                  (one(T) - P.rootlitfrac) * dz[c, j] / zi[c, F.nlev_soildecomp_standard]
                else
                    partition_z = rootfr_col[c, j] * P.rootlitfrac
                end
            else
                if (somhr[c] + lithr[c]) > zero(T)
                    partition_z = hr_vr[c, j] * dz[c, j] / (somhr[c] + lithr[c])
                else
                    partition_z = one(T)
                end
            end
        else
            partition_z = one(T)
        end

        # Adjust f_ch4 for temperature
        f_ch4_adj = one(T)
        if !F.lake
            t_fact_ch4 = P.q10ch4^((t_soisno[c, j] - P.q10ch4base) / T(10))
            f_ch4_adj = P.f_ch4 * t_fact_ch4
            if F.ch4rmcnlim
                jj = j > F.nlevdecomp ? 1 : j
                if fphr[c, jj] > zero(T)
                    f_ch4_adj /= fphr[c, jj]
                end
            end
        else
            f_ch4_adj = T(0.5)
        end

        # pH factor
        if !F.lake && F.usephfact
            if S.pH_col[c] > P.pHmin && S.pH_col[c] < P.pHmax
                pH_fact_ch4 = T(10)^(-T(0.2235) * S.pH_col[c]^2 + T(2.7727) * S.pH_col[c] - T(8.6))
                f_ch4_adj *= pH_fact_ch4
            end
        end

        # Redox factor
        if !F.lake && F.sat == 1 && S.finundated_lag_col[c] < S.finundated_col[c]
            f_ch4_adj *= S.finundated_lag_col[c] / S.finundated_col[c]
        elseif F.sat == 0 && j > jwt[c]
            f_ch4_adj *= S.layer_sat_lag_col[c, j]
        end

        f_ch4_adj = smooth_min(f_ch4_adj, T(0.5))

        # O2 decomposition demand
        S.o2_decomp_depth[c, j] = base_decomp * partition_z / dz[c, j]
        if F.anoxia
            if !F.lake && j > F.nlevdecomp
                if o_scalar[c, 1] > zero(T)
                    S.o2_decomp_depth[c, j] /= o_scalar[c, 1]
                end
            elseif !F.lake
                if o_scalar[c, j] > zero(T)
                    S.o2_decomp_depth[c, j] /= o_scalar[c, j]
                end
            end
        end

        # Add root respiration
        if !F.lake
            S.o2_decomp_depth[c, j] += rr_vr[c, j] / P.catomw / dz[c, j]
        end

        # Add nitrification O2 demand
        if F.use_nitrif_denitrif && !F.lake && j <= F.nlevdecomp_full
            S.o2_decomp_depth[c, j] += pot_f_nit_vr[c, j] * T(2) / T(14)
        end

        # CH4 production
        if j > jwt[c]
            S.ch4_prod_depth[c, j] = f_ch4_adj * base_decomp * partition_z / dz[c, j]
        else
            if F.anoxicmicrosites
                S.ch4_prod_depth[c, j] = f_ch4_adj * base_decomp * partition_z / dz[c, j] /
                                         (one(T) + P.oxinhib * S.conc_o2[c, j])
            else
                S.ch4_prod_depth[c, j] = zero(T)
            end
        end
    end
end

# ---------------------------------------------------------------------------
# ch4_prod! — Methane production
# ---------------------------------------------------------------------------

"""
    ch4_prod!(ch4, params, ch4vc, mask_soil, mask_soilp,
              patch_column, patch_itype, patch_wtcol, is_fates,
              crootfr, rootfr_col, watsat, h2osoi_vol, t_soisno,
              somhr, lithr, hr_vr, o_scalar, fphr, pot_f_nit_vr,
              rr, jwt, sat, lake,
              dz, z, zi, nlevsoi, nlevdecomp, nlevdecomp_full,
              nlev_soildecomp_standard, mino2lim,
              dtime, noveg, use_cn, use_nitrif_denitrif, anoxia)

Calculate CH4 production in each soil layer.
Ported from `ch4_prod` in `ch4Mod.F90`.
"""
function ch4_prod!(ch4::CH4Data,
                   params::CH4Params,
                   ch4vc::CH4VarCon,
                   mask_soil::AbstractVector{Bool},
                   mask_soilp::AbstractVector{Bool},
                   patch_column::AbstractVector{<:Integer},
                   patch_itype::AbstractVector{<:Integer},
                   patch_wtcol::AbstractVector{<:Real},
                   is_fates::AbstractVector{Bool},
                   crootfr::AbstractMatrix{<:Real},
                   rootfr_col::AbstractMatrix{<:Real},
                   watsat::AbstractMatrix{<:Real},
                   h2osoi_vol::AbstractMatrix{<:Real},
                   t_soisno::AbstractMatrix{<:Real},
                   somhr::AbstractVector{<:Real},
                   lithr::AbstractVector{<:Real},
                   hr_vr::AbstractMatrix{<:Real},
                   o_scalar::AbstractMatrix{<:Real},
                   fphr::AbstractMatrix{<:Real},
                   pot_f_nit_vr::AbstractMatrix{<:Real},
                   rr::AbstractVector{<:Real},
                   jwt::AbstractVector{<:Integer},
                   sat::Int,
                   lake::Bool,
                   dz::AbstractMatrix{<:Real},
                   z::AbstractMatrix{<:Real},
                   zi::AbstractMatrix{<:Real},
                   nlevsoi::Int,
                   nlevdecomp::Int,
                   nlevdecomp_full::Int,
                   nlev_soildecomp_standard::Int,
                   mino2lim::Real,
                   dtime::Real,
                   noveg::Int,
                   use_cn::Bool,
                   use_nitrif_denitrif::Bool,
                   anoxia::Bool)

    # Select sat/unsat arrays
    if sat == 0
        conc_o2 = ch4.conc_o2_unsat_col
        ch4_prod_depth = ch4.ch4_prod_depth_unsat_col
        o2_decomp_depth = ch4.o2_decomp_depth_unsat_col
        co2_decomp_depth = ch4.co2_decomp_depth_unsat_col
    else
        conc_o2 = ch4.conc_o2_sat_col
        ch4_prod_depth = ch4.ch4_prod_depth_sat_col
        o2_decomp_depth = ch4.o2_decomp_depth_sat_col
        co2_decomp_depth = ch4.co2_decomp_depth_sat_col
    end

    q10lake = params.q10ch4 * 1.5

    nc = length(mask_soil)
    FT = eltype(t_soisno)

    # Calculate vertically resolved column-averaged root respiration.
    # Zero rr_vr (device-resident scratch), then scatter the per-patch root
    # respiration into the column-resolved rr_vr[c, j] (many patches -> one column).
    rr_vr = fill!(similar(t_soisno, FT, nc, nlevsoi), zero(FT))
    if !lake
        _launch!(_ch4prod_zero_rrvr_kernel!, rr_vr, mask_soil;
                 ndrange = (nc, nlevsoi))
        _launch!(_ch4prod_rrvr_kernel!, rr_vr, mask_soilp,
                 patch_column, patch_itype, patch_wtcol, is_fates,
                 rr, crootfr, noveg, nlevsoi;
                 ndrange = length(mask_soilp))
    end

    # Bundle CH4Data views, scalar params (converted to FT), and flags.
    S = _Ch4ProdState{typeof(ch4.sif_col),typeof(ch4_prod_depth)}(
        sif_col            = ch4.sif_col,
        annavg_finrw_col   = ch4.annavg_finrw_col,
        finundated_col     = ch4.finundated_col,
        finundated_lag_col = ch4.finundated_lag_col,
        pH_col             = ch4.pH_col,
        lake_soilc_col     = ch4.lake_soilc_col,
        layer_sat_lag_col  = ch4.layer_sat_lag_col,
        conc_o2            = conc_o2,
        ch4_prod_depth     = ch4_prod_depth,
        o2_decomp_depth    = o2_decomp_depth)

    P = _Ch4ProdScal{FT}(
        q10ch4base       = FT(params.q10ch4base),
        f_ch4            = FT(params.f_ch4),
        rootlitfrac      = FT(params.rootlitfrac),
        cnscalefactor    = FT(params.cnscalefactor),
        lake_decomp_fact = FT(params.lake_decomp_fact),
        pHmax            = FT(params.pHmax),
        pHmin            = FT(params.pHmin),
        oxinhib          = FT(params.oxinhib),
        q10ch4           = FT(params.q10ch4),
        q10lake          = FT(q10lake),
        q10lakebase      = FT(params.q10lakebase),
        mino2lim         = FT(mino2lim),
        catomw           = FT(CATOMW),
        tfrz             = FT(TFRZ),
        spval            = FT(SPVAL))

    F = _Ch4ProdFlags(
        sat                = sat,
        lake               = lake,
        use_cn             = use_cn,
        use_nitrif_denitrif = use_nitrif_denitrif,
        anoxia             = anoxia,
        usephfact          = ch4vc.usephfact,
        ch4rmcnlim         = ch4vc.ch4rmcnlim,
        anoxicmicrosites   = ch4vc.anoxicmicrosites,
        noveg              = noveg,
        nlevdecomp         = nlevdecomp,
        nlevdecomp_full    = nlevdecomp_full,
        nlev_soildecomp_standard = nlev_soildecomp_standard)

    # Per-(column, level) production kernel (mask first so it owns the backend).
    _launch!(_ch4prod_main_kernel!, mask_soil, S, P, F,
             somhr, lithr, hr_vr, o_scalar, fphr, pot_f_nit_vr, rr_vr,
             rootfr_col, t_soisno, dz, zi, jwt;
             ndrange = (nc, nlevsoi))
    nothing
end

# ---------------------------------------------------------------------------
# ch4_oxid! — Methane oxidation
# ---------------------------------------------------------------------------

# Per-(column, soil layer) CH4 oxidation via double Michaelis-Menten kinetics.
# One thread per (c, j); each thread writes only ch4_oxid_depth[c,j] and
# o2_oxid_depth[c,j] from inputs at its own index — fully independent (no
# reduction, no scatter). The sat/unsat variant is selected by the `sat` Int
# and per-column water-table index jwt[c], passed as kernel args. Every Float64
# literal/param is carried at the working element type so the kernel emits no
# Float64 on a Float32-only backend (Metal); on Float64 this is byte-identical.
@kernel function _ch4oxid_kernel!(ch4_oxid_depth, o2_oxid_depth,
                                  @Const(mask), @Const(jwt),
                                  @Const(watsat), @Const(h2osoi_vol),
                                  @Const(smp_l), @Const(t_soisno),
                                  @Const(conc_ch4), @Const(conc_o2),
                                  sat::Int,
                                  vmax_ch4_oxid, k_m, q10_ch4oxid, smp_crit,
                                  k_m_o2, k_m_unsat, vmax_oxid_unsat,
                                  c_h_inv1, c_h_inv2, kh_theta1, kh_theta2,
                                  kh_tbase, rgaslatm, tfrz, t0)
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        T = eltype(ch4_oxid_depth)

        if sat == 1 || j > jwt[c]
            k_m_eff = k_m
            vmax_eff = vmax_ch4_oxid
        else
            k_m_eff = k_m_unsat
            vmax_eff = vmax_oxid_unsat
        end

        porevol = smooth_max(watsat[c, j] - h2osoi_vol[c, j], zero(T))
        h2osoi_vol_min = smooth_min(watsat[c, j], h2osoi_vol[c, j])

        if j <= jwt[c] && smp_l[c, j] < zero(T)
            smp_fact = exp(-smp_l[c, j] / smp_crit)
        else
            smp_fact = one(T)
        end

        if j <= jwt[c]
            k_h_inv = exp(-c_h_inv1 * (one(T) / t_soisno[c, j] - one(T) / kh_tbase) + log(kh_theta1))
            k_h_cc = t_soisno[c, j] / k_h_inv * rgaslatm
            conc_ch4_rel = conc_ch4[c, j] / (h2osoi_vol_min + porevol / k_h_cc)

            k_h_inv = exp(-c_h_inv2 * (one(T) / t_soisno[c, j] - one(T) / kh_tbase) + log(kh_theta2))
            k_h_cc = t_soisno[c, j] / k_h_inv * rgaslatm
            conc_o2_rel = conc_o2[c, j] / (h2osoi_vol_min + porevol / k_h_cc)
        else
            conc_ch4_rel = conc_ch4[c, j] / watsat[c, j]
            conc_o2_rel = conc_o2[c, j] / watsat[c, j]
        end

        oxid_a = vmax_eff * h2osoi_vol_min * conc_ch4_rel / (k_m_eff + conc_ch4_rel) *
                 conc_o2_rel / (k_m_o2 + conc_o2_rel) *
                 q10_ch4oxid^((t_soisno[c, j] - t0) / T(10.0)) * smp_fact

        if t_soisno[c, j] <= tfrz
            oxid_a = zero(T)
        end

        ch4_oxid_depth[c, j] = oxid_a
        o2_oxid_depth[c, j] = ch4_oxid_depth[c, j] * T(2.0)
    end
end

# Loose-array launcher: extracts no struct fields, takes the sat/unsat arrays
# directly so the harness can drive the device path with adapted arrays.
function meth_oxid!(ch4_oxid_depth, o2_oxid_depth, mask, jwt,
                    watsat, h2osoi_vol, smp_l, t_soisno, conc_ch4, conc_o2,
                    sat::Int, vmax_ch4_oxid, k_m, q10_ch4oxid, smp_crit,
                    k_m_o2, k_m_unsat, vmax_oxid_unsat, nc::Int, nlevsoi::Int)
    FT = eltype(ch4_oxid_depth)
    t0 = FT(TFRZ) + FT(12.0)
    _launch!(_ch4oxid_kernel!, ch4_oxid_depth, o2_oxid_depth, mask, jwt,
             watsat, h2osoi_vol, smp_l, t_soisno, conc_ch4, conc_o2,
             sat, FT(vmax_ch4_oxid), FT(k_m), FT(q10_ch4oxid), FT(smp_crit),
             FT(k_m_o2), FT(k_m_unsat), FT(vmax_oxid_unsat),
             FT(C_H_INV[1]), FT(C_H_INV[2]), FT(KH_THETA[1]), FT(KH_THETA[2]),
             FT(KH_TBASE), FT(rgasLatm), FT(TFRZ), t0;
             ndrange = (nc, nlevsoi))
end

"""
    ch4_oxid!(ch4, params, mask_soil, watsat, h2osoi_vol, smp_l, t_soisno,
              jwt, sat, lake, nlevsoi, dtime)

Calculate CH4 oxidation in each soil layer using double Michaelis-Menten kinetics.
Ported from `ch4_oxid` in `ch4Mod.F90`. Kernelized: a single per-(column, layer)
KernelAbstractions kernel (`_ch4oxid_kernel!`); backend-agnostic (CPU loop or GPU).
"""
function ch4_oxid!(ch4::CH4Data,
                   params::CH4Params,
                   mask_soil::AbstractVector{Bool},
                   watsat::AbstractMatrix{<:Real},
                   h2osoi_vol::AbstractMatrix{<:Real},
                   smp_l::AbstractMatrix{<:Real},
                   t_soisno::AbstractMatrix{<:Real},
                   jwt::AbstractVector{<:Integer},
                   sat::Int,
                   lake::Bool,
                   nlevsoi::Int,
                   dtime::Real)

    if sat == 0
        ch4_oxid_depth = ch4.ch4_oxid_depth_unsat_col
        o2_oxid_depth = ch4.o2_oxid_depth_unsat_col
        conc_ch4 = ch4.conc_ch4_unsat_col
        conc_o2 = ch4.conc_o2_unsat_col
    else
        ch4_oxid_depth = ch4.ch4_oxid_depth_sat_col
        o2_oxid_depth = ch4.o2_oxid_depth_sat_col
        conc_ch4 = ch4.conc_ch4_sat_col
        conc_o2 = ch4.conc_o2_sat_col
    end

    nc = length(mask_soil)
    meth_oxid!(ch4_oxid_depth, o2_oxid_depth, mask_soil, jwt,
               watsat, h2osoi_vol, smp_l, t_soisno, conc_ch4, conc_o2,
               sat, params.vmax_ch4_oxid, params.k_m, params.q10_ch4oxid,
               params.smp_crit, params.k_m_o2, params.k_m_unsat,
               params.vmax_oxid_unsat, nc, nlevsoi)
    nothing
end

# ---------------------------------------------------------------------------
# site_ox_aere! — Site-level aerenchyma helper
# ---------------------------------------------------------------------------

"""
    site_ox_aere!(tranloss, aere, oxaere,
                  is_vegetated, watsat, h2osoi_vol, t_soisno,
                  conc_ch4, rootr, qflx_tran_veg, jwt,
                  annsum_npp, annavg_agnpp, annavg_bgnpp,
                  elai, poros_tiller, rootfr, grnd_ch4_cond,
                  conc_o2, c_atm, z, dz, sat, nlevsoi,
                  params, ch4vc)

Site-level aerenchyma fluxes (O2 in, CH4 out, transpiration loss).
Ported from `SiteOxAere` in `ch4Mod.F90`.
"""
function site_ox_aere!(tranloss::Vector{<:Real},
                       aere::Vector{<:Real},
                       oxaere::Vector{<:Real},
                       is_vegetated::Bool,
                       watsat::AbstractVector{Float64},
                       h2osoi_vol::AbstractVector{Float64},
                       t_soisno::AbstractVector{Float64},
                       conc_ch4::AbstractVector{Float64},
                       rootr::AbstractVector{Float64},
                       qflx_tran_veg::Real,
                       jwt::Int,
                       annsum_npp::Real,
                       annavg_agnpp::Real,
                       annavg_bgnpp::Real,
                       elai::Real,
                       poros_tiller::Real,
                       rootfr::AbstractVector{Float64},
                       grnd_ch4_cond::Real,
                       conc_o2::AbstractVector{Float64},
                       c_atm::AbstractVector{Float64},
                       z::AbstractVector{Float64},
                       dz::AbstractVector{Float64},
                       sat::Int,
                       nlevsoi::Int,
                       params::CH4Params,
                       ch4vc::CH4VarCon)

    smallnumber = 1.0e-12
    diffus_aere = D_CON_G[1, 1] * 1.0e-4  # for CH4: m^2/s

    for j in 1:nlevsoi
        # Transpiration loss
        if ch4vc.transpirationloss && is_vegetated
            h2osoi_vol_min = smooth_min(watsat[j], h2osoi_vol[j])
            k_h_inv = exp(-C_H_INV[1] * (1.0 / t_soisno[j] - 1.0 / KH_TBASE) + log(KH_THETA[1]))
            k_h_cc = t_soisno[j] / k_h_inv * rgasLatm
            conc_ch4_wat = conc_ch4[j] / ((watsat[j] - h2osoi_vol_min) / k_h_cc + h2osoi_vol_min)
            tranloss[j] = conc_ch4_wat * rootr[j] * qflx_tran_veg / dz[j] / 1000.0
            tranloss[j] = smooth_max(tranloss[j], 0.0)
        else
            tranloss[j] = 0.0
        end

        # Aerenchyma diffusion
        if j > jwt && t_soisno[j] > TFRZ && is_vegetated
            anpp = smooth_max(annsum_npp, 0.0)

            if annavg_agnpp != SPVAL && annavg_bgnpp != SPVAL &&
               annavg_agnpp > 0.0 && annavg_bgnpp > 0.0
                nppratio = annavg_bgnpp / (annavg_agnpp + annavg_bgnpp)
            else
                nppratio = 0.5
            end

            m_tiller = anpp * nppratio * 4.0
            n_tiller = m_tiller / 0.22

            pt = poros_tiller
            if sat == 0
                pt *= params.unsat_aere_ratio
            end
            pt = smooth_max(pt, params.porosmin)

            area_tiller = params.scale_factor_aere * n_tiller * pt * RPI * 2.9e-3^2

            k_h_inv = exp(-C_H_INV[1] * (1.0 / t_soisno[j] - 1.0 / KH_TBASE) + log(KH_THETA[1]))
            k_h_cc = t_soisno[j] / k_h_inv * rgasLatm
            aerecond = area_tiller * rootfr[j] * diffus_aere / (z[j] * params.rob)
            aerecond = 1.0 / (1.0 / (aerecond + smallnumber) + 1.0 / (grnd_ch4_cond + smallnumber))

            aere[j] = aerecond * (conc_ch4[j] / watsat[j] / k_h_cc - c_atm[1]) / dz[j]
            aere[j] = smooth_max(aere[j], 0.0)

            # O2 diffusion
            k_h_inv = exp(-C_H_INV[2] * (1.0 / t_soisno[j] - 1.0 / KH_TBASE) + log(KH_THETA[2]))
            k_h_cc = t_soisno[j] / k_h_inv * rgasLatm
            oxdiffus = diffus_aere * D_CON_G[2, 1] / D_CON_G[1, 1]
            aerecond = area_tiller * rootfr[j] * oxdiffus / (z[j] * params.rob)
            aerecond = 1.0 / (1.0 / (aerecond + smallnumber) + 1.0 / (grnd_ch4_cond + smallnumber))
            oxaere[j] = -aerecond * (conc_o2[j] / watsat[j] / k_h_cc - c_atm[2]) / dz[j]
            oxaere[j] = smooth_max(oxaere[j], 0.0)

            if !ch4vc.use_aereoxid_prog
                oxaere[j] = 0.0
            end
        else
            aere[j] = 0.0
            oxaere[j] = 0.0
        end
    end
    nothing
end

# ---------------------------------------------------------------------------
# ch4_aere! — Aerenchyma transport  (kernelized: per-(c,j) zero + per-patch scatter)
# ---------------------------------------------------------------------------

# Masked per-(column, level) zeroing of the three patch→column scatter targets.
@kernel function _ch4aere_zero_kernel!(ch4_aere_depth, @Const(mask),
                                       ch4_tran_depth, o2_aere_depth)
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        T = eltype(ch4_aere_depth)
        ch4_aere_depth[c, j] = zero(T)
        ch4_tran_depth[c, j] = zero(T)
        o2_aere_depth[c, j]  = zero(T)
    end
end

# Device-view bundle of the array fields the per-patch aerenchyma kernel touches.
# Field names mirror the ch4_aere! locals so the body reads verbatim. Separate
# index/flag vectors (Int/Bool) from float arrays so adapt keeps their eltype.
Base.@kwdef struct _CH4AereDV{M,V,VI,VB}
    # column [c,j] matrices (read)
    watsat::M; h2osoi_vol::M; t_soisno::M; conc_ch4::M; conc_o2::M; z::M; dz::M
    ch4_prod_depth::M
    # column [c,j] matrices (scatter targets)
    ch4_aere_depth::M; ch4_tran_depth::M; o2_aere_depth::M
    # patch [p,j] matrices
    rootr::M; rootfr::M
    # gridcell [g,2] matrix
    c_atm_grc::M
    # patch [p] vectors
    qflx_tran_veg::V; annsum_npp::V; annavg_agnpp_patch::V; annavg_bgnpp_patch::V
    grnd_ch4_cond_patch::V; patch_wtcol::V
    # index/flag vectors
    patch_column::VI; col_gridcell::VI; patch_itype::VI; jwt::VI; is_fates::VB
end
Adapt.@adapt_structure _CH4AereDV

# Scalar parameter/flag bundle (kept off the loose-arg list; Metal caps ~31 args).
Base.@kwdef struct _CH4AereP{T}
    smallnumber::T; diffus_aere::T; dcong21::T; dcong11::T; rgaslatm::T
    unsat_aere_ratio::T; porosmin::T; scale_factor_aere::T; rob::T; nongrassporosratio::T
    c_h_inv1::T; c_h_inv2::T; kh_theta1::T; kh_theta2::T; kh_tbase::T
    tfrz::T; spval::T; rpi::T; poros_base::T
    transpirationloss::Bool; use_aereoxid_prog::Bool
    sat::Int; noveg::Int; nlevsoi::Int
end
Adapt.@adapt_structure _CH4AereP

# Per-PATCH aerenchyma kernel. Runs site_ox_aere!'s math inline as a sequential
# in-thread j-loop (own per-level scalars, no scratch arrays), fuses the second
# j-loop, and _scatter_add!s the per-patch contribution into the column arrays.
@kernel function _ch4aere_patch_kernel!(_out, d, @Const(mask), p_par, dtime)
    p = @index(Global)
    @inbounds if mask[p]
        T = typeof(dtime)
        c = d.patch_column[p]
        g = d.col_gridcell[c]
        wt = d.patch_wtcol[p]

        # is_vegetated / poros_tiller: FATES columns force vegetated (default poros).
        if d.is_fates[c]
            is_vegetated = true
        else
            is_vegetated = d.patch_itype[p] != p_par.noveg
        end
        poros_tiller = p_par.poros_base * p_par.nongrassporosratio

        jwt_c = d.jwt[c]
        for j in 1:p_par.nlevsoi
            # --- Transpiration loss ---
            if p_par.transpirationloss && is_vegetated
                h2osoi_vol_min = smooth_min(d.watsat[c, j], d.h2osoi_vol[c, j])
                k_h_inv = exp(-p_par.c_h_inv1 * (one(T) / d.t_soisno[c, j] - one(T) / p_par.kh_tbase) + log(p_par.kh_theta1))
                k_h_cc = d.t_soisno[c, j] / k_h_inv * p_par.rgaslatm
                conc_ch4_wat = d.conc_ch4[c, j] / ((d.watsat[c, j] - h2osoi_vol_min) / k_h_cc + h2osoi_vol_min)
                tranloss = conc_ch4_wat * d.rootr[p, j] * d.qflx_tran_veg[p] / d.dz[c, j] / T(1000.0)
                tranloss = smooth_max(tranloss, zero(T))
            else
                tranloss = zero(T)
            end

            # --- Aerenchyma diffusion ---
            if j > jwt_c && d.t_soisno[c, j] > p_par.tfrz && is_vegetated
                anpp = smooth_max(d.annsum_npp[p], zero(T))

                if d.annavg_agnpp_patch[p] != p_par.spval && d.annavg_bgnpp_patch[p] != p_par.spval &&
                   d.annavg_agnpp_patch[p] > zero(T) && d.annavg_bgnpp_patch[p] > zero(T)
                    nppratio = d.annavg_bgnpp_patch[p] / (d.annavg_agnpp_patch[p] + d.annavg_bgnpp_patch[p])
                else
                    nppratio = T(0.5)
                end

                m_tiller = anpp * nppratio * T(4.0)
                n_tiller = m_tiller / T(0.22)

                pt = poros_tiller
                if p_par.sat == 0
                    pt *= p_par.unsat_aere_ratio
                end
                pt = smooth_max(pt, p_par.porosmin)

                area_tiller = p_par.scale_factor_aere * n_tiller * pt * p_par.rpi * T(2.9e-3)^2

                k_h_inv = exp(-p_par.c_h_inv1 * (one(T) / d.t_soisno[c, j] - one(T) / p_par.kh_tbase) + log(p_par.kh_theta1))
                k_h_cc = d.t_soisno[c, j] / k_h_inv * p_par.rgaslatm
                aerecond = area_tiller * d.rootfr[p, j] * p_par.diffus_aere / (d.z[c, j] * p_par.rob)
                aerecond = one(T) / (one(T) / (aerecond + p_par.smallnumber) + one(T) / (d.grnd_ch4_cond_patch[p] + p_par.smallnumber))

                aere = aerecond * (d.conc_ch4[c, j] / d.watsat[c, j] / k_h_cc - d.c_atm_grc[g, 1]) / d.dz[c, j]
                aere = smooth_max(aere, zero(T))

                # O2 diffusion
                k_h_inv = exp(-p_par.c_h_inv2 * (one(T) / d.t_soisno[c, j] - one(T) / p_par.kh_tbase) + log(p_par.kh_theta2))
                k_h_cc = d.t_soisno[c, j] / k_h_inv * p_par.rgaslatm
                oxdiffus = p_par.diffus_aere * p_par.dcong21 / p_par.dcong11
                aerecond = area_tiller * d.rootfr[p, j] * oxdiffus / (d.z[c, j] * p_par.rob)
                aerecond = one(T) / (one(T) / (aerecond + p_par.smallnumber) + one(T) / (d.grnd_ch4_cond_patch[p] + p_par.smallnumber))
                oxaere = -aerecond * (d.conc_o2[c, j] / d.watsat[c, j] / k_h_cc - d.c_atm_grc[g, 2]) / d.dz[c, j]
                oxaere = smooth_max(oxaere, zero(T))

                if !p_par.use_aereoxid_prog
                    oxaere = zero(T)
                end
            else
                aere = zero(T)
                oxaere = zero(T)
            end

            # --- Scatter the per-patch contribution into the column arrays ---
            aeretran = smooth_min(aere + tranloss, d.conc_ch4[c, j] / dtime + d.ch4_prod_depth[c, j])
            _scatter_add!(d.ch4_aere_depth, c, j, aeretran * wt)
            _scatter_add!(d.ch4_tran_depth, c, j, smooth_min(tranloss, aeretran) * wt)
            _scatter_add!(d.o2_aere_depth,  c, j, oxaere * wt)
        end
    end
end

"""
    ch4_aere!(ch4, params, ch4vc, mask_soil, mask_soilp,
              patch_column, patch_itype, patch_wtcol, is_fates,
              watsat, h2osoi_vol, t_soisno, rootr, rootfr, elai,
              qflx_tran_veg, grnd_ch4_cond, annsum_npp,
              z, dz, jwt, sat, lake, nlevsoi, dtime, noveg)

Calculate CH4 loss and O2 gain via aerenchyma transport.
Ported from `ch4_aere` in `ch4Mod.F90`.
"""
function ch4_aere!(ch4::CH4Data,
                   params::CH4Params,
                   ch4vc::CH4VarCon,
                   mask_soil::AbstractVector{Bool},
                   mask_soilp::AbstractVector{Bool},
                   patch_column::AbstractVector{<:Integer},
                   patch_itype::AbstractVector{<:Integer},
                   patch_wtcol::AbstractVector{<:Real},
                   col_gridcell::AbstractVector{<:Integer},
                   is_fates::AbstractVector{Bool},
                   watsat::AbstractMatrix{<:Real},
                   h2osoi_vol::AbstractMatrix{<:Real},
                   t_soisno::AbstractMatrix{<:Real},
                   rootr::AbstractMatrix{<:Real},
                   rootfr::AbstractMatrix{<:Real},
                   elai::AbstractVector{<:Real},
                   qflx_tran_veg::AbstractVector{<:Real},
                   annsum_npp::AbstractVector{<:Real},
                   z::AbstractMatrix{<:Real},
                   dz::AbstractMatrix{<:Real},
                   jwt::AbstractVector{<:Integer},
                   sat::Int,
                   lake::Bool,
                   nlevsoi::Int,
                   dtime::Real,
                   noveg::Int)

    if sat == 0
        ch4_aere_depth = ch4.ch4_aere_depth_unsat_col
        ch4_tran_depth = ch4.ch4_tran_depth_unsat_col
        o2_aere_depth = ch4.o2_aere_depth_unsat_col
        conc_ch4 = ch4.conc_ch4_unsat_col
        conc_o2 = ch4.conc_o2_unsat_col
        ch4_prod_depth = ch4.ch4_prod_depth_unsat_col
    else
        ch4_aere_depth = ch4.ch4_aere_depth_sat_col
        ch4_tran_depth = ch4.ch4_tran_depth_sat_col
        o2_aere_depth = ch4.o2_aere_depth_sat_col
        conc_ch4 = ch4.conc_ch4_sat_col
        conc_o2 = ch4.conc_o2_sat_col
        ch4_prod_depth = ch4.ch4_prod_depth_sat_col
    end

    nc = length(mask_soil)

    # Initialize (masked per-(column, level) zeroing kernel)
    _launch!(_ch4aere_zero_kernel!, ch4_aere_depth, mask_soil,
             ch4_tran_depth, o2_aere_depth; ndrange = (nc, nlevsoi))

    if !lake
        FT_ae = eltype(t_soisno)
        T = FT_ae

        d = _CH4AereDV(;
            watsat = watsat, h2osoi_vol = h2osoi_vol, t_soisno = t_soisno,
            conc_ch4 = conc_ch4, conc_o2 = conc_o2, z = z, dz = dz,
            ch4_prod_depth = ch4_prod_depth,
            ch4_aere_depth = ch4_aere_depth, ch4_tran_depth = ch4_tran_depth,
            o2_aere_depth = o2_aere_depth,
            rootr = rootr, rootfr = rootfr, c_atm_grc = ch4.c_atm_grc,
            qflx_tran_veg = qflx_tran_veg, annsum_npp = annsum_npp,
            annavg_agnpp_patch = ch4.annavg_agnpp_patch,
            annavg_bgnpp_patch = ch4.annavg_bgnpp_patch,
            grnd_ch4_cond_patch = ch4.grnd_ch4_cond_patch,
            patch_wtcol = patch_wtcol,
            patch_column = patch_column, col_gridcell = col_gridcell,
            patch_itype = patch_itype, jwt = jwt, is_fates = is_fates)

        p_par = _CH4AereP(;
            smallnumber = T(1.0e-12),
            diffus_aere = T(D_CON_G[1, 1] * 1.0e-4),
            dcong21 = T(D_CON_G[2, 1]), dcong11 = T(D_CON_G[1, 1]),
            rgaslatm = T(rgasLatm),
            unsat_aere_ratio = T(params.unsat_aere_ratio),
            porosmin = T(params.porosmin),
            scale_factor_aere = T(params.scale_factor_aere),
            rob = T(params.rob),
            nongrassporosratio = T(params.nongrassporosratio),
            c_h_inv1 = T(C_H_INV[1]), c_h_inv2 = T(C_H_INV[2]),
            kh_theta1 = T(KH_THETA[1]), kh_theta2 = T(KH_THETA[2]),
            kh_tbase = T(KH_TBASE), tfrz = T(TFRZ), spval = T(SPVAL),
            rpi = T(RPI), poros_base = T(0.3),
            transpirationloss = ch4vc.transpirationloss,
            use_aereoxid_prog = ch4vc.use_aereoxid_prog,
            sat = sat, noveg = noveg, nlevsoi = nlevsoi)

        # Per-patch kernel: site_ox_aere! math inline (sequential j-loop) +
        # patch→column scatter into the (zeroed) column arrays. The first arg
        # `ch4_aere_depth` is the backend carrier (also reachable via `d`).
        np = length(mask_soilp)
        _launch!(_ch4aere_patch_kernel!, ch4_aere_depth, d, mask_soilp, p_par, T(dtime);
                 ndrange = np)
    end
    nothing
end

# ---------------------------------------------------------------------------
# ch4_ebul! — Ebullition
# ---------------------------------------------------------------------------

"""
    ch4_ebul!(ch4, params, mask_soil, watsat, h2osoi_vol, t_soisno,
              forc_pbot, h2osfc, frac_h2osfc, lake_icefrac, lakedepth,
              z, dz, zi, jwt, sat, lake, nlevsoi, dtime)

Calculate CH4 ebullition losses from saturated soil layers.
Ported from `ch4_ebul` in `ch4Mod.F90`.
"""
function ch4_ebul!(ch4::CH4Data,
                   params::CH4Params,
                   mask_soil::BitVector,
                   watsat::Matrix{<:Real},
                   h2osoi_vol::Matrix{<:Real},
                   t_soisno::Matrix{<:Real},
                   forc_pbot::Vector{<:Real},
                   h2osfc::Vector{<:Real},
                   frac_h2osfc::Vector{<:Real},
                   lake_icefrac::Matrix{<:Real},
                   lakedepth::Vector{<:Real},
                   z::Matrix{<:Real},
                   dz::Matrix{<:Real},
                   zi::AbstractMatrix{<:Real},
                   jwt::AbstractVector{<:Integer},
                   sat::Int,
                   lake::Bool,
                   nlevsoi::Int,
                   dtime::Real)

    if sat == 0
        ch4_ebul_depth = ch4.ch4_ebul_depth_unsat_col
        conc_ch4 = ch4.conc_ch4_unsat_col
    else
        ch4_ebul_depth = ch4.ch4_ebul_depth_sat_col
        conc_ch4 = ch4.conc_ch4_sat_col
    end

    vgc_max = params.vgc_max
    bubble_f = 0.57
    vgc_min = vgc_max
    ebul_timescale = dtime
    rgasm = RGAS / 1000.0

    nc = length(mask_soil)
    meth_ebul!(ch4_ebul_depth, mask_soil, jwt, t_soisno, conc_ch4, watsat,
               forc_pbot, h2osfc, frac_h2osfc, lake_icefrac, lakedepth, z, zi,
               sat, lake, vgc_max, vgc_min, bubble_f, ebul_timescale, rgasm, rgasLatm,
               nc, nlevsoi)
    nothing
end

# ---------------------------------------------------------------------------
# ch4_tran! — Reaction-diffusion transport solver
# ---------------------------------------------------------------------------

"""
    ch4_tran!(ch4, params, ch4vc, mask_soil, col_gridcell,
              watsat, h2osoi_vol, h2osoi_liq, h2osoi_ice, h2osfc,
              bsw, cellorg, smp_l, t_soisno, t_grnd, t_h2osfc,
              frac_h2osfc, snow_depth, snl,
              z, dz, zi, jwt, sat, lake,
              nlevsoi, nlevsno, dtime_ch4,
              organic_max, c_atm_local)

Solve reaction-diffusion equation for CH4 and O2 transport.
Ported from `ch4_tran` in `ch4Mod.F90`.
"""

# ===========================================================================
# ch4_tran! GPU kernelization. Every loop is per-(column,level) or per-column and
# the tridiagonal solves are per-column, so the whole routine is one big per-column
# kernel: thread c runs competition, the ebullition/aerenchyma reductions, the
# 2-species source/epsilon setup, then for each species s=1,2 builds the Patankar
# tridiagonal system and runs an IN-THREAD Thomas solve (inlined to match
# tridiagonal_solve! exactly), then the post-solve clamps + balance check. All
# writes target the thread's own column row -> race-free, byte-identical. Scratch is
# grouped into device-view bundles; meth_khcc! (existing 2D kernel) fills k_h_cc first.
# ===========================================================================
Base.@kwdef struct _Ch4TS{V,M}   # selected sat/unsat CH4 state arrays
    o2_decomp_depth::M; o2stress::M; ch4_oxid_depth::M; ch4_prod_depth::M
    ch4_aere_depth::M; ch4_ebul_depth::M; o2_oxid_depth::M; o2_aere_depth::M
    ch4stress::M; conc_ch4::M; conc_o2::M
    ch4_ebul_total::V; ch4_surf_aere::V; ch4_surf_ebul::V; ch4_surf_diff::V
end
Adapt.@adapt_structure _Ch4TS
Base.@kwdef struct _Ch4TF{V,M}   # forcing
    watsat::M; h2osoi_vol::M; h2osoi_liq::M; h2osoi_ice::M; bsw::M; cellorg::M
    t_soisno::M; dz::M
    h2osfc::V; t_grnd::V; t_h2osfc::V; frac_h2osfc::V
end
Adapt.@adapt_structure _Ch4TF
Base.@kwdef struct _Ch4Scr{M,A3}   # per-column scratch
    epsilon_t::A3; source::A3; k_h_cc::A3
    conc_ch4_bef::M; conc_ch4_rel::M; conc_o2_rel::M; conc_ch4_rel_old::M; conc_work::M
    h2osoi_vol_min::M; liqfrac::M; diffus::M; dp1_zp1::M; dm1_zm1::M
    at::M; bt::M; ct::M; rt::M; spec_grnd_cond::M; cp::M; dp::M
end
Adapt.@adapt_structure _Ch4Scr
Base.@kwdef struct _Ch4TP{T}   # isbits scalar params (at working precision)
    satpow::T; sfgas::T; sfliq::T; capthick::T; aereoxid::T; om_frac_sf::T
    dtime::T; organic_max::T; smallnumber::T
end

@kernel function _ch4tran_column_kernel!(ts, tf, scr, @Const(c_atm), @Const(jwt),
        @Const(col_gridcell), grnd_ch4_cond, @Const(dcong), @Const(dconw), @Const(mask),
        p::_Ch4TP, ch4frzout::Bool, use_aereoxid_prog::Bool, lake::Bool, sat::Int, nlevsoi::Int)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(scr.diffus)
        dtime = p.dtime
        sn = p.smallnumber
        g = col_gridcell[c]

        # --- Competition for O2 and CH4 ---
        for j in 1:nlevsoi
            o2demand = ts.o2_decomp_depth[c, j] + ts.o2_oxid_depth[c, j]
            if o2demand > zero(T)
                if (ts.conc_o2[c, j] / dtime + ts.o2_aere_depth[c, j]) > o2demand
                    ts.o2stress[c, j] = one(T)
                else
                    ts.o2stress[c, j] = (ts.conc_o2[c, j] / dtime + ts.o2_aere_depth[c, j]) / o2demand
                end
            else
                ts.o2stress[c, j] = one(T)
            end
            ch4demand = ts.ch4_oxid_depth[c, j] + ts.ch4_aere_depth[c, j] + ts.ch4_ebul_depth[c, j]
            if ch4demand > zero(T)
                ts.ch4stress[c, j] = smooth_min((ts.conc_ch4[c, j] / dtime + ts.ch4_prod_depth[c, j]) / ch4demand, one(T))
            else
                ts.ch4stress[c, j] = one(T)
            end
            if ts.o2stress[c, j] < one(T) || ts.ch4stress[c, j] < one(T)
                if ts.ch4stress[c, j] <= ts.o2stress[c, j]
                    if ts.o2stress[c, j] < one(T)
                        o2demand2 = ts.o2_decomp_depth[c, j]
                        if o2demand2 > zero(T)
                            ts.o2stress[c, j] = smooth_min((ts.conc_o2[c, j] / dtime + ts.o2_aere_depth[c, j] -
                                                   ts.ch4stress[c, j] * ts.o2_oxid_depth[c, j]) / o2demand2, one(T))
                        else
                            ts.o2stress[c, j] = one(T)
                        end
                    end
                    ts.ch4_oxid_depth[c, j] *= ts.ch4stress[c, j]
                    ts.o2_oxid_depth[c, j] *= ts.ch4stress[c, j]
                else
                    if ts.ch4stress[c, j] < one(T)
                        ch4demand2 = ts.ch4_aere_depth[c, j] + ts.ch4_ebul_depth[c, j]
                        if ch4demand2 > zero(T)
                            ts.ch4stress[c, j] = smooth_min((ts.conc_ch4[c, j] / dtime + ts.ch4_prod_depth[c, j] -
                                                    ts.o2stress[c, j] * ts.ch4_oxid_depth[c, j]) / ch4demand2, one(T))
                        else
                            ts.ch4stress[c, j] = one(T)
                        end
                    end
                    ts.ch4_oxid_depth[c, j] *= ts.o2stress[c, j]
                    ts.o2_oxid_depth[c, j] *= ts.o2stress[c, j]
                end
            end
            ts.ch4_aere_depth[c, j] *= ts.ch4stress[c, j]
            ts.ch4_ebul_depth[c, j] *= ts.ch4stress[c, j]
            ts.o2_decomp_depth[c, j] *= ts.o2stress[c, j]
        end

        # --- Accumulate ebullition (per-column reduction) ---
        ts.ch4_ebul_total[c] = zero(T)
        for j in 1:nlevsoi
            ts.ch4_ebul_total[c] += ts.ch4_ebul_depth[c, j] * tf.dz[c, j]
        end

        # --- Source terms + epsilon_t (k_h_cc already filled) ---
        for j in 1:nlevsoi
            scr.h2osoi_vol_min[c, j] = smooth_min(tf.watsat[c, j], tf.h2osoi_vol[c, j])
            if ch4frzout
                scr.liqfrac[c, j] = smooth_max(T(0.05), (tf.h2osoi_liq[c, j] / T(DENH2O) + sn) /
                                           (tf.h2osoi_liq[c, j] / T(DENH2O) + tf.h2osoi_ice[c, j] / T(DENICE) + sn))
            else
                scr.liqfrac[c, j] = one(T)
            end
            if j <= jwt[c]
                for s in 1:2
                    scr.epsilon_t[c, j, s] = tf.watsat[c, j] - (one(T) - scr.k_h_cc[c, j+1, s]) * scr.h2osoi_vol_min[c, j] * scr.liqfrac[c, j]
                end
            else
                for s in 1:2
                    scr.epsilon_t[c, j, s] = tf.watsat[c, j] * scr.liqfrac[c, j]
                end
            end
            if !use_aereoxid_prog
                ts.ch4_oxid_depth[c, j] += p.aereoxid * ts.ch4_aere_depth[c, j]
                ts.ch4_aere_depth[c, j] -= p.aereoxid * ts.ch4_aere_depth[c, j]
            end
            scr.source[c, j, 1] = ts.ch4_prod_depth[c, j] - ts.ch4_oxid_depth[c, j] -
                                  ts.ch4_aere_depth[c, j] - ts.ch4_ebul_depth[c, j]
            scr.source[c, j, 2] = -ts.o2_oxid_depth[c, j] - ts.o2_decomp_depth[c, j] + ts.o2_aere_depth[c, j]
            scr.conc_ch4_bef[c, j] = ts.conc_ch4[c, j]
        end

        # --- Accumulate aerenchyma surface flux (reduction) ---
        ts.ch4_surf_aere[c] = zero(T)
        for j in 1:nlevsoi
            ts.ch4_surf_aere[c] += ts.ch4_aere_depth[c, j] * tf.dz[c, j]
        end

        # --- Add ebullition to source at the water-table layer ---
        if jwt[c] != 0
            scr.source[c, jwt[c], 1] += ts.ch4_ebul_total[c] / tf.dz[c, jwt[c]]
        end

        # --- Relative concentrations (j=0 boundary -> atmosphere) ---
        for j in 0:nlevsoi
            if j == 0
                scr.conc_ch4_rel[c, 1] = c_atm[g, 1]
                scr.conc_o2_rel[c, 1] = c_atm[g, 2]
            else
                scr.conc_ch4_rel[c, j+1] = ts.conc_ch4[c, j] / scr.epsilon_t[c, j, 1]
                scr.conc_o2_rel[c, j+1] = ts.conc_o2[c, j] / scr.epsilon_t[c, j, 2]
            end
        end
        for jj in 1:(nlevsoi + 1)
            scr.conc_ch4_rel_old[c, jj] = scr.conc_ch4_rel[c, jj]
        end

        nlevs = nlevsoi + 1
        for s in 1:2
            if s == 1
                for jj in 1:nlevs; scr.conc_work[c, jj] = scr.conc_ch4_rel[c, jj]; end
            else
                for jj in 1:nlevs; scr.conc_work[c, jj] = scr.conc_o2_rel[c, jj]; end
            end

            # Snow/pond resistance + ground conductance
            if grnd_ch4_cond[c] < sn && s == 1
                grnd_ch4_cond[c] = sn
            end
            snowres = zero(T)
            pondres = zero(T)
            if !lake && sat == 1 && tf.frac_h2osfc[c] > zero(T)
                if tf.t_h2osfc[c] >= T(TFRZ)
                    tsc = tf.t_h2osfc[c] - T(TFRZ)
                    ponddiff = (dconw[s, 1] + dconw[s, 2] * tsc + dconw[s, 3] * tsc^2) * T(1.0e-9) * p.sfliq
                    pondz = tf.h2osfc[c] / T(1000.0) / tf.frac_h2osfc[c]
                    pondres = pondz / ponddiff
                elseif tf.h2osfc[c] / tf.frac_h2osfc[c] > p.capthick
                    pondres = one(T) / sn
                end
            end
            scr.spec_grnd_cond[c, s] = one(T) / (one(T) / grnd_ch4_cond[c] + snowres + pondres)

            # Gas/liquid diffusivity
            for j in 1:nlevsoi
                tsc = tf.t_soisno[c, j] - T(TFRZ)
                if j <= jwt[c]
                    f_a = one(T) - scr.h2osoi_vol_min[c, j] / tf.watsat[c, j]
                    eps = tf.watsat[c, j] - scr.h2osoi_vol_min[c, j]
                    if p.organic_max > zero(T)
                        om_frac = smooth_min(p.om_frac_sf * tf.cellorg[c, j] / p.organic_max, one(T))
                    else
                        om_frac = one(T)
                    end
                    scr.diffus[c, j] = (dcong[s, 1] + dcong[s, 2] * tsc) * T(1.0e-4) *
                                   (om_frac * f_a^(T(10.0) / T(3.0)) / tf.watsat[c, j]^2 +
                                    (one(T) - om_frac) * eps^2 * f_a^(T(3.0) / tf.bsw[c, j])) * p.sfgas
                else
                    eps = tf.watsat[c, j]
                    scr.diffus[c, j] = eps^p.satpow * (dconw[s, 1] + dconw[s, 2] * tsc + dconw[s, 3] * tsc^2) * T(1.0e-9) * p.sfliq
                    if tf.t_soisno[c, j] <= T(TFRZ)
                        scr.diffus[c, j] *= (tf.h2osoi_liq[c, j] / T(DENH2O) + sn) /
                                        (tf.h2osoi_liq[c, j] / T(DENH2O) + tf.h2osoi_ice[c, j] / T(DENICE) + sn)
                    end
                end
                scr.diffus[c, j] = smooth_max(scr.diffus[c, j], sn)
            end

            # Tridiagonal coefficients dm1_zm1 / dp1_zp1 (jwt-branchy)
            for j in 1:nlevsoi
                if j == 1 && j != jwt[c] && j != jwt[c] + 1
                    scr.dm1_zm1[c, j] = one(T) / (one(T) / scr.spec_grnd_cond[c, s] + tf.dz[c, j] / (scr.diffus[c, j] * T(2.0)))
                    scr.dp1_zp1[c, j] = j < nlevsoi ? T(2.0) / (tf.dz[c, j] / scr.diffus[c, j] + tf.dz[c, j+1] / scr.diffus[c, j+1]) : zero(T)
                elseif j == 1 && j == jwt[c]
                    scr.dm1_zm1[c, j] = one(T) / (one(T) / scr.spec_grnd_cond[c, s] + tf.dz[c, j] / (scr.diffus[c, j] * T(2.0)))
                    scr.dp1_zp1[c, j] = j < nlevsoi ? T(2.0) / (tf.dz[c, j] * scr.k_h_cc[c, j+1, s] / scr.diffus[c, j] + tf.dz[c, j+1] / scr.diffus[c, j+1]) : zero(T)
                elseif j == 1
                    scr.dm1_zm1[c, j] = one(T) / (scr.k_h_cc[c, j, s] / scr.spec_grnd_cond[c, s] + tf.dz[c, j] / (scr.diffus[c, j] * T(2.0)))
                    scr.dp1_zp1[c, j] = j < nlevsoi ? T(2.0) / (tf.dz[c, j] / scr.diffus[c, j] + tf.dz[c, j+1] / scr.diffus[c, j+1]) : zero(T)
                elseif j < nlevsoi && j != jwt[c] && j != jwt[c] + 1
                    scr.dm1_zm1[c, j] = T(2.0) / (tf.dz[c, j] / scr.diffus[c, j] + tf.dz[c, j-1] / scr.diffus[c, j-1])
                    scr.dp1_zp1[c, j] = T(2.0) / (tf.dz[c, j] / scr.diffus[c, j] + tf.dz[c, j+1] / scr.diffus[c, j+1])
                elseif j < nlevsoi && j == jwt[c]
                    scr.dm1_zm1[c, j] = T(2.0) / (tf.dz[c, j] / scr.diffus[c, j] + tf.dz[c, j-1] / scr.diffus[c, j-1])
                    scr.dp1_zp1[c, j] = T(2.0) / (tf.dz[c, j] * scr.k_h_cc[c, j+1, s] / scr.diffus[c, j] + tf.dz[c, j+1] / scr.diffus[c, j+1])
                elseif j < nlevsoi
                    scr.dm1_zm1[c, j] = T(2.0) / (tf.dz[c, j] / scr.diffus[c, j] + tf.dz[c, j-1] * scr.k_h_cc[c, j, s] / scr.diffus[c, j-1])
                    scr.dp1_zp1[c, j] = T(2.0) / (tf.dz[c, j] / scr.diffus[c, j] + tf.dz[c, j+1] / scr.diffus[c, j+1])
                elseif j != jwt[c] + 1
                    scr.dm1_zm1[c, j] = T(2.0) / (tf.dz[c, j] / scr.diffus[c, j] + tf.dz[c, j-1] / scr.diffus[c, j-1])
                else
                    scr.dm1_zm1[c, j] = T(2.0) / (tf.dz[c, j] / scr.diffus[c, j] + tf.dz[c, j-1] * scr.k_h_cc[c, j, s] / scr.diffus[c, j-1])
                end
            end

            # Build tridiagonal system (j=0 -> index 1)
            for j in 0:nlevsoi
                jj = j + 1
                if j == 0
                    scr.at[c, jj] = zero(T); scr.bt[c, jj] = one(T); scr.ct[c, jj] = zero(T)
                    scr.rt[c, jj] = c_atm[g, s]
                elseif j < nlevsoi && j == jwt[c]
                    dzj = tf.dz[c, j]
                    scr.at[c, jj] = -T(0.5) / dzj * scr.dm1_zm1[c, j]
                    scr.bt[c, jj] = scr.epsilon_t[c, j, s] / dtime + T(0.5) / dzj * (scr.dp1_zp1[c, j] * scr.k_h_cc[c, j+1, s] + scr.dm1_zm1[c, j])
                    scr.ct[c, jj] = -T(0.5) / dzj * scr.dp1_zp1[c, j]
                    scr.rt[c, jj] = scr.epsilon_t[c, j, s] / dtime * scr.conc_work[c, jj] +
                                    T(0.5) / dzj * (scr.dp1_zp1[c, j] * (scr.conc_work[c, jj+1] - scr.conc_work[c, jj] * scr.k_h_cc[c, j+1, s]) -
                                                 scr.dm1_zm1[c, j] * (scr.conc_work[c, jj] - scr.conc_work[c, jj-1])) + scr.source[c, j, s]
                elseif j < nlevsoi && j == jwt[c] + 1
                    dzj = tf.dz[c, j]
                    scr.at[c, jj] = -T(0.5) / dzj * scr.dm1_zm1[c, j] * scr.k_h_cc[c, j, s]
                    scr.bt[c, jj] = scr.epsilon_t[c, j, s] / dtime + T(0.5) / dzj * (scr.dp1_zp1[c, j] + scr.dm1_zm1[c, j])
                    scr.ct[c, jj] = -T(0.5) / dzj * scr.dp1_zp1[c, j]
                    scr.rt[c, jj] = scr.epsilon_t[c, j, s] / dtime * scr.conc_work[c, jj] +
                                    T(0.5) / dzj * (scr.dp1_zp1[c, j] * (scr.conc_work[c, jj+1] - scr.conc_work[c, jj]) -
                                                 scr.dm1_zm1[c, j] * (scr.conc_work[c, jj] - scr.conc_work[c, jj-1] * scr.k_h_cc[c, j, s])) + scr.source[c, j, s]
                elseif j < nlevsoi
                    dzj = tf.dz[c, j]
                    scr.at[c, jj] = -T(0.5) / dzj * scr.dm1_zm1[c, j]
                    scr.bt[c, jj] = scr.epsilon_t[c, j, s] / dtime + T(0.5) / dzj * (scr.dp1_zp1[c, j] + scr.dm1_zm1[c, j])
                    scr.ct[c, jj] = -T(0.5) / dzj * scr.dp1_zp1[c, j]
                    scr.rt[c, jj] = scr.epsilon_t[c, j, s] / dtime * scr.conc_work[c, jj] +
                                    T(0.5) / dzj * (scr.dp1_zp1[c, j] * (scr.conc_work[c, jj+1] - scr.conc_work[c, jj]) -
                                                 scr.dm1_zm1[c, j] * (scr.conc_work[c, jj] - scr.conc_work[c, jj-1])) + scr.source[c, j, s]
                elseif j == nlevsoi && j == jwt[c] + 1
                    dzj = tf.dz[c, j]
                    scr.at[c, jj] = -T(0.5) / dzj * scr.dm1_zm1[c, j] * scr.k_h_cc[c, j, s]
                    scr.bt[c, jj] = scr.epsilon_t[c, j, s] / dtime + T(0.5) / dzj * scr.dm1_zm1[c, j]
                    scr.ct[c, jj] = zero(T)
                    scr.rt[c, jj] = scr.epsilon_t[c, j, s] / dtime * scr.conc_work[c, jj] +
                                    T(0.5) / dzj * (-scr.dm1_zm1[c, j] * (scr.conc_work[c, jj] - scr.conc_work[c, jj-1] * scr.k_h_cc[c, j, s])) + scr.source[c, j, s]
                else
                    dzj = tf.dz[c, j]
                    scr.at[c, jj] = -T(0.5) / dzj * scr.dm1_zm1[c, j]
                    scr.bt[c, jj] = scr.epsilon_t[c, j, s] / dtime + T(0.5) / dzj * scr.dm1_zm1[c, j]
                    scr.ct[c, jj] = zero(T)
                    scr.rt[c, jj] = scr.epsilon_t[c, j, s] / dtime * scr.conc_work[c, jj] +
                                    T(0.5) / dzj * (-scr.dm1_zm1[c, j] * (scr.conc_work[c, jj] - scr.conc_work[c, jj-1])) + scr.source[c, j, s]
                end
            end

            # In-thread Thomas solve (jtop=1), matching tridiagonal_solve! exactly.
            scr.cp[c, 1] = scr.ct[c, 1] / scr.bt[c, 1]
            scr.dp[c, 1] = scr.rt[c, 1] / scr.bt[c, 1]
            for jj in 2:nlevs
                denom = scr.bt[c, jj] - scr.at[c, jj] * scr.cp[c, jj-1]
                scr.cp[c, jj] = scr.ct[c, jj] / denom
                scr.dp[c, jj] = (scr.rt[c, jj] - scr.at[c, jj] * scr.dp[c, jj-1]) / denom
            end
            scr.conc_work[c, nlevs] = scr.dp[c, nlevs]
            for jj in (nlevs-1):-1:1
                scr.conc_work[c, jj] = scr.dp[c, jj] - scr.cp[c, jj] * scr.conc_work[c, jj+1]
            end

            if s == 1
                if jwt[c] != 0
                    ts.ch4_surf_diff[c] = scr.dm1_zm1[c, 1] * ((scr.conc_work[c, 2] + scr.conc_ch4_rel_old[c, 2]) / T(2.0) - c_atm[g, s])
                    ts.ch4_surf_ebul[c] = zero(T)
                else
                    ts.ch4_surf_diff[c] = scr.dm1_zm1[c, 1] * ((scr.conc_work[c, 2] + scr.conc_ch4_rel_old[c, 2]) / T(2.0) - c_atm[g, s] * scr.k_h_cc[c, 1, s])
                    ts.ch4_surf_ebul[c] = ts.ch4_ebul_total[c]
                end
                for j in 1:nlevsoi
                    jj = j + 1
                    if scr.conc_work[c, jj] < zero(T)
                        deficit = -scr.conc_work[c, jj] * scr.epsilon_t[c, j, 1] * tf.dz[c, j]
                        scr.conc_work[c, jj] = zero(T)
                        ts.ch4_surf_diff[c] -= deficit / dtime
                    end
                end
                for jj in 1:nlevs; scr.conc_ch4_rel[c, jj] = scr.conc_work[c, jj]; end
            else
                for j in 1:nlevsoi
                    jj = j + 1
                    scr.conc_work[c, jj] = smooth_max(scr.conc_work[c, jj], T(1.0e-12))
                    scr.conc_work[c, jj] = smooth_min(scr.conc_work[c, jj], c_atm[g, 2] / scr.epsilon_t[c, j, 2])
                end
                for jj in 1:nlevs; scr.conc_o2_rel[c, jj] = scr.conc_work[c, jj]; end
            end
        end  # species loop

        # Update absolute concentrations
        for j in 1:nlevsoi
            ts.conc_ch4[c, j] = scr.conc_ch4_rel[c, j+1] * scr.epsilon_t[c, j, 1]
            ts.conc_o2[c, j] = scr.conc_o2_rel[c, j+1] * scr.epsilon_t[c, j, 2]
        end

        # Balance check
        errch4 = zero(T)
        for j in 1:nlevsoi
            errch4 += (ts.conc_ch4[c, j] - scr.conc_ch4_bef[c, j]) * tf.dz[c, j]
            errch4 -= ts.ch4_prod_depth[c, j] * tf.dz[c, j] * dtime
            errch4 += ts.ch4_oxid_depth[c, j] * tf.dz[c, j] * dtime
        end
        errch4 += (ts.ch4_surf_aere[c] + ts.ch4_surf_ebul[c] + ts.ch4_surf_diff[c]) * dtime
        if abs(errch4) < T(1.0e-8)
            ts.ch4_surf_diff[c] -= errch4 / dtime
        end
        grnd_ch4_cond[c] = scr.spec_grnd_cond[c, 1]
    end
end


function ch4_tran!(ch4::CH4Data,
                   params::CH4Params,
                   ch4vc::CH4VarCon,
                   mask_soil::AbstractVector{Bool},
                   col_gridcell::AbstractVector{<:Integer},
                   watsat::AbstractMatrix{<:Real},
                   h2osoi_vol::AbstractMatrix{<:Real},
                   h2osoi_liq::AbstractMatrix{<:Real},
                   h2osoi_ice::AbstractMatrix{<:Real},
                   h2osfc::AbstractVector{<:Real},
                   bsw::AbstractMatrix{<:Real},
                   cellorg::AbstractMatrix{<:Real},
                   t_soisno::AbstractMatrix{<:Real},
                   t_grnd::AbstractVector{<:Real},
                   t_h2osfc::AbstractVector{<:Real},
                   frac_h2osfc::AbstractVector{<:Real},
                   snow_depth::AbstractVector{<:Real},
                   snl::AbstractVector{<:Integer},
                   z::AbstractMatrix{<:Real},
                   dz::AbstractMatrix{<:Real},
                   zi::AbstractMatrix{<:Real},
                   jwt::AbstractVector{<:Integer},
                   sat::Int,
                   lake::Bool,
                   nlevsoi::Int,
                   nlevsno::Int,
                   dtime_ch4::Real,
                   organic_max::Real)

    dtime = dtime_ch4
    FT = eltype(t_soisno)
    nc = length(mask_soil)

    # Select the sat/unsat CH4 arrays into the device-view state bundle.
    ts = sat == 0 ?
        _Ch4TS(; o2_decomp_depth=ch4.o2_decomp_depth_unsat_col, o2stress=ch4.o2stress_unsat_col,
                 ch4_oxid_depth=ch4.ch4_oxid_depth_unsat_col, ch4_prod_depth=ch4.ch4_prod_depth_unsat_col,
                 ch4_aere_depth=ch4.ch4_aere_depth_unsat_col, ch4_ebul_depth=ch4.ch4_ebul_depth_unsat_col,
                 o2_oxid_depth=ch4.o2_oxid_depth_unsat_col, o2_aere_depth=ch4.o2_aere_depth_unsat_col,
                 ch4stress=ch4.ch4stress_unsat_col, conc_ch4=ch4.conc_ch4_unsat_col, conc_o2=ch4.conc_o2_unsat_col,
                 ch4_ebul_total=ch4.ch4_ebul_total_unsat_col, ch4_surf_aere=ch4.ch4_surf_aere_unsat_col,
                 ch4_surf_ebul=ch4.ch4_surf_ebul_unsat_col, ch4_surf_diff=ch4.ch4_surf_diff_unsat_col) :
        _Ch4TS(; o2_decomp_depth=ch4.o2_decomp_depth_sat_col, o2stress=ch4.o2stress_sat_col,
                 ch4_oxid_depth=ch4.ch4_oxid_depth_sat_col, ch4_prod_depth=ch4.ch4_prod_depth_sat_col,
                 ch4_aere_depth=ch4.ch4_aere_depth_sat_col, ch4_ebul_depth=ch4.ch4_ebul_depth_sat_col,
                 o2_oxid_depth=ch4.o2_oxid_depth_sat_col, o2_aere_depth=ch4.o2_aere_depth_sat_col,
                 ch4stress=ch4.ch4stress_sat_col, conc_ch4=ch4.conc_ch4_sat_col, conc_o2=ch4.conc_o2_sat_col,
                 ch4_ebul_total=ch4.ch4_ebul_total_sat_col, ch4_surf_aere=ch4.ch4_surf_aere_sat_col,
                 ch4_surf_ebul=ch4.ch4_surf_ebul_sat_col, ch4_surf_diff=ch4.ch4_surf_diff_sat_col)

    tf = _Ch4TF(; watsat, h2osoi_vol, h2osoi_liq, h2osoi_ice, bsw, cellorg, t_soisno, dz,
                  h2osfc, t_grnd, t_h2osfc, frac_h2osfc)

    ref = t_soisno
    _m(w) = fill!(similar(ref, FT, nc, w), zero(FT))
    _a3(w) = fill!(similar(ref, FT, nc, w, 2), zero(FT))
    scr = _Ch4Scr(; epsilon_t=_a3(nlevsoi), source=_a3(nlevsoi), k_h_cc=_a3(nlevsoi + 1),
                    conc_ch4_bef=_m(nlevsoi), conc_ch4_rel=_m(nlevsoi + 1), conc_o2_rel=_m(nlevsoi + 1),
                    conc_ch4_rel_old=_m(nlevsoi + 1), conc_work=_m(nlevsoi + 1),
                    h2osoi_vol_min=_m(nlevsoi), liqfrac=fill!(similar(ref, FT, nc, nlevsoi), one(FT)),
                    diffus=_m(nlevsoi), dp1_zp1=_m(nlevsoi), dm1_zm1=_m(nlevsoi),
                    at=_m(nlevsoi + 1), bt=_m(nlevsoi + 1), ct=_m(nlevsoi + 1), rt=_m(nlevsoi + 1),
                    spec_grnd_cond=_m(2), cp=_m(nlevsoi + 1), dp=_m(nlevsoi + 1))

    # mask + integer index vectors onto the state backend (BitVector/host non-bitstype on device).
    mask = similar(ref, Bool, length(mask_soil)); copyto!(mask, collect(Bool, mask_soil))
    jwt_d = similar(ref, Int, length(jwt)); copyto!(jwt_d, jwt)
    gc_d = similar(ref, Int, length(col_gridcell)); copyto!(gc_d, col_gridcell)
    _devm(v) = (d = similar(ref, FT, size(v)...); copyto!(d, FT.(v)); d)
    c_atm_d = _devm(ch4.c_atm_grc)
    dcong = _devm(D_CON_G); dconw = _devm(D_CON_W)

    # Henry's-law coefficients (existing 2D kernel) -> fills scr.k_h_cc.
    meth_khcc!(scr.k_h_cc, mask, t_grnd, t_soisno, nc, nlevsoi)

    pp = _Ch4TP(; satpow=FT(params.satpow), sfgas=FT(params.scale_factor_gasdiff),
                  sfliq=FT(params.scale_factor_liqdiff), capthick=FT(params.capthick),
                  aereoxid=FT(params.aereoxid), om_frac_sf=FT(params.om_frac_sf),
                  dtime=FT(dtime), organic_max=FT(organic_max), smallnumber=FT(1.0e-12))

    # Struct-first kernel: manual backend + synchronize (the bundle args carry no backend).
    backend = _kernel_backend(scr.diffus)
    _ch4tran_column_kernel!(backend)(ts, tf, scr, c_atm_d, jwt_d, gc_d,
             ch4.grnd_ch4_cond_col, dcong, dconw, mask, pp,
             ch4vc.ch4frzout, ch4vc.use_aereoxid_prog, lake, sat, nlevsoi; ndrange = nc)
    KA.synchronize(backend)
    nothing
end

# ---------------------------------------------------------------------------
# ch4_init_gridcell_balance_check! — Gridcell balance check init
# ---------------------------------------------------------------------------

"""
    ch4_init_gridcell_balance_check!(ch4, mask_nolake, mask_lake,
        col_gridcell, col_wtgcell, dz, nlevsoi, ng, allowlakeprod)

Calculate beginning gridcell-level CH4 balance for mass conservation check.
Computes total column CH4, then aggregates to gridcell level (c2g).
Ported from `ch4_init_gridcell_balance_check` in `ch4Mod.F90`.
"""
function ch4_init_gridcell_balance_check!(ch4::CH4Data,
                                           mask_nolake::BitVector,
                                           mask_lake::BitVector,
                                           col_gridcell::Vector{Int},
                                           col_wtgcell::Vector{<:Real},
                                           dz::Matrix{<:Real},
                                           nlevsoi::Int,
                                           ng::Int,
                                           allowlakeprod::Bool)
    nc = length(mask_nolake)
    FT = eltype(dz)
    totcolch4_bef_col = zeros(FT, nc)

    ch4_totcolch4!(totcolch4_bef_col, ch4, mask_nolake, mask_lake,
                   dz, nlevsoi, allowlakeprod)

    # c2g: column-to-gridcell aggregation (unity scaling)
    ch4.totcolch4_bef_grc .= 0.0
    for c in 1:nc
        (mask_nolake[c] || mask_lake[c]) || continue
        g = col_gridcell[c]
        ch4.totcolch4_bef_grc[g] += totcolch4_bef_col[c] * col_wtgcell[c]
    end
    nothing
end

# ---------------------------------------------------------------------------
# ch4! — Main methane driver
# ---------------------------------------------------------------------------

"""
    ch4!(ch4, params, ch4vc,
         mask_soil, mask_soilp, mask_lake, mask_nolake,
         col_gridcell, col_wtgcell, patch_column, patch_itype, patch_wtcol,
         is_fates, latdeg,
         forc_pbot, forc_t, forc_po2, forc_pco2, forc_pch4,
         watsat, h2osoi_vol, h2osoi_liq, h2osoi_ice, h2osfc,
         bsw, cellorg, smp_l, t_soisno, t_grnd, t_h2osfc,
         frac_h2osfc, snow_depth, snl, qflx_surf,
         rootfr, rootfr_col, crootfr, rootr, elai,
         qflx_tran_veg, annsum_npp, rr,
         somhr, lithr, hr_vr, o_scalar, fphr, pot_f_nit_vr,
         lake_icefrac, lakedepth,
         z, dz, zi,
         nlevsoi, nlevsno, nlevdecomp, nlevdecomp_full,
         nlev_soildecomp_standard, mino2lim, organic_max,
         ng, noveg, dtime,
         use_cn, use_nitrif_denitrif, anoxia,
         agnpp, bgnpp, secsperyear)

Main timestep driver for the methane emissions model.
Orchestrates production, oxidation, aerenchyma, ebullition, and transport
for saturated and unsaturated fractions, and optionally lakes.
Ported from `ch4` in `ch4Mod.F90`.
"""
function ch4!(ch4d::CH4Data,
              params::CH4Params,
              ch4vc::CH4VarCon,
              mask_soil::BitVector,
              mask_soilp::BitVector,
              mask_lake::BitVector,
              mask_nolake::BitVector,
              col_gridcell::Vector{Int},
              col_wtgcell::Vector{<:Real},
              patch_column::Vector{Int},
              patch_itype::Vector{Int},
              patch_wtcol::Vector{<:Real},
              is_fates::BitVector,
              latdeg::Vector{<:Real},
              forc_pbot::Vector{<:Real},
              forc_t::Vector{<:Real},
              forc_po2::Vector{<:Real},
              forc_pco2::Vector{<:Real},
              forc_pch4::Vector{<:Real},
              watsat::Matrix{<:Real},
              h2osoi_vol::Matrix{<:Real},
              h2osoi_liq::Matrix{<:Real},
              h2osoi_ice::Matrix{<:Real},
              h2osfc::Vector{<:Real},
              bsw::Matrix{<:Real},
              cellorg::Matrix{<:Real},
              smp_l::Matrix{<:Real},
              t_soisno::Matrix{<:Real},
              t_grnd::Vector{<:Real},
              t_h2osfc::Vector{<:Real},
              frac_h2osfc::Vector{<:Real},
              snow_depth::Vector{<:Real},
              snl::Vector{Int},
              qflx_surf::Vector{<:Real},
              rootfr::Matrix{<:Real},
              rootfr_col::Matrix{<:Real},
              crootfr::Matrix{<:Real},
              rootr::Matrix{<:Real},
              elai::Vector{<:Real},
              qflx_tran_veg::Vector{<:Real},
              annsum_npp::Vector{<:Real},
              rr::Vector{<:Real},
              somhr::Vector{<:Real},
              lithr::Vector{<:Real},
              hr_vr::Matrix{<:Real},
              o_scalar::Matrix{<:Real},
              fphr::Matrix{<:Real},
              pot_f_nit_vr::Matrix{<:Real},
              lake_icefrac::Matrix{<:Real},
              lakedepth::Vector{<:Real},
              z::Matrix{<:Real},
              dz::Matrix{<:Real},
              zi::AbstractMatrix{<:Real},
              nlevsoi::Int,
              nlevsno::Int,
              nlevdecomp::Int,
              nlevdecomp_full::Int,
              nlev_soildecomp_standard::Int,
              mino2lim::Real,
              organic_max::Real,
              ng::Int,
              noveg::Int,
              dtime::Real,
              use_cn::Bool,
              use_nitrif_denitrif::Bool,
              anoxia::Bool,
              agnpp::Vector{<:Real},
              bgnpp::Vector{<:Real},
              secsperyear::Real)

    nc = length(mask_soil)
    np = length(mask_soilp)

    rgasm = RGAS / 1000.0
    dtime_ch4 = dtime

    redoxlag = params.redoxlag
    redoxlag_vertical = params.redoxlag_vertical
    atmch4 = params.atmch4
    qflxlagd = params.qflxlagd
    highlatfact = params.highlatfact

    redoxlags = redoxlag * SECSPDAY
    redoxlags_vertical = redoxlag_vertical * SECSPDAY

    jwt = fill(typemax(Int) >> 1, nc)  # large sentinel

    ch4_surf_flux_tot = ch4d.ch4_surf_flux_tot_col
    FT = eltype(t_soisno)
    ch4_prod_tot = zeros(FT, nc)
    ch4_oxid_tot = zeros(FT, nc)
    nem_col = zeros(FT, nc)
    fsat_bef = zeros(FT, nc)

    ch4_surf_flux_tot .= 0.0

    # Map gridcell-level forcing to column level
    forc_pbot_col = zeros(FT, nc)
    for c in 1:nc
        forc_pbot_col[c] = forc_pbot[col_gridcell[c]]
    end

    # --- Atmospheric concentrations ---
    for g in 1:ng
        if ch4vc.ch4offline
            forc_pch4_g = atmch4 * forc_pbot[g]
        else
            forc_pch4_g = forc_pch4[g]
        end
        ch4d.c_atm_grc[g, 1] = forc_pch4_g / rgasm / forc_t[g]
        ch4d.c_atm_grc[g, 2] = forc_po2[g] / rgasm / forc_t[g]
        ch4d.c_atm_grc[g, 3] = forc_pco2[g] / rgasm / forc_t[g]
    end

    # --- Save finundated before, calculate lagged surface runoff ---
    for c in eachindex(mask_soil)
        mask_soil[c] || continue
        g = col_gridcell[c]

        fsat_bef[c] = ch4d.finundated_col[c]

        if latdeg[g] < 45.0
            qflxlags = qflxlagd * SECSPDAY
        else
            qflxlags = qflxlagd * SECSPDAY * highlatfact
        end
        ch4d.qflx_surf_lag_col[c] = ch4d.qflx_surf_lag_col[c] * exp(-dtime / qflxlags) +
                                     qflx_surf[c] * (1.0 - exp(-dtime / qflxlags))
    end

    # --- Calculate finundated before snow and lagged version ---
    for c in eachindex(mask_soil)
        mask_soil[c] || continue

        if snow_depth[c] <= 0.0
            ch4d.finundated_col[c] = smooth_max(smooth_min(ch4d.finundated_col[c], 1.0), 0.0)
            ch4d.finundated_pre_snow_col[c] = ch4d.finundated_col[c]
        else
            ch4d.finundated_col[c] = ch4d.finundated_pre_snow_col[c]
        end

        if redoxlags > 0.0
            ch4d.finundated_lag_col[c] = ch4d.finundated_lag_col[c] * exp(-dtime / redoxlags) +
                                          ch4d.finundated_col[c] * (1.0 - exp(-dtime / redoxlags))
        else
            ch4d.finundated_lag_col[c] = ch4d.finundated_col[c]
        end
    end

    # --- Check finundated changes → adjust conc_ch4_sat or add flux ---
    for j in 1:nlevsoi
        for c in eachindex(mask_soil)
            mask_soil[c] || continue
            g = col_gridcell[c]

            if j == 1
                ch4d.ch4_dfsat_flux_col[c] = 0.0
            end

            if !ch4d.ch4_first_time_grc[g]
                if ch4d.finundated_col[c] > fsat_bef[c]
                    dfsat = ch4d.finundated_col[c] - fsat_bef[c]
                    ch4d.conc_ch4_sat_col[c, j] = (fsat_bef[c] * ch4d.conc_ch4_sat_col[c, j] +
                                                    dfsat * ch4d.conc_ch4_unsat_col[c, j]) /
                                                   ch4d.finundated_col[c]
                elseif ch4d.finundated_col[c] < fsat_bef[c]
                    ch4d.ch4_dfsat_flux_col[c] += (fsat_bef[c] - ch4d.finundated_col[c]) *
                        (ch4d.conc_ch4_sat_col[c, j] - ch4d.conc_ch4_unsat_col[c, j]) *
                        dz[c, j] / dtime * CATOMW / 1000.0
                end
            end
        end
    end

    # --- Annual averages ---
    ch4_annualupdate!(ch4d, mask_soil, mask_soilp, patch_column,
                      is_fates, somhr, agnpp, bgnpp, dtime, secsperyear)

    # --- p2c for grnd_ch4_cond ---
    for c in eachindex(mask_soil)
        mask_soil[c] || continue
        ch4d.grnd_ch4_cond_col[c] = 0.0
    end
    for p in eachindex(mask_soilp)
        mask_soilp[p] || continue
        c = patch_column[p]
        ch4d.grnd_ch4_cond_col[c] += ch4d.grnd_ch4_cond_patch[p] * patch_wtcol[p]
    end

    # ===================================================================
    # Loop over saturated (sat=1) and unsaturated (sat=0) for soil
    # ===================================================================
    lake = false

    for sat in 0:1
        if sat == 0
            get_jwt!(jwt, mask_soil, watsat, h2osoi_vol, t_soisno, nlevsoi, params)
            for c in eachindex(mask_soil)
                mask_soil[c] || continue
                ch4d.zwt_ch4_unsat_col[c] = jwt[c] > 0 ? zi[c, jwt[c]] : 0.0
            end

            # Update lagged saturation status
            for j in 1:nlevsoi
                for c in eachindex(mask_soil)
                    mask_soil[c] || continue
                    if j > jwt[c] && redoxlags_vertical > 0.0
                        ch4d.layer_sat_lag_col[c, j] = ch4d.layer_sat_lag_col[c, j] * exp(-dtime / redoxlags_vertical) +
                                                        (1.0 - exp(-dtime / redoxlags_vertical))
                    elseif redoxlags_vertical > 0.0
                        ch4d.layer_sat_lag_col[c, j] = ch4d.layer_sat_lag_col[c, j] * exp(-dtime / redoxlags_vertical)
                    elseif j > jwt[c]
                        ch4d.layer_sat_lag_col[c, j] = 1.0
                    else
                        ch4d.layer_sat_lag_col[c, j] = 0.0
                    end
                end
            end
        else  # saturated
            for c in eachindex(mask_soil)
                mask_soil[c] || continue
                jwt[c] = 0
            end
        end

        # CH4 production
        ch4_prod!(ch4d, params, ch4vc, mask_soil, mask_soilp,
                  patch_column, patch_itype, patch_wtcol, is_fates,
                  crootfr, rootfr_col, watsat, h2osoi_vol, t_soisno,
                  somhr, lithr, hr_vr, o_scalar, fphr, pot_f_nit_vr,
                  rr, jwt, sat, lake,
                  dz, z, zi, nlevsoi, nlevdecomp, nlevdecomp_full,
                  nlev_soildecomp_standard, mino2lim,
                  dtime, noveg, use_cn, use_nitrif_denitrif, anoxia)

        # CH4 oxidation
        ch4_oxid!(ch4d, params, mask_soil, watsat, h2osoi_vol,
                  smp_l, t_soisno, jwt, sat, lake, nlevsoi, dtime)

        # Aerenchyma transport
        ch4_aere!(ch4d, params, ch4vc, mask_soil, mask_soilp,
                  patch_column, patch_itype, patch_wtcol,
                  col_gridcell, is_fates,
                  watsat, h2osoi_vol, t_soisno, rootr, rootfr,
                  elai, qflx_tran_veg, annsum_npp,
                  z, dz, jwt, sat, lake, nlevsoi, dtime, noveg)

        # Ebullition
        ch4_ebul!(ch4d, params, mask_soil, watsat, h2osoi_vol,
                  t_soisno, forc_pbot_col, h2osfc, frac_h2osfc, lake_icefrac,
                  lakedepth, z, dz, zi, jwt, sat, lake, nlevsoi, dtime)

        # Transport / diffusion solve
        ch4_tran!(ch4d, params, ch4vc, mask_soil, col_gridcell,
                  watsat, h2osoi_vol, h2osoi_liq, h2osoi_ice, h2osfc,
                  bsw, cellorg, t_soisno, t_grnd, t_h2osfc,
                  frac_h2osfc, snow_depth, snl,
                  z, dz, zi, jwt, sat, lake,
                  nlevsoi, nlevsno, dtime_ch4, organic_max)
    end

    # ===================================================================
    # Optionally do lakes
    # ===================================================================
    if ch4vc.allowlakeprod
        lake = true
        sat = 1
        for c in eachindex(mask_lake)
            mask_lake[c] || continue
            jwt[c] = 0
        end

        # Create empty patch mask for lakes
        mask_lakep = falses(np)

        ch4_prod!(ch4d, params, ch4vc, mask_lake, mask_lakep,
                  patch_column, patch_itype, patch_wtcol, is_fates,
                  crootfr, rootfr_col, watsat, h2osoi_vol, t_soisno,
                  somhr, lithr, hr_vr, o_scalar, fphr, pot_f_nit_vr,
                  rr, jwt, sat, lake,
                  dz, z, zi, nlevsoi, nlevdecomp, nlevdecomp_full,
                  nlev_soildecomp_standard, mino2lim,
                  dtime, noveg, use_cn, use_nitrif_denitrif, anoxia)

        ch4_oxid!(ch4d, params, mask_lake, watsat, h2osoi_vol,
                  smp_l, t_soisno, jwt, sat, lake, nlevsoi, dtime)

        ch4_aere!(ch4d, params, ch4vc, mask_lake, mask_lakep,
                  patch_column, patch_itype, patch_wtcol,
                  col_gridcell, is_fates,
                  watsat, h2osoi_vol, t_soisno, rootr, rootfr,
                  elai, qflx_tran_veg, annsum_npp,
                  z, dz, jwt, sat, lake, nlevsoi, dtime, noveg)

        ch4_ebul!(ch4d, params, mask_lake, watsat, h2osoi_vol,
                  t_soisno, forc_pbot_col, h2osfc, frac_h2osfc, lake_icefrac,
                  lakedepth, z, dz, zi, jwt, sat, lake, nlevsoi, dtime)

        ch4_tran!(ch4d, params, ch4vc, mask_lake, col_gridcell,
                  watsat, h2osoi_vol, h2osoi_liq, h2osoi_ice, h2osfc,
                  bsw, cellorg, t_soisno, t_grnd, t_h2osfc,
                  frac_h2osfc, snow_depth, snl,
                  z, dz, zi, jwt, sat, lake,
                  nlevsoi, nlevsno, dtime_ch4, organic_max)
    end

    # ===================================================================
    # Average up to gridcell flux and column oxidation/production rate
    # ===================================================================

    # Weight soil columns by finundated
    for j in 1:nlevsoi
        for c in eachindex(mask_soil)
            mask_soil[c] || continue

            if j == 1
                totalsat = ch4d.ch4_surf_diff_sat_col[c] + ch4d.ch4_surf_aere_sat_col[c] +
                           ch4d.ch4_surf_ebul_sat_col[c]
                totalunsat = ch4d.ch4_surf_diff_unsat_col[c] + ch4d.ch4_surf_aere_unsat_col[c] +
                             ch4d.ch4_surf_ebul_unsat_col[c]
                ch4_surf_flux_tot[c] = (ch4d.finundated_col[c] * totalsat +
                                        (1.0 - ch4d.finundated_col[c]) * totalunsat) *
                                       CATOMW / 1000.0
            end

            ch4_oxid_tot[c] += (ch4d.finundated_col[c] * ch4d.ch4_oxid_depth_sat_col[c, j] +
                                (1.0 - ch4d.finundated_col[c]) * ch4d.ch4_oxid_depth_unsat_col[c, j]) *
                               dz[c, j] * CATOMW

            ch4_prod_tot[c] += (ch4d.finundated_col[c] * ch4d.ch4_prod_depth_sat_col[c, j] +
                                (1.0 - ch4d.finundated_col[c]) * ch4d.ch4_prod_depth_unsat_col[c, j]) *
                               dz[c, j] * CATOMW

            if j == nlevsoi
                nem_col[c] -= ch4_prod_tot[c]
                nem_col[c] += ch4_oxid_tot[c]
            end
        end
    end

    # Add dfsat flux correction
    for c in eachindex(mask_soil)
        mask_soil[c] || continue
        ch4_surf_flux_tot[c] += ch4d.ch4_dfsat_flux_col[c]
    end

    # Lake flux and diagnostics
    if ch4vc.allowlakeprod
        for j in 1:nlevsoi
            for c in eachindex(mask_lake)
                mask_lake[c] || continue

                if j == 1
                    totalsat = ch4d.ch4_surf_diff_sat_col[c] + ch4d.ch4_surf_aere_sat_col[c] +
                               ch4d.ch4_surf_ebul_sat_col[c]
                    ch4_surf_flux_tot[c] = totalsat * CATOMW / 1000.0
                end

                ch4_oxid_tot[c] += ch4d.ch4_oxid_depth_sat_col[c, j] * dz[c, j] * CATOMW
                ch4_prod_tot[c] += ch4d.ch4_prod_depth_sat_col[c, j] * dz[c, j] * CATOMW

                if !ch4vc.replenishlakec
                    ch4d.lake_soilc_col[c, j] -= 2.0 * ch4d.ch4_prod_depth_sat_col[c, j] * dtime * CATOMW
                end

                if j == nlevsoi
                    if !ch4vc.replenishlakec
                        nem_col[c] += ch4_prod_tot[c]
                    else
                        nem_col[c] -= ch4_prod_tot[c]
                    end
                    nem_col[c] += ch4_oxid_tot[c]
                end

                # Lake diagnostic output
                ch4d.ch4_prod_depth_lake_col[c, j] = ch4d.ch4_prod_depth_sat_col[c, j]
                ch4d.conc_ch4_lake_col[c, j] = ch4d.conc_ch4_sat_col[c, j]
                ch4d.conc_o2_lake_col[c, j] = ch4d.conc_o2_sat_col[c, j]
                ch4d.ch4_oxid_depth_lake_col[c, j] = ch4d.ch4_oxid_depth_sat_col[c, j]
                if j == 1
                    ch4d.ch4_surf_diff_lake_col[c] = ch4d.ch4_surf_diff_sat_col[c]
                    ch4d.ch4_surf_ebul_lake_col[c] = ch4d.ch4_surf_ebul_sat_col[c]
                end
            end
        end
    end

    # --- Finalize CH4 balance ---
    ch4_totcolch4!(ch4d.totcolch4_col, ch4d, mask_nolake, mask_lake,
                   dz, nlevsoi, ch4vc.allowlakeprod)

    # Column-level balance check
    for c in eachindex(mask_soil)
        mask_soil[c] || continue
        g = col_gridcell[c]

        if !ch4d.ch4_first_time_grc[g]
            errch4 = ch4d.totcolch4_col[c] - ch4d.totcolch4_bef_col[c] -
                     dtime * (ch4_prod_tot[c] - ch4_oxid_tot[c] -
                              ch4_surf_flux_tot[c] * 1000.0)
            if abs(errch4) > 1.0e-7
                @warn "Column-level CH4 conservation error" c errch4
            end
        end
    end

    if ch4vc.allowlakeprod
        for c in eachindex(mask_lake)
            mask_lake[c] || continue
            g = col_gridcell[c]

            if !ch4d.ch4_first_time_grc[g]
                errch4 = ch4d.totcolch4_col[c] - ch4d.totcolch4_bef_col[c] -
                         dtime * (ch4_prod_tot[c] - ch4_oxid_tot[c] -
                                  ch4_surf_flux_tot[c] * 1000.0)
                if abs(errch4) > 1.0e-7
                    @warn "Column-level CH4 conservation error (lake)" c errch4
                end
            end
        end
    end

    # --- c2g aggregation ---
    ch4d.ch4co2f_grc .= 0.0
    ch4d.ch4prodg_grc .= 0.0
    ch4d.totcolch4_grc .= 0.0
    nem_grc = zeros(FT, ng)
    ch4_surf_flux_tot_grc = zeros(FT, ng)

    for c in 1:nc
        (mask_soil[c] || mask_lake[c]) || continue
        g = col_gridcell[c]
        w = col_wtgcell[c]
        ch4d.ch4co2f_grc[g] += ch4_oxid_tot[c] * w
        ch4d.ch4prodg_grc[g] += ch4_prod_tot[c] * w
        ch4d.totcolch4_grc[g] += ch4d.totcolch4_col[c] * w
        nem_grc[g] += nem_col[c] * w
        ch4_surf_flux_tot_grc[g] += ch4_surf_flux_tot[c] * w
    end

    # Gridcell-level balance check
    for g in 1:ng
        if !ch4d.ch4_first_time_grc[g]
            errch4 = ch4d.totcolch4_grc[g] - ch4d.totcolch4_bef_grc[g] +
                     dtime * (nem_grc[g] + ch4_surf_flux_tot_grc[g] * 1000.0)
            if abs(errch4) > 1.0e-7
                @warn "Gridcell-level CH4 conservation error" g errch4
            end
        end
    end

    ch4d.ch4_first_time_grc .= false

    nothing
end
