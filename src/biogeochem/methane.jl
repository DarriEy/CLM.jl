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
Base.@kwdef mutable struct CH4Data{FT<:Real}
    # Per-layer production/consumption/transport rates (nc × nlevsoi)
    ch4_prod_depth_sat_col      ::Matrix{FT} = Matrix{FT}(undef, 0, 0)
    ch4_prod_depth_unsat_col    ::Matrix{FT} = Matrix{FT}(undef, 0, 0)
    ch4_prod_depth_lake_col     ::Matrix{FT} = Matrix{FT}(undef, 0, 0)
    ch4_oxid_depth_sat_col      ::Matrix{FT} = Matrix{FT}(undef, 0, 0)
    ch4_oxid_depth_unsat_col    ::Matrix{FT} = Matrix{FT}(undef, 0, 0)
    ch4_oxid_depth_lake_col     ::Matrix{FT} = Matrix{FT}(undef, 0, 0)
    ch4_aere_depth_sat_col      ::Matrix{FT} = Matrix{FT}(undef, 0, 0)
    ch4_aere_depth_unsat_col    ::Matrix{FT} = Matrix{FT}(undef, 0, 0)
    ch4_tran_depth_sat_col      ::Matrix{FT} = Matrix{FT}(undef, 0, 0)
    ch4_tran_depth_unsat_col    ::Matrix{FT} = Matrix{FT}(undef, 0, 0)
    ch4_ebul_depth_sat_col      ::Matrix{FT} = Matrix{FT}(undef, 0, 0)
    ch4_ebul_depth_unsat_col    ::Matrix{FT} = Matrix{FT}(undef, 0, 0)
    o2_oxid_depth_sat_col       ::Matrix{FT} = Matrix{FT}(undef, 0, 0)
    o2_oxid_depth_unsat_col     ::Matrix{FT} = Matrix{FT}(undef, 0, 0)
    o2_aere_depth_sat_col       ::Matrix{FT} = Matrix{FT}(undef, 0, 0)
    o2_aere_depth_unsat_col     ::Matrix{FT} = Matrix{FT}(undef, 0, 0)
    co2_decomp_depth_sat_col    ::Matrix{FT} = Matrix{FT}(undef, 0, 0)
    co2_decomp_depth_unsat_col  ::Matrix{FT} = Matrix{FT}(undef, 0, 0)
    co2_oxid_depth_sat_col      ::Matrix{FT} = Matrix{FT}(undef, 0, 0)
    co2_oxid_depth_unsat_col    ::Matrix{FT} = Matrix{FT}(undef, 0, 0)
    co2_aere_depth_sat_col      ::Matrix{FT} = Matrix{FT}(undef, 0, 0)
    co2_aere_depth_unsat_col    ::Matrix{FT} = Matrix{FT}(undef, 0, 0)

    # Column-level surface fluxes (nc)
    ch4_ebul_total_sat_col      ::Vector{FT} = FT[]
    ch4_ebul_total_unsat_col    ::Vector{FT} = FT[]
    ch4_surf_aere_sat_col       ::Vector{FT} = FT[]
    ch4_surf_aere_unsat_col     ::Vector{FT} = FT[]
    ch4_surf_ebul_sat_col       ::Vector{FT} = FT[]
    ch4_surf_ebul_unsat_col     ::Vector{FT} = FT[]
    ch4_surf_ebul_lake_col      ::Vector{FT} = FT[]
    ch4_surf_diff_sat_col       ::Vector{FT} = FT[]
    ch4_surf_diff_unsat_col     ::Vector{FT} = FT[]
    ch4_surf_diff_lake_col      ::Vector{FT} = FT[]
    ch4_dfsat_flux_col          ::Vector{FT} = FT[]
    ch4_surf_flux_tot_col       ::Vector{FT} = FT[]

    # Concentrations (nc × nlevsoi)
    conc_ch4_sat_col            ::Matrix{FT} = Matrix{FT}(undef, 0, 0)
    conc_ch4_unsat_col          ::Matrix{FT} = Matrix{FT}(undef, 0, 0)
    conc_ch4_lake_col           ::Matrix{FT} = Matrix{FT}(undef, 0, 0)
    conc_o2_sat_col             ::Matrix{FT} = Matrix{FT}(undef, 0, 0)
    conc_o2_unsat_col           ::Matrix{FT} = Matrix{FT}(undef, 0, 0)
    conc_o2_lake_col            ::Matrix{FT} = Matrix{FT}(undef, 0, 0)
    o2_decomp_depth_sat_col     ::Matrix{FT} = Matrix{FT}(undef, 0, 0)
    o2_decomp_depth_unsat_col   ::Matrix{FT} = Matrix{FT}(undef, 0, 0)

    # Stress factors (nc × nlevsoi)
    o2stress_sat_col            ::Matrix{FT} = Matrix{FT}(undef, 0, 0)
    o2stress_unsat_col          ::Matrix{FT} = Matrix{FT}(undef, 0, 0)
    ch4stress_sat_col           ::Matrix{FT} = Matrix{FT}(undef, 0, 0)
    ch4stress_unsat_col         ::Matrix{FT} = Matrix{FT}(undef, 0, 0)

    # Column-level state (nc)
    zwt_ch4_unsat_col           ::Vector{FT} = FT[]
    lake_soilc_col              ::Matrix{FT} = Matrix{FT}(undef, 0, 0)
    totcolch4_col               ::Vector{FT} = FT[]
    totcolch4_bef_col           ::Vector{FT} = FT[]
    annsum_counter_col          ::Vector{FT} = FT[]
    tempavg_somhr_col           ::Vector{FT} = FT[]
    annavg_somhr_col            ::Vector{FT} = FT[]
    tempavg_finrw_col           ::Vector{FT} = FT[]
    annavg_finrw_col            ::Vector{FT} = FT[]
    sif_col                     ::Vector{FT} = FT[]
    qflx_surf_lag_col           ::Vector{FT} = FT[]
    finundated_col              ::Vector{FT} = FT[]
    finundated_pre_snow_col     ::Vector{FT} = FT[]
    finundated_lag_col          ::Vector{FT} = FT[]
    layer_sat_lag_col           ::Matrix{FT} = Matrix{FT}(undef, 0, 0)
    pH_col                      ::Vector{FT} = FT[]

    # Gridcell-level (ng)
    c_atm_grc                   ::Matrix{FT} = Matrix{FT}(undef, 0, 0)
    ch4co2f_grc                 ::Vector{FT} = FT[]
    ch4prodg_grc                ::Vector{FT} = FT[]
    totcolch4_grc               ::Vector{FT} = FT[]
    totcolch4_bef_grc           ::Vector{FT} = FT[]

    # Patch-level aerenchyma (np)
    annavg_agnpp_patch          ::Vector{FT} = FT[]
    annavg_bgnpp_patch          ::Vector{FT} = FT[]
    tempavg_agnpp_patch         ::Vector{FT} = FT[]
    tempavg_bgnpp_patch         ::Vector{FT} = FT[]

    # Boundary layer conductance
    grnd_ch4_cond_patch         ::Vector{FT} = FT[]
    grnd_ch4_cond_col           ::Vector{FT} = FT[]

    # First-time flag (ng)
    ch4_first_time_grc          ::Vector{Bool}    = Bool[]
end

# ---------------------------------------------------------------------------
# get_jwt! — Water table layer identification
# ---------------------------------------------------------------------------

"""
    get_jwt!(jwt, mask_soil, watsat, h2osoi_vol, t_soisno, nlevsoi, params)

Find the first unsaturated layer going up (layer right above water table).
Also allows perched water table over ice.
Ported from `get_jwt` in `ch4Mod.F90`.
"""
function get_jwt!(jwt::Vector{Int},
                  mask_soil::BitVector,
                  watsat::Matrix{<:Real},
                  h2osoi_vol::Matrix{<:Real},
                  t_soisno::Matrix{<:Real},
                  nlevsoi::Int,
                  params::CH4Params)

    f_sat = params.f_sat

    for c in eachindex(mask_soil)
        mask_soil[c] || continue

        # Check for frozen saturated layers → perched water table
        perch = nlevsoi
        for j in nlevsoi:-1:1
            if t_soisno[c, j] < TFRZ && h2osoi_vol[c, j] > f_sat * watsat[c, j]
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
    nothing
end

# ---------------------------------------------------------------------------
# ch4_annualupdate! — Annual average update
# ---------------------------------------------------------------------------

"""
    ch4_annualupdate!(ch4, mask_soil, mask_soilp, patch_column, is_fates,
                      somhr, finundated, agnpp, bgnpp,
                      dt, secsperyear, nlevsoi)

Update annual mean fields for methane-related variables.
Ported from `ch4_annualupdate` in `ch4Mod.F90`.
"""
function ch4_annualupdate!(ch4::CH4Data,
                           mask_soil::BitVector,
                           mask_soilp::BitVector,
                           patch_column::Vector{Int},
                           is_fates::BitVector,
                           somhr::Vector{<:Real},
                           agnpp::Vector{<:Real},
                           bgnpp::Vector{<:Real},
                           dt::Real,
                           secsperyear::Real)

    for c in eachindex(mask_soil)
        mask_soil[c] || continue
        ch4.annsum_counter_col[c] += dt
    end

    for c in eachindex(mask_soil)
        mask_soil[c] || continue
        if ch4.annsum_counter_col[c] >= secsperyear
            ch4.annavg_somhr_col[c] = ch4.tempavg_somhr_col[c]
            ch4.tempavg_somhr_col[c] = 0.0
            if ch4.annavg_somhr_col[c] > 0.0
                ch4.annavg_finrw_col[c] = ch4.tempavg_finrw_col[c] / ch4.annavg_somhr_col[c]
            else
                ch4.annavg_finrw_col[c] = 0.0
            end
            ch4.tempavg_finrw_col[c] = 0.0
        else
            ch4.tempavg_somhr_col[c] += dt / secsperyear * somhr[c]
            ch4.tempavg_finrw_col[c] += dt / secsperyear * ch4.finundated_col[c] * somhr[c]
        end
    end

    for p in eachindex(mask_soilp)
        mask_soilp[p] || continue
        c = patch_column[p]
        if !is_fates[c]
            if ch4.annsum_counter_col[c] >= secsperyear
                ch4.annavg_agnpp_patch[p] = ch4.tempavg_agnpp_patch[p]
                ch4.tempavg_agnpp_patch[p] = 0.0
                ch4.annavg_bgnpp_patch[p] = ch4.tempavg_bgnpp_patch[p]
                ch4.tempavg_bgnpp_patch[p] = 0.0
            else
                ch4.tempavg_agnpp_patch[p] += dt / secsperyear * agnpp[p]
                ch4.tempavg_bgnpp_patch[p] += dt / secsperyear * bgnpp[p]
            end
        end
    end

    for c in eachindex(mask_soil)
        mask_soil[c] || continue
        if ch4.annsum_counter_col[c] >= secsperyear
            ch4.annsum_counter_col[c] = 0.0
        end
    end
    nothing
end

# ---------------------------------------------------------------------------
# ch4_totcolch4! — Total column CH4 calculation
# ---------------------------------------------------------------------------

"""
    ch4_totcolch4!(totcolch4, ch4, mask_nolake, mask_lake,
                   dz, col_landunit, lun_itype, nlevsoi, allowlakeprod)

Compute total column CH4 by integrating concentrations across soil layers.
Ported from `ch4_totcolch4` in `ch4Mod.F90`.
"""
function ch4_totcolch4!(totcolch4::Vector{<:Real},
                        ch4::CH4Data,
                        mask_nolake::BitVector,
                        mask_lake::BitVector,
                        dz::Matrix{<:Real},
                        nlevsoi::Int,
                        allowlakeprod::Bool)

    for c in eachindex(mask_nolake)
        mask_nolake[c] || continue
        totcolch4[c] = 0.0
    end
    for c in eachindex(mask_lake)
        mask_lake[c] || continue
        totcolch4[c] = 0.0
    end

    for j in 1:nlevsoi
        for c in eachindex(mask_nolake)
            mask_nolake[c] || continue
            totcolch4[c] += (ch4.finundated_col[c] * ch4.conc_ch4_sat_col[c, j] +
                             (1.0 - ch4.finundated_col[c]) * ch4.conc_ch4_unsat_col[c, j]) *
                            dz[c, j] * CATOMW
        end

        if allowlakeprod
            for c in eachindex(mask_lake)
                mask_lake[c] || continue
                totcolch4[c] += ch4.conc_ch4_sat_col[c, j] * dz[c, j] * CATOMW
            end
        end
    end
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
                   mask_soil::BitVector,
                   mask_soilp::BitVector,
                   patch_column::Vector{Int},
                   patch_itype::Vector{Int},
                   patch_wtcol::Vector{<:Real},
                   is_fates::BitVector,
                   crootfr::Matrix{<:Real},
                   rootfr_col::Matrix{<:Real},
                   watsat::Matrix{<:Real},
                   h2osoi_vol::Matrix{<:Real},
                   t_soisno::Matrix{<:Real},
                   somhr::Vector{<:Real},
                   lithr::Vector{<:Real},
                   hr_vr::Matrix{<:Real},
                   o_scalar::Matrix{<:Real},
                   fphr::Matrix{<:Real},
                   pot_f_nit_vr::Matrix{<:Real},
                   rr::Vector{<:Real},
                   jwt::Vector{Int},
                   sat::Int,
                   lake::Bool,
                   dz::Matrix{<:Real},
                   z::Matrix{<:Real},
                   zi::Matrix{<:Real},
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

    q10ch4 = params.q10ch4
    q10ch4base = params.q10ch4base
    f_ch4 = params.f_ch4
    rootlitfrac = params.rootlitfrac
    cnscalefactor = params.cnscalefactor
    lake_decomp_fact = params.lake_decomp_fact
    pHmax = params.pHmax
    pHmin = params.pHmin
    oxinhib = params.oxinhib
    q10lakebase = params.q10lakebase
    q10lake = q10ch4 * 1.5

    # Calculate vertically resolved column-averaged root respiration
    nc = length(mask_soil)
    FT = eltype(t_soisno)
    rr_vr = zeros(FT, nc, nlevsoi)
    if !lake
        for c in eachindex(mask_soil)
            mask_soil[c] || continue
            rr_vr[c, :] .= 0.0
        end
        for p in eachindex(mask_soilp)
            mask_soilp[p] || continue
            c = patch_column[p]
            if !is_fates[c]
                if patch_wtcol[p] > 0.0 && patch_itype[p] != noveg
                    for j in 1:nlevsoi
                        rr_vr[c, j] += rr[p] * crootfr[p, j] * patch_wtcol[p]
                    end
                end
            end
        end
    end

    for j in 1:nlevsoi
        for c in eachindex(mask_soil)
            mask_soil[c] || continue

            base_decomp = 0.0
            partition_z = 1.0

            if !lake
                if use_cn
                    base_decomp = (somhr[c] + lithr[c]) / CATOMW
                    if sat == 1
                        ch4.sif_col[c] = 1.0
                        if !anoxia
                            if ch4.annavg_finrw_col[c] != SPVAL
                                seasonalfin = smooth_max(ch4.finundated_col[c] - ch4.annavg_finrw_col[c], 0.0)
                                if seasonalfin > 0.0
                                    ch4.sif_col[c] = (ch4.annavg_finrw_col[c] + mino2lim * seasonalfin) / ch4.finundated_col[c]
                                    base_decomp *= ch4.sif_col[c]
                                end
                            end
                        end
                    end
                end
                base_decomp *= cnscalefactor
            else
                base_decomp = lake_decomp_fact * ch4.lake_soilc_col[c, j] * dz[c, j] *
                              q10lake^((t_soisno[c, j] - q10lakebase) / 10.0) / CATOMW
            end

            if t_soisno[c, j] <= TFRZ && (nlevdecomp == 1 || lake)
                base_decomp = 0.0
            end

            # Depth dependence
            if !lake
                if nlevdecomp == 1
                    if j <= nlev_soildecomp_standard
                        partition_z = rootfr_col[c, j] * rootlitfrac +
                                      (1.0 - rootlitfrac) * dz[c, j] / zi[c, nlev_soildecomp_standard]
                    else
                        partition_z = rootfr_col[c, j] * rootlitfrac
                    end
                else
                    if (somhr[c] + lithr[c]) > 0.0
                        partition_z = hr_vr[c, j] * dz[c, j] / (somhr[c] + lithr[c])
                    else
                        partition_z = 1.0
                    end
                end
            else
                partition_z = 1.0
            end

            # Adjust f_ch4 for temperature
            f_ch4_adj = 1.0
            if !lake
                t_fact_ch4 = q10ch4^((t_soisno[c, j] - q10ch4base) / 10.0)
                f_ch4_adj = f_ch4 * t_fact_ch4
                if ch4vc.ch4rmcnlim
                    jj = j > nlevdecomp ? 1 : j
                    if fphr[c, jj] > 0.0
                        f_ch4_adj /= fphr[c, jj]
                    end
                end
            else
                f_ch4_adj = 0.5
            end

            # pH factor
            if !lake && ch4vc.usephfact
                if ch4.pH_col[c] > pHmin && ch4.pH_col[c] < pHmax
                    pH_fact_ch4 = 10.0^(-0.2235 * ch4.pH_col[c]^2 + 2.7727 * ch4.pH_col[c] - 8.6)
                    f_ch4_adj *= pH_fact_ch4
                end
            end

            # Redox factor
            if !lake && sat == 1 && ch4.finundated_lag_col[c] < ch4.finundated_col[c]
                f_ch4_adj *= ch4.finundated_lag_col[c] / ch4.finundated_col[c]
            elseif sat == 0 && j > jwt[c]
                f_ch4_adj *= ch4.layer_sat_lag_col[c, j]
            end

            f_ch4_adj = smooth_min(f_ch4_adj, 0.5)

            # O2 decomposition demand
            o2_decomp_depth[c, j] = base_decomp * partition_z / dz[c, j]
            if anoxia
                if !lake && j > nlevdecomp
                    if o_scalar[c, 1] > 0.0
                        o2_decomp_depth[c, j] /= o_scalar[c, 1]
                    end
                elseif !lake
                    if o_scalar[c, j] > 0.0
                        o2_decomp_depth[c, j] /= o_scalar[c, j]
                    end
                end
            end

            # Add root respiration
            if !lake
                o2_decomp_depth[c, j] += rr_vr[c, j] / CATOMW / dz[c, j]
            end

            # Add nitrification O2 demand
            if use_nitrif_denitrif && !lake && j <= nlevdecomp_full
                o2_decomp_depth[c, j] += pot_f_nit_vr[c, j] * 2.0 / 14.0
            end

            # CH4 production
            if j > jwt[c]
                ch4_prod_depth[c, j] = f_ch4_adj * base_decomp * partition_z / dz[c, j]
            else
                if ch4vc.anoxicmicrosites
                    ch4_prod_depth[c, j] = f_ch4_adj * base_decomp * partition_z / dz[c, j] /
                                           (1.0 + oxinhib * conc_o2[c, j])
                else
                    ch4_prod_depth[c, j] = 0.0
                end
            end
        end
    end
    nothing
end

# ---------------------------------------------------------------------------
# ch4_oxid! — Methane oxidation
# ---------------------------------------------------------------------------

"""
    ch4_oxid!(ch4, params, mask_soil, watsat, h2osoi_vol, smp_l, t_soisno,
              jwt, sat, lake, nlevsoi, dtime)

Calculate CH4 oxidation in each soil layer using double Michaelis-Menten kinetics.
Ported from `ch4_oxid` in `ch4Mod.F90`.
"""
function ch4_oxid!(ch4::CH4Data,
                   params::CH4Params,
                   mask_soil::BitVector,
                   watsat::Matrix{<:Real},
                   h2osoi_vol::Matrix{<:Real},
                   smp_l::Matrix{<:Real},
                   t_soisno::Matrix{<:Real},
                   jwt::Vector{Int},
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

    vmax_ch4_oxid = params.vmax_ch4_oxid
    k_m = params.k_m
    q10_ch4oxid = params.q10_ch4oxid
    smp_crit = params.smp_crit
    k_m_o2 = params.k_m_o2
    k_m_unsat = params.k_m_unsat
    vmax_oxid_unsat = params.vmax_oxid_unsat

    t0 = TFRZ + 12.0

    for j in 1:nlevsoi
        for c in eachindex(mask_soil)
            mask_soil[c] || continue

            if sat == 1 || j > jwt[c]
                k_m_eff = k_m
                vmax_eff = vmax_ch4_oxid
            else
                k_m_eff = k_m_unsat
                vmax_eff = vmax_oxid_unsat
            end

            porevol = smooth_max(watsat[c, j] - h2osoi_vol[c, j], 0.0)
            h2osoi_vol_min = smooth_min(watsat[c, j], h2osoi_vol[c, j])

            if j <= jwt[c] && smp_l[c, j] < 0.0
                smp_fact = exp(-smp_l[c, j] / smp_crit)
            else
                smp_fact = 1.0
            end

            if j <= jwt[c]
                k_h_inv = exp(-C_H_INV[1] * (1.0 / t_soisno[c, j] - 1.0 / KH_TBASE) + log(KH_THETA[1]))
                k_h_cc = t_soisno[c, j] / k_h_inv * rgasLatm
                conc_ch4_rel = conc_ch4[c, j] / (h2osoi_vol_min + porevol / k_h_cc)

                k_h_inv = exp(-C_H_INV[2] * (1.0 / t_soisno[c, j] - 1.0 / KH_TBASE) + log(KH_THETA[2]))
                k_h_cc = t_soisno[c, j] / k_h_inv * rgasLatm
                conc_o2_rel = conc_o2[c, j] / (h2osoi_vol_min + porevol / k_h_cc)
            else
                conc_ch4_rel = conc_ch4[c, j] / watsat[c, j]
                conc_o2_rel = conc_o2[c, j] / watsat[c, j]
            end

            oxid_a = vmax_eff * h2osoi_vol_min * conc_ch4_rel / (k_m_eff + conc_ch4_rel) *
                     conc_o2_rel / (k_m_o2 + conc_o2_rel) *
                     q10_ch4oxid^((t_soisno[c, j] - t0) / 10.0) * smp_fact

            if t_soisno[c, j] <= TFRZ
                oxid_a = 0.0
            end

            ch4_oxid_depth[c, j] = oxid_a
            o2_oxid_depth[c, j] = ch4_oxid_depth[c, j] * 2.0
        end
    end
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
# ch4_aere! — Aerenchyma transport
# ---------------------------------------------------------------------------

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
                   mask_soil::BitVector,
                   mask_soilp::BitVector,
                   patch_column::Vector{Int},
                   patch_itype::Vector{Int},
                   patch_wtcol::Vector{<:Real},
                   col_gridcell::Vector{Int},
                   is_fates::BitVector,
                   watsat::Matrix{<:Real},
                   h2osoi_vol::Matrix{<:Real},
                   t_soisno::Matrix{<:Real},
                   rootr::Matrix{<:Real},
                   rootfr::Matrix{<:Real},
                   elai::Vector{<:Real},
                   qflx_tran_veg::Vector{<:Real},
                   annsum_npp::Vector{<:Real},
                   z::Matrix{<:Real},
                   dz::Matrix{<:Real},
                   jwt::Vector{Int},
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

    # Initialize
    for j in 1:nlevsoi
        for c in eachindex(mask_soil)
            mask_soil[c] || continue
            ch4_aere_depth[c, j] = 0.0
            ch4_tran_depth[c, j] = 0.0
            o2_aere_depth[c, j] = 0.0
        end
    end

    if !lake
        FT_ae = eltype(t_soisno)
        tranloss = zeros(FT_ae, nlevsoi)
        aere = zeros(FT_ae, nlevsoi)
        oxaere = zeros(FT_ae, nlevsoi)

        for p in eachindex(mask_soilp)
            mask_soilp[p] || continue
            c = patch_column[p]
            g = col_gridcell[c]

            if !is_fates[c]
                is_vegetated = patch_itype[p] != noveg
                itype = patch_itype[p]
                # Grass porosity check (simplified: grasses get 0.3, others scaled)
                poros_tiller = 0.3 * params.nongrassporosratio  # default non-grass
                # Note: in full CLM, grasses (nc3_arctic_grass, etc.) would get 0.3
                rootfr_vr = rootfr[p, 1:nlevsoi]
            else
                is_vegetated = true
                poros_tiller = 0.3 * params.nongrassporosratio
                rootfr_vr = rootfr[p, 1:nlevsoi]
            end

            site_ox_aere!(tranloss, aere, oxaere,
                          is_vegetated,
                          view(watsat, c, 1:nlevsoi),
                          view(h2osoi_vol, c, 1:nlevsoi),
                          view(t_soisno, c, 1:nlevsoi),
                          view(conc_ch4, c, 1:nlevsoi),
                          view(rootr, p, 1:nlevsoi),
                          qflx_tran_veg[p],
                          jwt[c],
                          annsum_npp[p],
                          ch4.annavg_agnpp_patch[p],
                          ch4.annavg_bgnpp_patch[p],
                          elai[p],
                          poros_tiller,
                          rootfr_vr,
                          ch4.grnd_ch4_cond_patch[p],
                          view(conc_o2, c, 1:nlevsoi),
                          view(ch4.c_atm_grc, g, 1:2),
                          view(z, c, 1:nlevsoi),
                          view(dz, c, 1:nlevsoi),
                          sat, nlevsoi, params, ch4vc)

            for j in 1:nlevsoi
                aeretran = smooth_min(aere[j] + tranloss[j], conc_ch4[c, j] / dtime + ch4_prod_depth[c, j])
                ch4_aere_depth[c, j] += aeretran * patch_wtcol[p]
                ch4_tran_depth[c, j] += smooth_min(tranloss[j], aeretran) * patch_wtcol[p]
                o2_aere_depth[c, j] += oxaere[j] * patch_wtcol[p]
            end
        end
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
                   zi::Matrix{<:Real},
                   jwt::Vector{Int},
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

    for j in 1:nlevsoi
        for c in eachindex(mask_soil)
            mask_soil[c] || continue

            if j > jwt[c] && t_soisno[c, j] > TFRZ
                k_h_inv = exp(-C_H_INV[1] * (1.0 / t_soisno[c, j] - 1.0 / KH_TBASE) + log(KH_THETA[1]))
                k_h = 1.0 / k_h_inv
                k_h_cc = t_soisno[c, j] * k_h * rgasLatm

                if !lake
                    zi_jwt = jwt[c] > 0 ? zi[c, jwt[c]] : 0.0
                    pressure = forc_pbot[c] + DENH2O * GRAV * (z[c, j] - zi_jwt)
                    if sat == 1 && frac_h2osfc[c] > 0.0
                        pressure += DENH2O * GRAV * h2osfc[c] / 1000.0 / frac_h2osfc[c]
                    end
                else
                    pressure = forc_pbot[c] + DENH2O * GRAV * (z[c, j] + lakedepth[c])
                end

                vgc = conc_ch4[c, j] / watsat[c, j] / k_h_cc * rgasm * t_soisno[c, j] / pressure

                if vgc > vgc_max * bubble_f
                    ch4_ebul_depth[c, j] = (vgc - vgc_min * bubble_f) * conc_ch4[c, j] / ebul_timescale
                else
                    ch4_ebul_depth[c, j] = 0.0
                end
            else
                ch4_ebul_depth[c, j] = 0.0
            end

            if lake && lake_icefrac[c, 1] > 0.1
                ch4_ebul_depth[c, j] = 0.0
            end
        end
    end
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
function ch4_tran!(ch4::CH4Data,
                   params::CH4Params,
                   ch4vc::CH4VarCon,
                   mask_soil::BitVector,
                   col_gridcell::Vector{Int},
                   watsat::Matrix{<:Real},
                   h2osoi_vol::Matrix{<:Real},
                   h2osoi_liq::Matrix{<:Real},
                   h2osoi_ice::Matrix{<:Real},
                   h2osfc::Vector{<:Real},
                   bsw::Matrix{<:Real},
                   cellorg::Matrix{<:Real},
                   t_soisno::Matrix{<:Real},
                   t_grnd::Vector{<:Real},
                   t_h2osfc::Vector{<:Real},
                   frac_h2osfc::Vector{<:Real},
                   snow_depth::Vector{<:Real},
                   snl::Vector{Int},
                   z::Matrix{<:Real},
                   dz::Matrix{<:Real},
                   zi::Matrix{<:Real},
                   jwt::Vector{Int},
                   sat::Int,
                   lake::Bool,
                   nlevsoi::Int,
                   nlevsno::Int,
                   dtime_ch4::Real,
                   organic_max::Real)

    smallnumber = 1.0e-12
    dtime = dtime_ch4

    # Select sat/unsat arrays
    if sat == 0
        o2_decomp_depth = ch4.o2_decomp_depth_unsat_col
        o2stress = ch4.o2stress_unsat_col
        ch4_oxid_depth = ch4.ch4_oxid_depth_unsat_col
        ch4_prod_depth = ch4.ch4_prod_depth_unsat_col
        ch4_aere_depth = ch4.ch4_aere_depth_unsat_col
        ch4_surf_aere = ch4.ch4_surf_aere_unsat_col
        ch4_ebul_depth = ch4.ch4_ebul_depth_unsat_col
        ch4_ebul_total = ch4.ch4_ebul_total_unsat_col
        ch4_surf_ebul = ch4.ch4_surf_ebul_unsat_col
        ch4_surf_diff = ch4.ch4_surf_diff_unsat_col
        o2_oxid_depth = ch4.o2_oxid_depth_unsat_col
        o2_aere_depth = ch4.o2_aere_depth_unsat_col
        ch4stress = ch4.ch4stress_unsat_col
        co2_decomp_depth = ch4.co2_decomp_depth_unsat_col
        conc_ch4 = ch4.conc_ch4_unsat_col
        conc_o2 = ch4.conc_o2_unsat_col
    else
        o2_decomp_depth = ch4.o2_decomp_depth_sat_col
        o2stress = ch4.o2stress_sat_col
        ch4_oxid_depth = ch4.ch4_oxid_depth_sat_col
        ch4_prod_depth = ch4.ch4_prod_depth_sat_col
        ch4_aere_depth = ch4.ch4_aere_depth_sat_col
        ch4_surf_aere = ch4.ch4_surf_aere_sat_col
        ch4_ebul_depth = ch4.ch4_ebul_depth_sat_col
        ch4_ebul_total = ch4.ch4_ebul_total_sat_col
        ch4_surf_ebul = ch4.ch4_surf_ebul_sat_col
        ch4_surf_diff = ch4.ch4_surf_diff_sat_col
        o2_oxid_depth = ch4.o2_oxid_depth_sat_col
        o2_aere_depth = ch4.o2_aere_depth_sat_col
        ch4stress = ch4.ch4stress_sat_col
        co2_decomp_depth = ch4.co2_decomp_depth_sat_col
        conc_ch4 = ch4.conc_ch4_sat_col
        conc_o2 = ch4.conc_o2_sat_col
    end

    satpow = params.satpow
    scale_factor_gasdiff = params.scale_factor_gasdiff
    scale_factor_liqdiff = params.scale_factor_liqdiff
    capthick = params.capthick
    aereoxid = params.aereoxid

    nc = length(mask_soil)
    c_atm = ch4.c_atm_grc

    # --- Competition for O2 and CH4 ---
    for j in 1:nlevsoi
        for c in eachindex(mask_soil)
            mask_soil[c] || continue

            o2demand = o2_decomp_depth[c, j] + o2_oxid_depth[c, j]
            if o2demand > 0.0
                if (conc_o2[c, j] / dtime + o2_aere_depth[c, j]) > o2demand
                    o2stress[c, j] = 1.0
                else
                    o2stress[c, j] = (conc_o2[c, j] / dtime + o2_aere_depth[c, j]) / o2demand
                end
            else
                o2stress[c, j] = 1.0
            end

            ch4demand = ch4_oxid_depth[c, j] + ch4_aere_depth[c, j] + ch4_ebul_depth[c, j]
            if ch4demand > 0.0
                ch4stress[c, j] = smooth_min((conc_ch4[c, j] / dtime + ch4_prod_depth[c, j]) / ch4demand, 1.0)
            else
                ch4stress[c, j] = 1.0
            end

            # Resolve competition
            if o2stress[c, j] < 1.0 || ch4stress[c, j] < 1.0
                if ch4stress[c, j] <= o2stress[c, j]
                    if o2stress[c, j] < 1.0
                        o2demand2 = o2_decomp_depth[c, j]
                        if o2demand2 > 0.0
                            o2stress[c, j] = smooth_min((conc_o2[c, j] / dtime + o2_aere_depth[c, j] -
                                                   ch4stress[c, j] * o2_oxid_depth[c, j]) / o2demand2, 1.0)
                        else
                            o2stress[c, j] = 1.0
                        end
                    end
                    ch4_oxid_depth[c, j] *= ch4stress[c, j]
                    o2_oxid_depth[c, j] *= ch4stress[c, j]
                else
                    if ch4stress[c, j] < 1.0
                        ch4demand2 = ch4_aere_depth[c, j] + ch4_ebul_depth[c, j]
                        if ch4demand2 > 0.0
                            ch4stress[c, j] = smooth_min((conc_ch4[c, j] / dtime + ch4_prod_depth[c, j] -
                                                    o2stress[c, j] * ch4_oxid_depth[c, j]) / ch4demand2, 1.0)
                        else
                            ch4stress[c, j] = 1.0
                        end
                    end
                    ch4_oxid_depth[c, j] *= o2stress[c, j]
                    o2_oxid_depth[c, j] *= o2stress[c, j]
                end
            end

            ch4_aere_depth[c, j] *= ch4stress[c, j]
            ch4_ebul_depth[c, j] *= ch4stress[c, j]
            o2_decomp_depth[c, j] *= o2stress[c, j]
        end
    end

    # --- Accumulate ebullition ---
    for j in 1:nlevsoi
        for c in eachindex(mask_soil)
            mask_soil[c] || continue
            if j == 1
                ch4_ebul_total[c] = 0.0
            end
            ch4_ebul_total[c] += ch4_ebul_depth[c, j] * dz[c, j]
        end
    end

    # --- Set up and solve diffusion for each species ---
    # Allocate working arrays
    FT = eltype(t_soisno)
    k_h_cc_arr = zeros(FT, nc, nlevsoi + 1, 2)  # 0:nlevsoi → 1:nlevsoi+1
    epsilon_t = zeros(FT, nc, nlevsoi, 2)
    source = zeros(FT, nc, nlevsoi, 2)
    conc_ch4_bef = zeros(FT, nc, nlevsoi)
    conc_ch4_rel = zeros(FT, nc, nlevsoi + 1)
    conc_o2_rel = zeros(FT, nc, nlevsoi + 1)
    h2osoi_vol_min = zeros(FT, nc, nlevsoi)
    liqfrac = ones(FT, nc, nlevsoi)

    # Henry's law coefficients
    for j in 0:nlevsoi
        for c in eachindex(mask_soil)
            mask_soil[c] || continue
            g = col_gridcell[c]
            for s in 1:2
                if j == 0
                    k_h_inv = exp(-C_H_INV[s] * (1.0 / t_grnd[c] - 1.0 / KH_TBASE) + log(KH_THETA[s]))
                    k_h_cc_arr[c, j+1, s] = t_grnd[c] / k_h_inv * rgasLatm
                else
                    k_h_inv = exp(-C_H_INV[s] * (1.0 / t_soisno[c, j] - 1.0 / KH_TBASE) + log(KH_THETA[s]))
                    k_h_cc_arr[c, j+1, s] = t_soisno[c, j] / k_h_inv * rgasLatm
                end
            end
        end
    end

    # Source terms and epsilon_t
    for j in 1:nlevsoi
        for c in eachindex(mask_soil)
            mask_soil[c] || continue
            g = col_gridcell[c]

            h2osoi_vol_min[c, j] = smooth_min(watsat[c, j], h2osoi_vol[c, j])
            if ch4vc.ch4frzout
                liqfrac[c, j] = smooth_max(0.05, (h2osoi_liq[c, j] / DENH2O + smallnumber) /
                                           (h2osoi_liq[c, j] / DENH2O + h2osoi_ice[c, j] / DENICE + smallnumber))
            else
                liqfrac[c, j] = 1.0
            end

            if j <= jwt[c]
                for s in 1:2
                    epsilon_t[c, j, s] = watsat[c, j] - (1.0 - k_h_cc_arr[c, j+1, s]) * h2osoi_vol_min[c, j] * liqfrac[c, j]
                end
            else
                for s in 1:2
                    epsilon_t[c, j, s] = watsat[c, j] * liqfrac[c, j]
                end
            end

            if !ch4vc.use_aereoxid_prog
                ch4_oxid_depth[c, j] += aereoxid * ch4_aere_depth[c, j]
                ch4_aere_depth[c, j] -= aereoxid * ch4_aere_depth[c, j]
            end

            source[c, j, 1] = ch4_prod_depth[c, j] - ch4_oxid_depth[c, j] -
                              ch4_aere_depth[c, j] - ch4_ebul_depth[c, j]
            source[c, j, 2] = -o2_oxid_depth[c, j] - o2_decomp_depth[c, j] + o2_aere_depth[c, j]

            conc_ch4_bef[c, j] = conc_ch4[c, j]
        end
    end

    # Accumulate aerenchyma surface flux
    for j in 1:nlevsoi
        for c in eachindex(mask_soil)
            mask_soil[c] || continue
            if j == 1
                ch4_surf_aere[c] = 0.0
            end
            ch4_surf_aere[c] += ch4_aere_depth[c, j] * dz[c, j]
        end
    end

    # Add ebullition to source at jwt layer
    for c in eachindex(mask_soil)
        mask_soil[c] || continue
        if jwt[c] != 0
            source[c, jwt[c], 1] += ch4_ebul_total[c] / dz[c, jwt[c]]
        end
    end

    # Relative concentrations
    for j in 0:nlevsoi
        for c in eachindex(mask_soil)
            mask_soil[c] || continue
            g = col_gridcell[c]
            if j == 0
                conc_ch4_rel[c, 1] = c_atm[g, 1]
                conc_o2_rel[c, 1] = c_atm[g, 2]
            else
                conc_ch4_rel[c, j+1] = conc_ch4[c, j] / epsilon_t[c, j, 1]
                conc_o2_rel[c, j+1] = conc_o2[c, j] / epsilon_t[c, j, 2]
            end
        end
    end

    # Tridiagonal solver arrays
    diffus = zeros(FT, nc, nlevsoi)
    dp1_zp1 = zeros(FT, nc, nlevsoi)
    dm1_zm1 = zeros(FT, nc, nlevsoi)
    at_arr = zeros(FT, nc, nlevsoi + 1)
    bt_arr = zeros(FT, nc, nlevsoi + 1)
    ct_arr = zeros(FT, nc, nlevsoi + 1)
    rt_arr = zeros(FT, nc, nlevsoi + 1)
    spec_grnd_cond = zeros(FT, nc, 2)
    conc_ch4_rel_old = copy(conc_ch4_rel)

    for s in 1:2
        # Snow/pond resistance and ground conductance
        for c in eachindex(mask_soil)
            mask_soil[c] || continue
            if ch4.grnd_ch4_cond_col[c] < smallnumber && s == 1
                ch4.grnd_ch4_cond_col[c] = smallnumber
            end
            snowres = 0.0
            # Simplified: skip snow layers for now, just use ground conductance
            pondres = 0.0
            if !lake && sat == 1 && frac_h2osfc[c] > 0.0
                if t_h2osfc[c] >= TFRZ
                    t_soisno_c = t_h2osfc[c] - TFRZ
                    ponddiff = (D_CON_W[s, 1] + D_CON_W[s, 2] * t_soisno_c + D_CON_W[s, 3] * t_soisno_c^2) * 1.0e-9 *
                               scale_factor_liqdiff
                    pondz = h2osfc[c] / 1000.0 / frac_h2osfc[c]
                    pondres = pondz / ponddiff
                elseif h2osfc[c] / frac_h2osfc[c] > capthick
                    pondres = 1.0 / smallnumber
                end
            end
            spec_grnd_cond[c, s] = 1.0 / (1.0 / ch4.grnd_ch4_cond_col[c] + snowres + pondres)
        end

        # Gas diffusivity
        for j in 1:nlevsoi
            for c in eachindex(mask_soil)
                mask_soil[c] || continue

                t_soisno_c = t_soisno[c, j] - TFRZ

                if j <= jwt[c]
                    f_a = 1.0 - h2osoi_vol_min[c, j] / watsat[c, j]
                    eps = watsat[c, j] - h2osoi_vol_min[c, j]
                    if organic_max > 0.0
                        om_frac = smooth_min(params.om_frac_sf * cellorg[c, j] / organic_max, 1.0)
                    else
                        om_frac = 1.0
                    end
                    diffus[c, j] = (D_CON_G[s, 1] + D_CON_G[s, 2] * t_soisno_c) * 1.0e-4 *
                                   (om_frac * f_a^(10.0 / 3.0) / watsat[c, j]^2 +
                                    (1.0 - om_frac) * eps^2 * f_a^(3.0 / bsw[c, j])) *
                                   scale_factor_gasdiff
                else
                    eps = watsat[c, j]
                    diffus[c, j] = eps^satpow * (D_CON_W[s, 1] + D_CON_W[s, 2] * t_soisno_c +
                                                  D_CON_W[s, 3] * t_soisno_c^2) * 1.0e-9 * scale_factor_liqdiff
                    if t_soisno[c, j] <= TFRZ
                        diffus[c, j] *= (h2osoi_liq[c, j] / DENH2O + smallnumber) /
                                        (h2osoi_liq[c, j] / DENH2O + h2osoi_ice[c, j] / DENICE + smallnumber)
                    end
                end
                diffus[c, j] = smooth_max(diffus[c, j], smallnumber)
            end
        end

        # Tridiagonal coefficients dm1_zm1 and dp1_zp1
        for j in 1:nlevsoi
            for c in eachindex(mask_soil)
                mask_soil[c] || continue
                if j == 1 && j != jwt[c] && j != jwt[c] + 1
                    dm1_zm1[c, j] = 1.0 / (1.0 / spec_grnd_cond[c, s] + dz[c, j] / (diffus[c, j] * 2.0))
                    dp1_zp1[c, j] = j < nlevsoi ? 2.0 / (dz[c, j] / diffus[c, j] + dz[c, j+1] / diffus[c, j+1]) : 0.0
                elseif j == 1 && j == jwt[c]
                    dm1_zm1[c, j] = 1.0 / (1.0 / spec_grnd_cond[c, s] + dz[c, j] / (diffus[c, j] * 2.0))
                    dp1_zp1[c, j] = j < nlevsoi ? 2.0 / (dz[c, j] * k_h_cc_arr[c, j+1, s] / diffus[c, j] + dz[c, j+1] / diffus[c, j+1]) : 0.0
                elseif j == 1
                    dm1_zm1[c, j] = 1.0 / (k_h_cc_arr[c, j, s] / spec_grnd_cond[c, s] + dz[c, j] / (diffus[c, j] * 2.0))
                    dp1_zp1[c, j] = j < nlevsoi ? 2.0 / (dz[c, j] / diffus[c, j] + dz[c, j+1] / diffus[c, j+1]) : 0.0
                elseif j < nlevsoi && j != jwt[c] && j != jwt[c] + 1
                    dm1_zm1[c, j] = 2.0 / (dz[c, j] / diffus[c, j] + dz[c, j-1] / diffus[c, j-1])
                    dp1_zp1[c, j] = 2.0 / (dz[c, j] / diffus[c, j] + dz[c, j+1] / diffus[c, j+1])
                elseif j < nlevsoi && j == jwt[c]
                    dm1_zm1[c, j] = 2.0 / (dz[c, j] / diffus[c, j] + dz[c, j-1] / diffus[c, j-1])
                    dp1_zp1[c, j] = 2.0 / (dz[c, j] * k_h_cc_arr[c, j+1, s] / diffus[c, j] + dz[c, j+1] / diffus[c, j+1])
                elseif j < nlevsoi  # j == jwt+1
                    dm1_zm1[c, j] = 2.0 / (dz[c, j] / diffus[c, j] + dz[c, j-1] * k_h_cc_arr[c, j, s] / diffus[c, j-1])
                    dp1_zp1[c, j] = 2.0 / (dz[c, j] / diffus[c, j] + dz[c, j+1] / diffus[c, j+1])
                elseif j != jwt[c] + 1  # j == nlevsoi
                    dm1_zm1[c, j] = 2.0 / (dz[c, j] / diffus[c, j] + dz[c, j-1] / diffus[c, j-1])
                else  # jwt == nlevsoi-1
                    dm1_zm1[c, j] = 2.0 / (dz[c, j] / diffus[c, j] + dz[c, j-1] * k_h_cc_arr[c, j, s] / diffus[c, j-1])
                end
            end
        end

        # Select concentration array for this species
        conc_rel = s == 1 ? conc_ch4_rel : conc_o2_rel

        # Build tridiagonal system (j=0 is index 1 in Julia)
        for j in 0:nlevsoi
            for c in eachindex(mask_soil)
                mask_soil[c] || continue
                g = col_gridcell[c]

                jj = j + 1  # Julia 1-based index

                if j == 0
                    at_arr[c, jj] = 0.0
                    bt_arr[c, jj] = 1.0
                    ct_arr[c, jj] = 0.0
                    rt_arr[c, jj] = c_atm[g, s]
                elseif j < nlevsoi && j == jwt[c]
                    dzj = dz[c, j]
                    at_arr[c, jj] = -0.5 / dzj * dm1_zm1[c, j]
                    bt_arr[c, jj] = epsilon_t[c, j, s] / dtime_ch4 + 0.5 / dzj * (dp1_zp1[c, j] * k_h_cc_arr[c, j+1, s] + dm1_zm1[c, j])
                    ct_arr[c, jj] = -0.5 / dzj * dp1_zp1[c, j]
                    rt_arr[c, jj] = epsilon_t[c, j, s] / dtime_ch4 * conc_rel[c, jj] +
                                    0.5 / dzj * (dp1_zp1[c, j] * (conc_rel[c, jj+1] - conc_rel[c, jj] * k_h_cc_arr[c, j+1, s]) -
                                                 dm1_zm1[c, j] * (conc_rel[c, jj] - conc_rel[c, jj-1])) +
                                    source[c, j, s]
                elseif j < nlevsoi && j == jwt[c] + 1
                    dzj = dz[c, j]
                    at_arr[c, jj] = -0.5 / dzj * dm1_zm1[c, j] * k_h_cc_arr[c, j, s]
                    bt_arr[c, jj] = epsilon_t[c, j, s] / dtime_ch4 + 0.5 / dzj * (dp1_zp1[c, j] + dm1_zm1[c, j])
                    ct_arr[c, jj] = -0.5 / dzj * dp1_zp1[c, j]
                    rt_arr[c, jj] = epsilon_t[c, j, s] / dtime_ch4 * conc_rel[c, jj] +
                                    0.5 / dzj * (dp1_zp1[c, j] * (conc_rel[c, jj+1] - conc_rel[c, jj]) -
                                                 dm1_zm1[c, j] * (conc_rel[c, jj] - conc_rel[c, jj-1] * k_h_cc_arr[c, j, s])) +
                                    source[c, j, s]
                elseif j < nlevsoi
                    dzj = dz[c, j]
                    at_arr[c, jj] = -0.5 / dzj * dm1_zm1[c, j]
                    bt_arr[c, jj] = epsilon_t[c, j, s] / dtime_ch4 + 0.5 / dzj * (dp1_zp1[c, j] + dm1_zm1[c, j])
                    ct_arr[c, jj] = -0.5 / dzj * dp1_zp1[c, j]
                    rt_arr[c, jj] = epsilon_t[c, j, s] / dtime_ch4 * conc_rel[c, jj] +
                                    0.5 / dzj * (dp1_zp1[c, j] * (conc_rel[c, jj+1] - conc_rel[c, jj]) -
                                                 dm1_zm1[c, j] * (conc_rel[c, jj] - conc_rel[c, jj-1])) +
                                    source[c, j, s]
                elseif j == nlevsoi && j == jwt[c] + 1
                    dzj = dz[c, j]
                    at_arr[c, jj] = -0.5 / dzj * dm1_zm1[c, j] * k_h_cc_arr[c, j, s]
                    bt_arr[c, jj] = epsilon_t[c, j, s] / dtime_ch4 + 0.5 / dzj * dm1_zm1[c, j]
                    ct_arr[c, jj] = 0.0
                    rt_arr[c, jj] = epsilon_t[c, j, s] / dtime_ch4 * conc_rel[c, jj] +
                                    0.5 / dzj * (-dm1_zm1[c, j] * (conc_rel[c, jj] - conc_rel[c, jj-1] * k_h_cc_arr[c, j, s])) +
                                    source[c, j, s]
                else  # j == nlevsoi
                    dzj = dz[c, j]
                    at_arr[c, jj] = -0.5 / dzj * dm1_zm1[c, j]
                    bt_arr[c, jj] = epsilon_t[c, j, s] / dtime_ch4 + 0.5 / dzj * dm1_zm1[c, j]
                    ct_arr[c, jj] = 0.0
                    rt_arr[c, jj] = epsilon_t[c, j, s] / dtime_ch4 * conc_rel[c, jj] +
                                    0.5 / dzj * (-dm1_zm1[c, j] * (conc_rel[c, jj] - conc_rel[c, jj-1])) +
                                    source[c, j, s]
                end
            end
        end

        # Solve tridiagonal system for each column
        nlevs = nlevsoi + 1  # 0:nlevsoi → 1:nlevsoi+1
        for c in eachindex(mask_soil)
            mask_soil[c] || continue
            tridiagonal_solve!(
                view(conc_rel, c, :),
                view(at_arr, c, :),
                view(bt_arr, c, :),
                view(ct_arr, c, :),
                view(rt_arr, c, :),
                1, nlevs)
        end

        if s == 1  # CH4
            # Surface flux
            for c in eachindex(mask_soil)
                mask_soil[c] || continue
                g = col_gridcell[c]
                if jwt[c] != 0
                    ch4_surf_diff[c] = dm1_zm1[c, 1] * ((conc_rel[c, 2] + conc_ch4_rel_old[c, 2]) / 2.0 - c_atm[g, s])
                    ch4_surf_ebul[c] = 0.0
                else
                    ch4_surf_diff[c] = dm1_zm1[c, 1] * ((conc_rel[c, 2] + conc_ch4_rel_old[c, 2]) / 2.0 -
                                                          c_atm[g, s] * k_h_cc_arr[c, 1, s])
                    ch4_surf_ebul[c] = ch4_ebul_total[c]
                end
            end

            # Ensure non-negative concentrations
            for j in 1:nlevsoi
                for c in eachindex(mask_soil)
                    mask_soil[c] || continue
                    jj = j + 1
                    if conc_rel[c, jj] < 0.0
                        deficit = -conc_rel[c, jj] * epsilon_t[c, j, 1] * dz[c, j]
                        conc_rel[c, jj] = 0.0
                        ch4_surf_diff[c] -= deficit / dtime_ch4
                    end
                end
            end

            # Copy back to conc_ch4_rel for balance check
            conc_ch4_rel .= conc_rel

        elseif s == 2  # O2
            for j in 1:nlevsoi
                for c in eachindex(mask_soil)
                    mask_soil[c] || continue
                    g = col_gridcell[c]
                    jj = j + 1
                    conc_rel[c, jj] = smooth_max(conc_rel[c, jj], 1.0e-12)
                    conc_rel[c, jj] = smooth_min(conc_rel[c, jj], c_atm[g, 2] / epsilon_t[c, j, 2])
                end
            end
            conc_o2_rel .= conc_rel
        end
    end  # species loop

    # Update absolute concentrations
    for j in 1:nlevsoi
        for c in eachindex(mask_soil)
            mask_soil[c] || continue
            conc_ch4[c, j] = conc_ch4_rel[c, j+1] * epsilon_t[c, j, 1]
            conc_o2[c, j] = conc_o2_rel[c, j+1] * epsilon_t[c, j, 2]
        end
    end

    # Balance check
    for c in eachindex(mask_soil)
        mask_soil[c] || continue
        errch4 = 0.0
        for j in 1:nlevsoi
            errch4 += (conc_ch4[c, j] - conc_ch4_bef[c, j]) * dz[c, j]
            errch4 -= ch4_prod_depth[c, j] * dz[c, j] * dtime
            errch4 += ch4_oxid_depth[c, j] * dz[c, j] * dtime
        end
        errch4 += (ch4_surf_aere[c] + ch4_surf_ebul[c] + ch4_surf_diff[c]) * dtime
        if abs(errch4) < 1.0e-8
            ch4_surf_diff[c] -= errch4 / dtime
        end
        ch4.grnd_ch4_cond_col[c] = spec_grnd_cond[c, 1]
    end
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
              zi::Matrix{<:Real},
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
