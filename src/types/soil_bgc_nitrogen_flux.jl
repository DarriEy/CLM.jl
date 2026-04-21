# ==========================================================================
# Ported from: src/soilbiogeochem/SoilBiogeochemNitrogenFluxType.F90
# Soil biogeochemistry nitrogen flux data type allocation and initialization
# ==========================================================================

"""
    SoilBiogeochemNitrogenFluxData

Soil biogeochemistry nitrogen flux data structure. Holds decomposition N fluxes,
mineralization/immobilization, nitrification/denitrification, leaching, and
related diagnostic fields at column level.

Ported from `soilbiogeochem_nitrogenflux_type` in `SoilBiogeochemNitrogenFluxType.F90`.
"""
Base.@kwdef mutable struct SoilBiogeochemNitrogenFluxData{FT<:Real}
    # --- Deposition fluxes (1D col) ---
    ndep_to_sminn_col                        ::Vector{FT} = Float64[]  # (gN/m2/s) atmospheric N deposition to soil mineral N
    nfix_to_sminn_col                        ::Vector{FT} = Float64[]  # (gN/m2/s) symbiotic/asymbiotic N fixation
    ffix_to_sminn_col                        ::Vector{FT} = Float64[]  # (gN/m2/s) free living N fixation
    fert_to_sminn_col                        ::Vector{FT} = Float64[]  # (gN/m2/s) fertilizer N to soil mineral N
    soyfixn_to_sminn_col                     ::Vector{FT} = Float64[]  # (gN/m2/s) soybean fixation

    # --- Decomposition cascade fluxes ---
    decomp_cascade_ntransfer_vr_col          ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)  # (gN/m3/s) vr N transfer along decomp cascade
    decomp_cascade_ntransfer_col             ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (gN/m2/s) vert-int N transfer along decomp cascade
    decomp_cascade_sminn_flux_vr_col         ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)  # (gN/m3/s) vr mineral N flux along decomp cascade
    decomp_cascade_sminn_flux_col            ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (gN/m2/s) vert-int mineral N flux along decomp cascade

    # --- Immobilization / mineralization fluxes ---
    potential_immob_vr_col                   ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (gN/m3/s) potential N immobilization vr
    potential_immob_col                      ::Vector{FT} = Float64[]                        # (gN/m2/s) potential N immobilization
    actual_immob_vr_col                      ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (gN/m3/s) actual N immobilization vr
    actual_immob_col                         ::Vector{FT} = Float64[]                        # (gN/m2/s) actual N immobilization
    sminn_to_plant_vr_col                    ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (gN/m3/s) plant uptake of soil mineral N vr
    sminn_to_plant_col                       ::Vector{FT} = Float64[]                        # (gN/m2/s) plant uptake of soil mineral N
    supplement_to_sminn_vr_col               ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (gN/m3/s) supplemental N supply vr
    supplement_to_sminn_col                  ::Vector{FT} = Float64[]                        # (gN/m2/s) supplemental N supply
    gross_nmin_vr_col                        ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (gN/m3/s) gross N mineralization vr
    gross_nmin_col                           ::Vector{FT} = Float64[]                        # (gN/m2/s) gross N mineralization
    net_nmin_vr_col                          ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (gN/m3/s) net N mineralization vr
    net_nmin_col                             ::Vector{FT} = Float64[]                        # (gN/m2/s) net N mineralization
    sminn_to_plant_fun_col                   ::Vector{FT} = Float64[]                        # (gN/m2/s) total soil N uptake of FUN
    sminn_to_plant_fun_vr_col                ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (gN/m2/s) total layer soil N uptake of FUN
    sminn_to_plant_fun_no3_vr_col            ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (gN/m2/s) total layer no3 uptake of FUN
    sminn_to_plant_fun_nh4_vr_col            ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (gN/m2/s) total layer nh4 uptake of FUN

    # --- Nitrification / denitrification fluxes ---
    f_nit_vr_col                             ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (gN/m3/s) soil nitrification flux vr
    f_denit_vr_col                           ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (gN/m3/s) soil denitrification flux vr
    f_nit_col                                ::Vector{FT} = Float64[]                        # (gN/m2/s) soil nitrification flux
    f_denit_col                              ::Vector{FT} = Float64[]                        # (gN/m2/s) soil denitrification flux
    pot_f_nit_vr_col                         ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (gN/m3/s) potential nitrification vr
    pot_f_denit_vr_col                       ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (gN/m3/s) potential denitrification vr
    pot_f_nit_col                            ::Vector{FT} = Float64[]                        # (gN/m2/s) potential nitrification
    pot_f_denit_col                          ::Vector{FT} = Float64[]                        # (gN/m2/s) potential denitrification
    n2_n2o_ratio_denit_vr_col                ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (gN/gN) N2:N2O ratio from denitrification
    f_n2o_denit_vr_col                       ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (gN/m3/s) N2O from denitrification vr
    f_n2o_denit_col                          ::Vector{FT} = Float64[]                        # (gN/m2/s) N2O from denitrification
    f_n2o_nit_vr_col                         ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (gN/m3/s) N2O from nitrification vr
    f_n2o_nit_col                            ::Vector{FT} = Float64[]                        # (gN/m2/s) N2O from nitrification

    # --- NO3/NH4 immobilization / uptake fluxes (nitrif_denitrif) ---
    actual_immob_no3_vr_col                  ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (gN/m3/s) actual immobilization of NO3 vr
    actual_immob_nh4_vr_col                  ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (gN/m3/s) actual immobilization of NH4 vr
    smin_no3_to_plant_vr_col                 ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (gN/m3/s) plant uptake of NO3 vr
    smin_nh4_to_plant_vr_col                 ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (gN/m3/s) plant uptake of NH4 vr
    actual_immob_no3_col                     ::Vector{FT} = Float64[]                        # (gN/m2/s) actual immobilization of NO3
    actual_immob_nh4_col                     ::Vector{FT} = Float64[]                        # (gN/m2/s) actual immobilization of NH4
    smin_no3_to_plant_col                    ::Vector{FT} = Float64[]                        # (gN/m2/s) plant uptake of NO3
    smin_nh4_to_plant_col                    ::Vector{FT} = Float64[]                        # (gN/m2/s) plant uptake of NH4

    # --- NO3 leaching / runoff fluxes (nitrif_denitrif) ---
    smin_no3_leached_vr_col                  ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (gN/m3/s) NO3 loss to leaching vr
    smin_no3_leached_col                     ::Vector{FT} = Float64[]                        # (gN/m2/s) NO3 loss to leaching
    smin_no3_runoff_vr_col                   ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (gN/m3/s) NO3 loss to runoff vr
    smin_no3_runoff_col                      ::Vector{FT} = Float64[]                        # (gN/m2/s) NO3 loss to runoff

    # --- Nitrif/denitrif diagnostic quantities ---
    smin_no3_massdens_vr_col                 ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (ugN/g soil) soil nitrate concentration
    soil_bulkdensity_col                     ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (kg/m3) bulk density of soil
    k_nitr_t_vr_col                          ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    k_nitr_ph_vr_col                         ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    k_nitr_h2o_vr_col                        ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    k_nitr_vr_col                            ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    wfps_vr_col                              ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    fmax_denit_carbonsubstrate_vr_col        ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    fmax_denit_nitrate_vr_col                ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    f_denit_base_vr_col                      ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    diffus_col                               ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (m2/s) diffusivity
    ratio_k1_col                             ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    ratio_no3_co2_col                        ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    soil_co2_prod_col                        ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    fr_WFPS_col                              ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    r_psi_col                                ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    anaerobic_frac_col                       ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)

    # --- Non-nitrif_denitrif denitrification fluxes ---
    sminn_to_denit_decomp_cascade_vr_col     ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)  # (gN/m3/s) denitrif along decomp cascade vr
    sminn_to_denit_decomp_cascade_col        ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (gN/m2/s) denitrif along decomp cascade
    sminn_to_denit_excess_vr_col             ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (gN/m3/s) excess denitrification vr
    sminn_to_denit_excess_col                ::Vector{FT} = Float64[]                        # (gN/m2/s) excess denitrification

    # --- Non-nitrif_denitrif leaching fluxes ---
    sminn_leached_vr_col                     ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (gN/m3/s) sminn leaching vr
    sminn_leached_col                        ::Vector{FT} = Float64[]                        # (gN/m2/s) sminn leaching

    # --- Summary (diagnostic) flux variables ---
    denit_col                                ::Vector{FT} = Float64[]  # (gN/m2/s) total denitrification
    ninputs_col                              ::Vector{FT} = Float64[]  # (gN/m2/s) column N inputs
    noutputs_col                             ::Vector{FT} = Float64[]  # (gN/m2/s) column N outputs
    som_n_leached_col                        ::Vector{FT} = Float64[]  # (gN/m2/s) total SOM N loss from vert transport
    decomp_npools_leached_col                ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # (gN/m2/s) N loss from vert transport per pool
    decomp_npools_transport_tendency_col     ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)  # (gN/m3/s) N tendency due to vert transport

    # --- All n pools involved in decomposition ---
    decomp_npools_sourcesink_col             ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)  # (gN/m3) change in decomposing N pools
    fates_litter_flux                        ::Vector{FT} = Float64[]  # (gN/m2/s) total litter flux from FATES

    # --- Matrix-CN fields ---
    NE_AKallsoiln                            ::Int = 0
    RI_AKallsoiln                            ::Vector{Int} = Int[]
    CI_AKallsoiln                            ::Vector{Int} = Int[]
    RI_na                                    ::Vector{Int} = Int[]
    CI_na                                    ::Vector{Int} = Int[]
end

# ---------------------------------------------------------------------------
# Helper constructors (reuse from carbon state if already defined)
# ---------------------------------------------------------------------------
if !@isdefined(nanvec)
    nanvec(n) = fill(NaN, n)
end
if !@isdefined(nanmat)
    nanmat(r, c) = fill(NaN, r, c)
end
if !@isdefined(nan3d)
    nan3d(a, b, c) = fill(NaN, a, b, c)
end

"""
    soil_bgc_nitrogen_flux_init!(nf, nc, nlevdecomp_full, ndecomp_pools,
                                  ndecomp_cascade_transitions;
                                  nlevdecomp=nlevdecomp_full,
                                  ndecomp_cascade_outtransitions=0,
                                  ndecomp_pools_vr=ndecomp_pools*nlevdecomp,
                                  use_soil_matrixcn=false,
                                  use_fates=false,
                                  Ntri_setup=0)

Allocate all fields of `SoilBiogeochemNitrogenFluxData`.
Corresponds to `InitAllocate` in the Fortran source.
"""
function soil_bgc_nitrogen_flux_init!(nf::SoilBiogeochemNitrogenFluxData,
                                       nc::Int,
                                       nlevdecomp_full::Int,
                                       ndecomp_pools::Int,
                                       ndecomp_cascade_transitions::Int;
                                       nlevdecomp::Int=nlevdecomp_full,
                                       ndecomp_cascade_outtransitions::Int=0,
                                       ndecomp_pools_vr::Int=ndecomp_pools * nlevdecomp,
                                       use_soil_matrixcn::Bool=false,
                                       use_fates::Bool=false,
                                       Ntri_setup::Int=0)

    # --- Deposition fluxes (1D col) ---
    nf.ndep_to_sminn_col    = nanvec(nc)
    nf.nfix_to_sminn_col    = nanvec(nc)
    nf.ffix_to_sminn_col    = nanvec(nc)
    nf.fert_to_sminn_col    = nanvec(nc)
    nf.soyfixn_to_sminn_col = nanvec(nc)

    # --- Column-level 1D summary ---
    nf.sminn_to_plant_col     = nanvec(nc)
    nf.potential_immob_col    = nanvec(nc)
    nf.actual_immob_col       = nanvec(nc)
    nf.gross_nmin_col         = nanvec(nc)
    nf.net_nmin_col           = nanvec(nc)
    nf.denit_col              = nanvec(nc)
    nf.supplement_to_sminn_col = nanvec(nc)
    nf.ninputs_col            = nanvec(nc)
    nf.noutputs_col           = nanvec(nc)
    nf.som_n_leached_col      = nanvec(nc)
    nf.sminn_to_plant_fun_col = nanvec(nc)

    # --- Vertically-resolved (col × nlevdecomp_full) ---
    nf.potential_immob_vr_col       = nanmat(nc, nlevdecomp_full)
    nf.actual_immob_vr_col          = nanmat(nc, nlevdecomp_full)
    nf.sminn_to_plant_vr_col        = nanmat(nc, nlevdecomp_full)
    nf.supplement_to_sminn_vr_col   = nanmat(nc, nlevdecomp_full)
    nf.gross_nmin_vr_col            = nanmat(nc, nlevdecomp_full)
    nf.net_nmin_vr_col              = nanmat(nc, nlevdecomp_full)
    nf.sminn_to_plant_fun_vr_col    = nanmat(nc, nlevdecomp_full)
    nf.sminn_to_plant_fun_no3_vr_col = nanmat(nc, nlevdecomp_full)
    nf.sminn_to_plant_fun_nh4_vr_col = nanmat(nc, nlevdecomp_full)
    nf.r_psi_col                    = nanmat(nc, nlevdecomp_full)
    nf.anaerobic_frac_col           = nanmat(nc, nlevdecomp_full)

    # --- Nitrification / denitrification (col × nlevdecomp_full) ---
    nf.f_nit_vr_col                      = nanmat(nc, nlevdecomp_full)
    nf.f_denit_vr_col                    = nanmat(nc, nlevdecomp_full)
    nf.smin_no3_leached_vr_col           = nanmat(nc, nlevdecomp_full)
    nf.smin_no3_runoff_vr_col            = nanmat(nc, nlevdecomp_full)
    nf.pot_f_nit_vr_col                  = nanmat(nc, nlevdecomp_full)
    nf.pot_f_denit_vr_col                = nanmat(nc, nlevdecomp_full)
    nf.actual_immob_no3_vr_col           = nanmat(nc, nlevdecomp_full)
    nf.actual_immob_nh4_vr_col           = nanmat(nc, nlevdecomp_full)
    nf.smin_no3_to_plant_vr_col          = nanmat(nc, nlevdecomp_full)
    nf.smin_nh4_to_plant_vr_col          = nanmat(nc, nlevdecomp_full)
    nf.n2_n2o_ratio_denit_vr_col         = nanmat(nc, nlevdecomp_full)
    nf.f_n2o_denit_vr_col                = nanmat(nc, nlevdecomp_full)
    nf.f_n2o_nit_vr_col                  = nanmat(nc, nlevdecomp_full)
    nf.smin_no3_massdens_vr_col          = nanmat(nc, nlevdecomp_full)
    nf.soil_bulkdensity_col              = nanmat(nc, nlevdecomp_full)
    nf.k_nitr_t_vr_col                   = nanmat(nc, nlevdecomp_full)
    nf.k_nitr_ph_vr_col                  = nanmat(nc, nlevdecomp_full)
    nf.k_nitr_h2o_vr_col                 = nanmat(nc, nlevdecomp_full)
    nf.k_nitr_vr_col                     = nanmat(nc, nlevdecomp_full)
    nf.wfps_vr_col                       = nanmat(nc, nlevdecomp_full)
    nf.fmax_denit_carbonsubstrate_vr_col = nanmat(nc, nlevdecomp_full)
    nf.fmax_denit_nitrate_vr_col         = nanmat(nc, nlevdecomp_full)
    nf.f_denit_base_vr_col               = nanmat(nc, nlevdecomp_full)
    nf.diffus_col                        = nanmat(nc, nlevdecomp_full)
    nf.ratio_k1_col                      = nanmat(nc, nlevdecomp_full)
    nf.ratio_no3_co2_col                 = nanmat(nc, nlevdecomp_full)
    nf.soil_co2_prod_col                 = nanmat(nc, nlevdecomp_full)
    nf.fr_WFPS_col                       = nanmat(nc, nlevdecomp_full)

    # --- Nitrif/denitrif 1D col ---
    nf.smin_no3_leached_col = nanvec(nc)
    nf.smin_no3_runoff_col  = nanvec(nc)
    nf.pot_f_nit_col        = nanvec(nc)
    nf.pot_f_denit_col      = nanvec(nc)
    nf.f_nit_col            = nanvec(nc)
    nf.f_denit_col          = nanvec(nc)
    nf.f_n2o_denit_col      = nanvec(nc)
    nf.f_n2o_nit_col        = nanvec(nc)
    nf.actual_immob_no3_col = nanvec(nc)
    nf.actual_immob_nh4_col = nanvec(nc)
    nf.smin_no3_to_plant_col = nanvec(nc)
    nf.smin_nh4_to_plant_col = nanvec(nc)

    # --- Decomposition cascade fluxes (3D: col × nlev × ntrans) ---
    nf.decomp_cascade_ntransfer_vr_col  = nan3d(nc, nlevdecomp_full, ndecomp_cascade_transitions)
    nf.decomp_cascade_sminn_flux_vr_col = nan3d(nc, nlevdecomp_full, ndecomp_cascade_transitions)

    # --- Decomposition cascade fluxes (2D: col × ntrans) ---
    nf.decomp_cascade_ntransfer_col     = nanmat(nc, ndecomp_cascade_transitions)
    nf.decomp_cascade_sminn_flux_col    = nanmat(nc, ndecomp_cascade_transitions)

    # --- Non-nitrif_denitrif denitrification ---
    nf.sminn_to_denit_decomp_cascade_vr_col = nan3d(nc, nlevdecomp_full, ndecomp_cascade_transitions)
    nf.sminn_to_denit_decomp_cascade_col    = nanmat(nc, ndecomp_cascade_transitions)
    nf.sminn_to_denit_excess_vr_col         = nanmat(nc, nlevdecomp_full)
    nf.sminn_to_denit_excess_col            = nanvec(nc)
    nf.sminn_leached_vr_col                 = nanmat(nc, nlevdecomp_full)
    nf.sminn_leached_col                    = nanvec(nc)

    # --- Pool-level diagnostics ---
    nf.decomp_npools_leached_col            = nanmat(nc, ndecomp_pools)
    nf.decomp_npools_transport_tendency_col  = nan3d(nc, nlevdecomp_full, ndecomp_pools)

    # --- All n pools in decomposition (3D: col × nlev × npools) ---
    nf.decomp_npools_sourcesink_col = nan3d(nc, nlevdecomp_full, ndecomp_pools)

    # --- FATES litter flux ---
    if use_fates
        nf.fates_litter_flux = nanvec(nc)
    else
        nf.fates_litter_flux = fill(NaN, 1)
    end

    # --- Matrix-CN fields ---
    if use_soil_matrixcn
        Ntrans = (ndecomp_cascade_transitions - ndecomp_cascade_outtransitions) * nlevdecomp
        nf.NE_AKallsoiln = (Ntrans + nlevdecomp * ndecomp_pools) +
                           (Ntrans + Ntri_setup + nlevdecomp) +
                           (ndecomp_pools * nlevdecomp)
        nf.RI_AKallsoiln = fill(-9999, nf.NE_AKallsoiln)
        nf.CI_AKallsoiln = fill(-9999, nf.NE_AKallsoiln)

        Ntrans_diag = (ndecomp_cascade_transitions - ndecomp_cascade_outtransitions) * nlevdecomp + ndecomp_pools_vr
        nf.RI_na = fill(-9999, Ntrans_diag)
        nf.CI_na = fill(-9999, Ntrans_diag)
    end

    return nothing
end

"""
    soil_bgc_nitrogen_flux_set_values!(nf, mask_col, value_column;
                                        nlevdecomp_full=..., nlevdecomp=...,
                                        ndecomp_pools=..., ndecomp_cascade_transitions=...,
                                        use_nitrif_denitrif=false,
                                        use_soil_matrixcn=false)

Set nitrogen flux variables for columns matching `mask_col`.
Corresponds to `SetValues` in the Fortran source.
"""
function soil_bgc_nitrogen_flux_set_values!(nf::SoilBiogeochemNitrogenFluxData,
                                             mask_col::BitVector,
                                             value_column::Real;
                                             nlevdecomp_full::Int=size(nf.potential_immob_vr_col, 2),
                                             nlevdecomp::Int=nlevdecomp_full,
                                             ndecomp_pools::Int=size(nf.decomp_npools_leached_col, 2),
                                             ndecomp_cascade_transitions::Int=size(nf.decomp_cascade_ntransfer_col, 2),
                                             use_nitrif_denitrif::Bool=false,
                                             use_soil_matrixcn::Bool=false)

    # Vertically-resolved fields
    for j in 1:nlevdecomp_full
        for i in eachindex(mask_col)
            mask_col[i] || continue

            if !use_nitrif_denitrif
                nf.sminn_to_denit_excess_vr_col[i, j]      = value_column
                nf.sminn_leached_vr_col[i, j]              = value_column
                nf.sminn_to_plant_fun_vr_col[i, j]         = value_column
            else
                nf.f_nit_vr_col[i, j]                      = value_column
                nf.f_denit_vr_col[i, j]                    = value_column
                nf.smin_no3_leached_vr_col[i, j]           = value_column
                nf.smin_no3_runoff_vr_col[i, j]            = value_column
                nf.n2_n2o_ratio_denit_vr_col[i, j]         = value_column
                nf.pot_f_nit_vr_col[i, j]                  = value_column
                nf.pot_f_denit_vr_col[i, j]                = value_column
                nf.actual_immob_no3_vr_col[i, j]           = value_column
                nf.actual_immob_nh4_vr_col[i, j]           = value_column
                nf.smin_no3_to_plant_vr_col[i, j]          = value_column
                nf.smin_nh4_to_plant_vr_col[i, j]          = value_column
                nf.f_n2o_denit_vr_col[i, j]                = value_column
                nf.f_n2o_nit_vr_col[i, j]                  = value_column

                nf.smin_no3_massdens_vr_col[i, j]          = value_column
                nf.k_nitr_t_vr_col[i, j]                   = value_column
                nf.k_nitr_ph_vr_col[i, j]                  = value_column
                nf.k_nitr_h2o_vr_col[i, j]                 = value_column
                nf.k_nitr_vr_col[i, j]                     = value_column
                nf.wfps_vr_col[i, j]                       = value_column
                nf.fmax_denit_carbonsubstrate_vr_col[i, j] = value_column
                nf.fmax_denit_nitrate_vr_col[i, j]         = value_column
                nf.f_denit_base_vr_col[i, j]               = value_column

                nf.diffus_col[i, j]                        = value_column
                nf.ratio_k1_col[i, j]                      = value_column
                nf.ratio_no3_co2_col[i, j]                 = value_column
                nf.soil_co2_prod_col[i, j]                 = value_column
                nf.fr_WFPS_col[i, j]                       = value_column
                nf.soil_bulkdensity_col[i, j]              = value_column

                nf.r_psi_col[i, j]                         = value_column
                nf.anaerobic_frac_col[i, j]                = value_column
            end
            nf.potential_immob_vr_col[i, j]               = value_column
            nf.actual_immob_vr_col[i, j]                  = value_column
            nf.sminn_to_plant_vr_col[i, j]                = value_column
            nf.supplement_to_sminn_vr_col[i, j]           = value_column
            nf.gross_nmin_vr_col[i, j]                    = value_column
            nf.net_nmin_vr_col[i, j]                      = value_column
            nf.sminn_to_plant_fun_no3_vr_col[i, j]        = value_column
            nf.sminn_to_plant_fun_nh4_vr_col[i, j]        = value_column
        end
    end

    # Column-level 1D
    for i in eachindex(mask_col)
        mask_col[i] || continue

        nf.ndep_to_sminn_col[i]        = value_column
        nf.nfix_to_sminn_col[i]        = value_column
        nf.ffix_to_sminn_col[i]        = value_column
        nf.fert_to_sminn_col[i]        = value_column
        nf.soyfixn_to_sminn_col[i]     = value_column
        nf.potential_immob_col[i]      = value_column
        nf.actual_immob_col[i]         = value_column
        nf.sminn_to_plant_col[i]       = value_column
        nf.supplement_to_sminn_col[i]  = value_column
        nf.gross_nmin_col[i]           = value_column
        nf.net_nmin_col[i]             = value_column
        nf.denit_col[i]                = value_column
        nf.sminn_to_plant_fun_col[i]   = value_column
        if use_nitrif_denitrif
            nf.f_nit_col[i]            = value_column
            nf.pot_f_nit_col[i]        = value_column
            nf.f_denit_col[i]          = value_column
            nf.pot_f_denit_col[i]      = value_column
            nf.f_n2o_denit_col[i]      = value_column
            nf.f_n2o_nit_col[i]        = value_column
            nf.smin_no3_leached_col[i] = value_column
            nf.smin_no3_runoff_col[i]  = value_column
        else
            nf.sminn_to_denit_excess_col[i] = value_column
            nf.sminn_leached_col[i]         = value_column
        end
        nf.ninputs_col[i]              = value_column
        nf.noutputs_col[i]             = value_column
        nf.som_n_leached_col[i]        = value_column
    end

    # Pool-level fields
    for k in 1:ndecomp_pools
        for i in eachindex(mask_col)
            mask_col[i] || continue
            nf.decomp_npools_leached_col[i, k] = value_column
        end
    end

    # Pool × level fields
    for k in 1:ndecomp_pools
        for j in 1:nlevdecomp_full
            for i in eachindex(mask_col)
                mask_col[i] || continue
                nf.decomp_npools_transport_tendency_col[i, j, k] = value_column
            end
        end
    end

    # Decomposition cascade transitions (2D)
    for l in 1:ndecomp_cascade_transitions
        for i in eachindex(mask_col)
            mask_col[i] || continue
            nf.decomp_cascade_ntransfer_col[i, l]  = value_column
            nf.decomp_cascade_sminn_flux_col[i, l]  = value_column
            if !use_nitrif_denitrif
                nf.sminn_to_denit_decomp_cascade_col[i, l] = value_column
            end
        end
    end

    # Decomposition cascade transitions (3D: col × lev × trans)
    for l in 1:ndecomp_cascade_transitions
        for j in 1:nlevdecomp_full
            for i in eachindex(mask_col)
                mask_col[i] || continue
                nf.decomp_cascade_ntransfer_vr_col[i, j, l]  = value_column
                nf.decomp_cascade_sminn_flux_vr_col[i, j, l]  = value_column
                if !use_nitrif_denitrif
                    nf.sminn_to_denit_decomp_cascade_vr_col[i, j, l] = value_column
                end
            end
        end
    end

    # All decomposing N pools sourcesink
    for k in 1:ndecomp_pools
        for j in 1:nlevdecomp_full
            for i in eachindex(mask_col)
                mask_col[i] || continue
                nf.decomp_npools_sourcesink_col[i, j, k] = value_column
            end
        end
    end

    return nothing
end

"""
    soil_bgc_nitrogen_flux_init_cold!(nf, bounds_col; mask_special=nothing)

Cold-start initialization of soil N flux fields.
Corresponds to `InitCold` in the Fortran source.
Sets values to zero for special landunit columns.
"""
function soil_bgc_nitrogen_flux_init_cold!(nf::SoilBiogeochemNitrogenFluxData,
                                            bounds_col::UnitRange{Int};
                                            mask_special::Union{BitVector,Nothing}=nothing)
    if mask_special !== nothing
        soil_bgc_nitrogen_flux_set_values!(nf, mask_special, 0.0)
    end
    return nothing
end

"""
    soil_bgc_nitrogen_flux_summary!(nf, mask_bgc_soilc, bounds_col;
                                     ndecomp_cascade_transitions, nlevdecomp,
                                     ndecomp_pools, dzsoi_decomp_vals,
                                     use_nitrif_denitrif=false)

Nitrogen flux summary calculations: vertically integrate decomposition cascade
N transfer and mineral N fluxes, denitrification, leaching, and supplemental N.
Corresponds to `Summary` in the Fortran source.
"""
function soil_bgc_nitrogen_flux_summary!(nf::SoilBiogeochemNitrogenFluxData,
                                          mask_bgc_soilc::BitVector,
                                          bounds_col::UnitRange{Int};
                                          ndecomp_cascade_transitions::Int=size(nf.decomp_cascade_ntransfer_col, 2),
                                          nlevdecomp::Int=size(nf.potential_immob_vr_col, 2),
                                          ndecomp_pools::Int=size(nf.decomp_npools_leached_col, 2),
                                          dzsoi_decomp_vals::Vector{<:Real},
                                          use_nitrif_denitrif::Bool=false)

    # Initialize column-level accumulators
    for c in bounds_col
        mask_bgc_soilc[c] || continue
        nf.denit_col[c]              = 0.0
        nf.supplement_to_sminn_col[c] = 0.0
        nf.som_n_leached_col[c]      = 0.0
    end

    # Vertically integrate decomposing N cascade fluxes and mineral N fluxes
    for k in 1:ndecomp_cascade_transitions
        for j in 1:nlevdecomp
            for c in bounds_col
                mask_bgc_soilc[c] || continue

                nf.decomp_cascade_ntransfer_col[c, k] =
                    nf.decomp_cascade_ntransfer_col[c, k] +
                    nf.decomp_cascade_ntransfer_vr_col[c, j, k] * dzsoi_decomp_vals[j]

                nf.decomp_cascade_sminn_flux_col[c, k] =
                    nf.decomp_cascade_sminn_flux_col[c, k] +
                    nf.decomp_cascade_sminn_flux_vr_col[c, j, k] * dzsoi_decomp_vals[j]
            end
        end
    end

    if !use_nitrif_denitrif

        # Vertically integrate each denitrification flux
        for l in 1:ndecomp_cascade_transitions
            for j in 1:nlevdecomp
                for c in bounds_col
                    mask_bgc_soilc[c] || continue
                    nf.sminn_to_denit_decomp_cascade_col[c, l] =
                        nf.sminn_to_denit_decomp_cascade_col[c, l] +
                        nf.sminn_to_denit_decomp_cascade_vr_col[c, j, l] * dzsoi_decomp_vals[j]
                end
            end
        end

        # Vertically integrate bulk denitrification and leaching flux
        for j in 1:nlevdecomp
            for c in bounds_col
                mask_bgc_soilc[c] || continue
                nf.sminn_to_denit_excess_col[c] =
                    nf.sminn_to_denit_excess_col[c] +
                    nf.sminn_to_denit_excess_vr_col[c, j] * dzsoi_decomp_vals[j]

                nf.sminn_leached_col[c] =
                    nf.sminn_leached_col[c] +
                    nf.sminn_leached_vr_col[c, j] * dzsoi_decomp_vals[j]
            end
        end

        # Total N denitrification (DENIT)
        for l in 1:ndecomp_cascade_transitions
            for c in bounds_col
                mask_bgc_soilc[c] || continue
                nf.denit_col[c] =
                    nf.denit_col[c] +
                    nf.sminn_to_denit_decomp_cascade_col[c, l]
            end
        end

        for c in bounds_col
            mask_bgc_soilc[c] || continue
            nf.denit_col[c] =
                nf.denit_col[c] +
                nf.sminn_to_denit_excess_col[c]
        end

    else

        # Vertically integrate NO3 NH4 N2O fluxes and pools
        for j in 1:nlevdecomp
            for c in bounds_col
                mask_bgc_soilc[c] || continue

                # nitrification and denitrification fluxes
                nf.f_nit_col[c] =
                    nf.f_nit_col[c] +
                    nf.f_nit_vr_col[c, j] * dzsoi_decomp_vals[j]

                nf.f_denit_col[c] =
                    nf.f_denit_col[c] +
                    nf.f_denit_vr_col[c, j] * dzsoi_decomp_vals[j]

                nf.pot_f_nit_col[c] =
                    nf.pot_f_nit_col[c] +
                    nf.pot_f_nit_vr_col[c, j] * dzsoi_decomp_vals[j]

                nf.pot_f_denit_col[c] =
                    nf.pot_f_denit_col[c] +
                    nf.pot_f_denit_vr_col[c, j] * dzsoi_decomp_vals[j]

                nf.f_n2o_nit_col[c] =
                    nf.f_n2o_nit_col[c] +
                    nf.f_n2o_nit_vr_col[c, j] * dzsoi_decomp_vals[j]

                nf.f_n2o_denit_col[c] =
                    nf.f_n2o_denit_col[c] +
                    nf.f_n2o_denit_vr_col[c, j] * dzsoi_decomp_vals[j]

                # leaching/runoff flux
                nf.smin_no3_leached_col[c] =
                    nf.smin_no3_leached_col[c] +
                    nf.smin_no3_leached_vr_col[c, j] * dzsoi_decomp_vals[j]

                nf.smin_no3_runoff_col[c] =
                    nf.smin_no3_runoff_col[c] +
                    nf.smin_no3_runoff_vr_col[c, j] * dzsoi_decomp_vals[j]
            end
        end

        for c in bounds_col
            mask_bgc_soilc[c] || continue
            nf.denit_col[c] = nf.f_denit_col[c]
        end

    end

    # Supplementary N supply
    for j in 1:nlevdecomp
        for c in bounds_col
            mask_bgc_soilc[c] || continue
            nf.supplement_to_sminn_col[c] =
                nf.supplement_to_sminn_col[c] +
                nf.supplement_to_sminn_vr_col[c, j] * dzsoi_decomp_vals[j]
        end
    end

    # Add up all vertical transport tendency terms and calculate total SOM N leaching loss
    for l in 1:ndecomp_pools
        for c in bounds_col
            mask_bgc_soilc[c] || continue
            nf.decomp_npools_leached_col[c, l] = 0.0
        end

        for j in 1:nlevdecomp
            for c in bounds_col
                mask_bgc_soilc[c] || continue
                nf.decomp_npools_leached_col[c, l] =
                    nf.decomp_npools_leached_col[c, l] +
                    nf.decomp_npools_transport_tendency_col[c, j, l] * dzsoi_decomp_vals[j]
            end
        end

        for c in bounds_col
            mask_bgc_soilc[c] || continue
            nf.som_n_leached_col[c] =
                nf.som_n_leached_col[c] +
                nf.decomp_npools_leached_col[c, l]
        end
    end

    return nothing
end

"""
    soil_bgc_nitrogen_flux_init_history!(nf, bounds_col)

Stub for history field registration (no-op in Julia port).
Corresponds to `InitHistory` in the Fortran source.
"""
function soil_bgc_nitrogen_flux_init_history!(nf::SoilBiogeochemNitrogenFluxData,
                                               bounds_col::UnitRange{Int})
    return nothing
end
