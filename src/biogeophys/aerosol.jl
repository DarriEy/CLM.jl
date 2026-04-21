# ==========================================================================
# Ported from: src/biogeophys/AerosolMod.F90 (803 lines)
# Aerosol mass tracking and deposition fluxes for snow-impurity radiative
# transfer (SNICAR).
#
# Public functions:
#   aerosol_init!          — Allocate and initialize AerosolData
#   aerosol_init_cold!     — Cold-start initialization of aerosol masses
#   aerosol_reset!         — Reset aerosol masses for a single column
#   aerosol_reset_filter!  — Reset aerosol masses for masked columns
#   aerosol_masses!        — Compute column-integrated masses & concentrations
#   aerosol_fluxes!        — Compute aerosol deposition fluxes into snowpack
#   aerosol_restart!       — Restart read/write (stub)
#   aerosol_init_history!  — History field registration (stub)
# ==========================================================================

# --------------------------------------------------------------------------
# Data structure
# --------------------------------------------------------------------------

"""
    AerosolData

Aerosol mass and deposition flux data for snow impurities (SNICAR).

2D arrays are dimensioned (nc, nlevsno) where the second index maps
Fortran snow-layer index j ∈ [-nlevsno+1, 0] to Julia index j + nlevsno.

Ported from `aerosol_type` in `AerosolMod.F90`.
"""
Base.@kwdef mutable struct AerosolData{FT<:Real}
    # --- Mass of aerosol species in snow (col, lyr) [kg] ---
    mss_bcpho_col   ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # hydrophobic BC
    mss_bcphi_col   ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # hydrophilic BC
    mss_bctot_col   ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # total BC (pho+phi)
    mss_bc_col_col  ::Vector{FT} = Float64[]                      # column-integrated BC
    mss_bc_top_col  ::Vector{FT} = Float64[]                      # top-layer BC

    mss_ocpho_col   ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # hydrophobic OC
    mss_ocphi_col   ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # hydrophilic OC
    mss_octot_col   ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # total OC (pho+phi)
    mss_oc_col_col  ::Vector{FT} = Float64[]                      # column-integrated OC
    mss_oc_top_col  ::Vector{FT} = Float64[]                      # top-layer OC

    mss_dst1_col    ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # dust species 1
    mss_dst2_col    ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # dust species 2
    mss_dst3_col    ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # dust species 3
    mss_dst4_col    ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # dust species 4
    mss_dsttot_col  ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # total dust
    mss_dst_col_col ::Vector{FT} = Float64[]                      # column-integrated dust
    mss_dst_top_col ::Vector{FT} = Float64[]                      # top-layer dust

    # --- Mass concentrations in snow (col, lyr) [kg/kg] ---
    mss_cnc_bcphi_col ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    mss_cnc_bcpho_col ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    mss_cnc_ocphi_col ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    mss_cnc_ocpho_col ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    mss_cnc_dst1_col  ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    mss_cnc_dst2_col  ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    mss_cnc_dst3_col  ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    mss_cnc_dst4_col  ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)

    # --- Deposition fluxes (col) [kg/m2/s] ---
    flx_dst_dep_dry1_col ::Vector{FT} = Float64[]
    flx_dst_dep_wet1_col ::Vector{FT} = Float64[]
    flx_dst_dep_dry2_col ::Vector{FT} = Float64[]
    flx_dst_dep_wet2_col ::Vector{FT} = Float64[]
    flx_dst_dep_dry3_col ::Vector{FT} = Float64[]
    flx_dst_dep_wet3_col ::Vector{FT} = Float64[]
    flx_dst_dep_dry4_col ::Vector{FT} = Float64[]
    flx_dst_dep_wet4_col ::Vector{FT} = Float64[]
    flx_dst_dep_col      ::Vector{FT} = Float64[]

    flx_bc_dep_dry_col   ::Vector{FT} = Float64[]
    flx_bc_dep_wet_col   ::Vector{FT} = Float64[]
    flx_bc_dep_pho_col   ::Vector{FT} = Float64[]
    flx_bc_dep_phi_col   ::Vector{FT} = Float64[]
    flx_bc_dep_col       ::Vector{FT} = Float64[]

    flx_oc_dep_dry_col   ::Vector{FT} = Float64[]
    flx_oc_dep_wet_col   ::Vector{FT} = Float64[]
    flx_oc_dep_pho_col   ::Vector{FT} = Float64[]
    flx_oc_dep_phi_col   ::Vector{FT} = Float64[]
    flx_oc_dep_col       ::Vector{FT} = Float64[]
end

# --------------------------------------------------------------------------
# Initialization
# --------------------------------------------------------------------------

"""
    aerosol_init!(aer, nc)

Allocate all fields of an `AerosolData` instance for `nc` columns.

Ported from `InitAllocate` in `AerosolMod.F90`.
"""
function aerosol_init!(aer::AerosolData{FT}, nc::Int) where {FT}
    nlevsno = varpar.nlevsno

    # Deposition fluxes — 1D (col)
    aer.flx_dst_dep_dry1_col = fill(FT(NaN), nc)
    aer.flx_dst_dep_wet1_col = fill(FT(NaN), nc)
    aer.flx_dst_dep_dry2_col = fill(FT(NaN), nc)
    aer.flx_dst_dep_wet2_col = fill(FT(NaN), nc)
    aer.flx_dst_dep_dry3_col = fill(FT(NaN), nc)
    aer.flx_dst_dep_wet3_col = fill(FT(NaN), nc)
    aer.flx_dst_dep_dry4_col = fill(FT(NaN), nc)
    aer.flx_dst_dep_wet4_col = fill(FT(NaN), nc)
    aer.flx_dst_dep_col      = fill(FT(NaN), nc)

    aer.flx_bc_dep_dry_col   = fill(FT(NaN), nc)
    aer.flx_bc_dep_wet_col   = fill(FT(NaN), nc)
    aer.flx_bc_dep_pho_col   = fill(FT(NaN), nc)
    aer.flx_bc_dep_phi_col   = fill(FT(NaN), nc)
    aer.flx_bc_dep_col       = fill(FT(NaN), nc)

    aer.flx_oc_dep_dry_col   = fill(FT(NaN), nc)
    aer.flx_oc_dep_wet_col   = fill(FT(NaN), nc)
    aer.flx_oc_dep_pho_col   = fill(FT(NaN), nc)
    aer.flx_oc_dep_phi_col   = fill(FT(NaN), nc)
    aer.flx_oc_dep_col       = fill(FT(NaN), nc)

    # Mass in snow layers — 2D (col, lyr)  Fortran (-nlevsno+1:0) → Julia 1:nlevsno
    aer.mss_bcpho_col   = fill(FT(NaN), nc, nlevsno)
    aer.mss_bcphi_col   = fill(FT(NaN), nc, nlevsno)
    aer.mss_bctot_col   = fill(FT(NaN), nc, nlevsno)
    aer.mss_bc_col_col  = fill(FT(NaN), nc)
    aer.mss_bc_top_col  = fill(FT(NaN), nc)

    aer.mss_ocpho_col   = fill(FT(NaN), nc, nlevsno)
    aer.mss_ocphi_col   = fill(FT(NaN), nc, nlevsno)
    aer.mss_octot_col   = fill(FT(NaN), nc, nlevsno)
    aer.mss_oc_col_col  = fill(FT(NaN), nc)
    aer.mss_oc_top_col  = fill(FT(NaN), nc)

    aer.mss_dst1_col    = fill(FT(NaN), nc, nlevsno)
    aer.mss_dst2_col    = fill(FT(NaN), nc, nlevsno)
    aer.mss_dst3_col    = fill(FT(NaN), nc, nlevsno)
    aer.mss_dst4_col    = fill(FT(NaN), nc, nlevsno)
    aer.mss_dsttot_col  = fill(FT(NaN), nc, nlevsno)
    aer.mss_dst_col_col = fill(FT(NaN), nc)
    aer.mss_dst_top_col = fill(FT(NaN), nc)

    # Mass concentrations — 2D (col, lyr)
    aer.mss_cnc_bcphi_col = fill(FT(NaN), nc, nlevsno)
    aer.mss_cnc_bcpho_col = fill(FT(NaN), nc, nlevsno)
    aer.mss_cnc_ocphi_col = fill(FT(NaN), nc, nlevsno)
    aer.mss_cnc_ocpho_col = fill(FT(NaN), nc, nlevsno)
    aer.mss_cnc_dst1_col  = fill(FT(NaN), nc, nlevsno)
    aer.mss_cnc_dst2_col  = fill(FT(NaN), nc, nlevsno)
    aer.mss_cnc_dst3_col  = fill(FT(NaN), nc, nlevsno)
    aer.mss_cnc_dst4_col  = fill(FT(NaN), nc, nlevsno)

    return nothing
end

"""
    aerosol_clean!(aer)

Deallocate (reset to empty) all fields of an `AerosolData` instance.
"""
function aerosol_clean!(aer::AerosolData{FT}) where {FT}
    aer.flx_dst_dep_dry1_col = FT[]
    aer.flx_dst_dep_wet1_col = FT[]
    aer.flx_dst_dep_dry2_col = FT[]
    aer.flx_dst_dep_wet2_col = FT[]
    aer.flx_dst_dep_dry3_col = FT[]
    aer.flx_dst_dep_wet3_col = FT[]
    aer.flx_dst_dep_dry4_col = FT[]
    aer.flx_dst_dep_wet4_col = FT[]
    aer.flx_dst_dep_col      = FT[]

    aer.flx_bc_dep_dry_col   = FT[]
    aer.flx_bc_dep_wet_col   = FT[]
    aer.flx_bc_dep_pho_col   = FT[]
    aer.flx_bc_dep_phi_col   = FT[]
    aer.flx_bc_dep_col       = FT[]

    aer.flx_oc_dep_dry_col   = FT[]
    aer.flx_oc_dep_wet_col   = FT[]
    aer.flx_oc_dep_pho_col   = FT[]
    aer.flx_oc_dep_phi_col   = FT[]
    aer.flx_oc_dep_col       = FT[]

    aer.mss_bcpho_col   = Matrix{FT}(undef, 0, 0)
    aer.mss_bcphi_col   = Matrix{FT}(undef, 0, 0)
    aer.mss_bctot_col   = Matrix{FT}(undef, 0, 0)
    aer.mss_bc_col_col  = FT[]
    aer.mss_bc_top_col  = FT[]

    aer.mss_ocpho_col   = Matrix{FT}(undef, 0, 0)
    aer.mss_ocphi_col   = Matrix{FT}(undef, 0, 0)
    aer.mss_octot_col   = Matrix{FT}(undef, 0, 0)
    aer.mss_oc_col_col  = FT[]
    aer.mss_oc_top_col  = FT[]

    aer.mss_dst1_col    = Matrix{FT}(undef, 0, 0)
    aer.mss_dst2_col    = Matrix{FT}(undef, 0, 0)
    aer.mss_dst3_col    = Matrix{FT}(undef, 0, 0)
    aer.mss_dst4_col    = Matrix{FT}(undef, 0, 0)
    aer.mss_dsttot_col  = Matrix{FT}(undef, 0, 0)
    aer.mss_dst_col_col = FT[]
    aer.mss_dst_top_col = FT[]

    aer.mss_cnc_bcphi_col = Matrix{FT}(undef, 0, 0)
    aer.mss_cnc_bcpho_col = Matrix{FT}(undef, 0, 0)
    aer.mss_cnc_ocphi_col = Matrix{FT}(undef, 0, 0)
    aer.mss_cnc_ocpho_col = Matrix{FT}(undef, 0, 0)
    aer.mss_cnc_dst1_col  = Matrix{FT}(undef, 0, 0)
    aer.mss_cnc_dst2_col  = Matrix{FT}(undef, 0, 0)
    aer.mss_cnc_dst3_col  = Matrix{FT}(undef, 0, 0)
    aer.mss_cnc_dst4_col  = Matrix{FT}(undef, 0, 0)

    return nothing
end

"""
    aerosol_init_cold!(aer, bounds)

Cold-start initialization: zero all mass and concentration arrays.

Ported from `InitCold` in `AerosolMod.F90`.
"""
function aerosol_init_cold!(aer::AerosolData{FT}, bounds::UnitRange{Int}) where {FT}
    nlevsno = varpar.nlevsno

    for c in bounds
        for j in 1:nlevsno
            aer.mss_cnc_bcphi_col[c, j] = 0.0
            aer.mss_cnc_bcpho_col[c, j] = 0.0
            aer.mss_cnc_ocphi_col[c, j] = 0.0
            aer.mss_cnc_ocpho_col[c, j] = 0.0
            aer.mss_cnc_dst1_col[c, j]  = 0.0
            aer.mss_cnc_dst2_col[c, j]  = 0.0
            aer.mss_cnc_dst3_col[c, j]  = 0.0
            aer.mss_cnc_dst4_col[c, j]  = 0.0

            aer.mss_bctot_col[c, j] = 0.0
            aer.mss_bcpho_col[c, j] = 0.0
            aer.mss_bcphi_col[c, j] = 0.0

            aer.mss_octot_col[c, j] = 0.0
            aer.mss_ocpho_col[c, j] = 0.0
            aer.mss_ocphi_col[c, j] = 0.0

            aer.mss_dst1_col[c, j]   = 0.0
            aer.mss_dst2_col[c, j]   = 0.0
            aer.mss_dst3_col[c, j]   = 0.0
            aer.mss_dst4_col[c, j]   = 0.0
            aer.mss_dsttot_col[c, j] = 0.0
        end
    end

    return nothing
end

# --------------------------------------------------------------------------
# Reset
# --------------------------------------------------------------------------

"""
    aerosol_reset!(aer, c)

Reset all aerosol mass variables to zero for column `c`.

Ported from `Reset` in `AerosolMod.F90`.
"""
function aerosol_reset!(aer::AerosolData{FT}, c::Int) where {FT}
    nlevsno = varpar.nlevsno

    for j in 1:nlevsno
        aer.mss_bcpho_col[c, j] = 0.0
        aer.mss_bcphi_col[c, j] = 0.0
        aer.mss_bctot_col[c, j] = 0.0

        aer.mss_ocpho_col[c, j] = 0.0
        aer.mss_ocphi_col[c, j] = 0.0
        aer.mss_octot_col[c, j] = 0.0

        aer.mss_dst1_col[c, j]   = 0.0
        aer.mss_dst2_col[c, j]   = 0.0
        aer.mss_dst3_col[c, j]   = 0.0
        aer.mss_dst4_col[c, j]   = 0.0
        aer.mss_dsttot_col[c, j] = 0.0
    end

    aer.mss_bc_col_col[c]  = 0.0
    aer.mss_bc_top_col[c]  = 0.0
    aer.mss_oc_col_col[c]  = 0.0
    aer.mss_oc_top_col[c]  = 0.0
    aer.mss_dst_col_col[c] = 0.0
    aer.mss_dst_top_col[c] = 0.0

    return nothing
end

"""
    aerosol_reset_filter!(aer, mask, bounds)

Reset aerosol masses for all columns where `mask[c]` is true.

Ported from `ResetFilter` in `AerosolMod.F90`.
"""
function aerosol_reset_filter!(aer::AerosolData{FT}, mask::BitVector, bounds::UnitRange{Int}) where {FT}
    for c in bounds
        mask[c] || continue
        aerosol_reset!(aer, c)
    end
    return nothing
end

# --------------------------------------------------------------------------
# AerosolMasses
# --------------------------------------------------------------------------

"""
    aerosol_masses!(aer, mask_on, mask_off, bounds, snl, h2osoi_ice_col,
                    h2osoi_liq_col, h2osno_top_col, snw_rds_col)

Calculate column-integrated aerosol masses and mass concentrations for
radiative calculations and output.

Must be called after the snow filter is rebuilt in Hydrology2.

Index convention: Fortran snow layer j ∈ [-nlevsno+1, 0] maps to Julia
index j + nlevsno.  `snl[c]` is the (negative) number of snow layers,
identical to the Fortran convention.

## Arguments
- `aer`: aerosol data (modified in place)
- `mask_on`: columns with snow (replaces filter_on)
- `mask_off`: columns without snow (replaces filter_off)
- `bounds`: column index range
- `snl`: number of snow layers per column (negative)
- `h2osoi_ice_col`: ice lens per layer [kg/m2] (nc × (nlevsno+nlevmaxurbgrnd))
- `h2osoi_liq_col`: liquid water per layer [kg/m2]
- `h2osno_top_col`: top-layer snow mass [kg] (output, modified)
- `snw_rds_col`: snow grain radius [microns] (nc × nlevsno, modified)

Ported from `AerosolMasses` in `AerosolMod.F90`.
"""
function aerosol_masses!(aer::AerosolData,
                         mask_on::BitVector,
                         mask_off::BitVector,
                         bounds::UnitRange{Int},
                         snl::Vector{Int},
                         h2osoi_ice_col::Matrix{<:Real},
                         h2osoi_liq_col::Matrix{<:Real},
                         h2osno_top_col::Vector{<:Real},
                         snw_rds_col::Matrix{<:Real})
    nlevsno = varpar.nlevsno

    # --- Columns with snow ---
    for c in bounds
        mask_on[c] || continue

        # Zero column-integrated mass before summation
        aer.mss_bc_col_col[c]  = 0.0
        aer.mss_oc_col_col[c]  = 0.0
        aer.mss_dst_col_col[c] = 0.0

        for j in 1:nlevsno
            # j_fortran = j - nlevsno (Fortran snow layer index)
            # Layer mass of snow (h2osoi uses same offset)
            snowmass = h2osoi_ice_col[c, j] + h2osoi_liq_col[c, j]

            # Active snow layer check: Fortran j >= snl(c)+1
            # → Julia: j >= snl[c] + nlevsno + 1
            if j >= snl[c] + nlevsno + 1
                # BC
                aer.mss_bctot_col[c, j]     = aer.mss_bcpho_col[c, j] + aer.mss_bcphi_col[c, j]
                aer.mss_bc_col_col[c]      += aer.mss_bctot_col[c, j]
                aer.mss_cnc_bcphi_col[c, j] = aer.mss_bcphi_col[c, j] / snowmass
                aer.mss_cnc_bcpho_col[c, j] = aer.mss_bcpho_col[c, j] / snowmass

                # OC
                aer.mss_octot_col[c, j]     = aer.mss_ocpho_col[c, j] + aer.mss_ocphi_col[c, j]
                aer.mss_oc_col_col[c]      += aer.mss_octot_col[c, j]
                aer.mss_cnc_ocphi_col[c, j] = aer.mss_ocphi_col[c, j] / snowmass
                aer.mss_cnc_ocpho_col[c, j] = aer.mss_ocpho_col[c, j] / snowmass

                # Dust
                aer.mss_dsttot_col[c, j]    = aer.mss_dst1_col[c, j] + aer.mss_dst2_col[c, j] +
                                               aer.mss_dst3_col[c, j] + aer.mss_dst4_col[c, j]
                aer.mss_dst_col_col[c]     += aer.mss_dsttot_col[c, j]
                aer.mss_cnc_dst1_col[c, j]  = aer.mss_dst1_col[c, j] / snowmass
                aer.mss_cnc_dst2_col[c, j]  = aer.mss_dst2_col[c, j] / snowmass
                aer.mss_cnc_dst3_col[c, j]  = aer.mss_dst3_col[c, j] / snowmass
                aer.mss_cnc_dst4_col[c, j]  = aer.mss_dst4_col[c, j] / snowmass
            else
                # Empty (inactive) snow layer — zero everything
                snw_rds_col[c, j] = 0.0

                aer.mss_bcpho_col[c, j]     = 0.0
                aer.mss_bcphi_col[c, j]     = 0.0
                aer.mss_bctot_col[c, j]     = 0.0
                aer.mss_cnc_bcphi_col[c, j] = 0.0
                aer.mss_cnc_bcpho_col[c, j] = 0.0

                aer.mss_ocpho_col[c, j]     = 0.0
                aer.mss_ocphi_col[c, j]     = 0.0
                aer.mss_octot_col[c, j]     = 0.0
                aer.mss_cnc_ocphi_col[c, j] = 0.0
                aer.mss_cnc_ocpho_col[c, j] = 0.0

                aer.mss_dst1_col[c, j]      = 0.0
                aer.mss_dst2_col[c, j]      = 0.0
                aer.mss_dst3_col[c, j]      = 0.0
                aer.mss_dst4_col[c, j]      = 0.0
                aer.mss_dsttot_col[c, j]    = 0.0
                aer.mss_cnc_dst1_col[c, j]  = 0.0
                aer.mss_cnc_dst2_col[c, j]  = 0.0
                aer.mss_cnc_dst3_col[c, j]  = 0.0
                aer.mss_cnc_dst4_col[c, j]  = 0.0
            end
        end

        # Top-layer diagnostics
        # Fortran: snl(c)+1 → Julia: snl[c] + nlevsno + 1
        top_j = snl[c] + nlevsno + 1
        h2osno_top_col[c]      = h2osoi_ice_col[c, top_j] + h2osoi_liq_col[c, top_j]
        aer.mss_bc_top_col[c]  = aer.mss_bctot_col[c, top_j]
        aer.mss_oc_top_col[c]  = aer.mss_octot_col[c, top_j]
        aer.mss_dst_top_col[c] = aer.mss_dsttot_col[c, top_j]
    end

    # --- Columns without snow: zero all mass variables ---
    for c in bounds
        mask_off[c] || continue

        aer.mss_bc_top_col[c]  = 0.0
        aer.mss_bc_col_col[c]  = 0.0
        aer.mss_oc_top_col[c]  = 0.0
        aer.mss_oc_col_col[c]  = 0.0
        aer.mss_dst_top_col[c] = 0.0
        aer.mss_dst_col_col[c] = 0.0

        for j in 1:nlevsno
            aer.mss_bcpho_col[c, j]     = 0.0
            aer.mss_bcphi_col[c, j]     = 0.0
            aer.mss_bctot_col[c, j]     = 0.0
            aer.mss_cnc_bcphi_col[c, j] = 0.0
            aer.mss_cnc_bcpho_col[c, j] = 0.0

            aer.mss_ocpho_col[c, j]     = 0.0
            aer.mss_ocphi_col[c, j]     = 0.0
            aer.mss_octot_col[c, j]     = 0.0
            aer.mss_cnc_ocphi_col[c, j] = 0.0
            aer.mss_cnc_ocpho_col[c, j] = 0.0

            aer.mss_dst1_col[c, j]      = 0.0
            aer.mss_dst2_col[c, j]      = 0.0
            aer.mss_dst3_col[c, j]      = 0.0
            aer.mss_dst4_col[c, j]      = 0.0
            aer.mss_dsttot_col[c, j]    = 0.0
            aer.mss_cnc_dst1_col[c, j]  = 0.0
            aer.mss_cnc_dst2_col[c, j]  = 0.0
            aer.mss_cnc_dst3_col[c, j]  = 0.0
            aer.mss_cnc_dst4_col[c, j]  = 0.0
        end
    end

    return nothing
end

# --------------------------------------------------------------------------
# AerosolFluxes
# --------------------------------------------------------------------------

"""
    aerosol_fluxes!(aer, mask_snow, bounds, snl, col_gridcell, forc_aer, dtime;
                    snicar_use_aerosol=varctl.snicar_use_aerosol)

Compute aerosol deposition fluxes from the atmosphere and add deposited
aerosol mass to the top snow layer.

## Arguments
- `aer`: aerosol data (modified in place)
- `mask_snow`: columns with snow (replaces filter_snowc)
- `bounds`: column index range
- `snl`: number of snow layers per column (negative)
- `col_gridcell`: column → gridcell mapping
- `forc_aer`: aerosol forcing array (ng × 14) [kg/m2/s]
- `dtime`: model time step [s]
- `snicar_use_aerosol`: flag to include aerosol effects in snow (default from varctl)

Ported from `AerosolFluxes` in `AerosolMod.F90`.
"""
function aerosol_fluxes!(aer::AerosolData,
                         mask_snow::BitVector,
                         bounds::UnitRange{Int},
                         snl::Vector{Int},
                         col_gridcell::Vector{Int},
                         forc_aer::Matrix{<:Real},
                         dtime::Real;
                         snicar_use_aerosol::Bool = varctl.snicar_use_aerosol)
    nlevsno = varpar.nlevsno

    # Set deposition fluxes from forcing array for ALL columns
    for c in bounds
        g = col_gridcell[c]

        aer.flx_bc_dep_dry_col[c] = forc_aer[g, 1] + forc_aer[g, 2]
        aer.flx_bc_dep_wet_col[c] = forc_aer[g, 3]
        aer.flx_bc_dep_phi_col[c] = forc_aer[g, 1] + forc_aer[g, 3]
        aer.flx_bc_dep_pho_col[c] = forc_aer[g, 2]
        aer.flx_bc_dep_col[c]     = forc_aer[g, 1] + forc_aer[g, 2] + forc_aer[g, 3]

        aer.flx_oc_dep_dry_col[c] = forc_aer[g, 4] + forc_aer[g, 5]
        aer.flx_oc_dep_wet_col[c] = forc_aer[g, 6]
        aer.flx_oc_dep_phi_col[c] = forc_aer[g, 4] + forc_aer[g, 6]
        aer.flx_oc_dep_pho_col[c] = forc_aer[g, 5]
        aer.flx_oc_dep_col[c]     = forc_aer[g, 4] + forc_aer[g, 5] + forc_aer[g, 6]

        aer.flx_dst_dep_wet1_col[c] = forc_aer[g, 7]
        aer.flx_dst_dep_dry1_col[c] = forc_aer[g, 8]
        aer.flx_dst_dep_wet2_col[c] = forc_aer[g, 9]
        aer.flx_dst_dep_dry2_col[c] = forc_aer[g, 10]
        aer.flx_dst_dep_wet3_col[c] = forc_aer[g, 11]
        aer.flx_dst_dep_dry3_col[c] = forc_aer[g, 12]
        aer.flx_dst_dep_wet4_col[c] = forc_aer[g, 13]
        aer.flx_dst_dep_dry4_col[c] = forc_aer[g, 14]
        aer.flx_dst_dep_col[c]      = forc_aer[g, 7] + forc_aer[g, 8] +
                                       forc_aer[g, 9] + forc_aer[g, 10] +
                                       forc_aer[g, 11] + forc_aer[g, 12] +
                                       forc_aer[g, 13] + forc_aer[g, 14]
    end

    # If aerosol effect in snow is turned off, zero deposition fluxes
    if !snicar_use_aerosol
        for c in bounds
            aer.flx_bc_dep_dry_col[c]   = 0.0
            aer.flx_bc_dep_wet_col[c]   = 0.0
            aer.flx_bc_dep_phi_col[c]   = 0.0
            aer.flx_bc_dep_pho_col[c]   = 0.0
            aer.flx_bc_dep_col[c]       = 0.0
            aer.flx_oc_dep_dry_col[c]   = 0.0
            aer.flx_oc_dep_wet_col[c]   = 0.0
            aer.flx_oc_dep_phi_col[c]   = 0.0
            aer.flx_oc_dep_pho_col[c]   = 0.0
            aer.flx_oc_dep_col[c]       = 0.0
            aer.flx_dst_dep_wet1_col[c] = 0.0
            aer.flx_dst_dep_dry1_col[c] = 0.0
            aer.flx_dst_dep_wet2_col[c] = 0.0
            aer.flx_dst_dep_dry2_col[c] = 0.0
            aer.flx_dst_dep_wet3_col[c] = 0.0
            aer.flx_dst_dep_dry3_col[c] = 0.0
            aer.flx_dst_dep_wet4_col[c] = 0.0
            aer.flx_dst_dep_dry4_col[c] = 0.0
            aer.flx_dst_dep_col[c]      = 0.0
        end
    end

    # Add deposition fluxes into top snow layer
    # Done after inter-layer fluxes so aerosol is present in top layer
    # before radiative calculations
    for c in bounds
        mask_snow[c] || continue

        # Fortran: snl(c)+1 → Julia: snl[c] + nlevsno + 1
        top_j = snl[c] + nlevsno + 1

        aer.mss_bcphi_col[c, top_j] += aer.flx_bc_dep_phi_col[c] * dtime
        aer.mss_bcpho_col[c, top_j] += aer.flx_bc_dep_pho_col[c] * dtime
        aer.mss_ocphi_col[c, top_j] += aer.flx_oc_dep_phi_col[c] * dtime
        aer.mss_ocpho_col[c, top_j] += aer.flx_oc_dep_pho_col[c] * dtime

        aer.mss_dst1_col[c, top_j] += (aer.flx_dst_dep_dry1_col[c] + aer.flx_dst_dep_wet1_col[c]) * dtime
        aer.mss_dst2_col[c, top_j] += (aer.flx_dst_dep_dry2_col[c] + aer.flx_dst_dep_wet2_col[c]) * dtime
        aer.mss_dst3_col[c, top_j] += (aer.flx_dst_dep_dry3_col[c] + aer.flx_dst_dep_wet3_col[c]) * dtime
        aer.mss_dst4_col[c, top_j] += (aer.flx_dst_dep_dry4_col[c] + aer.flx_dst_dep_wet4_col[c]) * dtime
    end

    return nothing
end

# --------------------------------------------------------------------------
# Restart (stub)
# --------------------------------------------------------------------------

"""
    aerosol_restart!(aer, bounds; flag="read",
                     h2osoi_ice_col, h2osoi_liq_col)

Read/write aerosol fields from/to restart file.
Currently a stub; computes derived mass concentrations on read.

Ported from `Restart` in `AerosolMod.F90`.
"""
function aerosol_restart!(aer::AerosolData, bounds::UnitRange{Int};
                          flag::String = "read",
                          h2osoi_ice_col::Matrix{<:Real} = Matrix{Float64}(undef, 0, 0),
                          h2osoi_liq_col::Matrix{<:Real} = Matrix{Float64}(undef, 0, 0))
    nlevsno = varpar.nlevsno

    if flag == "read" && size(h2osoi_ice_col, 1) > 0
        for j in 1:nlevsno
            for c in bounds
                total_water = h2osoi_ice_col[c, j] + h2osoi_liq_col[c, j]
                if total_water > 0.0
                    aer.mss_cnc_bcpho_col[c, j] = aer.mss_bcpho_col[c, j] / total_water
                    aer.mss_cnc_bcphi_col[c, j] = aer.mss_bcphi_col[c, j] / total_water
                    aer.mss_cnc_ocpho_col[c, j] = aer.mss_ocpho_col[c, j] / total_water
                    aer.mss_cnc_ocphi_col[c, j] = aer.mss_ocphi_col[c, j] / total_water
                    aer.mss_cnc_dst1_col[c, j]  = aer.mss_dst1_col[c, j] / total_water
                    aer.mss_cnc_dst2_col[c, j]  = aer.mss_dst2_col[c, j] / total_water
                    aer.mss_cnc_dst3_col[c, j]  = aer.mss_dst3_col[c, j] / total_water
                    aer.mss_cnc_dst4_col[c, j]  = aer.mss_dst4_col[c, j] / total_water
                else
                    aer.mss_cnc_bcpho_col[c, j] = 0.0
                    aer.mss_cnc_bcphi_col[c, j] = 0.0
                    aer.mss_cnc_ocpho_col[c, j] = 0.0
                    aer.mss_cnc_ocphi_col[c, j] = 0.0
                    aer.mss_cnc_dst1_col[c, j]  = 0.0
                    aer.mss_cnc_dst2_col[c, j]  = 0.0
                    aer.mss_cnc_dst3_col[c, j]  = 0.0
                    aer.mss_cnc_dst4_col[c, j]  = 0.0
                end
            end
        end
    end

    return nothing
end

"""
    aerosol_init_history!(aer, bounds)

Register aerosol fields for history file output.
Stub — requires history infrastructure.

Ported from `InitHistory` in `AerosolMod.F90`.
"""
function aerosol_init_history!(aer::AerosolData{FT}, bounds::UnitRange{Int}) where {FT}
    return nothing
end
