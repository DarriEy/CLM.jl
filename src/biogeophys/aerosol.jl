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
Base.@kwdef mutable struct AerosolData{FT<:Real,
                           V<:AbstractVector{FT},
                           M<:AbstractMatrix{FT}}
    # --- Mass of aerosol species in snow (col, lyr) [kg] ---
    mss_bcpho_col   ::M = Matrix{Float64}(undef, 0, 0)  # hydrophobic BC
    mss_bcphi_col   ::M = Matrix{Float64}(undef, 0, 0)  # hydrophilic BC
    mss_bctot_col   ::M = Matrix{Float64}(undef, 0, 0)  # total BC (pho+phi)
    mss_bc_col_col  ::V = Float64[]                      # column-integrated BC
    mss_bc_top_col  ::V = Float64[]                      # top-layer BC

    mss_ocpho_col   ::M = Matrix{Float64}(undef, 0, 0)  # hydrophobic OC
    mss_ocphi_col   ::M = Matrix{Float64}(undef, 0, 0)  # hydrophilic OC
    mss_octot_col   ::M = Matrix{Float64}(undef, 0, 0)  # total OC (pho+phi)
    mss_oc_col_col  ::V = Float64[]                      # column-integrated OC
    mss_oc_top_col  ::V = Float64[]                      # top-layer OC

    mss_dst1_col    ::M = Matrix{Float64}(undef, 0, 0)  # dust species 1
    mss_dst2_col    ::M = Matrix{Float64}(undef, 0, 0)  # dust species 2
    mss_dst3_col    ::M = Matrix{Float64}(undef, 0, 0)  # dust species 3
    mss_dst4_col    ::M = Matrix{Float64}(undef, 0, 0)  # dust species 4
    mss_dsttot_col  ::M = Matrix{Float64}(undef, 0, 0)  # total dust
    mss_dst_col_col ::V = Float64[]                      # column-integrated dust
    mss_dst_top_col ::V = Float64[]                      # top-layer dust

    # --- Mass concentrations in snow (col, lyr) [kg/kg] ---
    mss_cnc_bcphi_col ::M = Matrix{Float64}(undef, 0, 0)
    mss_cnc_bcpho_col ::M = Matrix{Float64}(undef, 0, 0)
    mss_cnc_ocphi_col ::M = Matrix{Float64}(undef, 0, 0)
    mss_cnc_ocpho_col ::M = Matrix{Float64}(undef, 0, 0)
    mss_cnc_dst1_col  ::M = Matrix{Float64}(undef, 0, 0)
    mss_cnc_dst2_col  ::M = Matrix{Float64}(undef, 0, 0)
    mss_cnc_dst3_col  ::M = Matrix{Float64}(undef, 0, 0)
    mss_cnc_dst4_col  ::M = Matrix{Float64}(undef, 0, 0)

    # --- Deposition fluxes (col) [kg/m2/s] ---
    flx_dst_dep_dry1_col ::V = Float64[]
    flx_dst_dep_wet1_col ::V = Float64[]
    flx_dst_dep_dry2_col ::V = Float64[]
    flx_dst_dep_wet2_col ::V = Float64[]
    flx_dst_dep_dry3_col ::V = Float64[]
    flx_dst_dep_wet3_col ::V = Float64[]
    flx_dst_dep_dry4_col ::V = Float64[]
    flx_dst_dep_wet4_col ::V = Float64[]
    flx_dst_dep_col      ::V = Float64[]

    flx_bc_dep_dry_col   ::V = Float64[]
    flx_bc_dep_wet_col   ::V = Float64[]
    flx_bc_dep_pho_col   ::V = Float64[]
    flx_bc_dep_phi_col   ::V = Float64[]
    flx_bc_dep_col       ::V = Float64[]

    flx_oc_dep_dry_col   ::V = Float64[]
    flx_oc_dep_wet_col   ::V = Float64[]
    flx_oc_dep_pho_col   ::V = Float64[]
    flx_oc_dep_phi_col   ::V = Float64[]
    flx_oc_dep_col       ::V = Float64[]
end

AerosolData{FT}(; kwargs...) where {FT<:Real} =
    AerosolData{FT, Vector{FT}, Matrix{FT}}(; kwargs...)
Adapt.@adapt_structure AerosolData


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
function aerosol_reset_filter!(aer::AerosolData{FT}, mask::AbstractVector{Bool}, bounds::UnitRange{Int}) where {FT}
    for c in bounds
        mask[c] || continue
        aerosol_reset!(aer, c)
    end
    return nothing
end

# --------------------------------------------------------------------------
# AerosolMasses
# --------------------------------------------------------------------------

# Device-view bundles of AerosolData array fields used by the aerosol_masses!
# kernels. Structs are not isbits kernel args, so the per-column physics receives
# these small adapt-able bundles instead of the full AerosolData container.
struct _AerMassMass{M}
    mss_bcpho::M; mss_bcphi::M; mss_bctot::M
    mss_ocpho::M; mss_ocphi::M; mss_octot::M
    mss_dst1::M; mss_dst2::M; mss_dst3::M; mss_dst4::M; mss_dsttot::M
end
Adapt.@adapt_structure _AerMassMass

struct _AerMassCnc{M}
    bcphi::M; bcpho::M; ocphi::M; ocpho::M
    dst1::M; dst2::M; dst3::M; dst4::M
end
Adapt.@adapt_structure _AerMassCnc

struct _AerMassCol{V}
    bc_col::V; oc_col::V; dst_col::V
    bc_top::V; oc_top::V; dst_top::V
end
Adapt.@adapt_structure _AerMassCol

# Columns WITH snow: integrate masses & concentrations over snow layers.
@kernel function _aerosol_masses_on_kernel!(
    h2osno_top_col, snw_rds_col,
    m::_AerMassMass, cnc::_AerMassCnc, colv::_AerMassCol,
    @Const(mask_on), @Const(snl),
    @Const(h2osoi_ice_col), @Const(h2osoi_liq_col), nlevsno::Int)
    c = @index(Global)
    @inbounds if mask_on[c]
        bc_sum = zero(eltype(colv.bc_col))
        oc_sum = zero(eltype(colv.oc_col))
        dst_sum = zero(eltype(colv.dst_col))
        for j in 1:nlevsno
            snowmass = h2osoi_ice_col[c, j] + h2osoi_liq_col[c, j]
            if j >= snl[c] + nlevsno + 1
                m.mss_bctot[c, j]  = m.mss_bcpho[c, j] + m.mss_bcphi[c, j]
                bc_sum += m.mss_bctot[c, j]
                cnc.bcphi[c, j] = m.mss_bcphi[c, j] / snowmass
                cnc.bcpho[c, j] = m.mss_bcpho[c, j] / snowmass

                m.mss_octot[c, j]  = m.mss_ocpho[c, j] + m.mss_ocphi[c, j]
                oc_sum += m.mss_octot[c, j]
                cnc.ocphi[c, j] = m.mss_ocphi[c, j] / snowmass
                cnc.ocpho[c, j] = m.mss_ocpho[c, j] / snowmass

                m.mss_dsttot[c, j] = m.mss_dst1[c, j] + m.mss_dst2[c, j] +
                                     m.mss_dst3[c, j] + m.mss_dst4[c, j]
                dst_sum += m.mss_dsttot[c, j]
                cnc.dst1[c, j] = m.mss_dst1[c, j] / snowmass
                cnc.dst2[c, j] = m.mss_dst2[c, j] / snowmass
                cnc.dst3[c, j] = m.mss_dst3[c, j] / snowmass
                cnc.dst4[c, j] = m.mss_dst4[c, j] / snowmass
            else
                snw_rds_col[c, j] = 0.0

                m.mss_bcpho[c, j] = 0.0
                m.mss_bcphi[c, j] = 0.0
                m.mss_bctot[c, j] = 0.0
                cnc.bcphi[c, j] = 0.0
                cnc.bcpho[c, j] = 0.0

                m.mss_ocpho[c, j] = 0.0
                m.mss_ocphi[c, j] = 0.0
                m.mss_octot[c, j] = 0.0
                cnc.ocphi[c, j] = 0.0
                cnc.ocpho[c, j] = 0.0

                m.mss_dst1[c, j] = 0.0
                m.mss_dst2[c, j] = 0.0
                m.mss_dst3[c, j] = 0.0
                m.mss_dst4[c, j] = 0.0
                m.mss_dsttot[c, j] = 0.0
                cnc.dst1[c, j] = 0.0
                cnc.dst2[c, j] = 0.0
                cnc.dst3[c, j] = 0.0
                cnc.dst4[c, j] = 0.0
            end
        end
        colv.bc_col[c]  = bc_sum
        colv.oc_col[c]  = oc_sum
        colv.dst_col[c] = dst_sum

        top_j = snl[c] + nlevsno + 1
        h2osno_top_col[c] = h2osoi_ice_col[c, top_j] + h2osoi_liq_col[c, top_j]
        colv.bc_top[c]  = m.mss_bctot[c, top_j]
        colv.oc_top[c]  = m.mss_octot[c, top_j]
        colv.dst_top[c] = m.mss_dsttot[c, top_j]
    end
end

# Columns WITHOUT snow: zero all mass variables.
@kernel function _aerosol_masses_off_kernel!(
    mask_off, m::_AerMassMass, cnc::_AerMassCnc, colv::_AerMassCol,
    nlevsno::Int)
    c = @index(Global)
    @inbounds if mask_off[c]
        colv.bc_top[c]  = 0.0
        colv.bc_col[c]  = 0.0
        colv.oc_top[c]  = 0.0
        colv.oc_col[c]  = 0.0
        colv.dst_top[c] = 0.0
        colv.dst_col[c] = 0.0
        for j in 1:nlevsno
            m.mss_bcpho[c, j] = 0.0
            m.mss_bcphi[c, j] = 0.0
            m.mss_bctot[c, j] = 0.0
            cnc.bcphi[c, j] = 0.0
            cnc.bcpho[c, j] = 0.0

            m.mss_ocpho[c, j] = 0.0
            m.mss_ocphi[c, j] = 0.0
            m.mss_octot[c, j] = 0.0
            cnc.ocphi[c, j] = 0.0
            cnc.ocpho[c, j] = 0.0

            m.mss_dst1[c, j] = 0.0
            m.mss_dst2[c, j] = 0.0
            m.mss_dst3[c, j] = 0.0
            m.mss_dst4[c, j] = 0.0
            m.mss_dsttot[c, j] = 0.0
            cnc.dst1[c, j] = 0.0
            cnc.dst2[c, j] = 0.0
            cnc.dst3[c, j] = 0.0
            cnc.dst4[c, j] = 0.0
        end
    end
end

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
                         mask_on::AbstractVector{Bool},
                         mask_off::AbstractVector{Bool},
                         bounds::UnitRange{Int},
                         snl::AbstractVector{<:Integer},
                         h2osoi_ice_col::AbstractMatrix{<:Real},
                         h2osoi_liq_col::AbstractMatrix{<:Real},
                         h2osno_top_col::AbstractVector{<:Real},
                         snw_rds_col::AbstractMatrix{<:Real})
    nlevsno = varpar.nlevsno

    m = _AerMassMass(aer.mss_bcpho_col, aer.mss_bcphi_col, aer.mss_bctot_col,
                     aer.mss_ocpho_col, aer.mss_ocphi_col, aer.mss_octot_col,
                     aer.mss_dst1_col, aer.mss_dst2_col, aer.mss_dst3_col,
                     aer.mss_dst4_col, aer.mss_dsttot_col)
    cnc = _AerMassCnc(aer.mss_cnc_bcphi_col, aer.mss_cnc_bcpho_col,
                      aer.mss_cnc_ocphi_col, aer.mss_cnc_ocpho_col,
                      aer.mss_cnc_dst1_col, aer.mss_cnc_dst2_col,
                      aer.mss_cnc_dst3_col, aer.mss_cnc_dst4_col)
    colv = _AerMassCol(aer.mss_bc_col_col, aer.mss_oc_col_col, aer.mss_dst_col_col,
                       aer.mss_bc_top_col, aer.mss_oc_top_col, aer.mss_dst_top_col)

    # --- Columns with snow ---
    _launch!(_aerosol_masses_on_kernel!,
             h2osno_top_col, snw_rds_col, m, cnc, colv, mask_on, snl,
             h2osoi_ice_col, h2osoi_liq_col, nlevsno;
             ndrange = length(mask_on))

    # --- Columns without snow: zero all mass variables ---
    _launch!(_aerosol_masses_off_kernel!, mask_off, m, cnc, colv, nlevsno;
             ndrange = length(mask_off))

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
                         mask_snow::AbstractVector{Bool},
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
