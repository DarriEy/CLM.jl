# ==========================================================================
# Ported from: src/biogeochem/CNDVType.F90
# Dynamic Global Vegetation Simulator (DGVS) data types and initialization.
#
# Public types:
#   DGVSData        — patch-level DGVS state variables
#   DGVEcophysCon   — ecophysiological constants for DGVS (module-level singleton)
#
# Public functions:
#   dgvs_init!          — allocate all vectors
#   dgvs_init_cold!     — cold-start initialization
#   dgv_ecophyscon_init! — initialize ecophysiological constants from pftcon
# ==========================================================================

# --------------------------------------------------------------------------
# DGVSData — patch-level DGVS state
# --------------------------------------------------------------------------

"""
    DGVSData

Patch-level state variables for the Dynamic Global Vegetation Simulator (DGVS).
Used when `use_cndv = true`.

Ported from `dgvs_type` in `CNDVType.F90`.
"""
Base.@kwdef mutable struct DGVSData{FT<:Real}
    # --- Accumulated growing degree days / temperature ---
    agdd_patch          ::Vector{FT} = Float64[]   # accumulated GDD this year (degree-days)
    agddtw_patch        ::Vector{FT} = Float64[]   # accumulated GDD for twmax check (degree-days)
    agdd20_patch        ::Vector{FT} = Float64[]   # 20-year running mean annual GDD (degree-days)
    tmomin20_patch      ::Vector{FT} = Float64[]   # 20-year running mean of annual min monthly temp (K)

    # --- Vegetation state ---
    present_patch       ::Vector{Bool} = Bool[]     # whether PFT present in gridcell
    pftmayexist_patch   ::Vector{Bool} = Bool[]     # whether PFT may establish in gridcell
    nind_patch          ::Vector{FT} = Float64[]    # number of individuals (#/m2)
    lm_ind_patch        ::Vector{FT} = Float64[]    # individual leaf mass (gC/ind)
    lai_ind_patch       ::Vector{FT} = Float64[]    # individual LAI (m2/m2)
    fpcinc_patch        ::Vector{FT} = Float64[]    # foliar projective cover increment (fraction)
    fpcgrid_patch       ::Vector{FT} = Float64[]    # foliar projective cover on gridcell (fraction)
    fpcgridold_patch    ::Vector{FT} = Float64[]    # previous-year fpcgrid (fraction)
    crownarea_patch     ::Vector{FT} = Float64[]    # crown area per individual (m2)
    greffic_patch       ::Vector{FT} = Float64[]    # growth efficiency (gC/m2 leaf/yr)
    heatstress_patch    ::Vector{FT} = Float64[]    # heat stress mortality factor
end

# --------------------------------------------------------------------------
# dgvs_init! — allocate all vectors
# --------------------------------------------------------------------------

"""
    dgvs_init!(d::DGVSData, np::Int)

Allocate all vectors in DGVSData to size `np` (number of patches).
"""
function dgvs_init!(d::DGVSData{FT}, np::Int) where {FT}
    d.agdd_patch        = zeros(FT, np)
    d.agddtw_patch      = zeros(FT, np)
    d.agdd20_patch       = zeros(FT, np)
    d.tmomin20_patch     = zeros(FT, np)
    d.present_patch      = fill(false, np)
    d.pftmayexist_patch  = fill(false, np)
    d.nind_patch         = zeros(FT, np)
    d.lm_ind_patch       = zeros(FT, np)
    d.lai_ind_patch      = zeros(FT, np)
    d.fpcinc_patch       = zeros(FT, np)
    d.fpcgrid_patch      = zeros(FT, np)
    d.fpcgridold_patch   = zeros(FT, np)
    d.crownarea_patch    = zeros(FT, np)
    d.greffic_patch      = zeros(FT, np)
    d.heatstress_patch   = zeros(FT, np)
    return nothing
end

# --------------------------------------------------------------------------
# dgvs_init_cold! — cold-start initialization
# --------------------------------------------------------------------------

"""
    dgvs_init_cold!(d::DGVSData, np::Int)

Cold-start initialization for DGVS state. Sets default values:
- `present_patch` = false (no PFTs present initially)
- `crownarea_patch` = 0
- `nind_patch` = 0
- `agdd20_patch` = 0
- `tmomin20_patch` = TFRZ - 5 = 268.15 K
- All other fields = 0

Ported from `InitCold` in `CNDVType.F90`.
"""
function dgvs_init_cold!(d::DGVSData{FT}, np::Int) where {FT}
    # First allocate
    dgvs_init!(d, np)

    # Then set cold-start values
    fill!(d.present_patch, false)
    fill!(d.pftmayexist_patch, false)
    fill!(d.crownarea_patch, FT(0))
    fill!(d.nind_patch, FT(0))
    fill!(d.lm_ind_patch, FT(0))
    fill!(d.lai_ind_patch, FT(0))
    fill!(d.fpcinc_patch, FT(0))
    fill!(d.fpcgrid_patch, FT(0))
    fill!(d.fpcgridold_patch, FT(0))
    fill!(d.greffic_patch, FT(0))
    fill!(d.heatstress_patch, FT(0))
    fill!(d.agdd_patch, FT(0))
    fill!(d.agddtw_patch, FT(0))
    fill!(d.agdd20_patch, FT(0))
    fill!(d.tmomin20_patch, FT(TFRZ - 5.0))  # 268.15 K
    return nothing
end

# --------------------------------------------------------------------------
# DGVEcophysCon — ecophysiological constants for DGVS
# --------------------------------------------------------------------------

"""
    DGVEcophysCon

Ecophysiological constants for the Dynamic Global Vegetation Simulator.
Maps pftcon parameters (pftpar20/28/29/30/31) to semantic DGVS names, and
provides allometric constants.

Ported from `dgv_ecophyscon` in `CNDVEcosystemDynMod.F90`.
"""
Base.@kwdef mutable struct DGVEcophysCon
    crownarea_max ::Vector{Float64} = Float64[]  # max crown area (m2)         — from pftpar20
    tcmin         ::Vector{Float64} = Float64[]  # min coldest monthly mean T  — from pftpar28
    tcmax         ::Vector{Float64} = Float64[]  # max coldest monthly mean T  — from pftpar29
    gddmin        ::Vector{Float64} = Float64[]  # min GDD (>= 5 deg C)       — from pftpar30
    twmax         ::Vector{Float64} = Float64[]  # upper limit warmest month T — from pftpar31
    reinickerp    ::Vector{Float64} = Float64[]  # Reinicke's p (allometric)   — constant 1.6
    allom1        ::Vector{Float64} = Float64[]  # allometric param 1          — 100 (250 shrubs)
    allom2        ::Vector{Float64} = Float64[]  # allometric param 2          — 40 (8 shrubs)
    allom3        ::Vector{Float64} = Float64[]  # allometric param 3          — 0.5
end

"""
Module-level singleton holding DGVS ecophysiological constants.
Initialized by `dgv_ecophyscon_init!`.
"""
const dgv_ecophyscon = DGVEcophysCon()

# --------------------------------------------------------------------------
# dgv_ecophyscon_init! — populate from pftcon
# --------------------------------------------------------------------------

"""
    dgv_ecophyscon_init!(eco::DGVEcophysCon, pftcon::PftconType)

Initialize DGVS ecophysiological constants from PFT parameter data.

Reads pftpar20/28/29/30/31 from `pftcon` and sets allometric constants.
Shrub PFTs (identified by `pftcon.is_shrub`) get modified allom1 (250)
and allom2 (8) values.

Ported from `dgv_ecophys_init` in `CNDVEcosystemDynMod.F90`.
"""
function dgv_ecophyscon_init!(eco::DGVEcophysCon, pftcon_in::PftconType)
    n = length(pftcon_in.pftpar20)
    if n == 0
        error("dgv_ecophyscon_init!: pftcon must be allocated before calling this function")
    end

    eco.crownarea_max = copy(pftcon_in.pftpar20)
    eco.tcmin         = copy(pftcon_in.pftpar28)
    eco.tcmax         = copy(pftcon_in.pftpar29)
    eco.gddmin        = copy(pftcon_in.pftpar30)
    eco.twmax         = copy(pftcon_in.pftpar31)

    # Allometric constants: default values, modified for shrubs
    eco.reinickerp = fill(REINICKERP, n)
    eco.allom1     = fill(ALLOM1, n)
    eco.allom2     = fill(ALLOM2, n)
    eco.allom3     = fill(ALLOM3, n)

    # Override for shrubs
    for i in 1:n
        if pftcon_in.is_shrub[i]
            eco.allom1[i] = ALLOM1S
            eco.allom2[i] = ALLOM2S
        end
    end

    return nothing
end
