# ==========================================================================
# Ported from: src/biogeochem/CNPhenologyMod.F90
# Phenology routines for coupled carbon-nitrogen code (CN).
# Handles evergreen, seasonal deciduous, stress deciduous, and crop phenology.
#
# Public functions:
#   cn_phenology_set_params!        — Set parameter defaults (for unit tests)
#   cn_phenology_set_nml!           — Set namelist settings (for unit tests)
#   cn_phenology_init!              — Initialization
#   cn_phenology!                   — Main driver (two-phase)
#   crop_phase!                     — Get current crop phase
#   days_past_planting              — Days since planting
#   seasonal_decid_onset            — Seasonal deciduous onset test
#   seasonal_critical_daylength     — Critical daylength for offset
#   get_swindow                     — Get next sowing window
#   was_sown_in_this_window         — Check if crop sown in current window
#
# Private functions:
#   cn_phenology_climate!           — Climate averaging
#   cn_evergreen_phenology!         — Evergreen phenology
#   cn_season_decid_phenology!      — Seasonal deciduous phenology
#   cn_stress_decid_phenology!      — Stress deciduous phenology
#   crop_phenology!                 — Crop phenology
#   crop_phenology_init!            — Crop phenology init
#   plant_crop!                     — Initialize crop at planting
#   vernalization!                  — Vernalization for winter cereal
#   cn_onset_growth!                — Transfer → display during onset
#   cn_offset_litterfall!           — Display → litter during offset
#   cn_background_litterfall!       — Background litterfall
#   cn_livewood_turnover!           — Live wood → dead wood turnover
#   cn_crop_harvest_to_product_pools! — Crop harvest to product pools
#   cn_litter_to_column!            — Patch litter → column level
# ==========================================================================

# ---------------------------------------------------------------------------
# Parameters (replaces Fortran params_type)
# ---------------------------------------------------------------------------
Base.@kwdef mutable struct PhenologyParams
    crit_dayl             ::Float64 = 39200.0    # critical daylength for senescence (s)
    crit_dayl_at_high_lat ::Float64 = 54000.0    # critical daylength at high latitudes (s)
    crit_dayl_lat_slope   ::Float64 = 720.0      # slope of critical daylength with latitude (s/deg)
    ndays_off             ::Float64 = 30.0       # number of days to complete leaf offset
    fstor2tran            ::Float64 = 0.5        # fraction of storage to move to transfer
    crit_onset_fdd        ::Float64 = 15.0       # critical freezing days for onset GDD trigger
    crit_onset_swi        ::Float64 = 15.0       # critical soil water index for onset
    soilpsi_on            ::Float64 = -0.6       # soil water potential for onset (MPa)
    crit_offset_fdd       ::Float64 = 15.0       # critical freezing days for offset trigger
    crit_offset_swi       ::Float64 = 15.0       # critical soil water index for offset
    soilpsi_off           ::Float64 = -0.8       # soil water potential for offset (MPa)
    lwtop                 ::Float64 = 0.7        # live wood turnover proportion (annual)
    phenology_soil_depth  ::Float64 = 0.08       # soil depth for phenology triggers (m)
    snow5d_thresh_for_onset::Float64 = 0.2       # 5-day snow depth threshold for onset (m)
end

# Module-level instance (Fortran CNPhenologyMod params_inst). Populated from the
# parameter file by readParams_CNPhenology! (infrastructure/read_params.jl) and
# consumed by cn_driver_no_leaching!. The struct defaults above mirror Fortran's
# CNPhenologySetParams (its *unit-test* defaults) — the production values come
# from the params file, where e.g. ndays_off is 15 days, not 30.
const cn_phenology_params = PhenologyParams()

# ---------------------------------------------------------------------------
# Module-level state (replaces Fortran module variables)
# ---------------------------------------------------------------------------
Base.@kwdef mutable struct PhenologyState
    dt                    ::Float64 = 0.0        # time step (s)
    fracday               ::Float64 = 0.0        # fraction of day per timestep
    crit_dayl             ::Float64 = 0.0        # critical daylength for offset (s)
    ndays_off             ::Float64 = 0.0        # number of days to complete offset
    fstor2tran            ::Float64 = 0.0        # fraction of storage → transfer
    crit_onset_fdd        ::Float64 = 0.0        # critical freezing days for onset
    crit_onset_swi        ::Float64 = 0.0        # critical soil water index for onset
    soilpsi_on            ::Float64 = 0.0        # water potential for onset (MPa)
    crit_offset_fdd       ::Float64 = 0.0        # critical freezing days for offset
    crit_offset_swi       ::Float64 = 0.0        # critical soil water index for offset
    soilpsi_off           ::Float64 = 0.0        # water potential for offset (MPa)
    lwtop                 ::Float64 = 0.0        # live wood turnover rate (per second)
    phenology_soil_layer  ::Int     = 1          # soil layer index for phenology triggers
    # Crop constants
    p1d                   ::Float64 = 0.004      # photoperiod factor for vernalization
    p1v                   ::Float64 = 0.003      # vernalization factor constant
    hti                   ::Float64 = 1.0        # cold hardening index threshold
    tbase                 ::Float64 = 0.0        # base temperature for vernalization
    # Crop state arrays
    inhemi                ::Vector{Int}     = Int[]          # hemisphere per patch (1=NH, 2=SH)
    minplantjday          ::Matrix{Int}     = Matrix{Int}(undef, 0, 0)  # min planting jday(pft, hemi)
    maxplantjday          ::Matrix{Int}     = Matrix{Int}(undef, 0, 0)  # max planting jday(pft, hemi)
    jdayyrstart           ::Vector{Int}     = [1, 182]       # julian day start of year per hemisphere
end

# ---------------------------------------------------------------------------
# PFT constants needed by phenology (from pftcon)
# ---------------------------------------------------------------------------
Base.@kwdef mutable struct PftConPhenology{FT<:Real,
                                           V<:AbstractVector{FT},
                                           M<:AbstractMatrix{FT},
                                           VI<:AbstractVector{<:Integer},
                                           VB<:AbstractVector{Bool}}
    evergreen             ::V = Float64[]   # binary evergreen flag
    season_decid          ::V = Float64[]   # binary seasonal deciduous flag
    season_decid_temperate::V = Float64[]   # binary temperate seasonal deciduous flag
    stress_decid          ::V = Float64[]   # binary stress deciduous flag
    woody                 ::V = Float64[]   # binary woody flag
    leaf_long             ::V = Float64[]   # leaf longevity (yrs)
    leafcn                ::V = Float64[]   # leaf C:N (gC/gN)
    frootcn               ::V = Float64[]   # fine root C:N (gC/gN)
    lflitcn               ::V = Float64[]   # leaf litter C:N (gC/gN)
    livewdcn              ::V = Float64[]   # live wood C:N (gC/gN)
    deadwdcn              ::V = Float64[]   # dead wood C:N (gC/gN)
    ndays_on              ::V = Float64[]   # days to complete onset
    crit_onset_gdd_sf     ::V = Float64[]   # scale factor for crit_onset_gdd
    lf_f                  ::M = Matrix{Float64}(undef, 0, 0)  # leaf litter fractions (pft, litr)
    fr_f                  ::M = Matrix{Float64}(undef, 0, 0)  # fine root litter fractions (pft, litr)
    biofuel_harvfrac      ::V = Float64[]   # biofuel harvest fraction
    repr_structure_harvfrac::M = Matrix{Float64}(undef, 0, 0) # repr structure harvest frac
    # Crop-specific PFT parameters
    minplanttemp          ::V = Float64[]
    planttemp             ::V = Float64[]
    gddmin                ::V = Float64[]
    lfemerg               ::V = Float64[]
    grnfill               ::V = Float64[]
    hybgdd                ::V = Float64[]
    mxmat                 ::VI = Int[]
    manunitro             ::V = Float64[]
    is_pft_known_to_model ::VB = Bool[]
    mnNHplantdate         ::V = Float64[]
    mxNHplantdate         ::V = Float64[]
    mnSHplantdate         ::V = Float64[]
    mxSHplantdate         ::V = Float64[]
end
PftConPhenology{FT}(; kwargs...) where {FT<:Real} =
    PftConPhenology{FT, Vector{FT}, Matrix{FT}, Vector{Int}, Vector{Bool}}(; kwargs...)
Adapt.@adapt_structure PftConPhenology

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
const NOT_Planted   = 999
const NOT_Harvested = 999
const inNH = 1
const inSH = 2

# Critical daylight method constants
const critical_daylight_constant           = 1
const critical_daylight_depends_on_lat     = 2
const critical_daylight_depends_on_veg     = 3
const critical_daylight_depends_on_latnveg = 4

# Critical offset high latitude (degrees)
const critical_offset_high_lat = 65.0

# Harvest reason constants
const HARVEST_REASON_MATURE          = 1.0
const HARVEST_REASON_MAXSEASLENGTH   = 2.0
const HARVEST_REASON_SOWNBADDEC31    = 3.0
const HARVEST_REASON_SOWTODAY        = 4.0
const HARVEST_REASON_SOWTOMORROW     = 5.0
const HARVEST_REASON_IDOPTOMORROW    = 6.0
const HARVEST_REASON_VERNFREEZEKILL  = 7.0

# Namelist-controlled module variables
const _initial_seed_at_planting = Ref(3.0)
const _onset_thresh_depends_on_veg = Ref(false)
const _critical_daylight_method = Ref(critical_daylight_constant)
const _generate_crop_gdds = Ref(false)
const _use_mxmat = Ref(true)
const _min_gddmaturity = Ref(1.0)
const _min_gdd20_baseline = Ref(0.0)

# Minimum critical daylength for onset (s)
const _min_critical_daylength_onset = 39300.0 / 2.0

# ==========================================================================
# KernelAbstractions kernels for simple, fully-independent per-patch loops.
# Each thread handles one patch p, reading/writing only patch-p indices.
# Backend is taken from the written output array (CPU loop or GPU). The patch
# mask is a BitVector (CPU backend); read-only args are @Const.
# ==========================================================================

# --- cn_phenology_climate!: 10-day running-mean t2m accumulator ----------
@kernel function _phen_climate_kernel!(tempavg_t2m, @Const(mask),
                                       @Const(t_ref2m), fracday,
                                       spval)
    T = eltype(tempavg_t2m)
    p = @index(Global)
    @inbounds if mask[p]
        if t_ref2m[p] != spval
            tempavg_t2m[p] = tempavg_t2m[p] + fracday * (t_ref2m[p] - tempavg_t2m[p])
        end
    end
end

phen_climate!(tempavg_t2m, mask, t_ref2m, fracday, spval) =
    _launch!(_phen_climate_kernel!, tempavg_t2m, mask, t_ref2m,
             eltype(tempavg_t2m)(fracday), eltype(tempavg_t2m)(spval))

# --- cn_evergreen_phenology!: storage→transfer for evergreen patches -----
# --- device-view structs for _phen_evergreen_kernel! ---------------------
# 16 read+write per-patch output arrays (all per-patch Vectors, shared eltype V)
Base.@kwdef struct _PhenEvergreenState{V}
    bglfr::V; bgtr::V; lgsf::V
    leafc_stor2xfer::V; frootc_stor2xfer::V
    livestemc_stor2xfer::V; deadstemc_stor2xfer::V; livecrootc_stor2xfer::V
    deadcrootc_stor2xfer::V; gresp_stor2xfer::V
    leafn_stor2xfer::V; frootn_stor2xfer::V
    livestemn_stor2xfer::V; deadstemn_stor2xfer::V; livecrootn_stor2xfer::V
    deadcrootn_stor2xfer::V
end
Adapt.@adapt_structure _PhenEvergreenState

# 13 read-only per-patch storage inputs (@Const conceptually, all per-patch Vectors)
Base.@kwdef struct _PhenEvergreenInputs{V}
    leafc_storage::V; frootc_storage::V
    livestemc_storage::V; deadstemc_storage::V
    livecrootc_storage::V; deadcrootc_storage::V
    gresp_storage::V
    leafn_storage::V; frootn_storage::V
    livestemn_storage::V; deadstemn_storage::V
    livecrootn_storage::V; deadcrootn_storage::V
end
Adapt.@adapt_structure _PhenEvergreenInputs

# isbits scalar bundle at working precision T (routes the 4 Float64 scalars)
Base.@kwdef struct _PhenEvergreenScalars{T}
    dt::T; fstor2tran::T; avg_dayspyr::T; secspday::T
end

# --- cn_evergreen_phenology!: storage→transfer for evergreen patches -----
@kernel function _phen_evergreen_kernel!(
        st::_PhenEvergreenState, inp::_PhenEvergreenInputs,
        @Const(mask), @Const(itype), @Const(evergreen), @Const(woody),
        @Const(leaf_long),
        ps::_PhenEvergreenScalars)
    p = @index(Global)
    @inbounds if mask[p]
        T = eltype(st.bglfr)

        # --- alias scalar bundle fields to Fortran-named locals ---
        dt          = ps.dt
        fstor2tran  = ps.fstor2tran
        avg_dayspyr = ps.avg_dayspyr
        secspday    = ps.secspday

        # --- alias read+write per-patch output arrays ---
        bglfr                = st.bglfr
        bgtr                 = st.bgtr
        lgsf                 = st.lgsf
        leafc_stor2xfer      = st.leafc_stor2xfer
        frootc_stor2xfer     = st.frootc_stor2xfer
        livestemc_stor2xfer  = st.livestemc_stor2xfer
        deadstemc_stor2xfer  = st.deadstemc_stor2xfer
        livecrootc_stor2xfer = st.livecrootc_stor2xfer
        deadcrootc_stor2xfer = st.deadcrootc_stor2xfer
        gresp_stor2xfer      = st.gresp_stor2xfer
        leafn_stor2xfer      = st.leafn_stor2xfer
        frootn_stor2xfer     = st.frootn_stor2xfer
        livestemn_stor2xfer  = st.livestemn_stor2xfer
        deadstemn_stor2xfer  = st.deadstemn_stor2xfer
        livecrootn_stor2xfer = st.livecrootn_stor2xfer
        deadcrootn_stor2xfer = st.deadcrootn_stor2xfer

        # --- alias read-only per-patch storage inputs ---
        leafc_storage      = inp.leafc_storage
        frootc_storage     = inp.frootc_storage
        livestemc_storage  = inp.livestemc_storage
        deadstemc_storage  = inp.deadstemc_storage
        livecrootc_storage = inp.livecrootc_storage
        deadcrootc_storage = inp.deadcrootc_storage
        gresp_storage      = inp.gresp_storage
        leafn_storage      = inp.leafn_storage
        frootn_storage     = inp.frootn_storage
        livestemn_storage  = inp.livestemn_storage
        deadstemn_storage  = inp.deadstemn_storage
        livecrootn_storage = inp.livecrootn_storage
        deadcrootn_storage = inp.deadcrootn_storage

        # --- body verbatim ---
        ivt = itype[p] + 1
        if evergreen[ivt] == one(T)
            bglfr[p] = one(T) / (leaf_long[ivt] * avg_dayspyr * secspday)
            bgtr[p]  = zero(T)
            lgsf[p]  = zero(T)
            leafc_stor2xfer[p]  = fstor2tran * leafc_storage[p] / dt
            frootc_stor2xfer[p] = fstor2tran * frootc_storage[p] / dt
            if woody[ivt] == one(T)
                livestemc_stor2xfer[p]  = fstor2tran * livestemc_storage[p] / dt
                deadstemc_stor2xfer[p]  = fstor2tran * deadstemc_storage[p] / dt
                livecrootc_stor2xfer[p] = fstor2tran * livecrootc_storage[p] / dt
                deadcrootc_stor2xfer[p] = fstor2tran * deadcrootc_storage[p] / dt
                gresp_stor2xfer[p]      = fstor2tran * gresp_storage[p] / dt
            end
            leafn_stor2xfer[p]  = fstor2tran * leafn_storage[p] / dt
            frootn_stor2xfer[p] = fstor2tran * frootn_storage[p] / dt
            if woody[ivt] == one(T)
                livestemn_stor2xfer[p]  = fstor2tran * livestemn_storage[p] / dt
                deadstemn_stor2xfer[p]  = fstor2tran * deadstemn_storage[p] / dt
                livecrootn_stor2xfer[p] = fstor2tran * livecrootn_storage[p] / dt
                deadcrootn_stor2xfer[p] = fstor2tran * deadcrootn_storage[p] / dt
            end
        end
    end
end

function phen_evergreen!(bglfr, bgtr, lgsf,
        leafc_stor2xfer, frootc_stor2xfer,
        livestemc_stor2xfer, deadstemc_stor2xfer, livecrootc_stor2xfer,
        deadcrootc_stor2xfer, gresp_stor2xfer,
        leafn_stor2xfer, frootn_stor2xfer,
        livestemn_stor2xfer, deadstemn_stor2xfer, livecrootn_stor2xfer,
        deadcrootn_stor2xfer,
        mask, itype, evergreen, woody, leaf_long,
        leafc_storage, frootc_storage,
        livestemc_storage, deadstemc_storage, livecrootc_storage,
        deadcrootc_storage, gresp_storage,
        leafn_storage, frootn_storage,
        livestemn_storage, deadstemn_storage, livecrootn_storage,
        deadcrootn_storage,
        dt, fstor2tran, avg_dayspyr, secspday)
    # working precision = element type of the device output arrays
    T = eltype(bglfr)

    st = _PhenEvergreenState(;
        bglfr=bglfr, bgtr=bgtr, lgsf=lgsf,
        leafc_stor2xfer=leafc_stor2xfer, frootc_stor2xfer=frootc_stor2xfer,
        livestemc_stor2xfer=livestemc_stor2xfer, deadstemc_stor2xfer=deadstemc_stor2xfer,
        livecrootc_stor2xfer=livecrootc_stor2xfer, deadcrootc_stor2xfer=deadcrootc_stor2xfer,
        gresp_stor2xfer=gresp_stor2xfer,
        leafn_stor2xfer=leafn_stor2xfer, frootn_stor2xfer=frootn_stor2xfer,
        livestemn_stor2xfer=livestemn_stor2xfer, deadstemn_stor2xfer=deadstemn_stor2xfer,
        livecrootn_stor2xfer=livecrootn_stor2xfer, deadcrootn_stor2xfer=deadcrootn_stor2xfer)

    inp = _PhenEvergreenInputs(;
        leafc_storage=leafc_storage, frootc_storage=frootc_storage,
        livestemc_storage=livestemc_storage, deadstemc_storage=deadstemc_storage,
        livecrootc_storage=livecrootc_storage, deadcrootc_storage=deadcrootc_storage,
        gresp_storage=gresp_storage,
        leafn_storage=leafn_storage, frootn_storage=frootn_storage,
        livestemn_storage=livestemn_storage, deadstemn_storage=deadstemn_storage,
        livecrootn_storage=livecrootn_storage, deadcrootn_storage=deadcrootn_storage)

    ps = _PhenEvergreenScalars(;
        dt=T(dt), fstor2tran=T(fstor2tran),
        avg_dayspyr=T(avg_dayspyr), secspday=T(secspday))

    # Struct-first kernel: manual backend + synchronize (the bundle args carry no backend).
    backend = _kernel_backend(st.bglfr)
    _phen_evergreen_kernel!(backend)(
        st, inp,
        mask, itype, evergreen, woody, leaf_long,
        ps;
        ndrange = length(mask))
    KA.synchronize(backend)

    return nothing
end

# --- crop_phase!: derive current crop phase into an output array ----------
@kernel function _phen_crop_phase_kernel!(crop_phase_out, @Const(mask),
        @Const(croplive), @Const(gddtsoi), @Const(hui),
        @Const(huileaf), @Const(huigrain),
        planted, leafemerge, grainfill)
    T = eltype(crop_phase_out)
    p = @index(Global)
    @inbounds if mask[p]
        if croplive[p]
            crop_phase_out[p] = planted
            if gddtsoi[p] >= huileaf[p] && hui[p] < huigrain[p]
                crop_phase_out[p] = leafemerge
            elseif hui[p] >= huigrain[p]
                crop_phase_out[p] = grainfill
            end
        end
    end
end

phen_crop_phase!(crop_phase_out, mask, croplive, gddtsoi, hui, huileaf,
                 huigrain, planted, leafemerge, grainfill) =
    _launch!(_phen_crop_phase_kernel!, crop_phase_out, mask, croplive, gddtsoi,
             hui, huileaf, huigrain, eltype(crop_phase_out)(planted),
             eltype(crop_phase_out)(leafemerge), eltype(crop_phase_out)(grainfill))

# --- cn_onset_growth!: transfer→display fluxes during onset --------------
# ==========================================================================
# cn_onset_growth! device-view structs (group >31 kernel args for Metal)
# ==========================================================================
Base.@kwdef struct _PhenOnsetGrowthOut{V}   # the 12 write-only transfer→display growth-flux arrays (per-patch Vectors)
    leafc_xfer_to_leafc::V;     frootc_xfer_to_frootc::V
    leafn_xfer_to_leafn::V;     frootn_xfer_to_frootn::V
    livestemc_xfer_to_livestemc::V; deadstemc_xfer_to_deadstemc::V
    livecrootc_xfer_to_livecrootc::V; deadcrootc_xfer_to_deadcrootc::V
    livestemn_xfer_to_livestemn::V; deadstemn_xfer_to_deadstemn::V
    livecrootn_xfer_to_livecrootn::V; deadcrootn_xfer_to_deadcrootn::V
end
Adapt.@adapt_structure _PhenOnsetGrowthOut

Base.@kwdef struct _PhenOnsetGrowthIn{V}   # read-only per-patch inputs (@Const), all per-patch Vectors share eltype V
    onset_flag::V; onset_counter::V; bgtr::V
    leafc_xfer::V;  frootc_xfer::V
    leafn_xfer::V;  frootn_xfer::V
    livestemc_xfer::V; deadstemc_xfer::V
    livecrootc_xfer::V; deadcrootc_xfer::V
    livestemn_xfer::V; deadstemn_xfer::V
    livecrootn_xfer::V; deadcrootn_xfer::V
end
Adapt.@adapt_structure _PhenOnsetGrowthIn

Base.@kwdef struct _PhenOnsetScalars{T}   # isbits scalar bundle at working precision T
    dt::T
end

# --- cn_onset_growth!: transfer→display fluxes during onset --------------
@kernel function _phen_onset_growth_kernel!(
        out::_PhenOnsetGrowthOut, inp::_PhenOnsetGrowthIn,
        @Const(mask), @Const(itype), @Const(woody),
        ps::_PhenOnsetScalars)
    T = eltype(out.leafc_xfer_to_leafc)
    p = @index(Global)
    @inbounds if mask[p]
        # --- alias scalar bundle to Fortran-named local ---
        dt = ps.dt

        # --- alias write-only outputs to Fortran-named locals ---
        leafc_xfer_to_leafc         = out.leafc_xfer_to_leafc
        frootc_xfer_to_frootc       = out.frootc_xfer_to_frootc
        leafn_xfer_to_leafn         = out.leafn_xfer_to_leafn
        frootn_xfer_to_frootn       = out.frootn_xfer_to_frootn
        livestemc_xfer_to_livestemc = out.livestemc_xfer_to_livestemc
        deadstemc_xfer_to_deadstemc = out.deadstemc_xfer_to_deadstemc
        livecrootc_xfer_to_livecrootc = out.livecrootc_xfer_to_livecrootc
        deadcrootc_xfer_to_deadcrootc = out.deadcrootc_xfer_to_deadcrootc
        livestemn_xfer_to_livestemn = out.livestemn_xfer_to_livestemn
        deadstemn_xfer_to_deadstemn = out.deadstemn_xfer_to_deadstemn
        livecrootn_xfer_to_livecrootn = out.livecrootn_xfer_to_livecrootn
        deadcrootn_xfer_to_deadcrootn = out.deadcrootn_xfer_to_deadcrootn

        # --- alias read-only per-patch inputs to Fortran-named locals ---
        onset_flag      = inp.onset_flag
        onset_counter   = inp.onset_counter
        bgtr            = inp.bgtr
        leafc_xfer      = inp.leafc_xfer
        frootc_xfer     = inp.frootc_xfer
        leafn_xfer      = inp.leafn_xfer
        frootn_xfer     = inp.frootn_xfer
        livestemc_xfer  = inp.livestemc_xfer
        deadstemc_xfer  = inp.deadstemc_xfer
        livecrootc_xfer = inp.livecrootc_xfer
        deadcrootc_xfer = inp.deadcrootc_xfer
        livestemn_xfer  = inp.livestemn_xfer
        deadstemn_xfer  = inp.deadstemn_xfer
        livecrootn_xfer = inp.livecrootn_xfer
        deadcrootn_xfer = inp.deadcrootn_xfer

        # --- body verbatim ---
        ivt = itype[p] + 1
        if onset_flag[p] == one(T)
            if abs(onset_counter[p] - dt) <= dt / T(2)
                t1 = one(T) / dt
            else
                t1 = T(2) / onset_counter[p]
            end
            leafc_xfer_to_leafc[p]   = t1 * leafc_xfer[p]
            frootc_xfer_to_frootc[p] = t1 * frootc_xfer[p]
            leafn_xfer_to_leafn[p]   = t1 * leafn_xfer[p]
            frootn_xfer_to_frootn[p] = t1 * frootn_xfer[p]
            if woody[ivt] == one(T)
                livestemc_xfer_to_livestemc[p]   = t1 * livestemc_xfer[p]
                deadstemc_xfer_to_deadstemc[p]   = t1 * deadstemc_xfer[p]
                livecrootc_xfer_to_livecrootc[p] = t1 * livecrootc_xfer[p]
                deadcrootc_xfer_to_deadcrootc[p] = t1 * deadcrootc_xfer[p]
                livestemn_xfer_to_livestemn[p]   = t1 * livestemn_xfer[p]
                deadstemn_xfer_to_deadstemn[p]   = t1 * deadstemn_xfer[p]
                livecrootn_xfer_to_livecrootn[p] = t1 * livecrootn_xfer[p]
                deadcrootn_xfer_to_deadcrootn[p] = t1 * deadcrootn_xfer[p]
            end
        end
        if bgtr[p] > zero(T)
            leafc_xfer_to_leafc[p]   = leafc_xfer[p] / dt
            frootc_xfer_to_frootc[p] = frootc_xfer[p] / dt
            leafn_xfer_to_leafn[p]   = leafn_xfer[p] / dt
            frootn_xfer_to_frootn[p] = frootn_xfer[p] / dt
            if woody[ivt] == one(T)
                livestemc_xfer_to_livestemc[p]   = livestemc_xfer[p] / dt
                deadstemc_xfer_to_deadstemc[p]   = deadstemc_xfer[p] / dt
                livecrootc_xfer_to_livecrootc[p] = livecrootc_xfer[p] / dt
                deadcrootc_xfer_to_deadcrootc[p] = deadcrootc_xfer[p] / dt
                livestemn_xfer_to_livestemn[p]   = livestemn_xfer[p] / dt
                deadstemn_xfer_to_deadstemn[p]   = deadstemn_xfer[p] / dt
                livecrootn_xfer_to_livecrootn[p] = livecrootn_xfer[p] / dt
                deadcrootn_xfer_to_deadcrootn[p] = deadcrootn_xfer[p] / dt
            end
        end
    end
end

function phen_onset_growth!(
        leafc_xfer_to_leafc, frootc_xfer_to_frootc,
        leafn_xfer_to_leafn, frootn_xfer_to_frootn,
        livestemc_xfer_to_livestemc, deadstemc_xfer_to_deadstemc,
        livecrootc_xfer_to_livecrootc, deadcrootc_xfer_to_deadcrootc,
        livestemn_xfer_to_livestemn, deadstemn_xfer_to_deadstemn,
        livecrootn_xfer_to_livecrootn, deadcrootn_xfer_to_deadcrootn,
        mask, itype, woody, onset_flag, onset_counter, bgtr,
        leafc_xfer, frootc_xfer, leafn_xfer, frootn_xfer,
        livestemc_xfer, deadstemc_xfer, livecrootc_xfer, deadcrootc_xfer,
        livestemn_xfer, deadstemn_xfer, livecrootn_xfer, deadcrootn_xfer,
        dt)

    # working precision = element type of the device flux arrays
    T = eltype(leafc_xfer_to_leafc)

    out = _PhenOnsetGrowthOut(;
        leafc_xfer_to_leafc=leafc_xfer_to_leafc, frootc_xfer_to_frootc=frootc_xfer_to_frootc,
        leafn_xfer_to_leafn=leafn_xfer_to_leafn, frootn_xfer_to_frootn=frootn_xfer_to_frootn,
        livestemc_xfer_to_livestemc=livestemc_xfer_to_livestemc, deadstemc_xfer_to_deadstemc=deadstemc_xfer_to_deadstemc,
        livecrootc_xfer_to_livecrootc=livecrootc_xfer_to_livecrootc, deadcrootc_xfer_to_deadcrootc=deadcrootc_xfer_to_deadcrootc,
        livestemn_xfer_to_livestemn=livestemn_xfer_to_livestemn, deadstemn_xfer_to_deadstemn=deadstemn_xfer_to_deadstemn,
        livecrootn_xfer_to_livecrootn=livecrootn_xfer_to_livecrootn, deadcrootn_xfer_to_deadcrootn=deadcrootn_xfer_to_deadcrootn)

    inp = _PhenOnsetGrowthIn(;
        onset_flag=onset_flag, onset_counter=onset_counter, bgtr=bgtr,
        leafc_xfer=leafc_xfer, frootc_xfer=frootc_xfer,
        leafn_xfer=leafn_xfer, frootn_xfer=frootn_xfer,
        livestemc_xfer=livestemc_xfer, deadstemc_xfer=deadstemc_xfer,
        livecrootc_xfer=livecrootc_xfer, deadcrootc_xfer=deadcrootc_xfer,
        livestemn_xfer=livestemn_xfer, deadstemn_xfer=deadstemn_xfer,
        livecrootn_xfer=livecrootn_xfer, deadcrootn_xfer=deadcrootn_xfer)

    ps = _PhenOnsetScalars(; dt=T(dt))

    # Struct-first kernel: manual backend + synchronize (the bundle args carry no backend).
    backend = _kernel_backend(out.leafc_xfer_to_leafc)
    _phen_onset_growth_kernel!(backend)(
        out, inp,
        mask, itype, woody,
        ps;
        ndrange = length(mask))
    KA.synchronize(backend)

    return nothing
end

# --- cn_offset_litterfall!: display→litter fluxes during offset ----------
@kernel function _phen_offset_litterfall_kernel!(
        leafc_to_litter, frootc_to_litter, leafc_to_litter_fun,
        prev_leafc_to_litter, prev_frootc_to_litter,
        leafn_to_litter, leafn_to_retransn, frootn_to_litter,
        @Const(mask), @Const(itype),
        @Const(offset_flag), @Const(offset_counter),
        @Const(leafc), @Const(frootc), @Const(leafn), @Const(frootn),
        @Const(lflitcn), @Const(leafcn), @Const(frootcn),
        dt, CNratio_floating::Bool, use_fun::Bool)
    T = eltype(leafc_to_litter)
    p = @index(Global)
    @inbounds if mask[p]
        if offset_flag[p] == one(T)
            ivt = itype[p] + 1
            if abs(offset_counter[p] - dt) <= dt / T(2)
                t1 = one(T) / dt
                frootc_to_litter[p] = t1 * frootc[p]
                leafc_to_litter[p]  = t1 * leafc[p]
            else
                t1 = dt * T(2) / (offset_counter[p]^2)
                leafc_to_litter[p]  = prev_leafc_to_litter[p] +
                    t1 * (leafc[p] - prev_leafc_to_litter[p] * offset_counter[p])
                frootc_to_litter[p] = prev_frootc_to_litter[p] +
                    t1 * (frootc[p] - prev_frootc_to_litter[p] * offset_counter[p])
            end

            if CNratio_floating
                fr_leafn_to_litter = T(0.5)
                if leafc[p] == zero(T)
                    ntovr_leaf = zero(T)
                else
                    ntovr_leaf = leafc_to_litter[p] * (leafn[p] / leafc[p])
                end
                leafn_to_litter[p]   = fr_leafn_to_litter * ntovr_leaf
                leafn_to_retransn[p] = ntovr_leaf - leafn_to_litter[p]
                if frootc[p] == zero(T)
                    frootn_to_litter[p] = zero(T)
                else
                    frootn_to_litter[p] = frootc_to_litter[p] * (frootn[p] / frootc[p])
                end
            else
                leafn_to_litter[p]   = leafc_to_litter[p] / lflitcn[ivt]
                leafn_to_retransn[p] = (leafc_to_litter[p] / leafcn[ivt]) - leafn_to_litter[p]
                frootn_to_litter[p]  = frootc_to_litter[p] / frootcn[ivt]
            end

            # FUN uses the current step's leaf litterfall for retranslocation
            # accounting (CNPhenologyMod.F90:3749). Set inside the litterfall
            # branch so non-litterfall patches keep their prior value.
            if use_fun
                leafc_to_litter_fun[p] = leafc_to_litter[p]
            end

            prev_leafc_to_litter[p]  = leafc_to_litter[p]
            prev_frootc_to_litter[p] = frootc_to_litter[p]
        end
    end
end

function phen_offset_litterfall!(
        leafc_to_litter, frootc_to_litter, leafc_to_litter_fun,
        prev_leafc_to_litter, prev_frootc_to_litter,
        leafn_to_litter, leafn_to_retransn, frootn_to_litter,
        mask, itype, offset_flag, offset_counter,
        leafc, frootc, leafn, frootn, lflitcn, leafcn, frootcn,
        dt, CNratio_floating, use_fun)
    _launch!(_phen_offset_litterfall_kernel!,
        leafc_to_litter, frootc_to_litter, leafc_to_litter_fun,
        prev_leafc_to_litter, prev_frootc_to_litter,
        leafn_to_litter, leafn_to_retransn, frootn_to_litter,
        mask, itype, offset_flag, offset_counter,
        leafc, frootc, leafn, frootn, lflitcn, leafcn, frootcn,
        eltype(leafc_to_litter)(dt), CNratio_floating, use_fun)
end

# --- cn_background_litterfall!: background litter fluxes ------------------
@kernel function _phen_background_litterfall_kernel!(
        leafc_to_litter, frootc_to_litter, leafc_to_litter_fun,
        leafn_to_litter, leafn_to_retransn, frootn_to_litter,
        @Const(mask), @Const(itype), @Const(bglfr),
        @Const(leafc), @Const(frootc), @Const(leafn), @Const(frootn),
        @Const(lflitcn), @Const(leafcn), @Const(frootcn),
        CNratio_floating::Bool, use_fun::Bool)
    T = eltype(leafc_to_litter)
    p = @index(Global)
    @inbounds if mask[p]
        if bglfr[p] > zero(T)
            ivt = itype[p] + 1
            leafc_to_litter[p]  = bglfr[p] * leafc[p]
            frootc_to_litter[p] = bglfr[p] * frootc[p]
            # FUN uses the current step's leaf litterfall (CNPhenologyMod.F90:3983)
            if use_fun
                leafc_to_litter_fun[p] = leafc_to_litter[p]
            end
            if CNratio_floating
                fr_leafn_to_litter = T(0.5)
                if leafc[p] == zero(T)
                    ntovr_leaf = zero(T)
                else
                    ntovr_leaf = leafc_to_litter[p] * (leafn[p] / leafc[p])
                end
                leafn_to_litter[p]   = fr_leafn_to_litter * ntovr_leaf
                leafn_to_retransn[p] = ntovr_leaf - leafn_to_litter[p]
                if frootc[p] == zero(T)
                    frootn_to_litter[p] = zero(T)
                else
                    frootn_to_litter[p] = frootc_to_litter[p] * (frootn[p] / frootc[p])
                end
            else
                leafn_to_litter[p]   = leafc_to_litter[p] / lflitcn[ivt]
                leafn_to_retransn[p] = (leafc_to_litter[p] / leafcn[ivt]) - leafn_to_litter[p]
                frootn_to_litter[p]  = frootc_to_litter[p] / frootcn[ivt]
            end
        end
    end
end

function phen_background_litterfall!(
        leafc_to_litter, frootc_to_litter, leafc_to_litter_fun,
        leafn_to_litter, leafn_to_retransn, frootn_to_litter,
        mask, itype, bglfr, leafc, frootc, leafn, frootn,
        lflitcn, leafcn, frootcn, CNratio_floating, use_fun)
    _launch!(_phen_background_litterfall_kernel!,
        leafc_to_litter, frootc_to_litter, leafc_to_litter_fun,
        leafn_to_litter, leafn_to_retransn, frootn_to_litter,
        mask, itype, bglfr, leafc, frootc, leafn, frootn,
        lflitcn, leafcn, frootcn, CNratio_floating, use_fun)
end

# --- cn_livewood_turnover!: live wood → dead wood turnover ----------------
@kernel function _phen_livewood_turnover_kernel!(
        livestemc_to_deadstemc, livestemn_to_deadstemn, livestemn_to_retransn,
        livecrootc_to_deadcrootc, livecrootn_to_deadcrootn, livecrootn_to_retransn,
        @Const(mask), @Const(itype), @Const(woody),
        @Const(livewdcn), @Const(deadwdcn),
        @Const(livestemc), @Const(livestemn),
        @Const(livecrootc), @Const(livecrootn),
        lwtop, CNratio_floating::Bool)
    T = eltype(livestemc_to_deadstemc)
    p = @index(Global)
    @inbounds if mask[p]
        ivt = itype[p] + 1
        if woody[ivt] > zero(T)
            ctovr = livestemc[p] * lwtop
            ntovr = ctovr / livewdcn[ivt]
            livestemc_to_deadstemc[p] = ctovr
            livestemn_to_deadstemn[p] = ctovr / deadwdcn[ivt]
            if CNratio_floating
                if livestemc[p] == zero(T)
                    ntovr = zero(T)
                    livestemn_to_deadstemn[p] = zero(T)
                else
                    ntovr = ctovr * (livestemn[p] / livestemc[p])
                    livestemn_to_deadstemn[p] = ctovr / deadwdcn[ivt]
                end
            end
            livestemn_to_retransn[p] = ntovr - livestemn_to_deadstemn[p]

            ctovr = livecrootc[p] * lwtop
            ntovr = ctovr / livewdcn[ivt]
            livecrootc_to_deadcrootc[p] = ctovr
            livecrootn_to_deadcrootn[p] = ctovr / deadwdcn[ivt]
            if CNratio_floating
                if livecrootc[p] == zero(T)
                    ntovr = zero(T)
                    livecrootn_to_deadcrootn[p] = zero(T)
                else
                    ntovr = ctovr * (livecrootn[p] / livecrootc[p])
                    livecrootn_to_deadcrootn[p] = ctovr / deadwdcn[ivt]
                end
            end
            livecrootn_to_retransn[p] = ntovr - livecrootn_to_deadcrootn[p]
        end
    end
end

function phen_livewood_turnover!(
        livestemc_to_deadstemc, livestemn_to_deadstemn, livestemn_to_retransn,
        livecrootc_to_deadcrootc, livecrootn_to_deadcrootn, livecrootn_to_retransn,
        mask, itype, woody, livewdcn, deadwdcn,
        livestemc, livestemn, livecrootc, livecrootn,
        lwtop, CNratio_floating)
    _launch!(_phen_livewood_turnover_kernel!,
        livestemc_to_deadstemc, livestemn_to_deadstemn, livestemn_to_retransn,
        livecrootc_to_deadcrootc, livecrootn_to_deadcrootn, livecrootn_to_retransn,
        mask, itype, woody, livewdcn, deadwdcn,
        livestemc, livestemn, livecrootc, livecrootn,
        eltype(livestemc_to_deadstemc)(lwtop), CNratio_floating)
end

# --- cn_crop_harvest_to_product_pools!: per-patch harvest sums ------------
@kernel function _phen_crop_harvest_kernel!(
        harvestc_to_cropprodc, harvestn_to_cropprodn,
        @Const(mask),
        @Const(leafc_to_biofuelc), @Const(livestemc_to_biofuelc),
        @Const(leafc_to_removedresiduec), @Const(livestemc_to_removedresiduec),
        @Const(leafn_to_biofueln), @Const(livestemn_to_biofueln),
        @Const(leafn_to_removedresiduen), @Const(livestemn_to_removedresiduen))
    p = @index(Global)
    @inbounds if mask[p]
        harvestc_to_cropprodc[p] = leafc_to_biofuelc[p] + livestemc_to_biofuelc[p] +
            leafc_to_removedresiduec[p] + livestemc_to_removedresiduec[p]
        harvestn_to_cropprodn[p] = leafn_to_biofueln[p] + livestemn_to_biofueln[p] +
            leafn_to_removedresiduen[p] + livestemn_to_removedresiduen[p]
    end
end

phen_crop_harvest!(harvestc_to_cropprodc, harvestn_to_cropprodn, mask,
        leafc_to_biofuelc, livestemc_to_biofuelc,
        leafc_to_removedresiduec, livestemc_to_removedresiduec,
        leafn_to_biofueln, livestemn_to_biofueln,
        leafn_to_removedresiduen, livestemn_to_removedresiduen) =
    _launch!(_phen_crop_harvest_kernel!, harvestc_to_cropprodc,
        harvestn_to_cropprodn, mask,
        leafc_to_biofuelc, livestemc_to_biofuelc,
        leafc_to_removedresiduec, livestemc_to_removedresiduec,
        leafn_to_biofueln, livestemn_to_biofueln,
        leafn_to_removedresiduen, livestemn_to_removedresiduen)

# ==========================================================================
# cn_phenology_set_params!  — Set parameter defaults for unit testing
# ==========================================================================
function cn_phenology_set_params!(params::PhenologyParams)
    params.crit_dayl             = 39200.0
    params.crit_dayl_at_high_lat = 54000.0
    params.crit_dayl_lat_slope   = 720.0
    params.ndays_off             = 30.0
    params.fstor2tran            = 0.5
    params.crit_onset_fdd        = 15.0
    params.crit_onset_swi        = 15.0
    params.soilpsi_on            = -0.6
    params.crit_offset_fdd       = 15.0
    params.crit_offset_swi       = 15.0
    params.soilpsi_off           = -0.8
    params.lwtop                 = 0.7
    params.phenology_soil_depth  = 0.08
    params.snow5d_thresh_for_onset = 0.2
    return nothing
end

# ==========================================================================
# cn_phenology_set_nml!  — Set namelist items for unit testing
# ==========================================================================
function cn_phenology_set_nml!(;
        onset_thresh_depends_on_veg::Bool = false,
        critical_daylight_method_in::Int  = critical_daylight_constant)
    _onset_thresh_depends_on_veg[] = onset_thresh_depends_on_veg
    if critical_daylight_method_in < critical_daylight_constant ||
       critical_daylight_method_in > critical_daylight_depends_on_latnveg
        error("ERROR critical_daylight_method out of range")
    end
    _critical_daylight_method[] = critical_daylight_method_in
    return nothing
end

# ==========================================================================
# cn_phenology_init!  — Initialization (after time manager & pftcon ready)
# ==========================================================================
function cn_phenology_init!(pstate::PhenologyState, params::PhenologyParams,
                            dt_in::Real;
                            use_crop::Bool=false,
                            pftcon::PftConPhenology=PftConPhenology(),
                            patch_data::PatchData=PatchData(),
                            gridcell::GridcellData=GridcellData(),
                            npcropmin::Int=17, npcropmax::Int=78, maxveg::Int=78,
                            find_soil_layer_fn::Function = (depth) -> 1)

    pstate.dt      = dt_in
    pstate.fracday = dt_in / SECSPDAY

    # Set constants from parameters
    pstate.crit_dayl     = params.crit_dayl
    pstate.ndays_off     = params.ndays_off
    pstate.fstor2tran    = params.fstor2tran

    pstate.phenology_soil_layer = find_soil_layer_fn(params.phenology_soil_depth)

    pstate.crit_onset_fdd  = params.crit_onset_fdd
    pstate.crit_onset_swi  = params.crit_onset_swi
    pstate.soilpsi_on      = params.soilpsi_on

    pstate.crit_offset_fdd = params.crit_offset_fdd
    pstate.crit_offset_swi = params.crit_offset_swi
    pstate.soilpsi_off     = params.soilpsi_off

    # Live wood turnover: annual fraction → per second
    pstate.lwtop = params.lwtop / 31536000.0

    # Crop-specific initialization
    if use_crop
        crop_phenology_init!(pstate, pftcon, patch_data, gridcell,
                             npcropmin, npcropmax, maxveg)
    end

    # Error checking for daylength methods
    meth = _critical_daylight_method[]
    if meth == critical_daylight_depends_on_lat ||
       meth == critical_daylight_depends_on_veg ||
       meth == critical_daylight_depends_on_latnveg
        if params.crit_dayl_at_high_lat < params.crit_dayl
            error("ERROR crit_dayl_at_high_lat should be higher than crit_dayl")
        end
        if params.crit_dayl_at_high_lat >= SECSPDAY
            error("ERROR crit_dayl_at_high_lat >= seconds in a day")
        end
        if params.crit_dayl >= SECSPDAY
            error("ERROR crit_dayl >= seconds in a day")
        end
    end
    if meth == critical_daylight_depends_on_lat ||
       meth == critical_daylight_depends_on_latnveg
        if params.crit_dayl_lat_slope <= 0.0
            error("ERROR crit_dayl_lat_slope must be > 0")
        end
    end

    return nothing
end

# ==========================================================================
# cn_phenology!  — Main driver (two-phase)
# ==========================================================================
function cn_phenology!(pstate::PhenologyState, params::PhenologyParams,
                       pftcon::PftConPhenology,
                       mask_soilp::AbstractVector{Bool},   # soil patches
                       mask_pcropp::AbstractVector{Bool},   # prognostic crop patches
                       mask_soilc::AbstractVector{Bool},    # soil columns
                       temperature::TemperatureData,
                       water_diag::WaterDiagnosticBulkData,
                       canopy_state::CanopyStateData,
                       soil_state::SoilStateData,
                       cnveg_state::CNVegStateData,
                       cnveg_cs::CNVegCarbonStateData,
                       cnveg_cf::CNVegCarbonFluxData,
                       cnveg_ns::CNVegNitrogenStateData,
                       cnveg_nf::CNVegNitrogenFluxData,
                       crop::CropData,
                       patch_data::PatchData,
                       gridcell::GridcellData,
                       cn_params::CNSharedParamsData,
                       leaf_prof_patch::AbstractMatrix{<:Real},
                       froot_prof_patch::AbstractMatrix{<:Real},
                       phase::Int;
                       varctl::VarCtl=VarCtl(),
                       is_first_step::Bool=false,
                       avg_dayspyr::Real=365.0,
                       prec10_patch::AbstractVector{<:Real}=Float64[])

    if phase == 1
        cn_phenology_climate!(pstate, mask_soilp, temperature, cnveg_state, crop,
                              patch_data, pftcon)

        cn_evergreen_phenology!(pstate, mask_soilp, pftcon, cnveg_state,
                                cnveg_cs, cnveg_ns, cnveg_cf, cnveg_nf,
                                patch_data, avg_dayspyr)

        cn_season_decid_phenology!(pstate, params, mask_soilp, pftcon,
                                   temperature, water_diag, cnveg_state,
                                   cnveg_cs, cnveg_ns, cnveg_cf, cnveg_nf,
                                   patch_data, gridcell; use_cndv=varctl.use_cndv)

        cn_stress_decid_phenology!(pstate, mask_soilp, pftcon,
                                   soil_state, temperature,
                                   cnveg_state, cnveg_cs, cnveg_ns,
                                   cnveg_cf, cnveg_nf,
                                   patch_data, gridcell, cn_params;
                                   avg_dayspyr=avg_dayspyr,
                                   prec10_patch=prec10_patch)

        if any(mask_pcropp) && !is_first_step
            crop_phenology!(pstate, params, mask_pcropp, pftcon,
                            water_diag, temperature, crop, canopy_state,
                            cnveg_state, cnveg_cs, cnveg_ns, cnveg_cf, cnveg_nf,
                            patch_data, gridcell;
                            varctl=varctl, avg_dayspyr=avg_dayspyr)
        end

    elseif phase == 2
        cn_onset_growth!(pstate, mask_soilp, pftcon,
                         cnveg_state, cnveg_cs, cnveg_ns, cnveg_cf, cnveg_nf,
                         patch_data)

        cn_offset_litterfall!(pstate, mask_soilp, pftcon,
                              cnveg_state, cnveg_cs, cnveg_ns, cnveg_cf, cnveg_nf,
                              crop, patch_data; use_fun=cn_params.use_fun,
                              CNratio_floating=varctl.CNratio_floating,
                              for_testing_no_crop_seed_replenishment=varctl.for_testing_no_crop_seed_replenishment)

        cn_background_litterfall!(pstate, mask_soilp, pftcon,
                                  cnveg_state, cnveg_cs, cnveg_ns, cnveg_cf, cnveg_nf,
                                  patch_data; use_fun=cn_params.use_fun,
                                  CNratio_floating=varctl.CNratio_floating)

        cn_livewood_turnover!(pstate, mask_soilp, pftcon,
                              cnveg_cs, cnveg_ns, cnveg_cf, cnveg_nf,
                              patch_data; CNratio_floating=varctl.CNratio_floating)

        cn_crop_harvest_to_product_pools!(mask_soilp, mask_soilc,
                                          cnveg_cf, cnveg_nf, patch_data;
                                          use_crop=varctl.use_crop)

        cn_litter_to_column!(mask_soilp, pftcon,
                             cnveg_state, cnveg_cf, cnveg_nf,
                             patch_data, leaf_prof_patch, froot_prof_patch;
                             use_grainproduct=false)
    else
        error("bad phase: $phase")
    end

    return nothing
end

# ==========================================================================
# cn_phenology_climate!  — Climate averaging for phenology triggers
# ==========================================================================
function cn_phenology_climate!(pstate::PhenologyState,
                               mask_soilp::AbstractVector{Bool},
                               temperature::TemperatureData,
                               cnveg_state::CNVegStateData,
                               crop::CropData,
                               patch_data::PatchData,
                               pftcon::PftConPhenology)
    fracday = pstate.fracday

    # Update tempavg_t2m accumulator (10-day running mean) — per-patch independent.
    phen_climate!(cnveg_state.tempavg_t2m_patch, mask_soilp,
                  temperature.t_ref2m_patch, fracday, SPVAL)

    return nothing
end

# ==========================================================================
# cn_evergreen_phenology!
# ==========================================================================
function cn_evergreen_phenology!(pstate::PhenologyState,
                                 mask_soilp::AbstractVector{Bool},
                                 pftcon::PftConPhenology,
                                 cnveg_state::CNVegStateData,
                                 cnveg_cs::CNVegCarbonStateData,
                                 cnveg_ns::CNVegNitrogenStateData,
                                 cnveg_cf::CNVegCarbonFluxData,
                                 cnveg_nf::CNVegNitrogenFluxData,
                                 patch_data::PatchData,
                                 avg_dayspyr::Real)
    dt         = pstate.dt
    # The evergreen storage→transfer uses Fortran's FIXED tranr = 0.0002
    # (CNPhenologyMod CNEvergreenPhenology), NOT fstor2tran=0.5 — that 0.5 is the
    # DECIDUOUS onset fraction. Using fstor2tran here drained evergreen storage
    # 2500x too fast (halving leafc_storage every step).
    #
    # KNOWN DIVERGENCE (deliberate, not an oversight — do not "fix" casually):
    # Fortran applies this transfer ONLY when `cn_evergreen_phenology_opt == 1`;
    # with opt == 0 it skips the storage→transfer entirely. `varctl` HAS the flag
    # (`CN_evergreen_phenology_opt::Int = 0`, varctl.jl:112) and its default is 0 —
    # so a faithful port would skip. CLM.jl instead applies tranr = 0.0002
    # UNCONDITIONALLY. The rate is small, but "small" is not "zero", and this is a
    # real fidelity gap, not a rounding detail.
    #
    # It is NOT gated here because the whole 16-biome parity scorecard (1102/1104
    # vars within tol) was validated against the current always-on behaviour;
    # switching it off is a result-changing change that needs its own re-validation
    # against the Fortran reference, not a rider on an unrelated PR. Gating it =
    # `tranr = (varctl.CN_evergreen_phenology_opt == 1) ? 0.0002 : 0.0`, then re-run
    # scripts/parity_run_domain.jl across the scorecard.
    tranr = 0.0002

    # Per-patch independent: evergreen background rates + storage→transfer.
    phen_evergreen!(
        cnveg_state.bglfr_patch, cnveg_state.bgtr_patch, cnveg_state.lgsf_patch,
        cnveg_cf.leafc_storage_to_xfer_patch, cnveg_cf.frootc_storage_to_xfer_patch,
        cnveg_cf.livestemc_storage_to_xfer_patch, cnveg_cf.deadstemc_storage_to_xfer_patch,
        cnveg_cf.livecrootc_storage_to_xfer_patch, cnveg_cf.deadcrootc_storage_to_xfer_patch,
        cnveg_cf.gresp_storage_to_xfer_patch,
        cnveg_nf.leafn_storage_to_xfer_patch, cnveg_nf.frootn_storage_to_xfer_patch,
        cnveg_nf.livestemn_storage_to_xfer_patch, cnveg_nf.deadstemn_storage_to_xfer_patch,
        cnveg_nf.livecrootn_storage_to_xfer_patch, cnveg_nf.deadcrootn_storage_to_xfer_patch,
        mask_soilp, patch_data.itype, pftcon.evergreen, pftcon.woody, pftcon.leaf_long,
        cnveg_cs.leafc_storage_patch, cnveg_cs.frootc_storage_patch,
        cnveg_cs.livestemc_storage_patch, cnveg_cs.deadstemc_storage_patch,
        cnveg_cs.livecrootc_storage_patch, cnveg_cs.deadcrootc_storage_patch,
        cnveg_cs.gresp_storage_patch,
        cnveg_ns.leafn_storage_patch, cnveg_ns.frootn_storage_patch,
        cnveg_ns.livestemn_storage_patch, cnveg_ns.deadstemn_storage_patch,
        cnveg_ns.livecrootn_storage_patch, cnveg_ns.deadcrootn_storage_patch,
        dt, tranr, Float64(avg_dayspyr), SECSPDAY)

    return nothing
end

# ==========================================================================
# cn_season_decid_phenology!
# ==========================================================================
Base.@kwdef struct _PhenSeasonSDState{V}   # read+write per-patch state/flux arrays (all per-patch Vectors)
    bglfr_patch::V; bgtr_patch::V; lgsf_patch::V
    offset_flag_patch::V; offset_counter_patch::V
    dormant_flag_patch::V; days_active_patch::V
    onset_flag_patch::V; onset_counter_patch::V
    onset_gdd_patch::V; onset_gddflag_patch::V
    prev_leafc_to_litter_patch::V; prev_frootc_to_litter_patch::V
    leafc_xfer_to_leafc_patch::V; frootc_xfer_to_frootc_patch::V
    leafn_xfer_to_leafn_patch::V; frootn_xfer_to_frootn_patch::V
    livestemc_xfer_to_livestemc_patch::V; deadstemc_xfer_to_deadstemc_patch::V
    livecrootc_xfer_to_livecrootc_patch::V; deadcrootc_xfer_to_deadcrootc_patch::V
    livestemn_xfer_to_livestemn_patch::V; deadstemn_xfer_to_deadstemn_patch::V
    livecrootn_xfer_to_livecrootn_patch::V; deadcrootn_xfer_to_deadcrootn_patch::V
    leafc_xfer_patch::V; leafn_xfer_patch::V; frootc_xfer_patch::V; frootn_xfer_patch::V
    livestemc_xfer_patch::V; livestemn_xfer_patch::V; deadstemc_xfer_patch::V; deadstemn_xfer_patch::V
    livecrootc_xfer_patch::V; livecrootn_xfer_patch::V; deadcrootc_xfer_patch::V; deadcrootn_xfer_patch::V
    leafc_storage_to_xfer_patch::V; frootc_storage_to_xfer_patch::V
    livestemc_storage_to_xfer_patch::V; deadstemc_storage_to_xfer_patch::V
    livecrootc_storage_to_xfer_patch::V; deadcrootc_storage_to_xfer_patch::V
    gresp_storage_to_xfer_patch::V
    leafn_storage_to_xfer_patch::V; frootn_storage_to_xfer_patch::V
    livestemn_storage_to_xfer_patch::V; deadstemn_storage_to_xfer_patch::V
    livecrootn_storage_to_xfer_patch::V; deadcrootn_storage_to_xfer_patch::V
end
Adapt.@adapt_structure _PhenSeasonSDState

Base.@kwdef struct _PhenSeasonSDInputs{V}   # read-only per-patch inputs (@Const), all per-patch Vectors
    season_decid::V; woody::V; season_decid_temperate::V
    crit_onset_gdd_sf::V; ndays_on::V
    annavg_t2m_patch::V; t_a5min_patch::V
    leafc_storage_patch::V; frootc_storage_patch::V
    livestemc_storage_patch::V; deadstemc_storage_patch::V
    livecrootc_storage_patch::V; deadcrootc_storage_patch::V
    gresp_storage_patch::V
    leafn_storage_patch::V; frootn_storage_patch::V
    livestemn_storage_patch::V; deadstemn_storage_patch::V
    livecrootn_storage_patch::V; deadcrootn_storage_patch::V
end
Adapt.@adapt_structure _PhenSeasonSDInputs

Base.@kwdef struct _PhenSeasonSDColGrc{V,M}   # per-column matrix + per-column/per-gridcell vectors
    t_soisno_col::M
    snow_5day_col::V; soila10_col::V
    dayl::V; prev_dayl::V; latdeg::V
end
Adapt.@adapt_structure _PhenSeasonSDColGrc

Base.@kwdef struct _PhenSeasonSDScalars{T}   # isbits scalar bundle at working precision T
    dt::T; fracday::T; crit_dayl::T; ndays_off::T
    fstor2tran::T; snow5d_thresh::T
    crit_dayl_at_high_lat::T; crit_dayl_lat_slope::T
    secspday::T; tfrz::T
    crit_offset_high_lat::T; min_critical_daylength_onset::T
end

@kernel function _phen_season_decid_kernel!(
        st::_PhenSeasonSDState, inp::_PhenSeasonSDInputs, cg::_PhenSeasonSDColGrc,
        @Const(mask_soilp), @Const(itype), @Const(column), @Const(gridcell_idx),
        ps::_PhenSeasonSDScalars,
        soil_layer::Int, ts_layer::Int, use_cndv::Bool,
        crit_dayl_method::Int, onset_thresh_depends_on_veg::Bool)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        T = eltype(st.bglfr_patch)

        # --- alias scalar bundle fields to Fortran-named locals (verbatim body below) ---
        dt                          = ps.dt
        fracday                     = ps.fracday
        crit_dayl                   = ps.crit_dayl
        ndays_off                   = ps.ndays_off
        fstor2tran                  = ps.fstor2tran
        snow5d_thresh               = ps.snow5d_thresh
        crit_dayl_at_high_lat       = ps.crit_dayl_at_high_lat
        crit_dayl_lat_slope         = ps.crit_dayl_lat_slope
        secspday                    = ps.secspday
        tfrz                        = ps.tfrz
        crit_offset_high_lat        = ps.crit_offset_high_lat
        min_critical_daylength_onset = ps.min_critical_daylength_onset

        # --- alias read-only per-patch inputs ---
        season_decid           = inp.season_decid
        woody                  = inp.woody
        season_decid_temperate = inp.season_decid_temperate
        crit_onset_gdd_sf      = inp.crit_onset_gdd_sf
        ndays_on               = inp.ndays_on
        annavg_t2m_patch       = inp.annavg_t2m_patch
        t_a5min_patch          = inp.t_a5min_patch
        leafc_storage_patch      = inp.leafc_storage_patch
        frootc_storage_patch     = inp.frootc_storage_patch
        livestemc_storage_patch  = inp.livestemc_storage_patch
        deadstemc_storage_patch  = inp.deadstemc_storage_patch
        livecrootc_storage_patch = inp.livecrootc_storage_patch
        deadcrootc_storage_patch = inp.deadcrootc_storage_patch
        gresp_storage_patch      = inp.gresp_storage_patch
        leafn_storage_patch      = inp.leafn_storage_patch
        frootn_storage_patch     = inp.frootn_storage_patch
        livestemn_storage_patch  = inp.livestemn_storage_patch
        deadstemn_storage_patch  = inp.deadstemn_storage_patch
        livecrootn_storage_patch = inp.livecrootn_storage_patch
        deadcrootn_storage_patch = inp.deadcrootn_storage_patch

        # --- alias per-column / per-gridcell fields ---
        t_soisno_col = cg.t_soisno_col
        snow_5day_col = cg.snow_5day_col
        soila10_col   = cg.soila10_col
        dayl          = cg.dayl
        prev_dayl     = cg.prev_dayl
        latdeg        = cg.latdeg

        # --- alias read+write per-patch state arrays ---
        bglfr_patch          = st.bglfr_patch
        bgtr_patch           = st.bgtr_patch
        lgsf_patch           = st.lgsf_patch
        offset_flag_patch    = st.offset_flag_patch
        offset_counter_patch = st.offset_counter_patch
        dormant_flag_patch   = st.dormant_flag_patch
        days_active_patch    = st.days_active_patch
        onset_flag_patch     = st.onset_flag_patch
        onset_counter_patch  = st.onset_counter_patch
        onset_gdd_patch      = st.onset_gdd_patch
        onset_gddflag_patch  = st.onset_gddflag_patch
        prev_leafc_to_litter_patch  = st.prev_leafc_to_litter_patch
        prev_frootc_to_litter_patch = st.prev_frootc_to_litter_patch
        leafc_xfer_to_leafc_patch   = st.leafc_xfer_to_leafc_patch
        frootc_xfer_to_frootc_patch = st.frootc_xfer_to_frootc_patch
        leafn_xfer_to_leafn_patch   = st.leafn_xfer_to_leafn_patch
        frootn_xfer_to_frootn_patch = st.frootn_xfer_to_frootn_patch
        livestemc_xfer_to_livestemc_patch   = st.livestemc_xfer_to_livestemc_patch
        deadstemc_xfer_to_deadstemc_patch   = st.deadstemc_xfer_to_deadstemc_patch
        livecrootc_xfer_to_livecrootc_patch = st.livecrootc_xfer_to_livecrootc_patch
        deadcrootc_xfer_to_deadcrootc_patch = st.deadcrootc_xfer_to_deadcrootc_patch
        livestemn_xfer_to_livestemn_patch   = st.livestemn_xfer_to_livestemn_patch
        deadstemn_xfer_to_deadstemn_patch   = st.deadstemn_xfer_to_deadstemn_patch
        livecrootn_xfer_to_livecrootn_patch = st.livecrootn_xfer_to_livecrootn_patch
        deadcrootn_xfer_to_deadcrootn_patch = st.deadcrootn_xfer_to_deadcrootn_patch
        leafc_xfer_patch  = st.leafc_xfer_patch
        leafn_xfer_patch  = st.leafn_xfer_patch
        frootc_xfer_patch = st.frootc_xfer_patch
        frootn_xfer_patch = st.frootn_xfer_patch
        livestemc_xfer_patch  = st.livestemc_xfer_patch
        livestemn_xfer_patch  = st.livestemn_xfer_patch
        deadstemc_xfer_patch  = st.deadstemc_xfer_patch
        deadstemn_xfer_patch  = st.deadstemn_xfer_patch
        livecrootc_xfer_patch = st.livecrootc_xfer_patch
        livecrootn_xfer_patch = st.livecrootn_xfer_patch
        deadcrootc_xfer_patch = st.deadcrootc_xfer_patch
        deadcrootn_xfer_patch = st.deadcrootn_xfer_patch
        leafc_storage_to_xfer_patch  = st.leafc_storage_to_xfer_patch
        frootc_storage_to_xfer_patch = st.frootc_storage_to_xfer_patch
        livestemc_storage_to_xfer_patch  = st.livestemc_storage_to_xfer_patch
        deadstemc_storage_to_xfer_patch  = st.deadstemc_storage_to_xfer_patch
        livecrootc_storage_to_xfer_patch = st.livecrootc_storage_to_xfer_patch
        deadcrootc_storage_to_xfer_patch = st.deadcrootc_storage_to_xfer_patch
        gresp_storage_to_xfer_patch      = st.gresp_storage_to_xfer_patch
        leafn_storage_to_xfer_patch  = st.leafn_storage_to_xfer_patch
        frootn_storage_to_xfer_patch = st.frootn_storage_to_xfer_patch
        livestemn_storage_to_xfer_patch  = st.livestemn_storage_to_xfer_patch
        deadstemn_storage_to_xfer_patch  = st.deadstemn_storage_to_xfer_patch
        livecrootn_storage_to_xfer_patch = st.livecrootn_storage_to_xfer_patch
        deadcrootn_storage_to_xfer_patch = st.deadcrootn_storage_to_xfer_patch

        ivt = itype[p] + 1  # 0-based Fortran → 1-based Julia
        c   = column[p]
        g   = gridcell_idx[p]

        if season_decid[ivt] == one(T)
            # set background rates to 0
            bglfr_patch[p] = zero(T)
            bgtr_patch[p]  = zero(T)
            lgsf_patch[p]  = zero(T)

            # onset GDD sum
            crit_onset_gdd = crit_onset_gdd_sf[ivt] *
                exp(T(4.8) + T(0.13) * (annavg_t2m_patch[p] - tfrz))

            # winter→summer flag
            ws_flag = dayl[g] >= prev_dayl[g] ? one(T) : zero(T)

            # update offset_counter
            if offset_flag_patch[p] == one(T)
                offset_counter_patch[p] -= dt

                if offset_counter_patch[p] < dt / T(2)
                    offset_flag_patch[p]    = zero(T)
                    offset_counter_patch[p]  = zero(T)
                    dormant_flag_patch[p]    = one(T)
                    days_active_patch[p]     = zero(T)

                    prev_leafc_to_litter_patch[p]  = zero(T)
                    prev_frootc_to_litter_patch[p] = zero(T)
                end
            end

            # update onset_counter
            if onset_flag_patch[p] == one(T)
                onset_counter_patch[p] -= dt

                if onset_counter_patch[p] < dt / T(2)
                    onset_flag_patch[p]    = zero(T)
                    onset_counter_patch[p] = zero(T)
                    # zero transfer growth rates
                    leafc_xfer_to_leafc_patch[p]   = zero(T)
                    frootc_xfer_to_frootc_patch[p] = zero(T)
                    leafn_xfer_to_leafn_patch[p]   = zero(T)
                    frootn_xfer_to_frootn_patch[p] = zero(T)
                    if woody[ivt] == one(T)
                        livestemc_xfer_to_livestemc_patch[p]   = zero(T)
                        deadstemc_xfer_to_deadstemc_patch[p]   = zero(T)
                        livecrootc_xfer_to_livecrootc_patch[p] = zero(T)
                        deadcrootc_xfer_to_deadcrootc_patch[p] = zero(T)
                        livestemn_xfer_to_livestemn_patch[p]   = zero(T)
                        deadstemn_xfer_to_deadstemn_patch[p]   = zero(T)
                        livecrootn_xfer_to_livecrootn_patch[p] = zero(T)
                        deadcrootn_xfer_to_deadcrootn_patch[p] = zero(T)
                    end
                    # zero transfer pools
                    leafc_xfer_patch[p]  = zero(T)
                    leafn_xfer_patch[p]  = zero(T)
                    frootc_xfer_patch[p] = zero(T)
                    frootn_xfer_patch[p] = zero(T)
                    if woody[ivt] == one(T)
                        livestemc_xfer_patch[p]  = zero(T)
                        livestemn_xfer_patch[p]  = zero(T)
                        deadstemc_xfer_patch[p]  = zero(T)
                        deadstemn_xfer_patch[p]  = zero(T)
                        livecrootc_xfer_patch[p] = zero(T)
                        livecrootn_xfer_patch[p] = zero(T)
                        deadcrootc_xfer_patch[p] = zero(T)
                        deadcrootn_xfer_patch[p] = zero(T)
                    end
                end
            end

            # test dormant → growth
            if dormant_flag_patch[p] == one(T)
                # t_soisno_col is snow-padded (nlevsno + nlevgrnd); soil layer
                # `soil_layer` lives at `nlevsno + soil_layer` (= ts_layer). Reading at
                # `soil_layer` directly hit a snow slot (frozen year-round) → onset_gdd
                # never accumulated → season_decid leaf-out never fired.
                soilt = t_soisno_col[c, ts_layer]
                snow_5day = snow_5day_col[c]
                soila10 = soila10_col[c]
                t_a5min = t_a5min_patch[p]

                # --- INLINED seasonal_decid_onset ---
                og  = onset_gdd_patch[p]
                ogf = onset_gddflag_patch[p]
                do_result = false

                # switch on GDD sum at winter solstice
                if ogf == zero(T) && ws_flag == one(T)
                    ogf = one(T)
                    og  = zero(T)
                end

                # reset if past summer solstice without reaching threshold
                if ogf == one(T) && ws_flag == zero(T)
                    ogf = zero(T)
                    og  = zero(T)
                end

                # accumulate GDD (smoothed TFRZ branch for AD)
                if ogf == one(T)
                    og += smooth_max(soilt - tfrz, zero(T)) * fracday
                end

                if onset_thresh_depends_on_veg
                    if og > crit_onset_gdd && season_decid_temperate[ivt] == one(T)
                        do_result = true
                    elseif season_decid_temperate[ivt] == zero(T) && ogf == one(T) &&
                           soila10 > tfrz && t_a5min > tfrz && ws_flag == one(T) &&
                           dayl[g] > min_critical_daylength_onset &&
                           snow_5day < snow5d_thresh
                        do_result = true
                    end
                else
                    if og > crit_onset_gdd
                        do_result = true
                    end
                end
                # --- end INLINED seasonal_decid_onset ---

                # update in/out args
                onset_gdd_patch[p]     = og
                onset_gddflag_patch[p]  = ogf

                if do_result
                    onset_flag_patch[p]     = one(T)
                    dormant_flag_patch[p]    = zero(T)
                    onset_gddflag_patch[p]   = zero(T)
                    onset_gdd_patch[p]       = zero(T)
                    onset_counter_patch[p]   = ndays_on[ivt] * secspday

                    # storage → transfer (non-matrix). ONE-TIME-EVENT SMOOTHING: weight the GDD-
                    # threshold onset burst by a smooth step of (onset_gdd − crit_onset_gdd). Under
                    # :auto smooth_heaviside is the EXACT step and do_result ⇒ og>crit ⇒ ow=1, so
                    # this is BYTE-IDENTICAL; under :always ow RAMPS as og crosses crit, turning the
                    # delta-function d(leaf-out)/d(crit_onset_gdd) of the instant event into a finite
                    # gradient (the burst magnitude varies smoothly with how far og is past crit).
                    # The veg-dependent trigger path (daylength/temperature) keeps the hard event (ow=1).
                    ow = onset_thresh_depends_on_veg ? one(T) : smooth_heaviside(og - crit_onset_gdd)
                    fs = ow * fstor2tran
                    leafc_storage_to_xfer_patch[p]  = fs * leafc_storage_patch[p] / dt
                    frootc_storage_to_xfer_patch[p] = fs * frootc_storage_patch[p] / dt
                    if woody[ivt] == one(T)
                        livestemc_storage_to_xfer_patch[p]  = fs * livestemc_storage_patch[p] / dt
                        deadstemc_storage_to_xfer_patch[p]  = fs * deadstemc_storage_patch[p] / dt
                        livecrootc_storage_to_xfer_patch[p] = fs * livecrootc_storage_patch[p] / dt
                        deadcrootc_storage_to_xfer_patch[p] = fs * deadcrootc_storage_patch[p] / dt
                        gresp_storage_to_xfer_patch[p]      = fs * gresp_storage_patch[p] / dt
                    end
                    leafn_storage_to_xfer_patch[p]  = fs * leafn_storage_patch[p] / dt
                    frootn_storage_to_xfer_patch[p] = fs * frootn_storage_patch[p] / dt
                    if woody[ivt] == one(T)
                        livestemn_storage_to_xfer_patch[p]  = fs * livestemn_storage_patch[p] / dt
                        deadstemn_storage_to_xfer_patch[p]  = fs * deadstemn_storage_patch[p] / dt
                        livecrootn_storage_to_xfer_patch[p] = fs * livecrootn_storage_patch[p] / dt
                        deadcrootn_storage_to_xfer_patch[p] = fs * deadcrootn_storage_patch[p] / dt
                    end
                end

            # test growth → offset
            elseif offset_flag_patch[p] == zero(T)
                if use_cndv
                    days_active_patch[p] += fracday
                end

                # --- INLINED seasonal_critical_daylength ---
                crit_daylat = crit_dayl
                if crit_dayl_method == critical_daylight_depends_on_latnveg
                    if season_decid_temperate[ivt] == one(T)
                        crit_daylat = crit_dayl
                    else
                        cd = crit_dayl_at_high_lat -
                             crit_dayl_lat_slope * (crit_offset_high_lat - smooth_abs(latdeg[g]))
                        crit_daylat = smooth_max(cd, crit_dayl)
                    end
                elseif crit_dayl_method == critical_daylight_depends_on_veg
                    if season_decid_temperate[ivt] == one(T)
                        crit_daylat = crit_dayl
                    else
                        crit_daylat = crit_dayl_at_high_lat
                    end
                elseif crit_dayl_method == critical_daylight_depends_on_lat
                    cd = crit_dayl_at_high_lat -
                         crit_dayl_lat_slope * (crit_offset_high_lat - smooth_abs(latdeg[g]))
                    crit_daylat = smooth_max(cd, crit_dayl)
                else # critical_daylight_constant (host-validated in cn_phenology_set_nml!)
                    crit_daylat = crit_dayl
                end
                # --- end INLINED seasonal_critical_daylength ---

                if ws_flag == zero(T) && dayl[g] < crit_daylat
                    offset_flag_patch[p]    = one(T)
                    offset_counter_patch[p]  = ndays_off * secspday
                    prev_leafc_to_litter_patch[p]  = zero(T)
                    prev_frootc_to_litter_patch[p] = zero(T)
                end
            end
        end
    end
end

function cn_season_decid_phenology!(pstate::PhenologyState,
                                    params::PhenologyParams,
                                    mask_soilp::AbstractVector{Bool},
                                    pftcon::PftConPhenology,
                                    temperature::TemperatureData,
                                    water_diag::WaterDiagnosticBulkData,
                                    cnveg_state::CNVegStateData,
                                    cnveg_cs::CNVegCarbonStateData,
                                    cnveg_ns::CNVegNitrogenStateData,
                                    cnveg_cf::CNVegCarbonFluxData,
                                    cnveg_nf::CNVegNitrogenFluxData,
                                    patch_data::PatchData,
                                    gridcell::GridcellData;
                                    use_cndv::Bool=false)
    dt         = pstate.dt
    fracday    = pstate.fracday
    crit_dayl  = pstate.crit_dayl
    ndays_off  = pstate.ndays_off
    fstor2tran = pstate.fstor2tran
    soil_layer = pstate.phenology_soil_layer

    # host-resolved config flags (never read Refs / branch on String inside kernel)
    crit_dayl_method            = _critical_daylight_method[]
    onset_thresh_depends_on_veg = _onset_thresh_depends_on_veg[]

    # working precision = element type of the device state arrays
    T = eltype(cnveg_state.bglfr_patch)

    st = _PhenSeasonSDState(;
        bglfr_patch=cnveg_state.bglfr_patch, bgtr_patch=cnveg_state.bgtr_patch, lgsf_patch=cnveg_state.lgsf_patch,
        offset_flag_patch=cnveg_state.offset_flag_patch, offset_counter_patch=cnveg_state.offset_counter_patch,
        dormant_flag_patch=cnveg_state.dormant_flag_patch, days_active_patch=cnveg_state.days_active_patch,
        onset_flag_patch=cnveg_state.onset_flag_patch, onset_counter_patch=cnveg_state.onset_counter_patch,
        onset_gdd_patch=cnveg_state.onset_gdd_patch, onset_gddflag_patch=cnveg_state.onset_gddflag_patch,
        prev_leafc_to_litter_patch=cnveg_cf.prev_leafc_to_litter_patch,
        prev_frootc_to_litter_patch=cnveg_cf.prev_frootc_to_litter_patch,
        leafc_xfer_to_leafc_patch=cnveg_cf.leafc_xfer_to_leafc_patch,
        frootc_xfer_to_frootc_patch=cnveg_cf.frootc_xfer_to_frootc_patch,
        leafn_xfer_to_leafn_patch=cnveg_nf.leafn_xfer_to_leafn_patch,
        frootn_xfer_to_frootn_patch=cnveg_nf.frootn_xfer_to_frootn_patch,
        livestemc_xfer_to_livestemc_patch=cnveg_cf.livestemc_xfer_to_livestemc_patch,
        deadstemc_xfer_to_deadstemc_patch=cnveg_cf.deadstemc_xfer_to_deadstemc_patch,
        livecrootc_xfer_to_livecrootc_patch=cnveg_cf.livecrootc_xfer_to_livecrootc_patch,
        deadcrootc_xfer_to_deadcrootc_patch=cnveg_cf.deadcrootc_xfer_to_deadcrootc_patch,
        livestemn_xfer_to_livestemn_patch=cnveg_nf.livestemn_xfer_to_livestemn_patch,
        deadstemn_xfer_to_deadstemn_patch=cnveg_nf.deadstemn_xfer_to_deadstemn_patch,
        livecrootn_xfer_to_livecrootn_patch=cnveg_nf.livecrootn_xfer_to_livecrootn_patch,
        deadcrootn_xfer_to_deadcrootn_patch=cnveg_nf.deadcrootn_xfer_to_deadcrootn_patch,
        leafc_xfer_patch=cnveg_cs.leafc_xfer_patch, leafn_xfer_patch=cnveg_ns.leafn_xfer_patch,
        frootc_xfer_patch=cnveg_cs.frootc_xfer_patch, frootn_xfer_patch=cnveg_ns.frootn_xfer_patch,
        livestemc_xfer_patch=cnveg_cs.livestemc_xfer_patch, livestemn_xfer_patch=cnveg_ns.livestemn_xfer_patch,
        deadstemc_xfer_patch=cnveg_cs.deadstemc_xfer_patch, deadstemn_xfer_patch=cnveg_ns.deadstemn_xfer_patch,
        livecrootc_xfer_patch=cnveg_cs.livecrootc_xfer_patch, livecrootn_xfer_patch=cnveg_ns.livecrootn_xfer_patch,
        deadcrootc_xfer_patch=cnveg_cs.deadcrootc_xfer_patch, deadcrootn_xfer_patch=cnveg_ns.deadcrootn_xfer_patch,
        leafc_storage_to_xfer_patch=cnveg_cf.leafc_storage_to_xfer_patch,
        frootc_storage_to_xfer_patch=cnveg_cf.frootc_storage_to_xfer_patch,
        livestemc_storage_to_xfer_patch=cnveg_cf.livestemc_storage_to_xfer_patch,
        deadstemc_storage_to_xfer_patch=cnveg_cf.deadstemc_storage_to_xfer_patch,
        livecrootc_storage_to_xfer_patch=cnveg_cf.livecrootc_storage_to_xfer_patch,
        deadcrootc_storage_to_xfer_patch=cnveg_cf.deadcrootc_storage_to_xfer_patch,
        gresp_storage_to_xfer_patch=cnveg_cf.gresp_storage_to_xfer_patch,
        leafn_storage_to_xfer_patch=cnveg_nf.leafn_storage_to_xfer_patch,
        frootn_storage_to_xfer_patch=cnveg_nf.frootn_storage_to_xfer_patch,
        livestemn_storage_to_xfer_patch=cnveg_nf.livestemn_storage_to_xfer_patch,
        deadstemn_storage_to_xfer_patch=cnveg_nf.deadstemn_storage_to_xfer_patch,
        livecrootn_storage_to_xfer_patch=cnveg_nf.livecrootn_storage_to_xfer_patch,
        deadcrootn_storage_to_xfer_patch=cnveg_nf.deadcrootn_storage_to_xfer_patch)

    inp = _PhenSeasonSDInputs(;
        season_decid=pftcon.season_decid, woody=pftcon.woody,
        season_decid_temperate=pftcon.season_decid_temperate,
        crit_onset_gdd_sf=pftcon.crit_onset_gdd_sf, ndays_on=pftcon.ndays_on,
        annavg_t2m_patch=cnveg_state.annavg_t2m_patch, t_a5min_patch=temperature.t_a5min_patch,
        leafc_storage_patch=cnveg_cs.leafc_storage_patch, frootc_storage_patch=cnveg_cs.frootc_storage_patch,
        livestemc_storage_patch=cnveg_cs.livestemc_storage_patch, deadstemc_storage_patch=cnveg_cs.deadstemc_storage_patch,
        livecrootc_storage_patch=cnveg_cs.livecrootc_storage_patch, deadcrootc_storage_patch=cnveg_cs.deadcrootc_storage_patch,
        gresp_storage_patch=cnveg_cs.gresp_storage_patch,
        leafn_storage_patch=cnveg_ns.leafn_storage_patch, frootn_storage_patch=cnveg_ns.frootn_storage_patch,
        livestemn_storage_patch=cnveg_ns.livestemn_storage_patch, deadstemn_storage_patch=cnveg_ns.deadstemn_storage_patch,
        livecrootn_storage_patch=cnveg_ns.livecrootn_storage_patch, deadcrootn_storage_patch=cnveg_ns.deadcrootn_storage_patch)

    cg = _PhenSeasonSDColGrc(;
        t_soisno_col=temperature.t_soisno_col,
        snow_5day_col=water_diag.snow_5day_col, soila10_col=temperature.soila10_col,
        dayl=gridcell.dayl, prev_dayl=gridcell.prev_dayl, latdeg=gridcell.latdeg)

    ps = _PhenSeasonSDScalars(;
        dt=T(dt), fracday=T(fracday), crit_dayl=T(crit_dayl), ndays_off=T(ndays_off),
        fstor2tran=T(fstor2tran), snow5d_thresh=T(params.snow5d_thresh_for_onset),
        crit_dayl_at_high_lat=T(params.crit_dayl_at_high_lat),
        crit_dayl_lat_slope=T(params.crit_dayl_lat_slope),
        secspday=T(SECSPDAY), tfrz=T(TFRZ),
        crit_offset_high_lat=T(critical_offset_high_lat),
        min_critical_daylength_onset=T(_min_critical_daylength_onset))

    # Struct-first kernel: manual backend + synchronize (the bundle args carry no backend).
    backend = _kernel_backend(st.bglfr_patch)
    _phen_season_decid_kernel!(backend)(
        st, inp, cg,
        mask_soilp, patch_data.itype, patch_data.column, patch_data.gridcell,
        ps,
        soil_layer, varpar.nlevsno + soil_layer, use_cndv,
        crit_dayl_method, onset_thresh_depends_on_veg;
        ndrange = length(mask_soilp))
    KA.synchronize(backend)

    return nothing
end

# ==========================================================================
# seasonal_critical_daylength — Critical daylength for seasonal deciduous offset
# ==========================================================================
function seasonal_critical_daylength(g::Int, p::Int,
                                     crit_dayl_val::Real,
                                     params::PhenologyParams,
                                     pftcon::PftConPhenology,
                                     gridcell::GridcellData,
                                     patch_data::PatchData)
    method = _critical_daylight_method[]
    ivt = patch_data.itype[p] + 1  # 0-based Fortran → 1-based Julia

    if method == critical_daylight_depends_on_latnveg
        if pftcon.season_decid_temperate[ivt] == 1.0
            return crit_dayl_val
        else
            cd = params.crit_dayl_at_high_lat -
                 params.crit_dayl_lat_slope * (critical_offset_high_lat - smooth_abs(gridcell.latdeg[g]))
            return smooth_max(cd, crit_dayl_val)
        end
    elseif method == critical_daylight_depends_on_veg
        if pftcon.season_decid_temperate[ivt] == 1.0
            return crit_dayl_val
        else
            return params.crit_dayl_at_high_lat
        end
    elseif method == critical_daylight_depends_on_lat
        cd = params.crit_dayl_at_high_lat -
             params.crit_dayl_lat_slope * (critical_offset_high_lat - smooth_abs(gridcell.latdeg[g]))
        return smooth_max(cd, crit_dayl_val)
    elseif method == critical_daylight_constant
        return crit_dayl_val
    else
        error("ERROR: critical_daylight_method not implemented")
    end
end

# ==========================================================================
# seasonal_decid_onset — Determine if seasonal deciduous onset should happen
# ==========================================================================
function seasonal_decid_onset(onset_gdd::Real, onset_gddflag::Real,
                              soilt::Real, soila10::Real, t_a5min::Real,
                              dayl::Real, snow_5day::Real,
                              ws_flag::Real, crit_onset_gdd::Real,
                              season_decid_temperate::Real,
                              fracday::Real, snow5d_thresh::Real)
    og      = onset_gdd
    ogf     = onset_gddflag
    result  = false

    # switch on GDD sum at winter solstice
    if ogf == 0.0 && ws_flag == 1.0
        ogf = 1.0
        og  = 0.0
    end

    # reset if past summer solstice without reaching threshold
    if ogf == 1.0 && ws_flag == 0.0
        ogf = 0.0
        og  = 0.0
    end

    # accumulate GDD (smoothed TFRZ branch for AD)
    if ogf == 1.0
        og += smooth_max(soilt - TFRZ, 0.0) * fracday
    end

    if _onset_thresh_depends_on_veg[]
        if og > crit_onset_gdd && season_decid_temperate == 1.0
            result = true
        elseif season_decid_temperate == 0.0 && ogf == 1.0 &&
               soila10 > TFRZ && t_a5min > TFRZ && ws_flag == 1.0 &&
               dayl > _min_critical_daylength_onset &&
               snow_5day < snow5d_thresh
            result = true
        end
    else
        if og > crit_onset_gdd
            result = true
        end
    end

    return (result=result, onset_gdd=og, onset_gddflag=ogf)
end

# ==========================================================================
# cn_stress_decid_phenology!
# ==========================================================================
# ==========================================================================
# Device-view bundles for _phen_stress_decid_kernel! (Metal 31-arg limit).
# All per-patch float arrays share one element type -> single Vector param V.
# t_soisno_col / soilpsi_col are matrices -> param M in the Col bundle.
# ==========================================================================

# Read+written per-patch state — split into two bundles (<=~30 fields each).
Base.@kwdef struct _PhenStressSD_S1{V}
    offset_flag_patch::V; offset_counter_patch::V; dormant_flag_patch::V
    days_active_patch::V; prev_leafc_to_litter_patch::V; prev_frootc_to_litter_patch::V
    onset_flag_patch::V; onset_counter_patch::V
    onset_gddflag_patch::V; onset_fdd_patch::V; onset_swi_patch::V; onset_gdd_patch::V
    offset_swi_patch::V; offset_fdd_patch::V
    lgsf_patch::V; bglfr_patch::V; bgtr_patch::V
    leafc_xfer_to_leafc_patch::V; frootc_xfer_to_frootc_patch::V
    leafn_xfer_to_leafn_patch::V; frootn_xfer_to_frootn_patch::V
    livestemc_xfer_to_livestemc_patch::V; deadstemc_xfer_to_deadstemc_patch::V
    livecrootc_xfer_to_livecrootc_patch::V; deadcrootc_xfer_to_deadcrootc_patch::V
    livestemn_xfer_to_livestemn_patch::V; deadstemn_xfer_to_deadstemn_patch::V
    livecrootn_xfer_to_livecrootn_patch::V; deadcrootn_xfer_to_deadcrootn_patch::V
end
Adapt.@adapt_structure _PhenStressSD_S1

Base.@kwdef struct _PhenStressSD_S2{V}
    leafc_xfer_patch::V; leafn_xfer_patch::V; frootc_xfer_patch::V; frootn_xfer_patch::V
    livestemc_xfer_patch::V; livestemn_xfer_patch::V; deadstemc_xfer_patch::V
    deadstemn_xfer_patch::V; livecrootc_xfer_patch::V; livecrootn_xfer_patch::V
    deadcrootc_xfer_patch::V; deadcrootn_xfer_patch::V
    leafc_storage_to_xfer_patch::V; frootc_storage_to_xfer_patch::V
    livestemc_storage_to_xfer_patch::V; deadstemc_storage_to_xfer_patch::V
    livecrootc_storage_to_xfer_patch::V; deadcrootc_storage_to_xfer_patch::V
    gresp_storage_to_xfer_patch::V
    leafn_storage_to_xfer_patch::V; frootn_storage_to_xfer_patch::V
    livestemn_storage_to_xfer_patch::V; deadstemn_storage_to_xfer_patch::V
    livecrootn_storage_to_xfer_patch::V; deadcrootn_storage_to_xfer_patch::V
end
Adapt.@adapt_structure _PhenStressSD_S2

# Read-only per-patch inputs (the @Const float arrays + per-ivt pftcon vectors).
Base.@kwdef struct _PhenStressSD_In{V}
    stress_decid::V; woody::V; crit_onset_gdd_sf::V
    ndays_on::V; leaf_long::V
    annavg_t2m_patch::V
    leafc_storage_patch::V; frootc_storage_patch::V
    livestemc_storage_patch::V; deadstemc_storage_patch::V
    livecrootc_storage_patch::V; deadcrootc_storage_patch::V
    gresp_storage_patch::V
    leafn_storage_patch::V; frootn_storage_patch::V
    livestemn_storage_patch::V; deadstemn_storage_patch::V
    livecrootn_storage_patch::V; deadcrootn_storage_patch::V
    leafc_patch::V; frootc_patch::V
    # 10-day running mean of total precipitation (mm H2O/s). Feeds the
    # `constrain_stress_deciduous_onset` rain trigger below.
    prec10_patch::V
end
Adapt.@adapt_structure _PhenStressSD_In

# Per-column matrices + per-gridcell vector.
Base.@kwdef struct _PhenStressSD_Col{V,M}
    t_soisno_col::M; soilpsi_col::M
    dayl::V
end
Adapt.@adapt_structure _PhenStressSD_Col

# isbits scalar bundle at working precision T (no Float64 reaches the kernel).
# tfrz / secspday are the module Float64 consts routed through here.
Base.@kwdef struct _PhenStressSDScalars{T}
    dt::T; fracday::T; ndays_off_val::T; fstor2tran::T
    crit_onset_fdd::T; crit_onset_swi::T; soilpsi_on::T
    crit_offset_fdd::T; crit_offset_swi::T; soilpsi_off::T
    secspqtrday::T; avg_dayspyr::T
    tfrz::T; secspday::T
    # Rain trigger for stress-deciduous onset (CNPhenologyMod: rain_threshold =
    # 20 mm accumulated over the prec10 window). `constrain_onset` is
    # CNParamsShare's `constrain_stress_deciduous_onset`; when false the condition
    # is unconditionally true, which is the clm4_5 behaviour.
    rain_threshold::T
    constrain_onset::Bool
end

@kernel function _phen_stress_decid_kernel!(
        s1::_PhenStressSD_S1, s2::_PhenStressSD_S2,
        inp::_PhenStressSD_In, col::_PhenStressSD_Col,
        @Const(mask_soilp), @Const(itype), @Const(column), @Const(gridcell_idx),
        sc::_PhenStressSDScalars, soil_layer::Int, ts_layer::Int)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        # working precision from an output array
        T = eltype(s1.offset_flag_patch)

        # --- alias scalar bundle to Fortran-named locals ---
        dt              = sc.dt
        fracday         = sc.fracday
        ndays_off_val   = sc.ndays_off_val
        fstor2tran      = sc.fstor2tran
        crit_onset_fdd  = sc.crit_onset_fdd
        crit_onset_swi  = sc.crit_onset_swi
        soilpsi_on      = sc.soilpsi_on
        crit_offset_fdd = sc.crit_offset_fdd
        crit_offset_swi = sc.crit_offset_swi
        soilpsi_off     = sc.soilpsi_off
        secspqtrday     = sc.secspqtrday
        avg_dayspyr     = sc.avg_dayspyr
        TFRZ            = sc.tfrz
        SECSPDAY        = sc.secspday

        # --- alias S1 state arrays ---
        offset_flag_patch              = s1.offset_flag_patch
        offset_counter_patch           = s1.offset_counter_patch
        dormant_flag_patch             = s1.dormant_flag_patch
        days_active_patch              = s1.days_active_patch
        prev_leafc_to_litter_patch     = s1.prev_leafc_to_litter_patch
        prev_frootc_to_litter_patch    = s1.prev_frootc_to_litter_patch
        onset_flag_patch               = s1.onset_flag_patch
        onset_counter_patch            = s1.onset_counter_patch
        onset_gddflag_patch            = s1.onset_gddflag_patch
        onset_fdd_patch                = s1.onset_fdd_patch
        onset_swi_patch                = s1.onset_swi_patch
        onset_gdd_patch                = s1.onset_gdd_patch
        offset_swi_patch               = s1.offset_swi_patch
        offset_fdd_patch               = s1.offset_fdd_patch
        lgsf_patch                     = s1.lgsf_patch
        bglfr_patch                    = s1.bglfr_patch
        bgtr_patch                     = s1.bgtr_patch
        leafc_xfer_to_leafc_patch      = s1.leafc_xfer_to_leafc_patch
        frootc_xfer_to_frootc_patch    = s1.frootc_xfer_to_frootc_patch
        leafn_xfer_to_leafn_patch      = s1.leafn_xfer_to_leafn_patch
        frootn_xfer_to_frootn_patch    = s1.frootn_xfer_to_frootn_patch
        livestemc_xfer_to_livestemc_patch   = s1.livestemc_xfer_to_livestemc_patch
        deadstemc_xfer_to_deadstemc_patch   = s1.deadstemc_xfer_to_deadstemc_patch
        livecrootc_xfer_to_livecrootc_patch = s1.livecrootc_xfer_to_livecrootc_patch
        deadcrootc_xfer_to_deadcrootc_patch = s1.deadcrootc_xfer_to_deadcrootc_patch
        livestemn_xfer_to_livestemn_patch   = s1.livestemn_xfer_to_livestemn_patch
        deadstemn_xfer_to_deadstemn_patch   = s1.deadstemn_xfer_to_deadstemn_patch
        livecrootn_xfer_to_livecrootn_patch = s1.livecrootn_xfer_to_livecrootn_patch
        deadcrootn_xfer_to_deadcrootn_patch = s1.deadcrootn_xfer_to_deadcrootn_patch

        # --- alias S2 state arrays ---
        leafc_xfer_patch     = s2.leafc_xfer_patch
        leafn_xfer_patch     = s2.leafn_xfer_patch
        frootc_xfer_patch    = s2.frootc_xfer_patch
        frootn_xfer_patch    = s2.frootn_xfer_patch
        livestemc_xfer_patch = s2.livestemc_xfer_patch
        livestemn_xfer_patch = s2.livestemn_xfer_patch
        deadstemc_xfer_patch = s2.deadstemc_xfer_patch
        deadstemn_xfer_patch = s2.deadstemn_xfer_patch
        livecrootc_xfer_patch = s2.livecrootc_xfer_patch
        livecrootn_xfer_patch = s2.livecrootn_xfer_patch
        deadcrootc_xfer_patch = s2.deadcrootc_xfer_patch
        deadcrootn_xfer_patch = s2.deadcrootn_xfer_patch
        leafc_storage_to_xfer_patch      = s2.leafc_storage_to_xfer_patch
        frootc_storage_to_xfer_patch     = s2.frootc_storage_to_xfer_patch
        livestemc_storage_to_xfer_patch  = s2.livestemc_storage_to_xfer_patch
        deadstemc_storage_to_xfer_patch  = s2.deadstemc_storage_to_xfer_patch
        livecrootc_storage_to_xfer_patch = s2.livecrootc_storage_to_xfer_patch
        deadcrootc_storage_to_xfer_patch = s2.deadcrootc_storage_to_xfer_patch
        gresp_storage_to_xfer_patch      = s2.gresp_storage_to_xfer_patch
        leafn_storage_to_xfer_patch      = s2.leafn_storage_to_xfer_patch
        frootn_storage_to_xfer_patch     = s2.frootn_storage_to_xfer_patch
        livestemn_storage_to_xfer_patch  = s2.livestemn_storage_to_xfer_patch
        deadstemn_storage_to_xfer_patch  = s2.deadstemn_storage_to_xfer_patch
        livecrootn_storage_to_xfer_patch = s2.livecrootn_storage_to_xfer_patch
        deadcrootn_storage_to_xfer_patch = s2.deadcrootn_storage_to_xfer_patch

        # --- alias read-only inputs ---
        stress_decid      = inp.stress_decid
        woody             = inp.woody
        crit_onset_gdd_sf = inp.crit_onset_gdd_sf
        ndays_on          = inp.ndays_on
        leaf_long         = inp.leaf_long
        annavg_t2m_patch  = inp.annavg_t2m_patch
        leafc_storage_patch      = inp.leafc_storage_patch
        frootc_storage_patch     = inp.frootc_storage_patch
        livestemc_storage_patch  = inp.livestemc_storage_patch
        deadstemc_storage_patch  = inp.deadstemc_storage_patch
        livecrootc_storage_patch = inp.livecrootc_storage_patch
        deadcrootc_storage_patch = inp.deadcrootc_storage_patch
        gresp_storage_patch      = inp.gresp_storage_patch
        leafn_storage_patch      = inp.leafn_storage_patch
        frootn_storage_patch     = inp.frootn_storage_patch
        livestemn_storage_patch  = inp.livestemn_storage_patch
        deadstemn_storage_patch  = inp.deadstemn_storage_patch
        livecrootn_storage_patch = inp.livecrootn_storage_patch
        deadcrootn_storage_patch = inp.deadcrootn_storage_patch
        leafc_patch  = inp.leafc_patch
        frootc_patch = inp.frootc_patch

        # --- alias per-column / per-gridcell ---
        t_soisno_col = col.t_soisno_col
        soilpsi_col  = col.soilpsi_col
        dayl         = col.dayl

        ivt = itype[p] + 1  # 0-based Fortran → 1-based Julia
        c   = column[p]
        g   = gridcell_idx[p]

        if stress_decid[ivt] == one(T)
            soilt = t_soisno_col[c, ts_layer]    # snow-padded array → nlevsno offset
            psi   = soilpsi_col[c, soil_layer]   # soil-only array → no offset

            crit_onset_gdd = crit_onset_gdd_sf[ivt] *
                exp(T(4.8) + T(0.13) * (annavg_t2m_patch[p] - TFRZ))

            # offset counter
            if offset_flag_patch[p] == one(T)
                offset_counter_patch[p] -= dt
                if offset_counter_patch[p] < dt / T(2)
                    offset_flag_patch[p]    = zero(T)
                    offset_counter_patch[p]  = zero(T)
                    dormant_flag_patch[p]    = one(T)
                    days_active_patch[p]     = zero(T)
                    prev_leafc_to_litter_patch[p]  = zero(T)
                    prev_frootc_to_litter_patch[p] = zero(T)
                end
            end

            # onset counter
            if onset_flag_patch[p] == one(T)
                onset_counter_patch[p] -= dt
                if onset_counter_patch[p] < dt / T(2)
                    onset_flag_patch[p]    = zero(T)
                    onset_counter_patch[p] = zero(T)
                    leafc_xfer_to_leafc_patch[p]   = zero(T)
                    frootc_xfer_to_frootc_patch[p] = zero(T)
                    leafn_xfer_to_leafn_patch[p]   = zero(T)
                    frootn_xfer_to_frootn_patch[p] = zero(T)
                    if woody[ivt] == one(T)
                        livestemc_xfer_to_livestemc_patch[p]   = zero(T)
                        deadstemc_xfer_to_deadstemc_patch[p]   = zero(T)
                        livecrootc_xfer_to_livecrootc_patch[p] = zero(T)
                        deadcrootc_xfer_to_deadcrootc_patch[p] = zero(T)
                        livestemn_xfer_to_livestemn_patch[p]   = zero(T)
                        deadstemn_xfer_to_deadstemn_patch[p]   = zero(T)
                        livecrootn_xfer_to_livecrootn_patch[p] = zero(T)
                        deadcrootn_xfer_to_deadcrootn_patch[p] = zero(T)
                    end
                    leafc_xfer_patch[p]  = zero(T)
                    leafn_xfer_patch[p]  = zero(T)
                    frootc_xfer_patch[p] = zero(T)
                    frootn_xfer_patch[p] = zero(T)
                    if woody[ivt] == one(T)
                        livestemc_xfer_patch[p]  = zero(T)
                        livestemn_xfer_patch[p]  = zero(T)
                        deadstemc_xfer_patch[p]  = zero(T)
                        deadstemn_xfer_patch[p]  = zero(T)
                        livecrootc_xfer_patch[p] = zero(T)
                        livecrootn_xfer_patch[p] = zero(T)
                        deadcrootc_xfer_patch[p] = zero(T)
                        deadcrootn_xfer_patch[p] = zero(T)
                    end
                end
            end

            # test dormant → growth
            if dormant_flag_patch[p] == one(T)
                # freezing degree days (smoothed TFRZ branch for AD)
                if onset_gddflag_patch[p] == zero(T)
                    onset_fdd_patch[p] += smooth_heaviside(TFRZ - soilt) * fracday
                end
                if onset_fdd_patch[p] > crit_onset_fdd
                    onset_gddflag_patch[p] = one(T)
                    onset_fdd_patch[p]     = zero(T)
                    onset_swi_patch[p]     = zero(T)
                end
                # GDD accumulation (smoothed TFRZ branch for AD)
                if onset_gddflag_patch[p] == one(T)
                    onset_gdd_patch[p] += smooth_max(soilt - TFRZ, zero(T)) * fracday
                end
                # Rain trigger (CNPhenologyMod.F90:1638-1644). Was hardcoded
                # `true` because prec10 did not exist as a field; it does now (a
                # real 10-day running mean, accumulated in atm2lnd_update_acc_vars!).
                # Gated on constrain_stress_deciduous_onset, whose CLM.jl default is
                # false → this stays unconditionally true by default, exactly as
                # before. (NOTE: Fortran's CLM5/CLM6 namelist default for that flag
                # is .true.; only clm4_5 is .false. Flipping the CLM.jl default is a
                # separate, parity-validated change — not made here.)
                additional_onset_condition = true
                if sc.constrain_onset
                    if inp.prec10_patch[p] * (T(3600) * T(10) * T(24)) < sc.rain_threshold
                        additional_onset_condition = false
                    end
                end
                if psi >= soilpsi_on
                    onset_swi_patch[p] += fracday
                end
                if onset_swi_patch[p] > crit_onset_swi && additional_onset_condition
                    onset_flag_patch[p] = one(T)
                    if onset_gddflag_patch[p] == one(T) &&
                       onset_gdd_patch[p] < crit_onset_gdd
                        onset_flag_patch[p] = zero(T)
                    end
                end
                # only allow onset if dayl > 6hrs
                if onset_flag_patch[p] == one(T) && dayl[g] <= secspqtrday
                    onset_flag_patch[p] = zero(T)
                end
                # if onset triggered, reset and do storage→transfer
                if onset_flag_patch[p] == one(T)
                    dormant_flag_patch[p]    = zero(T)
                    days_active_patch[p]     = zero(T)
                    onset_gddflag_patch[p]   = zero(T)
                    onset_fdd_patch[p]       = zero(T)
                    onset_gdd_patch[p]       = zero(T)
                    onset_swi_patch[p]       = zero(T)
                    onset_counter_patch[p]   = ndays_on[ivt] * SECSPDAY

                    leafc_storage_to_xfer_patch[p]  = fstor2tran * leafc_storage_patch[p] / dt
                    frootc_storage_to_xfer_patch[p] = fstor2tran * frootc_storage_patch[p] / dt
                    if woody[ivt] == one(T)
                        livestemc_storage_to_xfer_patch[p]  = fstor2tran * livestemc_storage_patch[p] / dt
                        deadstemc_storage_to_xfer_patch[p]  = fstor2tran * deadstemc_storage_patch[p] / dt
                        livecrootc_storage_to_xfer_patch[p] = fstor2tran * livecrootc_storage_patch[p] / dt
                        deadcrootc_storage_to_xfer_patch[p] = fstor2tran * deadcrootc_storage_patch[p] / dt
                        gresp_storage_to_xfer_patch[p]      = fstor2tran * gresp_storage_patch[p] / dt
                    end
                    leafn_storage_to_xfer_patch[p]  = fstor2tran * leafn_storage_patch[p] / dt
                    frootn_storage_to_xfer_patch[p] = fstor2tran * frootn_storage_patch[p] / dt
                    if woody[ivt] == one(T)
                        livestemn_storage_to_xfer_patch[p]  = fstor2tran * livestemn_storage_patch[p] / dt
                        deadstemn_storage_to_xfer_patch[p]  = fstor2tran * deadstemn_storage_patch[p] / dt
                        livecrootn_storage_to_xfer_patch[p] = fstor2tran * livecrootn_storage_patch[p] / dt
                        deadcrootn_storage_to_xfer_patch[p] = fstor2tran * deadcrootn_storage_patch[p] / dt
                    end
                end

            # test growth → offset
            elseif offset_flag_patch[p] == zero(T)
                if psi <= soilpsi_off
                    offset_swi_patch[p] += fracday
                    if offset_swi_patch[p] >= crit_offset_swi &&
                       onset_flag_patch[p] == zero(T)
                        offset_flag_patch[p] = one(T)
                    end
                elseif psi >= soilpsi_on
                    offset_swi_patch[p] -= fracday
                    offset_swi_patch[p] = smooth_max(offset_swi_patch[p], zero(T))
                end
                # freezing day accumulation (smoothed TFRZ branch for AD)
                # Blend thaw-decay and freeze-accumulation using smooth_heaviside
                _hw_above = smooth_heaviside(soilt - TFRZ)   # ~1 when above freezing
                _hw_below = one(T) - _hw_above               # ~1 when below freezing
                if offset_fdd_patch[p] > zero(T)
                    offset_fdd_patch[p] -= _hw_above * fracday
                    offset_fdd_patch[p] = smooth_max(zero(T), offset_fdd_patch[p])
                end
                offset_fdd_patch[p] += _hw_below * fracday
                if offset_fdd_patch[p] > crit_offset_fdd &&
                   onset_flag_patch[p] == zero(T)
                    offset_flag_patch[p] = one(T)
                end
                # force offset if dayl < 6 hrs
                if dayl[g] <= secspqtrday
                    offset_flag_patch[p] = one(T)
                end
                # set offset params
                if offset_flag_patch[p] == one(T)
                    offset_fdd_patch[p]     = zero(T)
                    offset_swi_patch[p]     = zero(T)
                    offset_counter_patch[p]  = ndays_off_val * SECSPDAY
                    prev_leafc_to_litter_patch[p]  = zero(T)
                    prev_frootc_to_litter_patch[p] = zero(T)
                end
            end

            # days active tracking
            if dormant_flag_patch[p] == zero(T)
                days_active_patch[p] += fracday
            end

            # long growing season factor
            lgsf_patch[p] = smooth_clamp(
                T(3) * (days_active_patch[p] - leaf_long[ivt] * avg_dayspyr) / avg_dayspyr,
                zero(T), one(T))

            # background litterfall rate
            if offset_flag_patch[p] == one(T)
                bglfr_patch[p] = zero(T)
            else
                bglfr_patch[p] = (one(T) / (leaf_long[ivt] * avg_dayspyr * SECSPDAY)) *
                    lgsf_patch[p]
            end

            # background transfer rate
            if onset_flag_patch[p] == one(T)
                bgtr_patch[p] = zero(T)
            else
                bgtr_val = (one(T) / (avg_dayspyr * SECSPDAY)) * lgsf_patch[p]
                bgtr_patch[p] = bgtr_val

                # non-matrix storage → transfer
                leafc_storage_to_xfer_patch[p]  = smooth_max(zero(T), leafc_storage_patch[p] - leafc_patch[p]) * bgtr_val
                frootc_storage_to_xfer_patch[p] = smooth_max(zero(T), frootc_storage_patch[p] - frootc_patch[p]) * bgtr_val
                if woody[ivt] == one(T)
                    livestemc_storage_to_xfer_patch[p]  = livestemc_storage_patch[p] * bgtr_val
                    deadstemc_storage_to_xfer_patch[p]  = deadstemc_storage_patch[p] * bgtr_val
                    livecrootc_storage_to_xfer_patch[p] = livecrootc_storage_patch[p] * bgtr_val
                    deadcrootc_storage_to_xfer_patch[p] = deadcrootc_storage_patch[p] * bgtr_val
                    gresp_storage_to_xfer_patch[p]      = gresp_storage_patch[p] * bgtr_val
                end
                leafn_storage_to_xfer_patch[p]  = leafn_storage_patch[p] * bgtr_val
                frootn_storage_to_xfer_patch[p] = frootn_storage_patch[p] * bgtr_val
                if woody[ivt] == one(T)
                    livestemn_storage_to_xfer_patch[p]  = livestemn_storage_patch[p] * bgtr_val
                    deadstemn_storage_to_xfer_patch[p]  = deadstemn_storage_patch[p] * bgtr_val
                    livecrootn_storage_to_xfer_patch[p] = livecrootn_storage_patch[p] * bgtr_val
                    deadcrootn_storage_to_xfer_patch[p] = deadcrootn_storage_patch[p] * bgtr_val
                end
            end
        end
    end
end

function cn_stress_decid_phenology!(pstate::PhenologyState,
                                    mask_soilp::AbstractVector{Bool},
                                    pftcon::PftConPhenology,
                                    soil_state::SoilStateData,
                                    temperature::TemperatureData,
                                    cnveg_state::CNVegStateData,
                                    cnveg_cs::CNVegCarbonStateData,
                                    cnveg_ns::CNVegNitrogenStateData,
                                    cnveg_cf::CNVegCarbonFluxData,
                                    cnveg_nf::CNVegNitrogenFluxData,
                                    patch_data::PatchData,
                                    gridcell::GridcellData,
                                    cn_params::CNSharedParamsData;
                                    avg_dayspyr::Real=365.0,
                                    prec10_patch::AbstractVector{<:Real}=Float64[])
    dt              = pstate.dt
    fracday         = pstate.fracday
    ndays_off_val   = pstate.ndays_off
    fstor2tran      = pstate.fstor2tran
    soil_layer      = pstate.phenology_soil_layer
    crit_onset_fdd  = pstate.crit_onset_fdd
    crit_onset_swi  = pstate.crit_onset_swi
    soilpsi_on      = pstate.soilpsi_on
    crit_offset_fdd = pstate.crit_offset_fdd
    crit_offset_swi = pstate.crit_offset_swi
    soilpsi_off     = pstate.soilpsi_off
    secspqtrday     = SECSPDAY / 4.0

    # Working precision from a representative output array.
    T = eltype(cnveg_state.offset_flag_patch)

    s1 = _PhenStressSD_S1(;
        offset_flag_patch              = cnveg_state.offset_flag_patch,
        offset_counter_patch           = cnveg_state.offset_counter_patch,
        dormant_flag_patch             = cnveg_state.dormant_flag_patch,
        days_active_patch              = cnveg_state.days_active_patch,
        prev_leafc_to_litter_patch     = cnveg_cf.prev_leafc_to_litter_patch,
        prev_frootc_to_litter_patch    = cnveg_cf.prev_frootc_to_litter_patch,
        onset_flag_patch               = cnveg_state.onset_flag_patch,
        onset_counter_patch            = cnveg_state.onset_counter_patch,
        onset_gddflag_patch            = cnveg_state.onset_gddflag_patch,
        onset_fdd_patch                = cnveg_state.onset_fdd_patch,
        onset_swi_patch                = cnveg_state.onset_swi_patch,
        onset_gdd_patch                = cnveg_state.onset_gdd_patch,
        offset_swi_patch               = cnveg_state.offset_swi_patch,
        offset_fdd_patch               = cnveg_state.offset_fdd_patch,
        lgsf_patch                     = cnveg_state.lgsf_patch,
        bglfr_patch                    = cnveg_state.bglfr_patch,
        bgtr_patch                     = cnveg_state.bgtr_patch,
        leafc_xfer_to_leafc_patch      = cnveg_cf.leafc_xfer_to_leafc_patch,
        frootc_xfer_to_frootc_patch    = cnveg_cf.frootc_xfer_to_frootc_patch,
        leafn_xfer_to_leafn_patch      = cnveg_nf.leafn_xfer_to_leafn_patch,
        frootn_xfer_to_frootn_patch    = cnveg_nf.frootn_xfer_to_frootn_patch,
        livestemc_xfer_to_livestemc_patch   = cnveg_cf.livestemc_xfer_to_livestemc_patch,
        deadstemc_xfer_to_deadstemc_patch   = cnveg_cf.deadstemc_xfer_to_deadstemc_patch,
        livecrootc_xfer_to_livecrootc_patch = cnveg_cf.livecrootc_xfer_to_livecrootc_patch,
        deadcrootc_xfer_to_deadcrootc_patch = cnveg_cf.deadcrootc_xfer_to_deadcrootc_patch,
        livestemn_xfer_to_livestemn_patch   = cnveg_nf.livestemn_xfer_to_livestemn_patch,
        deadstemn_xfer_to_deadstemn_patch   = cnveg_nf.deadstemn_xfer_to_deadstemn_patch,
        livecrootn_xfer_to_livecrootn_patch = cnveg_nf.livecrootn_xfer_to_livecrootn_patch,
        deadcrootn_xfer_to_deadcrootn_patch = cnveg_nf.deadcrootn_xfer_to_deadcrootn_patch)

    s2 = _PhenStressSD_S2(;
        leafc_xfer_patch     = cnveg_cs.leafc_xfer_patch,
        leafn_xfer_patch     = cnveg_ns.leafn_xfer_patch,
        frootc_xfer_patch    = cnveg_cs.frootc_xfer_patch,
        frootn_xfer_patch    = cnveg_ns.frootn_xfer_patch,
        livestemc_xfer_patch = cnveg_cs.livestemc_xfer_patch,
        livestemn_xfer_patch = cnveg_ns.livestemn_xfer_patch,
        deadstemc_xfer_patch = cnveg_cs.deadstemc_xfer_patch,
        deadstemn_xfer_patch = cnveg_ns.deadstemn_xfer_patch,
        livecrootc_xfer_patch = cnveg_cs.livecrootc_xfer_patch,
        livecrootn_xfer_patch = cnveg_ns.livecrootn_xfer_patch,
        deadcrootc_xfer_patch = cnveg_cs.deadcrootc_xfer_patch,
        deadcrootn_xfer_patch = cnveg_ns.deadcrootn_xfer_patch,
        leafc_storage_to_xfer_patch      = cnveg_cf.leafc_storage_to_xfer_patch,
        frootc_storage_to_xfer_patch     = cnveg_cf.frootc_storage_to_xfer_patch,
        livestemc_storage_to_xfer_patch  = cnveg_cf.livestemc_storage_to_xfer_patch,
        deadstemc_storage_to_xfer_patch  = cnveg_cf.deadstemc_storage_to_xfer_patch,
        livecrootc_storage_to_xfer_patch = cnveg_cf.livecrootc_storage_to_xfer_patch,
        deadcrootc_storage_to_xfer_patch = cnveg_cf.deadcrootc_storage_to_xfer_patch,
        gresp_storage_to_xfer_patch      = cnveg_cf.gresp_storage_to_xfer_patch,
        leafn_storage_to_xfer_patch      = cnveg_nf.leafn_storage_to_xfer_patch,
        frootn_storage_to_xfer_patch     = cnveg_nf.frootn_storage_to_xfer_patch,
        livestemn_storage_to_xfer_patch  = cnveg_nf.livestemn_storage_to_xfer_patch,
        deadstemn_storage_to_xfer_patch  = cnveg_nf.deadstemn_storage_to_xfer_patch,
        livecrootn_storage_to_xfer_patch = cnveg_nf.livecrootn_storage_to_xfer_patch,
        deadcrootn_storage_to_xfer_patch = cnveg_nf.deadcrootn_storage_to_xfer_patch)

    inp = _PhenStressSD_In(;
        stress_decid      = pftcon.stress_decid,
        woody             = pftcon.woody,
        crit_onset_gdd_sf = pftcon.crit_onset_gdd_sf,
        ndays_on          = pftcon.ndays_on,
        leaf_long         = pftcon.leaf_long,
        annavg_t2m_patch  = cnveg_state.annavg_t2m_patch,
        leafc_storage_patch      = cnveg_cs.leafc_storage_patch,
        frootc_storage_patch     = cnveg_cs.frootc_storage_patch,
        livestemc_storage_patch  = cnveg_cs.livestemc_storage_patch,
        deadstemc_storage_patch  = cnveg_cs.deadstemc_storage_patch,
        livecrootc_storage_patch = cnveg_cs.livecrootc_storage_patch,
        deadcrootc_storage_patch = cnveg_cs.deadcrootc_storage_patch,
        gresp_storage_patch      = cnveg_cs.gresp_storage_patch,
        leafn_storage_patch      = cnveg_ns.leafn_storage_patch,
        frootn_storage_patch     = cnveg_ns.frootn_storage_patch,
        livestemn_storage_patch  = cnveg_ns.livestemn_storage_patch,
        deadstemn_storage_patch  = cnveg_ns.deadstemn_storage_patch,
        livecrootn_storage_patch = cnveg_ns.livecrootn_storage_patch,
        deadcrootn_storage_patch = cnveg_ns.deadcrootn_storage_patch,
        leafc_patch  = cnveg_cs.leafc_patch,
        frootc_patch = cnveg_cs.frootc_patch,
        # prec10 is only READ when constrain_stress_deciduous_onset; supply a
        # correctly-sized zero buffer if the accumulator was never allocated
        # (use_cn=false paths that still construct this bundle).
        prec10_patch = (length(prec10_patch) >= length(cnveg_cs.leafc_patch) ?
                        prec10_patch :
                        fill!(similar(cnveg_cs.leafc_patch), zero(eltype(cnveg_cs.leafc_patch)))))

    col = _PhenStressSD_Col(;
        t_soisno_col = temperature.t_soisno_col,
        soilpsi_col  = soil_state.soilpsi_col,
        dayl         = gridcell.dayl)

    sc = _PhenStressSDScalars{T}(;
        dt              = T(dt),
        fracday         = T(fracday),
        ndays_off_val   = T(ndays_off_val),
        fstor2tran      = T(fstor2tran),
        crit_onset_fdd  = T(crit_onset_fdd),
        crit_onset_swi  = T(crit_onset_swi),
        soilpsi_on      = T(soilpsi_on),
        crit_offset_fdd = T(crit_offset_fdd),
        crit_offset_swi = T(crit_offset_swi),
        soilpsi_off     = T(soilpsi_off),
        secspqtrday     = T(secspqtrday),
        avg_dayspyr     = T(avg_dayspyr),
        tfrz            = T(TFRZ),
        secspday        = T(SECSPDAY),
        rain_threshold  = T(20.0),
        constrain_onset = cn_params.constrain_stress_deciduous_onset)

    # Struct-first kernel: manual backend (struct args carry no backend) + synchronize.
    isempty(mask_soilp) && return nothing
    backend = _kernel_backend(s1.offset_flag_patch)
    _phen_stress_decid_kernel!(backend)(
        s1, s2, inp, col,
        mask_soilp, patch_data.itype, patch_data.column, patch_data.gridcell,
        sc, soil_layer, varpar.nlevsno + soil_layer; ndrange = length(mask_soilp))
    KA.synchronize(backend)

    return nothing
end

# ==========================================================================
# get_swindow — Determine next sowing window
# ==========================================================================
function get_swindow(jday::Int, rx_starts::AbstractVector{Int},
                     rx_ends::AbstractVector{Int},
                     param_start::Int, param_end::Int)
    # No prescribed sowing windows
    if maximum(rx_starts) < 1
        return (w=1, start_w=param_start, end_w=param_end)
    end
    # Today is after latest sowing window end
    if jday > maximum(rx_ends)
        return (w=1, start_w=rx_starts[1], end_w=rx_ends[1])
    end
    # Find first window whose end >= today
    for w in eachindex(rx_starts)
        if min(rx_starts[w], rx_ends[w]) < 1
            break
        end
        if jday <= rx_ends[w]
            return (w=w, start_w=rx_starts[w], end_w=rx_ends[w])
        end
    end
    error("get_swindow(): No sowing window found")
end

# ==========================================================================
# was_sown_in_this_window
# ==========================================================================
function was_sown_in_this_window(sowing_window_startdate::Int,
                                 sowing_window_enddate::Int,
                                 jday::Int, idop::Int,
                                 sown_in_this_window::Bool)
    result = sown_in_this_window

    # Check if in sowing window
    is_in_sw = _is_doy_in_interval(sowing_window_startdate, sowing_window_enddate, jday)
    if !is_in_sw
        return false
    end
    # Check if planting date is in the window
    idop_in_sw = _is_doy_in_interval(sowing_window_startdate, sowing_window_enddate, idop)
    if is_in_sw && !idop_in_sw
        return false
    end
    # Check for same window vs different occurrence
    if sowing_window_startdate < sowing_window_enddate && idop > jday
        result = false
    elseif sowing_window_startdate > sowing_window_enddate
        if jday <= sowing_window_enddate && idop <= sowing_window_enddate && idop > jday
            result = false
        elseif jday >= sowing_window_startdate && (idop > jday || idop <= sowing_window_enddate)
            result = false
        end
    end
    return result
end

# Helper: is day-of-year in interval (wrapping around year boundary)
function _is_doy_in_interval(start_doy::Int, end_doy::Int, doy::Int)
    if start_doy <= end_doy
        return doy >= start_doy && doy <= end_doy
    else
        return doy >= start_doy || doy <= end_doy
    end
end

# ==========================================================================
# crop_phenology_init! — Crop-specific initialization
# ==========================================================================
function crop_phenology_init!(pstate::PhenologyState,
                              pftcon::PftConPhenology,
                              patch_data::PatchData,
                              gridcell::GridcellData,
                              npcropmin::Int, npcropmax::Int, maxveg::Int)
    np = length(patch_data.itype)
    pstate.inhemi = zeros(Int, np)

    pstate.minplantjday = fill(typemax(Int), maxveg + 1, inSH)  # 0:maxveg mapped to 1:maxveg+1
    pstate.maxplantjday = fill(typemax(Int), maxveg + 1, inSH)

    pstate.jdayyrstart = [1, 182]

    # Convert planting dates
    for n in npcropmin:npcropmax
        if n <= length(pftcon.is_pft_known_to_model) && pftcon.is_pft_known_to_model[n]
            pstate.minplantjday[n, inNH] = round(Int, pftcon.mnNHplantdate[n])
            pstate.maxplantjday[n, inNH] = round(Int, pftcon.mxNHplantdate[n])
            pstate.minplantjday[n, inSH] = round(Int, pftcon.mnSHplantdate[n])
            pstate.maxplantjday[n, inSH] = round(Int, pftcon.mxSHplantdate[n])
        end
    end

    # Determine hemisphere for each patch
    for p in 1:np
        g = patch_data.gridcell[p]
        if gridcell.latdeg[g] > 0.0
            pstate.inhemi[p] = inNH
        else
            pstate.inhemi[p] = inSH
        end
    end

    # Vernalization constants
    pstate.p1d   = 0.004
    pstate.p1v   = 0.003
    pstate.hti   = 1.0
    pstate.tbase = 0.0

    return nothing
end

# ==========================================================================
# crop_phenology! — Crop lifecycle management (simplified)
# ==========================================================================
@kernel function _phen_crop_phenology_kernel!(
        cphase_patch, bglfr_patch, bgtr_patch, lgsf_patch,
        onset_flag_patch, offset_flag_patch, onset_counter_patch,
        @Const(mask), @Const(itype), @Const(croplive_patch),
        @Const(hui_patch), @Const(huigrain_patch), @Const(leaf_long),
        dt, avg_dayspyr, secspday,
        cphase_planted_in, cphase_grainfill_in)
    T = eltype(bglfr_patch)
    p = @index(Global)
    @inbounds if mask[p]
        ivt = itype[p] + 1  # 0-based Fortran → 1-based Julia

        # reset background rates
        bglfr_patch[p] = zero(T)
        bgtr_patch[p]  = zero(T)
        lgsf_patch[p]  = zero(T)

        onset_flag_patch[p]  = zero(T)
        offset_flag_patch[p] = zero(T)

        if croplive_patch[p]
            cphase_patch[p] = cphase_planted_in
            onset_counter_patch[p] -= dt

            # grain fill phase: background litterfall
            if hui_patch[p] >= huigrain_patch[p]
                cphase_patch[p] = cphase_grainfill_in
                bglfr_patch[p] = one(T) / (leaf_long[ivt] * avg_dayspyr * secspday)
            end
        end
    end
end

function crop_phenology!(pstate::PhenologyState, params::PhenologyParams,
                         mask_pcropp::AbstractVector{Bool},
                         pftcon::PftConPhenology,
                         water_diag::WaterDiagnosticBulkData,
                         temperature::TemperatureData,
                         crop::CropData,
                         canopy_state::CanopyStateData,
                         cnveg_state::CNVegStateData,
                         cnveg_cs::CNVegCarbonStateData,
                         cnveg_ns::CNVegNitrogenStateData,
                         cnveg_cf::CNVegCarbonFluxData,
                         cnveg_nf::CNVegNitrogenFluxData,
                         patch_data::PatchData,
                         gridcell::GridcellData;
                         varctl::VarCtl=VarCtl(),
                         avg_dayspyr::Real=365.0,
                         jday::Int=1, kyr::Int=1, dayspyr::Real=365.0,
                         use_fertilizer::Bool=false)
    dt = pstate.dt
    T = eltype(cnveg_state.bglfr_patch)   # working precision (Float32 on Metal)

    _launch!(_phen_crop_phenology_kernel!,
        crop.cphase_patch, cnveg_state.bglfr_patch, cnveg_state.bgtr_patch,
        cnveg_state.lgsf_patch, cnveg_state.onset_flag_patch,
        cnveg_state.offset_flag_patch, cnveg_state.onset_counter_patch,
        mask_pcropp, patch_data.itype, crop.croplive_patch,
        crop.hui_patch, cnveg_state.huigrain_patch, pftcon.leaf_long,
        T(dt), T(avg_dayspyr), T(SECSPDAY),
        T(cphase_planted), T(cphase_grainfill))

    return nothing
end

# ==========================================================================
# crop_phase! — Get current crop phase
# ==========================================================================
function crop_phase!(mask_pcropp::AbstractVector{Bool}, crop::CropData,
                     cnveg_state::CNVegStateData,
                     crop_phase_out::AbstractVector{<:Real})
    # Per-patch independent: derive crop phase into output array.
    phen_crop_phase!(crop_phase_out, mask_pcropp, crop.croplive_patch,
                     crop.gddtsoi_patch, crop.hui_patch,
                     cnveg_state.huileaf_patch, cnveg_state.huigrain_patch,
                     cphase_planted, cphase_leafemerge, cphase_grainfill)
    return nothing
end

# ==========================================================================
# plant_crop! — Initialize crop at planting
# ==========================================================================
function plant_crop!(p::Int, leafcn_in::Real, jday::Int, kyr::Int,
                     crop::CropData, cnveg_state::CNVegStateData,
                     cnveg_cs::CNVegCarbonStateData,
                     cnveg_ns::CNVegNitrogenStateData,
                     cnveg_cf::CNVegCarbonFluxData,
                     cnveg_nf::CNVegNitrogenFluxData,
                     dt::Real)
    crop.croplive_patch[p]      = true
    crop.sown_in_this_window[p] = true
    cnveg_state.idop_patch[p]   = jday
    cnveg_state.iyop_patch[p]   = kyr
    crop.harvdate_patch[p]      = NOT_Harvested
    crop.sowing_count[p]       += 1

    seed = _initial_seed_at_planting[]
    cnveg_cs.leafc_xfer_patch[p] = seed
    cnveg_ns.leafn_xfer_patch[p] = seed / leafcn_in
    cnveg_cf.crop_seedc_to_leaf_patch[p] = seed / dt
    cnveg_nf.crop_seedn_to_leaf_patch[p] = (seed / leafcn_in) / dt

    # Initialize allocation coefficients
    cnveg_state.aleaf_patch[p]  = 1.0
    cnveg_state.aleafi_patch[p] = 1.0
    cnveg_state.astem_patch[p]  = 0.0
    cnveg_state.astemi_patch[p] = 0.0
    cnveg_state.aroot_patch[p]  = 0.0

    return nothing
end

# ==========================================================================
# days_past_planting — Calculate days since planting
# ==========================================================================
function days_past_planting(idop::Int, jday::Int; dayspyr::Int=365)
    if jday >= idop
        return jday - idop
    else
        return jday - idop + dayspyr
    end
end

# ==========================================================================
# vernalization! — Vernalization for winter cereal
# ==========================================================================
function vernalization!(p::Int, pstate::PhenologyState,
                        canopy_state::CanopyStateData,
                        temperature::TemperatureData,
                        water_diag::WaterDiagnosticBulkData,
                        cnveg_state::CNVegStateData,
                        crop::CropData,
                        patch_data::PatchData)
    c = patch_data.column[p]
    p1v   = pstate.p1v
    hti   = pstate.hti
    tbase = pstate.tbase

    t_ref2m     = temperature.t_ref2m_patch[p]
    t_ref2m_min = temperature.t_ref2m_min_patch[p]
    t_ref2m_max = temperature.t_ref2m_max_patch[p]
    snow_depth  = water_diag.snow_depth_col[c]

    force_harvest = false

    # crown temperature (smoothed TFRZ branch for AD)
    _tcrown_frozen = 2.0 + (t_ref2m - TFRZ) * (0.4 + 0.0018 *
        (smooth_min(snow_depth * 100.0, 15.0) - 15.0)^2)
    _tcrown_unfrozen = t_ref2m - TFRZ
    tcrown = smooth_ifelse(t_ref2m - TFRZ, _tcrown_unfrozen, _tcrown_frozen)

    # vernalization factor (smoothed TFRZ branches for AD)
    _hw_tmax_above_frz = smooth_heaviside(t_ref2m_max - TFRZ)
    _hw_tmin_below_15  = smooth_heaviside(TFRZ + 15.0 - t_ref2m_min)
    begin
        vd1 = 1.4 - 0.0778 * tcrown
        vd2 = 0.5 + 13.44 / ((t_ref2m_max - t_ref2m_min + 3.0)^2) * tcrown
        vd  = smooth_clamp(smooth_min(vd1, vd2), 0.0, 1.0)
        cnveg_state.cumvd_patch[p] += vd * _hw_tmax_above_frz * _hw_tmin_below_15
    end
    if cnveg_state.cumvd_patch[p] < 10.0
        cnveg_state.cumvd_patch[p] -= 0.5 * smooth_max(t_ref2m_max - TFRZ - 30.0, 0.0) * _hw_tmax_above_frz
    end
    cnveg_state.cumvd_patch[p] = smooth_max(0.0, cnveg_state.cumvd_patch[p])

    _vf_new = smooth_clamp(1.0 - p1v * (50.0 - cnveg_state.cumvd_patch[p]), 0.0, 1.0)
    crop.vf_patch[p] = _vf_new * _hw_tmax_above_frz + crop.vf_patch[p] * (1.0 - _hw_tmax_above_frz)

    # cold hardening (smoothed TFRZ branches for AD)
    hdidx = cnveg_state.hdidx_patch[p]
    _hw_thaw_decay = smooth_heaviside(t_ref2m_max - tbase - TFRZ - 10.0)
    _thaw_delta = smooth_max(t_ref2m_max - tbase - TFRZ - 10.0, 0.0)
    if t_ref2m_min <= TFRZ - 3.0 || hdidx != 0.0
        if hdidx >= hti
            hdidx += 0.083
            hdidx = smooth_min(hdidx, hti * 2.0)
        end
        hdidx -= 0.02 * _thaw_delta * _hw_thaw_decay
        if hdidx > hti
            hdidx -= 0.02 * _thaw_delta * _hw_thaw_decay
        end
        hdidx = smooth_max(0.0, hdidx) * _hw_thaw_decay + hdidx * (1.0 - _hw_thaw_decay)
    elseif tcrown >= tbase - 1.0
        if tcrown <= tbase + 8.0
            hdidx += 0.1 - (tcrown - tbase + 3.5)^2 / 506.0
            if hdidx >= hti && tcrown <= tbase + 0.0
                hdidx += 0.083
                hdidx = smooth_min(hdidx, hti * 2.0)
            end
        end
        hdidx -= 0.02 * _thaw_delta * _hw_thaw_decay
        if hdidx > hti
            hdidx -= 0.02 * _thaw_delta * _hw_thaw_decay
        end
        hdidx = smooth_max(0.0, hdidx) * _hw_thaw_decay + hdidx * (1.0 - _hw_thaw_decay)
    end
    cnveg_state.hdidx_patch[p] = hdidx

    # killing temperature check
    if t_ref2m_min <= TFRZ - 6.0
        tkil = (tbase - 6.0) - 6.0 * hdidx
        if tkil >= tcrown
            if (0.95 - 0.02 * (tcrown - tkil)^2) < 0.02 && canopy_state.tlai_patch[p] > 0.0
                force_harvest = true
            end
        end
    end

    return force_harvest
end

# ==========================================================================
# cn_onset_growth! — Transfer → display during onset
# ==========================================================================
function cn_onset_growth!(pstate::PhenologyState,
                          mask_soilp::AbstractVector{Bool},
                          pftcon::PftConPhenology,
                          cnveg_state::CNVegStateData,
                          cnveg_cs::CNVegCarbonStateData,
                          cnveg_ns::CNVegNitrogenStateData,
                          cnveg_cf::CNVegCarbonFluxData,
                          cnveg_nf::CNVegNitrogenFluxData,
                          patch_data::PatchData)
    dt = pstate.dt

    # Per-patch independent: transfer→display growth fluxes.
    phen_onset_growth!(
        cnveg_cf.leafc_xfer_to_leafc_patch, cnveg_cf.frootc_xfer_to_frootc_patch,
        cnveg_nf.leafn_xfer_to_leafn_patch, cnveg_nf.frootn_xfer_to_frootn_patch,
        cnveg_cf.livestemc_xfer_to_livestemc_patch, cnveg_cf.deadstemc_xfer_to_deadstemc_patch,
        cnveg_cf.livecrootc_xfer_to_livecrootc_patch, cnveg_cf.deadcrootc_xfer_to_deadcrootc_patch,
        cnveg_nf.livestemn_xfer_to_livestemn_patch, cnveg_nf.deadstemn_xfer_to_deadstemn_patch,
        cnveg_nf.livecrootn_xfer_to_livecrootn_patch, cnveg_nf.deadcrootn_xfer_to_deadcrootn_patch,
        mask_soilp, patch_data.itype, pftcon.woody,
        cnveg_state.onset_flag_patch, cnveg_state.onset_counter_patch, cnveg_state.bgtr_patch,
        cnveg_cs.leafc_xfer_patch, cnveg_cs.frootc_xfer_patch,
        cnveg_ns.leafn_xfer_patch, cnveg_ns.frootn_xfer_patch,
        cnveg_cs.livestemc_xfer_patch, cnveg_cs.deadstemc_xfer_patch,
        cnveg_cs.livecrootc_xfer_patch, cnveg_cs.deadcrootc_xfer_patch,
        cnveg_ns.livestemn_xfer_patch, cnveg_ns.deadstemn_xfer_patch,
        cnveg_ns.livecrootn_xfer_patch, cnveg_ns.deadcrootn_xfer_patch,
        dt)

    return nothing
end

# ==========================================================================
# cn_offset_litterfall! — Display → litter during offset
# ==========================================================================
function cn_offset_litterfall!(pstate::PhenologyState,
                               mask_soilp::AbstractVector{Bool},
                               pftcon::PftConPhenology,
                               cnveg_state::CNVegStateData,
                               cnveg_cs::CNVegCarbonStateData,
                               cnveg_ns::CNVegNitrogenStateData,
                               cnveg_cf::CNVegCarbonFluxData,
                               cnveg_nf::CNVegNitrogenFluxData,
                               crop::CropData,
                               patch_data::PatchData;
                               use_fun::Bool=false,
                               CNratio_floating::Bool=false,
                               for_testing_no_crop_seed_replenishment::Bool=false)
    dt = pstate.dt

    # Per-patch independent: display→litter fluxes during offset.
    phen_offset_litterfall!(
        cnveg_cf.leafc_to_litter_patch, cnveg_cf.frootc_to_litter_patch,
        cnveg_cf.leafc_to_litter_fun_patch,
        cnveg_cf.prev_leafc_to_litter_patch, cnveg_cf.prev_frootc_to_litter_patch,
        cnveg_nf.leafn_to_litter_patch, cnveg_nf.leafn_to_retransn_patch,
        cnveg_nf.frootn_to_litter_patch,
        mask_soilp, patch_data.itype,
        cnveg_state.offset_flag_patch, cnveg_state.offset_counter_patch,
        cnveg_cs.leafc_patch, cnveg_cs.frootc_patch,
        cnveg_ns.leafn_patch, cnveg_ns.frootn_patch,
        pftcon.lflitcn, pftcon.leafcn, pftcon.frootcn,
        dt, CNratio_floating, use_fun)

    return nothing
end

# ==========================================================================
# cn_background_litterfall! — Background litterfall
# ==========================================================================
function cn_background_litterfall!(pstate::PhenologyState,
                                   mask_soilp::AbstractVector{Bool},
                                   pftcon::PftConPhenology,
                                   cnveg_state::CNVegStateData,
                                   cnveg_cs::CNVegCarbonStateData,
                                   cnveg_ns::CNVegNitrogenStateData,
                                   cnveg_cf::CNVegCarbonFluxData,
                                   cnveg_nf::CNVegNitrogenFluxData,
                                   patch_data::PatchData;
                                   use_fun::Bool=false,
                                   CNratio_floating::Bool=false)
    # Per-patch independent: background litterfall fluxes.
    phen_background_litterfall!(
        cnveg_cf.leafc_to_litter_patch, cnveg_cf.frootc_to_litter_patch,
        cnveg_cf.leafc_to_litter_fun_patch,
        cnveg_nf.leafn_to_litter_patch, cnveg_nf.leafn_to_retransn_patch,
        cnveg_nf.frootn_to_litter_patch,
        mask_soilp, patch_data.itype, cnveg_state.bglfr_patch,
        cnveg_cs.leafc_patch, cnveg_cs.frootc_patch,
        cnveg_ns.leafn_patch, cnveg_ns.frootn_patch,
        pftcon.lflitcn, pftcon.leafcn, pftcon.frootcn, CNratio_floating, use_fun)

    return nothing
end

# ==========================================================================
# cn_livewood_turnover! — Live wood → dead wood turnover
# ==========================================================================
function cn_livewood_turnover!(pstate::PhenologyState,
                               mask_soilp::AbstractVector{Bool},
                               pftcon::PftConPhenology,
                               cnveg_cs::CNVegCarbonStateData,
                               cnveg_ns::CNVegNitrogenStateData,
                               cnveg_cf::CNVegCarbonFluxData,
                               cnveg_nf::CNVegNitrogenFluxData,
                               patch_data::PatchData;
                               CNratio_floating::Bool=false)
    lwtop = pstate.lwtop

    # Per-patch independent: live wood → dead wood turnover.
    phen_livewood_turnover!(
        cnveg_cf.livestemc_to_deadstemc_patch, cnveg_nf.livestemn_to_deadstemn_patch,
        cnveg_nf.livestemn_to_retransn_patch,
        cnveg_cf.livecrootc_to_deadcrootc_patch, cnveg_nf.livecrootn_to_deadcrootn_patch,
        cnveg_nf.livecrootn_to_retransn_patch,
        mask_soilp, patch_data.itype, pftcon.woody,
        pftcon.livewdcn, pftcon.deadwdcn,
        cnveg_cs.livestemc_patch, cnveg_ns.livestemn_patch,
        cnveg_cs.livecrootc_patch, cnveg_ns.livecrootn_patch,
        lwtop, CNratio_floating)

    return nothing
end

# ==========================================================================
# cn_crop_harvest_to_product_pools! — Crop harvest to product pools
# ==========================================================================
function cn_crop_harvest_to_product_pools!(mask_soilp::AbstractVector{Bool},
                                           mask_soilc::AbstractVector{Bool},
                                           cnveg_cf::CNVegCarbonFluxData,
                                           cnveg_nf::CNVegNitrogenFluxData,
                                           patch_data::PatchData;
                                           use_crop::Bool=false)
    if !use_crop
        return nothing
    end

    # Per-patch independent: harvest C/N to crop-product pools.
    phen_crop_harvest!(
        cnveg_cf.crop_harvestc_to_cropprodc_patch,
        cnveg_nf.crop_harvestn_to_cropprodn_patch, mask_soilp,
        cnveg_cf.leafc_to_biofuelc_patch, cnveg_cf.livestemc_to_biofuelc_patch,
        cnveg_cf.leafc_to_removedresiduec_patch, cnveg_cf.livestemc_to_removedresiduec_patch,
        cnveg_nf.leafn_to_biofueln_patch, cnveg_nf.livestemn_to_biofueln_patch,
        cnveg_nf.leafn_to_removedresiduen_patch, cnveg_nf.livestemn_to_removedresiduen_patch)

    return nothing
end

# ==========================================================================
# cn_litter_to_column! — Aggregate patch litter to column level
# ==========================================================================
@kernel function _phen_litter_to_column_kernel!(
        phenology_c_to_litr_c_col, phenology_n_to_litr_n_col,
        @Const(mask), @Const(itype), @Const(column), @Const(wtcol),
        @Const(lf_f), @Const(fr_f),
        @Const(leafc_to_litter_patch), @Const(frootc_to_litter_patch),
        @Const(leafn_to_litter_patch), @Const(frootn_to_litter_patch),
        @Const(livestemc_to_litter_patch), @Const(livestemn_to_litter_patch),
        @Const(leaf_prof), @Const(froot_prof),
        nlevdecomp::Int, i_litr_min::Int, i_litr_max::Int, npcropmin::Int)
    p = @index(Global)
    @inbounds if mask[p]
        ivt = itype[p] + 1  # 0-based Fortran → 1-based Julia
        c   = column[p]
        wt  = wtcol[p]
        ncol = size(phenology_c_to_litr_c_col, 1)
        for j in 1:nlevdecomp
            for i in i_litr_min:i_litr_max
                # leaf litter C and N
                if ncol >= c
                    _scatter_add!(phenology_c_to_litr_c_col, c, j, i,
                        leafc_to_litter_patch[p] * lf_f[ivt, i] * wt * leaf_prof[p, j])

                    _scatter_add!(phenology_n_to_litr_n_col, c, j, i,
                        leafn_to_litter_patch[p] * lf_f[ivt, i] * wt * leaf_prof[p, j])

                    # fine root litter C and N
                    _scatter_add!(phenology_c_to_litr_c_col, c, j, i,
                        frootc_to_litter_patch[p] * fr_f[ivt, i] * wt * froot_prof[p, j])

                    _scatter_add!(phenology_n_to_litr_n_col, c, j, i,
                        frootn_to_litter_patch[p] * fr_f[ivt, i] * wt * froot_prof[p, j])
                end
            end

            # crop stem litter uses leaf litter fractions
            if ivt >= npcropmin
                for i in i_litr_min:i_litr_max
                    if ncol >= c
                        _scatter_add!(phenology_c_to_litr_c_col, c, j, i,
                            livestemc_to_litter_patch[p] * lf_f[ivt, i] * wt * leaf_prof[p, j])
                        _scatter_add!(phenology_n_to_litr_n_col, c, j, i,
                            livestemn_to_litter_patch[p] * lf_f[ivt, i] * wt * leaf_prof[p, j])
                    end
                end
            end
        end # decomp level loop
    end
end

function cn_litter_to_column!(mask_soilp::AbstractVector{Bool},
                              pftcon::PftConPhenology,
                              cnveg_state::CNVegStateData,
                              cnveg_cf::CNVegCarbonFluxData,
                              cnveg_nf::CNVegNitrogenFluxData,
                              patch_data::PatchData,
                              leaf_prof::AbstractMatrix{<:Real},
                              froot_prof::AbstractMatrix{<:Real};
                              nlevdecomp::Int=1,
                              i_litr_min::Int=1, i_litr_max::Int=3,
                              npcropmin::Int=17,
                              use_grainproduct::Bool=false)
    _launch!(_phen_litter_to_column_kernel!,
        cnveg_cf.phenology_c_to_litr_c_col, cnveg_nf.phenology_n_to_litr_n_col,
        mask_soilp, patch_data.itype, patch_data.column, patch_data.wtcol,
        pftcon.lf_f, pftcon.fr_f,
        cnveg_cf.leafc_to_litter_patch, cnveg_cf.frootc_to_litter_patch,
        cnveg_nf.leafn_to_litter_patch, cnveg_nf.frootn_to_litter_patch,
        cnveg_cf.livestemc_to_litter_patch, cnveg_nf.livestemn_to_litter_patch,
        leaf_prof, froot_prof,
        nlevdecomp, i_litr_min, i_litr_max, npcropmin;
        ndrange = length(mask_soilp))
    return nothing
end
