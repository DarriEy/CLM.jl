# ==========================================================================
# Ported from: src/main/clm_varpar.F90 (358 lines)
# Model dimension parameters
# ==========================================================================

# --- Fixed parameters (Fortran `parameter`) ---

const NLEV_EQUALSPACE  = 15    # number of equal-spaced levels
const TOPLEV_EQUALSPACE = 6    # top level for equal spacing
const NGASES           = 3     # CH4, O2, CO2
const NLEVCAN          = 1     # number of canopy layers
const NVEGWCS          = 4     # number of vegetation water conductance segments
const NUMWAT           = 5     # number of water types
const NUMRAD           = 2     # number of radiation bands (VIS, NIR)
const IVIS             = 1     # index for visible band
const INIR             = 2     # index for near-infrared band
const NUMSOLAR         = 2     # number of solar type bands
const NDST             = 4     # number of dust size classes
const DST_SRC_NBR      = 3     # number of source size distributions
const SZ_NBR           = 200   # number of sub-bin bins
const MXPFT            = 78    # maximum number of PFTs
const MXSOWINGS        = 1     # maximum number of crop sowings per year
const MXHARVESTS       = MXSOWINGS + 1  # maximum number of crop harvests per year
const NREPR            = 1     # number of crop reproductive pools (from CropReprPoolsMod, default=1)
const NLAYER           = 3     # number of VIC layers
const NVARIANTS        = 2     # number of soil moisture variants

# --- Vegetation pool indices ---
const NVEGPOOL_NATVEG = 18
const NVEGPOOL_CROP   = 3
const NVEG_RETRANSN   = 1

# Leaf pool indices
const ILEAF       = 1;  const ILEAF_ST    = 2;  const ILEAF_XF    = 3
const IFROOT      = 4;  const IFROOT_ST   = 5;  const IFROOT_XF   = 6
const ILIVESTEM   = 7;  const ILIVESTEM_ST = 8; const ILIVESTEM_XF = 9
const IDEADSTEM   = 10; const IDEADSTEM_ST = 11; const IDEADSTEM_XF = 12
const ILIVECROOT  = 13; const ILIVECROOT_ST = 14; const ILIVECROOT_XF = 15
const IDEADCROOT  = 16; const IDEADCROOT_ST = 17; const IDEADCROOT_XF = 18
const IGRAIN      = 19; const IGRAIN_ST   = 20; const IGRAIN_XF   = 21

# Litter pool indices (set in clm_varpar_init!)
const I_LITR1 = 1  # base litter index

"""
Runtime-computed model parameters.
These depend on configuration (use_crop, use_cn, etc.) and are set by `varpar_init!()`.
Mirrors variables computed in `clm_varpar_init` in Fortran.
"""
Base.@kwdef mutable struct VarPar
    nlevsoi::Int = 0                 # number of soil levels
    nlevsoifl::Int = 0               # number of soil levels (full)
    nlevgrnd::Int = 0                # number of ground levels (soil + bedrock)
    nlevurb::Int = 0                 # number of urban levels
    nlevmaxurbgrnd::Int = 0          # max urban ground levels
    nlevlak::Int = 0                 # number of lake levels
    nlevdecomp::Int = 0              # number of decomposition levels
    nlevdecomp_full::Int = 0         # full decomposition levels
    nlevsno::Int = -1                # number of snow levels (negative = computed)
    mxharvests::Int = 0              # max harvests
    nlayert::Int = 0                 # total number of layers
    nvegcpool::Int = 0               # number of veg carbon pools
    nvegnpool::Int = 0               # number of veg nitrogen pools
    maxveg::Int = 0                  # max vegetation types
    maxpatch_urb::Int = 5            # max urban patches
    maxsoil_patches::Int = 0         # max soil patches
    maxpatch_glc::Int = 0            # max glacier patches

    # Litter pool indices (runtime-dependent)
    i_litr2::Int = -9
    i_litr3::Int = -9
    i_litr_min::Int = -9
    i_litr_max::Int = -9
    i_met_lit::Int = -9
    i_cop_mic::Int = -9
    i_oli_mic::Int = -9
    i_cwd::Int = -9
    i_cwdl2::Int = -9

    # Decomposition pool counts
    ndecomp_pools_max::Int = 0
    ndecomp_pools::Int = 0
    ndecomp_cascade_transitions::Int = 0
    ndecomp_cascade_outtransitions::Int = 0
    ndecomp_pools_vr::Int = 0

    # PFT/CFT ranges
    natpft_lb::Int = 0
    natpft_ub::Int = 0
    natpft_size::Int = 0
    surfpft_lb::Int = 0
    surfpft_ub::Int = 0
    cft_lb::Int = 0
    cft_ub::Int = 0
    cft_size::Int = 0

    # Retranslocated N index
    iretransn::Int = 0
    ioutc::Int = 0
    ioutn::Int = 0
end

# Global parameters instance
const varpar = VarPar()

"""
    varpar_init!(vp::VarPar, actual_maxsoil_patches, surf_numpft, surf_numcft, actual_nlevurb)

Initialize runtime-dependent parameters. Called once during model setup.
Ported from: clm_varpar_init() in clm_varpar.F90
"""
function varpar_init!(vp::VarPar, actual_maxsoil_patches::Int,
                      surf_numpft::Int, surf_numcft::Int, actual_nlevurb::Int)
    vp.maxsoil_patches = actual_maxsoil_patches
    vp.nlevurb = actual_nlevurb

    # Soil levels depend on layer structure
    if varctl.soil_layerstruct_predefined == "10SL_3.5m"
        vp.nlevsoi = 10
        vp.nlevgrnd = 15
    elseif varctl.soil_layerstruct_predefined == "23SL_3.5m"
        vp.nlevsoi = 20
        vp.nlevgrnd = 25
    elseif varctl.soil_layerstruct_predefined == "49SL_10m"
        vp.nlevsoi = 49
        vp.nlevgrnd = 54
    elseif varctl.soil_layerstruct_predefined == "20SL_8.5m"
        vp.nlevsoi = 20
        vp.nlevgrnd = 25
    elseif varctl.soil_layerstruct_predefined == "4SL_2m"
        vp.nlevsoi = 4
        vp.nlevgrnd = 4
    elseif varctl.soil_layerstruct_predefined == "UNSET"
        # Default: standard 20-layer CLM structure
        vp.nlevsoi = 20
        vp.nlevgrnd = 25
    else
        error("Unknown soil_layerstruct_predefined: $(varctl.soil_layerstruct_predefined)")
    end

    vp.nlevsoifl = vp.nlevsoi
    vp.nlevmaxurbgrnd = max(vp.nlevurb, vp.nlevgrnd)

    # Lake levels
    vp.nlevlak = varctl.use_extralakelayers ? 25 : 10

    # Snow levels
    if vp.nlevsno == -1
        vp.nlevsno = 12  # default max snow layers
    end

    # Decomposition levels — vertically-resolved soil C/N (the `_vr` pools).
    # Matches Fortran clm_varpar.F90:302-310, which for the BGC-active branch
    # (use_cn .or. (use_fates .and. .not. use_fates_sp)) sets:
    #     nlevdecomp = nlevsoi (20), nlevdecomp_full = nlevgrnd (25)
    # The Fortran spinup restart confirms it (depth-decreasing _vr pools).
    #
    # The old "nlevdecomp=1 placeholder / NaNs pools 2-7" TODO here was resolved
    # by 1348de3 (which flipped these to the vertical values and fixed the three
    # unpopulated-param NaN bugs behind it) — the comment simply outlived the fix.
    # Verified: clm_drv! with use_cn=true runs the full vertical decomposition and
    # litter-transport chain over all 7 decomp pools × 20 layers with zero NaN.
    #
    # Known (numerically inert) deviation: Fortran's else-branch drops to
    # nlevdecomp = nlevdecomp_full = 1 when BGC is off. We size unconditionally to
    # the vertical values instead. With BGC off nothing reads the _vr pools, so this
    # is pure allocation headroom, not a physics difference; gating on use_cn would
    # silently collapse the BGC test suite to a single layer, since those tests call
    # varpar_init! against the varctl default (use_cn=false) and then read
    # varpar.nlevdecomp. Revisit only if non-BGC memory footprint becomes a concern.
    vp.nlevdecomp = vp.nlevsoi
    vp.nlevdecomp_full = vp.nlevgrnd

    # PFT/CFT ranges
    vp.natpft_lb = 0
    vp.natpft_ub = surf_numpft
    vp.natpft_size = vp.natpft_ub - vp.natpft_lb + 1

    # Fortran clm_varpar.F90:209 branches on create_crop_landunit, NOT use_crop,
    # and comments its else-branch "only true when FATES is active". Branching on
    # use_crop here gave a default non-crop, non-FATES run cft_size=0 where CTSM
    # gives surf_numcft. See docs/DRIVER_DEFAULTS_AUDIT.md M5.
    if varctl.create_crop_landunit
        vp.cft_lb = vp.natpft_ub + 1
        vp.cft_ub = vp.cft_lb + surf_numcft - 1
        vp.cft_size = vp.cft_ub - vp.cft_lb + 1
    else
        vp.cft_lb = vp.natpft_ub + 1
        vp.cft_ub = vp.cft_lb
        vp.cft_size = 0
    end

    vp.surfpft_lb = vp.natpft_lb
    vp.surfpft_ub = max(vp.natpft_ub, vp.cft_ub)
    vp.maxveg = vp.surfpft_ub

    # Vegetation pools
    if varctl.use_crop
        vp.nvegcpool = NVEGPOOL_NATVEG + NVEGPOOL_CROP
        vp.nvegnpool = NVEGPOOL_NATVEG + NVEGPOOL_CROP
    else
        vp.nvegcpool = NVEGPOOL_NATVEG
        vp.nvegnpool = NVEGPOOL_NATVEG
    end

    vp.iretransn = vp.nvegnpool + NVEG_RETRANSN
    vp.ioutc = vp.nvegcpool + 1
    vp.ioutn = vp.nvegnpool + NVEG_RETRANSN + 1

    # Glacier patches
    vp.maxpatch_glc = 10

    # Total layers
    vp.nlayert = vp.nlevsoi + vp.nlevlak + NLAYER

    nothing
end
