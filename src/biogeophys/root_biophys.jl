# ==========================================================================
# Ported from: src/biogeophys/RootBiophysMod.F90
# Root biophysics: initialization of plant root profiles
#
# Contains subroutines for computing root fraction distributions using
# three methods:
#   - Zeng (2001) two-parameter exponential profile
#   - Jackson (1996) beta-function profile
#   - Koven exponential profile
#
# HISTORY: Original Fortran by Jinyun Tang, Mar 1, 2014
# ==========================================================================

# --- Rooting profile method constants ---
const ZENG_2001_ROOT    = 0  # Zeng 2001 root profile function
const JACKSON_1996_ROOT = 1  # Jackson 1996 root profile function
const KOVEN_EXP_ROOT    = 2  # Koven exponential root profile function

# --- Module-level rooting profile configuration ---

"""
    RootingProfileConfig

Configuration for rooting profile parameterization.
Selects the method and variant index for water and carbon root profiles.

Ported from module-level variables in `RootBiophysMod.F90`.
"""
Base.@kwdef mutable struct RootingProfileConfig
    rooting_profile_method_water::Int    = ZENG_2001_ROOT  # method for water root profile
    rooting_profile_method_carbon::Int   = ZENG_2001_ROOT  # method for carbon root profile
    rooting_profile_varindex_water::Int  = 1               # variant index for water (Jackson only)
    rooting_profile_varindex_carbon::Int = 2               # variant index for carbon (Jackson only)
end

const rooting_profile_config = RootingProfileConfig()

# --------------------------------------------------------------------------
# init_rootprof!
# --------------------------------------------------------------------------

"""
    init_rootprof!(config::RootingProfileConfig;
                   method_water::Int=ZENG_2001_ROOT,
                   method_carbon::Int=ZENG_2001_ROOT,
                   varindex_water::Int=1,
                   varindex_carbon::Int=2)

Initialize methods for root profile calculation.
In the Fortran version this reads from a namelist; here we accept keyword arguments.

Ported from `init_rootprof` in `RootBiophysMod.F90`.
"""
function init_rootprof!(config::RootingProfileConfig;
                        method_water::Int=ZENG_2001_ROOT,
                        method_carbon::Int=ZENG_2001_ROOT,
                        varindex_water::Int=1,
                        varindex_carbon::Int=2)
    config.rooting_profile_method_water    = method_water
    config.rooting_profile_method_carbon   = method_carbon
    config.rooting_profile_varindex_water  = varindex_water
    config.rooting_profile_varindex_carbon = varindex_carbon
    return nothing
end

# --------------------------------------------------------------------------
# init_vegrootfr!
# --------------------------------------------------------------------------

"""
    init_vegrootfr!(rootfr, col_zi, col_z, col_dz, col_nbedrock,
                    patch_column, patch_itype, patch_is_fates,
                    pftcon, config,
                    bounds_p, nlevsoi, nlevgrnd, water_carbon)

Initialize plant root fraction profiles.

The root fractions are computed per-patch and per-layer using one of three
methods (Zeng 2001, Jackson 1996, or Koven exponential), then roots below
bedrock are redistributed equally among the layers above bedrock.

Arguments:
- `rootfr`:          output root fraction matrix (npatches x nlevgrnd), modified in-place
- `col_zi`:          column interface depths (ncols x nlevels), soil-layer indexed
                     col_zi[c, lev] = bottom interface of soil layer lev;
                     col_zi[c, 0] is surface (accessed as lev-1 when lev=1 needs the 0 interface,
                     so we pass a pre-offset matrix or handle the zero-index externally)
- `col_z`:           column layer center depths (ncols x nlevels), soil-layer indexed
- `col_dz`:          column layer thicknesses (ncols x nlevels), soil-layer indexed
- `col_nbedrock`:    bedrock layer index per column (ncols)
- `patch_column`:    column index for each patch (npatches)
- `patch_itype`:     PFT type index for each patch (npatches, 0-based Fortran convention)
- `patch_is_fates`:  FATES flag per patch (npatches)
- `pftcon`:          PFT constants (PftconType)
- `config`:          rooting profile configuration (RootingProfileConfig)
- `bounds_p`:        patch index range
- `nlevsoi`:         number of hydrologically active soil layers
- `nlevgrnd`:        number of ground layers (soil + bedrock)
- `water_carbon`:    "water" or "carbon" to select which profile parameters to use

Ported from `init_vegrootfr` in `RootBiophysMod.F90`.
"""
function init_vegrootfr!(rootfr::Matrix{<:Real},
                         col_zi::Matrix{<:Real},
                         col_z::Matrix{<:Real},
                         col_dz::Matrix{<:Real},
                         col_nbedrock::Vector{Int},
                         patch_column::Vector{Int},
                         patch_itype::Vector{Int},
                         patch_is_fates::Vector{Bool},
                         pftcon::PftconType,
                         config::RootingProfileConfig,
                         bounds_p::UnitRange{Int},
                         nlevsoi::Int,
                         nlevgrnd::Int,
                         water_carbon::String)

    # Select method and variant based on water vs carbon
    if water_carbon == "water"
        rooting_profile_method = config.rooting_profile_method_water
        rooting_profile_varidx = config.rooting_profile_varindex_water
    elseif water_carbon == "carbon"
        rooting_profile_method = config.rooting_profile_method_carbon
        rooting_profile_varidx = config.rooting_profile_varindex_carbon
    else
        error("init_vegrootfr!: input type can only be water or carbon, got: $water_carbon")
    end

    # Compute root fractions for layers 1:nlevsoi using selected method
    if rooting_profile_method == ZENG_2001_ROOT
        zeng2001_rootfr!(rootfr, col_zi, patch_column, patch_itype,
                         patch_is_fates, pftcon, bounds_p, nlevsoi)
    elseif rooting_profile_method == JACKSON_1996_ROOT
        jackson1996_rootfr!(rootfr, col_zi, patch_column, patch_itype,
                            patch_is_fates, pftcon, bounds_p, nlevsoi,
                            rooting_profile_varidx)
    elseif rooting_profile_method == KOVEN_EXP_ROOT
        exponential_rootfr!(rootfr, col_z, col_dz, patch_column, patch_itype,
                            patch_is_fates, bounds_p, nlevsoi)
    else
        error("init_vegrootfr!: a root fraction function must be specified! Got method=$rooting_profile_method")
    end

    # Zero out root fractions below nlevsoi (bedrock layers)
    if nlevgrnd > nlevsoi && !isempty(bounds_p)
        root_zero_bedrock!(rootfr, first(bounds_p), last(bounds_p),
                           nlevsoi, nlevgrnd)
    end

    # Shift roots above bedrock boundary: redistribute roots below bedrock
    # equally among the layers above bedrock
    for p in bounds_p
        c = patch_column[p]
        nb = col_nbedrock[c]
        if nb < nlevsoi
            # Sum of root fractions in layers below bedrock but within nlevsoi
            root_below = 0.0
            for lev in (nb + 1):nlevsoi
                root_below += rootfr[p, lev]
            end
            # Redistribute equally to layers 1:nb
            if nb > 0
                add_per_layer = root_below / nb
                for lev in 1:nb
                    rootfr[p, lev] += add_per_layer
                end
            end
            # Zero out layers below bedrock
            for lev in (nb + 1):nlevsoi
                rootfr[p, lev] = 0.0
            end
        end
    end

    return nothing
end

# --------------------------------------------------------------------------
# zeng2001_rootfr!
# --------------------------------------------------------------------------

"""
    zeng2001_rootfr!(rootfr, col_zi, patch_column, patch_itype,
                     patch_is_fates, pftcon, bounds_p, ubj)

Compute root profile for soil water uptake using the Zeng (2001) equation
from J. Hydrometeorology.

Formula (computing from surface, d is depth in meters):
  Y = 1 - 0.5*(exp(-a*d) + exp(-b*d))
under the constraint that Y(d=0.1m) = 1 - beta^(10cm) and Y(d=d_obs) = 0.99
with beta and d_obs from Zeng et al. (1998).

Root fraction in layer lev:
  rootfr(lev) = 0.5*(exp(-a*zi(lev-1)) + exp(-b*zi(lev-1))
                    - exp(-a*zi(lev))   - exp(-b*zi(lev)))

For the bottom layer (ubj):
  rootfr(ubj) = 0.5*(exp(-a*zi(ubj-1)) + exp(-b*zi(ubj-1)))

Note on indexing: col_zi is indexed with soil-layer convention where
col_zi[c, lev] is the bottom interface depth of soil layer lev, and
col_zi[c, 0] (the surface, depth=0) is handled by using depth 0.0
explicitly when lev-1 = 0.

Ported from `zeng2001_rootfr` in `RootBiophysMod.F90`.
"""
function zeng2001_rootfr!(rootfr::Matrix{<:Real},
                          col_zi::Matrix{<:Real},
                          patch_column::Vector{Int},
                          patch_itype::Vector{Int},
                          patch_is_fates::Vector{Bool},
                          pftcon::PftconType,
                          bounds_p::UnitRange{Int},
                          ubj::Int)

    if !isempty(bounds_p)
        root_zeng2001!(rootfr, first(bounds_p), last(bounds_p), col_zi,
                       patch_column, patch_itype, patch_is_fates,
                       pftcon.roota_par, pftcon.rootb_par, ubj)
    end

    return nothing
end

# --------------------------------------------------------------------------
# jackson1996_rootfr!
# --------------------------------------------------------------------------

"""
    jackson1996_rootfr!(rootfr, col_zi, patch_column, patch_itype,
                        patch_is_fates, pftcon, bounds_p, ubj, varindx)

Compute root profile for soil water uptake using the Jackson et al. (1996)
equation from Oecologia.

Formula (computing from surface, d is depth in centimeters):
  Y = 1 - beta^d
Root fraction in layer lev:
  rootfr(lev) = beta^(zi(lev-1)*100) - beta^(zi(lev)*100)

where beta is a PFT-specific shape parameter from pftcon.rootprof_beta.

Ported from `jackson1996_rootfr` in `RootBiophysMod.F90`.
"""
function jackson1996_rootfr!(rootfr::Matrix{<:Real},
                             col_zi::Matrix{<:Real},
                             patch_column::Vector{Int},
                             patch_itype::Vector{Int},
                             patch_is_fates::Vector{Bool},
                             pftcon::PftconType,
                             bounds_p::UnitRange{Int},
                             ubj::Int,
                             varindx::Int)

    if !isempty(bounds_p)
        root_jackson1996!(rootfr, first(bounds_p), last(bounds_p), col_zi,
                          patch_column, patch_itype, patch_is_fates,
                          pftcon.rootprof_beta, ubj, varindx)
    end

    return nothing
end

# --------------------------------------------------------------------------
# exponential_rootfr!
# --------------------------------------------------------------------------

"""
    exponential_rootfr!(rootfr, col_z, col_dz, patch_column, patch_itype,
                        patch_is_fates, bounds_p, ubj)

Compute root profile for soil water/carbon uptake using the Koven
exponential equation.

The raw root density at each layer is:
  rootfr(lev) = exp(-rootprof_exp * z(lev)) * dz(lev)

Then normalized by:
  norm = -1/rootprof_exp * (exp(-rootprof_exp * z(ubj)) - 1)

where rootprof_exp = 3.0 (1/m) is the e-folding depth parameter.

Ported from `exponential_rootfr` in `RootBiophysMod.F90`.
"""
function exponential_rootfr!(rootfr::Matrix{<:Real},
                             col_z::Matrix{<:Real},
                             col_dz::Matrix{<:Real},
                             patch_column::Vector{Int},
                             patch_itype::Vector{Int},
                             patch_is_fates::Vector{Bool},
                             bounds_p::UnitRange{Int},
                             ubj::Int)

    rootprof_exp = 3.0  # e-folding depth parameter (1/m)

    if !isempty(bounds_p)
        root_exponential!(rootfr, first(bounds_p), last(bounds_p), col_z, col_dz,
                          patch_column, patch_is_fates, ubj, rootprof_exp)
    end

    return nothing
end

# --------------------------------------------------------------------------
# GPU kernels for root-fraction profiles (KernelAbstractions).
# Each (patch, layer) output element is independent: profile values depend
# only on per-patch parameters and per-column depth inputs (no loop-carried
# / accumulation dependence within these routines). The per-patch sum-based
# bedrock redistribution in init_vegrootfr! is intentionally NOT kernelized.
# --------------------------------------------------------------------------

# Zero out root fractions in bedrock layers (nlevsoi+1 .. nlevgrnd).
@kernel function _root_zero_bedrock_kernel!(rootfr, pmin::Int, levoff::Int)
    p0, j = @index(Global, NTuple)
    @inbounds rootfr[pmin + p0 - 1, levoff + j] = 0.0
end

"""
    root_zero_bedrock!(rootfr, pmin, pmax, nlevsoi, nlevgrnd)

Zero root fractions in bedrock layers `(nlevsoi+1):nlevgrnd` over patches
`pmin:pmax`. Backend-agnostic; one thread per (patch, bedrock layer).
"""
function root_zero_bedrock!(rootfr, pmin::Int, pmax::Int, nlevsoi::Int, nlevgrnd::Int)
    np = pmax - pmin + 1
    nbed = nlevgrnd - nlevsoi
    _launch!(_root_zero_bedrock_kernel!, rootfr, pmin, nlevsoi;
             ndrange = (np, nbed))
end

# Zeng (2001) two-parameter exponential root profile.
@kernel function _root_zeng2001_kernel!(rootfr, pmin::Int, @Const(col_zi),
                                        @Const(patch_column), @Const(patch_itype),
                                        @Const(patch_is_fates), @Const(roota_par),
                                        @Const(rootb_par), ubj::Int)
    p0, lev = @index(Global, NTuple)
    p = pmin + p0 - 1
    @inbounds begin
        if patch_is_fates[p]
            rootfr[p, lev] = 0.0
        else
            c = patch_column[p]
            pft_idx = patch_itype[p] + 1
            a = roota_par[pft_idx]
            b = rootb_par[pft_idx]
            if lev == ubj
                zi_upper_bot = col_zi[c, ubj - 1]
                rootfr[p, lev] = 0.5 * (exp(-a * zi_upper_bot) + exp(-b * zi_upper_bot))
            else
                zi_upper = lev == 1 ? 0.0 : col_zi[c, lev - 1]
                zi_lower = col_zi[c, lev]
                rootfr[p, lev] = 0.5 * (
                    exp(-a * zi_upper) + exp(-b * zi_upper) -
                    exp(-a * zi_lower) - exp(-b * zi_lower))
            end
        end
    end
end

"""
    root_zeng2001!(rootfr, pmin, pmax, col_zi, patch_column, patch_itype,
                   patch_is_fates, roota_par, rootb_par, ubj)

Zeng (2001) root profile over patches `pmin:pmax`, layers `1:ubj`.
Backend-agnostic; one thread per (patch, layer).
"""
function root_zeng2001!(rootfr, pmin::Int, pmax::Int, col_zi, patch_column,
                        patch_itype, patch_is_fates, roota_par, rootb_par, ubj::Int)
    np = pmax - pmin + 1
    _launch!(_root_zeng2001_kernel!, rootfr, pmin, col_zi, patch_column,
             patch_itype, patch_is_fates, roota_par, rootb_par, ubj;
             ndrange = (np, ubj))
end

# Jackson (1996) beta-function root profile.
@kernel function _root_jackson1996_kernel!(rootfr, pmin::Int, @Const(col_zi),
                                           @Const(patch_column), @Const(patch_itype),
                                           @Const(patch_is_fates), @Const(rootprof_beta),
                                           ubj::Int, varindx::Int)
    p0, lev = @index(Global, NTuple)
    p = pmin + p0 - 1
    @inbounds begin
        if patch_is_fates[p]
            rootfr[p, lev] = 0.0
        else
            c = patch_column[p]
            pft_idx = patch_itype[p] + 1
            beta = rootprof_beta[pft_idx, varindx]
            m_to_cm = 100.0
            zi_upper = lev == 1 ? 0.0 : col_zi[c, lev - 1]
            zi_lower = col_zi[c, lev]
            rootfr[p, lev] = beta^(zi_upper * m_to_cm) - beta^(zi_lower * m_to_cm)
        end
    end
end

"""
    root_jackson1996!(rootfr, pmin, pmax, col_zi, patch_column, patch_itype,
                      patch_is_fates, rootprof_beta, ubj, varindx)

Jackson (1996) root profile over patches `pmin:pmax`, layers `1:ubj`.
Backend-agnostic; one thread per (patch, layer).
"""
function root_jackson1996!(rootfr, pmin::Int, pmax::Int, col_zi, patch_column,
                           patch_itype, patch_is_fates, rootprof_beta, ubj::Int,
                           varindx::Int)
    np = pmax - pmin + 1
    _launch!(_root_jackson1996_kernel!, rootfr, pmin, col_zi, patch_column,
             patch_itype, patch_is_fates, rootprof_beta, ubj, varindx;
             ndrange = (np, ubj))
end

# Koven exponential root profile (norm recomputed per element from col_z[c,ubj]).
@kernel function _root_exponential_kernel!(rootfr, pmin::Int, @Const(col_z),
                                           @Const(col_dz), @Const(patch_column),
                                           @Const(patch_is_fates), ubj::Int,
                                           rootprof_exp)
    p0, lev = @index(Global, NTuple)
    p = pmin + p0 - 1
    @inbounds begin
        c = patch_column[p]
        norm = -1.0 / rootprof_exp * (exp(-rootprof_exp * col_z[c, ubj]) - 1.0)
        raw = patch_is_fates[p] ? 0.0 :
              exp(-rootprof_exp * col_z[c, lev]) * col_dz[c, lev]
        rootfr[p, lev] = raw / norm
    end
end

"""
    root_exponential!(rootfr, pmin, pmax, col_z, col_dz, patch_column,
                      patch_is_fates, ubj, rootprof_exp)

Koven exponential root profile over patches `pmin:pmax`, layers `1:ubj`.
The per-patch normalization is recomputed per element (depends only on the
input `col_z[c, ubj]`), so iterations remain independent. Backend-agnostic;
one thread per (patch, layer).
"""
function root_exponential!(rootfr, pmin::Int, pmax::Int, col_z, col_dz,
                           patch_column, patch_is_fates, ubj::Int, rootprof_exp)
    np = pmax - pmin + 1
    _launch!(_root_exponential_kernel!, rootfr, pmin, col_z, col_dz,
             patch_column, patch_is_fates, ubj, rootprof_exp;
             ndrange = (np, ubj))
end
