# ==========================================================================
# Ported from: src/biogeophys/SoilMoistStressMod.F90
# Calculates soil moisture stress for plant GPP and transpiration.
# ==========================================================================

# --- Module-level control state ---

const MOIST_STRESS_CLM_DEFAULT = 0

Base.@kwdef mutable struct SoilMoistStressControl
    root_moist_stress_method::Int = MOIST_STRESS_CLM_DEFAULT
    perchroot::Bool = false       # true => btran is based only on unfrozen soil levels
    perchroot_alt::Bool = false   # true => btran is based on active layer (defined over two years)
end

const soil_moist_stress_ctrl = SoilMoistStressControl()

# --------------------------------------------------------------------------
# AD-smoothing of the btran (root soil-moisture stress) discontinuities.
#
# The CLM4.5-default btran kernel (`_smstress_rootr_btran_kernel!`) has a HARD
# per-layer gate that zeroes the layer's root-water contribution the instant the
# layer's liquid water hits 0 or its temperature drops to tfrz-2 K, plus several
# hard clamps (s_node floor/ceil, smp_node floor, rresis ceil, rootr floor). As
# T or soil-water sweep across these thresholds, btran[p] — the factor that
# multiplies stomatal conductance / photosynthesis — has KINKS, so its forward-AD
# derivative jumps discontinuously (and reverse-AD can NaN at the gate).
#
# The smoothing is gated PURELY by element type via the shared `_use_smooth`
# (smooth_ad.jl): plain Float64 evaluates the EXACT hard physics (byte-identical
# default; full forward suite unchanged), ForwardDiff.Dual evaluates the smooth
# surrogate. `BTRAN_SMOOTH_K` is the logistic/softmax sharpness; the transition
# width ε ≈ 1/k, and as k → ∞ (ε → 0) the surrogate → the hard physics. Separate
# widths for the two gate axes because they live in different units (volumetric
# water O(0.1) vs temperature O(1 K)); both default sharp enough to be physically
# negligible while keeping a finite, continuous derivative across the threshold.
const BTRAN_SMOOTH_K        = Ref(50.0)   # clamps (s_node / smp_node / rresis / rootr floor)
const BTRAN_GATE_K_WATER    = Ref(200.0)  # liquid-water gate: σ(k·h2osoi_liqvol)
const BTRAN_GATE_K_TEMP     = Ref(5.0)    # temperature gate: σ(k·(t_soisno-(tfrz-2)))

"""
    init_root_moist_stress!()

Specify the method to compute root soil moisture stress.
Ported from `init_root_moist_stress` in `SoilMoistStressMod.F90`.
"""
function init_root_moist_stress!()
    soil_moist_stress_ctrl.root_moist_stress_method = MOIST_STRESS_CLM_DEFAULT
    return nothing
end

"""
    set_perchroot_opt!(perchroot_global, perchroot_alt_global)

Set up local perchroot logical switches.
Ported from `set_perchroot_opt` in `SoilMoistStressMod.F90`.
"""
function set_perchroot_opt!(perchroot_global::Bool, perchroot_alt_global::Bool)
    soil_moist_stress_ctrl.perchroot = perchroot_global
    soil_moist_stress_ctrl.perchroot_alt = perchroot_alt_global
    return nothing
end

# --------------------------------------------------------------------------
# calc_effective_soilporosity!
# --------------------------------------------------------------------------

"""
    calc_effective_soilporosity!(watsat, h2osoi_ice, col_dz, eff_por,
                                 mask, bounds_col, nlevgrnd, nlevsno)

Compute the effective soil porosity.
Ported from `calc_effective_soilporosity` in `SoilMoistStressMod.F90`.

Arguments:
- `watsat`: soil porosity (nc × nlevgrnd), soil-only indexing
- `h2osoi_ice`: ice water content kg/m2 (nc × nlevsno+nlevgrnd), combined indexing
- `col_dz`: layer thickness (nc × nlevsno+nlevgrnd), combined indexing
- `eff_por`: effective porosity output (nc × nlevgrnd), soil-only indexing
- `mask`: BitVector mask for active columns
- `bounds_col`: column index range
- `nlevgrnd`: number of ground levels (ubj)
- `nlevsno`: number of snow levels (for index offset)
"""
# --- Kernel: effective soil porosity over (column, soil layer) ---
@kernel function _smstress_effsoilpor_kernel!(eff_por, @Const(mask), @Const(watsat),
                                              @Const(h2osoi_ice), @Const(col_dz),
                                              joff::Int, denice, cmin::Int, cmax::Int)
    c, j = @index(Global, NTuple)
    @inbounds if cmin <= c <= cmax && mask[c]
        # compute the volumetric ice content
        vol_ice = min(watsat[c, j], h2osoi_ice[c, j + joff] / (denice * col_dz[c, j + joff]))
        # compute the maximum soil space to fill liquid water and air
        eff_por[c, j] = watsat[c, j] - vol_ice
    end
end

smstress_effsoilpor!(eff_por, mask, watsat, h2osoi_ice, col_dz,
                     bounds_col::UnitRange{Int}, joff::Int, nlevgrnd::Int, denice) =
    _launch!(_smstress_effsoilpor_kernel!, eff_por, mask, watsat, h2osoi_ice, col_dz,
             joff, eltype(eff_por)(denice), first(bounds_col), last(bounds_col);
             ndrange = (length(mask), nlevgrnd))

function calc_effective_soilporosity!(watsat::AbstractMatrix{<:Real},
                                      h2osoi_ice::AbstractMatrix{<:Real},
                                      col_dz::AbstractMatrix{<:Real},
                                      eff_por::AbstractMatrix{<:Real},
                                      mask::AbstractVector{Bool},
                                      bounds_col::UnitRange{Int},
                                      nlevgrnd::Int,
                                      nlevsno::Int)
    joff = nlevsno  # offset for combined snow+soil indexing
    smstress_effsoilpor!(eff_por, mask, watsat, h2osoi_ice, col_dz,
                         bounds_col, joff, nlevgrnd, DENICE)
    return nothing
end

# --------------------------------------------------------------------------
# calc_effective_snowporosity!
# --------------------------------------------------------------------------

"""
    calc_effective_snowporosity!(h2osoi_ice, col_dz, jtop, eff_por,
                                 mask, bounds_col, lbj, nlevsno)

Compute the effective porosity of snow.
Ported from `calc_effective_snowporosity` in `SoilMoistStressMod.F90`.

Arguments:
- `h2osoi_ice`: ice water content kg/m2 (nc × nlevsno+nlevgrnd), combined indexing
- `col_dz`: layer thickness (nc × nlevsno+nlevgrnd), combined indexing
- `jtop`: top level for each column (Fortran snow layer index, e.g. -4 to 0)
- `eff_por`: effective porosity output (nc × nlevsno), snow-only indexing (1:nlevsno maps to lbj:0)
- `mask`: BitVector mask for active columns
- `bounds_col`: column index range
- `lbj`: lower bound snow layer index (Fortran, e.g. -nlevsno+1)
- `nlevsno`: number of snow levels
"""
# --- Kernel: effective snow porosity over (column, snow layer) ---
# jj is the Julia 1-based combined index (1..nlevsno); the Fortran snow index is
# j = jj - joff. Only layers at or below the column top (j >= jtop) are written.
@kernel function _smstress_effsnowpor_kernel!(eff_por, @Const(mask), @Const(jtop),
                                              @Const(h2osoi_ice), @Const(col_dz),
                                              joff::Int, denice)
    c, jj = @index(Global, NTuple)
    @inbounds if mask[c]
        T = eltype(eff_por)
        j = jj - joff  # Fortran snow layer index (negative)
        if j >= jtop[c]
            # compute the volumetric ice content
            vol_ice = min(one(T), h2osoi_ice[c, jj] / (denice * col_dz[c, jj]))
            # compute the maximum snow void space to fill liquid water and air
            eff_por[c, jj] = one(T) - vol_ice
        end
    end
end

smstress_effsnowpor!(eff_por, mask, jtop, h2osoi_ice, col_dz, joff::Int, nlevsno::Int, denice) =
    _launch!(_smstress_effsnowpor_kernel!, eff_por, mask, jtop, h2osoi_ice, col_dz,
             joff, eltype(eff_por)(denice); ndrange = (length(mask), nlevsno))

function calc_effective_snowporosity!(h2osoi_ice::AbstractMatrix{<:Real},
                                      col_dz::AbstractMatrix{<:Real},
                                      jtop::AbstractVector{<:Integer},
                                      eff_por::AbstractMatrix{<:Real},
                                      mask::AbstractVector{Bool},
                                      bounds_col::UnitRange{Int},
                                      lbj::Int,
                                      nlevsno::Int)
    joff = nlevsno  # offset: Fortran j → Julia index j + nlevsno
    # Fortran snow layers run lbj:0, i.e. Julia combined indices 1:nlevsno.
    smstress_effsnowpor!(eff_por, mask, jtop, h2osoi_ice, col_dz, joff, nlevsno, DENICE)
    return nothing
end

# --------------------------------------------------------------------------
# calc_volumetric_h2oliq!
# --------------------------------------------------------------------------

"""
    calc_volumetric_h2oliq!(eff_porosity, h2osoi_liq, col_dz, jtop, vol_liq,
                             mask, bounds_col, lbj, ubj, nlevsno)

Compute the volumetric liquid water content.
Ported from `calc_volumetric_h2oliq` in `SoilMoistStressMod.F90`.

Arguments:
- `eff_porosity`: effective soil porosity (nc × dim), indexing depends on lbj
- `h2osoi_liq`: liquid water content kg/m2 (nc × nlevsno+nlevgrnd), combined indexing
- `col_dz`: layer thickness (nc × nlevsno+nlevgrnd), combined indexing
- `jtop`: top level for each column (Fortran snow layer index)
- `vol_liq`: volumetric liquid water content output (nc × dim)
- `mask`: BitVector mask for active columns
- `bounds_col`: column index range
- `lbj`: lower bound layer index (Fortran)
- `ubj`: upper bound layer index (Fortran)
- `nlevsno`: number of snow levels (for index offset)
"""
# --- Kernel: volumetric liquid water content over (column, layer) ---
# jrel is 1-based over the layer span lbj:ubj; the combined Julia index is
# jj = jrel + lbj - 1 + joff and the Fortran layer index is j = jj - joff.
# Only layers at or below the column top (j >= jtop) are written.
@kernel function _smstress_volh2oliq_kernel!(vol_liq, @Const(mask), @Const(jtop),
                                             @Const(eff_porosity), @Const(h2osoi_liq),
                                             @Const(col_dz), lbj::Int, joff::Int, denh2o)
    c, jrel = @index(Global, NTuple)
    @inbounds if mask[c]
        j = lbj + (jrel - 1)            # Fortran layer index
        jj = j + joff                   # Julia 1-based combined index
        if j >= jtop[c]
            # volume of liquid is no greater than effective void space
            vol_liq[c, jj] = min(eff_porosity[c, jj], h2osoi_liq[c, jj] / (col_dz[c, jj] * denh2o))
        end
    end
end

function smstress_volh2oliq!(vol_liq, mask, jtop, eff_porosity, h2osoi_liq, col_dz,
                             lbj::Int, ubj::Int, joff::Int, denh2o)
    nlay = ubj - lbj + 1
    _launch!(_smstress_volh2oliq_kernel!, vol_liq, mask, jtop, eff_porosity, h2osoi_liq,
             col_dz, lbj, joff, eltype(vol_liq)(denh2o); ndrange = (length(mask), nlay))
end

function calc_volumetric_h2oliq!(eff_porosity::AbstractMatrix{<:Real},
                                  h2osoi_liq::AbstractMatrix{<:Real},
                                  col_dz::AbstractMatrix{<:Real},
                                  jtop::AbstractVector{<:Integer},
                                  vol_liq::AbstractMatrix{<:Real},
                                  mask::AbstractVector{Bool},
                                  bounds_col::UnitRange{Int},
                                  lbj::Int,
                                  ubj::Int,
                                  nlevsno::Int)
    joff = nlevsno  # offset: Fortran j → Julia index j + nlevsno
    smstress_volh2oliq!(vol_liq, mask, jtop, eff_porosity, h2osoi_liq, col_dz,
                        lbj, ubj, joff, DENH2O)
    return nothing
end

# --------------------------------------------------------------------------
# array_normalization! (from SimpleMathMod)
# --------------------------------------------------------------------------

# --- Kernel: per-row (patch) sum-to-one normalization with internal layer loop ---
# One thread per patch; the row sum and division are an inter-layer reduction, so
# they stay sequential inside the thread (cf. the soil_temperature per-column
# kernels with internal level loops). Literals carried as the array eltype.
@kernel function _smstress_array_norm_kernel!(arr, @Const(mask), nlevgrnd::Int)
    p = @index(Global)
    @inbounds if mask[p]
        T = eltype(arr)
        rootsum = zero(T)
        for j in 1:nlevgrnd
            rootsum += arr[p, j]
        end
        if rootsum > zero(T)
            for j in 1:nlevgrnd
                arr[p, j] = arr[p, j] / rootsum
            end
        end
    end
end

"""
    array_normalization!(arr, mask, bounds, nlevgrnd)

Normalize each row of `arr` (patch × nlevgrnd) so that its sum equals 1.
If the sum is zero, the row is left as zeros.
Ported from `array_normalization` in `SimpleMathMod.F90`.
"""
function array_normalization!(arr::AbstractMatrix{<:Real},
                               mask::AbstractVector{Bool},
                               bounds::UnitRange{Int},
                               nlevgrnd::Int)
    _launch!(_smstress_array_norm_kernel!, arr, mask, nlevgrnd; ndrange = length(mask))
    return nothing
end

# --------------------------------------------------------------------------
# normalize_unfrozen_rootfr!
# --------------------------------------------------------------------------

"""
    normalize_unfrozen_rootfr!(rootfr, t_soisno, altmax_lastyear_indx, altmax_indx,
                                patch_column, rootfr_unf,
                                mask_patch, bounds_patch, nlevgrnd, nlevsno)

Normalize root fraction for total unfrozen depth.
Ported from `normalize_unfrozen_rootfr` in `SoilMoistStressMod.F90`.
"""
# --- Kernel: unfrozen rootfr from active layer depth, over (patch, layer) ---
@kernel function _smstress_rootfr_unf_alt_kernel!(rootfr_unf, @Const(mask_patch),
                                                  @Const(patch_column), @Const(rootfr),
                                                  @Const(altmax_lastyear_indx),
                                                  @Const(altmax_indx))
    p, j = @index(Global, NTuple)
    @inbounds if mask_patch[p]
        T = eltype(rootfr_unf)
        c = patch_column[p]
        # active-layer index threshold (Real comparison; literal as the row eltype)
        thresh = max(altmax_lastyear_indx[c], altmax_indx[c], one(eltype(altmax_indx)))
        if j <= thresh
            rootfr_unf[p, j] = rootfr[p, j]
        else
            rootfr_unf[p, j] = zero(T)
        end
    end
end

smstress_rootfr_unf_alt!(rootfr_unf, mask_patch, patch_column, rootfr,
                         altmax_lastyear_indx, altmax_indx, nlevgrnd::Int) =
    _launch!(_smstress_rootfr_unf_alt_kernel!, rootfr_unf, mask_patch, patch_column,
             rootfr, altmax_lastyear_indx, altmax_indx;
             ndrange = (length(mask_patch), nlevgrnd))

# --- Kernel: unfrozen rootfr from instantaneous temperature, over (patch, layer) ---
@kernel function _smstress_rootfr_unf_temp_kernel!(rootfr_unf, @Const(mask_patch),
                                                   @Const(patch_column), @Const(rootfr),
                                                   @Const(t_soisno), joff::Int, tfrz)
    p, j = @index(Global, NTuple)
    @inbounds if mask_patch[p]
        T = eltype(rootfr_unf)
        c = patch_column[p]
        if t_soisno[c, j + joff] >= tfrz
            rootfr_unf[p, j] = rootfr[p, j]
        else
            rootfr_unf[p, j] = zero(T)
        end
    end
end

smstress_rootfr_unf_temp!(rootfr_unf, mask_patch, patch_column, rootfr, t_soisno,
                          joff::Int, nlevgrnd::Int, tfrz) =
    _launch!(_smstress_rootfr_unf_temp_kernel!, rootfr_unf, mask_patch, patch_column,
             rootfr, t_soisno, joff, eltype(rootfr_unf)(tfrz);
             ndrange = (length(mask_patch), nlevgrnd))

function normalize_unfrozen_rootfr!(rootfr::AbstractMatrix{<:Real},
                                     t_soisno::AbstractMatrix{<:Real},
                                     altmax_lastyear_indx::AbstractVector{<:Real},
                                     altmax_indx::AbstractVector{<:Real},
                                     patch_column::AbstractVector{<:Integer},
                                     rootfr_unf::AbstractMatrix{<:Real},
                                     mask_patch::AbstractVector{Bool},
                                     bounds_patch::UnitRange{Int},
                                     nlevgrnd::Int,
                                     nlevsno::Int)
    joff = nlevsno
    perchroot = soil_moist_stress_ctrl.perchroot
    perchroot_alt = soil_moist_stress_ctrl.perchroot_alt

    if perchroot || perchroot_alt
        if perchroot_alt
            # use total active layer (max thaw depth for current and prior year)
            smstress_rootfr_unf_alt!(rootfr_unf, mask_patch, patch_column, rootfr,
                                     altmax_lastyear_indx, altmax_indx, nlevgrnd)
        else
            # use instantaneous temperature
            smstress_rootfr_unf_temp!(rootfr_unf, mask_patch, patch_column, rootfr,
                                      t_soisno, joff, nlevgrnd, TFRZ)
        end
    end

    # normalize the root fraction for each pft
    array_normalization!(rootfr_unf, mask_patch, bounds_patch, nlevgrnd)

    return nothing
end

# --------------------------------------------------------------------------
# soil_suction (inline replacement for SoilWaterRetentionCurveMod)
# --------------------------------------------------------------------------

"""
    soil_suction_clapp_hornberger(sucsat, s_node, bsw)

Compute soil matric potential using the Clapp-Hornberger parameterization.
Returns smp_node (matric potential in mm).
Ported from `soil_suction` in `SoilWaterRetentionCurveMod.F90`.
"""
function soil_suction_clapp_hornberger(sucsat::Real, s_node::Real, bsw::Real)
    return -sucsat * s_node^(-bsw)
end

# --------------------------------------------------------------------------
# calc_root_moist_stress_clm45default!
# --------------------------------------------------------------------------

# --- Kernel: normalize root resistances by btran, over (patch, layer) ---
@kernel function _smstress_normalize_rootr_kernel!(rootr, @Const(mask_patch),
                                                   @Const(btran), btran0)
    p, j = @index(Global, NTuple)
    @inbounds if mask_patch[p]
        T = eltype(rootr)
        if btran[p] > btran0
            rootr[p, j] = rootr[p, j] / btran[p]
        else
            rootr[p, j] = zero(T)
        end
    end
end

smstress_normalize_rootr!(rootr, mask_patch, btran, btran0, nlevgrnd::Int) =
    _launch!(_smstress_normalize_rootr_kernel!, rootr, mask_patch, btran,
             eltype(rootr)(btran0); ndrange = (length(mask_patch), nlevgrnd))

# --- Kernel: per-patch root resistance + btran accumulation (CLM4.5 default) ---
# One thread per patch. btran[p] is accumulated across layers (an inter-layer
# reduction), so the layer loop stays sequential inside the thread; rresis/rootr
# are written per layer. All literals/physical consts carried as the output
# eltype so no Float64 reaches a Float32-only backend.
@kernel function _smstress_rootr_btran_kernel!(rootr, rresis, btran,
                                               @Const(mask_patch), @Const(patch_column),
                                               @Const(patch_itype), @Const(rootfr),
                                               @Const(rootfr_unf), @Const(smpso),
                                               @Const(smpsc), @Const(t_soisno),
                                               @Const(watsat), @Const(sucsat),
                                               @Const(bsw), @Const(eff_porosity),
                                               @Const(h2osoi_liqvol),
                                               nlevgrnd::Int, joff::Int,
                                               tfrz, use_unf::Bool,
                                               smk, gkw, gkt)
    p = @index(Global)
    @inbounds if mask_patch[p]
        T = eltype(rootr)
        c = patch_column[p]
        itype = patch_itype[p] + 1  # 0-based Fortran PFT -> 1-based Julia
        acc = zero(T)
        for j in 1:nlevgrnd
            # --- HARD GATE (byte-identical for Float64 via smooth_heaviside) ---
            # Original: zero the layer if it is dry (h2osoi_liqvol<=0) or near-frozen
            # (t_soisno<=tfrz-2). w_gate = σ_w · σ_t is the smooth product of the two
            # one-sided gates; for plain Float64 each σ is an EXACT 0/1 step, so
            # w_gate ∈ {0,1} and the branch below reduces to the original bit-for-bit.
            #
            # The hard condition is STRICT ( >0 / >tfrz-2 ), so the gate must be OFF at
            # the threshold itself. smooth_heaviside uses the `>=` convention (ON at 0),
            # so we form the strict ">" gate as 1-σ(-x): at x=0 → 1-σ(0)=1-1=0 (OFF, like
            # the hard <=0 branch); x>0 → 1-σ(-x)=1 (ON); x<0 → 0 (OFF). Identical 0/1 on
            # Float64, smooth ramp on Dual.
            w_water = one(T) - smooth_heaviside(-(h2osoi_liqvol[c, j + joff]); k = gkw)
            w_temp  = one(T) - smooth_heaviside((tfrz - T(2.0)) - t_soisno[c, j + joff]; k = gkt)
            w_gate  = w_water * w_temp

            # Wetness-branch quantities. Computed unconditionally (the gate-off branch
            # is blended out by w_gate below) — the smooth_max floor on s_node keeps
            # s_node^(-bsw) finite even where the layer is dry/frozen, so the smooth
            # path never evaluates a singular power on a zeroed layer.
            s_node = smooth_max(h2osoi_liqvol[c, j + joff] / eff_porosity[c, j], T(0.01); k = smk)
            s_node = smooth_min(s_node, one(T); k = smk)

            smp_node = soil_suction_clapp_hornberger(sucsat[c, j], s_node, bsw[c, j])
            smp_node = smooth_max(smpsc[itype], smp_node; k = smk)

            rresis_j = min((eff_porosity[c, j] / watsat[c, j]) *
                (smp_node - smpsc[itype]) / (smpso[itype] - smpsc[itype]), one(T))   # HARD: the cap is the CONSTANT 1 (a fully unstressed layer). Its derivative there is zero, so smoothing recovers nothing and only cost a flat log(2)/50 = 1.4% off btran on EVERY well-watered rooted layer -> -1.4% on gs/GPP/transpiration. The s_node and smp_node transitions above stay smooth (BTRAN_SMOOTH_K) — those are the real wilting-point physics.

            rootfr_j = use_unf ? rootfr_unf[p, j] : rootfr[p, j]
            rootr_wet = rootfr_j * rresis_j

            # Blend across the gate. For Float64 the gate is {0,1}, so this is
            # byte-identical to the original branch:
            #   gate=1 → rresis_j / rootr_wet (the original active-layer writes);
            #   gate=0 → rresis[p,j] UNCHANGED (original never wrote it) and
            #            rootr[p,j]=0 (original `rootr[p,j]=zero(T)`).
            rresis[p, j] = w_gate * rresis_j + (one(T) - w_gate) * rresis[p, j]
            rootr[p, j]  = w_gate * rootr_wet

            # btran is the NORMALIZER for rootr (`rootr /= btran` below), so the
            # partition of transpiration over the soil column sums to one — and the
            # water the soil loses equals qflx_tran_veg — IF AND ONLY IF btran is the
            # exact sum of the SAME rootr values that get normalized.
            #
            # Fortran writes `btran += max(rootr(p,j), 0)` (SoilMoistStressMod.F90:413),
            # where the max is a pure guard: rootr = rootfr*rresis >= 0 always
            # (rootfr >= 0; rresis = min(x,1) with x >= 0 because eff_porosity is floored
            # at 0.01 and smp_node is floored at smpsc). So in the EXACT physics the max
            # never bites and btran == sum(rootr) identically — which is what makes
            # transpiration conserve.
            #
            # A SMOOTHED max does bite: smooth_max(r, 0) overshoots max(r, 0) by up to
            # log(2)/k per layer (0.014 at k=50, and rootr itself is O(0.01-0.5)), and it
            # overshoots the ZEROED layers too. That inflated btran WITHOUT changing
            # sum(rootr), so after normalization sum_j rootr[p,j] < 1: the soil gave up
            # less water than the qflx_tran_veg the balance check debits, creating water
            # every step (~9.4e-5 mm/step at Bow).
            #
            # Summing rootr directly keeps storage and flux on the SAME expression, so
            # sum(rootr)/btran == 1 by construction in both paths. Identical to Fortran
            # (and byte-identical for Float64) because rootr >= 0 there.
            acc += rootr[p, j]
        end
        btran[p] = btran[p] + acc
    end
end

function smstress_rootr_btran!(rootr, rresis, btran, mask_patch, patch_column,
                               patch_itype, rootfr, rootfr_unf, smpso, smpsc,
                               t_soisno, watsat, sucsat, bsw, eff_porosity,
                               h2osoi_liqvol, nlevgrnd::Int, joff::Int, tfrz,
                               use_unf::Bool)
    T = eltype(rootr)
    _launch!(_smstress_rootr_btran_kernel!, rootr, rresis, btran, mask_patch,
             patch_column, patch_itype, rootfr, rootfr_unf, smpso, smpsc, t_soisno,
             watsat, sucsat, bsw, eff_porosity, h2osoi_liqvol, nlevgrnd, joff,
             T(tfrz), use_unf, T(BTRAN_SMOOTH_K[]), T(BTRAN_GATE_K_WATER[]),
             T(BTRAN_GATE_K_TEMP[]); ndrange = length(mask_patch))
end

# --- Kernel: zero btran for active patches ---
@kernel function _smstress_zero_btran_kernel!(btran, @Const(mask_patch))
    p = @index(Global)
    @inbounds if mask_patch[p]
        btran[p] = zero(eltype(btran))
    end
end

smstress_zero_btran!(btran, mask_patch) =
    _launch!(_smstress_zero_btran_kernel!, btran, mask_patch; ndrange = length(mask_patch))

"""
    calc_root_moist_stress_clm45default!(...)

Compute root water stress using the default CLM4.5 approach.
Ported from `calc_root_moist_stress_clm45default` in `SoilMoistStressMod.F90`.
"""
function calc_root_moist_stress_clm45default!(rootfr_unf::AbstractMatrix{<:Real},
                                               rootfr::AbstractMatrix{<:Real},
                                               rootr::AbstractMatrix{<:Real},
                                               btran::AbstractVector{<:Real},
                                               rresis::AbstractMatrix{<:Real},
                                               smpso::AbstractVector{<:Real},
                                               smpsc::AbstractVector{<:Real},
                                               t_soisno::AbstractMatrix{<:Real},
                                               watsat::AbstractMatrix{<:Real},
                                               sucsat::AbstractMatrix{<:Real},
                                               bsw::AbstractMatrix{<:Real},
                                               eff_porosity::AbstractMatrix{<:Real},
                                               h2osoi_liqvol::AbstractMatrix{<:Real},
                                               patch_column::AbstractVector{<:Integer},
                                               patch_itype::AbstractVector{<:Integer},
                                               mask_patch::AbstractVector{Bool},
                                               bounds_patch::UnitRange{Int},
                                               nlevgrnd::Int,
                                               nlevsno::Int)
    btran0 = 0.0
    joff = nlevsno
    perchroot = soil_moist_stress_ctrl.perchroot
    perchroot_alt = soil_moist_stress_ctrl.perchroot_alt
    use_unf = perchroot || perchroot_alt

    # Per-patch root resistance + btran accumulation (sequential layer loop in-kernel).
    smstress_rootr_btran!(rootr, rresis, btran, mask_patch, patch_column,
                          patch_itype, rootfr, rootfr_unf, smpso, smpsc, t_soisno,
                          watsat, sucsat, bsw, eff_porosity, h2osoi_liqvol,
                          nlevgrnd, joff, TFRZ, use_unf)

    # Normalize root resistances to get layer contribution to ET.
    # btran[p] is fully accumulated above; here it is read-only, so each
    # (patch, layer) is independent.
    smstress_normalize_rootr!(rootr, mask_patch, btran, btran0, nlevgrnd)

    return nothing
end

# --------------------------------------------------------------------------
# calc_root_moist_stress! (top-level dispatcher)
# --------------------------------------------------------------------------

"""
    calc_root_moist_stress!(...)

Compute root water stress using different approaches.
Ported from `calc_root_moist_stress` in `SoilMoistStressMod.F90`.
"""
function calc_root_moist_stress!(soilstate::SoilStateData,
                                  energyflux::EnergyFluxData,
                                  temperature::TemperatureData,
                                  waterstatebulk::WaterStateBulkData,
                                  waterdiagbulk::WaterDiagnosticBulkData,
                                  col::ColumnData,
                                  patchdata::PatchData,
                                  smpso::AbstractVector{<:Real},
                                  smpsc::AbstractVector{<:Real},
                                  altmax_lastyear_indx::AbstractVector{<:Real},
                                  altmax_indx::AbstractVector{<:Real},
                                  mask_patch::AbstractVector{Bool},
                                  bounds_patch::UnitRange{Int},
                                  nlevgrnd::Int,
                                  nlevsno::Int)

    np = length(bounds_patch)

    # Initialize rootfr_unf to zero. Device-resident scratch (matches the backend
    # of an existing device array) so the whole routine runs on the GPU.
    rootfr_unf = similar(temperature.t_soisno_col, length(mask_patch), nlevgrnd)
    fill!(rootfr_unf, zero(eltype(rootfr_unf)))

    # Initialize btran to zero for accumulation (masked, on-device).
    smstress_zero_btran!(energyflux.btran_patch, mask_patch)

    # Define normalized rootfraction for unfrozen soil
    normalize_unfrozen_rootfr!(soilstate.rootfr_patch,
                                temperature.t_soisno_col,
                                altmax_lastyear_indx,
                                altmax_indx,
                                patchdata.column,
                                rootfr_unf,
                                mask_patch,
                                bounds_patch,
                                nlevgrnd,
                                nlevsno)

    method = soil_moist_stress_ctrl.root_moist_stress_method

    if method == MOIST_STRESS_CLM_DEFAULT
        calc_root_moist_stress_clm45default!(rootfr_unf,
                                              soilstate.rootfr_patch,
                                              soilstate.rootr_patch,
                                              energyflux.btran_patch,
                                              energyflux.rresis_patch,
                                              smpso,
                                              smpsc,
                                              temperature.t_soisno_col,
                                              soilstate.watsat_col,
                                              soilstate.sucsat_col,
                                              soilstate.bsw_col,
                                              soilstate.eff_porosity_col,
                                              waterdiagbulk.h2osoi_liqvol_col,
                                              patchdata.column,
                                              patchdata.itype,
                                              mask_patch,
                                              bounds_patch,
                                              nlevgrnd,
                                              nlevsno)
    else
        error("calc_root_moist_stress: a root moisture stress function must be specified!")
    end

    return nothing
end
