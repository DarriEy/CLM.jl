# ==========================================================================
# Ported from: src/main/surfrdUtilsMod.F90
# Weight normalization and collapse utilities for surface data reading
# ==========================================================================

"""
    check_sums_equal_1!(arr, name; tol=SUM_TO_1_TOL, sumto=nothing)

Confirm that sum(arr[n,:]) ≈ target for all rows n. `arr` is a Matrix.
If `sumto` is provided it should be a Vector of target sums per row (default 1.0).
Throws an error if any row's sum deviates by more than `tol`.
"""
function check_sums_equal_1!(arr::AbstractMatrix{Float64}, name::String;
                              tol::Float64=SUM_TO_1_TOL,
                              sumto::Union{Nothing,AbstractVector{Float64}}=nothing)
    nrows = size(arr, 1)
    for n in 1:nrows
        target = sumto === nothing ? 1.0 : sumto[n]
        s = sum(@view arr[n, :])
        if abs(s - target) > tol
            error("check_sums_equal_1!: sum of $name not $target at row $n, sum=$s")
        end
    end
    nothing
end

"""
    renormalize!(arr; target=1.0)

Renormalize each row of `arr` so that it sums to `target`.
Rows that sum to zero are left unchanged.
"""
function renormalize!(arr::AbstractMatrix{Float64}; target::Float64=1.0)
    for n in 1:size(arr, 1)
        s = sum(@view arr[n, :])
        if s != 0.0
            ratio = target / s
            arr[n, :] .*= ratio
        end
    end
    nothing
end

"""
    apply_convert_ocean_to_land!(wt_lunit)

Convert ocean points to land by moving ocean weight to natveg (ISTSOIL).
"""
function apply_convert_ocean_to_land!(wt_lunit::AbstractMatrix{Float64})
    for g in 1:size(wt_lunit, 1)
        wt_lunit[g, ISTSOIL] += wt_lunit[g, ISTOCN]
        wt_lunit[g, ISTOCN] = 0.0
    end
    nothing
end

"""
    collapse_crop_types!(wt_cft, fert_cft, begg, endg; sumto=nothing)

Collapse unused crop types into types used in this run.
- If !irrigate, merge irrigated CFTs into rainfed CFTs.
- Merge CFTs via pftcon.mergetoclmpft mapping.
"""
function collapse_crop_types!(wt_cft::AbstractMatrix{Float64},
                               fert_cft::AbstractMatrix{Float64},
                               begg::Int, endg::Int;
                               sumto::Union{Nothing,Vector{Float64}}=nothing)
    cftsize = size(wt_cft, 2)
    cftsize > 0 || return nothing

    total_sum = sumto === nothing ? ones(endg - begg + 1) : sumto

    # If not using irrigation, merge irrigated into rainfed
    if !varctl.irrigate
        for g in begg:endg
            # Stride-2 merge: rainfed at odd offsets, irrigated at even
            for m in nc3crop:2:(varpar.maxveg - 1)
                irr = m + 1
                if irr <= varpar.maxveg
                    gi = g - begg + 1
                    m_idx = m - varpar.cft_lb + 1
                    irr_idx = irr - varpar.cft_lb + 1
                    if m_idx >= 1 && m_idx <= cftsize && irr_idx >= 1 && irr_idx <= cftsize
                        wt_cft[gi, m_idx] += wt_cft[gi, irr_idx]
                        wt_cft[gi, irr_idx] = 0.0
                    end
                end
            end
        end
    end

    # Merge CFTs via mergetoclmpft mapping
    for g in begg:endg
        gi = g - begg + 1
        for m in 1:varpar.maxveg
            merge_to = pftcon.mergetoclmpft[m + 1]  # +1 for Julia 1-based
            if m != merge_to
                m_idx = m - varpar.cft_lb + 1
                to_idx = merge_to - varpar.cft_lb + 1
                if m_idx >= 1 && m_idx <= cftsize && to_idx >= 1 && to_idx <= cftsize
                    wt_to = wt_cft[gi, to_idx]
                    wt_from = wt_cft[gi, m_idx]
                    wt_merge = wt_to + wt_from
                    wt_cft[gi, to_idx] = wt_merge
                    wt_cft[gi, m_idx] = 0.0
                    if wt_merge > 0.0
                        fert_cft[gi, to_idx] = (wt_to * fert_cft[gi, to_idx] +
                                                  wt_from * fert_cft[gi, m_idx]) / wt_merge
                    end
                    pftcon.is_pft_known_to_model[m + 1] = false
                end
            end
        end
    end

    nothing
end

"""
    collapse_individual_lunits!(wt_lunit; toosmall_soil, toosmall_crop,
        toosmall_glacier, toosmall_lake, toosmall_wetland, toosmall_urban)

Remove landunits below user-defined percentage thresholds.
If all landunits get removed for a gridcell, keep the largest.
"""
function collapse_individual_lunits!(wt_lunit::AbstractMatrix{Float64};
                                      toosmall_soil::Float64=-1.0,
                                      toosmall_crop::Float64=-1.0,
                                      toosmall_glacier::Float64=-1.0,
                                      toosmall_lake::Float64=-1.0,
                                      toosmall_wetland::Float64=-1.0,
                                      toosmall_urban::Float64=-1.0)
    toosmall_any = toosmall_soil + toosmall_crop + toosmall_glacier +
                   toosmall_lake + toosmall_wetland + toosmall_urban
    toosmall_any > 0.0 || return nothing

    # Build threshold array (fraction, not percent)
    toosmall = zeros(MAX_LUNIT)
    toosmall[ISTSOIL] = toosmall_soil / 100.0
    toosmall[ISTCROP] = toosmall_crop / 100.0
    toosmall[ISTICE] = toosmall_glacier / 100.0
    toosmall[ISTDLAK] = toosmall_lake / 100.0
    toosmall[ISTWET] = toosmall_wetland / 100.0
    toosmall[ISTURB_TBD] = toosmall_urban / 100.0
    toosmall[ISTURB_HD] = toosmall_urban / 100.0
    toosmall[ISTURB_MD] = toosmall_urban / 100.0

    ng = size(wt_lunit, 1)
    residual = zeros(MAX_LUNIT)

    for g in 1:ng
        residual .= 0.0
        for m in 1:MAX_LUNIT
            if wt_lunit[g, m] > 0.0 && wt_lunit[g, m] <= toosmall[m]
                residual[m] = wt_lunit[g, m]
                wt_lunit[g, m] = 0.0
            end
        end
        # If all removed, keep largest
        if sum(@view wt_lunit[g, :]) == 0.0
            max_lu = argmax(residual)
            wt_lunit[g, max_lu] = residual[max_lu]
        end
    end

    renormalize!(wt_lunit; target=1.0)
    nothing
end

"""
    collapse_to_dominant!(weight, n_dominant)

Collapse to top N dominant PFTs or landunits per gridcell.
If `n_dominant <= 0` or `n_dominant >= size(weight,2)`, does nothing.
"""
function collapse_to_dominant!(weight::AbstractMatrix{Float64}, n_dominant::Int)
    ncols = size(weight, 2)
    (n_dominant > 0 && n_dominant < ncols) || return nothing

    ng = size(weight, 1)
    for g in 1:ng
        wt_sum = sum(@view weight[g, :])

        # Find top-N indices
        row = collect(@view weight[g, :])
        top_indices = partialsortperm(row, 1:n_dominant; rev=true)

        wt_dom_sum = sum(weight[g, idx] for idx in top_indices)

        # Normalize dominants and zero others
        if wt_sum > 0.0 && wt_dom_sum > 0.0
            for idx in top_indices
                weight[g, idx] = weight[g, idx] * wt_sum / wt_dom_sum
            end
        end
        for m in 1:ncols
            if !(m in top_indices)
                weight[g, m] = 0.0
            end
        end
    end
    nothing
end
