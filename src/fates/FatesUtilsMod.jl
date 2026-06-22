# FatesUtilsMod.jl
# Julia port of FATES src/fates/main/FatesUtilsMod.F90
#
# General helper functions: HLM-list checking, real-value sanity checks,
# great-circle distance, array search, and quadratic-root solvers.
# Deps: FatesConstantsMod (nearzero, earth_radius_eq, rad_per_deg),
#       FatesGlobals (fates_log, fates_endrun).

"""
    check_hlm_list(hlms, hlm_name) -> Bool

Compare a string of HLM tags to see if any match the active HLM name. Mirrors
the Fortran `scan(trim(hlms), trim(hlm_name)) > 0` semantics (a character-set
scan: true if any character of `hlm_name` appears in `hlms`).
"""
function check_hlm_list(hlms::AbstractString, hlm_name::AbstractString)
    h = strip(hlms)
    name = strip(hlm_name)
    # Fortran SCAN returns the position of the first char of `name` found in `h`.
    for ch in name
        if occursin(ch, h)
            return true
        end
    end
    return false
end

# =====================================================================================

"""
    check_var_real(r8_var, var_name) -> return_code

Sanity-check a real value. Returns 0 if clean; adds 1 (NaN), 10 (near overflow),
100 (near underflow) to the return code, and logs each issue.
"""
function check_var_real(r8_var::Real, var_name::AbstractString)
    overflow = floatmax(Float64)
    underflow = floatmin(Float64)

    return_code = 0

    # NaN check
    if r8_var != r8_var
        @warn "NaN detected, $(strip(var_name)): $r8_var"
        return_code = 1
    end

    # Overflow check (within 100th of max precision)
    if abs(r8_var) > 0.01 * overflow
        @warn "Nigh overflow detected, $(strip(var_name)): $r8_var"
        return_code += 10
    end

    # Underflow check (within 100x of min precision)
    if abs(r8_var) < 100.0 * underflow
        @warn "Nigh underflow detected, $(strip(var_name)): $r8_var"
        return_code += 100
    end

    return return_code
end

# =====================================================================================

"""
    GreatCircleDist(slons, slonf, slats, slatf) -> distance [m]

Great-circle distance between two points (source `s`, forepoint `f`), given in
degrees lon/lat. Accurate for small and large distances.
"""
function GreatCircleDist(slons::Real, slonf::Real, slats::Real, slatf::Real)
    # Convert co-ordinates to radians.
    lons = slons * rad_per_deg
    lonf = slonf * rad_per_deg
    lats = slats * rad_per_deg
    latf = slatf * rad_per_deg
    dlon = lonf - lons

    # Find the arcs.
    x = sin(lats) * sin(latf) + cos(lats) * cos(latf) * cos(dlon)
    y = sqrt((cos(latf) * sin(dlon)) * (cos(latf) * sin(dlon)) +
             (cos(lats) * sin(latf) - sin(lats) * cos(latf) * cos(dlon)) *
             (cos(lats) * sin(latf) - sin(lats) * cos(latf) * cos(dlon)))

    # Convert the arcs to actual distance.
    return earth_radius_eq * atan(y, x)
end

# =====================================================================================

"""
    GetNeighborDistance(gi, gj, latc, lonc) -> great-circle distance [m]

Distance between gridcells `gi` and `gj` given lat/lon vectors.
"""
function GetNeighborDistance(gi::Integer, gj::Integer,
                             latc::AbstractVector, lonc::AbstractVector)
    return GreatCircleDist(lonc[gi], lonc[gj], latc[gi], latc[gj])
end

# =====================================================================================

"""
    FindIndex(input_string_array, string_to_match) -> index

Standin for the intrinsic FINDLOC. Returns the (1-based) index of the last
matching string, or 0 if no match is found.
"""
function FindIndex(input_string_array::AbstractVector{<:AbstractString},
                   string_to_match::AbstractString)
    array_index = 0
    for i in 1:length(input_string_array)
        if strip(input_string_array[i]) == strip(string_to_match)
            array_index = i
        end
    end
    return array_index
end

# =====================================================================================

"""
    QuadraticRootsNSWC(a, b, c) -> (root1, root2)

Real roots of a x^2 + b x + c, based on the NSWC Mathematics Subroutine Library
(overflow-avoiding discriminant). Aborts if only imaginary roots are generated.
"""
function QuadraticRootsNSWC(a::Real, b::Real, c::Real)
    if abs(a) < nearzero
        root2 = 0.0
        if b != 0.0
            root2 = -c / b
        end
        root1 = 0.0
    elseif abs(c) > nearzero
        # compute discriminant avoiding overflow
        b1 = b / 2.0
        if abs(b1) < abs(c)
            e = a
            if c < 0.0
                e = -a
            end
            e = b1 * (b1 / abs(c)) - e
            d = sqrt(abs(e)) * sqrt(abs(c))
        else
            e = 1.0 - (a / b1) * (c / b1)
            d = sqrt(abs(e)) * abs(b1)
        end
        if e < 0.0
            # complex conjugate zeros
            @warn "error, imaginary roots detected in quadratic solve"
            fates_endrun("imaginary roots detected in QuadraticRootsNSWC")
        else
            # real zeros
            if b1 >= 0.0
                d = -d
            end
            root1 = (-b1 + d) / a
            root2 = 0.0
            if root1 != 0.0
                root2 = (c / root1) / a
            end
        end
    else
        root2 = 0.0
        root1 = -b / a
    end

    return root1, root2
end

# =====================================================================================

"""
    QuadraticRootsSridharachary(a, b, c) -> (root1, root2)

Real roots of a x^2 + b x + c via the Sridharacharya (quadratic) formula.
Aborts if only imaginary roots are generated.
"""
function QuadraticRootsSridharachary(a::Real, b::Real, c::Real)
    # If a is 0, equation is linear, not quadratic.
    if abs(a) < nearzero
        root2 = 0.0
        if abs(b) > nearzero
            root2 = -c / b
        end
        root1 = 0.0
        return root1, root2
    end

    d = b * b - 4.0 * a * c
    das = sqrt(abs(d))

    if d > nearzero
        root1 = (-b + das) / (2.0 * a)
        root2 = (-b - das) / (2.0 * a)
    elseif abs(d) <= nearzero
        root1 = -b / (2.0 * a)
        root2 = root1
    else
        @warn "error, imaginary roots detected in quadratic solve"
        fates_endrun("imaginary roots detected in QuadraticRootsSridharachary")
        root1 = NaN
        root2 = NaN
    end

    return root1, root2
end
