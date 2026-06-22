# DamageMainMod.jl
# Julia port of FATES src/fates/biogeochem/DamageMainMod.F90
#
# The FATES tree-damage module (e.g. wind / biotic structural damage). It
# provides:
#   * IsItDamageTime!      — calendar logic that decides, daily, whether a damage
#     event occurs this time step (sets the module-global `damage_time`).
#   * GetDamageFrac        — fraction of a cohort transitioning between damage
#     classes, looked up from the param_derived transition matrix.
#   * GetCrownReduction    — fraction of crown lost for a given damage class.
#   * GetDamageMortality!  — damage-class-dependent (crown-loss) mortality rate.
#
# Translation notes:
#   * fates_r8 -> Float64; fates_int -> Int.
#   * The Fortran `logical, protected :: damage_time` module variable -> a const
#     `Ref{Bool}` (`damage_time`), mirroring the FatesInterfaceTypesMod hlm_*
#     globals and the FatesDispersalMod dispersal-flag pattern. Access with
#     `damage_time[]`.
#   * `damage_event_code` is already stored as an `Int` in EDParamsMod (the
#     Fortran code does `icode = int(damage_event_code)`), so no extra
#     truncation is needed — read it via `ed_params().damage_event_code`.
#   * Calendar values (hlm_current_day/month/year, hlm_model_day, hlm_day_of_year)
#     are read from the FatesInterfaceTypesMod module-global Refs.
#   * The transition matrix is `param_derived().damage_transitions[cc_cd, nc_cd, pft]`
#     (nlevdamage × nlevdamage × numpft), and the crown bin edges are
#     `ed_params().ED_val_history_damage_bin_edges`.
#   * `damage_mort_p1`/`damage_mort_p2` come from `edpftvarcon_inst()`.
#   * Crown-reduction / damage-fraction / mortality formulas preserved exactly.
#
# Standalone: nothing here is added to CLMInstances or any dual-copied struct.

# =====================================================================================
# Module globals / constants
# =====================================================================================

"""
    damage_time

Module-global flag (Fortran `logical, protected :: damage_time`): `true` when a
damage event occurs during the current time step. Held in a `Ref{Bool}` so it can
be set by [`IsItDamageTime!`](@ref) and read elsewhere via `damage_time[]`.
"""
const damage_time = Ref{Bool}(false)

# Module-level debug flag (Fortran `logical :: debug = .false.`).
const fates_damage_debug = Ref{Bool}(false)

"""
    undamaged_class

Special damage class for undamaged plants (Fortran `undamaged_class = 1`). Used
wherever a cohort's `damageclass` is referenced, flagging that an undamaged plant
is assumed.
"""
const undamaged_class = 1

# =====================================================================================
# IsItDamageTime — daily calendar test for whether damage occurs this step.
# =====================================================================================

"""
    IsItDamageTime!(is_master::Integer)

Determine whether a damage event should occur during the current time step and
store the result in the module-global [`damage_time`](@ref). Called daily. This
is almost an exact replica of the FATES `IsItLoggingTime` logic.

The event is driven by the integer `damage_event_code` (`ed_params()`):
- `1`   — damage turned off.
- `2`   — damage on the first model day (`hlm_model_day == 1`).
- `3`   — damage every day (not recommended).
- `4`   — damage once a month (`hlm_current_day == 1`).
- `< 0 and > -366` — damage every year on day-of-year `abs(icode)`.
- `> 10000` — a specific `YYYYMMDD` event date.
- anything else — invalid; raises (Fortran `endrun`).

If a damage event is enacted and `is_master == itrue`, the event date is logged.
Mirrors Fortran `IsItDamageTime`.
"""
function IsItDamageTime!(is_master::Integer)
    damage_time[] = false
    icode = ed_params().damage_event_code   # already Int (Fortran: int(damage_event_code))

    model_day_int = round(Int, hlm_model_day[])

    if icode == 1
        # Damage is turned off
        damage_time[] = false

    elseif icode == 2
        # Damage event on first time step
        if model_day_int == 1
            damage_time[] = true
        end

    elseif icode == 3
        # Damage event every day - this is not recommended as it will result in a
        # very large number of cohorts which will likely be terminated
        damage_time[] = true

    elseif icode == 4
        # Damage event once a month
        if hlm_current_day[] == 1
            damage_time[] = true
        end

    elseif icode < 0 && icode > -366
        # Damage event every year on a specific day of the year,
        # specified as negative day of year
        if hlm_day_of_year[] == abs(icode)
            damage_time[] = true
        end

    elseif icode > 10000
        # Specific Event: YYYYMMDD
        damage_date  = icode - Int(100 * floor(icode / 100.0))
        damage_year  = Int(floor(icode / 10000.0))
        damage_month = Int(floor(icode / 100.0)) - damage_year * 100

        if hlm_current_day[]   == damage_date  &&
           hlm_current_month[] == damage_month &&
           hlm_current_year[]  == damage_year
            damage_time[] = true
        end

    else
        # Bad damage event flag
        msg = "An invalid damage code was specified in fates_params\n" *
              "Check DamageMainMod.jl:IsItDamageTime!()\n" *
              "for a breakdown of the valid codes and change\n" *
              "fates_damage_event_code in the file accordingly.\nexiting"
        fates_endrun(msg)
    end

    if damage_time[] && (is_master == itrue)
        @info string("Damage Event Enacted on date: ",
                      hlm_current_month[], "-", hlm_current_day[], "-", hlm_current_year[])
    end

    return nothing
end

# =====================================================================================
# GetDamageFrac — transition fraction between damage classes.
# =====================================================================================

"""
    GetDamageFrac(cc_cd::Integer, nc_cd::Integer, pft::Integer) -> Float64

Given the current cohort damage class `cc_cd` and a new damage class `nc_cd` for
a given `pft`, return the fraction of individuals transitioning from `cc_cd` to
`nc_cd`. Consults the `param_derived().damage_transitions` lookup table.
Mirrors Fortran `GetDamageFrac` (returns `dist_frac`).
"""
function GetDamageFrac(cc_cd::Integer, nc_cd::Integer, pft::Integer)
    dist_frac = param_derived().damage_transitions[cc_cd, nc_cd, pft]
    return dist_frac
end

# =====================================================================================
# GetCrownReduction — fraction of crown lost for a damage class.
# =====================================================================================

"""
    GetCrownReduction(crowndamage::Integer) -> Float64

Take the crown damage class `crowndamage` of a cohort and return the fraction of
the crown that is lost: `ED_val_history_damage_bin_edges[crowndamage] / 100`.
Mirrors Fortran `GetCrownReduction` (returns `crown_reduction`).
"""
function GetCrownReduction(crowndamage::Integer)
    crown_reduction = ed_params().ED_val_history_damage_bin_edges[crowndamage] / 100.0
    return crown_reduction
end

# =====================================================================================
# GetDamageMortality — damage-class (crown-loss) dependent mortality.
# =====================================================================================

"""
    GetDamageMortality(crowndamage::Integer, pft::Integer) -> Float64

Compute damage-dependent mortality (the `dgmort` rate) for a cohort with crown
damage class `crowndamage` and plant functional type `pft`. Captures mortality
from mechanisms not otherwise represented in FATES (e.g. pathogens, increased
windthrow vulnerability).

Mortality is a logistic function of crown loss (not of the damage class index,
so it need not be re-parameterised if the number of damage classes changes):

    crown_loss = ED_val_history_damage_bin_edges[crowndamage] / 100
    dgmort     = 1 / (1 + exp(-damage_mort_p2 * (crown_loss - damage_mort_p1)))

except that the undamaged class (`crowndamage == 1`) has `dgmort == 0`.
`damage_mort_p1` (inflection point) and `damage_mort_p2` (rate) are read per-PFT
from `edpftvarcon_inst()`. Mirrors Fortran `GetDamageMortality` (returns
`dgmort`).
"""
function GetDamageMortality(crowndamage::Integer, pft::Integer)
    damage_mort_p1 = edpftvarcon_inst().damage_mort_p1[pft]
    damage_mort_p2 = edpftvarcon_inst().damage_mort_p2[pft]

    # make damage mortality a function of crownloss and not crowndamage class so
    # that it doesn't need to be re-parameterised if the number of damage classes
    # change.
    crown_loss = ed_params().ED_val_history_damage_bin_edges[crowndamage] / 100.0

    if crowndamage == 1
        dgmort = 0.0
    else
        dgmort = 1.0 / (1.0 + exp(-1.0 * damage_mort_p2 * (crown_loss - damage_mort_p1)))
    end

    return dgmort
end
