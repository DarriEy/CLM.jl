# FatesFuelClassesMod.jl
# Julia port of FATES src/fates/fire/FatesFuelClassesMod.F90
#
# The fuel-class enumeration used by SPITFIRE. There are six fuel classes:
#   1) twigs, 2) small branches, 3) large branches, 4) trunks,
#   5) dead leaves, 6) live grass
# The Fortran defines a derived type `fuel_classes_type` with *private* integer
# fields holding the pool array indices, plus public accessor type-bound
# procedures (twigs/small_branches/large_branches/trunks/dead_leaves/live_grass)
# and a single module-level instance `fuel_classes` passed around.
#
# Julia port: the array indices become module-level `const`s (the canonical
# source of truth), and we reproduce the Fortran accessor pattern with a
# singleton `fuel_classes_type` and accessor methods that read those indices, so
# downstream code can call `twigs(fuel_classes)` exactly as the Fortran does.
# `fates_r8` is unused here. Deps: FatesLitterMod (ncwd, for the count check).

# ---------------------------------------------------------------------------
# Number of fuel classes
# ---------------------------------------------------------------------------
const num_fuel_classes = 6  # number of total fuel classes

# ---------------------------------------------------------------------------
# Fuel-class array indices (Fortran private fields of fuel_classes_type)
# ---------------------------------------------------------------------------
const twigs_i = 1           # array index for twigs pool
const small_branches_i = 2  # array index for small branches pool
const large_branches_i = 3  # array index for large branches pool
const trunks_i = 4          # array index for trunks pool
const dead_leaves_i = 5     # array index for dead leaves pool
const live_grass_i = 6      # array index for live grass pool

"""
    fuel_classes_type

Singleton accessor type for the SPITFIRE fuel-class indices (Fortran
`fuel_classes_type`). Carries no data; the indices live as module `const`s. Use
the module instance [`fuel_classes`](@ref) together with the accessor functions
([`twigs`](@ref) etc.) to retrieve a pool index, mirroring the Fortran
`fuel_classes%twigs()` calls.
"""
struct fuel_classes_type end

"""
    fuel_classes

The actual fuel-class accessor instance we can pass around (Fortran
`type(fuel_classes_type), public :: fuel_classes`).
"""
const fuel_classes = fuel_classes_type()

"""
    twigs(this::fuel_classes_type) -> Int

Array index for the twigs pool.
"""
twigs(this::fuel_classes_type) = twigs_i

"""
    small_branches(this::fuel_classes_type) -> Int

Array index for the small branches pool.
"""
small_branches(this::fuel_classes_type) = small_branches_i

"""
    large_branches(this::fuel_classes_type) -> Int

Array index for the large branches pool.
"""
large_branches(this::fuel_classes_type) = large_branches_i

"""
    trunks(this::fuel_classes_type) -> Int

Array index for the trunks pool.
"""
trunks(this::fuel_classes_type) = trunks_i

"""
    dead_leaves(this::fuel_classes_type) -> Int

Array index for the dead leaves pool.
"""
dead_leaves(this::fuel_classes_type) = dead_leaves_i

"""
    live_grass(this::fuel_classes_type) -> Int

Array index for the live grass pool.
"""
live_grass(this::fuel_classes_type) = live_grass_i
