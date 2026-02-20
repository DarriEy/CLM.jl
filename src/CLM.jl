module CLM

using LinearAlgebra

# ===========================================================================
# Tier 1: Constants & Parameters
# ===========================================================================
include("constants/varctl.jl")
include("constants/landunit_varcon.jl")
include("constants/column_varcon.jl")
include("constants/varpar.jl")
include("constants/varcon.jl")

# ===========================================================================
# Tier 1: Infrastructure (solvers, utilities)
# ===========================================================================
include("infrastructure/tridiagonal.jl")
include("infrastructure/band_diagonal.jl")

# ===========================================================================
# Tier 1: Biogeophysics (pure math)
# ===========================================================================
include("biogeophys/qsat.jl")

end # module CLM
