using Test
using CLM

@testset "CLM.jl" begin
    include("test_constants.jl")
    include("test_qsat.jl")
    include("test_tridiagonal.jl")
    include("test_decomp.jl")
    include("test_gridcell.jl")
    include("test_landunit.jl")
    include("test_column.jl")
    include("test_patch.jl")
    include("test_filters.jl")
    include("test_temperature.jl")
    include("test_energy_flux.jl")
    include("test_soil_state.jl")
    include("test_soil_hydrology.jl")
    include("test_canopy_state.jl")
    include("test_lake_state.jl")
    include("test_surface_albedo.jl")
    include("test_solar_absorbed.jl")
    include("test_urban_params.jl")
end
