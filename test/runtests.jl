using Test
using CLM

@testset "CLM.jl" begin
    include("test_constants.jl")
    include("test_qsat.jl")
    include("test_tridiagonal.jl")
    include("test_gridcell.jl")
    include("test_landunit.jl")
    include("test_column.jl")
    include("test_patch.jl")
end
