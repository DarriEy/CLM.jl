using Test
using CLM

@testset "CLM.jl" begin
    include("test_constants.jl")
    include("test_qsat.jl")
    include("test_tridiagonal.jl")
    include("test_gridcell.jl")
end
