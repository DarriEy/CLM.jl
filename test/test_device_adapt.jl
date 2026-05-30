# ==========================================================================
# test_device_adapt.jl — state structs and the CLMInstances tree are
# device-movable via Adapt.
#
# Uses a stand-in AbstractArray ("device") type since the CI/dev machine has no
# Float64 GPU; the same mechanism moves arrays to CuArray/ROCArray on real GPUs.
# ==========================================================================

using Test
using CLM
using Adapt

# Minimal stand-in "device" array type.
struct FakeDev{T,N} <: AbstractArray{T,N}
    d::Array{T,N}
end
Base.size(a::FakeDev) = size(a.d)
Base.getindex(a::FakeDev, i...) = getindex(a.d, i...)
Base.setindex!(a::FakeDev, v, i...) = setindex!(a.d, v, i...)
Adapt.adapt_storage(::Type{<:FakeDev}, x::Array) = FakeDev(x)

@testset "Device-movable structs (Adapt)" begin
    # Leaf state structs round-trip to the device array type.
    t = CLM.TemperatureData{Float64}()
    t.t_grnd_col = rand(4)
    t.t_soisno_col = rand(4, 3)
    td = adapt(FakeDev, t)
    @test td.t_grnd_col isa FakeDev
    @test td.t_soisno_col isa FakeDev
    @test td isa CLM.TemperatureData              # still a TemperatureData (different params)
    @test td.t_grnd_col == t.t_grnd_col           # values preserved

    c = CLM.ColumnData{Float64}()
    c.dz = rand(2, 5)
    @test adapt(FakeDev, c).dz isa FakeDev

    # Whole-instance tree: registered sub-structs move, the rest pass through.
    inst = CLM.CLMInstances()
    inst.temperature.t_grnd_col = rand(3)
    inst.column.dz = rand(3, 4)
    inst.energyflux.eflx_lh_tot_patch = rand(3)
    inst.atm2lnd.forc_t_not_downscaled_grc = rand(3)
    idev = adapt(FakeDev, inst)
    @test idev.temperature.t_grnd_col isa FakeDev
    @test idev.column.dz isa FakeDev
    @test idev.energyflux.eflx_lh_tot_patch isa FakeDev
    @test idev.atm2lnd.forc_t_not_downscaled_grc isa FakeDev   # forcing moves too
    # WaterData is intentionally not yet device-movable (bulk_and_tracers
    # reference-aliasing) — it passes through unchanged.
    @test typeof(idev.water) === typeof(inst.water)
    println("  device-adapt: CLMInstances tree moves registered sub-structs to device")
end
