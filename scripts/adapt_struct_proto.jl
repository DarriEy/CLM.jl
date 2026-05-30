# Prototype: make a CLM-style @kwdef mutable struct device-movable via Adapt,
# WITHOUT breaking the existing `T{FT}()` and `T(; kwargs)` construction patterns
# (used in cold_start and the Dual make_dual_copy path).
#
# Pattern: parametrize over FT + array types (M, V), default them to CPU Array in
# a convenience constructor, and register Adapt.@adapt_structure. Adapt then
# rebuilds the struct with the device array types (e.g. CuArray) automatically.

using Adapt

# Current CLM style is: mutable struct TemperatureData{FT<:Real}
#   t_soisno_col::Matrix{FT} = Matrix{FT}(undef,0,0); t_grnd_col::Vector{FT} = FT[]
# New device-movable style:
Base.@kwdef mutable struct TemperatureData{FT<:Real,
                                           M<:AbstractMatrix{FT},
                                           V<:AbstractVector{FT}}
    t_soisno_col::M = Matrix{FT}(undef, 0, 0)
    t_grnd_col::V = FT[]
    nlevsno::Int = 5            # non-array field (stays as-is)
end

# Convenience constructor preserving the existing call site `TemperatureData{Float64}()`:
# fills CPU Array defaults so M, V are inferred.
function TemperatureData{FT}(; kwargs...) where {FT<:Real}
    return TemperatureData{FT, Matrix{FT}, Vector{FT}}(; kwargs...)
end

# Adapt: rebuild with adapted (device) arrays. Non-array fields pass through.
Adapt.@adapt_structure TemperatureData

# ---------------- validation ----------------
# 1. Existing construction patterns still work
t0 = TemperatureData{Float64}()
println("T{Float64}() ok: ", typeof(t0).name.name, " eltype=", eltype(t0.t_grnd_col))

t1 = TemperatureData{Float64}(t_soisno_col = rand(3, 4), t_grnd_col = rand(3))
println("kwargs construction ok: size t_soisno=", size(t1.t_soisno_col))

# 2. Adapt round-trip. Use a custom AbstractArray wrapper as a stand-in "device"
#    array type to prove the struct genuinely holds a different array type
#    (real GPU would be CuArray; this machine is Metal/Float32-only so we mock).
struct FakeDevArray{T,N} <: AbstractArray{T,N}
    data::Array{T,N}
end
Base.size(a::FakeDevArray) = size(a.data)
Base.getindex(a::FakeDevArray, i...) = getindex(a.data, i...)
Base.setindex!(a::FakeDevArray, v, i...) = setindex!(a.data, v, i...)
struct ToFakeDev end
Adapt.adapt_storage(::ToFakeDev, x::Array) = FakeDevArray(x)

t1_dev = adapt(ToFakeDev(), t1)
println("after adapt: t_soisno type = ", typeof(t1_dev.t_soisno_col))
println("after adapt: struct type params = ", typeof(t1_dev))

ok = (t0 isa TemperatureData) &&
     (t1.t_soisno_col isa Matrix) &&
     (t1_dev.t_soisno_col isa FakeDevArray) &&
     (t1_dev.t_grnd_col isa FakeDevArray) &&
     (t1_dev.nlevsno == t1.nlevsno)
println(ok ? "\nADAPT STRUCT PATTERN PASSED ✓" : "\nFAILED ✗")
