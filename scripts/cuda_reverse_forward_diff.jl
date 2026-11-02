# Forward-side companion of scripts/cuda_reverse_phase_diff.jl (diagnostic):
#  1. aliasing report: Float arrays shared between two paths of the driver bundle on the
#     host but separate after adapt to CUDA (or vice versa) — a shadow seeded per path
#     and a write through one alias then behave differently on the two backends;
#  2. forward state after each device phase, CPU (device list forced on host) vs CUDA,
#     every Float leaf, max |Δ| over entries finite on both, plus finite-pattern
#     mismatches.
# Run from the repo root: julia --project=scripts scripts/cuda_reverse_forward_diff.jl
let path = joinpath(@__DIR__, "gpu_ad_reverse_driver_validate.jl"), src = read(path, String)
    include_string(Main, replace(src, "exit(main())" => ""), path)
end

function fleaves!(out, x, path, seen)
    if x isa AbstractArray{<:AbstractFloat}
        push!(out, path => x)
    elseif x isa AbstractArray || x isa Number || x isa Symbol || x isa String ||
           x isa Function || x === nothing
    else
        objectid(x) in seen && return out
        ismutable(x) && push!(seen, objectid(x))
        for f in fieldnames(typeof(x))
            isdefined(x, f) || continue
            fleaves!(out, getfield(x, f), string(path, ".", f), seen)
        end
    end
    out
end
fleaves(x) = fleaves!(Pair{String,Any}[], x, "", Set{UInt}())

# path => canonical path (first path holding the same array object)
function alias_classes(x)
    first_of = Dict{UInt,String}(); canon = Dict{String,String}()
    for (p, a) in fleaves(x)
        id = objectid(a)
        canon[p] = get!(first_of, id, p)
    end
    canon
end

function main_fwd()
    inst64, bounds, filt, config = warmed_state()
    inst_gpu = CLM.Adapt.adapt(CUDA.CuArray, CLM.Adapt.adapt(_CUDAHost(), deepcopy(inst64)))
    filt_gpu = CLM.Adapt.adapt(CUDA.CuArray, CLM.Adapt.adapt(_CUDAHost(), deepcopy(filt)))
    bc = CLM.driver_rev_bundle(deepcopy(inst64)); bg = CLM.driver_rev_bundle(inst_gpu)
    ac = alias_classes(bc); ag = alias_classes(bg)
    nmis = 0
    for p in sort(collect(keys(ac)))
        haskey(ag, p) || (println("  path only on CPU: ", p); continue)
        ac[p] == ag[p] && continue
        nmis += 1
        nmis <= 40 && println("  alias mismatch: ", p, "  cpu->", ac[p], "  gpu->", ag[p])
    end
    println("alias mismatches: ", nmis, " of ", length(ac), " Float leaves")
    kw = (; n_canopy = 2, dtime = 1800.0, include_surface = true)
    pc = CLM.driver_rev_phases(bounds, filt, config; canopy_aux = CLM.canopy_rev_aux(bc.inst, bounds, filt),
                               inst = bc.inst, device = true, kw...)
    pg = CLM.driver_rev_phases(bounds, filt_gpu, config; canopy_aux = CLM.canopy_rev_aux(inst_gpu, bounds, filt_gpu),
                               inst = inst_gpu, kw...)
    e3(entry) = length(entry) == 3 ? entry[3] : identity
    for k in eachindex(pc)
        pc[k][1](e3(pc[k])(bc), pc[k][2]...); pg[k][1](e3(pg[k])(bg), pg[k][2]...)
        CUDA.synchronize()
        lg = Dict(fleaves(bg))
        rows = Tuple{Float64,Float64,Int,String}[]
        for (p, a) in fleaves(bc)
            haskey(lg, p) || continue
            ah = Array(a); gh = Array(lg[p]); size(ah) == size(gh) || continue
            m = isfinite.(ah) .& isfinite.(gh)
            nf = count(isfinite.(ah) .!= isfinite.(gh))
            d = maximum(abs.(ah[m] .- gh[m]); init = 0.0)
            sc = max(maximum(abs, ah[m]; init = 0.0), 1e-300)
            (d > 0 || nf > 0) && push!(rows, (d / sc, d, nf, p))
        end
        sort!(rows; rev = true)
        @printf("FWD %2d %-44s leaves differing: %d\n", k, string(pc[k][1]), length(rows))
        for r in first(rows, 8)
            @printf("    rel=%.3e abs=%.3e nonfinite-mismatch=%d  %s\n", r...)
        end
        flush(stdout)
    end
    println("forward diff done")
end
main_fwd()
