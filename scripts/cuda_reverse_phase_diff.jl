# Per-phase CUDA-vs-CPU reverse comparison for the device phase list (diagnostic for
# scripts/gpu_ad_reverse_driver_validate.jl; run from the repo root:
#   julia --project=scripts scripts/cuda_reverse_phase_diff.jl).
# For each selected phase: same checkpoint, shadow seeded with ones on every Float
# array, one Enzyme reverse call on CPU (device list forced on host) and on CUDA,
# then compare every Float shadow leaf.
#   PHASES=comma,separated,substrings  (default: a cheap diagnostic subset)
let src = read(joinpath(pwd(), "scripts", "gpu_ad_reverse_driver_validate.jl"), String)
    include_string(Main, replace(src, "exit(main())" => ""))
end

function leaves!(out, x, path, seen)
    if x isa AbstractArray{<:AbstractFloat}
        push!(out, path => x)
    elseif x isa AbstractArray || x isa Number || x isa Symbol || x isa String ||
           x isa Function || x === nothing
        return out
    else
        objectid(x) in seen && return out
        ismutable(x) && push!(seen, objectid(x))
        for f in fieldnames(typeof(x))
            isdefined(x, f) || continue
            leaves!(out, getfield(x, f), string(path, ".", f), seen)
        end
    end
    out
end
leaves(x) = leaves!(Pair{String,Any}[], x, "", Set{UInt}())

function seed_ones!(db)
    seen = Set{UInt}()
    for (_, a) in leaves(db)
        objectid(a) in seen && continue
        push!(seen, objectid(a)); fill!(a, 1.0)
    end
end

function run_phase(entry, b, revmode)
    bk = deepcopy(b); db = Enzyme.make_zero(bk); seed_ones!(e3(entry)(db))
    Enzyme.autodiff(revmode, entry[1], Enzyme.Const,
        Enzyme.Duplicated(e3(entry)(bk), e3(entry)(db)), map(Enzyme.Const, entry[2])...)
    db
end
e3(entry) = length(entry) == 3 ? entry[3] : identity

function main_diff()
    sel = split(get(ENV, "PHASES", "tk_interface,lwrad,heatdiff,infil_rev,setqflx"), ",")
    inst64, bounds, filt, config = warmed_state()
    inst_gpu = CLM.Adapt.adapt(CUDA.CuArray, CLM.Adapt.adapt(_CUDAHost(), deepcopy(inst64)))
    filt_gpu = CLM.Adapt.adapt(CUDA.CuArray, CLM.Adapt.adapt(_CUDAHost(), deepcopy(filt)))
    kw = (; n_canopy = 2, dtime = 1800.0, include_surface = true)
    bc = CLM.driver_rev_bundle(deepcopy(inst64)); bg = CLM.driver_rev_bundle(inst_gpu)
    pc = CLM.driver_rev_phases(bounds, filt, config; canopy_aux = CLM.canopy_rev_aux(bc.inst, bounds, filt),
                               inst = bc.inst, device = true, kw...)
    pg = CLM.driver_rev_phases(bounds, filt_gpu, config; canopy_aux = CLM.canopy_rev_aux(inst_gpu, bounds, filt_gpu),
                               inst = inst_gpu, kw...)
    Enzyme.API.strictAliasing!(false)
    revmode = Enzyme.set_runtime_activity(Enzyme.Reverse)
    for k in eachindex(pc)
        name = string(pc[k][1])
        if any(s -> occursin(s, name), sel)
            t0 = time()
            dc = run_phase(pc[k], bc, revmode); dg = run_phase(pg[k], bg, revmode)
            CUDA.synchronize()
            lc = leaves(e3(pc[k])(dc)); lg = leaves(e3(pg[k])(dg))
            rows = Tuple{Float64,String,Float64,Float64}[]; nnz_c = 0; nnz_g = 0
            for ((p, a), (_, g)) in zip(lc, lg)
                ah = Array(a); gh = Array(g)
                nnz_c += count(!iszero, ah); nnz_g += count(!iszero, gh)
                m = isfinite.(ah) .& isfinite.(gh)
                sc = max(maximum(abs, ah[m]; init = 0.0), 1.0)
                d = maximum(abs.(ah[m] .- gh[m]); init = 0.0) / sc
                push!(rows, (d, p, maximum(abs, ah[m]; init = 0.0), maximum(abs, gh[m]; init = 0.0)))
            end
            sort!(rows; rev = true)
            @printf("PHASE %2d %-44s nnz cpu=%d gpu=%d  (%.0f s)\n", k, name, nnz_c, nnz_g, time() - t0)
            for r in first(rows, 8)
                r[1] > 1e-12 && @printf("    rel=%.3e  %-60s max|cpu|=%.3e max|gpu|=%.3e\n", r...)
            end
            flush(stdout)
        end
        # advance both bundles forward through phase k
        pc[k][1](e3(pc[k])(bc), pc[k][2]...); pg[k][1](e3(pg[k])(bg), pg[k][2]...)
    end
    println("diff done")
end
main_diff()
