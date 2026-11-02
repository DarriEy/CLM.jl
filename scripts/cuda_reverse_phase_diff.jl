# Per-phase CUDA-vs-CPU reverse comparison for the device phase list (diagnostic for
# scripts/gpu_ad_reverse_driver_validate.jl; run from the repo root:
#   julia --project=scripts scripts/cuda_reverse_phase_diff.jl).
# For each selected phase: same checkpoint, the shadow of every Float array seeded
# (so every output of the phase carries a cotangent and every input receives a
# gradient), one Enzyme reverse call on CPU (device list forced on host) and on CUDA,
# then compare every Float shadow leaf.
#   PHASES=comma,separated,substrings  (default: a cheap diagnostic subset; "all" = every phase)
#   PROFILE=1 / PROFILE_PHASES=a,b  sample the first CUDA call of every / the named
#                    phases into profile_phase<k>.txt (flat) and _tree.txt
#   ENZYME_MAXTYPEDEPTH / ENZYME_MAXTYPEOFFSET  Enzyme type-lattice size for the CUDA
#                    calls (experiment only: 3/64 gives WRONG canopy gradients)
#   SEED=rand|ones   rand (default): independent uniform(0.5,1.5) seed per entry, the
#                    same on both backends, so cancelling adjoints cannot hide behind
#                    a uniform seed; ones: the original all-ones seed
let path = joinpath(@__DIR__, "gpu_ad_reverse_driver_validate.jl"), src = read(path, String)
    # pass the path so the gate's @__DIR__ (input defaults) resolves to scripts/
    include_string(Main, replace(src, "exit(main())" => ""), path)
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

using Random, Profile   # Profile: stdlib, used only with PROFILE=1
function seed_leaves!(db, k)
    seen = Set{UInt}()
    for (p, a) in leaves(db)
        objectid(a) in seen && continue
        push!(seen, objectid(a))
        # keyed on (phase, leaf path) so both backends draw the same seed per leaf
        rng = Random.MersenneTwister(hash((k, p)))
        if get(ENV, "SEED", "rand") == "ones"
            fill!(a, 1.0)
        else
            copyto!(a, 0.5 .+ rand(rng, eltype(a), size(a)))
        end
    end
end

function run_phase(entry, b, revmode, k)
    bk = deepcopy(b); db = Enzyme.make_zero(bk); seed_leaves!(e3(entry)(db), k)
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
    # compile-time experiments: Enzyme type-analysis lattice size (defaults 6 / 512),
    # applied to the CUDA calls only; the CPU reference always compiles with defaults.
    tdepth = parse(Int, get(ENV, "ENZYME_MAXTYPEDEPTH", "6"))
    toffset = parse(Int, get(ENV, "ENZYME_MAXTYPEOFFSET", "512"))
    knobs!(on) = (Enzyme.API.maxtypedepth!(on ? tdepth : 6);
                  Enzyme.API.maxtypeoffset!(on ? toffset : 512))
    CLM._GPU_LAUNCH_PROFILE[] = true
    Base.cumulative_compile_timing(true)
    revmode = Enzyme.set_runtime_activity(Enzyme.Reverse)
    for k in eachindex(pc)
        name = string(pc[k][1])
        if sel == ["all"] || any(s -> occursin(s, name), sel)
            t0 = time()
            knobs!(false)
            dc = run_phase(pc[k], bc, revmode, k); t1 = time()
            jc0 = Base.cumulative_compile_time_ns()[1]
            knobs!(true)
            dg = try
            if get(ENV, "PROFILE", "0") == "1" ||
               any(s -> occursin(s, name), split(get(ENV, "PROFILE_PHASES", ""), ","; keepempty = false))
                # sampling profile of the first (compiling) CUDA call, C frames included
                Profile.init(; n = 10^7, delay = 0.05)
                r = Profile.@profile run_phase(pg[k], bg, revmode, k)
                open("profile_phase$(k).txt", "w") do io
                    Profile.print(IOContext(io, :displaysize => (100000, 400)); format = :flat,
                               sortedby = :count, C = true, mincount = 20)
                end
                open("profile_phase$(k)_tree.txt", "w") do io
                    Profile.print(IOContext(io, :displaysize => (100000, 400)); format = :tree,
                               C = true, mincount = 200, maxdepth = 60, noisefloor = 2.0)
                end
                Profile.clear(); r
            else
                run_phase(pg[k], bg, revmode, k)
            end
            catch err
                @printf("PHASE %2d %-44s CUDA REVERSE ERROR after %.0f s: %s\n", k, name,
                        time() - t1, first(sprint(showerror, err), 600))
                flush(stdout)
                nothing
            end
            if dg === nothing
                pc[k][1](e3(pc[k])(bc), pc[k][2]...); pg[k][1](e3(pg[k])(bg), pg[k][2]...)
                continue
            end
            CUDA.synchronize(); t2 = time()
            jc = (Base.cumulative_compile_time_ns()[1] - jc0) / 1e9  # Julia inference+codegen
            # second CUDA call: compiled, so t2 - t1 minus this is the compile time
            run_phase(pg[k], bg, revmode, k); CUDA.synchronize(); t3 = time()
            lc = leaves(e3(pc[k])(dc)); lg = leaves(e3(pg[k])(dg))
            rows = Tuple{Float64,String,Float64,Float64,Float64}[]; nnz_c = 0; nnz_g = 0
            amax = 0.0; rmax = 0.0
            for ((p, a), (_, g)) in zip(lc, lg)
                ah = Array(a); gh = Array(g)
                nnz_c += count(!iszero, ah); nnz_g += count(!iszero, gh)
                m = isfinite.(ah) .& isfinite.(gh)
                sc = max(maximum(abs, ah[m]; init = 0.0), 1.0)
                ad = maximum(abs.(ah[m] .- gh[m]); init = 0.0)
                d = ad / sc
                # finite-pattern mismatch counts as a full-size difference
                any(isfinite.(ah) .!= isfinite.(gh)) && (d = Inf)
                amax = max(amax, ad)
                rmax = max(rmax, ad / max(maximum(abs, ah[m]; init = 0.0), floatmin()))
                push!(rows, (d, p, ad, maximum(abs, ah[m]; init = 0.0), maximum(abs, gh[m]; init = 0.0)))
            end
            sort!(rows; rev = true)
            @printf("PHASE %2d %-44s nnz cpu=%d gpu=%d maxabs=%.3e maxrel(leaf)=%.3e maxscaled=%.3e  cpu %.0f s, cuda %.0f s (julia-compile %.0f s, run %.1f s)\n",
                    k, name, nnz_c, nnz_g, amax, rmax, rows[1][1], t1 - t0, t2 - t1, jc, t3 - t2)
            for r in first(rows, 8)
                r[1] > 1e-12 && @printf("    scaled=%.3e  %-60s abs=%.3e max|cpu|=%.3e max|gpu|=%.3e\n", r...)
            end
            for (kn, ta, tr, na, nr, ns) in first(CLM.gpu_launch_profile(), 6)
                @printf("      kernel %-48s aug1=%7.1f s rev1=%7.1f s launches=%d/%d sigs=%d\n",
                        kn, ta, tr, na, nr, ns)
            end
            flush(stdout)
        end
        # advance both bundles forward through phase k
        pc[k][1](e3(pc[k])(bc), pc[k][2]...); pg[k][1](e3(pg[k])(bg), pg[k][2]...)
    end
    println("diff done")
end
main_diff()
