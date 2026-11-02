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
#   ENZYME_MAXTYPEDEPTH / ENZYME_MAXTYPEOFFSET  override CLM.GPU_REVERSE_TYPE_LATTICE
#                    (the CUDA calls' Enzyme type lattice); CLM_GPU_REVERSE_TYPE_LATTICE=default
#                    uses Enzyme's defaults
#   DETAIL=1         also: CUDA rerun determinism, and for the worst leaves a central
#                    FD (CPU forward) of the seeded loss at the 3 worst entries
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
    # the library's per-phase call (device type lattice on CUDA, Enzyme defaults on CPU)
    CLM._phase_autodiff!(revmode, entry[1], e3(entry)(bk), e3(entry)(db), entry[2])
    db
end
e3(entry) = length(entry) == 3 ? entry[3] : identity

# The scalar the seeded reverse differentiates: L = Σ_leaves Σ_i seed_i · leaf_i over
# the phase OUTPUT (deduplicated by array object, like the seed). Central FD of it
# w.r.t. one input entry, on the CPU device-list phase from checkpoint bc, summed as
# Σ seed·(out₊ − out₋) over entries finite in both runs (no large-sum cancellation).
function phase_fd(entry, bc, k, path, idx)
    function run(δ)
        bk = deepcopy(bc); v = e3(entry)(bk)
        a = Dict(leaves(v))[path]; a[idx] += δ
        entry[1](v, entry[2]...)
        v
    end
    x0 = Dict(leaves(e3(entry)(bc)))[path][idx]
    h = 1e-6 * max(abs(x0), 1.0)
    lp = leaves(run(h)); lm = Dict(leaves(run(-h)))
    seen = Set{UInt}(); s = 0.0
    for (q, a) in lp
        objectid(a) in seen && continue
        push!(seen, objectid(a))
        rng = Random.MersenneTwister(hash((k, q)))
        w = get(ENV, "SEED", "rand") == "ones" ? ones(eltype(a), size(a)) :
            0.5 .+ rand(rng, eltype(a), size(a))
        d = Array(a) .- Array(lm[q]); m = isfinite.(d)
        s += sum(w[m] .* d[m])
    end
    s / (2h), x0, h
end

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
    # Enzyme type-analysis lattice for the CUDA calls (CLM.GPU_REVERSE_TYPE_LATTICE,
    # default (3, 64); the CPU reference always compiles with Enzyme's 6 / 512).
    if haskey(ENV, "ENZYME_MAXTYPEDEPTH") || haskey(ENV, "ENZYME_MAXTYPEOFFSET")
        CLM.GPU_REVERSE_TYPE_LATTICE[] = (parse(Int, get(ENV, "ENZYME_MAXTYPEDEPTH", "6")),
                                          parse(Int, get(ENV, "ENZYME_MAXTYPEOFFSET", "512")))
    end
    println("CUDA type lattice: ", get(ENV, "CLM_GPU_REVERSE_TYPE_LATTICE", "") == "default" ?
            "Enzyme defaults" : CLM.GPU_REVERSE_TYPE_LATTICE[])
    CLM._GPU_LAUNCH_PROFILE[] = true
    Base.cumulative_compile_timing(true)
    revmode = Enzyme.set_runtime_activity(Enzyme.Reverse)
    for k in eachindex(pc)
        name = string(pc[k][1])
        if sel == ["all"] || any(s -> occursin(s, name), sel)
            t0 = time()
            dc = run_phase(pc[k], bc, revmode, k); t1 = time()
            jc0 = Base.cumulative_compile_time_ns()[1]
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
            dg2 = run_phase(pg[k], bg, revmode, k); CUDA.synchronize(); t3 = time()
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
            if get(ENV, "DETAIL", "0") == "1"
                # determinism of the CUDA reverse (same input, same seed, twice)
                nd = 0.0
                for ((_, a), (_, b2)) in zip(lg, leaves(e3(pg[k])(dg2)))
                    ah = Array(a); bh = Array(b2); m = isfinite.(ah) .& isfinite.(bh)
                    nd = max(nd, maximum(abs.(ah[m] .- bh[m]); init = 0.0))
                end
                @printf("    CUDA rerun max|Δ| = %.3e\n", nd)
                # arbitrate the worst leaves with a central FD of the seeded loss
                dcl = Dict(lc); dgl = Dict(lg)
                for r in first(rows, 4)
                    r[1] > 1e-8 || continue
                    p = r[2]; ah = Array(dcl[p]); gh = Array(dgl[p])
                    m = isfinite.(ah) .& isfinite.(gh)
                    dd = ifelse.(m, abs.(ah .- gh), 0.0)
                    for idx in first(sortperm(vec(dd); rev = true), 3)
                        ci = CartesianIndices(dd)[idx]
                        fd, x0, h = phase_fd(pc[k], bc, k, p, ci)
                        @printf("      %s%s x=%.6g: cpu=% .10e cuda=% .10e FD(h=%.1e)=% .10e  -> %s\n",
                                p, string(Tuple(ci)), x0, ah[ci], gh[ci], h, fd,
                                abs(fd - ah[ci]) < abs(fd - gh[ci]) ? "CPU closer" : "CUDA closer")
                    end
                end
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
