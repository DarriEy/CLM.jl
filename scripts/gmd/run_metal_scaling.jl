# ==========================================================================
# run_metal_scaling.jl — E9 accelerator performance scaling (protocol form).
#
# CPU-vs-Metal wall time of one whole clm_drv! biogeophys (SP) step over a
# logarithmic column-count sweep, reported per the experiment protocol §E9:
# median and IQR over R repetitions after a discarded warm-up/compile step,
# compilation excluded, host-device transfers excluded (state is resident on
# each side; the step is timed with a full device synchronize), Float64 host
# vs Float32 device (the production Metal precision), thread count and device
# recorded. Emits JSON with full provenance; the crossover N is reported, not
# only the maximum speedup.
#
#   julia --project=scripts scripts/gmd/run_metal_scaling.jl [out.json]
#
# Unlike scripts/gpu_scaling_bench.jl (min-of-trials, exploratory), this is
# the qualification runner: run it on an otherwise idle machine.
# ==========================================================================
using CLM, Printf, Statistics
include(joinpath(@__DIR__, "..", "gpu_backends.jl"))
include(joinpath(@__DIR__, "..", "gpu_adapt.jl"))
include(joinpath(@__DIR__, "..", "clmdrv_make_data.jl"))
mf(x) = mf(device_array_type(), x)

const DRV = (true, 1.0, 0.0, 0.0, 0.4091, false, false, "20260101", false)
const CFG = CLM.CLMDriverConfig(use_cn=false)   # biogeophys (SP) step
step!(inst, filt, fia, bounds, ps) = CLM.clm_drv!(CFG, inst, filt, fia, bounds, DRV...;
    nstep=1, is_first_step=false, is_beg_curr_day=false, is_beg_curr_year=false, photosyns=ps)

"R timed repetitions of `f` (each fully synchronized), AFTER one discarded warmup."
function timings(f, R)
    f()                                   # warm-up / compile — discarded
    return [(@elapsed f()) for _ in 1:R]
end

med(x) = median(x)
iqr(x) = quantile(x, 0.75) - quantile(x, 0.25)

function bench_N(N; R=11)
    inst, bounds, filt, fia, _c, ps = make_driver_data_physical(ng=1, nl=3, nc=N, np=N)
    inst.canopystate.frac_veg_nosno_alb_patch .= 1
    tc = timings(() -> step!(inst, filt, fia, bounds, ps), R)
    inst_d = mf(inst); ps_d = mf(ps); filt_d = mf(filt); fia_d = mf(fia)
    inst_d.temperature.t_soisno_col isa device_array_type() ||
        return (; N, cpu_s = tc, gpu_s = Float64[])
    tg = timings(() -> (step!(inst_d, filt_d, fia_d, bounds, ps_d); device_synchronize()), R)
    return (; N, cpu_s = tc, gpu_s = tg)
end

function main()
    out = length(ARGS) >= 1 ? ARGS[1] : "metal_scaling.json"
    R = parse(Int, get(ENV, "SCALING_REPS", "11"))
    rows = []
    println("  N cols    CPU med (ms)  ±IQR      GPU med (ms)  ±IQR      speedup")
    for N in (1, 4, 16, 64, 256, 1024, 4096, 16384, 65536)
        r = bench_N(N; R=R)
        isempty(r.gpu_s) && (println("  $(r.N): device adapt failed"); continue)
        cm, ci, gm, gi = med(r.cpu_s), iqr(r.cpu_s), med(r.gpu_s), iqr(r.gpu_s)
        @printf("  %7d  %10.2f  %8.2f  %10.2f  %8.2f  %8.2fx\n",
                r.N, 1e3cm, 1e3ci, 1e3gm, 1e3gi, cm / gm)
        push!(rows, (; r.N, cpu_median_s = cm, cpu_iqr_s = ci,
                       gpu_median_s = gm, gpu_iqr_s = gi,
                       speedup_median = cm / gm, reps = R,
                       cpu_all_s = r.cpu_s, gpu_all_s = r.gpu_s))
    end
    crossover = nothing
    for row in rows
        if row.speedup_median > 1.0
            crossover = row.N; break
        end
    end
    dev = try string(Metal.device()) catch; "unknown" end
    record = (;
        kind = "metal-scaling-e9",
        julia = string(VERSION),
        threads = Threads.nthreads(),
        device = dev,
        host_precision = "Float64",
        device_precision = "Float32",
        transfers_included = false,
        statistic = "median over reps after 1 discarded warmup; IQR dispersion",
        workload = "one clm_drv! biogeophys (SP) step, synthetic physical fixture, 1 patch/col",
        crossover_N = crossover,
        rows,
    )
    open(out, "w") do io
        # hand-rolled JSON to avoid adding a dependency
        function j(io, x)
            if x isa NamedTuple
                print(io, "{"); first_ = true
                for (k, v) in pairs(x)
                    first_ || print(io, ","); first_ = false
                    print(io, "\"", k, "\":"); j(io, v)
                end
                print(io, "}")
            elseif x isa AbstractVector
                print(io, "["); for (i, v) in enumerate(x)
                    i > 1 && print(io, ","); j(io, v)
                end; print(io, "]")
            elseif x isa AbstractString || x isa Symbol
                print(io, "\"", x, "\"")
            elseif x === nothing
                print(io, "null")
            else
                print(io, x)
            end
        end
        j(io, record); println(io)
    end
    println("wrote $out (crossover N = $crossover)")
end

main()
