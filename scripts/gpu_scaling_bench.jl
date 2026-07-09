# ==========================================================================
# gpu_scaling_bench.jl — CPU-vs-Metal scaling of a whole clm_drv! biogeophys
# step as the column count N grows. Single-column is GPU-hostile (all launch +
# marshaling overhead, no parallelism); this finds the crossover where the GPU
# overtakes the CPU and the asymptotic speedup, to extrapolate continental scale.
#
#   julia +1.12 --project=scripts scripts/gpu_scaling_bench.jl
#
# Timings are min-of-trials (min approximates the uncontended time, so this stays
# meaningful even while other Metal work runs in the background).
# ==========================================================================
using CLM, Printf
import Metal
include(joinpath(@__DIR__, "gpu_backends.jl"))
include(joinpath(@__DIR__, "gpu_adapt.jl"))
include(joinpath(@__DIR__, "clmdrv_make_data.jl"))
mf(x) = mf(Metal.MtlArray, x)

const DRV = (true, 1.0, 0.0, 0.0, 0.4091, false, false, "20260101", false)
const CFG = CLM.CLMDriverConfig(use_cn=false)   # biogeophys (SP) step — the parity path's bulk
step!(inst, filt, fia, bounds, ps) = CLM.clm_drv!(CFG, inst, filt, fia, bounds, DRV...;
    nstep=1, is_first_step=false, is_beg_curr_day=false, is_beg_curr_year=false, photosyns=ps)

# min-of-`k` elapsed seconds of `f` (f must fully finish, incl. GPU sync)
bestof(f, k) = minimum(_ -> (@elapsed f()), 1:k)

function bench_N(N; ktrial=5)
    # scale the fixture to N columns (1 patch/col, single gridcell/landunit column-set)
    inst, bounds, filt, fia, _c, ps = make_driver_data_physical(ng=1, nl=3, nc=N, np=N)
    inst.canopystate.frac_veg_nosno_alb_patch .= 1
    # ---- CPU ----
    step!(inst, filt, fia, bounds, ps)                       # warmup / compile
    tcpu = bestof(() -> step!(inst, filt, fia, bounds, ps), ktrial)
    # ---- Metal ----
    inst_d = mf(inst); ps_d = mf(ps); filt_d = mf(filt); fia_d = mf(fia)
    inst_d.temperature.t_soisno_col isa Metal.MtlArray || return (N, tcpu, NaN, NaN)
    step!(inst_d, filt_d, fia_d, bounds, ps_d); Metal.synchronize()   # warmup / compile
    tgpu = bestof(() -> (step!(inst_d, filt_d, fia_d, bounds, ps_d); Metal.synchronize()), ktrial)
    return (N, tcpu, tgpu, tcpu / tgpu)
end

function main()
    println("="^68)
    println("  clm_drv! (biogeophys) CPU-vs-Metal scaling — one step, min-of-trials")
    println("="^68)
    @printf("  %10s %12s %12s %10s\n", "N cols", "CPU (ms)", "GPU (ms)", "speedup")
    println("  " * "-"^46)
    for N in (1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144)
        try
            n, tc, tg, sp = bench_N(N)
            if isnan(tg)
                @printf("  %10d %12.3f %12s %10s  (no GPU)\n", n, tc*1e3, "-", "-")
            else
                @printf("  %10d %12.3f %12.3f %9.2f×%s\n", n, tc*1e3, tg*1e3, sp,
                        sp >= 1 ? "  <-- GPU wins" : "")
            end
            flush(stdout)
        catch e
            @printf("  %10d   FAILED: %s\n", N, first(split(sprint(showerror, e), "\n")))
            flush(stdout); (N >= 4096) && break
        end
    end
end
main()
