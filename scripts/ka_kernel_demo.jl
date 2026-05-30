# First KernelAbstractions kernel over the batched (ncopies=N) column layout.
#
# Kernelizes a real per-column physics op — specific humidity from vapor pressure
# (clm_driver.jl:530) — indexed by column with a per-column gridcell lookup, the
# exact data-parallel pattern the GPU port needs. Runs on KA's CPU backend here;
# the identical kernel runs on CUDA/AMD by passing a GPU backend + device arrays
# (this Apple-Silicon machine is Metal/Float32-only, so GPU execution is elsewhere).
#
#   julia --project=. scripts/ka_kernel_demo.jl [ncopies]

using CLM
using KernelAbstractions
using Printf

const FS = get(ENV, "CLM_FSURDAT",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc")
const PF = get(ENV, "CLM_PARAMFILE",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc")

# Scalar reference (the current driver loop).
function forc_q_scalar!(forc_q_col, gridcell, forc_vp_grc, forc_pbot_col)
    for c in eachindex(forc_q_col)
        g = gridcell[c]
        vp = forc_vp_grc[g]
        pbot = forc_pbot_col[c]
        forc_q_col[c] = 0.622 * vp / max(pbot - 0.378 * vp, 1.0)
    end
end

# KernelAbstractions kernel — one thread per column. Backend-agnostic: the same
# code compiles for CPU(), CUDABackend(), ROCBackend(), etc.
@kernel function forc_q_kernel!(forc_q_col, @Const(gridcell), @Const(forc_vp_grc), @Const(forc_pbot_col))
    c = @index(Global)
    g = gridcell[c]
    vp = forc_vp_grc[g]
    pbot = forc_pbot_col[c]
    forc_q_col[c] = 0.622 * vp / max(pbot - 0.378 * vp, 1.0)
end

ncopies = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 500
@printf("Initializing ncopies=%d ...\n", ncopies)
(inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat=FS, paramfile=PF, ncopies=ncopies)
nc = bounds.endc
ng = bounds.endg

# Forcing inputs (per-gridcell vp, per-column pbot).
for g in 1:ng
    inst.atm2lnd.forc_vp_grc[g] = 800.0
end
for c in 1:nc
    inst.atm2lnd.forc_pbot_downscaled_col[c] = 85000.0
end

gridcell = inst.column.gridcell
vp_grc = inst.atm2lnd.forc_vp_grc
pbot_col = inst.atm2lnd.forc_pbot_downscaled_col

# Reference
q_ref = zeros(nc)
forc_q_scalar!(q_ref, gridcell, vp_grc, pbot_col)

# KA kernel on the CPU backend
backend = CPU()
q_ka = zeros(nc)
kern = forc_q_kernel!(backend)
kern(q_ka, gridcell, vp_grc, pbot_col; ndrange=nc)
KernelAbstractions.synchronize(backend)

err = maximum(abs.(q_ka .- q_ref))
@printf("columns = %d  (ng=%d)\n", nc, ng)
@printf("max |KA - scalar| = %.3e\n", err)

# quick timing (warm)
kern(q_ka, gridcell, vp_grc, pbot_col; ndrange=nc); KernelAbstractions.synchronize(backend)
t_ka = @elapsed (for _ in 1:100; kern(q_ka, gridcell, vp_grc, pbot_col; ndrange=nc); end; KernelAbstractions.synchronize(backend))
t_sc = @elapsed (for _ in 1:100; forc_q_scalar!(q_ref, gridcell, vp_grc, pbot_col); end)
@printf("100x over %d cols: KA(CPU backend)=%.3f ms  scalar=%.3f ms\n", nc, 1e3 * t_ka / 100, 1e3 * t_sc / 100)

println(err < 1e-12 ? "\nKA KERNEL DEMO PASSED ✓ (matches scalar; same kernel runs on GPU backends)" : "\nFAILED ✗")
