# =============================================================================
# REVERSE-AD ON GPU (Metal) PROBE — (c) of a→b→c. The deferred frontier: the forward CLM kernels
# run on Metal (Phase A/B), and the reverse chain works on CPU (Enzyme); this probes whether
# Enzyme can reverse-differentiate Metal (Apple-Silicon) GPU code. Enzyme's GPU reverse support is
# primarily CUDA — Metal reverse is unproven — so this is an honest escalating probe (each test
# guarded), reporting exactly what works and where it breaks.
#   julia +1.12 --project=/tmp/gpu_rev_env scripts/enzyme_metal_reverse_probe.jl
# =============================================================================
using Metal, Enzyme, KernelAbstractions, Printf
const KA = KernelAbstractions

report(name, fn) = try
    r = fn(); @printf("  [PASS] %-46s → %s\n", name, string(r)); true
catch e
    msg = sprint(showerror, e); @printf("  [FAIL] %-46s → %s\n", name, first(split(msg, '\n')))
    false
end

println("="^80); println("REVERSE-AD ON METAL probe (Apple M-series)"); println("="^80)

# --- Test 0: Metal functions at all (forward) ---
println("\n[0] Metal forward sanity:")
report("MtlArray sum(x.^2), expect 30.0", () -> begin
    x = MtlArray(Float32[1,2,3,4]); Float64(sum(x.^2))
end)

# --- Test 1: Enzyme REVERSE over a GPU broadcast (no explicit kernel) ---
println("\n[1] Enzyme.Reverse over a Metal broadcast  f(x)=sum(x.^2),  expect dx=2x=[2,4,6,8]:")
report("autodiff Reverse, Duplicated MtlArray", () -> begin
    f(x) = sum(x .^ 2)
    x  = MtlArray(Float32[1,2,3,4]); dx = MtlArray(zeros(Float32,4))
    Enzyme.autodiff(Enzyme.Reverse, f, Enzyme.Active, Enzyme.Duplicated(x, dx))
    Array(dx)
end)

# --- Test 2: Enzyme REVERSE over a KernelAbstractions kernel launched on Metal ---
@kernel function sq_kernel!(out, @Const(inp))
    i = @index(Global)
    @inbounds out[i] = inp[i]^2
end
function run_sq!(out, inp, backend)
    sq_kernel!(backend, 64)(out, inp, ndrange=length(inp))
    KA.synchronize(backend)
    return nothing
end
println("\n[2] Enzyme.Reverse over a KA kernel on Metal  out=inp²,  seed dout=1 ⇒ dinp=2·inp=[2,4,6,8]:")
report("autodiff Reverse, MetalBackend KA kernel", () -> begin
    backend = MetalBackend()
    inp  = MtlArray(Float32[1,2,3,4]); out  = MtlArray(zeros(Float32,4))
    dinp = MtlArray(zeros(Float32,4)); dout = MtlArray(ones(Float32,4))
    Enzyme.autodiff(Enzyme.Reverse, run_sq!, Enzyme.Const,
        Enzyme.Duplicated(out, dout), Enzyme.Duplicated(inp, dinp), Enzyme.Const(backend))
    Array(dinp)
end)

# --- Test 3: CPU-KA control (same kernel on CPU backend) — proves the kernel/Enzyme path itself ---
println("\n[3] CONTROL: same KA kernel on the CPU backend (isolates Metal-specific failures):")
report("autodiff Reverse, CPU() KA kernel", () -> begin
    backend = CPU()
    inp  = Float32[1,2,3,4]; out  = zeros(Float32,4)
    dinp = zeros(Float32,4); dout = ones(Float32,4)
    Enzyme.autodiff(Enzyme.Reverse, run_sq!, Enzyme.Const,
        Enzyme.Duplicated(out, dout), Enzyme.Duplicated(inp, dinp), Enzyme.Const(backend))
    dinp
end)

println("\n", "="^80)
println("Summary: [0] = Metal alive; [1] = Enzyme reverse over Metal broadcast; [2] = Enzyme reverse")
println("over a Metal KA kernel (the real target); [3] = CPU-KA control. A [2] PASS = reverse-AD on")
println("Metal works; [2] FAIL + [3] PASS = the kernel/Enzyme path is fine, the gap is Metal-specific.")
println("="^80)
