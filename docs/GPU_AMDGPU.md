# AMD GPUs (ROCm) — the `:amdgpu` backend

CLM.jl's physics runs as KernelAbstractions kernels, so an AMD GPU is one more
backend behind the same registry as CPU, CUDA and Metal
(`src/infrastructure/backend.jl`). AMDGPU.jl is a **weak dependency**: `using CLM`
never loads it, and the default test suite does not need it.

> **Status: not yet run on AMD hardware.** Everything below the "Verified" line
> was checked on machines without an AMD GPU. Treat the `:amdgpu` path as
> untested until `scripts/gpu_validate_amdgpu.jl` passes on a ROCm box.

## Usage

```julia
using CLM
import AMDGPU                     # loads ext/CLMAMDGPUExt.jl, registers :amdgpu

CLM.clm_set_backend(:amdgpu; require_functional = true)   # errors if no ROCm device
inst_d = CLM.clm_device_array(inst)    # Adapt.adapt(ROCArray, inst) — Float64
CLM.clm_ka_backend()                   # AMDGPU.ROCBackend()
# ... run clm_drv! on inst_d exactly as on CPU/CUDA ...
CLM.clm_set_backend(:cpu)
```

- Precision is **Float64** (`clm_float_type(backend = :amdgpu) == Float64`), the
  same as CPU and CUDA. No Float32 down-adaptor is needed (that is Metal-only).
- BitVector masks must become `Vector{Bool}` before moving to the device, as on
  CUDA (see the `_HostPrep` adaptor in `scripts/gpu_adapt.jl`).
- Multi-GPU: `clm_bind_local_gpu!` passes CLM's 0-based rank index to the bind
  hook, which calls `AMDGPU.device!(AMDGPU.devices()[CLM._device_list_slot(idx, n)])`
  (slot = `idx + 1`; an out-of-range index raises `ArgumentError`, like
  `CUDA.device!`). The Metal extension uses the same helper.
- **Compat:** `AMDGPU = "1, 2"`. Only AMDGPU.jl **2.9 and later** need
  KernelAbstractions 0.10. AMDGPU.jl 2.0–2.8 still use KA `0.9.2 - 0.9`, so they
  install next to CLM's KA 0.9.45. With CUDA.jl 5 and Metal.jl in the same
  environment, the resolver picks **AMDGPU.jl 2.7.0** (2.8 needs GPUCompiler 2,
  which CUDA.jl 5 does not allow). AMDGPU.jl 2.9+ waits on the KA 0.10 move;
  see "KernelAbstractions 0.10" below.

## What the extension does

`ext/CLMAMDGPUExt.jl` mirrors `ext/CLMCUDAExt.jl`:

1. Registers `:amdgpu` with `adapt_to = adapt(ROCArray, ·)`,
   `ka_backend = ROCBackend()`, `functional = AMDGPU.functional()`, and the
   device-count / device-bind hooks.
2. Installs the device-only overlay `CLM._smooth_f64() = false`. Without it,
   Float64 `smooth_*` kernels read a host `Ref` from the device, which on CUDA
   showed up as `ERROR_ILLEGAL_ADDRESS`. The macro is `AMDGPU.Device.@device_override`.
   An earlier version checked for `@device_override` on the top-level `AMDGPU`
   module, where it is not defined, so the overlay was skipped without any
   warning. That is now fixed and checked by a test (`ext._SMOOTH_OVERRIDE_INSTALLED`).

No shared code path is CUDA-specific. The audit (Oct 2026) found no `CuArray`
type checks, no `CUDA.@allowscalar` (the driver uses `GPUArraysCore.@allowscalar`),
and no Enzyme rules tied to a backend: the `scatter_add_1d!` reverse rule launches
a KA kernel on whatever backend its arrays use. The 64-bit-atomics check in
`src/infrastructure/kernels.jl` excludes only Metal, so ROCm uses `Int64`
counters, as CUDA does.

## Verified (no AMD GPU)

Run on 2026-10-30 on the NVIDIA/Windows 11 server (RTX 5070 Ti, no ROCm, Julia
1.12.6). AMDGPU.jl was added to a scratch environment next to a checkout of this
branch. The resolver installed **AMDGPU.jl 1.3.6** with KernelAbstractions 0.9.45.

- `CLM → CLMAMDGPUExt` precompiles. `Base.get_extension(CLM, :CLMAMDGPUExt)` is loaded.
  `AMDGPU.functional() == false`. AMDGPU prints its normal "ROCm discovery failed /
  library unavailable" messages on a machine without ROCm.
- `isdefined(AMDGPU, Symbol("@device_override")) == false`, which confirms the old
  guard could never succeed. `ext._SMOOTH_OVERRIDE_INSTALLED == true` with the fix.
- `test/test_amdgpu_backend.jl`: the registration tier passes 9/9. Its checks:
  `:amdgpu` is registered, `clm_ka_backend(backend = :amdgpu) isa ROCBackend`,
  `clm_backend_functional(:amdgpu) == false`, the Float64 precision policy,
  `clm_set_backend(:amdgpu)` succeeds, and `require_functional = true` raises an
  error. The device tier is skipped, as expected without a device.
- The device tier's comparison code (`_finite_reldiff`, the Bool-mask host prep,
  and the `clm_drv!` call) was run once with the CPU standing in for the device.
  The result was bit-identical, as it must be, because both sides ran on the CPU.
  That run showed the 1-column fixture leaves the energy fields NaN, so the
  device tier uses the 4-column fixture instead. This run checks only the test
  code and is **not** evidence of AMD parity.

Without AMDGPU.jl installed (the default), the test file shows two `@test_skip`
and nothing else changes.

Re-checked on 2026-11-01 on an Apple Silicon Mac (no ROCm, Julia 1.12.6) with
`AMDGPU = "1, 2"`. The scratch environment held CLM (dev), AMDGPU, CUDA and Metal.
The resolver installed AMDGPU.jl 2.7.0, CUDA.jl 5.11.3, Metal.jl 1.9.3,
KernelAbstractions 0.9.45, Enzyme 0.13.216, GPUCompiler 1.15.2 and LLVM.jl 9.13.2.
`CLM → CLMAMDGPUExt` precompiles and the registration tier passes 9/9. The device
tier is skipped. `AMDGPU.devices`, `AMDGPU.device!(::HIPDevice)`, `ROCBackend`
and `AMDGPU.Device.@device_override` all still exist in 2.7.

## KernelAbstractions 0.10 — not installable with Enzyme yet (checked 2026-11-01)

Widening the compat to `KernelAbstractions = "0.9.45, 0.10"` and asking Pkg for
0.10 fails to resolve:

```
ERROR: Unsatisfiable requirements detected for package LLVM [929cbde3]:
 ├─restricted by compatibility requirements with KernelAbstractions [63c18a36] to versions: 10.0.0 - 10.4.2
 │   └─restricted to versions 0.10 by an explicit requirement, leaving only versions: 0.10.0 - 0.10.4
 └─restricted by compatibility requirements with Enzyme [7da242da] to versions: 9.8.1 - 9.13.2 — no versions left
   └─restricted to versions 0.13.152 - 0.13 by CLM [4245f2e0], leaving only versions: 0.13.152 - 0.13.216
```

- Every KA 0.10.x requires LLVM.jl 10. Every released Enzyme.jl, up to 0.13.216
  (and Enzyme.jl `main`, `LLVM = "9.10"`), requires LLVM.jl 9. Enzyme is a hard
  dependency of CLM, so no environment can hold both. The Enzyme upgrade is open as
  [EnzymeAD/Enzyme.jl#3741](https://github.com/EnzymeAD/Enzyme.jl/pull/3741)
  ("Upgrade to LLVM.jl 10").
- KA 0.10 also **removed its Enzyme extension**
  ([JuliaGPU/KernelAbstractions.jl#824](https://github.com/JuliaGPU/KernelAbstractions.jl/pull/824),
  "Temporarily remove Enzyme support"). Restoring it is
  [#826](https://github.com/JuliaGPU/KernelAbstractions.jl/issues/826), which is still
  open (milestone 0.10.1, but 0.10.4 shipped without it). CLM's reverse AD depends on
  that extension. The `_gpu_launch!` rule in `src/infrastructure/enzyme_rules.jl`
  forwards to KA's own kernel reverse rule, and Enzyme on CPU differentiates KA
  kernel launches through the same extension. Even once LLVM.jl resolves, CLM
  cannot move until #826 lands.
- Other KA 0.10 breaking changes that CLM would have to review then (from the
  v0.10.0 release notes): the `CPU` backend is now an OpenCL backend running on
  PoCL, so every CPU kernel and the CPU Enzyme path go through a new compiler.
  `CPU(; static=true)`, `isgpu`, `@ka_code_llvm` and `@ka_code_typed` are removed.
  `GPU` is deprecated (CLM's `_kcall!` dispatches on `backend isa KA.GPU`).
  `@private` arrays become a `PrivateArray`. `@Const` arguments are constified at
  launch. Backends must implement KernelInterface 0.4.
- KA 0.10 also means moving every GPU package at once, because each pins one KA
  minor. KA 0.10 needs CUDA.jl 6.5+ (CUDACore 6.5, so compat `CUDA = "5, 6"`),
  Metal.jl 1.12+, and AMDGPU.jl 2.9+.

Not tested, because the environment cannot be built: the CPU kernel tests and
CPU/CUDA/Metal reverse AD under KA 0.10. Re-check when Enzyme.jl releases with
LLVM.jl 10 support and KA has merged #826.

## Untested (needs a ROCm device)

- Any kernel launch on `ROCBackend`, including the one-step device-vs-CPU
  `clm_drv!` test (tier 3 of `test/test_amdgpu_backend.jl`).
- The kernel, forward-AD and whole-driver parity suites. Run them on a ROCm box with:
  ```bash
  julia --project=scripts -e 'using Pkg; Pkg.develop(path="."); Pkg.add("AMDGPU"); Pkg.instantiate()'
  julia --project=scripts scripts/gpu_validate_amdgpu.jl     # exit 0 pass / 1 fail / 2 no device
  ```
  The script sets `CLM_GPU_BACKEND=AMDGPU`, which makes `scripts/gpu_backends.jl`
  use AMDGPU even if another GPU is present.
- Reverse-mode Enzyme on ROCm, and multi-GPU binding on real AMD devices.
