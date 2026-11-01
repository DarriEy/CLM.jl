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
- Multi-GPU: `clm_bind_local_gpu!` calls `AMDGPU.device!(AMDGPU.devices()[idx + 1])`
  with CLM's 0-based rank index.
- **Compat:** `AMDGPU = "1"`. AMDGPU.jl 2.x requires KernelAbstractions 0.10,
  and CLM is on KernelAbstractions 0.9, so the resolver installs AMDGPU.jl 1.x
  next to CLM. Moving to AMDGPU 2.x means moving CLM to KA 0.10 first.

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
