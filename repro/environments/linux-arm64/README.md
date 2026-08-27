# Locked Linux ARM64 sensitivity environment

This container supplies the prespecified Linux compiler-sensitivity cells for the Bow
BTRAN oracle. It is a genuine Linux ARM64 environment, not a Darwin cross-compile. Both
the release (`-O`) and no-optimization (`-O0`) cases must be created from independent
build directories using the same image digest and input bundle.

Build from the repository root:

```sh
docker build -t clm-jl-gmd-linux-arm64:2026-08-27 \
  -f repro/environments/linux-arm64/Dockerfile .
```

The image contains Julia 1.12.6 and a single conda-forge ABI stack: GNU 14.3, MPICH
4.3.2, NetCDF-C 4.9.3, NetCDF-Fortran 4.6.2, PnetCDF 1.14.1, ESMF 8.9.1, CMake 3.31.8,
BLAS/LAPACK 3.11.0, and Python 3.11.16. The completed matrix used image digest
`sha256:fcecc1f80c4423ce2688bbf7033968f1dccbecbe04bfcbb693e1bed17690316c`.
CIME case creation must use `--machine gmd-linux-arm64 --extra-machines-dir` pointing at
the adjacent `machines` directory.

CTSM's bundled PIO2 must be built with the archived
`repro/patches/ctsm5.3.012-linux-build-compat.patch`. It disables parallel netCDF filters
that are absent from this locked ABI stack; it does not modify CLM scientific source.
The `-O0` cell additionally appends `-O0` to `FFLAGS` in its private case macro, after
all CIME defaults, and uses a distinct case and build tree.

The image does not contain CTSM source or scientific inputs. Mount the exact audited CTSM
checkout, CLM.jl checkout, and staged SYMFLUENCE data read-only; use separate writable
case, build, run, Julia depot, and result directories. Do not publish an image containing
input data until the redistribution audit is complete.
