# build_patched_cesm.jl — relink the existing CTSM cesm.exe with patched Fortran sources.
#
# Recompiles each given source with the exact command recorded in the case's lnd
# build log, replaces its object in a PRIVATE copy of libclm.a, and relinks a
# private cesm.exe with the recorded link line. The case's own build tree is only
# read, never written, so concurrent users of that case are unaffected.
#
#   julia build_patched_cesm.jl <outdir> <patched.F90>...
#
# Each patched file replaces the source of the same basename in the original build.
# Environment:
#   CTSM_CASE_BLD  bld dir of the case that built $SYMFLUENCE_DATA/installs/clm/bin/cesm.exe
#                  (default: $SYMFLUENCE_DATA/installs/clm/cases/bow_btran_macos/bld)
#   CTSM_TOOLCHAIN toolchain bin dir with mpif90 (default: .../installs/clm/toolchain/bin)
#
# With no patched files this reproduces the reference cesm.exe's results bit-for-bit
# (checked on the fates_parity_1pt planthydro case).
D = get(ENV, "SYMFLUENCE_DATA", "/Users/darri/repos/SYMFLUENCE_data")
B = get(ENV, "CTSM_CASE_BLD", "$D/installs/clm/cases/bow_btran_macos/bld")
ENV["PATH"] = get(ENV, "CTSM_TOOLCHAIN", "$D/installs/clm/toolchain/bin") * ":" * ENV["PATH"]
out = abspath(ARGS[1]); srcs = abspath.(ARGS[2:end])
mkpath("$out/build"); mkpath("$out/cpl")

gz(f) = read(pipeline(`gzcat $f`), String)
lndlog = gz(only(filter(f -> startswith(f, "lnd.bldlog"), readdir(B; join=false))) |> f -> "$B/$f")
cesmlog = gz(only(filter(f -> startswith(f, "cesm.bldlog"), readdir(B; join=false))) |> f -> "$B/$f")
libdir = "$B/gnu/mpich/nodebug/nothreads/nuopc/nuopc/esmf/lib/"
cp("$libdir/libclm.a", "$out/build/libclm.a"; force=true)
for f in readdir("$B/cpl/obj"); endswith(f, ".o") && cp("$B/cpl/obj/$f", "$out/cpl/$f"; force=true); end

for src in srcs
    base = basename(src)
    line = only(filter(l -> startswith(l, "mpif90 -c") && endswith(strip(l), "/" * base),
                       split(lndlog, '\n')))
    tok = String.(split(strip(line)))
    tok[end] = src                         # compile the patched file instead
    cd(() -> run(Cmd(tok)), "$out/build")  # .o/.mod land in the private build dir
    obj = replace(base, r"\.F90$" => ".o")
    run(`ar r $out/build/libclm.a $out/build/$obj`)
end

link = String.(split(strip(only(filter(l -> startswith(l, "mpif90 -o ../../cesm.exe"),
                                         split(cesmlog, '\n'))))))
link[findfirst(==("../../cesm.exe"), link)] = "$out/cesm.exe"
link[findfirst(==("-L$libdir"), link)] = "-L$out/build/"
cd(() -> run(Cmd(link)), "$out/cpl")
println("built $out/cesm.exe with ", isempty(srcs) ? "no patches" : join(basename.(srcs), ", "))
