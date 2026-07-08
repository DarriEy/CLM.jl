# ==========================================================================
# gpu_validate_matrixcn_e2e.jl — REAL-Metal parity for the matrix-CN sparse
# operator layer (sparse_matrix_multiply.jl). The matrix-CN veg + soil solvers
# are built entirely on these per-unit KA kernels; this runs each one on the CPU
# (Float64 golden) and on the GPU (Float32) over the SAME operands and compares.
# It is the "golden CPU vs real Metal" proof the CPU-only test_matrixcn_gpu_ready
# suite cannot give (that runs the KA CPU backend, not the device).
#
# Covers every kernelized op:
#   set_value_dm! / set_value_v! / set_value_v_scaler!   (scalar->eltype fills)
#   spmm_ak!  (A <- A*K)      spmm_ax!  (X <- X + A*X, the pool advance)
#   spmp_b_acc! (B <- B + A)  set_value_a! memoized value-fill (device M/RI/CI)
# plus the veg-solver GLUE kernels that replaced the host scalar-index loops in
# the solve body: _veg_aoned_kernel! / _veg_binput_kernel! / _veg_add_vec_kernel! /
# _veg_load_c_kernel! / _veg_writeback_c_kernel!.
# RI/CI/filter are carried to the device via Adapt (item: index arrays on device).
#
#   julia --project=scripts scripts/gpu_validate_matrixcn_e2e.jl
# ==========================================================================
using CLM
using Printf
import Metal
include(joinpath(@__DIR__, "gpu_backends.jl"))

# Adaptor that moves every array field of a sparse struct to the device via the
# backend's array converter (from detect_backend). Structs are already built in the
# device precision (FT), so this is a pure host->device move (RI/CI Int arrays too).
struct DevAdaptor{F}; f::F; end
CLM.Adapt.adapt_storage(d::DevAdaptor, x::AbstractArray) = d.f(x)

function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + max(abs(Float64(A[i])), abs(Float64(B[i])))))
    end
    return m
end

# Shared toy problem: 3 units, 4-pool state, 3 off-diagonal transfer entries.
const BEGU, ENDU = 1, 3
const NUNIT, SM, NUM = 3, 4, 3
const FILT = [1, 2, 3]
const RI, CI, NE = [2, 3, 4], [1, 2, 3], 3          # transfers 2<-1, 3<-2, 4<-3
const AI, AJ, NENON = [2, 3, 4], [1, 2, 3], 3       # set_value_a! off-diagonal (row,col)

# ---- CPU golden builders (Float64) + their GPU twins (Float32 via `dev`) ------
mkA(FT) = (A = CLM.SparseMatrixType(); CLM.init_sm!(A, SM, BEGU, ENDU; ref = nothing, FT = FT);
           A.RI[1:NE] .= RI; A.CI[1:NE] .= CI; A.NE = NE;
           for u in 1:NUNIT, k in 1:NE; A.M[u, k] = FT(0.1k + 0.05u); end; A)
mkK(FT) = (K = CLM.DiagMatrixType(); CLM.init_dm!(K, SM, BEGU, ENDU; ref = nothing, FT = FT);
           K)  # values set via set_value_dm!
mkX(FT) = (X = CLM.VectorType(); CLM.init_v!(X, SM, BEGU, ENDU; ref = nothing, FT = FT); X)
Ksrc(FT) = FT[0.5 + 0.1i for u in 1:NUNIT, i in 1:SM]
Xsrc(FT) = FT[100.0 + 10i + u for u in 1:NUNIT, i in 1:SM]
Bsrc(FT) = FT[1.0 + 0.2k + u for u in 1:NUNIT, k in 1:NE]
Asrc(FT) = FT[0.1k + 0.05u for u in 1:NUNIT, k in 1:NENON]

function main(backend)
    println("=" ^ 68); println("REAL-GPU parity for the matrix-CN sparse operator layer"); println("=" ^ 68)
    if backend === nothing; println("  No GPU backend available — nothing to prove."); return 0; end
    name, dev, FT = backend
    @printf("  Backend: %s   (precision: %s)\n", name, FT)
    a(x) = CLM.Adapt.adapt(DevAdaptor(dev), x)   # move a whole sparse struct to the device
    checks = Tuple{String,Any,Any}[]

    # (1) set_value_dm! : diagonal K fill  ------------------------------------
    K64 = mkK(Float64); CLM.set_value_dm!(K64, BEGU, ENDU, NUM, FILT, Ksrc(Float64))
    Km  = a(mkK(FT));   CLM.set_value_dm!(Km, BEGU, ENDU, NUM, dev(FILT), dev(Ksrc(FT)))
    push!(checks, ("set_value_dm!", K64.DM, Km.DM))

    # (2) set_value_v! + set_value_v_scaler! : vector fill then in-place scale --
    X64 = mkX(Float64); CLM.set_value_v!(X64, BEGU, ENDU, NUM, FILT, Xsrc(Float64))
    CLM.set_value_v_scaler!(X64, NUM, FILT, 2.0)
    Xm  = a(mkX(FT));   CLM.set_value_v!(Xm, BEGU, ENDU, NUM, dev(FILT), dev(Xsrc(FT)))
    CLM.set_value_v_scaler!(Xm, NUM, dev(FILT), 2.0)
    push!(checks, ("set_value_v!+scaler", X64.V, Xm.V))

    # (3) spmm_ak! then spmm_ax! : the full A*K*X pool-advance operator ---------
    A64 = mkA(Float64); Kk64 = mkK(Float64); CLM.set_value_dm!(Kk64, BEGU, ENDU, NUM, FILT, Ksrc(Float64))
    Xx64 = mkX(Float64); CLM.set_value_v!(Xx64, BEGU, ENDU, NUM, FILT, Xsrc(Float64))
    CLM.spmm_ak!(A64, NUM, FILT, Kk64); CLM.spmm_ax!(Xx64, NUM, FILT, A64)
    Am = a(mkA(FT)); Kkm = a(mkK(FT)); CLM.set_value_dm!(Kkm, BEGU, ENDU, NUM, dev(FILT), dev(Ksrc(FT)))
    Xxm = a(mkX(FT)); CLM.set_value_v!(Xxm, BEGU, ENDU, NUM, dev(FILT), dev(Xsrc(FT)))
    CLM.spmm_ak!(Am, NUM, dev(FILT), Kkm); CLM.spmm_ax!(Xxm, NUM, dev(FILT), Am)
    push!(checks, ("spmm_ak!*A", A64.M[:, 1:NE], Am.M[:, 1:NE]))
    push!(checks, ("spmm_ax! (X<-X+A*X)", Xx64.V, Xxm.V))

    # (4) spmp_b_acc! : same-structure sparse accumulate B <- B + A -------------
    Ba64 = mkA(Float64); Bb64 = mkA(Float64)
    for u in 1:NUNIT, k in 1:NE; Bb64.M[u, k] = Bsrc(Float64)[u, k]; end
    CLM.spmp_b_acc!(Bb64, NUM, FILT, Ba64)
    Bam = a(mkA(FT)); Bbm = a(mkA(FT))
    Bfill = dev(Bsrc(FT)); copyto!(view(Bbm.M, :, 1:NE), Bfill)
    CLM.spmp_b_acc!(Bbm, NUM, dev(FILT), Bam)
    push!(checks, ("spmp_b_acc!", Bb64.M[:, 1:NE], Bbm.M[:, 1:NE]))

    # (5) set_value_a! memoized value-fill : build structure on CPU, then run the
    #     kernelized (init_ready=true) fill on a DEVICE matrix w/ device list/RI/CI
    Aa = CLM.SparseMatrixType(); CLM.init_sm!(Aa, SM, BEGU, ENDU)
    list = fill(0, SM + NENON); RIa = fill(0, SM + NENON); CIa = fill(0, SM + NENON)
    ir = CLM.set_value_a!(Aa, BEGU, ENDU, NUM, FILT, Asrc(Float64), AI, AJ, NENON, false;
                          list = list, RI_A = RIa, CI_A = CIa)
    CLM.set_value_a!(Aa, BEGU, ENDU, NUM, FILT, Asrc(Float64), AI, AJ, NENON, ir;
                     list = list, RI_A = RIa, CI_A = CIa)
    Ad = CLM.SparseMatrixType(); CLM.init_sm!(Ad, SM, BEGU, ENDU; ref = nothing, FT = FT); Adm = a(Ad)
    CLM.set_value_a!(Adm, BEGU, ENDU, NUM, dev(FILT), dev(Asrc(FT)), AI, AJ, NENON, true;
                     list = dev(list), RI_A = dev(RIa), CI_A = dev(CIa))
    push!(checks, ("set_value_a! (memoized fill)", Aa.M[:, 1:Aa.NE], Adm.M[:, 1:Aa.NE]))

    # (6) veg-solver GLUE kernels — the per-patch load / B-input / Aoned / add /
    #     write-back loops that replaced the host scalar-index loops in the solve
    #     body. Each runs on CPU (golden) and device over the same operands.
    np = 3; nveg = 18; nnon = 5; filt = [1, 2, 3]; nump = 3; dt = 1800.0; doner = [1, 4, 7, 10, 13]
    # aoned: Aoned[p,k] = transfer·dt / turnover[doner[k]]
    let tr = FT[0.001(k + p) for p in 1:np, k in 1:nnon], tu = FT[0.05(d + p) for p in 1:np, d in 1:nveg]
        Ah = zeros(np, nnon); CLM._launch!(CLM._veg_aoned_kernel!, Ah, Float64.(tr), Float64.(tu), doner, filt, nnon, 1, dt; ndrange = nump)
        Ag = dev(zeros(FT, np, nnon)); CLM._launch!(CLM._veg_aoned_kernel!, Ag, dev(tr), dev(tu), dev(doner), dev(filt), nnon, 1, FT(dt); ndrange = nump)
        push!(checks, ("glue _veg_aoned_kernel!", Ah, Ag))
    end
    # binput: B[p,i] = alloc·input·dt
    let al = FT[0.1(i + p) for p in 1:np, i in 1:nveg], inp = FT[2.0 + p for p in 1:np]
        Bh = zeros(np, nveg); CLM._launch!(CLM._veg_binput_kernel!, Bh, Float64.(al), Float64.(inp), filt, nveg, 1, dt; ndrange = nump)
        Bg = dev(zeros(FT, np, nveg)); CLM._launch!(CLM._veg_binput_kernel!, Bg, dev(al), dev(inp), dev(filt), nveg, 1, FT(dt); ndrange = nump)
        push!(checks, ("glue _veg_binput_kernel!", Bh, Bg))
    end
    # add: X += B
    let X0 = FT[100.0 + i + p for p in 1:np, i in 1:nveg], Bv = FT[i + 0.5p for p in 1:np, i in 1:nveg]
        Xh = Float64.(X0); CLM._launch!(CLM._veg_add_vec_kernel!, Xh, Float64.(Bv), filt, nveg, 1; ndrange = nump)
        Xg = dev(copy(X0)); CLM._launch!(CLM._veg_add_vec_kernel!, Xg, dev(Bv), dev(filt), nveg, 1; ndrange = nump)
        push!(checks, ("glue _veg_add_vec_kernel!", Xh, Xg))
    end
    # load_c / write-back_c round-trip (crop inactive): fields -> Xvegc.V -> fields
    let fields = [FT[10.0k + p for p in 1:np] for k in 1:18], repro = [FT[0.0 for p in 1:np, r in 1:1] for _ in 1:3], ivt = [1, 2, 3]
        Xh = fill(-9999.0, np, nveg)
        CLM._launch!(CLM._veg_load_c_kernel!, Xh, (Float64.(f) for f in fields)..., (Float64.(r) for r in repro)..., ivt, filt, 1, 15, false, 1; ndrange = nump)
        Xg = dev(fill(FT(-9999), np, nveg))
        CLM._launch!(CLM._veg_load_c_kernel!, Xg, (dev(f) for f in fields)..., (dev(r) for r in repro)..., dev(ivt), dev(filt), 1, 15, false, 1; ndrange = nump)
        push!(checks, ("glue _veg_load_c_kernel!", Xh, Xg))
        # write-back into fresh field arrays and compare the 18 native pools
        oh = [zeros(np) for _ in 1:18]; ohr = [zeros(np, 1) for _ in 1:3]
        CLM._launch!(CLM._veg_writeback_c_kernel!, oh..., ohr..., Xh, ivt, filt, 1, 15, false, 1; ndrange = nump)
        og = [dev(zeros(FT, np)) for _ in 1:18]; ogr = [dev(zeros(FT, np, 1)) for _ in 1:3]
        CLM._launch!(CLM._veg_writeback_c_kernel!, og..., ogr..., Xg, dev(ivt), dev(filt), 1, 15, false, 1; ndrange = nump)
        push!(checks, ("glue _veg_writeback_c_kernel!", reduce(hcat, oh), reduce(hcat, [Array(x) for x in og])))
    end

    nfail = 0
    for (nm, cpu, gpu) in checks
        dd = reldiff(cpu, gpu); ok = dd < 1f-3
        @printf("  [%s] %-30s rel = %.3e\n", ok ? "PASS" : "FAIL", nm, dd); ok || (nfail += 1)
    end
    println()
    println(nfail == 0 ? "  matrix-CN sparse ops MATCH CPU ON $name ($FT)" : "  DIVERGENCE ($nfail op(s)).")
    return nfail == 0 ? 0 : 1
end

exit(main(detect_backend()))
