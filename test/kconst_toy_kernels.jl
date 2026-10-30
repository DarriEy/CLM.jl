# Toy kernels for test_kconst_alias.jl: each breaks (or keeps) the @Const contract in
# a known way. Parsed by CLM.kconst_kernel_signatures / kconst_static_writes and
# included (evaluated) by the test.
using KernelAbstractions: @kernel, @index, @Const

@kernel function _kct_copy!(out, @Const(a))            # contract kept (if not aliased)
    i = @index(Global)
    @inbounds out[i] = 2 * a[i]
end

@kernel function _kct_write_const!(out, @Const(a))     # writes its @Const argument
    i = @index(Global)
    @inbounds a[i] += 1
    @inbounds out[i] = a[i]
end

@kernel function _kct_write_view!(out, @Const(A))      # writes through a view of it
    i = @index(Global)
    v = @view A[:, 1]
    @inbounds v[i] = 0
    @inbounds out[i] = A[i, 2]
end
