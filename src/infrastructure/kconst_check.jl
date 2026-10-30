# ==========================================================================
# Opt-in @Const aliasing checker for `_launch!` (debug only).
#
# KernelAbstractions' `@Const(x)` promises that `x` is read-only AND does not alias
# any array the kernel writes; on the CPU backend (KA < 0.9.45, or any KA on Julia
# 1.10) it is lowered through `Base.Experimental.Const` + `@aliasscope`, so a broken
# promise is undefined behaviour (LLVM may reorder loads past stores).
#
# With the checker on, every `_launch!` on the CPU backend:
#   * snapshots each @Const array argument, runs the kernel, and reports any @Const
#     argument whose contents CHANGED (written by the kernel, directly or through an
#     alias) -- kind `:const_written`;
#   * reports every (@Const, non-@Const) array-argument pair with
#     `Base.mightalias` true -- kind `:alias_written` when the non-@Const one changed
#     during the launch, `:alias_readonly` when it did not (harmless today, fragile).
# The @Const positions come from parsing the `@kernel` definitions in `src/`.
# Kernels launched directly (`k!(backend)(...)`, not through `_launch!`) are not
# covered; test/test_kconst_alias.jl audits every `@kernel` body statically.
#
# Enable: `CLM.set_kconst_check!(:warn)` (record + warn once per kernel/kind),
# `:error` (throw on the first hazard), `:off` (default). Or set the environment
# variable CLM_CHECK_KCONST=warn|error before `using CLM`. The mode is a method
# constant, so the default (`:off`) path compiles to exactly the unchecked launch.
# Findings accumulate in `CLM.KCONST_FINDINGS` (cleared by `clear_kconst_findings!`).
# ==========================================================================

kconst_check_mode() = :off

"""
    set_kconst_check!(mode::Symbol)

Turn the `@Const` aliasing checker in `_launch!` on (`:warn`, `:error`) or off
(`:off`). Redefines a method, so call it at top level (not from inside a function
that then launches kernels in the same world age).
"""
function set_kconst_check!(mode::Symbol)
    mode in (:off, :warn, :error) || throw(ArgumentError("mode must be :off, :warn or :error"))
    @eval CLM kconst_check_mode() = $(QuoteNode(mode))
    return mode
end

const KCONST_FINDINGS = NamedTuple{(:kernel, :kind, :const_arg, :other_arg),
                                   Tuple{Symbol,Symbol,Symbol,Symbol}}[]
const _KCONST_SEEN = Set{Tuple{Symbol,Symbol,Symbol,Symbol}}()
clear_kconst_findings!() = (empty!(KCONST_FINDINGS); empty!(_KCONST_SEEN); nothing)

# kernel name => (argument names, @Const flags), positional, excluding __ctx__.
const _KCONST_REGISTRY = Ref{Union{Nothing,Dict{Symbol,Tuple{Vector{Symbol},Vector{Bool}}}}}(nothing)

_kc_argname(a::Symbol) = a
_kc_argname(a::Expr) = a.head === :(::) ? (length(a.args) == 2 ? _kc_argname(a.args[1]) : :_) :
                       a.head === :kw ? _kc_argname(a.args[1]) :
                       (a.head === :macrocall ? _kc_argname(a.args[end]) : :_)
_kc_isconst(a) = a isa Expr && a.head === :macrocall && a.args[1] === Symbol("@Const")

"Parse every `@kernel function` signature in `path` (a file, or every .jl file under a
directory; default `src/`): name => (argnames, isconst)."
function kconst_kernel_signatures(path::AbstractString = joinpath(@__DIR__, ".."))
    reg = Dict{Symbol,Tuple{Vector{Symbol},Vector{Bool}}}()
    files = isfile(path) ? [path] :
        [joinpath(r, f) for (r, _, fs) in walkdir(path) for f in fs if endswith(f, ".jl")]
    for file in files
        ex = Meta.parseall(read(file, String); filename = file)
        _kc_walk!(ex) do def
            sig = def.args[1]
            while sig isa Expr && sig.head === :where
                sig = sig.args[1]
            end
            (sig isa Expr && sig.head === :call && sig.args[1] isa Symbol) || return
            args = [a for a in sig.args[2:end] if !(a isa Expr && a.head === :parameters)]
            reg[sig.args[1]] = (Symbol[_kc_argname(a) for a in args], Bool[_kc_isconst(a) for a in args])
        end
    end
    return reg
end

# Call `f(funcdef)` for the function definition of every `@kernel` macro call.
function _kc_walk!(f, ex)
    ex isa Expr || return
    if ex.head === :macrocall && ex.args[1] === Symbol("@kernel")
        def = ex.args[end]
        def isa Expr && def.head === :function && f(def)
        return
    end
    foreach(a -> _kc_walk!(f, a), ex.args)
end

function _kconst_registry()
    r = _KCONST_REGISTRY[]
    r === nothing || return r
    return (_KCONST_REGISTRY[] = kconst_kernel_signatures())
end

function _kconst_report!(kname, kind, ca, oa)
    key = (kname, kind, ca, oa)
    push!(KCONST_FINDINGS, (kernel = kname, kind = kind, const_arg = ca, other_arg = oa))
    if kconst_check_mode() === :error && kind !== :alias_readonly
        error("@Const hazard in $kname: $kind ($ca vs $oa)")
    end
    key in _KCONST_SEEN && return
    push!(_KCONST_SEEN, key)
    kind === :alias_readonly ||
        @warn "@Const hazard in kernel $kname: $kind" const_arg = ca other_arg = oa
    return
end

_kc_addr(x, I) = UInt(pointer(x)) + sizeof(eltype(x)) * sum((Tuple(I) .- 1) .* strides(x); init = 0)
_kc_addrs(x) = Set(_kc_addr(x, I) for I in CartesianIndices(x))
_kc_lo(x) = minimum(_kc_addr(x, I) for I in (first(CartesianIndices(x)), last(CartesianIndices(x))))
_kc_hi(x) = maximum(_kc_addr(x, I) for I in (first(CartesianIndices(x)), last(CartesianIndices(x)))) +
            sizeof(eltype(x)) - 1

# Memory overlap of two kernel arguments. Strided bits arrays: compare their byte
# extents (so disjoint column views of one packed scratch array do NOT count);
# anything else: `Base.mightalias` (conservative).
function _kc_overlap(a, b)
    Base.mightalias(a, b) || return false
    if a isa StridedArray && b isa StridedArray && !(a isa BitArray) && !(b isa BitArray) &&
       isbitstype(eltype(a)) && isbitstype(eltype(b)) && !isempty(a) && !isempty(b)
        (_kc_hi(a) < _kc_lo(b) || _kc_hi(b) < _kc_lo(a)) && return false
        # extents cross: strided slices can interleave without sharing an element
        return !isdisjoint(_kc_addrs(a), _kc_addrs(b))
    end
    return true
end

# Array leaves of a kernel argument: the array itself, or the array fields of a
# bundle struct / NamedTuple / Tuple (recursively, a few levels), with dotted names.
function _kc_leaves!(acc, x, name::Symbol, depth::Int = 0)
    if x isa AbstractArray
        push!(acc, (name, x))
    elseif depth < 3 && (x isa Tuple || x isa NamedTuple || (isstructtype(typeof(x)) &&
            !(x isa Number || x isa Symbol || x isa AbstractString || x isa Function ||
              x isa Type || x isa Module)))
        for f in fieldnames(typeof(x))
            isdefined(x, f) || continue
            _kc_leaves!(acc, getfield(x, f), Symbol(name, ".", f), depth + 1)
        end
    end
    return acc
end

@noinline function _kconst_checked_launch!(kernel, backend, ndrange, kargs::Tuple)
    entry = backend isa KA.CPU ? get(_kconst_registry(), nameof(kernel), nothing) : nothing
    if entry === nothing || length(entry[1]) != length(kargs)
        kernel(backend)(kargs...; ndrange = ndrange)   # GPU, unknown or mismatched kernel
        return
    end
    names, isconst = entry
    cl = Tuple{Symbol,Any}[]; wl = Tuple{Symbol,Any}[]
    for i in eachindex(kargs)
        _kc_leaves!(isconst[i] ? cl : wl, kargs[i], names[i])
    end
    pairs = [(i, j) for i in eachindex(cl) for j in eachindex(wl) if _kc_overlap(cl[i][2], wl[j][2])]
    csnap = [copy(x) for (_, x) in cl]
    wsnap = Dict(j => copy(wl[j][2]) for j in unique(last.(pairs)))
    kernel(backend)(kargs...; ndrange = ndrange)
    KA.synchronize(backend)
    kname = nameof(kernel)
    for i in eachindex(cl)
        isequal(csnap[i], cl[i][2]) || _kconst_report!(kname, :const_written, cl[i][1], :_)
    end
    for (i, j) in pairs
        _kconst_report!(kname, isequal(wsnap[j], wl[j][2]) ? :alias_readonly : :alias_written,
                        cl[i][1], wl[j][1])
    end
    return
end

function __init__()
    m = Symbol(lowercase(get(ENV, "CLM_CHECK_KCONST", "off")))
    m === Symbol("1") && (m = :error)
    m === :off || set_kconst_check!(m)
end

# ---- static audit: is a @Const argument (or a local view/alias of one) written? ----
const _KC_UPDATE_OPS = (:(=), :(+=), :(-=), :(*=), :(/=), :(.=), :(.+=), :(.-=), :(.*=),
                        :(./=), :(^=), :(|=), :(&=))
const _KC_MUTATORS = (:_scatter_add!, :fill!, :copyto!, :setindex!, :push!, :empty!,
                      :resize!, :sort!, :reverse!, :map!, :broadcast!)

_kc_base(x::Symbol) = x
_kc_base(x::Expr) = x.head in (:ref, :., :macrocall) && !isempty(x.args) ?
    (x.head === :macrocall ? _kc_base(x.args[end]) : _kc_base(x.args[1])) : nothing
_kc_base(x) = nothing
_kc_viewsrc(rhs) = rhs isa Symbol ? rhs :
    (rhs isa Expr && rhs.head === :macrocall && rhs.args[1] in (Symbol("@view"), Symbol("@views"))) ?
        _kc_base(rhs.args[end]) :
    (rhs isa Expr && rhs.head === :call && rhs.args[1] in (:view, :reshape, :vec, :reinterpret, :PermutedDimsArray)) ?
        _kc_base(rhs.args[2]) : nothing

"""
    kconst_static_writes(path = src/) -> Vector of (kernel, const_arg, via, line)

For every `@kernel` under `path`, report each statement that writes to a `@Const`
argument -- directly (`A[i] = …`, `A[i] += …`, `_scatter_add!(A, …)`, `fill!(A, …)`,
`@atomic A[i] += …`) or through a local alias bound to it (`v = @view A[:, 1]`).
"""
function kconst_static_writes(path::AbstractString = joinpath(@__DIR__, ".."))
    out = NamedTuple{(:kernel, :const_arg, :via, :file, :line),Tuple{Symbol,Symbol,Symbol,String,Int}}[]
    files = isfile(path) ? [path] :
        [joinpath(r, f) for (r, _, fs) in walkdir(path) for f in fs if endswith(f, ".jl")]
    for file in files
        ex = Meta.parseall(read(file, String); filename = file)
        _kc_walk!(ex) do def
            sig = def.args[1]
            while sig isa Expr && sig.head === :where
                sig = sig.args[1]
            end
            (sig isa Expr && sig.head === :call) || return
            kname = sig.args[1] isa Symbol ? sig.args[1] : :_
            alias = Dict{Symbol,Symbol}(_kc_argname(a) => _kc_argname(a)
                                        for a in sig.args[2:end] if _kc_isconst(a))
            isempty(alias) && return
            # propagate local aliases (two passes cover chained views)
            for _ in 1:2
                _kc_foreach(def.args[2]) do e, _
                    if e.head === :(=) && e.args[1] isa Symbol
                        src = _kc_viewsrc(e.args[2])
                        src !== nothing && haskey(alias, src) && (alias[e.args[1]] = alias[src])
                    end
                end
            end
            _kc_foreach(def.args[2]) do e, line
                tgt = nothing
                if e.head in _KC_UPDATE_OPS && !(e.head === :(=) && e.args[1] isa Symbol)
                    tgt = _kc_base(e.args[1])
                elseif e.head === :call && e.args[1] in _KC_MUTATORS && length(e.args) >= 2
                    tgt = _kc_base(e.args[2])
                end
                if tgt !== nothing && haskey(alias, tgt)
                    push!(out, (kernel = kname, const_arg = alias[tgt], via = tgt, file = file, line = line))
                end
            end
        end
    end
    return out
end

function _kc_foreach(f, ex, line = 0)
    ex isa Expr || return line
    for a in ex.args
        a isa LineNumberNode && (line = a.line)
    end
    f(ex, line)
    l = line
    for a in ex.args
        a isa LineNumberNode ? (l = a.line) : _kc_foreach(f, a, l)
    end
    return line
end
