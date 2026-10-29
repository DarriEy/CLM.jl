# Test-file isolation for the in-process suite (`runtests.jl`).
#
# CLM.jl mirrors CTSM/FATES module state as module-level globals: parameter
# structs (`pftcon`, `varctl`, `varpar`, `EDPftvarcon_inst`, ...), control flags
# held in `Ref`s, namelist knobs such as `cnvegcstate_const.initial_vegC`, the
# backend registry, and a handful of non-`const` bindings. Many test files write
# these to set up a scenario and do not put them back, so a file's result used to
# depend on which files ran before it.
#
# Two mechanisms make every file see the same starting state, whatever the order:
#
#   1. `snapshot_clm_globals()` records every module-level binding of `CLM`
#      right after `using CLM` (≈0.3 MB). `reset_clm_globals!()` puts the live
#      objects back to that state *in place*: a mutable struct gets its fields
#      re-assigned, a `Ref` its value, a container its contents, so any code that
#      holds a reference to the global keeps seeing the restored object. Non-const
#      bindings are re-bound.
#
#   2. `run_test_file(path)` evaluates each file in a fresh anonymous module that
#      only has `Test` and `CLM` loaded, the way SafeTestsets does. Helpers, consts
#      and methods a file defines therefore cannot bleed into, or be overridden
#      by, another file's (e.g. a more specific method of a same-named helper).
#
# `run_test_file` resets the CLM globals before each file. A file that needs
# state from setup it doesn't do itself fails in every order, not only in some,
# so the dependency shows up and gets fixed in that file.
#
# Physics code is untouched: this file only lives under test/.

module CLMTestIsolation

using Test
import CLM

export snapshot_clm_globals, reset_clm_globals!, run_test_file, ordered_test_files

# Values that are never state worth restoring (or cannot be restored safely).
_skip_value(v) = v isa Module || v isa Function || v isa Type || v isa UnionAll ||
                 v isa IO || v isa Base.AbstractLock || v isa Task ||
                 v isa Channel || v isa Base.Condition || v isa Base.GenericCondition ||
                 v isa Core.MethodInstance || v isa Core.CodeInstance

# Non-mutable leaves: nothing below them can change.
_is_leaf(v) = isbits(v) || v isa Symbol || v isa String || v isa Nothing ||
              v isa Missing || v isa BigInt || v isa BigFloat || v isa Regex

struct GlobalSnapshot
    mod::Module
    consts::Vector{Tuple{Symbol,Any}}      # (name, deepcopy of value)
    nonconsts::Vector{Tuple{Symbol,Any}}
end

"""
    snapshot_clm_globals(mod = CLM) -> GlobalSnapshot

Record every module-level binding of `mod` that can change after load: mutable
`const` objects (structs, `Ref`s, arrays, dicts, sets, and immutable wrappers
around them) and every non-`const` global.
"""
function snapshot_clm_globals(mod::Module = CLM)
    consts = Tuple{Symbol,Any}[]
    nonconsts = Tuple{Symbol,Any}[]
    for n in names(mod; all = true, imported = false)
        s = String(n)
        startswith(s, '#') && continue            # docstring metadata, gensyms
        isdefined(mod, n) || continue
        v = getfield(mod, n)
        _skip_value(v) && continue
        if isconst(mod, n)
            _is_leaf(v) && continue
            push!(consts, (n, deepcopy(v)))
        else
            push!(nonconsts, (n, deepcopy(v)))
        end
    end
    return GlobalSnapshot(mod, consts, nonconsts)
end

# Restore `dst` (a live object) to `src` (a private copy, consumed here).
function _restore!(dst, src)
    dst === src && return dst
    _is_leaf(dst) && return dst
    _skip_value(dst) && return dst
    if dst isa Threads.Atomic
        dst[] = src[]
    elseif dst isa Base.RefValue
        dst[] = src[]
    elseif dst isa AbstractDict
        empty!(dst); merge!(dst, src)
    elseif dst isa AbstractSet
        empty!(dst); union!(dst, src)
    elseif dst isa AbstractVector && applicable(resize!, dst, 0)
        resize!(dst, length(src)); copyto!(dst, src)
    elseif dst isa AbstractArray
        if size(dst) == size(src)
            copyto!(dst, src)
        else
            error("cannot restore a global array whose shape changed: $(size(src)) → $(size(dst))")
        end
    elseif ismutable(dst)
        T = typeof(dst)
        for i in 1:fieldcount(T)
            isdefined(src, i) || continue
            sv = getfield(src, i)
            if isconst(T, i)
                isdefined(dst, i) && _restore!(getfield(dst, i), sv)
            elseif Base.isfieldatomic(T, i)
                setfield!(dst, i, sv, :sequentially_consistent)
            else
                setfield!(dst, i, sv)
            end
        end
    else
        # Immutable wrapper (Tuple, NamedTuple, immutable struct): restore the
        # mutable things it holds.
        for i in 1:nfields(dst)
            isdefined(dst, i) && isdefined(src, i) || continue
            _restore!(getfield(dst, i), getfield(src, i))
        end
    end
    return dst
end

const _COUNTER = Ref(0)
const _PRISTINE = Ref{Union{Nothing,GlobalSnapshot}}(nothing)

"""
    reset_clm_globals!(snap = <snapshot taken when this file was loaded>)

Put every recorded CLM global back to its snapshot state, in place.
"""
function reset_clm_globals!(snap::GlobalSnapshot = _PRISTINE[])
    mod = snap.mod
    for (n, v) in snap.consts
        _restore!(getfield(mod, n), deepcopy(v))
    end
    for (n, v) in snap.nonconsts
        Core.eval(mod, Expr(:(=), n, QuoteNode(deepcopy(v))))
    end
    return nothing
end

"""
    run_test_file(path; reset = true)

`include` the test file at `path` into a fresh module that has only `Test` and
`CLM` loaded, after putting the CLM globals back to their pristine state.
"""
function run_test_file(path::AbstractString; reset::Bool = true)
    reset && reset_clm_globals!()
    _COUNTER[] += 1
    name = Symbol("Sandbox", _COUNTER[], "_", splitext(basename(path))[1])
    m = Core.eval(Main, Expr(:module, true, name, Expr(:block, :(using Test, CLM))))
    Base.include(m, path)
    return nothing
end

"""
    ordered_test_files(files) -> Vector{String}

Apply the order requested by the `CLM_TEST_ORDER` environment variable:
`normal` (default), `reverse`, or `random:<seed>`. `CLM_TEST_ONLY` (a
comma-separated list of file names) restricts the run to those files.
"""
function ordered_test_files(files::AbstractVector{<:AbstractString})
    fs = collect(String, files)
    only = strip(get(ENV, "CLM_TEST_ONLY", ""))
    if !isempty(only)
        keep = Set(strip.(split(only, ',')))
        fs = filter(f -> f in keep, fs)
    end
    order = strip(get(ENV, "CLM_TEST_ORDER", "normal"))
    if order == "reverse"
        reverse!(fs)
    elseif startswith(order, "random:")
        seed = parse(UInt64, order[length("random:")+1:end])
        _shuffle!(fs, seed)
    elseif order != "normal"
        error("CLM_TEST_ORDER must be normal, reverse or random:<seed>; got $(repr(order))")
    end
    return fs
end

# Fisher–Yates with a self-contained SplitMix64 stream, so a seed gives the same
# order on every Julia version (Random's default RNG changed across releases).
function _shuffle!(v::Vector, seed::UInt64)
    s = seed
    for i in length(v):-1:2
        s += 0x9e3779b97f4a7c15
        z = s
        z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
        z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
        z = z ⊻ (z >> 31)
        j = Int(z % UInt64(i)) + 1
        v[i], v[j] = v[j], v[i]
    end
    return v
end

end # module CLMTestIsolation

using .CLMTestIsolation
CLMTestIsolation._PRISTINE[] === nothing &&
    (CLMTestIsolation._PRISTINE[] = snapshot_clm_globals())
