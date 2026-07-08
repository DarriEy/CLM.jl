# ==========================================================================
# gpu_host_loop_census.jl — the Phase-1 census tool for the GPU port. Parses
# every src file's AST and reports each `for`-loop that is NOT inside an @kernel
# function — i.e. the residual HOST loops that still need kernelizing for a
# whole-driver device run. AST-based (not grep) so @kernel bodies, comprehensions,
# and nested functions are classified correctly.
#
# Output: per-file, per-function host-loop counts (functions with 0 are omitted),
# with 1-based line numbers to jump to. A grand total + per-directory rollup.
#
#   julia scripts/gpu_host_loop_census.jl [srcdir1 srcdir2 ...]
#   (default dirs: biogeophys biogeochem infrastructure types driver)
# ==========================================================================

using Printf
const SRC = joinpath(dirname(@__DIR__), "src")
default_dirs = ["biogeophys", "biogeochem", "infrastructure", "types", "driver"]
dirs = isempty(ARGS) ? default_dirs : ARGS

# Is this expr a `@kernel function ...` macrocall?
is_kernel_macro(e) = e isa Expr && e.head === :macrocall && !isempty(e.args) &&
    (e.args[1] === Symbol("@kernel") || (e.args[1] isa Expr && e.args[1].head === :. &&
        e.args[1].args[end] === QuoteNode(Symbol("@kernel"))))

# Count `for` loops (and `while` loops) anywhere under `x`, NOT descending into
# any @kernel function it contains. Returns (count, first_lines::Vector{Int}).
function count_host_loops(x, lines)
    n = 0
    if x isa Expr
        if is_kernel_macro(x)
            return 0            # kernel body: its loops run on-device, skip
        end
        if x.head === :for || x.head === :while
            n += 1
        end
        for a in x.args
            if a isa LineNumberNode
                continue
            end
            c = count_host_loops(a, lines)
            n += c
        end
    end
    return n
end

# Walk toplevel, tracking the current source line, and report per function.
function census_file(path)
    src = read(path, String)
    local ast
    try
        ast = Meta.parseall(src; filename=path)
    catch err
        return [(name="<parse error>", line=0, loops=0)]
    end
    results = Tuple{String,Int,Int}[]
    curline = 0
    function funcname(sig)
        # sig is the function signature expr; dig to the call name
        s = sig
        while s isa Expr && s.head in (:where, :(::)); s = s.args[1]; end
        if s isa Expr && s.head === :call; s = s.args[1]; end
        if s isa Expr && s.head === :curly; s = s.args[1]; end
        return string(s)
    end
    function handle(e)
        if e isa LineNumberNode
            curline = e.line
            return
        end
        e isa Expr || return
        # unwrap @kernel function → skip entirely (device)
        if is_kernel_macro(e); return; end
        if e.head === :function || (e.head === :(=) && e.args[1] isa Expr && e.args[1].head === :call)
            body = e.args[end]
            nm = funcname(e.args[1])
            ln = curline
            nloops = count_host_loops(body, nothing)
            if nloops > 0
                push!(results, (nm, ln, nloops))
            end
            return
        end
        # descend into module/other toplevel wrappers
        for a in e.args; handle(a); end
    end
    if ast isa Expr && ast.head === :toplevel
        for a in ast.args; handle(a); end
    else
        handle(ast)
    end
    return [(name=r[1], line=r[2], loops=r[3]) for r in results]
end

total = 0
rollup = Dict{String,Int}()
for d in dirs
    dpath = joinpath(SRC, d)
    isdir(dpath) || continue
    files = sort(filter(f -> endswith(f, ".jl"), readdir(dpath)))
    dsum = 0
    printed_dir = false
    for f in files
        rows = census_file(joinpath(dpath, f))
        fsum = sum(r.loops for r in rows; init=0)
        fsum == 0 && continue
        if !printed_dir
            println("\n", "="^68, "\n  ", uppercase(d), "\n", "="^68)
            printed_dir = true
        end
        println("\n  $f  —  $fsum host loop(s)")
        for r in sort(rows; by=x -> -x.loops)
            @printf("    L%-5d %-42s %2d\n", r.line, r.name, r.loops)
        end
        dsum += fsum
    end
    rollup[d] = dsum
    global total += dsum
end

println("\n", "="^68, "\n  ROLLUP\n", "="^68)
for d in dirs
    haskey(rollup, d) && @printf("  %-18s %4d host loops\n", d, rollup[d])
end
@printf("  %-18s %4d host loops (candidates for kernelization)\n", "TOTAL", total)
