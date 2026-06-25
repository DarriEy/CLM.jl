# ==========================================================================
# Tests for on-node intra-process multi-clump parallelism — CTSM's
# OpenMP-over-clumps (`do nc = 1, nclumps`) model ported to the CLM.jl driver.
#
# The clump loop is correct iff:
#   (1) decompInit!(numg, clump_pproc=N) partitions the proc domain into N
#       disjoint, contiguous, complete clumps (get_clump_bounds(nc));
#   (2) a clump-safe per-column physics kernel run threaded over those clumps
#       produces EXACTLY the same state as the single-clump (whole-proc) run —
#       i.e. the answer is independent of clump count and threading
#       (DETERMINISM / race-freedom);
#   (3) each clump touches only its disjoint slice;
#   (4) the nclumps == 1 path is byte-identical to the proc-bounds path.
# ==========================================================================

using Test
using CLM
const _C = CLM

@testset "threaded multi-clump driver loop" begin

    # ----------------------------------------------------------------------
    # (1) Decomposition partitions the proc domain into disjoint clumps.
    # ----------------------------------------------------------------------
    @testset "decomposition partitions proc into disjoint clumps" begin
        numg = 12
        d = _C.DecompData()
        _C.decompInit!(numg; clump_pproc = 4, decomp_data = d)

        @test d.nclumps == 4
        @test _C.get_proc_clumps(decomp_data = d) == 4

        bp = _C.get_proc_bounds(decomp_data = d)
        covc = Int[]; covp = Int[]; covg = Int[]
        for nc in 1:_C.get_proc_clumps(decomp_data = d)
            b = _C.get_clump_bounds(nc; decomp_data = d)
            @test b.level == _C.BOUNDS_LEVEL_CLUMP
            @test b.clump_index == nc
            @test b.begc <= b.endc
            append!(covc, b.begc:b.endc)
            append!(covp, b.begp:b.endp)
            append!(covg, b.begg:b.endg)
        end
        # disjoint (no repeats) AND complete (covers the whole proc)
        @test length(unique(covc)) == length(covc)          # disjoint columns
        @test sort(covc) == collect(bp.begc:bp.endc)         # complete columns
        @test sort(covp) == collect(bp.begp:bp.endp)         # complete patches
        @test sort(covg) == collect(bp.begg:bp.endg)         # complete gridcells
    end

    # ----------------------------------------------------------------------
    # A clump-safe per-column physics kernel, written in the EXACT pattern the
    # CLM driver uses: iterate bounds_clump.begc:endc, gate on a proc-length
    # (absolute-indexed) filter mask, and write only into the disjoint proc-level
    # state slice for those columns. No cross-clump reads/writes, no shared
    # global mutation -> the contract clm_run_clump_physics! requires.
    # ----------------------------------------------------------------------
    function build_state(nc::Int)
        # proc-length state + filter (absolute-indexed, exactly like the driver)
        x      = collect(Float64, 1:nc) .* 0.5            # input
        y      = fill(-999.0, nc)                          # output (state slice)
        mask   = BitVector(isodd(c) for c in 1:nc)        # filt.soilc-style mask
        touch  = fill(-1, nc)                              # which clump wrote col c
        return (x = x, y = y, mask = mask, touch = touch)
    end

    function make_phys!(st, clump_tag::Dict{_C.BoundsType,Int})
        # closure: per-clump column physics. Reads st.x, writes st.y/st.touch only
        # for the columns in bounds_clump that pass the mask. clump_tag lets us
        # record which clump processed each column (proves disjoint coverage).
        return bc -> begin
            tag = clump_tag[bc]
            for c in bc.begc:bc.endc
                st.mask[c] || continue
                # a nonlinear-ish per-column update (order-independent across columns)
                st.y[c]     = sqrt(st.x[c]) + st.x[c]^2 - 3.0 * st.x[c]
                st.touch[c] = tag
            end
        end
    end

    # ----------------------------------------------------------------------
    # (2) + (4) multi-clump (threaded) == single-clump == proc-bounds reference.
    # ----------------------------------------------------------------------
    @testset "multi-clump result == single-clump result (determinism)" begin
        numg = 12
        d = _C.DecompData()
        _C.decompInit!(numg; clump_pproc = 4, decomp_data = d)
        nc = d.numc
        bp = _C.get_proc_bounds(decomp_data = d)

        # --- reference: single whole-proc clump, serial ---
        ref = build_state(nc)
        ref_tag = Dict(bp => 0)
        _C.clm_run_clump_physics!(make_phys!(ref, ref_tag), [bp]; threaded = false)

        # --- nclumps == 1 path is byte-identical to the proc-bounds path ---
        d1 = _C.DecompData()
        _C.decompInit!(numg; clump_pproc = 1, decomp_data = d1)
        b1 = _C.get_clump_bounds(1; decomp_data = d1)
        @test (b1.begc, b1.endc) == (bp.begc, bp.endc)
        @test (b1.begp, b1.endp) == (bp.begp, bp.endp)
        @test (b1.begg, b1.endg) == (bp.begg, bp.endg)
        single = build_state(nc)
        single_tag = Dict(b1 => 1)
        _C.clm_run_clump_physics!(make_phys!(single, single_tag), [b1]; threaded = false)
        @test single.y == ref.y                              # byte-identical

        # --- multi-clump (4 clumps), THREADED ---
        clump_bounds = [_C.get_clump_bounds(n; decomp_data = d) for n in 1:_C.get_proc_clumps(decomp_data = d)]
        multi = build_state(nc)
        tags = Dict(clump_bounds[i] => i for i in eachindex(clump_bounds))
        _C.clm_run_clump_physics!(make_phys!(multi, tags), clump_bounds; threaded = true)

        # same answer regardless of clump count / threading
        @test multi.y == ref.y

        # (3) each clump processed ONLY its disjoint slice: a column written by
        # clump i must lie inside clump i's begc:endc.
        for c in 1:nc
            multi.mask[c] || continue
            t = multi.touch[c]
            @test t >= 1
            b = clump_bounds[t]
            @test b.begc <= c <= b.endc
        end
    end

    # ----------------------------------------------------------------------
    # (2') No race: run the threaded multi-clump loop many times, every run must
    # produce the identical result (a data race would surface as nondeterminism).
    # ----------------------------------------------------------------------
    @testset "threaded loop is race-free (repeated runs identical)" begin
        numg = 60
        d = _C.DecompData()
        _C.decompInit!(numg; clump_pproc = 6, decomp_data = d)
        nc = d.numc
        bp = _C.get_proc_bounds(decomp_data = d)

        ref = build_state(nc)
        _C.clm_run_clump_physics!(make_phys!(ref, Dict(bp => 0)), [bp]; threaded = false)

        clump_bounds = [_C.get_clump_bounds(n; decomp_data = d) for n in 1:_C.get_proc_clumps(decomp_data = d)]
        for trial in 1:50
            st = build_state(nc)
            tags = Dict(clump_bounds[i] => i for i in eachindex(clump_bounds))
            _C.clm_run_clump_physics!(make_phys!(st, tags), clump_bounds; threaded = true)
            @test st.y == ref.y          # identical every trial -> no race
        end
    end

    # ----------------------------------------------------------------------
    # Threaded vs serial dispatch of the SAME clump loop agree (the threaded=
    # true/false branches of clm_run_clump_physics! are equivalent).
    # ----------------------------------------------------------------------
    @testset "threaded == serial dispatch" begin
        numg = 24
        d = _C.DecompData()
        _C.decompInit!(numg; clump_pproc = 3, decomp_data = d)
        nc = d.numc
        clump_bounds = [_C.get_clump_bounds(n; decomp_data = d) for n in 1:_C.get_proc_clumps(decomp_data = d)]

        a = build_state(nc); ta = Dict(clump_bounds[i] => i for i in eachindex(clump_bounds))
        b = build_state(nc); tb = Dict(clump_bounds[i] => i for i in eachindex(clump_bounds))
        _C.clm_run_clump_physics!(make_phys!(a, ta), clump_bounds; threaded = false)
        _C.clm_run_clump_physics!(make_phys!(b, tb), clump_bounds; threaded = true)
        @test a.y == b.y
        @test a.touch == b.touch
    end

    # ----------------------------------------------------------------------
    # Driver config plumbing: use_threaded_clumps defaults off (byte-identical
    # single-clump path), is settable, and the driver's clump-bounds list
    # collapses to [bounds_proc] when off / single-clump.
    # ----------------------------------------------------------------------
    @testset "config use_threaded_clumps default + plumbing" begin
        cfg_default = _C.CLMDriverConfig()
        @test cfg_default.use_threaded_clumps == false

        cfg_on = _C.CLMDriverConfig(use_threaded_clumps = true)
        @test cfg_on.use_threaded_clumps == true

        # the driver's per-clump bounds selection (mirrors clm_drv_core!):
        # off  -> [bounds_proc]; on + nclumps>1 -> the clump partition.
        numg = 12
        d = _C.DecompData()
        _C.decompInit!(numg; clump_pproc = 4, decomp_data = d)
        bp = _C.get_proc_bounds(decomp_data = d)

        sel(use_threaded, nclumps) =
            (use_threaded && nclumps > 1) ?
                [_C.get_clump_bounds(n; decomp_data = d) for n in 1:nclumps] : [bp]

        @test sel(false, _C.get_proc_clumps(decomp_data = d)) == [bp]      # default -> single
        @test length(sel(true, _C.get_proc_clumps(decomp_data = d))) == 4  # threaded -> 4 clumps
        # the threaded selection still partitions the proc (column ranges sum to proc size)
        partition = sel(true, _C.get_proc_clumps(decomp_data = d))
        @test sum(b -> b.endc - b.begc + 1, partition) == (bp.endc - bp.begc + 1)
    end

end
