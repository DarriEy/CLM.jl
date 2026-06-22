# Tests for the minimal CLM-style history (h0) writer in
# src/infrastructure/history_io.jl. Uses a lightweight stand-in `inst` with
# mutable per-element vectors so we can feed known values over N steps,
# write the NetCDF, re-read it, and assert the time-aggregated values.

using Test
using NCDatasets

# A tiny mutable stand-in for CLMInstances: the getters just read named fields.
Base.@kwdef mutable struct _HistFakeInst
    a::Vector{Float64} = Float64[]
    b::Vector{Float64} = Float64[]
end

@testset "history_io" begin
    @testset "avgflag accumulation (A/I/X/M)" begin
        inst = _HistFakeInst(a = [0.0], b = [0.0])
        tape = CLM.HistoryTape()
        CLM.hist_addfld!(tape, "AVG", "u", i -> i.a; avgflag="A")
        CLM.hist_addfld!(tape, "INS", "u", i -> i.a; avgflag="I")
        CLM.hist_addfld!(tape, "MAX", "u", i -> i.a; avgflag="X")
        CLM.hist_addfld!(tape, "MIN", "u", i -> i.a; avgflag="M")

        samples = [10.0, 30.0, 20.0, 0.0]
        for s in samples
            inst.a[1] = s
            CLM.hist_accumulate!(tape, inst)
        end
        @test tape.nsamples == 4

        vals = Dict(f.name => CLM.hist_field_value(f)[1] for f in tape.fields)
        @test vals["AVG"] ≈ sum(samples) / length(samples)   # 15.0
        @test vals["INS"] ≈ samples[end]                      # 0.0 (last)
        @test vals["MAX"] ≈ maximum(samples)                  # 30.0
        @test vals["MIN"] ≈ minimum(samples)                  # 0.0
    end

    @testset "register guards" begin
        tape = CLM.HistoryTape()
        CLM.hist_addfld!(tape, "X1", "u", i -> i.a)
        @test_throws ErrorException CLM.hist_addfld!(tape, "X1", "u", i -> i.a)   # dup
        @test_throws ErrorException CLM.hist_addfld!(tape, "X2", "u", i -> i.a; avgflag="Q")  # bad flag
    end

    @testset "SPVAL / non-finite samples skipped" begin
        inst = _HistFakeInst(a = [0.0])
        tape = CLM.HistoryTape()
        CLM.hist_addfld!(tape, "AVG", "u", i -> i.a; avgflag="A")
        # valid 10, then SPVAL, then NaN, then valid 20 → average over 2 valid = 15
        for s in (10.0, CLM.SPVAL, NaN, 20.0)
            inst.a[1] = s
            CLM.hist_accumulate!(tape, inst)
        end
        f = tape.fields[1]
        @test f.nacs[1] == 2
        @test CLM.hist_field_value(f)[1] ≈ 15.0
    end

    @testset "write + re-read NetCDF round trip" begin
        # 3 spatial elements, multi-field, known per-step values.
        inst = _HistFakeInst(a = zeros(3), b = zeros(3))
        tape = CLM.HistoryTape()
        CLM.hist_addfld!(tape, "FIELD_A", "W/m^2", i -> i.a; avgflag="A",
                         long_name="test average field")
        CLM.hist_addfld!(tape, "FIELD_X", "K", i -> i.b; avgflag="X",
                         long_name="test max field")

        # 4 steps. a[k] increments per step; b varies so max is deterministic.
        a_steps = [[1.0, 2.0, 3.0],
                   [3.0, 4.0, 5.0],
                   [5.0, 6.0, 7.0],
                   [7.0, 8.0, 9.0]]
        b_steps = [[1.0, 9.0, 2.0],
                   [4.0, 1.0, 8.0],
                   [2.0, 3.0, 5.0],
                   [9.0, 2.0, 1.0]]
        for k in 1:4
            inst.a .= a_steps[k]
            inst.b .= b_steps[k]
            CLM.hist_accumulate!(tape, inst)
        end

        expected_a = sum(a_steps) ./ 4
        expected_x = [maximum(getindex.(b_steps, j)) for j in 1:3]

        mktempdir() do dir
            fn = joinpath(dir, "test_h0.nc")
            CLM.hist_write!(tape, fn; time=12.5)

            @test isfile(fn)
            # After write with reset=true, accumulators cleared.
            @test tape.nsamples == 0

            NCDataset(fn, "r") do ds
                @test haskey(ds, "FIELD_A")
                @test haskey(ds, "FIELD_X")
                @test ds.dim["subgrid"] == 3
                @test ds.dim["time"] == 1
                # ds["time"].var = raw stored values (skip CF datetime decoding)
                @test Array(ds["time"].var)[1] ≈ 12.5
                @test ds["FIELD_A"].attrib["units"] == "W/m^2"
                @test ds["FIELD_A"].attrib["cell_methods"] == "time: mean"
                @test ds["FIELD_X"].attrib["cell_methods"] == "time: maximum"

                ra = Array(ds["FIELD_A"])[:, 1]
                rx = Array(ds["FIELD_X"])[:, 1]
                @test ra ≈ expected_a
                @test rx ≈ expected_x
            end
        end
    end

    @testset "append second record" begin
        inst = _HistFakeInst(a = zeros(2))
        tape = CLM.HistoryTape()
        CLM.hist_addfld!(tape, "FA", "u", i -> i.a; avgflag="A")

        mktempdir() do dir
            fn = joinpath(dir, "test_append.nc")

            inst.a .= [2.0, 4.0]; CLM.hist_accumulate!(tape, inst)
            inst.a .= [4.0, 8.0]; CLM.hist_accumulate!(tape, inst)
            CLM.hist_write!(tape, fn; time=1.0, append=true)   # record 1, avg [3,6]

            inst.a .= [10.0, 20.0]; CLM.hist_accumulate!(tape, inst)
            CLM.hist_write!(tape, fn; time=2.0, append=true)   # record 2, avg [10,20]

            NCDataset(fn, "r") do ds
                @test ds.dim["time"] == 2
                @test Array(ds["time"].var) ≈ [1.0, 2.0]
                fa = Array(ds["FA"])
                @test fa[:, 1] ≈ [3.0, 6.0]
                @test fa[:, 2] ≈ [10.0, 20.0]
            end
        end
    end

    @testset "error on empty write" begin
        tape = CLM.HistoryTape()
        @test_throws ErrorException CLM.hist_write!(tape, "nope.nc")  # no fields
        CLM.hist_addfld!(tape, "FA", "u", i -> i.a)
        @test_throws ErrorException CLM.hist_write!(tape, "nope.nc")  # no samples
    end

    @testset "name:flag parsing (getname/getflag)" begin
        @test CLM.hist_getname("TG") == "TG"
        @test CLM.hist_getflag("TG") == ""
        @test CLM.hist_getname("TG:I") == "TG"
        @test CLM.hist_getflag("TG:I") == "I"
        @test CLM.hist_getname(" TG : X ") == "TG"
        @test CLM.hist_getflag(" TG : X ") == "X"
    end

    @testset "set_hist_filename (h0/h1…)" begin
        @test CLM.set_hist_filename("BowCase", 1) == "BowCase.clm2.h0.nc"
        @test CLM.set_hist_filename("BowCase", 2) == "BowCase.clm2.h1.nc"
        @test CLM.set_hist_filename("BowCase", 3; ext_date="2000-01") ==
              "BowCase.clm2.h2.2000-01.nc"
        @test endswith(CLM.set_hist_filename("c", 4; ext_date="2000-01-01-00000"), ".nc")
        @test_throws ErrorException CLM.set_hist_filename("c", 0)
    end

    @testset "multi-tape build: fincl/fexcl + :flag override + nhtfrq" begin
        # A small master list with one on-by-default field and one off-by-default.
        master = CLM.HistFieldSpec[
            CLM.hist_master_field("AA", "u", i -> i.a; avgflag="A"),
            CLM.hist_master_field("BB", "u", i -> i.b; avgflag="A"),
            CLM.hist_master_field("CC", "u", i -> i.a; avgflag="A",
                                  on_by_default=false),
        ]

        # h0: defaults minus BB (fexcl).  h1: defaults + CC via fincl, with
        # tape-wide "X" default but AA overridden to "I". fincl is ADDITIVE
        # in Fortran (it does not replace the on-by-default set unless
        # hist_empty_htapes is set), so BB still appears on h1.
        set = CLM.build_history_tapes(;
            master = master,
            fincl  = [String[],            ["CC", "AA:I"]],
            fexcl  = [["BB"],              String[]],
            nhtfrq = [0,                   -24],
            avgflag_pertape = ["",         "X"],
            case   = "MultiCase",
        )

        @test length(set.tapes) == 2
        @test set.nhtfrq == [0, -24]

        # h0: AA only — BB excluded via fexcl, CC off-by-default & not in fincl.
        h0names = [f.name for f in set.tapes[1].fields]
        @test h0names == ["AA"]
        @test set.tapes[1].fields[1].avgflag == "A"

        # h1: AA (per-field "I") + BB (on-by-default, tape-wide "X") + CC
        # (fincl, tape-wide "X").
        h1 = Dict(f.name => f.avgflag for f in set.tapes[2].fields)
        @test Set(keys(h1)) == Set(["AA", "BB", "CC"])
        @test h1["AA"] == "I"                   # per-field override wins
        @test h1["BB"] == "X"                   # tape-wide pertape flag
        @test h1["CC"] == "X"                   # tape-wide pertape flag

        # hist_empty_htapes=true → tape contains ONLY fincl fields.
        set_empty = CLM.build_history_tapes(;
            master = master,
            fincl  = [["CC:I", "AA"]],
            hist_empty_htapes = true,
        )
        en = Dict(f.name => f.avgflag for f in set_empty.tapes[1].fields)
        @test Set(keys(en)) == Set(["AA", "CC"])   # BB NOT included
        @test en["CC"] == "I"
        @test en["AA"] == "A"                       # falls back to field default

        # Unknown fincl/fexcl name errors.
        @test_throws ErrorException CLM.build_history_tapes(;
            master=master, fincl=[["NOPE"]])
        @test_throws ErrorException CLM.build_history_tapes(;
            master=master, fexcl=[["NOPE"]])
    end

    @testset "multi-tape accumulate + write + re-read h0 & h1" begin
        master = CLM.HistFieldSpec[
            CLM.hist_master_field("AA", "W/m^2", i -> i.a; avgflag="A"),
            CLM.hist_master_field("BB", "K",     i -> i.b; avgflag="A"),
        ]
        # h0: both fields, A-average (defaults). h1: only BB — AA fexcl'd, and
        # BB set instantaneous via the :I fincl override.
        set = CLM.build_history_tapes(;
            master = master,
            fincl  = [String[],   ["BB:I"]],
            fexcl  = [String[],   ["AA"]],
            nhtfrq = [0,          -1],
            case   = "TwoTape",
        )

        inst = _HistFakeInst(a = zeros(2), b = zeros(2))
        a_steps = [[2.0, 4.0], [4.0, 8.0]]
        b_steps = [[1.0, 5.0], [3.0, 9.0]]
        for k in 1:2
            inst.a .= a_steps[k]
            inst.b .= b_steps[k]
            CLM.hist_accumulate!(set, inst)   # accumulates into BOTH tapes
        end

        mktempdir() do dir
            written = CLM.hist_write!(set; dir=dir, time=1.0)
            @test length(written) == 2
            fn_h0 = joinpath(dir, "TwoTape.clm2.h0.nc")
            fn_h1 = joinpath(dir, "TwoTape.clm2.h1.nc")
            @test isfile(fn_h0)
            @test isfile(fn_h1)

            NCDataset(fn_h0, "r") do ds
                @test haskey(ds, "AA")
                @test haskey(ds, "BB")
                @test ds["AA"].attrib["cell_methods"] == "time: mean"
                @test Array(ds["AA"])[:, 1] ≈ [3.0, 6.0]    # avg of a_steps
                @test Array(ds["BB"])[:, 1] ≈ [2.0, 7.0]    # avg of b_steps
            end

            NCDataset(fn_h1, "r") do ds
                @test haskey(ds, "BB")
                @test !haskey(ds, "AA")                     # only BB on h1
                @test ds["BB"].attrib["cell_methods"] == "time: point"  # "I"
                @test Array(ds["BB"])[:, 1] ≈ b_steps[end]  # last sample [3,9]
            end
        end

        # After write, both tapes reset.
        @test all(t -> t.nsamples == 0, set.tapes)
    end

    @testset "default_history_tape via master list" begin
        tape = CLM.default_history_tape()
        names = Set(f.name for f in tape.fields)
        @test "EFLX_LH_TOT" in names
        @test "TG" in names
        @test "H2OSNO" in names
        @test length(tape.fields) == length(CLM.default_master_field_list())
    end
end
