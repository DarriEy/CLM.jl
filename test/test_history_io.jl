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
end
