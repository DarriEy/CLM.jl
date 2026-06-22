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
                @test ds.dim["patch"] == 3   # default subgrid level
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
        @test "TSOI" in names                       # 2d column field
        @test length(tape.fields) == length(CLM.default_master_field_list())
        # native subgrid levels are tagged correctly
        bylvl = Dict(f.name => f.level for f in tape.fields)
        @test bylvl["EFLX_LH_TOT"] == CLM.HIST_PATCH
        @test bylvl["TG"] == CLM.HIST_COLUMN
        tsoi = tape.fields[findfirst(f -> f.name == "TSOI", tape.fields)]
        @test tsoi.level == CLM.HIST_COLUMN
        @test tsoi.nlev > 1                          # 2d (soil-layer) field
        @test tsoi.levdim == "levsoi"
    end

    # ------------------------------------------------------------------
    # Subgrid levels + 2D fields (Tier B: hist_dov2xy / hist_type1d_pertape
    # + hist_addfld2d). A richer fake `inst` carries patch / column /
    # landunit subgrid index+weight maps (for gridcell remap).
    # ------------------------------------------------------------------
    Base.@kwdef mutable struct _HistSubInst
        p_vals::Vector{Float64} = Float64[]   # patch level (np)
        c_vals::Vector{Float64} = Float64[]   # column level (nc)
        g_vals::Vector{Float64} = Float64[]   # gridcell level (ng)
        c2d::Matrix{Float64}    = Matrix{Float64}(undef, 0, 0)  # 2D column field (nc,nlev)
        patch    = nothing
        column   = nothing
        landunit = nothing
    end
    Base.@kwdef mutable struct _HistMap
        gridcell::Vector{Int}     = Int[]
        wtgcell::Vector{Float64}  = Float64[]
    end

    @testset "mixed subgrid levels coexist (no DimensionMismatch)" begin
        # Bow-like: np=4, nc=2, ng=1 — different lengths previously crashed
        # default_history_tape() with DimensionMismatch.
        inst = _HistSubInst(
            p_vals = zeros(4), c_vals = zeros(2), g_vals = zeros(1),
            patch  = _HistMap(gridcell=[1,1,1,1], wtgcell=fill(0.25, 4)),
            column = _HistMap(gridcell=[1,1],     wtgcell=[0.6, 0.4]),
            landunit = _HistMap(gridcell=Int[],   wtgcell=Float64[]),
        )
        tape = CLM.HistoryTape()
        CLM.hist_addfld!(tape, "PFLD", "u", i -> i.p_vals; level="patch")
        CLM.hist_addfld!(tape, "CFLD", "u", i -> i.c_vals; level="column")
        CLM.hist_addfld!(tape, "GFLD", "u", i -> i.g_vals; level="gridcell")

        p_steps = [[1.0,2.0,3.0,4.0], [2.0,3.0,4.0,5.0], [3.0,4.0,5.0,6.0]]
        c_steps = [[10.0,20.0], [20.0,30.0], [30.0,40.0]]
        g_steps = [[100.0], [200.0], [300.0]]
        for k in 1:3
            inst.p_vals .= p_steps[k]; inst.c_vals .= c_steps[k]; inst.g_vals .= g_steps[k]
            CLM.hist_accumulate!(tape, inst)   # must NOT raise DimensionMismatch
        end
        @test tape.nsamples == 3
        exp_p = sum(p_steps) ./ 3; exp_c = sum(c_steps) ./ 3; exp_g = sum(g_steps) ./ 3

        mktempdir() do dir
            fn = joinpath(dir, "mixed.nc")
            CLM.hist_write!(tape, fn; time=1.0)   # native-level write
            NCDataset(fn, "r") do ds
                @test ds.dim["patch"]    == 4   # separate dim per level
                @test ds.dim["column"]   == 2
                @test ds.dim["gridcell"] == 1
                @test Array(ds["PFLD"])[:, 1] ≈ exp_p
                @test Array(ds["CFLD"])[:, 1] ≈ exp_c
                @test Array(ds["GFLD"])[:, 1] ≈ exp_g
                @test ds["PFLD"].attrib["subgrid_level"] == "patch"
                @test ds["CFLD"].attrib["subgrid_level"] == "column"
            end
        end
    end

    @testset "2D column field (TSOI-style) round trip" begin
        nc, nlev = 2, 3
        inst = _HistSubInst(c2d = zeros(nc, nlev))
        tape = CLM.HistoryTape()
        CLM.hist_addfld!(tape, "TSOI", "K", i -> i.c2d;
                         level="column", nlev=nlev, levdim="levsoi",
                         long_name="soil temperature")

        s1 = [280.0 281.0 282.0; 270.0 271.0 272.0]
        s2 = [282.0 283.0 284.0; 272.0 273.0 274.0]
        for s in (s1, s2)
            inst.c2d .= s
            CLM.hist_accumulate!(tape, inst)
        end
        expected = (s1 .+ s2) ./ 2

        fv = CLM.hist_field_value(tape.fields[1])
        @test size(fv) == (nc, nlev)          # finalize -> (nc, nlev) matrix
        @test fv ≈ expected

        mktempdir() do dir
            fn = joinpath(dir, "tsoi.nc")
            CLM.hist_write!(tape, fn; time=5.0)
            NCDataset(fn, "r") do ds
                @test ds.dim["column"] == nc
                @test ds.dim["levsoi"] == nlev
                @test ds.dim["time"]   == 1
                rd = Array(ds["TSOI"])         # (column, levsoi, time)
                @test size(rd) == (nc, nlev, 1)
                @test rd[:, :, 1] ≈ expected
            end
        end
    end

    @testset "mixed 1D + 2D fields in one tape" begin
        nc, nlev = 2, 3
        inst = _HistSubInst(p_vals = zeros(4), c_vals = zeros(2), c2d = zeros(nc, nlev))
        tape = CLM.HistoryTape()
        CLM.hist_addfld!(tape, "PFLD", "u", i -> i.p_vals; level="patch")
        CLM.hist_addfld!(tape, "CFLD", "u", i -> i.c_vals; level="column")
        CLM.hist_addfld!(tape, "TSOI", "K", i -> i.c2d;
                         level="column", nlev=nlev, levdim="levsoi")
        inst.p_vals .= [1.0,2.0,3.0,4.0]
        inst.c_vals .= [10.0,20.0]
        inst.c2d    .= [280.0 281.0 282.0; 270.0 271.0 272.0]
        CLM.hist_accumulate!(tape, inst)
        @test tape.nsamples == 1
        mktempdir() do dir
            fn = joinpath(dir, "mixed12.nc")
            CLM.hist_write!(tape, fn; time=0.0)
            NCDataset(fn, "r") do ds
                @test ds.dim["patch"]  == 4
                @test ds.dim["column"] == 2
                @test ds.dim["levsoi"] == nlev
                @test Array(ds["PFLD"])[:, 1]   ≈ [1.0,2.0,3.0,4.0]
                @test Array(ds["CFLD"])[:, 1]   ≈ [10.0,20.0]
                @test Array(ds["TSOI"])[:, :, 1] ≈ inst.c2d
            end
        end
    end

    @testset "gridcell remap (dov2xy) area-weighted aggregate" begin
        # ng=2 gridcells. patch->gridcell and column->gridcell weights.
        inst = _HistSubInst(
            p_vals = zeros(4), c_vals = zeros(3),
            patch  = _HistMap(gridcell=[1,1,2,2], wtgcell=[0.3,0.7,0.5,0.5]),
            column = _HistMap(gridcell=[1,2,2],   wtgcell=[1.0,0.4,0.6]),
            landunit = _HistMap(gridcell=Int[],   wtgcell=Float64[]),
        )
        tape = CLM.HistoryTape(dov2xy=true)
        CLM.hist_addfld!(tape, "PFLD", "u", i -> i.p_vals; level="patch")
        CLM.hist_addfld!(tape, "CFLD", "u", i -> i.c_vals; level="column")

        pv = [10.0, 20.0, 30.0, 40.0]
        cv = [100.0, 200.0, 300.0]
        inst.p_vals .= pv; inst.c_vals .= cv
        CLM.hist_accumulate!(tape, inst)

        # patch g1 = (10*.3+20*.7)/(.3+.7)=17 ; g2 = (30*.5+40*.5)/1=35
        # col   g1 = 100*1/1=100 ; g2 = (200*.4+300*.6)/1=260
        exp_pg = [17.0, 35.0]
        exp_cg = [100.0, 260.0]

        ng = CLM._hist_ngridcell(inst)
        @test ng == 2
        pf = CLM.hist_field_value(tape.fields[1])
        cf = CLM.hist_field_value(tape.fields[2])
        @test CLM.hist_dov2xy(tape.fields[1], pf, ng, inst) ≈ exp_pg
        @test CLM.hist_dov2xy(tape.fields[2], cf, ng, inst) ≈ exp_cg

        mktempdir() do dir
            fn = joinpath(dir, "v2xy.nc")
            CLM.hist_write!(tape, fn; time=1.0, inst=inst)
            NCDataset(fn, "r") do ds
                @test ds.dim["gridcell"] == 2         # all fields share one dim
                @test !haskey(ds.dim, "patch")
                @test !haskey(ds.dim, "column")
                @test Array(ds["PFLD"])[:, 1] ≈ exp_pg
                @test Array(ds["CFLD"])[:, 1] ≈ exp_cg
                @test ds["PFLD"].attrib["subgrid_level"] == "gridcell"
            end
        end

        # dov2xy without inst must error.
        tape2 = CLM.HistoryTape(dov2xy=true)
        CLM.hist_addfld!(tape2, "PFLD", "u", i -> i.p_vals; level="patch")
        inst.p_vals .= pv
        CLM.hist_accumulate!(tape2, inst)
        @test_throws ErrorException CLM.hist_write!(tape2, "nope.nc")
    end

    @testset "2D field gridcell remap (dov2xy) per level" begin
        # nc=2 columns -> ng=1 gridcell, nlev=2.
        nlev = 2
        inst = _HistSubInst(
            c2d = zeros(2, nlev),
            column = _HistMap(gridcell=[1,1], wtgcell=[0.25, 0.75]),
            patch  = _HistMap(gridcell=Int[], wtgcell=Float64[]),
            landunit = _HistMap(gridcell=Int[], wtgcell=Float64[]),
        )
        tape = CLM.HistoryTape(dov2xy=true)
        CLM.hist_addfld!(tape, "TSOI", "K", i -> i.c2d;
                         level="column", nlev=nlev, levdim="levsoi")
        inst.c2d .= [280.0 300.0; 284.0 320.0]   # (col, lev)
        CLM.hist_accumulate!(tape, inst)

        # g1 lev1 = 280*.25+284*.75 = 283 ; lev2 = 300*.25+320*.75 = 315
        exp = reshape([283.0, 315.0], 1, 2)
        ng = CLM._hist_ngridcell(inst)
        fv = CLM.hist_field_value(tape.fields[1])
        @test CLM.hist_dov2xy(tape.fields[1], fv, ng, inst) ≈ exp
        mktempdir() do dir
            fn = joinpath(dir, "tsoi_v2xy.nc")
            CLM.hist_write!(tape, fn; time=0.0, inst=inst)
            NCDataset(fn, "r") do ds
                @test ds.dim["gridcell"] == 1
                @test ds.dim["levsoi"]   == nlev
                @test Array(ds["TSOI"])[:, :, 1] ≈ exp
            end
        end
    end

    @testset "bad subgrid level rejected" begin
        tape = CLM.HistoryTape()
        @test_throws ErrorException CLM.hist_addfld!(tape, "Z", "u", i -> i.a; level="bogus")
    end
end
