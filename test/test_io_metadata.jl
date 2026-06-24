# Tests for the CTSM history + restart NetCDF metadata added to
# src/infrastructure/history_io.jl + src/infrastructure/restart_io.jl:
#   - history time-bounds vars (time_bounds / mcdate / mcsec / mdcur / mscur)
#   - history mfilt rollover (multiple time slices per file, then a new file)
#   - history ndens output precision (single vs double per tape)
#   - restart global-attribute header (title/case/version/provenance)
# Mirrors the existing test_history_io.jl / test_restart_io.jl style.

using Test
using NCDatasets

Base.@kwdef mutable struct _MetaFakeInst
    a::Vector{Float64} = Float64[]
end

@testset "io_metadata" begin

    @testset "history time_bounds + mcdate/mcsec for an averaged field" begin
        inst = _MetaFakeInst(a = zeros(2))
        tape = CLM.HistoryTape()
        CLM.hist_addfld!(tape, "AVG", "u", i -> i.a; avgflag="A")

        # Accumulate two samples → time-mean field.
        inst.a .= [10.0, 20.0]; CLM.hist_accumulate!(tape, inst)
        inst.a .= [30.0, 40.0]; CLM.hist_accumulate!(tape, inst)

        mktempdir() do dir
            fn = joinpath(dir, "h0_bounds.nc")
            # Interval [0, 2) days since 2000-01-01 → end date 2000-01-03.
            CLM.hist_write!(tape, fn; time=2.0, begtime=0.0,
                            time_units="days since 2000-01-01")

            NCDataset(fn, "r") do ds
                # bounds wiring on the time coordinate
                @test haskey(ds, "time_bounds")
                @test ds["time"].attrib["bounds"] == "time_bounds"
                @test ds.dim["hist_interval"] == 2

                # .var: raw endpoints (skip CF datetime decoding of the bounds)
                tb = Array(ds["time_bounds"].var)     # (hist_interval, time)
                @test size(tb) == (2, 1)
                @test tb[1, 1] ≈ 0.0    # interval start (begtime)
                @test tb[2, 1] ≈ 2.0    # interval end   (time)
                @test ds["time_bounds"].attrib["long_name"] ==
                      "history time interval endpoints"

                # calendar bookkeeping vars present + correct
                for nm in ("mcdate", "mcsec", "mdcur", "mscur")
                    @test haskey(ds, nm)
                end
                # 2 days since 2000-01-01 00:00 → 2000-01-03, 0 sec of day
                # (.var bypasses CF fill-value masking for the Int slices)
                @test Array(ds["mcdate"].var)[1] == 20000103
                @test Array(ds["mcsec"].var)[1]  == 0
                @test Array(ds["mdcur"].var)[1]  == 2
                @test Array(ds["mscur"].var)[1]  == 0

                # the averaged field is still correct
                @test Array(ds["AVG"])[:, 1] ≈ [20.0, 30.0]
            end
        end
    end

    @testset "mcsec for a sub-day interval end" begin
        inst = _MetaFakeInst(a = zeros(1))
        tape = CLM.HistoryTape()
        CLM.hist_addfld!(tape, "F", "u", i -> i.a)
        inst.a .= [1.0]; CLM.hist_accumulate!(tape, inst)
        mktempdir() do dir
            fn = joinpath(dir, "h0_subday.nc")
            # 1.5 days → 2000-01-02 12:00:00 → mcsec = 12*3600
            CLM.hist_write!(tape, fn; time=1.5, begtime=1.0,
                            time_units="days since 2000-01-01")
            NCDataset(fn, "r") do ds
                @test Array(ds["mcdate"].var)[1] == 20000102
                @test Array(ds["mcsec"].var)[1]  == 12 * 3600
                @test Array(ds["mscur"].var)[1]  == 12 * 3600
                @test Array(ds["mdcur"].var)[1]  == 1
            end
        end
    end

    @testset "mfilt rollover: N slices per file then a new file" begin
        # mfilt=2 → two records land in h0, the third opens a new file.
        master = [CLM.hist_master_field("AVG", "u", i -> i.a; level="patch")]
        set = CLM.build_history_tapes(master=master, hist_mfilt=[2],
                                      case="rolltest")

        mktempdir() do dir
            written = String[]
            for k in 1:3
                inst = (a = [Float64(k)],)
                CLM.hist_accumulate!(set, inst)
                w = CLM.hist_write!(set; dir=dir, time=Float64(k))
                append!(written, w)
            end

            files = unique(written)
            # Three writes, mfilt=2 → first two share one file, third a new one.
            @test length(files) == 2
            # First file holds 2 time records, second holds 1.
            NCDataset(files[1], "r") do ds
                @test ds.dim["time"] == 2
                @test Array(ds["AVG"])[1, :] ≈ [1.0, 2.0]
            end
            NCDataset(files[2], "r") do ds
                @test ds.dim["time"] == 1
                @test Array(ds["AVG"])[1, 1] ≈ 3.0
            end
            # The rolled-over file carries the file-index suffix.
            @test occursin("rolltest.clm2.h0.nc", basename(files[1]))
            @test occursin("h0-00001.nc", basename(files[2]))
        end
    end

    @testset "ndens output precision: single vs double" begin
        a_steps = [[1.0, 2.0], [3.0, 4.0]]

        function write_with_ndens(ndens, fn)
            inst = _MetaFakeInst(a = zeros(2))
            tape = CLM.HistoryTape(ndens=ndens)
            CLM.hist_addfld!(tape, "FA", "u", i -> i.a; avgflag="A")
            for s in a_steps
                inst.a .= s; CLM.hist_accumulate!(tape, inst)
            end
            CLM.hist_write!(tape, fn; time=1.0)
        end

        mktempdir() do dir
            fdouble = joinpath(dir, "h0_double.nc")
            fsingle = joinpath(dir, "h0_single.nc")
            write_with_ndens(1, fdouble)   # 1 = double
            write_with_ndens(2, fsingle)   # 2 = single

            NCDataset(fdouble, "r") do ds
                @test eltype(ds["FA"].var) == Float64
            end
            NCDataset(fsingle, "r") do ds
                @test eltype(ds["FA"].var) == Float32
            end
        end
    end

    @testset "restart global-attribute header" begin
        ng, nl, nc, np = 2, 3, 6, 12
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
        inst = CLM.CLMInstances()
        CLM.clm_instInit!(inst; ng=ng, nl=nl, nc=nc, np=np)
        bounds = CLM.BoundsType(begc=1, endc=nc, begp=1, endp=np,
                                begg=1, endg=ng, begl=1, endl=nl)

        ctl = CLM.VarCtl(caseid="mycase", ctitle="My Test Case",
                         version="CLM.jl-0.1", source="CTSM-port",
                         conventions="CF-1.8", username="tester",
                         hostname="testhost", fsurdat="/data/surface.nc")

        mktempdir() do dir
            file = joinpath(dir, "clm_restart_hdr.nc")
            CLM.write_restart(inst, file; bounds=bounds, ctl=ctl)
            @test isfile(file)

            NCDataset(file, "r") do ds
                ga = ds.attrib
                @test ga["title"]       == "CLM Restart information"
                @test ga["case_id"]     == "mycase"
                @test ga["case_title"]  == "My Test Case"
                @test ga["version"]     == "CLM.jl-0.1"
                @test ga["source"]      == "CTSM-port"
                @test ga["Conventions"] == "CF-1.8"
                @test ga["username"]    == "tester"
                @test ga["host"]        == "testhost"
                @test ga["surface_dataset"] == "/data/surface.nc"
                # provenance "history" attr present + creation-stamp format
                @test haskey(ga, "history")
                @test occursin("created on", ga["history"])
                # boolean flag metadata serialized as 'true'/'false'
                @test ga["irrigate"] in ("true", "false")
                @test ga["created_glacier_mec_landunits"] == "true"
            end
        end
    end
end
