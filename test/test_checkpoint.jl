# Tests for the complete bit-exact-resumable checkpoint (checkpoint.jl).
# CI-safe: builds an in-memory instance via clm_instInit! (no external data).

using Dates

@testset "Checkpoint (complete, bit-exact resume)" begin
    ng, nl, nc, np = 2, 3, 6, 12
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

    inst = CLM.CLMInstances()
    CLM.clm_instInit!(inst; ng=ng, nl=nl, nc=nc, np=np)

    # Seed distinct finite state across several sub-instances (incl. ones the curated
    # NetCDF restart omits — accumulators, snow grain radius — which the checkpoint must
    # still round-trip).
    T = inst.temperature
    T.t_grnd_col   .= 270.0 .+ (1:nc)
    T.t_a10_patch  .= 280.0 .+ (1:np)        # accumulator (not in NetCDF restart registry)
    T.t_veg240_patch .= 285.0 .+ (1:np)
    wd = inst.water.waterdiagnosticbulk_inst
    wd.snw_rds_col .= 100.0 .+ reshape(1:length(wd.snw_rds_col), size(wd.snw_rds_col))
    ws = inst.water.waterstatebulk_inst.ws
    ws.h2osoi_liq_col .= 10.0 .+ reshape(1:length(ws.h2osoi_liq_col), size(ws.h2osoi_liq_col))

    tm = CLM.TimeManager(; nstep=137, current_date=DateTime(2003, 7, 15, 3),
                          start_date=DateTime(2003, 1, 1), dtime=1800)

    dir = mktempdir()
    path = joinpath(dir, "state.jls")

    @testset "round-trips the complete state + cursor" begin
        CLM.write_checkpoint(path, inst, tm)
        @test isfile(path)
        (inst2, nstep2, date2) = CLM.read_checkpoint(path)

        @test nstep2 == 137
        @test date2 == DateTime(2003, 7, 15, 3)
        # prognostic field the NetCDF restart also carries
        @test inst2.temperature.t_grnd_col == T.t_grnd_col
        @test inst2.water.waterstatebulk_inst.ws.h2osoi_liq_col == ws.h2osoi_liq_col
        # forward-feeding state the NetCDF restart OMITS — the checkpoint must carry it
        @test inst2.temperature.t_a10_patch == T.t_a10_patch
        @test inst2.temperature.t_veg240_patch == T.t_veg240_patch
        @test inst2.water.waterdiagnosticbulk_inst.snw_rds_col == wd.snw_rds_col
    end

    @testset "restored instance is independent (deep copy)" begin
        (inst2, _, _) = CLM.read_checkpoint(path)
        inst2.temperature.t_grnd_col[1] = -999.0
        @test T.t_grnd_col[1] != -999.0          # mutating the restored copy doesn't alias
    end

    @testset "a corrupt / non-checkpoint file is rejected" begin
        badpath = joinpath(dir, "bad.jls")
        write(badpath, "not a checkpoint")
        @test_throws Exception CLM.read_checkpoint(badpath)
    end
end
