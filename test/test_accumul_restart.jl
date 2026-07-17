@testset "Accumulator restart restore" begin
    using NCDatasets

    np = 3

    # Build a minimal atm2lnd with the runmean accumulator patch vectors allocated
    # (as InitAccBuffer does under use_cn) and seeded to a cold value.
    a2l = CLM.Atm2LndData{Float64}()
    a2l.rh30_patch   = fill(0.0, np)
    a2l.prec10_patch = fill(0.0, np)
    a2l.prec30_patch = fill(0.0, np)
    a2l.prec60_patch = fill(0.0, np)

    # A dump carrying the persisted VALUE + NSTEPS (dtime=3600 -> 30d window = 720 steps).
    tmp = tempname() * ".nc"
    NCDataset(tmp, "c") do ds
        defDim(ds, "pft", np)
        for (name, val) in (("RH30_VALUE", 67.63), ("PREC10_VALUE", 1.89e-5),
                            ("PREC30_VALUE", 1.90e-5), ("PREC60_VALUE", 1.87e-5))
            v = defVar(ds, name, Float64, ("pft",)); v[:] = fill(val, np)
        end
        for (name, ns) in (("RH30_NSTEPS", 720), ("PREC10_NSTEPS", 240),
                           ("PREC30_NSTEPS", 720), ("PREC60_NSTEPS", 1440))
            v = defVar(ds, name, Int32, ("pft",)); v[:] = fill(Int32(ns), np)
        end
    end

    bounds = (; endp = np)
    r = CLM.restore_atm2lnd_runmean_accum!(a2l, tmp, bounds; dtime = 3600)

    # VALUES restored into the patch vectors.
    @test all(a2l.rh30_patch   .== 67.63)
    @test all(a2l.prec10_patch .== 1.89e-5)
    @test all(a2l.prec30_patch .== 1.90e-5)
    @test all(a2l.prec60_patch .== 1.87e-5)

    # Reported restored fields + warm counter.
    @test Set(r.restored) == Set([:rh30_patch, :prec10_patch, :prec30_patch, :prec60_patch])
    @test r.nsteps[:prec60_patch] == 1440
    @test r.warm_nstep == 1440   # max persisted nsteps == the 60-day window in steps

    # A cold dump (no *_VALUE) leaves every field untouched (default byte-identical).
    a2l2 = CLM.Atm2LndData{Float64}()
    a2l2.rh30_patch = fill(-9.0, np)
    tmp2 = tempname() * ".nc"
    NCDataset(tmp2, "c") do ds
        defDim(ds, "pft", np)
        v = defVar(ds, "SOMETHING_ELSE", Float64, ("pft",)); v[:] = fill(1.0, np)
    end
    r2 = CLM.restore_atm2lnd_runmean_accum!(a2l2, tmp2, (; endp = np); dtime = 3600)
    @test isempty(r2.restored)
    @test r2.warm_nstep == 0
    @test all(a2l2.rh30_patch .== -9.0)

    # An unallocated target field (gate off) is skipped, not an error.
    a2l3 = CLM.Atm2LndData{Float64}()   # rh30_patch etc. stay empty
    r3 = CLM.restore_atm2lnd_runmean_accum!(a2l3, tmp, (; endp = np); dtime = 3600)
    @test isempty(r3.restored)

    rm(tmp; force = true); rm(tmp2; force = true)
end
