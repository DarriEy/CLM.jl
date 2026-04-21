#!/usr/bin/env julia
# ==========================================================================
# AD (Automatic Differentiation) Smoke Test for CLM.jl
# Tests whether ForwardDiff.Dual numbers flow through physics kernels
# ==========================================================================
using Test
using ForwardDiff
using CLM

println("=" ^ 70)
println("AD SMOKE TEST for CLM.jl")
println("=" ^ 70)

# =========================================================================
# Level 1: Pure Math Functions
# =========================================================================
@testset "Level 1 -- Pure Math" begin

    @testset "Tridiagonal solver" begin
        d = ForwardDiff.derivative(2.0) do x
            n = 5; T = typeof(x)
            a = zeros(T, n); b = fill(T(2), n); c = zeros(T, n)
            r = ones(T, n); u = zeros(T, n)
            b[3] = x
            for i in 2:n; a[i] = T(-0.5); end
            for i in 1:n-1; c[i] = T(-0.5); end
            CLM.tridiagonal_solve!(u, a, b, c, r, 1, n)
            u[3]
        end
        @test isfinite(d)
        println("  Tridiagonal solver: PASS (d=$d)")
    end

    @testset "Tridiagonal multi" begin
        d = ForwardDiff.derivative(2.0) do x
            nc, nl = 2, 4; T = typeof(x)
            a = zeros(T, nc, nl); b = fill(T(2), nc, nl); c = zeros(T, nc, nl)
            r = ones(T, nc, nl); u = zeros(T, nc, nl)
            b[1, 2] = x
            for j in 2:nl; a[:, j] .= T(-0.5); end
            for j in 1:nl-1; c[:, j] .= T(-0.5); end
            CLM.tridiagonal_multi!(u, a, b, c, r, [1,1], BitVector([true,true]), nc, nl)
            u[1, 2]
        end
        @test isfinite(d)
        println("  Tridiagonal multi: PASS (d=$d)")
    end

    @testset "_poly8" begin
        d = ForwardDiff.derivative(x -> CLM._poly8(x, CLM.QSAT_A), 10.0)
        @test isfinite(d)
        println("  _poly8: PASS (d=$d)")
    end

    @testset "qsat" begin
        d = ForwardDiff.derivative(x -> CLM.qsat(x, 101325.0)[1], 280.0)
        @test isfinite(d)
        println("  qsat: PASS (d=$d)")
    end

    @testset "qsat_no_derivs" begin
        d = ForwardDiff.derivative(x -> CLM.qsat_no_derivs(x, 101325.0)[1], 280.0)
        @test isfinite(d)
        println("  qsat_no_derivs: PASS (d=$d)")
    end

    @testset "erf" begin
        d = ForwardDiff.derivative(x -> CLM.erf(x), 1.0)
        @test isfinite(d)
        println("  erf: PASS (d=$d)")
    end

    @testset "daylength" begin
        d = ForwardDiff.derivative(x -> CLM.daylength(0.785, x), 0.3)
        @test isfinite(d)
        println("  daylength: PASS (d=$d)")
    end

    @testset "Band diagonal (pure-Julia fallback for AD)" begin
        d = ForwardDiff.derivative(2.0) do x
            n = 4; nband = 3; T = typeof(x)
            bm = zeros(T, nband, n)
            for j in 1:n; bm[2,j] = T(2); end
            bm[2,2] = x
            for j in 1:n-1; bm[1,j+1] = T(-0.5); bm[3,j] = T(-0.5); end
            r = ones(T, n); u = zeros(T, n)
            CLM.band_diagonal_column!(u, bm, r, 1, n, nband)
            u[2]
        end
        @test isfinite(d)
        println("  Band diagonal: PASS (d=$d)")
    end

    @testset "ft_photo" begin
        d = ForwardDiff.derivative(x -> CLM.ft_photo(x, 60000.0), 300.0)
        @test isfinite(d)
        println("  ft_photo: PASS (d=$d)")
    end

    @testset "fth_photo" begin
        d = ForwardDiff.derivative(x -> CLM.fth_photo(x, 150000.0, 490.0, 1.0), 300.0)
        @test isfinite(d)
        println("  fth_photo: PASS (d=$d)")
    end

    @testset "fth25_photo" begin
        d = ForwardDiff.derivative(x -> CLM.fth25_photo(150000.0, x), 490.0)
        @test isfinite(d)
        println("  fth25_photo: PASS (d=$d)")
    end

    @testset "quadratic_solve" begin
        d = ForwardDiff.derivative(x -> CLM.quadratic_solve(1.0, x, 1.0)[1], -3.0)
        @test isfinite(d)
        println("  quadratic_solve: PASS (d=$d)")
    end

    @testset "combo_scalar" begin
        d = ForwardDiff.derivative(t -> CLM.combo_scalar(0.1, 1.0, 0.5, t, 0.05, 0.5, 0.2, 270.0)[4], 271.0)
        @test isfinite(d)
        println("  combo_scalar: PASS (d=$d)")
    end

    @testset "fresh_snow_radius" begin
        d = ForwardDiff.derivative(t -> CLM.fresh_snow_radius(t), 260.0)
        @test isfinite(d)
        println("  fresh_snow_radius: PASS (d=$d)")
    end

    @testset "liquid_water_heat" begin
        d = ForwardDiff.derivative(t -> CLM.liquid_water_heat(t, 5.0), 280.0)
        @test isfinite(d)
        println("  liquid_water_heat: PASS (d=$d)")
    end

    @testset "overburden_compaction_anderson1976" begin
        d = ForwardDiff.derivative(b -> CLM.overburden_compaction_anderson1976(b, 1.0, -5.0, 250.0), 50.0)
        @test isfinite(d)
        println("  overburden_compaction_anderson1976: PASS (d=$d)")
    end

    @testset "overburden_compaction_vionnet2012" begin
        d = ForwardDiff.derivative(h -> CLM.overburden_compaction_vionnet2012(h, 0.1, 50.0, 1.0, -5.0, 250.0), 1.0)
        @test isfinite(d)
        println("  overburden_compaction_vionnet2012: PASS (d=$d)")
    end

    @testset "mass_weighted_snow_radius" begin
        d = ForwardDiff.derivative(r1 -> CLM.mass_weighted_snow_radius(r1, 100.0, 5.0, 3.0), 80.0)
        @test isfinite(d)
        println("  mass_weighted_snow_radius: PASS (d=$d)")
    end

    @testset "soil_suction_clapp_hornberger" begin
        d = ForwardDiff.derivative(s -> CLM.soil_suction_clapp_hornberger(100.0, s, 5.0), 0.5)
        @test isfinite(d)
        println("  soil_suction_clapp_hornberger: PASS (d=$d)")
    end

    @testset "qsat_water (lake)" begin
        d = ForwardDiff.derivative(t -> CLM.qsat_water(t, 101325.0), 280.0)
        @test isfinite(d)
        println("  qsat_water: PASS (d=$d)")
    end

    @testset "qsat_water_dT (lake)" begin
        d = ForwardDiff.derivative(t -> CLM.qsat_water_dT(t, 101325.0), 280.0)
        @test isfinite(d)
        println("  qsat_water_dT: PASS (d=$d)")
    end

    @testset "piecewise_linear_interp1d" begin
        xd = [1.0, 2.0, 3.0, 4.0]
        yd = [10.0, 20.0, 15.0, 25.0]
        d = ForwardDiff.derivative(x -> CLM.piecewise_linear_interp1d(xd, yd, x), 2.5)
        @test isfinite(d)
        println("  piecewise_linear_interp1d: PASS (d=$d)")
    end
end

# =========================================================================
# Level 2: Photosynthesis solver chain (P1 — FIXED)
# =========================================================================
println()
println("=" ^ 70)
println("Level 2: Photosynthesis & LUNA (widened to ::Real)")
println("=" ^ 70)

@testset "Level 2 -- Photosynthesis & LUNA" begin

    @testset "LUNA vcmx_t_kattge" begin
        d = ForwardDiff.derivative(t -> CLM.vcmx_t_kattge(20.0, t), 25.0)
        @test isfinite(d)
        println("  vcmx_t_kattge: PASS (d=$d)")
    end

    @testset "LUNA jmx_t_kattge" begin
        d = ForwardDiff.derivative(t -> CLM.jmx_t_kattge(20.0, t), 25.0)
        @test isfinite(d)
        println("  jmx_t_kattge: PASS (d=$d)")
    end

    @testset "LUNA resp_t_bernacchi" begin
        d = ForwardDiff.derivative(t -> CLM.resp_t_bernacchi(t), 25.0)
        @test isfinite(d)
        println("  resp_t_bernacchi: PASS (d=$d)")
    end

    @testset "LUNA quadratic_luna" begin
        d = ForwardDiff.derivative(x -> CLM.quadratic_luna(1.0, x, 1.0)[1], -3.0)
        @test isfinite(d)
        println("  quadratic_luna: PASS (d=$d)")
    end

    @testset "plc (vulnerability curve)" begin
        # plc(x, ivt, level, plc_method, params) — x is Real, rest are Int/struct
        params = CLM.PhotoParamsData()
        CLM.photo_params_init!(params)
        params.psi50[2, 1] = -200000.0
        params.ck[2, 1] = 2.0
        d = ForwardDiff.derivative(s -> CLM.plc(s, 2, 1, CLM.VEGETATION_WEIBULL, params), -50000.0)
        @test isfinite(d)
        println("  plc: PASS (d=$d)")
    end

    @testset "d1plc (vulnerability curve derivative)" begin
        params = CLM.PhotoParamsData()
        CLM.photo_params_init!(params)
        params.psi50[2, 1] = -200000.0
        params.ck[2, 1] = 2.0
        d = ForwardDiff.derivative(s -> CLM.d1plc(s, 2, 1, CLM.VEGETATION_WEIBULL, params), -50000.0)
        @test isfinite(d)
        println("  d1plc: PASS (d=$d)")
    end

end

# =========================================================================
# AD Status Summary
# =========================================================================
println()
println("=" ^ 70)
println("AD Coverage Summary")
println("=" ^ 70)
println("""
COMPLETED:
  P0: Struct init functions — FT-generic (zero(FT), fill(FT(x), ...))
  P1: Photosynthesis solver chain — all 15 functions widened to ::Real
  P2: Band diagonal solver — pure-Julia fallback for non-Float64 types
  P3: LUNA module — all functions widened to ::Real, LunaParamsData{FT}

REMAINING (low priority):
  P4: Utility functions (wind_drift_compaction!, dry_dep, irrigation helpers)
  P5: Full end-to-end AD through driver (needs struct arrays to store Dual)
""")
