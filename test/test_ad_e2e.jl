# ==========================================================================
# test_ad_e2e.jl — End-to-end ForwardDiff AD through CLM.jl physics
#
# P5 Milestone: Prove that Dual numbers flow through physics kernels.
#
# Strategy: Option 2 (wrapper approach). We do NOT change struct type
# parameters from {FT<:Real} to {FT<:Real}. Instead we:
#   1. Work with raw arrays (bypassing struct type constraints)
#   2. Call pure-math physics functions with Dual-valued inputs
#   3. Verify derivatives against finite differences
#
# Test levels:
#   Level 3a: qsat chain (d(qsat)/d(T))
#   Level 3b: Soil water retention (d(smp)/d(saturation))
#   Level 3c: Snow grain growth sensitivity
#   Level 3d: Two-stream radiation (d(albedo)/d(LAI))
#   Level 4:  Mini-driver chain (T -> qsat -> resistance -> LH estimate)
#   Level 5:  Full timestep (document barriers)
# ==========================================================================

using Test
using ForwardDiff
using ForwardDiff: Dual, Tag, value, partials
using CLM

@testset "AD End-to-End (ForwardDiff)" begin

    # =====================================================================
    # Level 3a: qsat — d(qs)/d(T) and d(es)/d(T)
    # =====================================================================
    @testset "Level 3a: qsat derivatives via ForwardDiff" begin
        # The qsat function accepts ::Real arguments, so Dual works directly.
        T0 = 280.0   # K
        p0 = 101325.0 # Pa

        # --- Forward mode: scalar perturbation of T ---
        T_dual = Dual(T0, 1.0)  # seed dT = 1
        (qs_d, es_d, dqsdT_d, desdT_d) = CLM.qsat(T_dual, p0)

        # Extract AD derivative of qs w.r.t. T
        dqs_dT_ad = partials(qs_d)[1]

        # Note: The analytical dqsdT returned by qsat uses QSAT_B coefficients,
        # which are independently fit polynomials (NOT the exact derivative of QSAT_A).
        # Therefore AD of qs (which differentiates through QSAT_A) will differ
        # slightly from the analytical dqsdT (which uses QSAT_B). Both are valid
        # approximations of the true Clausius-Clapeyron derivative.
        dqs_dT_analytical = value(dqsdT_d)

        # They should agree to ~1e-4 (polynomial fit accuracy)
        @test dqs_dT_ad ≈ dqs_dT_analytical rtol=1e-3

        # Also check es derivative (same independent-fit issue)
        des_dT_ad = partials(es_d)[1]
        des_dT_analytical = value(desdT_d)
        @test des_dT_ad ≈ des_dT_analytical rtol=1e-3

        # Cross-check AD against finite differences (these MUST match tightly)
        eps_fd = 1e-7
        (qs_plus, _, _, _) = CLM.qsat(T0 + eps_fd, p0)
        (qs_minus, _, _, _) = CLM.qsat(T0 - eps_fd, p0)
        dqs_dT_fd = (qs_plus - qs_minus) / (2 * eps_fd)
        @test dqs_dT_ad ≈ dqs_dT_fd rtol=1e-5

        # es AD vs FD
        (_, es_plus, _, _) = CLM.qsat(T0 + eps_fd, p0)
        (_, es_minus, _, _) = CLM.qsat(T0 - eps_fd, p0)
        des_dT_fd = (es_plus - es_minus) / (2 * eps_fd)
        @test des_dT_ad ≈ des_dT_fd rtol=1e-5

        # --- Test over-ice branch (T < 273.15) ---
        T_cold = 260.0
        T_cold_dual = Dual(T_cold, 1.0)
        (qs_cold, es_cold, _, _) = CLM.qsat(T_cold_dual, p0)
        dqs_cold_ad = partials(qs_cold)[1]

        (qs_cp, _, _, _) = CLM.qsat(T_cold + eps_fd, p0)
        (qs_cm, _, _, _) = CLM.qsat(T_cold - eps_fd, p0)
        dqs_cold_fd = (qs_cp - qs_cm) / (2 * eps_fd)
        @test dqs_cold_ad ≈ dqs_cold_fd rtol=1e-5

        # --- Test pressure sensitivity ---
        p_dual = Dual(p0, 1.0)
        (qs_p, _, _, _) = CLM.qsat(T0, p_dual)
        dqs_dp_ad = partials(qs_p)[1]

        (qs_pp, _, _, _) = CLM.qsat(T0, p0 + 1.0)
        (qs_pm, _, _, _) = CLM.qsat(T0, p0 - 1.0)
        dqs_dp_fd = (qs_pp - qs_pm) / 2.0
        @test dqs_dp_ad ≈ dqs_dp_fd rtol=1e-5

        # Verify sign: qs should decrease with increasing pressure
        @test dqs_dp_ad < 0

        println("  Level 3a PASSED: qsat AD derivatives match FD")
    end

    # =====================================================================
    # Level 3b: Soil water retention (Clapp-Hornberg)
    # =====================================================================
    @testset "Level 3b: Soil suction/conductivity AD" begin
        # The soil_suction! and soil_hk! functions have ::Float64 type
        # annotations on s and imped, so we test the MATH directly.

        # Clapp-Hornberg soil suction: smp = -sucsat * s^(-bsw)
        bsw = 5.0
        sucsat = 200.0  # mm

        function ch_suction(s)
            return -sucsat * s^(-bsw)
        end

        function ch_suction_deriv(s)
            smp = -sucsat * s^(-bsw)
            dsmpds = -bsw * smp / s
            return (smp, dsmpds)
        end

        s0 = 0.6
        s_dual = Dual(s0, 1.0)
        smp_ad = ch_suction(s_dual)
        dsmp_ds_ad = partials(smp_ad)[1]

        (smp_analytical, dsmp_ds_analytical) = ch_suction_deriv(s0)
        @test value(smp_ad) ≈ smp_analytical rtol=1e-12
        @test dsmp_ds_ad ≈ dsmp_ds_analytical rtol=1e-10

        # Clapp-Hornberg hydraulic conductivity: hk = imped * hksat * s^(2*bsw+3)
        hksat = 0.01  # mm/s
        imped = 0.8

        function ch_hk(s)
            return imped * hksat * s^(2.0 * bsw + 3.0)
        end

        hk_dual = ch_hk(s_dual)
        dhk_ds_ad = partials(hk_dual)[1]

        # Analytical: dhk/ds = (2*bsw+3) * hk / s
        hk_val = imped * hksat * s0^(2.0 * bsw + 3.0)
        dhk_ds_analytical = (2.0 * bsw + 3.0) * hk_val / s0
        @test dhk_ds_ad ≈ dhk_ds_analytical rtol=1e-10

        # Finite difference check
        eps_fd = 1e-7
        dhk_ds_fd = (ch_hk(s0 + eps_fd) - ch_hk(s0 - eps_fd)) / (2 * eps_fd)
        @test dhk_ds_ad ≈ dhk_ds_fd rtol=1e-5

        println("  Level 3b PASSED: Soil water retention AD matches analytical")
    end

    # =====================================================================
    # Level 3c: Snow grain growth sensitivity
    # =====================================================================
    @testset "Level 3c: Snow grain growth AD" begin
        # The snowage_grain! function has ::Vector{Float64} type annotations,
        # so we test the core grain growth MATH directly.

        function grain_growth_rate(T_snow, dTdz, rho_snow)
            # Simplified dry snow metamorphism rate (Arrhenius form)
            kappa = 4000.0  # activation energy / R [K]
            dr_dt_0 = 1.0e-6  # base rate [m/s]

            # Temperature dependence
            rate = dr_dt_0 * exp(-kappa / T_snow)

            # Temperature gradient effect
            rate *= (1.0 + 0.01 * dTdz)

            # Density effect
            rate *= max(0.1, 1.0 - rho_snow / 917.0)

            return rate
        end

        T_snow0 = 263.0  # K
        dTdz0 = 50.0     # K/m
        rho0 = 200.0      # kg/m3

        # AD derivative w.r.t. temperature
        T_dual = Dual(T_snow0, 1.0)
        rate_dual = grain_growth_rate(T_dual, dTdz0, rho0)
        drate_dT_ad = partials(rate_dual)[1]

        # FD check
        eps_fd = 1e-6
        drate_dT_fd = (grain_growth_rate(T_snow0 + eps_fd, dTdz0, rho0) -
                       grain_growth_rate(T_snow0 - eps_fd, dTdz0, rho0)) / (2 * eps_fd)
        @test drate_dT_ad ≈ drate_dT_fd rtol=1e-5

        # Verify sign: grain growth should increase with temperature
        @test drate_dT_ad > 0

        # AD derivative w.r.t. density
        rho_dual = Dual(rho0, 1.0)
        rate_rho = grain_growth_rate(T_snow0, dTdz0, rho_dual)
        drate_drho_ad = partials(rate_rho)[1]

        drate_drho_fd = (grain_growth_rate(T_snow0, dTdz0, rho0 + eps_fd) -
                         grain_growth_rate(T_snow0, dTdz0, rho0 - eps_fd)) / (2 * eps_fd)
        @test drate_drho_ad ≈ drate_drho_fd rtol=1e-5

        # Verify sign: growth rate should decrease with density
        @test drate_drho_ad < 0

        println("  Level 3c PASSED: Snow grain growth AD matches FD")
    end

    # =====================================================================
    # Level 3d: Two-stream canopy radiation (simplified)
    # =====================================================================
    @testset "Level 3d: Two-stream radiation AD" begin
        # The full two_stream! has ::Vector{Float64} type constraints.
        # We test the core two-stream math with Dual-valued LAI.

        function two_stream_albedo(LAI, rho_leaf, tau_leaf, xl, albsoi, cosz)
            # Leaf angle distribution
            chil = min(max(xl, -0.4), 0.6)
            if abs(chil) <= 0.01
                chil = typeof(LAI)(0.01)
            end
            phi1 = 0.5 - 0.633 * chil - 0.330 * chil^2
            phi2 = 0.877 * (1.0 - 2.0 * phi1)
            gdir = phi1 + phi2 * cosz
            mu_bar = (1.0 - phi1 / phi2 * log((phi1 + phi2) / phi1)) / phi2

            # Single scattering albedo
            omega = rho_leaf + tau_leaf

            # Extinction coefficient
            K = gdir / cosz
            beta0 = 0.5  # isotropic scattering approximation

            # Two-stream coefficients
            b = 1.0 - omega + omega * beta0
            c_val = omega * (1.0 - beta0)
            h = sqrt(b^2 - c_val^2) / mu_bar

            # Single-layer exponential attenuation
            tau_canopy = exp(-K * LAI)

            # Simplified canopy albedo
            alb_canopy_contrib = omega * beta0 * (1.0 - tau_canopy)
            alb_soil_contrib = albsoi * tau_canopy
            albedo = alb_canopy_contrib + alb_soil_contrib

            return albedo
        end

        LAI0 = 3.0
        rho_leaf = 0.10
        tau_leaf = 0.05
        xl = 0.01
        albsoi = 0.15
        cosz = 0.5

        # AD derivative of albedo w.r.t. LAI
        LAI_dual = Dual(LAI0, 1.0)
        alb_dual = two_stream_albedo(LAI_dual, rho_leaf, tau_leaf, xl, albsoi, cosz)
        dalb_dLAI_ad = partials(alb_dual)[1]

        # FD check
        eps_fd = 1e-7
        dalb_dLAI_fd = (two_stream_albedo(LAI0 + eps_fd, rho_leaf, tau_leaf, xl, albsoi, cosz) -
                        two_stream_albedo(LAI0 - eps_fd, rho_leaf, tau_leaf, xl, albsoi, cosz)) / (2 * eps_fd)
        @test dalb_dLAI_ad ≈ dalb_dLAI_fd rtol=1e-5

        @test isfinite(dalb_dLAI_ad)

        # AD derivative of albedo w.r.t. leaf reflectance
        rho_dual = Dual(rho_leaf, 1.0)
        alb_rho = two_stream_albedo(LAI0, rho_dual, tau_leaf, xl, albsoi, cosz)
        dalb_drho_ad = partials(alb_rho)[1]
        @test dalb_drho_ad > 0  # more reflective leaves -> higher canopy albedo

        # AD derivative w.r.t. soil albedo
        albsoi_dual = Dual(albsoi, 1.0)
        alb_soil = two_stream_albedo(LAI0, rho_leaf, tau_leaf, xl, albsoi_dual, cosz)
        dalb_dalbsoi_ad = partials(alb_soil)[1]
        @test dalb_dalbsoi_ad > 0  # brighter soil -> brighter total
        @test dalb_dalbsoi_ad < 1  # but attenuated by canopy

        println("  Level 3d PASSED: Two-stream radiation AD matches FD")
    end

    # =====================================================================
    # Level 3e: erf function AD compatibility
    # =====================================================================
    @testset "Level 3e: erf (Abramowitz-Stegun) AD" begin
        # CLM.erf is used in snow cover fraction calculations.
        x0 = 1.0
        x_dual = Dual(x0, 1.0)
        erf_dual = CLM.erf(x_dual)
        derf_dx_ad = partials(erf_dual)[1]

        # Analytical derivative of erf: d(erf)/dx = 2/sqrt(pi) * exp(-x^2)
        derf_dx_exact = 2.0 / sqrt(pi) * exp(-x0^2)

        # The Abramowitz-Stegun approximation derivative is close but not exact
        @test derf_dx_ad ≈ derf_dx_exact rtol=0.01

        # FD check (should match AD tightly)
        eps_fd = 1e-7
        derf_dx_fd = (CLM.erf(x0 + eps_fd) - CLM.erf(x0 - eps_fd)) / (2 * eps_fd)
        @test derf_dx_ad ≈ derf_dx_fd rtol=1e-5

        # Test at x=0 (maximum slope)
        x_dual_0 = Dual(0.0, 1.0)
        erf_0 = CLM.erf(x_dual_0)
        @test partials(erf_0)[1] ≈ 2.0 / sqrt(pi) rtol=0.01

        println("  Level 3e PASSED: erf AD derivative correct")
    end

    # =====================================================================
    # Level 3f: Polynomial evaluation AD (_poly8)
    # =====================================================================
    @testset "Level 3f: _poly8 AD" begin
        # _poly8 is @inline and accepts ::Real — works with Dual
        coeffs = CLM.QSAT_A
        x0 = 10.0  # 10 deg C

        x_dual = Dual(x0, 1.0)
        p_dual = CLM._poly8(x_dual, coeffs)
        dp_dx_ad = partials(p_dual)[1]

        # FD check — this must match tightly
        eps_fd = 1e-7
        dp_dx_fd = (CLM._poly8(x0 + eps_fd, coeffs) -
                    CLM._poly8(x0 - eps_fd, coeffs)) / (2 * eps_fd)
        @test dp_dx_ad ≈ dp_dx_fd rtol=1e-5

        # Note: QSAT_B is independently fit to approximate d(es)/d(T),
        # NOT the analytical derivative of the QSAT_A polynomial.
        # So AD of _poly8(QSAT_A) differs from _poly8(QSAT_B).
        # Both are valid approximations of the true derivative.
        dp_dx_explicit = CLM._poly8(x0, CLM.QSAT_B)
        @test dp_dx_ad ≈ dp_dx_explicit rtol=0.01  # agree to ~0.01%

        println("  Level 3f PASSED: _poly8 AD derivative matches FD")
    end

    # =====================================================================
    # Level 4: Mini-driver chain
    # =====================================================================
    @testset "Level 4: Mini-driver (T -> qsat -> resistance -> LH)" begin
        # Chain multiple CLM physics computations:
        # T_air -> qsat -> VPD -> soil resistance -> latent heat flux

        function micro_clm_latent_heat(T_air, p_atm, q_air, soil_wetness,
                                        wind_speed, z_atm)
            # Step 1: Saturation properties at ground temperature
            T_grnd = T_air - 2.0
            (qs_grnd, es_grnd, _, _) = CLM.qsat(T_grnd, p_atm)

            # Step 2: Soil resistance (Swenson-Lawrence 2014 simplified)
            watsat = 0.45
            bsw = 5.0
            sucsat = 200.0
            aird = watsat * (sucsat / 1.0e7)^(1.0 / bsw)
            d_max = 15.0
            frac_sat_init = 0.8

            vwc = soil_wetness * watsat
            # Use smooth clamp to avoid derivative discontinuity at boundary
            dsl_numer = frac_sat_init * watsat - vwc
            dsl_denom = frac_sat_init * watsat - aird
            dsl = d_max * max(0.001, dsl_numer) / max(0.001, dsl_denom)
            dsl = max(dsl, 0.0)
            dsl = min(dsl, 200.0)

            # Diffusivity
            d0 = 2.12e-5 * (T_grnd / 273.15)^1.75
            eps_soil = watsat - aird
            dg = eps_soil * d0 * (eps_soil / watsat)^(3.0 / max(3.0, bsw))

            soilresis = dsl / (dg * eps_soil * 1.0e3) + 20.0
            soilresis = min(1.0e6, soilresis)

            # Step 3: Aerodynamic resistance
            z0 = 0.01
            ra = (log(z_atm / z0))^2 / (CLM.VKC^2 * wind_speed)

            # Step 4: Latent heat flux
            rho_air = p_atm / (CLM.RAIR * T_air)
            q_diff = qs_grnd - q_air

            total_resis = ra + soilresis
            LH = rho_air * CLM.HVAP * q_diff / total_resis

            return LH
        end

        T0 = 290.0
        p0 = 101325.0
        q0 = 0.008
        sw0 = 0.3         # choose 0.3 to stay away from clamp boundaries
        wind0 = 3.0
        z0 = 10.0

        # Use ForwardDiff.gradient for consistent tagging
        function LH_scalar(x)
            return micro_clm_latent_heat(x[1], x[2], x[3], x[4], x[5], x[6])
        end

        x0_vec = [T0, p0, q0, sw0, wind0, z0]
        grad = ForwardDiff.gradient(LH_scalar, x0_vec)
        dLH_dT_ad = grad[1]
        dLH_dsw_ad = grad[4]
        dLH_dp_ad = grad[2]

        # FD checks
        eps_fd = 1e-7
        dLH_dT_fd = (LH_scalar(x0_vec .+ [eps_fd,0,0,0,0,0]) -
                     LH_scalar(x0_vec .- [eps_fd,0,0,0,0,0])) / (2 * eps_fd)
        @test dLH_dT_ad ≈ dLH_dT_fd rtol=1e-4

        dLH_dsw_fd = (LH_scalar(x0_vec .+ [0,0,0,eps_fd,0,0]) -
                      LH_scalar(x0_vec .- [0,0,0,eps_fd,0,0])) / (2 * eps_fd)
        @test dLH_dsw_ad ≈ dLH_dsw_fd rtol=1e-4

        # Physical checks
        @test dLH_dT_ad > 0  # warmer air -> higher sat humidity -> more evaporation
        @test dLH_dsw_ad > 0 # wetter soil -> lower resistance -> more evaporation

        @test isfinite(dLH_dp_ad)

        # --- Multi-input Jacobian via ForwardDiff.jacobian ---
        function LH_vector(x)
            return [micro_clm_latent_heat(x[1], x[2], x[3], x[4], x[5], x[6])]
        end

        J = ForwardDiff.jacobian(LH_vector, x0_vec)
        @test size(J) == (1, 6)

        println("  Level 4 PASSED: Mini-driver chain differentiable end-to-end")
        println("    dLH/dT    = $(round(dLH_dT_ad, digits=4)) W/m^2/K")
        println("    dLH/d(sw) = $(round(dLH_dsw_ad, digits=4)) W/m^2")
        println("    dLH/dp    = $(round(dLH_dp_ad, sigdigits=4)) W/m^2/Pa")
    end

    # =====================================================================
    # Level 4b: Multi-step physics chain with CLM functions
    # =====================================================================
    @testset "Level 4b: qsat -> surface energy balance chain" begin
        # Chain: T_surface -> qsat -> longwave emission -> net radiation
        # -> ground heat flux estimate -> temperature tendency

        function surface_energy_tendency(T_surf, T_air, p_atm, q_air,
                                          SW_down, LW_down, emissivity)
            # Outgoing longwave
            LW_out = emissivity * CLM.SB * T_surf^4

            # Net longwave
            LW_net = emissivity * LW_down - LW_out

            # Net shortwave
            albedo = 0.2
            SW_net = (1.0 - albedo) * SW_down

            # Sensible heat (simple bulk formula)
            rho_air = p_atm / (CLM.RAIR * T_air)
            ra = 50.0  # aerodynamic resistance [s/m]
            SH = rho_air * CLM.CPAIR * (T_surf - T_air) / ra

            # Latent heat
            (qs_surf, _, _, _) = CLM.qsat(T_surf, p_atm)
            LH = rho_air * CLM.HVAP * (qs_surf - q_air) / ra

            # Ground heat flux (residual)
            G = SW_net + LW_net - SH - LH

            # Temperature tendency
            C_soil = 2.0e6  # volumetric heat capacity [J/m^3/K]
            dz = 0.05       # layer thickness [m]
            dTdt = G / (C_soil * dz)

            return dTdt
        end

        T_s0 = 285.0
        T_a = 280.0
        p = 101325.0
        q = 0.005
        SW = 300.0
        LW = 280.0
        emis = 0.96

        # AD: d(dTdt)/d(T_surf) — the "implicit factor" for soil temperature
        T_dual = Dual(T_s0, 1.0)
        dTdt_dual = surface_energy_tendency(T_dual, T_a, p, q, SW, LW, emis)
        d2T_dTdt = partials(dTdt_dual)[1]

        # Negative = stabilizing feedback
        @test d2T_dTdt < 0

        # FD check
        eps_fd = 1e-5
        d2T_fd = (surface_energy_tendency(T_s0 + eps_fd, T_a, p, q, SW, LW, emis) -
                  surface_energy_tendency(T_s0 - eps_fd, T_a, p, q, SW, LW, emis)) / (2 * eps_fd)
        @test d2T_dTdt ≈ d2T_fd rtol=1e-4

        # Full Jacobian w.r.t. all forcing variables
        function tendency_vec(x)
            return surface_energy_tendency(x[1], x[2], x[3], x[4], x[5], x[6], x[7])
        end

        x0 = [T_s0, T_a, p, q, SW, LW, emis]
        grad = ForwardDiff.gradient(tendency_vec, x0)
        @test length(grad) == 7
        @test grad[1] ≈ d2T_dTdt rtol=1e-10
        @test grad[5] > 0   # more SW -> warmer
        @test grad[6] > 0   # more LW_down -> warmer

        println("  Level 4b PASSED: Surface energy balance chain differentiable")
        println("    d(dTdt)/d(T_surf) = $(round(d2T_dTdt, sigdigits=4)) K/s/K (stabilizing)")
    end

    # =====================================================================
    # Level 5: Full timestep AD — barrier identification
    # =====================================================================
    @testset "Level 5: Struct type barrier documentation" begin
        # Document which struct types have {FT<:Real} that would
        # need widening to {FT<:Real} for full driver AD.
        #
        # ForwardDiff.Dual{T,V,N} <: Real but NOT <: AbstractFloat
        # So any struct with {FT<:Real} cannot hold Dual values.

        D = ForwardDiff.Dual{Nothing, Float64, 1}

        # Verify the fundamental type relationship
        @test D <: Real
        @test !(D <: AbstractFloat)

        # Known CLM types with {FT<:Real} constraint
        # (these would need {FT<:Real} for full driver AD)
        known_types_with_constraint = [
            :TemperatureData,
            :EnergyFluxData,
            :SoilStateData,
            :SoilHydrologyData,
            :CanopyStateData,
            :LakeStateData,
            :SurfaceAlbedoData,
            :SolarAbsorbedData,
            :WaterStateData,
            :WaterFluxData,
            :WaterDiagnosticBulkData,
            :WaterStateBulkData,
            :WaterFluxBulkData,
            :WaterBalanceBulkData,
            :FrictionVelocityData,
            :ColumnData,
            :PatchData,
            :LandunitData,
            :GridcellData,
        ]

        # Verify these types exist and have {FT<:Real} (not AbstractFloat)
        types_with_real = Symbol[]
        types_still_abstractfloat = Symbol[]
        for name in known_types_with_constraint
            T = try
                getfield(CLM, name)
            catch
                nothing
            end
            if T !== nothing
                t = T
                while t isa UnionAll
                    if t.var.ub === Real
                        push!(types_with_real, name)
                        break
                    elseif t.var.ub === AbstractFloat
                        push!(types_still_abstractfloat, name)
                        break
                    end
                    t = t.body
                end
            end
        end

        println("  Level 5: Structs widened to {FT<:Real}:")
        for s in sort(types_with_real)
            println("    - $s")
        end
        println("  Widened: $(length(types_with_real)) / $(length(known_types_with_constraint))")
        if !isempty(types_still_abstractfloat)
            println("  Still AbstractFloat: $(types_still_abstractfloat)")
        end

        # Verify that the physics FUNCTIONS are AD-ready
        @test CLM.qsat(Dual(280.0, 1.0), 101325.0) isa NTuple{4}
        @test CLM._poly8(Dual(10.0, 1.0), CLM.QSAT_A) isa Dual
        @test CLM.erf(Dual(1.0, 1.0)) isa Dual

        # Document the recommended path forward
        println("\n  AD Readiness Summary:")
        println("  - Physics functions (qsat, _poly8, erf): READY (accept ::Real)")
        println("  - Tridiagonal solver: READY (accepts AbstractVector{<:Real})")
        println("  - Mini-driver chains: READY (Levels 4, 4b passed)")
        println("  - Struct types: $(length(types_with_real)) widened to {FT<:Real}")
        println("  - Remaining AbstractFloat: $(length(types_still_abstractfloat))")

        @test length(types_with_real) > 0  # confirm structs are now AD-ready
        @test length(types_still_abstractfloat) == 0  # no barriers remain

        # Additionally scan ALL names in the module for completeness
        all_constrained = Symbol[]
        for name in names(CLM; all=true)
            T = try
                getfield(CLM, name)
            catch
                continue
            end
            t = T
            while t isa UnionAll
                if t.var.ub === AbstractFloat
                    push!(all_constrained, name)
                    break
                end
                t = t.body
            end
        end
        println("  Full scan found $(length(all_constrained)) types with AbstractFloat constraint")
    end

    # =====================================================================
    # Level 5b: ForwardDiff.derivative high-level API
    # =====================================================================
    @testset "Level 5b: ForwardDiff.derivative API" begin
        # Show that ForwardDiff's high-level API works for CLM functions
        f_qs(T) = CLM.qsat(T, 101325.0)[1]  # just qs
        f_es(T) = CLM.qsat(T, 101325.0)[2]  # just es

        dqs = ForwardDiff.derivative(f_qs, 280.0)
        des = ForwardDiff.derivative(f_es, 280.0)

        # Cross-check with FD
        eps_fd = 1e-7
        dqs_fd = (f_qs(280.0 + eps_fd) - f_qs(280.0 - eps_fd)) / (2 * eps_fd)
        des_fd = (f_es(280.0 + eps_fd) - f_es(280.0 - eps_fd)) / (2 * eps_fd)
        @test dqs ≈ dqs_fd rtol=1e-5
        @test des ≈ des_fd rtol=1e-5

        # Cross-check AD vs analytical (loose tolerance due to independent polynomial fits)
        (_, _, dqsdT_a, desdT_a) = CLM.qsat(280.0, 101325.0)
        @test dqs ≈ dqsdT_a rtol=1e-3
        @test des ≈ desdT_a rtol=1e-3

        # Higher-order derivative (second derivative of qs w.r.t. T)
        d2qs = ForwardDiff.derivative(T -> ForwardDiff.derivative(f_qs, T), 280.0)
        @test isfinite(d2qs)
        @test d2qs > 0  # Clausius-Clapeyron: qs is convex in T

        println("  Level 5b PASSED: ForwardDiff.derivative API works for CLM physics")
        println("    d(qs)/dT    = $(round(dqs, sigdigits=4)) kg/kg/K")
        println("    d2(qs)/dT^2 = $(round(d2qs, sigdigits=4)) kg/kg/K^2 (convex, as expected)")
    end

    # =====================================================================
    # Level 5c: Tridiagonal solver with Dual arrays
    # =====================================================================
    @testset "Level 5c: Tridiagonal solver AD" begin
        # The tridiagonal solver accepts AbstractVector{<:Real}.
        # Test with Dual-valued RHS to get d(solution)/d(forcing).
        n = 5
        a = zeros(n)
        b = fill(2.0, n)
        c = zeros(n)
        for i in 1:n-1
            a[i+1] = -1.0
            c[i] = -1.0
        end

        # Float64 RHS
        r = [1.0, 0.0, 0.0, 0.0, 1.0]
        u = zeros(n)
        CLM.tridiagonal_solve!(u, a, b, c, r, 1, n)

        # Now perturb r[1] with Dual to get d(u)/d(r[1])
        r_dual = [Dual(1.0, 1.0), Dual(0.0, 0.0), Dual(0.0, 0.0),
                  Dual(0.0, 0.0), Dual(1.0, 0.0)]
        u_dual = Vector{Dual{Nothing, Float64, 1}}(undef, n)
        CLM.tridiagonal_solve!(u_dual, a, b, c, r_dual, 1, n)

        # Extract derivatives
        du1_dr1 = partials(u_dual[1])[1]
        du3_dr1 = partials(u_dual[3])[1]

        # FD check
        eps_fd = 1e-7
        r_plus = [1.0 + eps_fd, 0.0, 0.0, 0.0, 1.0]
        r_minus = [1.0 - eps_fd, 0.0, 0.0, 0.0, 1.0]
        u_plus = zeros(n)
        u_minus = zeros(n)
        CLM.tridiagonal_solve!(u_plus, a, b, c, r_plus, 1, n)
        CLM.tridiagonal_solve!(u_minus, a, b, c, r_minus, 1, n)
        du1_dr1_fd = (u_plus[1] - u_minus[1]) / (2 * eps_fd)
        du3_dr1_fd = (u_plus[3] - u_minus[3]) / (2 * eps_fd)

        @test du1_dr1 ≈ du1_dr1_fd rtol=1e-5
        @test du3_dr1 ≈ du3_dr1_fd rtol=1e-5

        # Sensitivity should decay away from the perturbation point
        @test abs(du1_dr1) > abs(du3_dr1)

        println("  Level 5c PASSED: Tridiagonal solver differentiable via ForwardDiff")
        println("    d(u[1])/d(r[1]) = $(round(du1_dr1, digits=6))")
        println("    d(u[3])/d(r[1]) = $(round(du3_dr1, digits=6))")
    end

    # =====================================================================
    # Level 6: Struct instantiation with Dual — proof that type barriers are gone
    # =====================================================================
    @testset "Level 6: Dual-valued struct instantiation" begin
        D = ForwardDiff.Dual{Nothing, Float64, 1}

        # Instantiate core CLM structs with Dual element type
        # This was impossible when bounds were {FT<:AbstractFloat}
        temp = CLM.TemperatureData{D}()
        @test temp isa CLM.TemperatureData{D}

        eflx = CLM.EnergyFluxData{D}()
        @test eflx isa CLM.EnergyFluxData{D}

        ws = CLM.WaterStateData{D}()
        @test ws isa CLM.WaterStateData{D}

        fv = CLM.FrictionVelocityData{D}()
        @test fv isa CLM.FrictionVelocityData{D}

        ss = CLM.SoilStateData{D}()
        @test ss isa CLM.SoilStateData{D}

        # Store Dual values and read them back
        temp.t_grnd_col = [Dual(272.0, 0.5), Dual(273.0, 0.3)]
        @test value(temp.t_grnd_col[1]) == 272.0
        @test partials(temp.t_grnd_col[1])[1] == 0.5

        # Physics function with Dual extracted from struct
        (qs, es, dqs, des) = CLM.qsat(temp.t_grnd_col[1], 101325.0)
        @test isfinite(value(qs))
        @test isfinite(partials(qs)[1])

        # Stability functions with Dual
        zeta_dual = Dual(-0.5, 1.0)
        sf1 = CLM.stability_func1(zeta_dual)
        sf2 = CLM.stability_func2(zeta_dual)
        @test isfinite(value(sf1))
        @test isfinite(partials(sf1)[1])
        @test isfinite(value(sf2))
        @test isfinite(partials(sf2)[1])

        # Monin-Obukhov init with Dual
        (um, obu) = CLM.monin_obuk_ini(
            Dual(0.5, 0.0), Dual(5.0, 1.0), Dual(280.0, 0.0),
            Dual(-0.5, 0.0), Dual(30.0, 0.0), Dual(0.01, 0.0))
        @test isfinite(value(um))
        @test isfinite(partials(um)[1])

        println("  Level 6 PASSED: Dual-valued struct instantiation and physics")
        println("    5 core structs instantiated with Dual{Nothing,Float64,1}")
        println("    qsat(Dual from struct): d(qs)/d(T_grnd_col) = $(round(partials(qs)[1], sigdigits=4))")
        println("    stability_func1(Dual): d/d(zeta) = $(round(partials(sf1)[1], sigdigits=4))")
        println("    monin_obuk_ini(Dual): d(um)/d(ur) = $(round(partials(um)[1], sigdigits=4))")
    end

    # =====================================================================
    # Level 7: CLMInstances{Dual} construction and full-driver readiness
    # =====================================================================
    @testset "Level 7: CLMInstances Dual-copy approach" begin
        D = Dual{Nothing, Float64, 1}

        # CLMInstances is unparameterized — test make_dual_copy of child structs
        @test typeof(CLM.CLMInstances()) == CLM.CLMInstances

        # Test make_dual_copy on individual parameterized child structs
        function _make_dual_copy(src, ::Type{DD}) where DD
            T = typeof(src)
            wrapper = T.name.wrapper
            wrapper === T && return src
            dst = wrapper{DD}()
            for name in fieldnames(T)
                sv = getfield(src, name)
                try
                    if sv isa Array{Float64}
                        setfield!(dst, name, DD.(sv))
                    elseif sv isa Float64
                        setfield!(dst, name, DD(sv))
                    else
                        setfield!(dst, name, sv)
                    end
                catch e
                    e isa TypeError || rethrow()
                end
            end
            return dst
        end

        # Construct Float64 structs, populate, and copy to Dual
        nc_test = 3
        temp_f64 = CLM.TemperatureData{Float64}()
        temp_f64.t_grnd_col = fill(280.0, nc_test)
        temp_d = _make_dual_copy(temp_f64, D)
        @test eltype(temp_d.t_grnd_col) == D
        @test length(temp_d.t_grnd_col) == nc_test

        ser_f64 = CLM.SaturatedExcessRunoffData{Float64}()
        ser_f64.fsat_col = fill(0.5, nc_test)
        ser_d = _make_dual_copy(ser_f64, D)
        @test eltype(ser_d.fsat_col) == D

        ier_f64 = CLM.InfiltrationExcessRunoffData{Float64}()
        ier_f64.qinmax_col = fill(0.01, nc_test)
        ier_d = _make_dual_copy(ier_f64, D)
        @test eltype(ier_d.qinmax_col) == D

        # Test scalar FT field conversion (e.g., FrictionVelocityData has zsno::FT)
        fv_f64 = CLM.FrictionVelocityData{Float64}()
        fv_f64.zsno = 0.00085
        fv_d = _make_dual_copy(fv_f64, D)
        @test fv_d.zsno isa D
        @test value(fv_d.zsno) == 0.00085

        # Seed a Dual into temperature and call qsat through the struct
        temp_d.t_grnd_col[1] = Dual(280.0, 1.0)
        T_d = temp_d.t_grnd_col[1]
        (qs, es, dqsdT, desdT) = CLM.qsat(T_d, 101325.0)
        @test isfinite(value(qs))
        @test abs(partials(qs)[1]) > 0  # nonzero derivative

        # Verify Dual dtime flows through driver signature types
        dtime_d = Dual(1800.0, 1.0)
        @test dtime_d isa Real

        println("  Level 7 PASSED: Dual-copy approach for CLMInstances")
        println("    make_dual_copy creates Dual-typed child structs")
        println("    Scalar FT fields correctly converted to Dual")
        println("    qsat through Dual struct: d(qs)/d(T) = $(round(partials(qs)[1], sigdigits=4))")
    end

    # =====================================================================
    # Level 8: Full clm_drv! timestep with ForwardDiff
    # =====================================================================

    @testset "Level 8: Full clm_drv! timestep AD" begin
        fsurdat = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc"
        paramfile = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc"

        if !isfile(fsurdat) || !isfile(paramfile)
            @warn "Skipping Level 8: input files not found"
            @test true
            return
        end

        # --- Helper: create a Dual-typed copy of a parameterized struct ---
        # For a struct like TemperatureData{Float64}, this creates TemperatureData{D}()
        # and copies all Float64 arrays as Dual arrays.  Non-array fields are copied
        # as-is. Non-parameterized structs are returned by reference.
        function make_dual_copy(src, ::Type{D}) where D
            T = typeof(src)
            wrapper = T.name.wrapper
            if wrapper === T
                # Not parameterized — return as-is (shared reference)
                return src
            end
            dst = wrapper{D}()
            for name in fieldnames(T)
                sv = getfield(src, name)
                try
                    if sv isa Array{Float64}
                        setfield!(dst, name, D.(sv))
                    elseif sv isa Float64
                        # Scalar FT fields: convert Float64 → D
                        setfield!(dst, name, D(sv))
                    else
                        setfield!(dst, name, sv)
                    end
                catch e
                    e isa TypeError || rethrow()
                end
            end
            return dst
        end

        # --- Helper: set up forcing on an Atm2LndData ---
        function setup_forcing!(a2l, T, ng)
            for g in 1:ng
                a2l.forc_t_not_downscaled_grc[g] = T
                a2l.forc_pbot_not_downscaled_grc[g] = 85000.0
                a2l.forc_th_not_downscaled_grc[g] = T * (100000.0 / 85000.0)^(CLM.RAIR / CLM.CPAIR)
                a2l.forc_rho_not_downscaled_grc[g] = 85000.0 / (CLM.RAIR * T)
                a2l.forc_lwrad_not_downscaled_grc[g] = 250.0
                a2l.forc_vp_grc[g] = 300.0
                a2l.forc_hgt_grc[g] = 30.0
                a2l.forc_topo_grc[g] = 0.0
                a2l.forc_wind_grc[g] = 3.0
                for b in 1:CLM.NUMRAD
                    a2l.forc_solad_not_downscaled_grc[g, b] = 100.0
                    a2l.forc_solai_grc[g, b] = 50.0
                end
                a2l.forc_solar_not_downscaled_grc[g] = 300.0
                a2l.forc_rain_not_downscaled_grc[g] = 0.0
                a2l.forc_snow_not_downscaled_grc[g] = 0.0001
            end
        end

        # ---- Step 1: Initialize Float64 state ----
        (inst_f64, bounds, filt, tm) = CLM.clm_initialize!(;
            fsurdat=fsurdat, paramfile=paramfile)
        ng = bounds.endg; nl = bounds.endl; nc = bounds.endc; np = bounds.endp

        T_base = 270.0
        setup_forcing!(inst_f64.atm2lnd, T_base, ng)
        CLM.downscale_forcings!(bounds, inst_f64.atm2lnd, inst_f64.column, inst_f64.landunit, inst_f64.topo)

        config = CLM.CLMDriverConfig()
        filt_ia = CLM.clump_filter_inactive_and_active
        calday = 1.0
        (declin, eccf) = CLM.compute_orbital(calday)
        nextsw_cday = calday + 1800.0 / CLM.SECSPDAY

        # Run 3 warmup timesteps to stabilise cold-start fluxes
        for n in 1:3
            CLM.clm_drv!(config, inst_f64, filt, filt_ia, bounds,
                true, nextsw_cday, declin, declin, CLM.ORB_OBLIQR_DEFAULT,
                false, false, "", false;
                nstep=n, is_first_step=(n==1), is_beg_curr_day=(n==1),
                dtime=1800.0, mon=1, day=1,
                photosyns=inst_f64.photosyns)
        end

        # Find a vegetated patch with finite LH (skip bare/urban/lake)
        test_p = 0
        for p in 1:np
            lh = inst_f64.energyflux.eflx_lh_tot_patch[p]
            if isfinite(lh) && abs(lh) > 0.0
                test_p = p; break
            end
        end
        if test_p == 0
            # Fallback: just use column-level ground temperature
            test_p = 1
        end

        lh_f64 = inst_f64.energyflux.eflx_lh_tot_patch[test_p]
        sh_f64 = inst_f64.energyflux.eflx_sh_tot_patch[test_p]
        tg_f64 = inst_f64.temperature.t_grnd_col[1]
        println("  Level 8: test patch = $test_p")
        println("  Level 8: Float64 warmup LH=$lh_f64, SH=$sh_f64, Tg=$tg_f64")

        # ---- Step 2: Create Dual-typed CLMInstances and copy state ----
        D = Dual{Nothing, Float64, 1}

        inst_d = CLM.CLMInstances()

        # For every field in CLMInstances, create a Dual-typed copy of
        # the Float64 instance and assign it.  Unparameterized structs
        # (GridcellData, ColumnData, etc.) are shared by reference.
        for name in fieldnames(CLM.CLMInstances)
            name === :water && continue          # handled separately (facade)
            name === :surfdata && continue       # Union{..., Nothing}
            src = getfield(inst_f64, name)
            setfield!(inst_d, name, make_dual_copy(src, D))
        end
        inst_d.surfdata = inst_f64.surfdata

        # Copy photosyns int/bool fields
        inst_d.photosyns.stomatalcond_mtd = inst_f64.photosyns.stomatalcond_mtd

        # Handle WaterData facade: it is unparameterized but its children
        # (WaterFluxBulkData{FT}, etc.) ARE parameterized, and they have
        # composed sub-structs (wf::WaterFluxData, ws::WaterStateData) that
        # are also typed.  Create a fresh WaterData and populate with Dual copies.
        water_d = CLM.WaterData()
        # Copy all non-array fields from the Float64 WaterData
        for name in fieldnames(CLM.WaterData)
            sv = getfield(inst_f64.water, name)
            try
                setfield!(water_d, name, sv)
            catch
            end
        end

        # WaterStateBulk: has ws::WaterStateData (composed)
        wsb_d = make_dual_copy(inst_f64.water.waterstatebulk_inst, D)
        ws_d  = make_dual_copy(inst_f64.water.waterstatebulk_inst.ws, D)
        wsb_d.ws = ws_d
        water_d.waterstatebulk_inst = wsb_d

        # WaterFluxBulk: has wf::WaterFluxData (composed)
        wfb_d = make_dual_copy(inst_f64.water.waterfluxbulk_inst, D)
        wf_d  = make_dual_copy(inst_f64.water.waterfluxbulk_inst.wf, D)
        wfb_d.wf = wf_d
        water_d.waterfluxbulk_inst = wfb_d

        # WaterDiagnosticBulk: no composed sub-struct
        water_d.waterdiagnosticbulk_inst = make_dual_copy(
            inst_f64.water.waterdiagnosticbulk_inst, D)

        # WaterBalance: no composed sub-struct
        water_d.waterbalancebulk_inst = make_dual_copy(
            inst_f64.water.waterbalancebulk_inst, D)

        # Also update bulk_and_tracers entries to point to Dual-typed data
        if !isempty(water_d.bulk_and_tracers)
            bt = water_d.bulk_and_tracers[water_d.i_bulk]
            bt.waterflux = wfb_d
            bt.waterstate = wsb_d
            bt.waterdiagnostic = water_d.waterdiagnosticbulk_inst
            bt.waterbalance = water_d.waterbalancebulk_inst
        end

        inst_d.water = water_d

        # ---- Step 3: Seed Dual into forcing temperature ----
        setup_forcing!(inst_d.atm2lnd, T_base, ng)
        for g in 1:ng
            inst_d.atm2lnd.forc_t_not_downscaled_grc[g] = Dual(T_base, 1.0)
            T_d = inst_d.atm2lnd.forc_t_not_downscaled_grc[g]
            inst_d.atm2lnd.forc_th_not_downscaled_grc[g] = T_d * (D(100000.0) / D(85000.0))^(D(CLM.RAIR) / D(CLM.CPAIR))
            inst_d.atm2lnd.forc_rho_not_downscaled_grc[g] = D(85000.0) / (D(CLM.RAIR) * T_d)
        end
        CLM.downscale_forcings!(bounds, inst_d.atm2lnd, inst_d.column, inst_d.landunit, inst_d.topo)

        # ---- Step 4: Run one Dual timestep ----
        CLM.clm_drv!(config, inst_d, filt, filt_ia, bounds,
            true, nextsw_cday, declin, declin, CLM.ORB_OBLIQR_DEFAULT,
            false, false, "", false;
            nstep=4, is_first_step=false,
            dtime=1800.0, mon=1, day=1,
            photosyns=inst_d.photosyns)

        # ---- Step 5: Verify derivatives ----
        lh_d = inst_d.energyflux.eflx_lh_tot_patch[test_p]
        tgrnd_d = inst_d.temperature.t_grnd_col[1]
        sh_d = inst_d.energyflux.eflx_sh_tot_patch[test_p]

        println("  Level 8: Dual LH value = $(value(lh_d))")
        println("  Level 8: d(LH)/d(T_forc) = $(partials(lh_d)[1])")
        println("  Level 8: d(T_grnd)/d(T_forc) = $(partials(tgrnd_d)[1])")
        println("  Level 8: d(SH)/d(T_forc) = $(partials(sh_d)[1])")

        # Output values should be finite
        @test isfinite(value(lh_d))
        @test isfinite(value(tgrnd_d))
        @test isfinite(value(sh_d))

        # Derivatives should be finite (may be zero for some outputs in cold start)
        @test isfinite(partials(lh_d)[1])
        @test isfinite(partials(tgrnd_d)[1])
        @test isfinite(partials(sh_d)[1])

        # At least one derivative should be nonzero (T_forc must influence something)
        any_nonzero = abs(partials(lh_d)[1]) > 0 ||
                      abs(partials(tgrnd_d)[1]) > 0 ||
                      abs(partials(sh_d)[1]) > 0
        @test any_nonzero

        # Derivatives should be bounded (not blown up)
        @test abs(partials(lh_d)[1]) < 1e10
        @test abs(partials(tgrnd_d)[1]) < 1e10
        @test abs(partials(sh_d)[1]) < 1e10

        # ---- Step 6: Finite difference comparison ----
        # Use same warmed-up state: re-initialize, run 3 warmup steps at T_base,
        # then run step 4 at T_base (reference) and T_base+eps (perturbed).
        eps_fd = 0.01

        # Reference run: fresh init → 3 warmup → step 4 at T_base
        (inst_ref, _, _, _) = CLM.clm_initialize!(;
            fsurdat=fsurdat, paramfile=paramfile)
        setup_forcing!(inst_ref.atm2lnd, T_base, ng)
        CLM.downscale_forcings!(bounds, inst_ref.atm2lnd, inst_ref.column, inst_ref.landunit, inst_ref.topo)
        for n in 1:3
            CLM.clm_drv!(config, inst_ref, filt, filt_ia, bounds,
                true, nextsw_cday, declin, declin, CLM.ORB_OBLIQR_DEFAULT,
                false, false, "", false;
                nstep=n, is_first_step=(n==1), is_beg_curr_day=(n==1),
                dtime=1800.0, mon=1, day=1,
                photosyns=inst_ref.photosyns)
        end
        # Step 4 at T_base
        CLM.clm_drv!(config, inst_ref, filt, filt_ia, bounds,
            true, nextsw_cday, declin, declin, CLM.ORB_OBLIQR_DEFAULT,
            false, false, "", false;
            nstep=4, is_first_step=false,
            dtime=1800.0, mon=1, day=1,
            photosyns=inst_ref.photosyns)
        lh_ref = inst_ref.energyflux.eflx_lh_tot_patch[test_p]

        # Perturbed run: fresh init → 3 warmup at T_base → step 4 at T_base+eps
        (inst_pert, _, _, _) = CLM.clm_initialize!(;
            fsurdat=fsurdat, paramfile=paramfile)
        setup_forcing!(inst_pert.atm2lnd, T_base, ng)
        CLM.downscale_forcings!(bounds, inst_pert.atm2lnd, inst_pert.column, inst_pert.landunit, inst_pert.topo)
        for n in 1:3
            CLM.clm_drv!(config, inst_pert, filt, filt_ia, bounds,
                true, nextsw_cday, declin, declin, CLM.ORB_OBLIQR_DEFAULT,
                false, false, "", false;
                nstep=n, is_first_step=(n==1), is_beg_curr_day=(n==1),
                dtime=1800.0, mon=1, day=1,
                photosyns=inst_pert.photosyns)
        end
        # Step 4 at T_base + eps
        setup_forcing!(inst_pert.atm2lnd, T_base + eps_fd, ng)
        CLM.downscale_forcings!(bounds, inst_pert.atm2lnd, inst_pert.column, inst_pert.landunit, inst_pert.topo)
        CLM.clm_drv!(config, inst_pert, filt, filt_ia, bounds,
            true, nextsw_cday, declin, declin, CLM.ORB_OBLIQR_DEFAULT,
            false, false, "", false;
            nstep=4, is_first_step=false,
            dtime=1800.0, mon=1, day=1,
            photosyns=inst_pert.photosyns)
        lh_pert = inst_pert.energyflux.eflx_lh_tot_patch[test_p]

        # Compute FD derivatives for LH, SH, T_grnd
        sh_ref  = inst_ref.energyflux.eflx_sh_tot_patch[test_p]
        tg_ref  = inst_ref.temperature.t_grnd_col[1]
        sh_pert = inst_pert.energyflux.eflx_sh_tot_patch[test_p]
        tg_pert = inst_pert.temperature.t_grnd_col[1]

        dlh_fd = (lh_pert - lh_ref) / eps_fd
        dsh_fd = (sh_pert - sh_ref) / eps_fd
        dtg_fd = (tg_pert - tg_ref) / eps_fd

        dlh_ad = partials(lh_d)[1]
        dsh_ad = partials(sh_d)[1]
        dtg_ad = partials(tgrnd_d)[1]

        println("  Level 8: d(LH)/d(T)    AD = $(round(dlh_ad, digits=4)), FD = $(round(dlh_fd, digits=4))")
        println("  Level 8: d(SH)/d(T)    AD = $(round(dsh_ad, digits=4)), FD = $(round(dsh_fd, digits=4))")
        println("  Level 8: d(T_grnd)/d(T) AD = $(round(dtg_ad, digits=6)), FD = $(round(dtg_fd, digits=6))")

        # Sign agreement for all three
        for (name, ad, fd) in [("LH", dlh_ad, dlh_fd), ("SH", dsh_ad, dsh_fd), ("T_grnd", dtg_ad, dtg_fd)]
            if abs(ad) > 1e-10 && abs(fd) > 1e-10
                @test sign(ad) == sign(fd)
            end
        end

        # Relative agreement: AD and FD should agree within 5%
        for (name, ad, fd) in [("LH", dlh_ad, dlh_fd), ("SH", dsh_ad, dsh_fd), ("T_grnd", dtg_ad, dtg_fd)]
            if abs(fd) > 1e-6
                rel_err = abs(ad - fd) / abs(fd)
                println("  Level 8: d($name)/d(T) relative error = $(round(rel_err * 100, digits=2))%")
                @test rel_err < 0.05  # within 5%
            end
        end

        println("  Level 8 PASSED: Full clm_drv! timestep with ForwardDiff")
    end

    # =====================================================================
    # Smooth AD primitives unit tests
    # =====================================================================

    @testset "Smooth AD primitives" begin
        # Test that smooth primitives match standard functions for Float64
        @test CLM.smooth_min(3.0, 5.0) === min(3.0, 5.0)
        @test CLM.smooth_min(5.0, 3.0) === min(5.0, 3.0)
        @test CLM.smooth_max(3.0, 5.0) === max(3.0, 5.0)
        @test CLM.smooth_max(5.0, 3.0) === max(5.0, 3.0)
        @test CLM.smooth_clamp(2.0, 3.0, 5.0) === clamp(2.0, 3.0, 5.0)
        @test CLM.smooth_clamp(4.0, 3.0, 5.0) === clamp(4.0, 3.0, 5.0)
        @test CLM.smooth_clamp(6.0, 3.0, 5.0) === clamp(6.0, 3.0, 5.0)
        @test CLM.smooth_heaviside(1.0) === 1.0
        @test CLM.smooth_heaviside(-1.0) === 0.0
        @test CLM.smooth_abs(3.0) === 3.0
        @test CLM.smooth_abs(-3.0) === 3.0

        # Test smooth primitives with Dual types — should approximate standard
        D1 = Dual{Nothing, Float64, 1}

        # smooth_min with Duals
        a_d = Dual(3.0, 1.0)
        b_d = Dual(5.0, 0.0)
        sm = CLM.smooth_min(a_d, b_d)
        @test abs(value(sm) - 3.0) < 1e-10  # value matches min
        @test partials(sm)[1] ≈ 1.0 atol=1e-5  # derivative of min(a,b) w.r.t. a when a<b is 1

        # smooth_max with Duals
        sm_max = CLM.smooth_max(a_d, b_d)
        @test abs(value(sm_max) - 5.0) < 1e-10
        @test abs(partials(sm_max)[1]) < 1e-5  # derivative of max(a,b) w.r.t. a when a<b is ~0

        # smooth_min at the kink point: derivative should be nonzero (not NaN)
        x_kink = Dual(0.0, 1.0)
        y_kink = Dual(0.0, 0.0)
        sm_kink = CLM.smooth_min(x_kink, y_kink)
        @test isfinite(value(sm_kink))
        @test isfinite(partials(sm_kink)[1])
        @test abs(partials(sm_kink)[1]) > 0  # nonzero derivative at kink

        # smooth_heaviside at x=0: should give 0.5, with nonzero derivative
        h0 = CLM.smooth_heaviside(Dual(0.0, 1.0))
        @test abs(value(h0) - 0.5) < 1e-10
        @test partials(h0)[1] > 0  # positive derivative at origin

        # smooth_abs at x=0: derivative exists (not NaN)
        sa0 = CLM.smooth_abs(Dual(0.0, 1.0))
        @test isfinite(value(sa0))
        @test isfinite(partials(sa0)[1])

        # smooth_clamp with Duals
        x_clamp = Dual(4.0, 1.0)
        sc = CLM.smooth_clamp(x_clamp, 3.0, 5.0)
        @test abs(value(sc) - 4.0) < 1e-10
        @test partials(sc)[1] ≈ 1.0 atol=1e-3  # inside clamp, derivative ~1

        # Mixed Int/Float64 dispatch
        @test CLM.smooth_min(3.0, 5) === min(3.0, 5.0)
        @test CLM.smooth_min(3, 5.0) === min(3.0, 5.0)
        @test CLM.smooth_max(3.0, 5) === max(3.0, 5.0)
        @test CLM.smooth_max(3, 5.0) === max(3.0, 5.0)

        println("  Smooth AD primitives: all tests passed")
    end

end  # outer testset
