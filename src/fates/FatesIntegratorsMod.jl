# FatesIntegratorsMod.jl
# Julia port of FATES src/fates/main/FatesIntegratorsMod.F90
#
# Generic ODE integrators: Euler and Runge-Kutta-Fehlberg 4/5 adaptive.
# The derivative function has the signature
#     DerivFunction(Y, Ymask, x, param_array) -> dYdx
# where Y, dYdx, param_array are vectors and Ymask is a Bool vector.
# Deps: FatesConstantsMod.

const max_states = 20

"""
    RKF45(DerivFunction, Y, Ymask, dx, x, max_err, param_array, Yout)

Runge-Kutta-Fehlberg 4/5 order adaptive explicit integration.

`Yout` is updated in place with the 5th-order solution. Returns
`(opt_dx, l_pass)`: the suggested optimal next step size and whether the step
was accepted (estimated error within `max_err`).
"""
function RKF45(DerivFunction, Y::AbstractVector, Ymask::AbstractVector{Bool},
               dx::Real, x::Real, max_err::Real, param_array::AbstractVector,
               Yout::AbstractVector)

    min_step_fraction = 0.25

    t1 = 1.0 / 4.0
    f1_0 = 1.0 / 4.0

    t2 = 3.0 / 8.0
    f2_0 = 3.0 / 32.0
    f2_1 = 9.0 / 32.0

    t3 = 12.0 / 13.0
    f3_0 = 1932.0 / 2197.0
    f3_1 = -7200.0 / 2197.0
    f3_2 = 7296.0 / 2197.0

    t4 = 1.0
    f4_0 = 439.0 / 216.0
    f4_1 = -8.0
    f4_2 = 3680.0 / 513.0
    f4_3 = -845.0 / 4104.0

    t5 = 0.5
    f5_0 = -8.0 / 27.0
    f5_1 = 2.0
    f5_2 = -3544.0 / 2565.0
    f5_3 = 1859.0 / 4104.0
    f5_4 = -11.0 / 40.0

    y_0 = 25.0 / 216.0
    y_2 = 1408.0 / 2565.0
    y_3 = 2197.0 / 4104.0
    y_4 = -1.0 / 5.0

    z_0 = 16.0 / 135.0
    z_2 = 6656.0 / 12825.0
    z_3 = 28561.0 / 56430.0
    z_4 = -9.0 / 50.0
    z_5 = 2.0 / 55.0

    nY = length(Y)

    # 0th Step
    Ytemp = Y[1:nY]
    xtemp = x
    K0 = DerivFunction(Ytemp, Ymask, xtemp, param_array)

    # 1st Step
    Ytemp = Y[1:nY] .+ dx .* (f1_0 .* K0)
    xtemp = x + t1 * dx
    K1 = DerivFunction(Ytemp, Ymask, xtemp, param_array)

    # 2nd Step
    Ytemp = Y[1:nY] .+ dx .* (f2_0 .* K0 .+ f2_1 .* K1)
    xtemp = x + t2 * dx
    K2 = DerivFunction(Ytemp, Ymask, xtemp, param_array)

    # 3rd Step
    Ytemp = Y[1:nY] .+ dx .* (f3_0 .* K0 .+ f3_1 .* K1 .+ f3_2 .* K2)
    xtemp = x + t3 * dx
    K3 = DerivFunction(Ytemp, Ymask, xtemp, param_array)

    # 4th Step
    Ytemp = Y[1:nY] .+ dx .* (f4_0 .* K0 .+ f4_1 .* K1 .+ f4_2 .* K2 .+ f4_3 .* K3)
    xtemp = x + t4 * dx
    K4 = DerivFunction(Ytemp, Ymask, xtemp, param_array)

    # 5th Step
    Ytemp = Y[1:nY] .+ dx .* (f5_0 .* K0 .+ f5_1 .* K1 .+ f5_2 .* K2 .+
                              f5_3 .* K3 .+ f5_4 .* K4)
    xtemp = x + t5 * dx
    K5 = DerivFunction(Ytemp, Ymask, xtemp, param_array)

    # Evaluate error on the 4/5 steps
    # 4th order
    Y4 = Y[1:nY] .+ dx .* (y_0 .* K0 .+ y_2 .* K2 .+ y_3 .* K3 .+ y_4 .* K4)
    # 5th order
    Yout[1:nY] .= Y[1:nY] .+ dx .* (z_0 .* K0 .+ z_2 .* K2 .+ z_3 .* K3 .+
                                    z_4 .* K4 .+ z_5 .* K5)

    # Take the maximum absolute error across all variables
    err45 = maximum(abs.(Yout[1:nY] .- Y4))

    # Update our estimate of the optimal time-step.
    opt_dx = dx * max(min_step_fraction,
                      0.840896 * (max_err / max(err45, 0.00001 * max_err))^0.25)

    l_pass = !(err45 > max_err)

    return opt_dx, l_pass
end

# ===================================================================================

"""
    Euler(DerivFunction, Y, Ymask, dx, x, param_array, Yout)

Simple Euler integration. `Yout` is updated in place.
"""
function Euler(DerivFunction, Y::AbstractVector, Ymask::AbstractVector{Bool},
               dx::Real, x::Real, param_array::AbstractVector, Yout::AbstractVector)

    nY = length(Y)
    dYdx = DerivFunction(Y[1:nY], Ymask, x, param_array)
    Yout[1:nY] .= Y[1:nY] .+ dx .* dYdx
    return nothing
end
