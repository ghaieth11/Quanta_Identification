using LinearAlgebra
using DifferentialEquations: ODEProblem, solve, Tsit5
using Plots


include("QubitModule.jl")
include("ControlModule.jl")
include("UtilsModule.jl")


mutable struct QuantumSystem
    x::Float64
    y::Float64
    z::Float64
    x0::Float64
    y0::Float64
    z0::Float64
end

function QuantumSystem(x::Float64, y::Float64, z::Float64)
    return QuantumSystem(x, y, z, x, y, z)
end

function copy(sys::QuantumSystem)
    return QuantumSystem(sys.x, sys.y, sys.z, sys.x0, sys.y0, sys.z0)
end

get_x(s::QuantumSystem) = s.x
get_y(s::QuantumSystem) = s.y
get_z(s::QuantumSystem) = s.z
get_coordinates(s::QuantumSystem) = [s.x, s.y, s.z]

function initialize!(s::QuantumSystem)
    s.x = s.x0
    s.y = s.y0
    s.z = s.z0
end

function evolve!(s::QuantumSystem, q::Qubit, c::Control, t::Float64)
    ω, κ, γ1, γ2 = get_param(q)
    u = get_u(c)
    v0 = get_coordinates(s)

    J = [
        -γ1 / 2 - 2 * γ2   -ω                0;
         ω                -γ1 / 2 - 2 * γ2  -u * κ;
         0                 u * κ           -γ1
    ]
    b = [0.0, 0.0, γ1]

    function f!(dv, v, p, t)
        dv[:] = J * v .+ b
    end

    prob = ODEProblem(f!, v0, (0.0, t))
    sol = solve(prob, Tsit5(), saveat=t)
    vt = sol.u[end]

    s.x, s.y, s.z = vt
end

function draw(s::QuantumSystem, q::Qubit, c::Control, t::Float64)
    times = range(0, t, length=100)
    ω, κ, γ1, γ2 = get_param(q)
    u = get_u(c)
    v0 = get_coordinates(s)

    J = [
        -γ1 / 2 - 2 * γ2   -ω                0;
         ω                -γ1 / 2 - 2 * γ2  -u * κ;
         0                 u * κ           -γ1
    ]
    b = [0.0, 0.0, γ1]

    function f!(dv, v, p, t)
        dv[:] = J * v .+ b
    end

    prob = ODEProblem(f!, v0, (0.0, t))
    sol = solve(prob, Tsit5(), saveat=times)

    x_vals = [u[1] for u in sol.u]

    plot(times, x_vals, xlabel="Time", ylabel="x(t)", label="x(t) evolution", lw=2)
end
