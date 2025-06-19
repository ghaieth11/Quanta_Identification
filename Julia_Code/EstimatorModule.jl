
using LinearAlgebra
using Random
using StatsFuns 

include("QubitModule.jl")
include("ControlModule.jl")
include("ObserverModule.jl")
include("SystemModule.jl")
include("UtilsModule.jl")



mutable struct Estimator
    system::QuantumSystem
end

# -- Estimation de gamma1 --
function estimate_gamma1(est::Estimator, qubit::Qubit, n::Int64)
    control = Control(0.0)
    observer = Observer(E1)
    t = 2.0
    p = measure(observer, est.system, qubit, control, t, n)
    return -log(p) / t
end

# -- Estimation de kappa --
function estimate_kappa(est::Estimator, qubit::Qubit, n::Int)
    control = Control(1000.0)
    observer = Observer(UtilsModule.E1)
    t = 1 / get_u(control)
    p = measure(observer, est.system, qubit, control, t, n)
    
    theta = acos(2p - 1)
    KAPPA = [θ for k in -5:4 for θ in (theta + k * π, -theta + k * π) if θ > 0]

    ω, κ, γ1, γ2 = get_param(qubit)
    params = [ω, get_u(control), κ, γ1, γ2]
    v0 = get_coordinates(est.system)

    return EliminationAlgorithmKappa(KAPPA, n, t, params, v0)
end

# -- Estimation de gamma2 --
function estimate_gamma2(est::Estimator, qubit::Qubit, n::Int)
    control = Control(0.0)
    v_init = get_coordinates(est.system)
    γ1 = get_param(qubit)[3]

    t = 1.0
    s1 = _simulate_rotated(est, qubit, control, t, n)
    initialize!(est.system)
    s2 = _simulate_rotated(est, qubit, control, 2t, n)

    ln_arg = 2s2 - 1 + 2(1 - 2s1)^2
    if ln_arg <= 0
        error("log argument non-positive: $ln_arg")
    end

    ln_term = -log(ln_arg) / (4t)
    return ln_term - 0.25 * γ1
end

# -- Estimation de omega --
function estimate_omega(est::Estimator, qubit::Qubit, n::Int)
    control = Control(0.0)
    ω, κ, γ1, γ2 = get_param(qubit)
    params = [ω, get_u(control), κ, γ1, γ2]
    v0 = get_coordinates(est.system)
    t = 1.0

    s1 = _simulate_rotated(est, qubit, control, t, n)
    γ = 0.5 * γ1 + 2 * γ2
    θ = acos((1 - 2s1) * exp(γ * t))
    OMEGA = [θ_k / t for k in -30:29 for θ_k in (θ + k * π, -θ + k * π) if θ_k > 0]

    return EliminationAlgorithmOmega2(OMEGA, n, t, params, v0)
end

# -- Simulation avec rotation UX --
function _simulate_rotated(est::Estimator, qubit::Qubit, control::Control, t::Float64, n::Int)
    evolve!(est.system, qubit, control, t)
    v = get_coordinates(est.system)

    θ = π / 4
    UX = exp(1im * θ * UtilsModule.sx)
    ρ = 0.5 * (I + v[1] * UtilsModule.sx + v[2] * UtilsModule.sy + v[3] * UtilsModule.sz)
    ρ_rot = UX' * ρ * UX

    p = real(tr(UtilsModule.E1 * ρ_rot))
    return sum(rand(Binomial(1, p), n)) / n
end
