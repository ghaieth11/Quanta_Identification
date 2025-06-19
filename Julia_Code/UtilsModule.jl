module UtilsModule

using LinearAlgebra
using Random
using DifferentialEquations
using Statistics

export sx, sy, sz, E1

# Définition des matrices de Pauli
const sx = [0.0 1.0; 1.0 0.0]
const sy = [0.0 -im; im 0.0]
const sz = [1.0 0.0; 0.0 -1.0]
const paulis = [sx, sy, sz]
const E1 = [0.0 0.0; 0.0 1.0]

# États sur la sphère de Bloch
const BlochStates = Dict(
    "|0⟩" => [0.0, 0.0, 1.0],
    "|1⟩" => [0.0, 0.0, -1.0],
    "|+⟩" => [1.0, 0.0, 0.0],
    "|-⟩" => [-1.0, 0.0, 0.0],
    "|i⟩" => [0.0, 1.0, 0.0],
    "|-i⟩" => [0.0, -1.0, 0.0],
)

# --- Simulation du qubit ---
function QubitSimulation(parameters::Vector{Float64}, v0::Vector{Float64}, t::Float64, n::Int)
    ω, u, κ, γ1, γ2 = parameters

    J = [
        -γ1/2 - 2*γ2   -ω           0;
         ω           -γ1/2 - 2*γ2  -u*κ;
         0            u*κ        -γ1
    ]
    b = [0.0, 0.0, γ1]

    function dyn!(dv, v, p, t)
        dv .= J * v + b
    end

    prob = ODEProblem(dyn!, v0, (0.0, t))
    sol = solve(prob, Tsit5(), saveat=t)
    vt = sol.u[end]

    # Construction de la matrice densité
    ρ = 0.5 * (Matrix{ComplexF64}(I, 2, 2) + vt[1]*sx + vt[2]*sy + vt[3]*sz)

    # Probabilité de mesurer l'état excité |1⟩
    O = E1
    p = real(tr(O * ρ))

    # Simulation de mesures binomiales
    return mean(rand(Binomial(1, p), n))
end

# --- Résolution de l'équation différentielle seule ---
function solve_time(params::Vector{Float64}, v0::Vector{Float64}, t::Float64)
    ω, u, κ, γ1, γ2 = params

    J = [
        -γ1/2 - 2γ2   -ω           0;
         ω           -γ1/2 - 2γ2  -u*κ;
         0            u*κ        -γ1
    ]
    b = [0.0, 0.0, γ1]

    function dyn!(dv, v, p, t)
        dv .= J * v + b
    end

    prob = ODEProblem(dyn!, v0, (0.0, t))
    sol = solve(prob, Tsit5(), saveat=t)
    return sol.u[end]
end

# --- Simulation d’une mesure avec rotation ---
function simulate_measurement(params::Vector{Float64}, v_init::Vector{Float64}, t::Float64, n::Int)
    v = solve_time(params, v_init, t)

    θ = π / 4
    UX = exp(1im * θ * sx)
    ρ = 0.5 * (Matrix{ComplexF64}(I, 2, 2) + v[1]*sx + v[2]*sy + v[3]*sz)
    ρ_rot = UX' * ρ * UX

    O = E1
    p = real(tr(O * ρ_rot))
    return mean(rand(Binomial(1, p), n))
end

# --- Algorithme d'élimination sur κ ---
function EliminationAlgorithmKappa(kappa_list::Vector{Float64}, n::Int, t1::Float64, parameters::Vector{Float64}, v0::Vector{Float64})
    tol = 0.15
    while length(kappa_list) > 1
        t_new = rand() * (t1 / 2) + (t1 / 2)

        results = Float64[]
        for κ in kappa_list
            new_params = copy(parameters)
            new_params[3] = κ  # param[3] == κ
            push!(results, QubitSimulation(new_params, v0, t_new, n))
        end

        true_result = QubitSimulation(parameters, v0, t_new, n)

        # Filtrage des candidats
        filtered = [κ for (κ, r) in zip(kappa_list, results) if abs(r - true_result) <= tol]
        kappa_list = filtered
        if isempty(kappa_list)
            break
        end

        tol /= 2
        t1 = t_new
    end

    if length(kappa_list) == 1 && abs(kappa_list[1] - parameters[3]) <= 0.012
        return kappa_list[1]
    else
        return nothing
    end
end

# --- Algorithme d'élimination sur ω (EliminationAlgorithmOmega2) ---
function EliminationAlgorithmOmega2(omega_list::Vector{Float64}, n::Int, t1::Float64, parameters::Vector{Float64}, v0::Vector{Float64})
    tol = 0.15
    while length(omega_list) > 1
        t_new = rand() * (t1 / 2) + (t1 / 2)

        results = Float64[]
        for ω in omega_list
            new_params = copy(parameters)
            new_params[1] = ω  # param[1] == ω
            push!(results, simulate_measurement(new_params, v0, t_new, n))
        end

        true_result = simulate_measurement(parameters, v0, t_new, n)

        filtered = [ω for (ω, r) in zip(omega_list, results) if abs(r - true_result) <= tol]
        omega_list = filtered
        if isempty(omega_list)
            break
        end

        tol /= 2
        t1 = t_new
    end

    if length(omega_list) == 1
        return omega_list[1]
    else
        return nothing
    end
end

end 