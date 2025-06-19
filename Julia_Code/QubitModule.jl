
using Random
using Distributions

export Qubit, get_param, set_param

mutable struct Qubit
    omega::Float64
    kappa::Float64
    gamma1::Float64
    gamma2::Float64
end

# Constructor from values (déjà généré automatiquement par Julia)


# Constructor from another Qubit
function Base.copy(q::Qubit)
    return Qubit(q.omega, q.kappa, q.gamma1, q.gamma2)
end

# Random constructor
function Qubit_random()
    omega = rand(Uniform(0.5, 5))
    kappa = rand(Uniform(0.1, 2))
    gamma1 = rand(Uniform(0.1, 1))
    gamma2 = rand(Uniform(0.01, 0.5))
    return Qubit(omega, kappa, gamma1, gamma2)
end

# String representation
Base.show(io::IO, q::Qubit) = print(io,
    "Qubit(omega=$(round(q.omega, digits=5)), kappa=$(round(q.kappa, digits=5)), " *
    "gamma1=$(round(q.gamma1, digits=5)), gamma2=$(round(q.gamma2, digits=5)))"
)

# Getters
get_omega(q::Qubit) = q.omega
get_kappa(q::Qubit) = q.kappa
get_gamma1(q::Qubit) = q.gamma1
get_gamma2(q::Qubit) = q.gamma2
get_param(q::Qubit) = [q.omega, q.kappa, q.gamma1, q.gamma2]

# Setter (with optional arguments)
function set_param!(q::Qubit; omega=nothing, kappa=nothing, gamma1=nothing, gamma2=nothing)
    if omega !== nothing; q.omega = omega; end
    if kappa !== nothing; q.kappa = kappa; end
    if gamma1 !== nothing; q.gamma1 = gamma1; end
    if gamma2 !== nothing; q.gamma2 = gamma2; end
    return q
end
