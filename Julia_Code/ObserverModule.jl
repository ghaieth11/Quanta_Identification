
using LinearAlgebra
using Random
using Distributions

include("ControlModule.jl")
include("QubitModule.jl")
include("SystemModule.jl")
include("UtilsModule.jl")

export Observer, measure, get_obs, set_obs

mutable struct Observer
    obs::Matrix{ComplexF64}
end

get_obs(observer::Observer) = observer.obs

function set_obs!(observer::Observer, new_obs::Matrix{ComplexF64})
    observer.obs = new_obs
end

function measure(observer::Observer, system::QuantumSystem, qubit::Qubit, control::Control, t::Float64, n::Int)::Float64
    evolve!(system, qubit, control, t)
    x, y, z = get_coordinates(system)

    ρ = 0.5 * (Matrix{ComplexF64}(I, 2, 2) + x * sx + y * sy + z * sz)
    p = real(tr(observer.obs * ρ))

    results = rand(Binomial(1, p), n)
    return sum(results) / n
end
