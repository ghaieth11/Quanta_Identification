# Load your project files
include("UtilsModule.jl")
include("QubitModule.jl")
include("ControlModule.jl")
include("SystemModule.jl")
include("ObserverModule.jl")
include("EstimatorModule.jl")

# Create instances
q = Qubit_random()
c = Control(2.0)

# Gamma1 Estimation
S1 = QuantumSystem(0.0, 0.0, -1.0)
E = Estimator(S1)
gam = estimate_gamma1(E, q, 100_000)
println("estimated γ₁ = $gam, true γ₁ = $(get_gamma1(q))")

# Kappa Estimation
initialize!(S1)
E = Estimator(S1)
kappa = estimate_kappa(E, q, 100_000)
println("estimated κ = $kappa, true κ = $(get_kappa(q))")

# Gamma2 Estimation
S2 = QuantumSystem(0.0, 1.0, 0.0)
E2 = Estimator(S2)
gamma2 = estimate_gamma2(E2, q, 100_000)
println("estimated γ₂ = $gamma2, true γ₂ = $(get_gamma2(q))")

# Omega Estimation
initialize!(S2)
E = Estimator(S2)
omega = estimate_omega(E, q, 100_000)
println("estimated ω = $omega, true ω = $(get_omega(q))")