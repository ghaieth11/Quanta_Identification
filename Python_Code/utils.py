import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm

sx = np.array([[0, 1], [1, 0]])
sy = np.array([[0, -1j], [1j, 0]])
sz = np.array([[1, 0], [0, -1]])
paulis = [sx, sy, sz]
E1 = np.array([[0, 0], [0, 1]])

BlochStates = {
    "|0⟩": np.array([0, 0, 1]),
    "|1⟩": np.array([0, 0, -1]),
    "|+⟩": np.array([1, 0, 0]),
    "|-⟩": np.array([-1, 0, 0]),
    "|i⟩": np.array([0, 1, 0]),
    "|-i⟩": np.array([0, -1, 0]),
}

def EliminationAlgorithmKappa(kappa_list: list, n: int, t1: float, parameters: list, v0: np.array) -> float:
    tolerance = 0.15

    while len(kappa_list) > 1:
        # Random new measurement time
        t_new = np.random.uniform(t1 / 2, t1)

        # Evaluate simulation results for all kappa candidates
        results = []
        for kappa_potential in kappa_list:
            new_params = parameters.copy()
            new_params[2] = kappa_potential
            res = QubitSimulation(new_params, v0, t_new, n)
            results.append(res)

        # True result using actual kappa (assumed known from original parameters)
        true_result = QubitSimulation(parameters, v0, t_new, n)

        # Filter candidates based on how close they are to the true result
        final_list = []
        for candidate_kappa, res_candidate in zip(kappa_list, results):
            if abs(res_candidate - true_result) <= tolerance:
                final_list.append(candidate_kappa)

        kappa_list = np.array(final_list)
        if len(kappa_list) == 0:
            break

        tolerance /= 2  # Reduce tolerance
        t1 = t_new

    # Final check
    if len(kappa_list) == 1 and abs(kappa_list[0] - parameters[2]) <= 0.012:
        return kappa_list[0]
    else:
        return "Null"
    
    
    
    
def QubitSimulation(parameters: list, v0: np.array, t: float, n: int) -> float:
    omega, u, kappa, gamma1, gamma2 = parameters

    # Define the matrix J
    J = np.array([
        [-gamma1 / 2 - 2 * gamma2, -omega,           0],
        [ omega,                  -gamma1 / 2 - 2 * gamma2, -u * kappa],
        [    0,                    u * kappa,       -gamma1]
    ])

    # Define vector b
    b = np.array([0, 0, gamma1])

    # Solve the differential equation dv/dt = Jv + b
    sol = solve_ivp(lambda tau, v: J @ v + b, [0, t], v0, t_eval=[t])
    v_t = sol.y[:, -1]
    
    # Construct the density matrix using the Pauli basis
    rho = 0.5 * (np.eye(2) + v_t[0]*sx + v_t[1]*sy + v_t[2]*sz)

    # Projector onto excited state |1><1|
    O = np.array([[0, 0], [0, 1]])

    # Compute the probability of measuring the excited state
    p = np.real(np.trace(O @ rho))

    # Simulate n measurements using a binomial distribution
    results = np.random.binomial(n=1, p=p, size=n)
    return np.sum(results) / n



def EliminationAlgorithmOmega2(kappa_list: list, n: int, t1: float, parameters: list, v0: np.array) -> float:
    tolerance = 0.15
    while len(kappa_list) > 1:
        # Random new measurement time
        t_new = np.random.uniform(t1 / 2, t1)

        # Evaluate simulation results for all kappa candidates
        results = []
        for kappa_potential in kappa_list:
            new_params = parameters.copy()
            new_params[0] = kappa_potential
            res = simulate_measurement(new_params, v0, t_new, n)
            results.append(res)

        # True result using actual kappa (assumed known from original parameters)
        true_result = simulate_measurement(parameters, v0, t_new, n)

        # Filter candidates based on how close they are to the true result
        final_list = []
        for candidate_kappa, res_candidate in zip(kappa_list, results):
            if abs(res_candidate - true_result) <= tolerance:
                final_list.append(candidate_kappa)

        kappa_list = np.array(final_list)
        if len(kappa_list) == 0:
            break

        tolerance /= 2  # Reduce tolerance
        t1 = t_new

    # Final check
    if len(kappa_list) == 1:
        return kappa_list[0]
    else:
        return "Null"
    
    
    
def solve_time(params, v0, t):
    omega, u, kappa, gamma1, gamma2 = params
    J = np.array([
        [-gamma1 / 2 - 2 * gamma2, -omega, 0],
        [ omega, -gamma1 / 2 - 2 * gamma2, -u * kappa],
        [0, u * kappa, -gamma1]
    ])
    b = np.array([0, 0, gamma1])
    sol = solve_ivp(lambda tau, v: J @ v + b, [0, t], v0, t_eval=[t])
    return sol.y[:, -1]  # Return final state

# --- Helper: simulate one measurement ---
def simulate_measurement(params, v_init, t,n):
    v = solve_time(params, v_init, t)

    # Rotate around X by π/4
    theta = np.pi / 4
    UX = expm(1j * theta * sx)
    rho = 0.5 * (np.eye(2) + v[0] * sx + v[1] * sy + v[2] * sz)
    rho_rot = UX.conj().T @ rho @ UX

    # Measure projector onto |1⟩
    O = np.array([[0, 0], [0, 1]])
    p = np.real(np.trace(O @ rho_rot))
    results = np.random.binomial(n=1, p=p, size=n)
    return np.sum(results) / n  # empirical mean
