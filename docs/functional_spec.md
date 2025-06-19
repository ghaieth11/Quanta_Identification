# Functional Specification Document

**Project:** Qubit Parameter Identification\
**Author:** *ALOUI Ghaieth*\
**Date:** 2025-05-28\
**Version:** 1.0

------------------------------------------------------------------------

## 1. Purpose

This system aims to estimate unknown physical parameters of a quantum two-level system (qubit) by simulating its dynamics and statistically analyzing measurement outcomes. The simulation mimics the behavior of a controlled open quantum system and enables extraction of key physical constants such as:

-   Natural frequency (ω)
-   Control interaction strength (κ)
-   Relaxation rate (γ₁)
-   Dephasing rate (γ₂)

------------------------------------------------------------------------

## 2. Functional Overview

### Key Functional Goals:

| Feature               | Description                                                               |
|-----------------------|---------------------------------------------------------------------------|
| Qubit modeling        | Create and manage physical parameters of a qubit                          |
| Control interface     | Apply control fields $u(t)$                                               |
| System dynamics       | Evolve the qubit’s Bloch vector in time using differential equations      |
| Observation process   | Simulate projective measurements using quantum observables                |
| Parameter estimation  | Statistically infer γ₁, γ₂, κ, and ω from repeated simulations            |
| Elimination algorithm | Reduce ambiguity in parameter space using simulated statistical filtering |

------------------------------------------------------------------------

## 3. Workflow / System Behavior

### Step-by-step usage:

1.  **Define a Qubit**
    -   Random or user-specified physical parameters.
2.  **Initialize a System**
    -   Specify Bloch vector coordinates (x, y, z).
3.  **Apply Control**
    -   Input control signal value `u` affecting Hamiltonian dynamics.
4.  **Simulate Evolution**
    -   Use the Bloch equation $\frac{dv}{dt} = Jv + b$ to evolve the qubit.
5.  **Measure**
    -   Project onto observable $E_1 = |1\rangle \langle 1|$ and simulate measurement outcomes.
6.  **Estimate Parameters**
    -   Use custom estimator methods to infer unknown values.

------------------------------------------------------------------------

## 4. Functional Requirements

### 4.1 Qubit Management

| Function         | Description                                     |
|------------------|-------------------------------------------------|
| `Qubit.random()` | Returns a Qubit with random physical parameters |
| `Qubit.copy()`   | Creates a new identical instance                |
| `get_param()`    | Returns list of `[ω, κ, γ₁, γ₂]`                |
| `set_param(...)` | Updates parameters individually or in groups    |

------------------------------------------------------------------------

### 4.2 Control Management

| Function   | Description                    |
|------------|--------------------------------|
| `get_u()`  | Returns current control signal |
| `set_u(u)` | Updates control signal         |

------------------------------------------------------------------------

### 4.3 System Dynamics

| Function            | Description                                                          |
|---------------------|----------------------------------------------------------------------|
| `evolve()`          | Evolves the Bloch state under a given control and qubit for time `t` |
| `get_coordinates()` | Returns current state [x, y, z]                                      |
| `initialize()`      | Resets system to original state                                      |

------------------------------------------------------------------------

### 4.4 Observation

| Function    | Description                                                        |
|-------------|--------------------------------------------------------------------|
| `measure()` | Simulates `n` measurements using the density matrix and observable |

------------------------------------------------------------------------

### 4.5 Estimation

| Method              | Parameter Estimated | Logic Description                                |
|---------------------|---------------------|--------------------------------------------------|
| `estimate_gamma1()` | γ₁                  | Based on decay of excited state                  |
| `estimate_gamma2()` | γ₂                  | Using two-point measurements with log expression |
| `estimate_kappa()`  | κ                   | Uses elimination algorithm with candidate list   |
| `estimate_omega()`  | ω                   | Uses rotation-based elimination algorithm        |

------------------------------------------------------------------------

##  5. Simulation Behavior

The quantum state is represented by a **Bloch vector** $v = (x, y, z)$ and evolves according to:

$$
\frac{dv}{dt} = Jv + b
$$

Where $J$ and $b$ depend on system parameters and control:

-   $J$: 3×3 real matrix based on ω, u, κ, γ₁, γ₂\
-   $b = [0, 0, γ₁]$: external driving vector

The solution is computed numerically via `scipy.integrate.solve_ivp`.

------------------------------------------------------------------------

##  6. Measurement Simulation

The quantum density matrix is:

$$
\rho = \frac{1}{2} (I + x \sigma_x + y \sigma_y + z \sigma_z)
$$

The probability of measurement outcome `1` is:

$$
p = \text{Tr}(E_1 \rho), \quad E_1 = |1\rangle \langle 1|
$$

Measurements are simulated using:

$$
\texttt{np.random.binomial(n=1, p=p, size=n)}
$$

------------------------------------------------------------------------

## 7. Estimation via Elimination Algorithms

**Why needed:** Ambiguity in inverse trigonometric functions during estimation (due to ±θ + 2kπ).

### Elimination Logic:

-   Generate candidate values (e.g., for κ or ω)
-   Simulate system under each candidate
-   Compare output with "true" outcome
-   Eliminate far-off candidates
-   Repeat until convergence or declare null

------------------------------------------------------------------------

# 8. Configuration Parameters

| Parameter | Description                   | Default Range         |
|-----------|-------------------------------|-----------------------|
| `ω`       | Natural frequency             | [0.5, 5.0]            |
| `κ`       | Control interaction strength  | [0.1, 2.0]            |
| `γ₁`      | Relaxation rate               | [0.1, 1.0]            |
| `γ₂`      | Dephasing rate                | [0.01, 0.5]           |
| `u`       | Control signal                | User-defined or fixed |
| `t`       | Simulation duration           | Typically 1.0 or 2.0  |
| `n`       | Number of measurement samples | 1000+ recommended     |

------------------------------------------------------------------------
