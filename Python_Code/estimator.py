from utils import *
from qubit import Qubit
from control import Control
from system import System
from observer import Observer
import numpy as np
from scipy.linalg import expm


class Estimator:
    def __init__(self, system: System):
        self.system = system
    
    def estimate_gamma1(self, qubit: Qubit, n):
        control = Control(u=0)
        observer = Observer(E1)
        t = 2.0
        p = observer.measure(self.system, qubit, control, t, n)
        gamma1 = -np.log(p) / t
        return gamma1
    
    def estimate_kappa(self, qubit:Qubit, n) : 
        control = Control(u=1000)
        observer = Observer(E1)
        t = 1/control.get_u()
        p = observer.measure(self.system, qubit, control, t, n)
        
        theta = np.arccos(2 * p - 1)
        KAPPA = np.array([val for k in range(-5, 5) for val in [theta + k * np.pi, -theta + k * np.pi] if val > 0])

        # Full parameters list to pass to simulation: omega, u, kappa, gamma1, gamma2
        omega, kappa, gamma1, gamma2 = qubit.get_param()
        params = [omega, control.get_u(), kappa, gamma1, gamma2]

        # Initial state vector
        v0 = self.system.get_coordinates()

    # Call elimination algorithm
        return EliminationAlgorithmKappa(KAPPA, n, t, params, v0)
    
    def estimate_gamma2(self, qubit: Qubit, n: int):
        control = Control(u=0)
        observer = Observer(E1)
        v_init = self.system.get_coordinates()  # Y-axis state
        qubit_params = qubit.get_param()
        
        t = 1.0
        # Simulate measurements at t and 2t
        s1 = self._simulate_rotated(qubit,control,t,n)
        self.system.initialize()
        s2 = self._simulate_rotated(qubit,control, 2*t,n)

        ln_argument = 2 * s2 - 1 + 2 * (1 - 2 * s1) ** 2
        if ln_argument <= 0:
            raise ValueError(f"log argument non-positive: {ln_argument:.5f}")
        
        ln_term = -np.log(ln_argument) / (4 * t)
        gamma2 = ln_term - 0.25 * qubit_params[2]  # subtract gamma1/4
        return gamma2
        
    def estimate_omega(self, qubit: Qubit, n: int):
        control = Control(u=0)
        observer = Observer(E1)
        omega, kappa, gamma1, gamma2 = qubit.get_param()
        params = [omega, control.get_u(), kappa, gamma1, gamma2]
        v0 = self.system.get_coordinates()
        
        t = 1.0

        s1 = self._simulate_rotated(qubit,control,t, n)
        gamma = 0.5 * params[3] + 2 * params[4]

        theta = np.arccos((1 - 2 * s1) * np.exp(gamma * t))
        OMEGA = np.array([
            val / t
            for k in range(-30, 30)
            for val in [theta + k * np.pi, -theta + k * np.pi]
            if val > 0
        ])
        

        return EliminationAlgorithmOmega2(OMEGA, n, t, params, v0)   
        
        
            
    def _simulate_rotated(self,qubit,control,t,n):
        self.system.evolve(qubit,control,t)
        v = self.system.get_coordinates()

        theta = np.pi / 4
        UX = expm(1j * theta * sx)
        rho = 0.5 * (np.eye(2) + v[0]*sx + v[1]*sy + v[2]*sz)
        rho_rot = UX.conj().T @ rho @ UX

        p = np.real(np.trace(E1 @ rho_rot))
        return np.sum(np.random.binomial(1, p, size=n)) / n
    
