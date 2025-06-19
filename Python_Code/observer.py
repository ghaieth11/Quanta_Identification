from qubit import Qubit
from control import Control
from system import System
from utils import sx, sy, sz
import numpy as np

class Observer:
    def __init__(self, obs):
        self.obs = obs  # Observable, e.g., np.array([[0, 0], [0, 1]])

    def get_obs(self):
        return self.obs

    def set_obs(self, o):
        self.obs = o

    def measure(self, system: System, qubit: Qubit, control: Control, t: float, n: int) -> float:
        # Evolve the system's state to time t
        system.evolve(qubit, control, t)
        x, y, z = system.get_coordinates()

        # Construct the density matrix using the Pauli basis
        rho = 0.5 * (np.eye(2) + x * sx + y * sy + z * sz)

        # Compute the probability of measuring the observable
        p = np.real(np.trace(self.obs @ rho))

        # Simulate n measurements using a binomial distribution
        results = np.random.binomial(n=1, p=p, size=n)
        return np.sum(results) / n