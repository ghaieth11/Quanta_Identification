import numpy as np
from scipy.integrate import solve_ivp
from utils import sx, sy, sz
from qubit import Qubit
from control import Control
import matplotlib.pyplot as plt

class System:
    # CONSTRUCTOR
    def __init__(self, x: float, y: float, z: float):
        self.x = x
        self.y = y
        self.z = z
        self._initial_state = (x, y, z)

    # COPY CONSTRUCTOR
    
    @classmethod
    def copy(cls, other_system):
        if not isinstance(other_system, cls):
            raise TypeError("Can only copy from another System instance.")
        return cls(other_system.x, other_system.y, other_system.z)

    # GETTERS
    def get_x(self):
        return self.x

    def get_y(self):
        return self.y

    def get_z(self):
        return self.z

    def get_coordinates(self):
        return [self.x, self.y, self.z]

    # EVOLUTION METHOD
    def evolve(self, q: Qubit, c: Control, t: float):
        omega, kappa, gamma1, gamma2 = q.get_param()
        u = c.get_u()
        v0 = self.get_coordinates()

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
        v_t = sol.y[:,-1]

        # Update internal state
        self.x, self.y, self.z = v_t.tolist()
        
    def initialize(self):
        """Reset the system's coordinates to their initial values."""
        self.x, self.y, self.z = self._initial_state
        
    def draw(self,q,c,t) :
        X = []
        times = np.linspace(0,t,100)
        for a in times : 
            self.evolve(q,c,float(a))
            X.append(self.x)
            self.initialize()
        X = np.array(X)    
        plt.plot(times,X)
         
        
            