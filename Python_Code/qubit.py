import numpy as np


class Qubit : 
    # CONSTRUCTORS
    def __init__(self,omega,kappa,gamma1,gamma2) : 
        self.omega = omega
        self.kappa = kappa
        self.gamma1 = gamma1
        self.gamma2 = gamma2
    
    @classmethod
    def random(cls):
        omega = np.random.uniform(0.5, 5)
        kappa = np.random.uniform(0.1, 2)
        gamma1 = np.random.uniform(0.1, 1)
        gamma2 = np.random.uniform(0.01, 0.5)
        return cls(omega, kappa, gamma1, gamma2)

    @classmethod
    def copy(cls, other_qubit):
        if not isinstance(other_qubit, cls):
            raise TypeError("Can only copy from another Qubit instance.")
        return cls(other_qubit.omega, other_qubit.kappa, other_qubit.gamma1, other_qubit.gamma2)
    
    # PRINT REPRESENTATION
    def __repr__(self):
        return (f"Qubit(omega={self.omega:.5f}, kappa={self.kappa:.5f}, "
                f"gamma1={self.gamma1:.5f}, gamma2={self.gamma2:.5f})")
    
    # GETTERS    
    def get_omega(self) : 
        return self.omega
    def get_kappa(self) : 
        return self.kappa
    def get_gamma1(self) : 
        return self.gamma1
    def get_gamma2(self) : 
        return self.gamma2
    def get_param(self) : 
        return [self.omega, self.kappa, self.gamma1, self.gamma2]
    
    # SETTERS 
    
    def set_param(self, omega = None, kappa = None, gamma1 = None, gamma2 = None) : 
        if omega is not None : self.omega = omega
        if kappa is not None : self.kappa = kappa
        if gamma1 is not None : self.gamma1 = gamma1
        if gamma2 is not None : self.gamma2 = gamma2
    
    