from utils import E1
from qubit import Qubit 
from control import Control
from system import System
from observer import Observer
from estimator import Estimator


q = Qubit.random()
c = Control(2)

#gamma1
S1 = System(0,0,-1)
E = Estimator(S1)
gam = E.estimate_gamma1(q,n=1000000)
print(f"estimated {gam} true {q.get_gamma1()}")

#kappa
S1.initialize()
E = Estimator(S1)
kappa = E.estimate_kappa(q,n=100000)
print(f"estimated {kappa} true {q.get_kappa()}")

#gamma2 
S2 = System(0,1,0)
E2 = Estimator(S2)
gamma2 = E2.estimate_gamma2(q,n=100000)
print(f"estimated {gamma2} true {q.get_gamma2()}")

#omega 
S2.initialize()
E = Estimator(S2)
omega = E2.estimate_omega(q,n=100000)
print(f"estimated {omega} true {q.get_omega()}")