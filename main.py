import numpy as np
from Projekt1 import DifferenceScheeme


if __name__ == '__main__':
    
    f = "function_string"
    time_samples = 30
    T = 1/2  
    dt = T/time_samples
    S0 = 100
    r = 0.05
    sigma = 0.05     
    K = 100


    scheeme = DifferenceScheeme.Scheeme(time_samples,time_samples,dt)
    scheeme.solve_PDE(S0,r,sigma,K,T)    
    scheeme.plot(T)
    
    # scheeme.run(iterations)Z[self._zdim - 1]