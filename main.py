import numpy as np
from Projekt1 import DifferenceScheeme


if __name__ == '__main__':
    
    f = "function_string"
    time_samples = 20
    T = 1/2  
    dt = T/time_samples

    scheeme = DifferenceScheeme.Scheeme(time_samples,time_samples,dt)
    scheeme.solve_PDE(100,0.05,0.5,100,T)    
    scheeme.plot(T)
    
    # scheeme.run(iterations)Z[self._zdim - 1]