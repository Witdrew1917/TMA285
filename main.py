import numpy as np
from Projekt1 import DifferenceScheeme


if __name__ == '__main__':
    
    f = "function_string"
    time_samples = 20
    T = 20   
    dt = T/time_samples

    scheeme = DifferenceScheeme.Scheeme(f,time_samples,time_samples,dt)
    scheeme.solve_PDE(1,1,1,1,T)    
    scheeme.plot(T)
    
    # scheeme.run(iterations)