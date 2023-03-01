import numpy as np
from Projekt1 import DifferenceScheeme


if __name__ == '__main__':
    
    f = "function_string"
    time_samples = 5
    T = 1
    dt = T/time_samples

    scheeme = DifferenceScheeme.Scheeme(f,time_samples,5,dt)
    scheeme.solve_PDE(1,1,1,1,1)
    scheeme.plot(1,1)
        
    # scheeme.run(iterations)