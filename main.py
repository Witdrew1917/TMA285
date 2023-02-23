import numpy as np
from Projekt1 import DifferenceScheeme


if __name__ == '__main__':
    
    f = "function_string"

    arr = np.zeros(shape=(3,2))
    print(arr.ndim)
    iterations = 100
    scheeme = DifferenceScheeme.Scheeme(f,5,5)
    # scheeme.run(iterations)