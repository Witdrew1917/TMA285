import numpy as np
from Projekt1 import DifferenceScheeme


if __name__ == '__main__':
    
    f = "function_string"

    scheeme = DifferenceScheeme.Scheeme(f,5,5)
    scheeme.calculate(0,1,1,0,1)
    # scheeme.run(iterations)