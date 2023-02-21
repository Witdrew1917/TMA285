import numpy as np
from Projekt1 import DifferenceScheeme


if __name__ == '__main__':
    
    f = "function_string"
    iterations = 100
    scheeme = DifferenceScheeme.Scheeme(f)
    scheeme.run(iterations)