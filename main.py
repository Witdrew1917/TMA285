import numpy as np
from Projekt1 import DifferenceScheeme


if __name__ == '__main__':
    
    f = "function string"
    iterations = 100
    scheeme = DifferenceScheeme(f)
    scheeme.run(iterations)