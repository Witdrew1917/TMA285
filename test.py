import numpy as np
from Projekt1 import DifferenceScheeme

def test_populate_sys():
    
    print("test populate_sys():")
    scheeme = DifferenceScheeme.Scheeme(6,6,1)
    scheeme._d = np.ones(5)
    scheeme._a = np.cumsum(np.ones(5))
    scheeme.populate_sys()
    ref1 = np.array([[ 2.,  -0.5,  0.,   0.,   0. ], \
                [-1.,   3.,  -1.,   0.,   0., ], \
                [ 0.,  -1.5,  4.,  -1.5,  0. ], \
                [ 0.,   0.,  -2.,   5. , -2. ], \
                [ 0.,   0.,   0.,  -2.5,  6. ]])
    
    ref2 = np.array([[ 0.,   0.5,  0.,   0.,   0. ], \
                [ 1.,  -1.,   1.,   0.,   0. ],\
                [ 0.,   1.5, -2.,   1.5,  0. ],\
                [ 0.,   0.,   2.,  -3.,   2. ],\
                [ 0.,   0.,   0.,   2.5, -4. ]])
    
    assert(np.sum(ref1 - scheeme._A_plus) <= 0.0001)
    assert(np.sum(ref2 - scheeme._A_minus) <= 0.0001)
    print("---------------------------------------- [Passed]")
    

if __name__ == '__main__':
    
    test_populate_sys()