import numpy as np
import math

class Scheeme:

    def __init__(self, function: str, time_length : int, space_width : int):

        assert time_length == space_width, \
            f"Expected equal sizes of the space partition and the time partition, got \
                {time_length}, {space_width} "

        self._zdim = space_width
        self._tdim = time_length
        self._U_list = []
        self._A_plus = np.zeros(shape=(self._zdim, self._zdim))
        self._A_minus = np.zeros(shape=(self._zdim, self._zdim))
        self._a = np.zeros(shape=(1,self._zdim))[0]
        self._u = np.zeros(shape=(self._zdim, 1))
        self._d = 0 

    def populate_sys(self):
        
        # TODO 
        # increase definition of A to encompase A(-1/2*d)^-1A(1/2*d)
        # there should also be a matrix for the coefficients a_ij

        temp = np.zeros(shape=(self._zdim, self._zdim))
        temp += np.diag(-1*self._a, 0)
        temp += np.diag(self._a[1:], -1)
        temp += np.diag(-1*self._a[:(len(self._a)-1)], 1)
        temp *= self._d
        self._A_plus = temp.copy()
        self._A_minus = -temp.copy()
        self._A_plus += np.diag(np.ones(self._zdim),0) 
        self._A_minus += np.diag(np.ones(self._zdim),0) 

        return 0
    
    def create_boundry(self):
        

        
        return 0


    def calculate(self, S:float, r:float, sigma:float, K:float, T:float):
        
        time = 0
        for t in range(self._tdim):
            
            if t == 0:
                continue

            # TODO 
            # provide calculation for aquiring Z
            # Z = ...
            Z = np.random.rand(self._zdim) # temporary assignment
            self._d = T * self._zdim**2/ (Z**2 * self._tdim)
            gamma = (1 - math.exp(-r*time))/(r*T)*np.ones(self._zdim)
            self._a = (sigma**2 / 2) * (gamma - Z)
            self.populate_sys()
            # self._u[1:,t+1] = self._A*self._u[1:,t]   
            time += (T / self._tdim)
        
        # print(Z)
        # print(gamma)
        # print(gamma - Z)
        # print(self._d)
        # print(self._a)
        # print(self._A_plus)
        # print(self._A_minus)



