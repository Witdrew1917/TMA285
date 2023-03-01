from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np
import math

class Scheeme:

    def __init__(self, function: str, time_length : int, space_width : int, dt:float):

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
        # self._U_list.append(self._u.copy())
        self._d = np.zeros(self._zdim)
        self._dt = dt

    def populate_sys(self):
        
        # TODO 
        # increase definition of A to encompase A(-1/2*d)^-1A(1/2*d)
        # there should also be a matrix for the coefficients a_ij

        temp = np.zeros(shape=(self._zdim, self._zdim))
        temp += np.diag(-1*self._a, 0)
        temp += np.diag(self._a[1:], -1)
        temp += np.diag(-1*self._a[:(len(self._a)-1)], 1)
        temp *= (np.ones(shape=(len(self._d),len(self._d))) * self._d).T
        self._A_plus = temp.copy()
        self._A_minus = -temp.copy()
        self._A_plus += np.diag(np.ones(self._zdim),0) 
        self._A_minus += np.diag(np.ones(self._zdim),0) 
        self._S = []

        return 0
    
    def price_calulation(self, S0:float, r:float,t:float):
        
        sigma = self._dt
        W = np.random.rand(len(t))
        W = 1/(np.sqrt(2*np.pi*(sigma**2)))*np.exp(-1/2*(W/sigma)**2) # this might not be correct
        
        self._S = S0*np.exp((r - sigma**2/2)*t + sigma*W)

        return 0

    def partition_z(self, S0:float, r:float, K:float, T:float):
        
        # partition t uniformly
        t = np.linspace(0, T, self._tdim + 1)
        self.price_calulation(S0, r, t)
        Q = np.zeros(len(t)).T
        for i in range(len(t)):
            Q[i] = sum(self._dt*self._S[:(i+1)])
        
        Z = []
        Z = 1/(r*T)*(1 - np.exp(-r*(T - t))) + np.matmul(np.exp(-r*(T - t))/self._S,(Q/T - K)) 
        Z.sort()
        dz = np.zeros(len(Z) - 1)
        for i in range(len(Z) - 1):
            dz[i] = Z[i + 1] - Z[i]

        return dz

    def solve_PDE(self, S0:float, r:float, sigma:float, K:float, T:float):
        
        time = 0
        self._u = np.random.rand(self._zdim).T
        self._U_list.append(self._u.copy())
        dz = self.partition_z(S0, r, K, T)
        for t in range(self._tdim):
            
            if t == 0:
                continue

            # TODO 
            # provide calculation for aquiring Z
            # Z = ...
            Z = np.random.rand(self._zdim) # temporary assignment
            
            # should be vector, see equation in chat
            self._d = self._dt*np.ones(len(dz))/(dz**2)
            gamma = (1 - math.exp(-r*time))/(r*T)*np.ones(self._zdim)
            self._a = (sigma**2 / 2) * (gamma - Z)
            self.populate_sys()
            self._u = np.matmul(np.linalg.inv(self._A_plus),np.matmul(self._A_minus,self._u))
            self._U_list.append(self._u.copy())   
            time += self._dt

    def calc_call(self):
        
        i = 0
        U = np.zeros(shape=(self._tdim, self._zdim))
        for u in self._U_list:
            U[i][:] = u
            i += 1

        # what is X in equation 6.43?
        # X = ...
        return self._S*U
       
    def plot(self, T:float, Zmax:float):

        fig = plt.figure()
        ax = plt.axes(projection = '3d')
        x = np.linspace(0,T, self._tdim)
        y = np.linspace(0,Zmax, self._zdim)
        X, Y = np.meshgrid(x,y)
        i = 0
        Z = np.zeros(shape=(self._tdim, self._zdim))
        for u in self._U_list:
            Z[i][:] = u
            i += 1

        print(Z)
        print(X)
        print(Y)
        ax.plot_surface(X, Y, Z, rstride = 1, cstride = 1, cmap='viridis', edgecolor = 'none')
        fig.show()
        input("Press Enter to continue...")
        return 0

