from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np
import math

class Scheeme:

    def __init__(self, time_length : int, space_width : int, dt:float):

        assert time_length == space_width, \
            f"Expected equal sizes of the space partition and the time partition, got \
                {time_length}, {space_width} "

        self._zdim = space_width
        self._tdim = time_length
        self._U_list = []
        self._A_plus = np.zeros(shape=(self._zdim - 1, self._zdim - 1))
        self._A_minus = np.zeros(shape=(self._zdim - 1, self._zdim - 1))
        self._a = np.zeros(shape=(1,self._zdim - 1))[0]
        self._u = np.zeros(shape=(self._zdim - 1, 1))
        self._d = np.zeros(self._zdim - 1)
        self._dt = dt
        self._sort_order_call = []
        self._sort_order_put = []
        self._S = []
        self._Q = []
        self._call_put = []
        self._parity = []


    def populate_sys(self):
        
        # TODO 
        # increase definition of A to encompase A(-1/2*d)^-1A(1/2*d)
        # there should also be a matrix for the coefficients a_ij

        temp = np.zeros(shape=(self._zdim - 1, self._zdim - 1))
        temp += np.diag(self._a, 0)
        temp += np.diag(-1/2*self._a[1:], -1)
        temp += np.diag(-1/2*self._a[:(len(self._a)-1)], 1)
        temp *= (np.ones(shape=(len(self._d),len(self._d))) * self._d).T
        self._A_plus = temp.copy()
        self._A_minus = -temp.copy()
        self._A_plus += np.diag(np.ones(self._zdim - 1),0) 
        self._A_minus += np.diag(np.ones(self._zdim - 1),0) 

        return 0
    
    
    def calc_gamma(self, time:float, r:float, T:float) -> tuple[float, float]:

        gamma = (1 - math.exp(-r*time))/(r*T)*np.ones(self._zdim)
        gamma_next = (1 - math.exp(-r*(time + self._dt)))/(r*T)*np.ones(self._zdim)

        return [gamma, gamma_next]
    

    def price_calulation(self, S0:float, r:float, sigma:float, t:float):
        
        W = np.random.normal(0, sigma, len(t))
        q = np.ones(len(t))
        self._S = S0*np.exp((r - sigma**2/2)*self._dt*np.cumsum(q) \
                            + sigma*np.sqrt(self._dt)*np.cumsum(W))

        return 0


    def partition_z(self, S0:float, r:float, sigma:float, K:float, T:float) -> tuple[np.ndarray, np.ndarray]:
        
        # partition t uniformly
        t = np.linspace(0, T, self._tdim)
        self.price_calulation(S0, r, sigma, t)
        Q = np.zeros(len(t)).T
        for i in range(len(t)):
            Q[i] = sum(self._dt*self._S[:(i+1)])

        self._Q = Q

        Z_call = []
        Z_call = 1/(r*T)*(1 - np.exp(-r*(T - t))) + (np.exp(-r*(T - t))/self._S *(Q/T - K)) 
        self._sort_order_call = np.argsort(Z_call)
        Z_call.sort()
        dz_call = np.zeros(len(Z_call) - 1)
        for i in range(len(Z_call) - 1):
            dz_call[i] = Z_call[i + 1] - Z_call[i]

        Z_put = []
        Z_put = 1/(r*T)*(1 - np.exp(-r*(T - t))) + (np.exp(-r*(T - t))/self._S * (K - Q/T)) 
        self._sort_order_put = np.argsort(Z_put)
        Z_put.sort()
        dz_put = np.zeros(len(Z_put) - 1)
        for i in range(len(Z_put) - 1):
            dz_put[i] = Z_put[i + 1] - Z_put[i]

        print(Z_call)
        print(Z_put)
        Z = []
        Z.append(Z_call)
        Z.append(Z_put)
        
        dz = []
        dz.append(dz_call)
        dz.append(dz_put)

        return [Z, dz]


    def solve_PDE(self, S0:float, r:float, sigma:float, K:float, T:float):
        
        time = 0
        Z_arr ,dz_arr = self.partition_z(S0, r, sigma, K, T)
        args = ["call", "put"]
        for Z, dz, arg in zip(Z_arr, dz_arr, args):
            
            self._u = zero_max(Z[:(len(Z)-1)])
            self._U_list.clear()
            self._U_list.append(self._u.copy())
            
            for t in np.linspace(0,T,self._tdim):
                
                if t == 0:
                    continue

                # TODO 
                # provide calculation for aquiring Z -> you've done that 
                # Z = ...
                # should be vector, see equation in chat
                self._d = self._dt*np.ones(len(dz))/(dz**2)
                gamma, gamma_next = self.calc_gamma(time, r, T)
                
                temp = ((sigma**2 / 2) * (gamma - Z))
                self._a = temp[:(self._zdim - 1)]
                a_n = ((sigma**2 / 2) * (gamma_next - Z))[self._zdim - 1]
                
                self.populate_sys()

                correction = np.zeros(self._zdim - 1).T
                correction[self._zdim - 2] = \
                    (1/2)*self._d[len(dz) - 1]*Z[self._zdim - 1]*(a_n* + temp[self._zdim - 1]) 
                
                self._u = np.matmul(np.linalg.inv(self._A_plus),np.matmul(self._A_minus,self._u) + correction) 
                self._U_list.append(self._u.copy())   
                time += self._dt
            

            self._call_put.append(self.calc_call_put(arg))

        self.calc_parity(r, K, T, np.linspace(0,T,self._tdim))

    def calc_call_put(self, arg:str) -> np.ndarray:
        
        i = 0
        U = np.zeros(shape=(self._tdim, self._zdim - 1))
        for u in self._U_list:
            U[i][:] = u
            i += 1
        
        temp = np.zeros(self._tdim)
        for i in range(len(temp) - 1):
            temp[i] = U[i][i]

        unsorted_u = np.zeros(self._tdim)
        it = 0
        if arg == "call":
            for j in self._sort_order_call:
                unsorted_u[it] = temp[j]
                it += 1
        
        if arg == "put":
            for j in self._sort_order_put:
                unsorted_u[it] = temp[j]
                it += 1
        
        
        return self._S*unsorted_u

    def calc_parity(self, r:float, K: float, T:float, t):

        self._parity = (1/T)*self._Q*np.exp(-r*(T - t))  \
            + self._S/(r*T)*(1 - np.exp(-r*(T - t)) - K*np.exp(-r*(T - t)))
        return 0


    def plot(self, T:float, ):

        x = self._call_put[0]
        y = self._call_put[1]
        t = np.linspace(0,T, self._tdim)
        # plt.plot(t,(x - y) - self._parity)
        plt.plot(t,x)
        # plt.plot(t, y)
        plt.show()
        return 0

def zero_max(arr:np.ndarray) -> np.ndarray:
    """
        Input: vector
        Output: vector

        Replaces values < 0 with 0 in input-vector
    """
    new_arr = np.zeros(len(arr))

    for i in range(len(arr)):
        
        if arr[i] > 0:
            new_arr[i] = arr[i]
        
        else:
            new_arr[i] = 0

    return new_arr