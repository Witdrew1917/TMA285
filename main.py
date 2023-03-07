from mpl_toolkits import mplot3d
import sys
import numpy as np
import matplotlib.pyplot as plt
from Projekt1 import DifferenceScheeme


def plot(sigma, call, put, price, num, subcap):
    
    plt.figure(num)
    plt.plot(sigma, call)
    plt.plot(sigma, put)
    plt.xlabel("sigma")
    plt.ylabel("call/put")
    plt.legend(['call', 'put'])
    plt.title("Initial price of the call/put option given S0 = {}".format(price))
    plt.show()
    return 0

if __name__ == '__main__':
    
    batch = ""
    if len(sys.argv) > 1:
        batch = sys.argv[1]

    f = "function_string"
    time_samples = 20   
    sigma_samples = 10
    price_samples = 5
    T = 1/2 
    dt = T/time_samples
    S0 = 80
    r = 0.05
    sigma = 0.3   
    K = 100
    if batch == "batch":

        num_sim = 10     

        PRICE = np.linspace(S0,S0 + 10*price_samples, price_samples)
        num = 0
        for price in PRICE: 
            num += 1

            init_price_call = []
            init_price_put = []
            SIGMA = np.linspace(0.01,sigma,sigma_samples)
            for sigma in SIGMA:
                scheeme = DifferenceScheeme.Scheeme(time_samples,time_samples,dt)
                scheeme.solve_PDE(price,r,sigma,K,T)
                init_price_call.append(scheeme._call_put[0][0])
                init_price_put.append(scheeme._call_put[1][0])
                
            plot(SIGMA, init_price_call, init_price_put, price, num, price_samples)
        
    else:
        scheeme = DifferenceScheeme.Scheeme(time_samples,time_samples,dt)
        scheeme.solve_PDE(S0,r,sigma,K,T)
        scheeme.plot(T,K)
        
            


