"""
LIF nonlocal.
F.I. Giasemis.

"""

import numpy as np
import time
from matplotlib import pyplot as plt 
from numpy import asarray
from numpy import savetxt
from numpy import loadtxt

N = 20 # Number of neurons.
mu = 1
R = 6
sigma = 0.7
urest = 0
uth = 0.98
dt = 0.001

u = np.random.rand(N)*uth
circles = np.zeros(N)

start_time = time.time()


def update_u(u_in):
    u_next = u_in + dt*(mu - u_in)

    k = np.random.randint(N)

    for i in range(k,k+1):
        coupling = 0 # R <= N/2.
        
        if i-R < 0: 
            for j in range(i+R+1):
                coupling = coupling - sigma*(u_in[j] - u_in[i])/(2*R)
                
            u_next[i] = u_next[i] + dt*coupling
                
            for j in range((i-R)%N, N):
                coupling = coupling - sigma*(u_in[j] - u_in[i])/(2*R)
                
            u_next[i] = u_next[i] + dt*coupling   
                
        elif i+R > N-1:
            for j in range(i-R, N):
                coupling = coupling - sigma*(u_in[j] - u_in[i])/(2*R)
                
            u_next[i] = u_next[i] + dt*coupling
                
            for j in range((i+R)%N+1):
                coupling = coupling - sigma*(u_in[j] - u_in[i])/(2*R)
                
            u_next[i] = u_next[i] + dt*coupling   
                
        else:
            for j in range(i-R, i+R+1):
                coupling = coupling - sigma*(u_in[j] - u_in[i])/(2*R)
                
            u_next[i] = u_next[i] + dt*coupling
        
    return u_next


def run(t, timemax, u_in, plotfreq = 1000, plotting = True):
    plotnum = 0
    print('--------------------------')
    print('Starting calculation')
    print('--------------------------')
    while t < timemax:
        u_in = update_u(u_in)
            
        t = t + dt
        plotnum = plotnum + 1 # Live plotting.
        # print('t =', t)
            
        if plotnum%plotfreq == 0 and plotting == True: # Live plotting.
            print('t =', round(t))
            plt.scatter(range(N), u_in)
            plt.ylim([0, 1])
            plt.ylabel("ui")
            plt.show()
            
        if plotnum%1000 == 0: # Save data.
            udata = asarray(u)
            
            timearray = np.array([t])
            timedata = asarray(timearray)
            circdata = asarray(circles)
            savetxt('data/udata.csv', udata, delimiter=',')
            savetxt('data/timedata.csv', timedata, delimiter=',')
            savetxt('data/circledata.csv', circdata, delimiter=',')
    
    print('--------------------------')
    print('End of calculation')
    print('--------------------------')
    
    return u_in


u = run(0, 500, u, 5000)

print("Time: %s seconds." % (round(time.time() - start_time, 2)))
print('--------------------------')

""" 
To restore data:
    
u = loadtxt('data/udata.csv', delimiter=',')
timearray = loadtxt('data/timedata.csv', delimiter=',')

"""

        


    

