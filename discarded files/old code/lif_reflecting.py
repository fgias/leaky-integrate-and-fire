import numpy as np
from matplotlib import pyplot as plt 

from numpy import asarray
from numpy import savetxt
#from numpy import loadtxt

N = 10 # Number of neurons.
mu = 1
R = 4
sigma = 0.7
urest = 0
uth = 0.98
dt = 0.01

u = np.random.rand(N)*uth
circles = np.zeros(N)
unext = np.zeros(N)


def check_uth(i):
    if unext[i] > uth:
        unext[i] = 0
        circles[i] = circles[i] + 1


def update_u():
    for i in range(N):
        unext[i] = u[i] + dt*(mu - u[i])
        
        coupling = 0 # Assuming R < N/2.
        
        if N-i-R < 0: 
            for j in range(N-i+R+1):
                coupling = coupling - sigma*(u[j] - u[i])/(2*R)
                
            unext[i] = unext[i] + dt*coupling
                
            for j in range((N-i-R)%N, N):
                coupling = coupling - sigma*(u[j] - u[i])/(2*R)
                
            unext[i] = unext[i] + dt*coupling   
                
        elif N-i+R > N-1:
            for j in range(N-i-R, N):
                coupling = coupling - sigma*(u[j] - u[i])/(2*R)
                
            unext[i] = unext[i] + dt*coupling
                
            for j in range((N-i+R)%N+1):
                coupling = coupling - sigma*(u[j] - u[i])/(2*R)
                
            unext[i] = unext[i] + dt*coupling   
                
        else:
            for j in range(N-i-R, N-i+R+1):
                coupling = coupling - sigma*(u[j] - u[i])/(2*R)
                
            unext[i] = unext[i] + dt*coupling
        
        check_uth(i)


def run(time, timemax, plotfreq = 1000, plotting = True):
    plotnum = 0
    print('--------------------------')
    print('Starting calculation')
    print('--------------------------')
    while time < timemax:
        update_u()
        
        for i in range(N):
            u[i] = unext[i]

            
        time = time + dt
        plotnum = plotnum + 1 # Live plotting.
        print('Time =', time)
            
        if plotnum%plotfreq == 0 and plotting == True: # Live plotting.
            plt.scatter(range(N), u)
            plt.ylim([0, 1])
            plt.show()
            
        if plotnum%1000 == 0: # Save data.
            udata = asarray(u)
            
            timearray = np.array([time])
            timedata = asarray(timearray)
            circdata = asarray(circles)
            savetxt('data/udata.csv', udata, delimiter=',')
            savetxt('data/timedata.csv', timedata, delimiter=',')
            savetxt('data/circledata.csv', circdata, delimiter=',')
    
    print('--------------------------')
    print('End of calculation')
    print('--------------------------')


run(0, 100, 100)

""" 
To restore data:
    
u = loadtxt('data/udata.csv', delimiter=',')
timearray = loadtxt('data/timedata.csv', delimiter=',')

"""

        


    

