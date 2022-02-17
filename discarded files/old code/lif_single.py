import numpy as np
from matplotlib import pyplot as plt 

N = 1
dt = 0.01
mu = 1
urest = 0
uth = 0.98

u = np.random.rand(N)*0
unext = np.zeros(N)
T=np.log(50)+1

def update_u():
    for i in range(N):
        unext[i] = u[i] + dt*(mu - u[i])
        
        if unext[i] > uth:
            unext[i] = 0

def run(time, timemax):
    step = 0
    sol = np.zeros(round(timemax/dt))
    while time < timemax-1:
        sol[step] = u[0]
        update_u()
        for i in range(N):
            u[i] = unext[i]
            
        time = time + dt
        step = step + 1 # Live plotting.
        print('Time =', time)
            
    return sol


solution = run(0, T)

plt.plot(np.arange(0, T-0.01,0.01), solution)
plt.ylim([0, 1])
plt.yticks([0, 0.5, 1])
plt.ylabel("u")
plt.xticks([])
plt.xlabel("time")
#plt.savefig('lif_single.png', dpi=800)

def update_u():
    for i in range(N):
        unext[i] = u[i] + dt*(mu - 0.8*u[i])
        
        if unext[i] > uth:
            unext[i] = 0
            
solution = run(0, T)

plt.plot(np.arange(0, T-0.01,0.01), solution)
plt.ylim([0, 1])
plt.yticks([0, 0.5, 1])
plt.ylabel("u")
plt.xticks([])
plt.xlabel("time")
#plt.savefig('lif_single.png', dpi=800)