import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode

N = 20 # Number of neurons.
mu = 1
R = 5
sigma = 0.7
urest = 0
uth = 0.98
circles = np.zeros(N)
u = np.zeros(N)

def fun(t, z):
    """
    Right hand side of the differential equations
      dx/dt = -omega * y
      dy/dt = omega * x
    """
    for i in range(N):
        u[i] = z[i]
    
    f = np.zeros(N)
    for i in range(N):
        f[i] = mu - u[i]
        
        coupling = 0 # Assuming R < N/2.
        if i-R < 0: 
            for j in range(i+R+1):
                coupling = coupling - sigma*(u[j] - u[i])/(2*R)

            f[i] = f[i] + coupling
                
            for j in range((i-R)%N, N):
                coupling = coupling - sigma*(u[j] - u[i])/(2*R)

            f[i] = f[i] + coupling       
                
        elif i+R > N-1:
            for j in range(i-R, N):
                coupling = coupling - sigma*(u[j] - u[i])/(2*R)

            f[i] = f[i] + coupling
                
            for j in range((i+R)%N+1):
                coupling = coupling - sigma*(u[j] - u[i])/(2*R)

            f[i] = f[i] + coupling   
                
        else:
            for j in range(i-R, i+R+1):
                coupling = coupling - sigma*(u[j] - u[i])/(2*R)

            f[i] = f[i] + coupling
                
        if z[i] > uth:
            z[i] = 0
            circles[i] = circles[i] + 1
     
    return f

# Create an `ode` instance to solve the system of differential
# equations defined by `fun`, and set the solver method to 'dop853'.
solver = ode(fun)
solver.set_integrator('dopri5')

# Set the initial value z(0) = z0.
t0 = 0.0
z0 = np.random.rand(N)*uth
solver.set_initial_value(z0, t0)

# Create the array `t` of time values at which to compute
# the solution, and create an array to hold the solution.
# Put the initial value in the solution array.
t1 = 3000
dt = 0.01
steps = round((t1-t0)/dt)
t = np.linspace(t0, t1, steps)
sol = np.empty((steps, N))
sol[0] = z0

# Repeatedly call the `integrate` method to advance the
# solution to time t[k], and save the solution in sol[k].
k = 1
while solver.successful() and solver.t < t1:
    solver.integrate(t[k])
    sol[k] = solver.y
    if k%1000 == 0:
        print('step =', k,'/', steps)
        plt.scatter(range(N), sol[k])
        plt.ylim([0, 1])
        plt.show()
        
    k += 1

