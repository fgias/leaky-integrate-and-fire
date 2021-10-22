import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
import time

mu = 1
u_th = 0.98
N = 60
R = round(0.35*N)
sigma = 1.9
time_max = 3000
dt = 0.1
circles = np.zeros(N)
start_time = time.time()

plotfreq = 10

def fun(t, z):
    """
    Right hand side of the differential equations
      dx/dt = -omega * y
      dy/dt = omega * x
    """
    u = z
    f = mu - u

    for i in range(N):
        coupling = 0 # R <= N/2.
        if i-R < 0:
            for j in range((i-R)%N,N):
                coupling = coupling + sigma*(u[i]-u[j])/(2*R)

            for j in range(i+R+1):
                coupling = coupling + sigma*(u[i]-u[j])/(2*R)

            f[i] = f[i] + coupling
        elif i+R > N-1:
            for j in range(i-R,N):
                coupling = coupling + sigma*(u[i]-u[j])/(2*R)

            for j in range((i+R)%N+1):
                coupling = coupling + sigma*(u[i]-u[j])/(2*R)

            f[i] = f[i] + coupling
        else: 
            for j in range(i-R,i+R+1):
                coupling = coupling + sigma*(u[i]-u[j])/(2*R)

            f[i] = f[i] + coupling

        if u[i] > u_th:
            z[i] = 0
            circles[i] += 1

    return f

# Create an `ode` instance to solve the system of differential
# equations defined by `fun`, and set the solver method to 'dop853'.
solver = ode(fun)
solver.set_integrator('dopri5')

# Give the value of omega to the solver. This is passed to
# `fun` when the solver calls it.


# Set the initial value z(0) = z0.
seed_num = 567632
np.random.seed(seed_num)
z0 = np.random.rand(N)*u_th

t0 = 0.0
solver.set_initial_value(z0, t0)

# Create the array `t` of time values at which to compute
# the solution, and create an array to hold the solution.
# Put the initial value in the solution array.
t1 = time_max
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
    if k%plotfreq == 0:
        print('step =', k,'/', steps)
        plt.scatter(range(N), sol[k])
        plt.ylim([0, 1])
        plt.show()

    k += 1
    
time_calc = time.time() - start_time
print(time_calc)

w = 2*np.pi*circles/time_max
plt.scatter(range(1, N+1), w)
plt.ylabel("ωi")
plt.yticks([round(min(w),2), round((min(w)+max(w))/2,2), round(max(w),2)])
plt.xticks([1, N/2, N])
plt.xlabel("index i")
plt.show()