import scipy.integrate as integrator
import numpy as np
import matplotlib.pyplot as plt

x_left = -25e-9
x_right = 25e-9
numPts = 10000

xGrid = np.linspace(x_left, x_right, num=numPts)
rho_init = np.ones(numPts) / (x_right - x_left)

vCoeff = 1e3
dCoeff = 1e3
tScale = 1e0

plt.figure(1)
plt.clf()
fig, ax = plt.subplots()
ax.plot(xGrid, rho_init)
##ax.legend(loc='best')
input('Press <Enter> to continue...')

'''
def velocity_func(t):
    velocity = -1.0 * vCoeff * xGrid
    return velocity

def flow_func(t, y):
    drift = velocity_func(t)*y
    diffusion = dCoeff * gradient(y)
    flow = drift - diffusion
    return flow

def fpe_rhs(t, y):
    drho_dt = -tScale*gradient(flow_func(t, y))
    return drho_dt

end_time = 1e-15
sol = integrator.solve_ivp(fpe_rhs, [0, end_time], rho_init, method='BDF', t_eval=end_time, vectorized=true, first_step=1e-21, max_step=1e-18, atol=1e-9, rtol=1e-3)

plt.figure(1)
plt.clf()
fig, ax = plt.subplots(num=1)
ax.plot(xGrid, sol.y[0], 'k--', label='Result')
ax.legend(loc='best')
input('Press <Enter> to continue...')
'''
