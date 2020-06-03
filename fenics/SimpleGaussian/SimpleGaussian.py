'''
The Fokker-Planck equation can be used to model the time evolution of a probability
density in a potential well. At steady-state, the diffusion process is exactly
balanced by the forces due to the potential well. The integration of the
probability density function yields the cumulative density function, which can be
compared to experimental data.

In this script, we model the effect of a 1D quadratic potential well.
Theoretically, the distribution will be approximately Gaussian if the potential
well is sufficiently deep as compared to the diffusion coefficient
'''
from fenics import *
from mshr import *
import numpy as np
import matplotlib.pyplot as plt

'''
Generate mesh in the interval from 0 to 1
This can be rescaled to give us the desired results
'''
mesh=UnitIntervalMesh(10001)

'''
The potential well and point velocities on the mesh are related as the relationship
between electric-field and the electrostatic potential:
velocity = -grad(potential)
For a quadratic potential well, the velocity will have an expression
velocity = -1.0 * beta * (x - x0)
where beta is a positive scalar factor, and x0 is the point at which the potential
well is at the minimum.
'''
beta = Constant(1000.0)
x0 = Constant(0.5)
velocity = Expression(('coeff*(offset-x[0])'), coeff=beta, offset=x0, domain=mesh, degree=1)

'''
The diffusion constant, D, is related to the Brownian motion that causes the
probability density to spread out evenly across the state-space. If rho is the
probability density function, then the diffusion process is described as:
d(rho)/dt = div(-1.0 * D * grad(rho)) = -1.0 * D * Laplacian(rho)
'''
D = Constant(2.0*beta)

'''
Initialize parameters for time-stepping
'''
num_steps=301
t=0
t_next=0
errTol=1e-4
dt=1e-18
delT=Constant(dt)
'''
Set up the variational form of the Fokker-Planck equation. This is formally written as:
v*d_rho*dx + dt*(div(u*rho)*dx) = dt*D*dot(grad(rho),grad(v))*dx
where u is the velocity and rho is the probability density function. If steady-state
results are desired, then the time-dependence of rho is zero, and the variational form
is instead written as:
div(velocity*rho)*dx = D*dot(grad(rho),grad(v))*dx
'''

#### Define basis function and test function
V=FunctionSpace(mesh,'CG',1)
drho_trial=TrialFunction(V)      # For forward Euler step
rho_pred=Function(V)        # For forward Euler step
rho_next=TrialFunction(V)        # For backward Euler corrector
rho_curr=Function(V)
v=TestFunction(V)

'''
Define solution at time = 0 (i.e., initial condition)
We will use a quadratic function with peak at the middle of the interval
rho0 = -omega*x*(x-1) = -omega*x**2 + omega*x
To find omega, we need to ensure the integration over the interval
gives result equal to one
Difference between -(omega/3)*x**3 + (omega/2)*x**2 evaluated at x = 1
and at x = 0 must be equal to 1
-omega/3 + omega/2 = omega*(1/2 - 1/3) = omega/6 = 1
i.e. omega = 6
'''
rho0=Expression('-6.0*x[0]*(x[0]-1.0)', degree=1)
rho_curr.interpolate(rho0)

#### Give the variational form of the Fokker-Planck equation
a1= drho_trial*v*dx                                                                               # For forward Euler
L1= -rho_curr*dot(velocity,v.dx(0))*dx - dot(D*grad(rho_curr),grad(v))*dx                         # For forward Euler
a2= rho_next*v*dx                                                                                 # For backward Euler
L2= rho_curr*v*dx - dt*rho_pred*dot(velocity,v.dx(0))*dx - dt*dot(D*grad(rho_pred),grad(v))*dx    # For backward Euler

#### Define temporary variables for time-stepper
drho_dt=Function(V)
rho_next=Function(V)

#### Define the problem and solver inside time-stepper
prob1=LinearVariationalProblem(a1, L1, drho_dt)
solv1=LinearVariationalSolver(prob1)
solv1.parameters["linear_solver"] = "gmres"
solv1.parameters["preconditioner"] = "ilu"
prob2=LinearVariationalProblem(a2, L2, rho_next)
solv2=LinearVariationalSolver(prob2)
solv2.parameters["linear_solver"] = "gmres"
solv2.parameters["preconditioner"] = "ilu"

#### Define the boundary conditions
#### Neumann boundary conditions are automatically defined
bc=[]

# Generate internal progess output to screen
progress = Progress('Time-stepping', num_steps)
set_log_level(LogLevel.PROGRESS)

# Create VTK file for saving solution
vtkfile = File('FP_Fenics_result_files/solution.pvd')
#vtkfile << (rho_curr, t)

for n in range(num_steps):
        # Perform forward Euler step
        solv1.solve()

        print('Calculating next time from t = ' + str(t) + ' sec')

        err_total = 10.0*errTol
        while errTol < err_total:
                print('    using dt = ' + str(dt))
                # Update temporary time step (this will be iterated until error meets tolerance)
                t_next = t + dt
                rho_pred.assign(rho_curr + dt*drho_dt)

                # Perform backward Euler correction
                solv2.solve()

                # Determine error
                err_total = abs(np.array(rho_next.vector()) - np.array(rho_pred.vector())).max()

                # Adapt time step
                scaleFac = 0.8*errTol/err_total
                if (scaleFac > 9.0):
                        scaleFac = 9.0 # Restrict expansion of dt
                dt = scaleFac*dt

        # Update current data point once error meets tolerance
        print('    ...accepted dt = ' + str(dt))
        rho_curr.assign(rho_next)
        t = t_next
        if (n%10 == 0):
                print('Saving VTK file...')
                vtkfile << (rho_curr, t)
                print('VTK file saved!')

plot(rho_curr)
plot(mesh)
vtkfile << (rho_curr, t)
plt.show()
