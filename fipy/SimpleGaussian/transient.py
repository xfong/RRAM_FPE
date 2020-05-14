from fipy import *
import scipy.integrate as integrate
import fipy.tools.dump as dumpster

## Set up simulation grid
nx = 1e5;
dx = 1.0e-3;
L = nx * dx
midPoint = L / 2
mesh = Grid1D(nx=nx, dx=dx)
dumpster.write(mesh, filename="grid_info.dat")

## The variable of interest
phi = CellVariable(name="\Phi",
                   mesh=mesh,
                   value=1/L)
dumpster.write(phi,filename="initial_phi.dat")

## Definition for scaling in time
timeScaleFac = 1e6

## Definition for diffusion
D = 1.0 * timeScaleFac

## Definition for drift term
cScaleValue = 1.0e0 * timeScaleFac
cCoeff = FaceVariable(mesh=mesh, value=[1.0])
Xgrid = mesh.faceCenters[0]
## Drift term is linear in x-axis (gap length)
cCoeff.setValue(cScaleValue * (Xgrid - midPoint))

## Set boundary conditions (if required)
#cCoeff.setValue(0., where=(Xgrid < dx) | (Xgrid > L - dx))
#phi.faceGrad.constrain([0.], mesh.facesLeft)
#phi.faceGrad.constrain([0.], mesh.facesRight)

# Define an explicit form of the Fokker-Planck equation
# and an implicit form. Due to the initial conditions,
# the implicit form will have numerical issues to start
# solving. The explicit form does not have this issue
# but it will have stability issues, giving rise to
# spurious oscillations in the result, which are
# numerical artifacts. Thus, we will solve the explicit
# form for a few initial time steps to get to a point
# that the implicit solver likes, and then switch over
# to the implicit solver to give the rest of the
# solution

## Explicit form of Fokker-Planck Equation
eqX = TransientTerm() == (ExplicitDiffusionTerm(coeff=D) + PowerLawConvectionTerm(coeff=cCoeff))
## Implicit form of Fokker-Planck Equation
eqI = TransientTerm() == (DiffusionTerm(coeff=D) + PowerLawConvectionTerm(coeff=cCoeff))

## Use extremely small time steps for the explicit
## solver to minimize spurious oscillations in the
## initial few time steps
timeStepDuration = 1.0e-18
steps = 5
viewer1 = Matplotlib1DViewer(vars=(phi),
                datamin=0., datamax=1.5)
viewer2 = Matplotlib1DViewer(vars=(phi),
                datamin=0., datamax=1.5)
t_i = 0.0

if __name__ == '__main__':
    viewer1.plot()
    input("Initial phi. Press <return> to proceed...")

## Loop to step through time
for step in range(steps-1):
    print("Loop 1, Step Count")
    print(step)
    eqX.solve(var=phi, dt=timeStepDuration)
    t_i = t_i + timeStepDuration

## Switch over to implicit solver using
## larger time steps
timeStepDuration = 1.0e-6
steps = 20

## Loop to step through time
for step in range(steps-1):
    dest_file = 'phi_transient_' + str(t_i) + 'sec.vtk'
    if __name__ == '__main__':
        viewer2.plot()
    print("Loop2 Step Count")
    print(step)
    eqI.solve(var=phi, dt=timeStepDuration)
    t_i = t_i + timeStepDuration

dumpster.write(phi,filename="final_phi.dat")

if __name__ == '__main__':
    viewer2.plot()
    input("Final phi. Press <return> to proceed...")

