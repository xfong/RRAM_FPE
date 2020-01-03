from fipy import *
nx = 100e3;
dx = 1.0e-3;
L = nx * dx
midPoint = L / 2
mesh = Grid1D(nx=nx, dx=dx)
phi = CellVariable(name="solution variable",
                   mesh=mesh,
                   value=1/L)

timeScaleFac = 1e6
D = 1.0 * timeScaleFac

cScaleValue = 1.0e0 * timeScaleFac
cCoeff = FaceVariable(mesh=mesh, value=[1.0])
Xgrid = mesh.faceCenters[0]
cCoeff.setValue(cScaleValue * (Xgrid - midPoint))
#cCoeff.setValue(0., where=(Xgrid < dx) | (Xgrid > L - dx))
#phi.faceGrad.constrain([0.], mesh.facesLeft)
#phi.faceGrad.constrain([0.], mesh.facesRight)

eqX = TransientTerm() == (ExplicitDiffusionTerm(coeff=D) + PowerLawConvectionTerm(coeff=cCoeff))
eqI = TransientTerm() == (DiffusionTerm(coeff=D) + PowerLawConvectionTerm(coeff=cCoeff))
#timeStepDuration = 0.09 * dx**2 / (2 * D)
timeStepDuration = 1.0e-18
steps = 5
viewer1 = Viewer(vars=(phi),
                datamin=0., datamax=1)
viewer2 = Viewer(vars=(phi),
                datamin=0., datamax=10)
if __name__ == '__main__':
    viewer1.plot()

for step in range(steps-1):
    print("Loop 1, Step Count")
    print(step)
    eqX.solve(var=phi, dt=timeStepDuration)
    if __name__ == '__main__':
        viewer1.plot()

timeStepDuration = 1.0e-6
steps = 200
for step in range(steps-1):
    print("Loop2 Step Count")
    print(step)
    eqI.solve(var=phi, dt=timeStepDuration)
    if __name__ == '__main__':
        viewer2.plot()

if __name__ == '__main__':
    input("Steady-state drift-diffusion. Press <return> to proceed...")
