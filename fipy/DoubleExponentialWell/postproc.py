from fipy import *
import fipy.tools.dump as dumpster

def getPhiVal(xcoor, phi):
  return phi([[xcoor]], order=1)

origGrid = dumpster.read("grid_info.dat")
origPhi = dumpster.read("final_phi.dat")

for idx in range(100):
  xcoor = idx * (1.0 / 1.0)
  print(getPhiVal(xcoor, origPhi))

