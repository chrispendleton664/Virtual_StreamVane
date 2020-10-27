from ... import swirlGenerator as sg
from ...createSU2 import createInlet

# Uniform streamwise velocity
axialVel = 1

# Mesh variables - should match the mesh being used in CFD solver
xSide = 1
ySide = 1
xNumCells = 10
yNumCells = 10

# Intialise flow field object
flowField = sg.FlowField([xSide,ySide],[xNumCells,yNumCells])

# Intialise domain configuration object
vortexDefs = sg.Vortices(1, [[0,0]], [2])

# Calculate velocity field
flowField.defineVortices(vortexDefs, axialVel=axialVel)

# Write SU2 boundary condition
createInlet(flowField, 'small_centered_isentropic.dat')

# Show flow fields
flowField.plotAll()