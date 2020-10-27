import swirlGenerator as sg
from createSU2 import createInlet

# Uniform streamwise velocity
axialVel = 1

# Mesh variables
xSide = 10                 # Side length of grid in x direction
ySide = 10                 # Side length of grid in y direction
xNumCells = 100            # Num of cells in x direction
yNumCells = 100            # Num of cells in y direction


# Initialise flow field
flowField = sg.FlowField([xSide,ySide],[xNumCells,yNumCells])

# Initialise domain configuration
oneVortex    = sg.Vortices(1, [[0,0]], [5])

# Calculate velocity field
flowField.defineVortices(oneVortex, axialVel=axialVel)

# Create SU2 boundary condition
createInlet(flowField, 'test.dat')

# Plot fields
flowField.plotAll()
