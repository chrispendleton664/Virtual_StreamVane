# --------------------------------------------------------------------------------
#
# Script to generate a test boundary condition
# -small Lamb-Oseen vortex in the middle of the domain
#
# --------------------------------------------------------------------------------

import swirlgenerator.core as sg
from swirlgenerator.writeBC import writeSU2
import sys

# Uniform streamwise velocity
axialVel = 1

# Mesh variables - should match the mesh being used in CFD solver
xSide = 1
ySide = 1
xNumCells = 20
yNumCells = 20

# Intialise flow field object
flowField = sg.FlowField([xSide,ySide],[xNumCells,yNumCells])

# Intialise domain configuration object
vortexDefs = sg.Vortices(2, [[0,0]], [2])

# Calculate velocity field
flowField.defineVortices(vortexDefs, axialVel=axialVel)

# Write SU2 boundary condition
writeSU2(flowField, 'small_centered_LO.dat')

# Show flow fields
if '-showfields' in sys.argv:
    flowField.plotAll()