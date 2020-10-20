import swirlGenerator as sg

# Uniform streamwise velocity
axialVel = 1

# Vortex Types
vortTypes = [1,2,3]
# Vortex centers
vortO = [[-2,0],[2,0],[0,0]]
# Vortex strengths
vortS = [-5,5,16]

# Mesh variables
xSide = 10                 # Side length of grid in x direction
ySide = 10                 # Side length of grid in y direction
xNumCells = 1000           # Num of cells in x direction
yNumCells = 1000           # Num of cells in y direction

# Pdf to output plots
outputPdf = '___.pdf'


# Setup coordinates grid
coordGrids = sg.makeGrid([xSide,ySide],[xNumCells,yNumCells])

# Initialise object to store data about multiple vortices of different types
VortexDefs = sg.Vortices(vortTypes, vortO, vortS, [5], [-1], [axialVel])

# Place vortices in domain
velGrids, rho, p = sg.defineVortices(VortexDefs, coordGrids, axialVel)

# Plot and save fields
sg.plotAll(coordGrids, velGrids, rho, p)