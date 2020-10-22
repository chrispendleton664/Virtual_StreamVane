import swirlGenerator as sg

# Uniform streamwise velocity
axialVel = 1

# Mesh variables
xSide = 10                 # Side length of grid in x direction
ySide = 10                 # Side length of grid in y direction
xNumCells = 1000           # Num of cells in x direction
yNumCells = 1000           # Num of cells in y direction


# Setup coordinates grid
coordGrids = sg.makeGrid([xSide,ySide],[xNumCells,yNumCells])

# Initialise domain configurations
bulkSwirl       = sg.Vortices(3, [[0,0]], [18], [5], [-1], [axialVel])
twinSwirlLO     = sg.Vortices(2, [[-1.25,0],[1.25,0]], [5,-5])
twinSwirlIso    = sg.Vortices(1, [[-1.25,0],[1.25,0]], [5,-5])
randomSwirl     = sg.Vortices(1, [[2,0],[3,3],[-3,-1]], [-5,2,3])

# Stack objects into array for iteration
swirlDefs = [bulkSwirl, twinSwirlLO, twinSwirlIso, randomSwirl]
swirlFiles = ["__BulkSwirl.pdf", "__LambOseen_TwinSwirl.pdf", "__Isentropic_TwinSwirl.pdf", "__Random_Mixed_Swirl.pdf"]

# Loop through swirl description
for i in range(len(swirlDefs)):
    # Calculate domain
    print(f"Calculating {swirlFiles[i][0:-4]}...")
    velGrids, rho, p = sg.defineVortices(swirlDefs[i], coordGrids, axialVel)

    # Save fields to pdf
    sg.plotAll(coordGrids, velGrids, rho, p, f"flow_fields/{swirlFiles[i]}")
