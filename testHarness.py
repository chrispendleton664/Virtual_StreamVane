import swirlGenerator as sg

# Uniform streamwise velocity
axialVel = 1

# Mesh variables
xSide = 10                 # Side length of grid in x direction
ySide = 10                 # Side length of grid in y direction
xNumCells = 1000           # Num of cells in x direction
yNumCells = 1000           # Num of cells in y direction


# Initialise flow fields
field1 = sg.FlowField([xSide,ySide],[xNumCells,yNumCells])
field2 = field1.copy()
field3 = field1.copy()
field4 = field1.copy()

# Initialise domain configurations
bulkSwirl       = sg.Vortices(3, [[0,0]], [18], [5], [-1], [axialVel])
twinSwirlLO     = sg.Vortices(2, [[-1.25,0],[1.25,0]], [5,-5])
twinSwirlIso    = sg.Vortices(1, [[-1.25,0],[1.25,0]], [5,-5])
randomSwirl     = sg.Vortices(1, [[2,0],[3,3],[-3,-1]], [-5,2,3])

# Stack objects into array for iteration
flowFields = [field1,field2,field3,field4]
swirlDefs = [bulkSwirl, twinSwirlLO, twinSwirlIso, randomSwirl]
swirlFiles = ["iBulkSwirl.pdf", "iLO_TwinSwirl.pdf", "iIso_TwinSwirl.pdf", "iRandomSwirl.pdf"]

# Loop through swirl description
for i in range(len(flowFields)):
    # Calculate domain
    print(f"Calculating {swirlFiles[i][0:-4]}...")
    velGrids = flowFields[i].defineVortices(swirlDefs[i], axialVel=axialVel)

    # Save fields to pdf
    flowFields[i].plotAll(f"flow_fields/{swirlFiles[i]}")
