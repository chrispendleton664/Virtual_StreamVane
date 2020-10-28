from matplotlib.pyplot import plot
import swirlgenerator.core as sg
import swirlgenerator.writeBC as bc
import swirlgenerator.maketestdomain as domain
import sys

'''
Check command line arguements
'''
# Configuration file name
if (len(sys.argv) < 2) or (sys.argv[1].find('-') == 0):
    raise RuntimeError('Configuration file missing')
else:
    configFile = sys.argv[1]

# Make meshed test domain
if '-makemesh' in sys.argv:
    # Get index of this command line arguement
    idx = sys.argv.index('-makemesh')

    # Check if there is a pair
    try:
        if (sys.argv[idx+1].find('-') == 0):
            raise RuntimeError
        else:
            # Use as filename
            meshFile = sys.argv[idx+1]
            makemesh = True
    except:
        raise RuntimeError("-makemesh arguement defined but no mesh filename given")
else:
    makemesh = False
    meshFile = None

# Show created test mesh
showmesh = (True if '-showmesh' in sys.argv else False)

# Show plots
showFields = (True if '-show' in sys.argv else False)

# Save plots
if '-saveplots' in sys.argv:
    # Get index of this command line arguement
    idx = sys.argv.index('-saveplots')

    # Check if there is a pair
    try:
        if (sys.argv[idx+1].find('-') == 0):
            raise RuntimeError
        else:
            # Use as filename
            plotsFile = sys.argv[idx+1]
            saveplots = True
    except:
        raise RuntimeError("-saveplots arguement defined but no pdf filename given")
else:
    saveplots = False
    plotsFile = None


'''
Call all necessary functions to create boundary condition from input data
'''
# Initialise Input object and read config file
inputData = sg.Input()
inputData.read(configFile)

# Intialise flow field object with coordinate system
flowField = sg.FlowField([inputData.xSide,inputData.ySide], [inputData.xNumCells, inputData.yNumCells])

# Initialise domain configuration object with vortex definitions
vortexDefs = sg.Vortices(inputData.vortModel, inputData.vortCoords, inputData.vortStrengths)

# Calculate velocity field
if inputData.axialVel != None:
    # Define the uniform axial velocity
    flowField.defineVortices(vortexDefs, axialVel=inputData.axialVel)
else:
    # Leave as default
    flowField.defineVortices(vortexDefs)

# Write inlet boundary condition file
bc.writeInlet(InputObject=inputData, flowField=flowField)


'''
Optional functionalities
'''
# Create a test mesh compatible with the boundary condition that's been generated
if makemesh:
    domain.simpleBox(inputData, meshFile, showmesh)

# Save flow fields in pdf if requested
if saveplots:
    flowField.plotAll(plotsFile)

# Show flow fields if requested
if showFields:
    flowField.plotAll()
