from matplotlib.pyplot import plot
import swirlgenerator.core as sg
import swirlgenerator.writeBC as bc
import sys

'''
Check command line arguements
'''
# Configuration file name
if (len(sys.argv) < 2) or (sys.argv[1].find('-') == 0):
    raise RuntimeError('Configuration file missing')
else:
    configFile = sys.argv[1]

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
inputData.read('example.config')

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
# Save flow fields in pdf if requested
if saveplots:
    flowField.plotAll(plotsFile)
    
# Show flow fields if requested
if showFields:
    flowField.plotAll()
