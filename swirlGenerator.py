import numpy as np
from matplotlib import pyplot as plt
import warnings
from typing import Union

'''
Heavy use of numpy for fast matrix and vector operations
Matplotlib for doing visualisations
'''

# Plots output pdf file
outputPdf = '___.pdf'

# Fluid parameters
GAMMA = 1.4
kin_visc = 1.81e-5

# Uniform streamwise velocity
axialVel = 1

# Vortex Types
vortType = [1,2,3]
# Vortex centers
vortO = [[-2,0],[2,0],[0,0]]
# Vortex strengths
vortS = [-5,5,16]

# Mesh variables
xSide = 10                 # Side length of grid in x direction
ySide = 10                 # Side length of grid in y direction
xNumCells = 1000           # Num of cells in x direction
yNumCells = 1000           # Num of cells in y direction

'''
For storing information about the vortices which have been defined for the domain
Each vortex can be defined as a different type
- Vortex types: 1 - Isentropic; 2 - Lamb-Oseen; 3 - Forced/Solid
- Positive vortex strength is defined as anti-clockwise rotation
'''
class Vortices: 
    # Object constructor accepts lists also for convenience, then later converts to numpy arrays
    def __init__(self, types: Union[list, np.ndarray], centres: Union[list, np.ndarray], strengths: Union[list, np.ndarray], radius: Union[list, np.ndarray] = [], clockwise: Union[list, np.ndarray] = []):
        self.numVortices    = len(types)

        # Make sure these are all numpy arrays not just lists
        self.centres        = (centres      if isinstance(centres,np.ndarray)   else np.array(centres))       # Vortex centre
        self.strengths      = (strengths    if isinstance(strengths,np.ndarray) else np.array(strengths))     # Vortex strength
        self.types          = (types        if isinstance(types,np.ndarray)     else np.array(types))         # Vortex type - which mathematical model to use
        self.radius         = (radius       if isinstance(radius,np.ndarray)    else np.array(radius))        # Vortex radius - can define a strong edge to the vortex, it has no effect on the flow outside this
        self.clockwise      = (clockwise    if isinstance(clockwise,np.ndarray) else np.array(clockwise))     # Rotation direction - only needed for some types of vortices

        self.vortNum        = 0             # For keeping track during iteration

        # If some but not all defined are solid vortices, and radius has not been defined for all of them;
        # create sparse matrices so that radius and clockwise values line up with correct vortex.
        # So that you don't have to manually define dummy radii/directions for other vortex types
        if (any(self.types == 3) and any(self.types != 3) and len(self.types) != len(self.radius)):
            self.radius = np.zeros(self.strengths.shape)
            self.clockwise = self.radius.copy()

            forced = self.types == 3
            
            self.radius[forced]     = radius
            self.clockwise[forced]  = clockwise

    # Return data for next vortex as tuple, for iteration
    def getNextVortex(self):
        # Check if at end of vortex list
        if (self.vortNum == self.numVortices):
            raise IndexError(f"Index {self.vortNum} is out of bounds of vortex list with size {self.numVortices}")
        else:
            # Output correct tuple format depending on vortex types
            if (any(self.types == 3)):
                data = (self.centres[self.vortNum], self.strengths[self.vortNum], self.radius[self.vortNum], self.clockwise[self.vortNum])
            else:
                data = (self.centres[self.vortNum], self.strengths[self.vortNum])

            self.vortNum += 1

        return data 

def main():
    # Setup coordinates grid
    coordGrids = makeGrid([xSide,ySide],[xNumCells,yNumCells])

    # Initialise object to store data about multiple vortices of different types
    VortexDefs = Vortices(vortType, vortO, vortS, [5], [-1])

    # Place vortices in domain
    velGrids, rho, p = defineVortices(VortexDefs, coordGrids)

    # Plot only vortex fields
    plotVelocity(coordGrids, velGrids)
    plt.show()

    # Plot and save fields
    #plotAll(coordGrids, velGrids, rho, p, 'Big_test.pdf')

'''
Generic multiple vortices function
Outputs velocity components, density and pressure as meshgrids
vortDefs - Vortices object; coordGrids - meshgrid for coordinates;  axialVel - uniform axial velocity to be applied
'''
def defineVortices(vortDefs, coordGrids, axialVel=1):
    # Intialise 3D arrays to store multiple meshgrids - one for the component effect of each vortex
    tComps = np.zeros(np.append(coordGrids[:,:,0].shape, vortDefs.strengths.shape[0]))    # Is temperature effect generic? Haven't implemented general method for doing pressure and density effects
    uComps = tComps.copy()
    vComps = tComps.copy()

    # Dictionary mapping for functions - will be faster then multiple if/else statements - also more readable code
    vortexType = {1:isoVortex, 2:loVortex, 3:solidVortex}

    # Loop through given vortices and calculate their effect on each cell of the grid
    for i in range(vortDefs.strengths.shape[0]):
        # Get function for this vortex type
        func = vortexType.get(vortDefs.types[i])
        # Call vortex function to fill component arrays - with data for a single vortex
        tComps[:,:,i], uComps[:,:,i], vComps[:,:,i] = func(vortDefs.getNextVortex(), coordGrids)

    # Collate effects of each vortex
    T = 1 - np.sum(tComps,axis=2)
    U = np.sum(uComps,axis=2)
    V = np.sum(vComps,axis=2)
    rho = np.power(T, 1/(GAMMA-1)) # Simple isentropic relations for rho and p
    p = rho*T

    # Add uniform axial velocity field? Or have some other equation for it
    W = np.ones(U.shape)

    # Stack velocity grids, and do other necessary cleanup before output
    velGrids = np.dstack([U,V,W])

    return velGrids, rho, p

'''
Function for outputting the effect of a simple isentropic vortex on the domain
- edge of grid currently not a bounding wall, ie vortex is acting like it's in an infinite domain and grid is just a smple of this
vortData - tuple produced by getNextVortex() function of Vortices class
'''
def isoVortex(vortData,coordGrids):
    # Get radius of each cell from centre of this vortex
    r = np.sqrt((coordGrids[:,:,0]-vortData[0][0])**2 + (coordGrids[:,:,1] - vortData[0][1])**2)

    # Take into account temperature effect of this vortex - and consequently, its effect on pressure and density
    tComp = (GAMMA-1) * vortData[1]**2 * np.exp(1-r**2) / (8*GAMMA*np.pi**2)

    # Velocity components due to this vortex
    uComp = (vortData[1]/(2*np.pi)) * np.exp(0.5*(1-r**2)) * (coordGrids[:,:,1] - vortData[0][1])
    vComp = (vortData[1]/(2*np.pi)) * np.exp(0.5*(1-r**2)) * (vortData[0][0] - coordGrids[:,:,0])

    return tComp, uComp, vComp

'''
Function for outputting the effect of a Lamb-Oseen vortex
- using adapted equations from StreamVane paper, generalised for an arbitrary number of vortices at arbitrary positions
------- Equations from StreamVane paper has slight error for v velocity component, correct here
- Generic twin swirl profile can be created with by placing two vortices with parameters from paper
- currently has dummy effect on temp

vortData - tuple produced by getNextVortex() function of Vortices class
'''
def loVortex(vortData,coordGrids):
    # Some parameter - a0=0.1 apparently creates maximum swirl angle of 16deg when other parameters same as paper
    a0 = 0.1

    # Get radius of each cell from centre of this vortex
    r = np.sqrt((coordGrids[:,:,0]-vortData[0][0])**2 + (coordGrids[:,:,1] - vortData[0][1])**2)
    rsquare = r**2

    # Get omega (actual vortex strength?)
    omega = vortData[1]/(np.pi * a0**2)

    # Velocity components due to this vortex
    uComp = 0.5 * (a0**2 * omega * (coordGrids[:,:,1] - vortData[0][1]) / r**2) * (1 - np.exp(-r**2/a0**2))
    vComp = 0.5 * (a0**2 * omega * (vortData[0][0] - coordGrids[:,:,0]) / r**2) * (1 - np.exp(-r**2/a0**2))

    # Dummy temp grid, no temperature effect currently
    tComp = np.zeros(uComp.shape)

    return tComp, uComp, vComp

'''
Function for outputting the effect of a forced vortex
- currently has dummy field for both density and pressure, only velocity vector fields at the moment
- linear increase in swirl angle from center to outer edge
- swirl angle defined as angle between resultant vector and axial velocity component (cause by tangential component velocity)
- solid/forced vortex - not realistic; ie instantaneously created vortex, no effect on cells outside it's radius
- also assumes no radial velocity

vortData - tuple produced by getNextVortex() function of Vortices class
'''
def solidVortex(vortData, coordGrids):
    # Get swirl angle and convert it to radians
    maxSwirlAngle = np.deg2rad(vortData[1])

    # Get axial coordinates
    r = np.sqrt((coordGrids[:,:,0]-vortData[0][0])**2 + (coordGrids[:,:,1] - vortData[0][1])**2)
    theta = np.arctan(coordGrids[:,:,1]/coordGrids[:,:,0])

    # Normalise radius for straightforward angle calculation and set cells outside vortex size to 0
    rNorm = r/vortData[2]
    rNorm[(rNorm > 1)] = 0

    # Get swirl angle distribution
    swirlAngles = maxSwirlAngle*rNorm

    # Transform so swirl is coherent (either clockwise or anticlockwise) - without this, the swirl profile produced is mirrored about the y axis
    swirlAngles[(coordGrids[:,:,0] * -vortData[3] < 0)] = swirlAngles[(coordGrids[:,:,0] * -vortData[3] < 0)] * -1

    # Get tangential velocity at each cell
    tangentVel = axialVel*np.tan(swirlAngles)

    # Get velocity vector components, in-plane cartesian (assume no radial velocity)
    uComp = -r*tangentVel*np.sin(theta)
    vComp =  r*tangentVel*np.cos(theta)

    # Dummy temp grid, no temperature effect modelled currently
    tComp = np.zeros(uComp.shape)

    return tComp, uComp, vComp

'''
Make meshgrids to store coordinate system; as a result, all variable fields will be meshgrids also, good performance since using numpy matrix operations
'''
def makeGrid(sideLengths=[10,10],numCells=[100,100]):
    # Get the resulting side length of the grid cells
    cellSides = np.divide(sideLengths,numCells)

    # Create coordinate system from mesh info
    x = np.linspace(-sideLengths[0]/2+cellSides[0]/2, sideLengths[0]/2-cellSides[0]/2, numCells[0])
    y = np.linspace(-sideLengths[1]/2+cellSides[1]/2, sideLengths[1]/2-cellSides[1]/2, numCells[1])

    # Create meshgrid to store coordinate grid - useful in this form for plotting later and reduces the amount of for loops since we can use numpy matrix operations instead
    X, Y = np.meshgrid(x,y, indexing='xy')        # Use familiar ij matrix indexing

    # Stack grids into a 3D array, for convenience when passing between functions - not sure about performance effect, better or worse or negligible?
    coordGrids = np.dstack([X,Y])

    return coordGrids

'''
Utility for showing and saving all plots
'''
def plotAll(coordGrids, velGrids, rho, p, pdfName):
    plotVelocity(coordGrids, velGrids)

    plotThermos(coordGrids,rho,p)

    plotSwirl(coordGrids, getSwirl(coordGrids,velGrids))

    saveFigsToPdf(pdfName)

    plt.show()

'''
Create plots for the swirling velocity profile as a quiver plot and a streamlines plot
'''
def plotVelocity(coordGrids, velGrids):
    # Get individual grids, more convenient for pruning etc. for plots
    X = coordGrids[:,:,0]
    Y = coordGrids[:,:,1]

    # For making quiver plot sparser, but independent of grid density
    quiverEvery = int(X.size**(1/2) / 20)
    quiverEvery = (1 if (quiverEvery == 0) else quiverEvery)

    # Make quiver plot
    plt.figure()
    plt.gca().set_aspect('equal', adjustable='box')
    plt.title("Quiver")
    skip = (slice(None,None,quiverEvery), slice(None,None,quiverEvery))       # Prune output so quiver plot is not so dense
    plt.quiver(coordGrids[:,:,0][skip], coordGrids[:,:,1][skip], velGrids[:,:,0][skip], velGrids[:,:,1][skip], scale = 5)

    # Make streamlines plot
    plt.figure()
    plt.gca().set_aspect('equal', adjustable='box')
    plt.title("Streamlines")
    plt.streamplot(coordGrids[1,:,0], coordGrids[:,1,1], velGrids[:,:,0], velGrids[:,:,1], density=2)            # streamplot uses vector axis for xy instead of meshgrid for some reason?

'''
Create contour plots for density and pressure field
'''
def plotThermos(coordGrids,rho,p):
    plt.figure()
    plt.title('Density')
    plt.contourf(coordGrids[:,:,0],coordGrids[:,:,1],rho,100)
    plt.colorbar()

    plt.figure()
    plt.title('Pressure')
    plt.contourf(coordGrids[:,:,0],coordGrids[:,:,1],p,100)
    plt.colorbar()

'''
Get swirl angles
'''
def getSwirl(coordGrids, velGrids):
    # Get tangential velocity
    velTheta = (coordGrids[:,:,0]*velGrids[:,:,1] - velGrids[:,:,0]*coordGrids[:,:,1]) / (coordGrids[:,:,0]**2 + coordGrids[:,:,1]**2)
    # Get swirl angle - as defined in literature
    swirlAngle = np.arctan(velTheta/velGrids[:,:,2])
    # Convert to degrees
    swirlAngle = np.rad2deg(swirlAngle)

    return swirlAngle

'''
Calculate Root Mean Square error between two swirl angle profiles, for comparisons
'''
def getError(swirlAngle, desiredSwirl):
    RMSE = np.sqrt((1/np.size(swirlAngle))*np.sum((swirlAngle-desiredSwirl)**2))

    return RMSE

'''
Create contour plot for swirl angle
'''
def plotSwirl(coordGrids, swirlAngle):
    # Make contour plot
    plt.figure()
    plt.title('Swirl angle')
    plt.contourf(coordGrids[:,:,0],coordGrids[:,:,1],swirlAngle,100)
    plt.colorbar()

'''
Save all current figures into a multi-page pdf
'''
def saveFigsToPdf(outputFile):
    from matplotlib.backends.backend_pdf import PdfPages

    with PdfPages(outputFile) as pdf:
        # Go through all active figures and save to a separate pdf page
        for fig in range(1, plt.gcf().number+1):
            pdf.savefig(fig)


if __name__ == '__main__':
    main()