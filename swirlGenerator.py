import numpy as np
from matplotlib import pyplot as plt
import warnings

'''
Heavy use of numpy for fast matrix and vector operations
Matplotlib for doing visualisations
'''

# Plots output pdf file
outputPdf = 'generic_twin_swirl.pdf'

# Fluid parameters
GAMMA = 1.4
kin_visc = 1.81e-5

# Uniform streamwise velocity
axialVel = 5

# Vortex centers
vortO = np.array([[-1.25,0],[1.25,0]])

# Vortex strengths
vortS = np.array([1,-1])

# Mesh variables
xSide = 10                 # Side length of grid in x direction
ySide = 10                 # Side length of grid in y direction
xNumCells = 1000           # Num of cells in x direction
yNumCells = 1000           # Num of cells in y direction


def main():

    # Setup coordinates grid
    coordGrids = makeGrid([xSide,ySide],[xNumCells,yNumCells])

    # Calculate flow fields -- all have the same interface, probs could make a generic function instead of 3 separate ones
    #velGrids, rho, p = isoVortices(vortO,vortS,coordGrids,axialVel)
    #velGrids, rho, p = bulkSwirl(5, 18,coordGrids)
    velGrids, rho, p = loVortices(vortO,vortS,coordGrids,axialVel)
    
    # Create swirled profiles
    plotVelocity(coordGrids,velGrids)

    # Create contour for density and pressure
    plotThermos(coordGrids,rho,p)

    # Calculate swirl angle of flow, for comparisons
    swirlAngle = getSwirl(coordGrids, velGrids)

    # Create contour of swirl angle
    plotSwirl(coordGrids,swirlAngle)

    # Save figures to pdf - this needs to be done before the plots are shown, does plt.show() clear the figures object??
    saveFigsToPdf(outputPdf)

    # Draw to screen
    plt.show()

'''
Generic multiple vortices function
'''
def defineVortices(centres, strengths, coordGrids, axialVel=1):
    pass
    # Intialise 3D arrays to store multiple meshgrids - one for the component effect of each vortex

    # Loop through given vortices and calculate their effect on each cell of the grid

    '''
    Some method of specifying the effect of this vortex on each cell of the grid
    Probably here we would make it modular and have the option of using different vortex types? (Lamb-Oseen, isentropic, forced)
    '''

    # Collate effects of each vortex - this probably needs to be modular too, different for each vortex type?

    # Add uniform axial velocity field? Or have some other equation for it

    # Stack velocity grids, and do other necessary cleanup before output


'''
For multiple simple isentropic vortices
- with uniform streamwise velocity
- edge of grid currently not a bounding wall, ie vortex is acting like it's in an infinite domain and grid is just sample of this
'''
def isoVortices(vortCentres, vortStrengths, coordGrids, axialVel=1):
    # Initialise 3D arrays to store multiple meshgrids
    tComps = np.zeros(np.append(coordGrids[:,:,0].shape, vortStrengths.size))
    uComps = tComps.copy()
    vComps = tComps.copy()

    # Loop through specified vortices and calculate their effect on each cell of the grid
    for i in range(vortStrengths.size):
        # Set radius of each cell from centre of vortex
        r = np.sqrt(np.square(coordGrids[:,:,0]-vortCentres[i,0]) + np.square(coordGrids[:,:,1] - vortCentres[i,1]))

        # Dummy parameter for taking into account temperature effect of multiple vortices on each cell (and consequently the rho and pi)
        tComps[:,:,i] = vortStrengths[i]**2 * np.exp(1-np.square(r))

        # Velocity components due to this vortex
        uComps[:,:,i] = (vortStrengths[i]/(2*np.pi)) * np.multiply(np.exp(0.5*(1-np.square(r))), vortCentres[i,1]-coordGrids[:,:,1])
        vComps[:,:,i] = (vortStrengths[i]/(2*np.pi)) * np.multiply(np.exp(0.5*(1-np.square(r))), coordGrids[:,:,0]-vortCentres[i,0])

    # Collate effect of all vortices on temp
    T = 1 - (((GAMMA-1)/(8*GAMMA*np.pi**2)) * np.sum(tComps,axis=2))
    # Check if any of T is negative - unrealistic, suggests that one of the betas specified for the vortices is too high
    if np.any(T < 0):
        warnings.warn('Vortex strength values specified may be too high resulting in negative temperatures, pressure and density fields produced may be invalid', RuntimeWarning)
    
    # Now calc density and pressure using temperature
    rho = np.power(T, 1/(GAMMA-1))
    p = np.multiply(rho,T)
    
    # Sum velocity effects of each vortex on the velocities
    U = np.sum(uComps,axis=2)
    V = np.sum(vComps,axis=2)

    # Apply uniform axial velocity profile - is this reasonable??
    W = np.ones(U.shape)*axialVel

    # Stack velocity vector components grids, similar to coords
    velGrids = np.dstack([U,V,W])

    return velGrids,rho,p

'''
For Bulk swirl profile
- currently has dummy field for both density and pressure, only velocity vector fields at the moment
- linear increase in swirl angle from center to outer edge
- swirl angle defined as angle between resultant vector and axial velocity component (cause by tangential component velocity)
- solid/forced vortex - not realistic; ie instantaneously created vortex, no effect on cells outside it's radius
- also assumes no radial velocity
'''
def bulkSwirl(vortR, maxSwirlAngleDegrees, coordGrids, axialVel=1, clockwise = False):
    # Convert to radians
    maxSwirlAngle = np.deg2rad(maxSwirlAngleDegrees)

    # Get axial coordinates
    r = np.sqrt(np.square(coordGrids[:,:,0]) + np.square(coordGrids[:,:,1]))
    theta = np.arctan(coordGrids[:,:,1]/coordGrids[:,:,0])
    
    # Normalise radius for straightforward swirl angle calculation and set cells outside vortex size to 0
    rNorm = r/vortR
    rNorm[(rNorm > 1)] = 0

    # Get swirl angle
    swirlAngles = maxSwirlAngle*rNorm

    # Transform so swirl is coherent (either clockwise or anticlockwise) - without this, the swirl profile produced is mirrored about the y axis
    swirlAngles[(coordGrids[:,:,0] * (-1 if clockwise else 1) < 0)] = swirlAngles[(coordGrids[:,:,0] * (-1 if clockwise else 1) < 0)] * -1

    # Get tangential velocity at each cell
    tangentVel = axialVel*np.tan(swirlAngles)

    # Get velocity vector components, in-plane cartesian (assume no radial velocity)
    U = -r*tangentVel*np.sin(theta)
    V =  r*tangentVel*np.cos(theta)
    # Apply uniform axial vel distribution
    W = np.ones(np.shape(U))*axialVel

    # Stack velocity vector components grids, similar to coords
    velGrids = np.dstack([U,V,W])

    # Create dummy fields for density and pressure
    rho = np.ones(np.shape(U))
    p = np.ones(np.shape(U))

    return velGrids,rho,p

'''
For multiple Lamb-Oseen vortices
- using adapted equations from StreamVane paper, generalised for an arbitrary number of vortices at arbitrary positions
------- Equations from StreamVane paper has slight error for v velocity component, correct here
- Generic twin swirl profile can be created with by placing two vortices with parameters from paper
- currently has dummy grids for the density and pressure fields
'''
def loVortices(vortCentres, vortStrengths, coordGrids, axialVel=1):
    # Some parameter - 0.1 creates maximum swirl angle of 16deg when other parameters same as paper
    a0 = 0.1
    a0square = a0**2

    # Initialise 3D arrays to store multiple meshgrids (one for each vortex)
    uComps = np.zeros(np.append(coordGrids[:,:,0].shape, vortStrengths.size))
    vComps = uComps.copy()

    # Loop through given vortices and calculate their effect on each cell of the grid
    for i in range(vortStrengths.size):
        # Get radius of each cell from centre of this vortex
        r = np.sqrt(np.square(coordGrids[:,:,0]-vortCentres[i,0]) + np.square(coordGrids[:,:,1] - vortCentres[i,1]))

        # Get omega (actual vortex strength?)
        omega = vortStrengths[i]/(np.pi * a0square)

        # Velocity components due to this vortex
        uComps[:,:,i] = 0.5 * (a0square * omega * (coordGrids[:,:,1] - vortCentres[i,1]) / np.square(r)) * (1 - np.exp(-np.square(r)/a0square))
        vComps[:,:,i] = 0.5 * (a0square * omega * (vortCentres[i,0] - coordGrids[:,:,0]) / np.square(r)) * (1 - np.exp(-np.square(r)/a0square))

    # Sum velocity effects of each vortex on the velocities
    U = np.sum(uComps,axis=2)
    V = np.sum(vComps,axis=2)

    # Apply uniform axial velocity profiles
    W = np.ones(U.shape)*axialVel
    
    # Stack velocity vector components grids, similar to coords
    velGrids = np.dstack([U,V,W])

    # Create dummy fields for density and pressure
    rho = np.ones(np.shape(U))
    p = np.ones(np.shape(U))

    return velGrids,rho,p
    
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
Calculate Root Mean Square error between two swirl angle profiles, for comparisons - not tested yet
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