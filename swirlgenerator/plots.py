# -------------------------------------------------------------------------------------------------------------------
#
# Module for plotting the variables contained in the flow field object
#
# -------------------------------------------------------------------------------------------------------------------

import swirlgenerator.core as sg
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

'''
Utility for showing and saving all plots
'''
def plotAll(flowfield: sg.FlowField, pdfName=None):
    plotVelocity(flowfield)

    #plotThermos(flowfield)

    plotSwirl(flowfield)

    # If saving, don't show plots
    if (pdfName != None):
        __saveFigsToPdf__(pdfName)
    else:
        plt.show()

    # Clear figures when done
    plt.close('all')

'''
Create plots for the swirling velocity profile as a quiver plot and a streamlines plot
'''
def plotVelocity(flowfield, maxNumArrows=30):
    # Reduced indices of grids - so only a maximum of n^2 arrows will be plotted on the quiver plot
    n = maxNumArrows
    gridDims = flowfield.coordGrids.shape
    step = [int(gridDims[0]/n),int(gridDims[1]/n)]
    step = [1 if s == 0 else s for s in step]               # Protection for when number of actual data points is less than maxNumArrows
    reduced = np.mgrid[0:gridDims[0]:step[0],0:gridDims[1]:step[1]].reshape(2,-1).T

    # Make quiver plot
    plt.figure()
    plt.gca().set_aspect('equal', adjustable='box')
    plt.title("Quiver")
    plt.quiver(flowfield.coordGrids[reduced[:,0],reduced[:,1],0], flowfield.coordGrids[reduced[:,0],reduced[:,1],1], flowfield.velGrids[reduced[:,0],reduced[:,1],0], flowfield.velGrids[reduced[:,0],reduced[:,1],1], units='dots', width=2,headwidth=5,headlength=5,headaxislength=2.5)

    # Correct velocity grids for circular domains - only needed since streamplot function takes axis as input rather than meshgrid coordinates
    correctedVel = flowfield.velGrids[:,:,0:2]
    correctedVel[np.dstack([flowfield.outside,flowfield.outside])] = 0

    # Make streamlines plot
    plt.figure()
    plt.gca().set_aspect('equal', adjustable='box')
    plt.title("Streamlines")
    # Need to take axis points from middle of grids in case the domain is circular
    plt.streamplot(flowfield.coordGrids[int(gridDims[0]/2),:,0], flowfield.coordGrids[:,int(gridDims[1]/2),1], correctedVel[:,:,0], correctedVel[:,:,1], density=2)            # streamplot uses vector axis for xy instead of meshgrid for some reason?

'''
Create contour plots for density and pressure field
'''
def plotThermos(flowfield):
    plt.figure()
    plt.title('Density')
    plt.contourf(flowfield.coordGrids[:,:,0],flowfield.coordGrids[:,:,1],flowfield.rho,100,cmap='jet')
    plt.colorbar()

    plt.figure()
    plt.title('Pressure')
    plt.contourf(flowfield.coordGrids[:,:,0],flowfield.coordGrids[:,:,1],flowfield.pressure,100,cmap='jet')
    plt.colorbar()

'''
Create contour plot for swirl angle
'''
def plotSwirl(flowfield):
    # Get maximum magnitude of swirl angle
    maxMag = max([abs(flowfield.swirlAngle.min()), flowfield.swirlAngle.max()])
    # Choose appropriate 'round to the nearest'
    if maxMag < 5:
        rounding = 1
    elif maxMag < 10:
        rounding = 5
    else:
        rounding = 10
    # Round max/min values to create range of swirl angles
    minVal = np.floor(flowfield.swirlAngle.min() / rounding) * rounding
    maxVal = np.ceil(flowfield.swirlAngle.max()  / rounding) * rounding

    # Make ticks for colormap
    #ticks = np.arange(minVal,maxVal,rounding)
    ticks = np.linspace(minVal,maxVal,11)

    # Make contour plot
    plt.figure()
    plt.title('Swirl angle')
    plt.contourf(flowfield.coordGrids[:,:,0],flowfield.coordGrids[:,:,1],flowfield.swirlAngle,100,cmap='jet',vmin=minVal,vmax=maxVal)
    plt.colorbar(ticks=ticks)

'''
Save all current figures into a multi-page pdf
'''
def __saveFigsToPdf__(outputFile):
    with PdfPages(outputFile) as pdf:
        # Go through all active figures and save to a separate pdf page
        for fig in range(1, plt.gcf().number+1):
            pdf.savefig(fig)