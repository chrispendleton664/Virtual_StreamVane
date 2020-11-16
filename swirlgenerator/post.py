# -------------------------------------------------------------------------------------------------------------------
#
# Module for plotting the variables contained in the flow field object
#
# -------------------------------------------------------------------------------------------------------------------

import core as sg
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np


def plotAll(flowfield: sg.FlowField, pdfName=None):
    '''
    Utility for showing and saving all plots
    - flowfield - FlowField object containing data to be plotted
    - pdfName - filename to save the plots to, should include extension
    '''

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


def plotVelocity(flowfield, maxNumArrows=40000):
    '''
    Create plots for the swirling velocity profile as a quiver plot and a streamlines plot
    - flowfield - FlowField object containing data to be plotted
    - maxNumArrows - maximum number of arrows to draw on the quiver plot
    '''

    # Reduced indices - so only a maximum of n arrows will be plotted on the quiver plot
    n = maxNumArrows
    numPoints = flowfield.coords.size
    step = int(numPoints/n)
    step = (1 if step == 0 else step)              # Protection for when number of actual data points is less than maxNumArrows
    reducedIdxs = np.arange(start=0,stop=numPoints,step=step)

    # Make quiver plot
    plt.figure()
    plt.gca().set_aspect('equal', adjustable='box')
    plt.title("Quiver")
    plt.quiver(flowfield.coords[reducedIdxs].real, flowfield.coords[reducedIdxs].imag, flowfield.velocity[reducedIdxs,0], flowfield.velocity[reducedIdxs,1], units='dots', width=2,headwidth=5,headlength=5,headaxislength=2.5)
    plt.axis('off')
    # Draw boundary
    plt.plot(flowfield.boundaryCurve.real, flowfield.boundaryCurve.imag,'k-')

    # Need to format list of velocity vectors into a meshgrid for plotting with streamplot
    size = [int(np.sqrt(flowfield.velocity[:,0].size))]*2
    uGrid = flowfield.velocity[:,0].reshape(size)
    vGrid = flowfield.velocity[:,1].reshape(size)

    # Make streamlines plot
    plt.figure()
    plt.gca().set_aspect('equal', adjustable='box')
    plt.title("Streamlines")
    plt.streamplot(flowfield.axis[0], flowfield.axis[1], uGrid, vGrid, density=2)            # streamplot uses vector axis for xy instead of meshgrid for some reason?
    plt.axis('off')
    # Draw boundary
    plt.plot(flowfield.boundaryCurve.real, flowfield.boundaryCurve.imag,'k-')


def plotThermos(flowfield):
    '''
    Create contour plots for density and pressure field
    - flowfield - FlowField object containing data to be plotted
    '''

    plt.figure()
    plt.title('Density')
    plt.contourf(flowfield.coords.real,flowfield.coords.imag,flowfield.rho,100,cmap='jet')
    plt.colorbar()

    plt.figure()
    plt.title('Pressure')
    plt.contourf(flowfield.coords.real,flowfield.coords.imag,flowfield.pressure,100,cmap='jet')
    plt.colorbar()


def plotSwirl(flowfield):
    '''
    Create contour plot for swirl angle
    - flowfield - FlowField object containing data to be plotted
    '''

    # Convert nans to zero so that max/min operations don't result in nan
    swirlAngle = np.nan_to_num(flowfield.swirlAngle)

    # Get maximum magnitude of swirl angle
    maxMag = max([abs(swirlAngle.min()), swirlAngle.max()])
    # Choose appropriate 'round to the nearest'
    if maxMag < 5:
        rounding = 1
    elif maxMag < 10:
        rounding = 5
    else:
        rounding = 10
    # Round max/min values to create range of swirl angles
    minVal = np.floor(swirlAngle.min() / rounding) * rounding
    maxVal = np.ceil(swirlAngle.max()  / rounding) * rounding

    # Make ticks for colormap
    ticks = np.linspace(minVal,maxVal,11)

    # Reshape swirl angle to a meshgrids for contour plot input
    size = [int(np.sqrt(flowfield.swirlAngle.size))]*2
    swirlAngle = flowfield.swirlAngle.reshape(size)

    # Make contour plot
    plt.figure()
    plt.gca().set_aspect('equal', adjustable='box')
    plt.title('Swirl angle')
    # For some reason contourf doesn't like when the coordinate grids have nans in them, so using zero instead of nan versions of array
    plt.contourf(flowfield.axis[0],flowfield.axis[1],swirlAngle,100,cmap='jet',vmin=minVal,vmax=maxVal)
    plt.colorbar(ticks=ticks)
    plt.axis('off')
    # Draw boundary
    plt.plot(flowfield.boundaryCurve.real, flowfield.boundaryCurve.imag,'k-')


def __saveFigsToPdf__(outputFile):
    '''
    Save all current figures into a multi-page pdf
    - Internal function, should not be used outside plots.py
    -outputFile - name of file to save to, including extension
    '''

    with PdfPages(outputFile) as pdf:
        # Go through all active figures and save to a separate pdf page
        for fig in range(1, plt.gcf().number+1):
            pdf.savefig(fig)