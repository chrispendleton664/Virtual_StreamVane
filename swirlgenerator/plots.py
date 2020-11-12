# -------------------------------------------------------------------------------------------------------------------
#
# Module for plotting the variables contained in the flow field object
#
# -------------------------------------------------------------------------------------------------------------------

import swirlgenerator.core as sg
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import plotly.graph_objects as go
import plotly.figure_factory as ff
import plotly.subplots as subplt
import numpy as np

'''
Utility for showing and saving all plots
'''
def plotAll(flowfield: sg.FlowField, pdfName=None):
    # Initialise array to store the figure objects
    figs = []

    # Plot velocity profiles as quiver and streamplots
    figs = plotVelocity(flowfield, figs)

    #plotThermos(flowfield)

    # Plot swirl angle
    figs = plotSwirl(flowfield, figs)

    # If saving, don't show plots
    if (pdfName != None):
        __saveFigsToPdf__(pdfName)
    #else:
        # plt.show()
        #fig.show()

    # Write figures to html file
    __figs_to_html__(figs,'test.html')

'''
Create plots for the swirling velocity profile as a quiver plot and a streamlines plot
'''
def plotVelocity(flowfield, figs, maxNumArrows=30):
    # Reduced indices of grids - so only a maximum of n^2 arrows will be plotted on the quiver plot
    n = maxNumArrows
    gridDims = flowfield.coordGrids.shape
    step = [int(gridDims[0]/n),int(gridDims[1]/n)]
    step = [1 if s == 0 else s for s in step]               # Protection for when number of actual data points is less than maxNumArrows
    reduced = np.mgrid[0:gridDims[0]:step[0],0:gridDims[1]:step[1]].reshape(2,-1).T

    # Make quiver plot
    quiver = ff.create_quiver(x=flowfield.coordGrids[reduced[:,0],reduced[:,1],0], y=flowfield.coordGrids[reduced[:,0],reduced[:,1],1], 
                            u=flowfield.velGrids[reduced[:,0],reduced[:,1],0], v=flowfield.velGrids[reduced[:,0],reduced[:,1],1],
                            scale=10)
    # Add to list of figures
    figs.append(quiver)

    # Flatten mesh grids into lists ------------- temporary solution
    #x = flowfield.coordGrids[:,:]
    # Make streamlines plot
    streamlines = ff.create_streamline(x=flowfield.axis[0], y=flowfield.axis[1], 
                                        u=flowfield.velGrids[:,:,0], v=flowfield.velGrids[:,:,1])
    # Add to list of figures
    figs.append(streamlines)

    return figs

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
def plotSwirl(flowfield, figs):
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

    # Make contour plot   
    figs.append(go.Figure(data=
                        go.Contour(
                            z=flowfield.swirlAngle, x=flowfield.axis[0], y=flowfield.axis[1],
                            colorscale='Jet',
                            contours=dict(start=minVal, end=maxVal, size=0.01, showlines=False))
                        ))

    return figs

'''
Save all current figures into a multi-page pdf
'''
def __saveFigsToPdf__(outputFile):
    with PdfPages(outputFile) as pdf:
        # Go through all active figures and save to a separate pdf page
        for fig in range(1, plt.gcf().number+1):
            pdf.savefig(fig)

'''
Saves list of plotly figures in a html file
Adapted from https://stackoverflow.com/a/59265030, answer by 'ttekampe'
'''
def __figs_to_html__(figs, filename):
    import plotly.offline as pyo

    with open(filename, 'w') as dashboard:
        # Page headings
        dashboard.write('<html><head></head><body>' + '\n')

        # Plotly.js code needs to be included, but only once, to make the plots interactive
        add_js = True

        for fig in figs:
            # Get the html version of the plotly fig, as just a div and not a full html file
            inner_html = pyo.plot(fig, include_plotlyjs=add_js, output_type='div')

            # Write to output html
            dashboard.write(inner_html)
            
            # Change flag to false so that plotly.js source code is only added once per file
            add_js = False

        # Close off html file
        dashboard.write('</body></html>' + '\n')
