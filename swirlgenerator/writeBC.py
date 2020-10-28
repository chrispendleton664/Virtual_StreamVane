# -------------------------------------------------------------------------------------------------------------------
#
# Module for writing the raw flow field data into a file format which can be uses as an inlet boundary condition
# Takes a FlowField object as input
#
# Assumes that the bottom corner of the inlet surface is at [0,0,0]
#
# Only has functionality for SU2 format at the moment
#
# -------------------------------------------------------------------------------------------------------------------

import swirlgenerator.core as sg
import numpy as np

'''
For writing into an SU2 useable format
Column format for SU2 inlet boundary condition (3D domain):
x, y, z, temperature, velocity magnitude,  x-component of flow direction unit vector, y-component of flow direction unit vector, z-component of flow direction unit vector
'''
def writeSU2(flowField: sg.FlowField, filename):
    # Flatten coordinate grids to get coordinates of every cell
    x = np.reshape(flowField.coordGrids[:,:,0],[flowField.coordGrids[:,:,0].size,1])
    y = np.reshape(flowField.coordGrids[:,:,1],[flowField.coordGrids[:,:,1].size,1])

    # Transform domain coordinates from 0,0 centered to match mesh
    x = x + flowField.sideLengths[0]/2
    y = y + flowField.sideLengths[1]/2
    
    # Z coords at end of mesh
    z = np.zeros(np.shape(x))

    # Flatten velocity grids to get velocity vector of every cell
    u = np.reshape(flowField.velGrids[:,:,0],[flowField.velGrids[:,:,0].size,1])
    v = np.reshape(flowField.velGrids[:,:,1],[flowField.velGrids[:,:,1].size,1])
    w = np.reshape(flowField.velGrids[:,:,2],[flowField.velGrids[:,:,2].size,1])
    velVec = np.column_stack((u,v,w))

    # Get velocity magnitude at each cell
    velMag = np.linalg.norm(velVec, axis=1)

    # Get direction components of unit vector
    flowDirection = velVec/ np.column_stack((velMag,velMag,velMag))

    # Dummy temperature array
    temp = np.zeros(np.shape(x))

    # Collate into one matrix for writing into file
    boundaryConditions = np.column_stack((x,y,z,temp,velMag,flowDirection))

    # Write to file
    with open(filename, 'w') as f:
        # Write metadata
        f.write(f"NMARK= 1\nMARKER_TAG= inlet\nNROW= {np.shape(boundaryConditions)[0]}\nNCOL= 8\n")

        # Write each row of matrix into separate line with 6 decimal places
        np.savetxt(f, boundaryConditions, fmt='%.6f')

'''
Wrapper function for calling the correct function to write the boundary condition in the requested format
'''
def writeInlet(InputObject: sg.Input, flowField: sg.FlowField):
    # Dictionary mapping for different boundary condition file formats and functions which create them
    formats = {'su2':writeSU2}

    # Write boundary condition in correct format
    func = formats.get(InputObject.format)
    func(flowField, InputObject.filename)