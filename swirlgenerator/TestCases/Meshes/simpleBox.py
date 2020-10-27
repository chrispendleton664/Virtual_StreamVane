# --------------------------------------------------------------------------------
#
# Script to generate a simple 3D rectangular geometry and mesh
# for testing the development of swirling profiles created by swirlGenerator
# using the python API of GMSH
# https://gitlab.onelab.info/gmsh/gmsh/-/blob/master/tutorial/python
#
# --------------------------------------------------------------------------------

import gmsh
import sys

gmsh.initialize()

# Output file - Gmsh guesses the correct format to write the mesh in with the file extension
outputFile = "box.su2"

# Box dimensions
HEIGHT  = 1
WIDTH   = 1
LENGTH  = 5

# Mesh definition
NUM_CELLS_HEIGHT    = 100
NUM_CELLS_WIDTH     = 100
NUM_CELLS_LENGTH    = 500

'''
GEOMETRY DEFINITION
'''
# Get vertices of box
vertices = [[0,0,0], [WIDTH,0,0], [WIDTH,HEIGHT,0], [0,HEIGHT,0], [0,0,LENGTH], [WIDTH,0,LENGTH], [WIDTH,HEIGHT,LENGTH], [0,HEIGHT,LENGTH]] 

# Add a new geometry model
gmsh.model.add("box")

# Define the points of the box
points = [gmsh.model.geo.addPoint(point[0],point[1],point[2]) for point in vertices]

# Define the lines between the points - Needs to be done in a certain order for later
lines = []
for i in range(4):
    # Lines in first square face
    lines.append(gmsh.model.geo.addLine(points[i], points[i+1]) if (i != 3) else gmsh.model.geo.addLine(points[i], points[i-3]))  
    
for i in range(4,8):
    # Lines in second square face
    lines.append(gmsh.model.geo.addLine(points[i], points[i+1]) if (i != 7) else gmsh.model.geo.addLine(points[i], points[i-3]) )

for i in range(4):
    # Lines connecting two square faces
    lines.append(gmsh.model.geo.addLine(points[i], points[i+4]))

# Define the curve loops of lines 
loops = []
# Loops for square faces
loops.append(gmsh.model.geo.addCurveLoop(lines[0:4]))
loops.append(gmsh.model.geo.addCurveLoop(lines[4:8]))
# Loops for faces aloing box length - order of points within lines are important, so need negative signs when going 'backwards' along a line
for i in range(4):
    loops.append(gmsh.model.geo.addCurveLoop([lines[i], (lines[i+9] if (i != 3) else lines[8]), -lines[i+4], -lines[i+8]]))

# Define plane surfaces with these curve loops
surfaces = [gmsh.model.geo.addPlaneSurface([loop]) for loop in loops]

# Define a 'surface loop' then use this to define a volume
surfLoop = gmsh.model.geo.addSurfaceLoop(surfaces)
vol = gmsh.model.geo.addVolume([surfLoop])

'''
MESH DEFINITION
'''
# Define a transfinite(uniform) mesh on the lines
for i in range(len(lines)):
    if (i > 7):
        numCells = NUM_CELLS_LENGTH
    elif (i % 2 == 0):
        numCells = NUM_CELLS_WIDTH
    else:
        numCells = NUM_CELLS_HEIGHT

    gmsh.model.geo.mesh.setTransfiniteCurve(lines[i], numCells+1)   # Gmsh defined meshes with points rather than cells

# Define transfinite mesh on surfaces
[gmsh.model.geo.mesh.setTransfiniteSurface(surface) for surface in surfaces]

# Recombine triangular mesh into quadrilaterals
[gmsh.model.geo.mesh.setRecombine(2, surface) for surface in surfaces]

# Define the volume as a transfinite mesh also
gmsh.model.geo.mesh.setTransfiniteVolume(vol)

# Before meshing, the CAD entities created needs to e synchronised with the GMsh model
gmsh.model.geo.synchronize()


# Set (2D) Physical groups on the geometry faces and name them so boundary conditions can be defined on them
pSurfaces = [gmsh.model.addPhysicalGroup(2, [surfaces[0]])]
pSurfaces.append(gmsh.model.addPhysicalGroup(2, [surfaces[1]]))
pSurfaces.append(gmsh.model.addPhysicalGroup(2, surfaces[2:]))
gmsh.model.setPhysicalName(2, pSurfaces[0], "Inlet")
gmsh.model.setPhysicalName(2, pSurfaces[1], "Outlet")
gmsh.model.setPhysicalName(2, pSurfaces[2], "Walls")

# Set 3D Physical group for domain volume
pVol = gmsh.model.addPhysicalGroup(3, [vol])
gmsh.model.setPhysicalName(3, pVol, "Domain")

# Generate 3D mesh
gmsh.model.mesh.generate(3)

# Write mesh to file
gmsh.write(outputFile)

# Visualise model in gui
if '-showmesh' in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()




