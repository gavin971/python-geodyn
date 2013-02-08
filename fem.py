#!/usr/bin/env python
import numpy
import scipy.spatial

# Create some random points to create the grid from
N_e = 30 # number of elements
N_d= 2 # number of dimensions
points = numpy.random.rand(N_e, N_d) # uniform distributed points

# Create Delaunay triangles with the points inside (using the Qhull lib)
tri = scipy.spatial.Delaunay(points)

# Triangle corner (node) coordinates
p = tri.points[tri.vertices]

# Number of nodes in the mesh
N = tri.nsimplex

# Construct the COORD matrix
COORD = numpy.array([p[:,0,0], p)   # allocate the matrix

# Plot mesh
import matplotlib.pyplot as plt
plt.hold(1)
plt.plot(points[:,0], points[:,1], '.')
#plt.plot(cc[0], cc[1], '*')
#plt.gca().add_collection(lines)
plt.axis('equal')
plt.xlim(-0.1, 1.1)
plt.ylim(-0.1, 1.1)
#plt.show()
