#!/usr/bin/env python
import numpy
import scipy.spatial
import matplotlib.pyplot
import matplotlib.collections

class mesh:
    """ Finite Element Method mesh class """

    def __init__(self):
        """ Class initializer """
        pass

    def randomRect(self, N = 10, x_min = 0.0, x_max = 1.0, y_min = 0.0, y_max = 1.0):
        """ Create N random mesh nodes in the rectangular area delimited by x and y limits
        @param N: Number of node points (int)
        @param x_min: Lower boundary of node positions along x (float)
        @param x_max: Upper boundary of node positions along x (float)
        @param y_min: Lower boundary of node positions along y (float)
        @param y_max: Upper boundary of node positions along y (float)
        """

        # Store parameters in object
        self.N = N
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max

        # Uniformly distributed nodes
        nodes = numpy.random.rand(self.N, 2)

        # Create Delaunay triangles with the points as centres (using the Qhull lib)
        self.delaunay(nodes)

    def delaunay(self, nodes):
        """ Perform a Delaunay triangulation with the nodes as triangle corners """

        # Create Delaunay triangles with the points as centres (using the Qhull lib)
        self.tri = scipy.spatial.Delaunay(nodes)

        # Number of elements (triangles) in the mesh
        self.N_e = self.tri.nsimplex

        # Construct the node coordinate matrix
        self.COORD = self.tri.points

        # Construct the topology matrix
        self.TOPO = self.tri.vertices


    def plot(self):
        """ Plot the mesh nodes and element borders """

        edge_points = []
        edges = set()
         
        def add_edge(i, j):
            """Add a line between the i-th and j-th points, if not in the list already"""
            if (i, j) in edges or (j, i) in edges:
                # already added
                return
            edges.add( (i, j) )
            edge_points.append(self.COORD[ [i, j] ])
         
        # loop over triangles:
        # ia, ib, ic = indices of corner points of the triangle
        for ia, ib, ic in self.TOPO:
            add_edge(ia, ib)
            add_edge(ib, ic)
            add_edge(ic, ia)

        lines = matplotlib.collections.LineCollection(edge_points)
        matplotlib.pyplot.hold(1)
        matplotlib.pyplot.plot(self.COORD[:,0], self.COORD[:,1], 'ko')
        matplotlib.pyplot.gca().add_collection(lines)
        matplotlib.pyplot.axis('equal')
        matplotlib.pyplot.xlim(self.x_min, self.x_max)
        matplotlib.pyplot.ylim(self.y_min, self.y_max)
        matplotlib.pyplot.show()



testgrid = mesh()
testgrid.randomRect()
testgrid.plot()
