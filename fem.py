#!/usr/bin/env python
import numpy
import scipy.spatial
import matplotlib.pyplot
import matplotlib.collections

class mesh:
    """ Finite Element Method mesh class """

    def __init__(meshobj, k = 1.0, A = 0.0):
        """ Class initializer
        @param k: Conductivity of each element (float)
        @param A: Source term for each element (float)
        """
        meshobj.k = k
        meshobj.A = A
        meshobj.elements = []


    def randomRect(meshobj, N = 10, x_min = 0.0, x_max = 1.0, y_min = 0.0, y_max = 1.0):
        """ Create N random mesh nodes in the rectangular area delimited by x and y limits
        @param N: Number of node points (int)
        @param x_min: Lower boundary of node positions along x (float)
        @param x_max: Upper boundary of node positions along x (float)
        @param y_min: Lower boundary of node positions along y (float)
        @param y_max: Upper boundary of node positions along y (float)
        """

        # Store parameters in object
        meshobj.N = N
        meshobj.x_min = x_min
        meshobj.x_max = x_max
        meshobj.y_min = y_min
        meshobj.y_max = y_max

        # Uniformly distributed nodes
        nodes = numpy.random.rand(meshobj.N, 2)

        # Create grid
        meshobj.delaunay(nodes)


    def delaunay(meshobj, nodes):
        """ Perform a Delaunay triangulation with the nodes as triangle corners """

        # Create Delaunay triangles with the points as centres (using the Qhull lib)
        meshobj.tri = scipy.spatial.Delaunay(nodes)

        # Number of elements (triangles) in the mesh
        meshobj.N_e = meshobj.tri.nsimplex

        # Construct the node coordinate matrix
        meshobj.COORD = meshobj.tri.points

        # Construct the topology matrix
        meshobj.TOPO = meshobj.tri.vertices

        # Create element objects
        for i in numpy.arange(meshobj.N_e):

            # Add element to mesh list of elements
            meshobj.elements.append(mesh.element(i, meshobj.TOPO[i,:], meshobj.COORD[meshobj.TOPO[i,:]], meshobj.k, meshobj.A))


    def findKf(meshobj):
        """ Construct the global stiffness matrix K and load vector f 
        from volumetric integration of the elements """
        meshobj.K = numpy.zeros((meshobj.N, meshobj.N))
        meshobj.f = numpy.zeros((meshobj.N, 1))
        for element in meshobj.elements:
            element.loadandstiffness()
            nodes = meshobj.TOPO[element.index,:]
            #print(nodes)
            #meshobj.K[meshobj.TOPO[element.index,:], meshobj.TOPO[element.index,:]] += element.K_e
            meshobj.f[nodes] += element.f_e
            i = 0
            for j in nodes:
                for k in nodes:
                    meshobj.K[j,k] += element.K_e.reshape(1,element.K_e.size)[0][i]
                    i += 1


    def ebc(meshobj, inodes, val):
        """ Enforce essential (Dirichlet) boundary condition on nodes with index i
        @inodes: Node indexes to enforce the essential boundary condition on (int array)
        @val: Value the nodes are fixed against (float)
        Using algorithm displayed in fig. 5, p. 12.
        """
        meshobj.f -= val * meshobj.K[inodes,:]
        meshobj.f[inodes] = val
        meshobj.K[inodes,:] = 0.0
        meshobj.K[:,inodes] = 0.0
        meshobj.K[inodes,inodes] = 1.0


    def plot(meshobj):
        """ Plot the mesh nodes and element borders """

        edge_points = []
        edges = set()
         
        def add_edge(i, j):
            """Add a line between the i-th and j-th points, if not in the list already"""
            if (i, j) in edges or (j, i) in edges:
                # already added
                return
            edges.add( (i, j) )
            edge_points.append(meshobj.COORD[ [i, j] ])
         
        # loop over triangles:
        # ia, ib, ic = indices of corner points of the triangle
        for ia, ib, ic in meshobj.TOPO:
            add_edge(ia, ib)
            add_edge(ib, ic)
            add_edge(ic, ia)

        lines = matplotlib.collections.LineCollection(edge_points)
        matplotlib.pyplot.hold(1)
        matplotlib.pyplot.plot(meshobj.COORD[:,0], meshobj.COORD[:,1], 'ko')
        matplotlib.pyplot.gca().add_collection(lines)
        matplotlib.pyplot.axis('equal')
        matplotlib.pyplot.xlim(meshobj.x_min, meshobj.x_max)
        matplotlib.pyplot.ylim(meshobj.y_min, meshobj.y_max)
        matplotlib.pyplot.show()


    class element:
        """ Finite Element Method element class """

        def __init__(elementobj, index, nodes, coord, k, A):
            """ Element object initializer 
            @param index: Index of this element (int)
            @param nodes: List of node indexes in this element (numpy.array)
            """
            elementobj.index = index
            elementobj.nodes = nodes
            elementobj.coord = coord
            elementobj.k = k
            elementobj.A = A


        def gradlt3(elementobj):
            """ Sets the local gradient matrix (eq. 17) """
            elementobj.dphi_l = numpy.array([[-1.0, 1.0, 0.0],[-1.0, 0.0, 1.0]])

        def jacobian(elementobj):
            """ Calculates the Jacobian matrix (eq. 26) and it's determinant """
            elementobj.J = numpy.dot(elementobj.dphi_l, elementobj.coord)
            elementobj.detJ = numpy.linalg.det(elementobj.J)

        def gradg(elementobj):
            """ Sets the global gradient matrix (eq. 22) """
            elementobj.jacobian()
            elementobj.dphi_g = numpy.dot(numpy.linalg.inv(elementobj.J), elementobj.dphi_l)

        def gausst3(elementobj, N_ip):
            """ Returns Gauss data for 3 node triangle (fig. 3, p. 9)
            @param N_ip: Number of integration points [int]
            Returns r_i as 1st col, s_i as 2nd col, w_i as 3rd col,
            and one row per integration point.
            """
            if (N_ip == 1):
                return numpy.array([[1.0/3.0, 1.0/3.0, 1.0/2.0]])

            elif (N_ip == 3):
                return numpy.array([
                    [1.0/6.0, 1.0/6.0, 1.0/6.0],
                    [2.0/3.0, 1.0/6.0, 1.0/6.0],
                    [1.0/6.0, 2.0/3.0, 1.0/6.0]])

            elif (N_ip == 4):
                return numpy.array([
                    [1.0/3.0, 1.0/3.0, -9.0/32.0],
                    [3.0/5.0, 1.0/5.0, 25.0/96.0],
                    [1.0/5.0, 3.0/5.0, 25.0/96.0],
                    [1.0/5.0, 1.0/5.0, 25.0/96.0]])

            elif (N_ip == 7):
                return numpy.array([
                    [0.0, 0.0, 1.0/40.0],
                    [0.5, 0.0, 1.0/15.0],
                    [1.0, 0.0, 1.0/40.0],
                    [0.5, 0.5, 1.0/15.0],
                    [0.0, 1.0, 1.0/40.0],
                    [0.0, 0.5, 1.0/15.0],
                    [1.0/3.0, 1.0/3.0, 9.0/40.0]])

            else :
                raise Exception("gausst3 can only return data for 1, 3, 4, or 7 integration points")

        def loadandstiffness(elementobj, N_ip = 1):
            """ Calculate the element load vector
            @param N_ip: Number of integration points to evaluate in each element (int)
            """
            # Get quadrature data [r_i, s_i, w_i], one row per N_ip
            g_i = elementobj.gausst3(N_ip)

            # Find local derivatives of the shape function (stored as elementobj.dphi_l)
            elementobj.gradlt3()

            # Get interpolation function values
            phi = elementobj.shplint3(g_i[:,0], g_i[:,1])

            # Find global gradients of the Jacobian (stored as elementobj.dphi_g)
            elementobj.gradg()

            elementobj.K_e = numpy.zeros((3,3)) # Element stiffness matrix
            elementobj.f_e = numpy.zeros((3,1)) # Element load vector

            # Loop over integration points
            for i in numpy.arange(N_ip):
                
                elementobj.K_e += g_i[i,2] * elementobj.k * numpy.dot(elementobj.dphi_l.T, elementobj.dphi_l * elementobj.detJ)    # eq. 38, K^e_i
                elementobj.f_e += g_i[i,2] * phi * elementobj.A * elementobj.detJ  # eq. 39, f^e_i
            

        # Integration point level
        def shplint3(elementobj, r, s):
            """ Returns the linear shape function vector """
            return numpy.array([1.0 - r - s, r, s])



testgrid = mesh()
testgrid.randomRect(N=100)
testgrid.findKf()
#testgrid.plot()
