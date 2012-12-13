#!/usr/bin/env python

class world:
    """
    World class containing grids (topography, etc.)
    """

    def __init__(self, x = 100e3, y = 100e3, time = 0.0, g = 9.81):
        """
        Class initializer
        @param x: World width, [m]
        @param y: World length, [m]
        @param time: Current time, [s]
        """

        self.x = x
        self.y = y
        self.time = time
        self.g = g

    def getSize(self):
        """
        Returns world size
        """

        return numpy.array([self.x, self.y], dtype=floattype)


    class material:
        """
        Class for physical material properties
        """

        def __init__(self, rho = 1.0, c = 1.0, k = 2.0, A = 1.0e-16, n = 3.0):
            """
            Initializes material with values of the physical properties.
            @param rho: Material denssity, [kg/m^3]
            @param c: Material heat capacity, [J/(kg*K)]
            @param k: Material thermal conductivity, [W/(m*K)]
            @param A: Material softness factor (Glen's flow law)
            @param n: Material exponent (Glen's flow law)
            """

            self.rho = rho
            self.c = c
            self.k = k
            self.A = A
            self.n = n

    class materials:
        """
        Material class, containing physical properties of the materials defined
        """
        pass


