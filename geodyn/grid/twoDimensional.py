#!/usr/bin/env python

import numpy
import matplotlib.pyplot

class uniformCubic:
    def __init__(self, world, resolution = 1e3, title = "Unnamed grid"):
        """
        Initialize homogeneous, cubic grid, fit to world dimensions.
        @param resolution: Cell length relative to world
        """

        self.grid = numpy.zeros([int(world.x/resolution), \
                                 int(world.y/resolution)])

        self.resolution = resolution
        self.title = title



    def addGauss(self, mu = 0.0, sigma = 1.0):
        """
        Add random numbers from the standard normal distribution to grid.
        @param mu: Mean value
        @param sigma: Variance
        """
        self.grid += numpy.sqrt(sigma) \
                     * numpy.random.randn(self.grid.shape[0], self.grid.shape[1]) \
                     + mu



    def coordinates(self, ix, iy):
        """
        Return the world coordinates of grid point
        @param ix: Grid point in x dimension
        @param iy: Grid point in y dimension
        """
        return numpy.array( \
                [ix*self.resolution, \
                 iy*self.resolution])


    def visualize(self, save = False, time = 0):
        """
        Plot grid using Matplotlib's imshow function.
        """
        fig = matplotlib.pyplot.figure(1)
        matplotlib.pyplot.imshow(self.grid)
        matplotlib.pyplot.colorbar()
        matplotlib.pyplot.title(self.title + ", t = {} s".format(time))
        if (save == True):
            fig.savefig(self.title + "_t={:10}s.png".format(time))
            fig.clf()
        else :
            fig.show()


    def gradient(self):
        """
        Returns the grid gradient, computed using central differences
        """
        return self.grid.gradient()


