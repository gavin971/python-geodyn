#!/usr/bin/env python

import numpy


def explTimeStep(H, world, material):
    """
    Updates the ice thickness grid using an explicit time step.
    @param H: Object of grid class, ice height
    @param world: Object of world class
    @param material: Object of material class
    Returns the used time step length, [s]
    """

    # Cell size
    dx2 = T.resolution**2.0
    dy2 = T.resolution**2.0

    secondsPerYr = 60 * 60 * 24 * 365.25

    # Ice rheology constant
    Gamma = 2.0 * material.A/secondsPerYr * (material.rho * world.g)**material.n / (material.n+2.0)

    # Find timestep length
    dt = 0.0002 * secondsPerYr



    return dt

