#!/usr/bin/env python
import numpy

def explTimeStep(T, material, dt = "max", f = 0.0):
    """
    Uses an explicit finite difference scheme with fixed
    boundary values to solve the diffusion equation.
    @param T: Object of grid class
    @param material: Object of material class
    @param dt: Time step length, [s]. If "max", the function will 
    determine the max. stable value, and return it.
    @param f: Local heat production, [W]. Can be a scalar or ndarray
    in the same dimensions as T.
    Returns the used time step length, [s] if dt is set to "max".
    """

    # Cell size
    dx2 = T.resolution**2.0
    dy2 = T.resolution**2.0

    # Material density, [kg/m^3]
    rho = 1.0

    # Material heat capacity, [J/(kg*K)]
    #c = 1.0
    c = material.c

    # Material thermal conductivity, [W/(m*K)]
    #k = 2.0
    k = material.k

    # Material diffusivity coefficient
    D = k / (rho * c)

    # Local heat sorce, [W]
    #f = 0.0

    # Source term
    F = f / (rho * c)

    # Set timestep to largest value possible if specified
    dt_max = dx2 * dy2 / (2.0 * D * (dx2 + dy2))
    returndt = False
    if (dt == "max"):
        dt = dt_max
        returndt = True
    elif (dt > dt_max):
        raise ValueError("Time step length (dt = {} s) is larger" \
                + " than the max. stable length (dt_max = {} s)"\
                .format(dt, dt_max))

    # Explicit 2D finite difference scheme
    T.grid[1:-1, 1:-1] += F * dt \
            + D * dt * ( \
                ( T.grid[2:, 1:-1] \
                  - 2.0 * T.grid[1:-1, 1:-1] \
                  + T.grid[:-2, 1:-1] )/dx2 \
              + ( T.grid[1:-1, 2:] \
                  - 2.0 * T.grid[1:-1, 1:-1] \
                  + T.grid[1:-1, :-2] ) /dy2 )

    if (returndt == True):
        return dt


