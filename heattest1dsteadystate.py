#!/usr/bin/env python

# 1D steady-state heat equation
# There are no thermal transients and no advection. k is constant.
# Solution to eq. 19
# Compares it to an analytical solution

# Import modules
import numpy
import matplotlib.pyplot

# Init grids
N = 30           # no of z-cells
z_max = 10.0     # max depth of z
dz = z_max/N     # cell width
z = numpy.arange(N)*dz # cell depths [m]
T = numpy.empty(N)     # cell temperatures [K]

# Initialize figure
fig = matplotlib.pyplot.figure(1)
ax = fig.add_subplot(111)
ax.grid(True)

# thermal conductivity
for k in [0.5, 1.5, 2.0]:
    # Physical parameters
    A = numpy.ones(N)*1e-6     # cell heat production values [W/m^3], try to change this

    # Boundary conditions
    T_s = 273.15 + 3.0 # surface temperature
    q_b = -0.065       # heat flux (neg. when upwards)

    # Construct diagonals to the coefficient matrix
    M = numpy.diag(numpy.ones(N)*-2, 0)     # main diagonal
    M = M + numpy.diag(numpy.ones(N-1), 1)  # upper diagonal
    M = M + numpy.diag(numpy.ones(N-1), -1) # lower diagonal
    # enforce BCs
    M[0,0] = 1; M[0,1] = 0          # const. val
    M[-1,-1] = 1; M[-1,-2] = -1     # const. flux

    # Construct the vector of known terms
    p = -A * dz**2 / k
    p[0] = T_s
    p[-1] = -q_b * dz/k

    # The temperature (T) is the vector of unknowns, and can be solved
    T = numpy.dot(numpy.linalg.inv(M),p)    # numpy does element-wise multiplication with * operator

    # Calculate the analytical solution for T to check the result
    # The analytical solution is equal to:
    #   import sympy
    #   T, z, A, k, c1, c2 = sympy.symbols('T z A k c1 c2')
    #   sympy.integrate(sympy.integrate(-A/k, z) + c1, z) + c2
    #   => -A*z**2/(2*k) + c1*z + c2
    # c2 is the temperature at z=0, i.e. T_s.
    # c1 is equal to the gradient at z=z_max, i.e. q_b/-k + A*z_max/k
    T_a = -A/(2*k) * z**2 + (-q_b/k + A*numpy.max(z)/k)*z + T_s

    # plot z and T
    ax.plot(T-273.15, z, '-', label='Finite difference, k=' + str(k))
    ax.plot(T_a-273.15, z, '+--', label='Analytical, k=' + str(k))

ax.legend()
ax.set_ylim(ax.get_ylim()[::-1]) # flip y axis
#ax.set_xlabel('Temperature [K]')
ax.set_xlabel('Temperature [C]')
ax.set_ylabel('Depth [m]')
ax.set_title('Steady state temperature')
fig.show()

