#!/usr/bin/env python

# 1D steady-state heat equation
# There are no thermal transients and no advection. k is constant.
# Solution to eq. 19
# Compares it to an analytical solution

# Import modules
import numpy
import matplotlib.pyplot

# Init grids
N = 30          # no of z-cells
z_max = 10.0     # max depth of z
dz = z_max/N     # cell width
z = numpy.arange(N)*dz # cell depths [m]
T = numpy.empty(N)     # cell temperatures [K]

# Physical parameters
k = 1.5     # thermal conductivity
A = numpy.ones(N)*1e-6     # cell heat production values [W/m^3], try to change this

# Boundary conditions
T_s = 273.15 + 3.0 # surface temperature
q_b = 0.065        # heat flux from below

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

# Calculate the analytical solution to check the result
#T_a = 

# plot z and T
fig = matplotlib.pyplot.figure(1)
ax = fig.add_subplot(111)
ax.grid(True)
ax.plot(z,T, label='Finite difference')
ax.legend()
ax.set_xlabel('Depth [m]')
ax.set_ylabel('Temperature [K]')
ax.set_title('Steady state temperature')
fig.show()



