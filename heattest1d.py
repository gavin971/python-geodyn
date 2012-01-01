#!/usr/bin/env python

# 1D steady-state heat equation
# There are no thermal transients and no advection. k is constant.
# Solution to eq. 19
# Compares it to an analytical solution

# Import modules
import numpy
import matplotlib.pyplot

# Init grids
N = 50           # no of z-cells
z_max = 20.0     # max depth of z
dz = z_max/N     # cell width
z = numpy.arange(N)*dz # cell depths [m]
T = numpy.empty(N)     # cell temperatures [K]

# Physical parameters
k = 1.5    # thermal conductivity [W/(kg K)]
rho = 2000 # material density [kg/m^3]
c_p = 1000 # material specific heat capacity [J/(kg K)]
A = numpy.ones(N)*1e-6     # cell heat production values [W/m^3], try to change this
kappa = k/(rho * c_p)      # thermal diffusivity [m^2/s]

# Temporal parameters
dt = 24.0*60.0*60.0  # time step length [s]
t = 0.0              # current time
t_end = 10.0*365.0*24.0*60.0*60.0 # end time [s]
dt_plot = 7.0*24.0*60.0*60.0 # time between plots [s]

# Boundary conditions
T0mean = 273.15 - 3.0 # mean temperature [K]
T_s = T0mean    # surface temperature [K]
deltaT = 12.0   # amplitude of surface temperature fluctuations [K]
q_b = 0.065     # heat flux from below [W/m^2]

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
T = numpy.dot(numpy.linalg.inv(M), p)    # numpy does element-wise multiplication with * operator
# The steady-state solution is used as the initial values to the transient solver

# Find the new coefficient matrix
M = numpy.diag(numpy.ones(N) + 2.0*dt * kappa / (dz**2), 0)       # main diagonal
M = M + numpy.diag(numpy.ones(N-1)*-1.0 * dt*kappa / (dz**2), 1)  # upper diagonal
M = M + numpy.diag(numpy.ones(N-1)*-1.0 * dt*kappa / (dz**2), -1) # lower diagonal
# enforce BCs
M[0,0] = 1; M[0,1] = 0          # const. val
M[-1,-1] = -(rho*c_p/dt + 0.5*(2*k)/(dz**2))
M[-1,-2] = 0.5*(2*k)/(dz**2)

# Create plot
fig = matplotlib.pyplot.figure(1)
ax = fig.add_subplot(111)
ax.grid(True)

# Iterate through time
t_plot = 0
for i in range(int((t_end-t)/dt)):

    # Find surface temperature at time t
    T[0] = deltaT * numpy.sin(2.0 * numpy.pi/365.0 * t/(24.0*60.0*60.0)) + T0mean

    # Construct the vector of known terms
    p = T + dt * A / (rho*c_p)
    p[-1] = -rho * c_p / dt * T[-1] - A[-1] + q_b/dz

    # Find the solution to t+dt
    T_new = numpy.dot(numpy.linalg.inv(M), p)

    # Increment current time, save new values as current
    t += dt
    T = T_new
    t_plot += dt

    if (t_plot >= dt_plot):
        # plot z and T
        ax.plot(T-273.15, z, '-', label='w=' + str(int(t/dt_plot)))
        t_plot = 0

ax.legend()
ax.set_ylim(ax.get_ylim()[::-1]) # flip y axis
ax.set_title('Transient temperature')
#ax.set_xlabel('Temperature [K]')
ax.set_xlabel('Temperature [C]')
ax.set_ylabel('Depth [m]')
fig.show()

