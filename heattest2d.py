#!/usr/bin/env python

# Import modules
import geodyn

# Create new world
world = geodyn.world()

# Initialize topographic grid
#world.topography = geodyn.grid.twoDimensional.uniformCubic(world, title = "Topography")

# Add random values to topography
#world.topography.addGauss(mu = 0.0, sigma = 5.0)

# Add a center layer of ice
world.heat = geodyn.grid.twoDimensional.uniformCubic(world, title = "Temperature")
world.heat.grid[8:13, 8:13] = 20.0
world.materials.layer1 = world.material()

# Plot topography
#world.topography.visualize()

world.heat.visualize(save = True, time = world.time)

for i in range(10):
    world.time += geodyn.surfaceProcesses.heat2D.explTimeStep(world.heat, world.materials.layer1)

    world.heat.visualize(save = True, time = world.time)


