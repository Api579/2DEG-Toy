import numpy
import matplotlib.pyplot as pyplot


X_grid_points_side = 10
Y_grid_points_side = 10
X_grid_spacing = 1.0
Y_grid_spacing = 1.0

x_positions = numpy.arange(X_grid_points_side) * X_grid_spacing
y_positions = numpy.arange(Y_grid_points_side) * Y_grid_spacing

#creates grid of points of x and y dimension,
#and coordinates for each point in grid
x_grid, y_grid = numpy.meshgrid(x_positions, y_positions)

# Draw a dot at each grid point
pyplot.scatter(x_grid, y_grid)
# Make one unit in x equal one unit in y
pyplot.gca().set_aspect("equal")

# Label the horizontal and vertical axes
pyplot.xlabel("x")
pyplot.ylabel("y")

# Save the plot as an image file
pyplot.savefig("electron.png")
