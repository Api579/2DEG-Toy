import numpy
import matplotlib.pyplot as pyplot



#+++++++++++++++++++++++++++++++++++++++++++
#GaAs parameters and physical constants
#+++++++++++++++++++++++++++++++++++++++++++
#Usually 2DEG forms within GaAs layer when AlGaAs is added
#Electron goes where energy is lower
#GaAs conduction band < ALGaAs conduction band
#Gallium Arsenide (GaAs), Aluminium Gallium Arsenide

# charge of electron, (coulombs)
electron_charge = 1.602176634e-19
# mass of a free electron, (kilograms)
free_electron_mass = 9.1093837e-31
# reduced Planck constant, (joule per second)
reduced_planck_constant = 1.054571817e-34


# effective mass in GaAs, mass in GaAs
effective_mass = 0.067 * free_electron_mass
# how easily electrons move in GaAs 300K, in square metres per volt seconds,
electron_mobility = 0.8
# num of electrons per square metre of the 2DEG typically
electron_density = 3e15
# The voltage applied on top of left edge and right edge of AlGaAs/GaAs sandwhich
# realistic value for experiments
applied_voltage = 1e-3


#+++++++++++++++++++++++++++++++++++++++++++
#Initialisation of Grid and Grid points
#+++++++++++++++++++++++++++++++++++++++++++

grid_points_side = 10

grid_spacing = 10e-9


x_positions = numpy.arange(grid_points_side) * grid_spacing
y_positions = numpy.arange(grid_points_side) * grid_spacing

#creates grid of points of x and y dimension,
#and coordinates for each point in grid
x_grid, y_grid = numpy.meshgrid(x_positions, y_positions)

#+++++++++++++++++++++++++++++++++++++++++++
#Electron Initialisation
#+++++++++++++++++++++++++++++++++++++++++++

number_of_electrons = 100
left_edge = x_positions[0]
right_edge = x_positions[-1]
bottom_edge = y_positions[0]
top_edge = y_positions[-1]

grid_length = right_edge-left_edge

random_number = numpy.random.default_rng()

electron_x_positions = numpy.full(number_of_electrons, left_edge)
electron_y_positions = random_number.uniform(bottom_edge, top_edge, number_of_electrons)

# femto second scale, average scattering is 0.3 pico so need smaller to "see" whats happening
# but not too small                              
time_step = 1e-15
number_of_steps = 2000

#init electron count leaving grid on right side
electrons_out_right = 0

#+++++++++++++++++++++++++++++++++++++++++++
#Calculated values
#++++++++++++++++++++++++++++++++++++++++++

# E = V/L
electric_field = applied_voltage / grid_length
                               
# newton's 2nd law and f=qe, a= (qe)/m
push_from_voltage = electron_charge * electric_field / effective_mass
                               
# drude model defination, tau = mobility * m/q (seconds)
time_between_collisions = electron_mobility * effective_mass / electron_charge

# for each tinme step, P=1-e^(-"small change in t"/tau), poisson process
scattering_probability = 1.0 - numpy.exp(-time_step / time_between_collisions)

# fermi velocity
electron_speed = reduced_planck_constant * numpy.sqrt(2.0 * numpy.pi * electron_density) / effective_mass

# Randomise direction for each electron, -num.py/2.0 is -90 degrees
rand_angles = random_number.uniform(-numpy.pi / 2.0, numpy.pi / 2.0, number_of_electrons)

# electron velocities x & y direction
electron_x_v = electron_speed * numpy.cos(rand_angles)
electron_y_v = electron_speed * numpy.sin(rand_angles)

#+++++++++++++++++++++++++++++++++++++++++++
# Poisson solver
#+++++++++++++++++++++++++++++++++++++++++++
# solve for e field from electron postion
# solve for potential
# average with grid points, nearest neighbours

#  (farads per metre)
vacuum_permittivity = 8.8541878128e-12
gaas_permittivity = 12.9 * vacuum_permittivity

# Poisson needs charge per volume, so sheet charge (per area) needs to be divided by 2DEG thickness
two_deg_thickness = 10e-9

# real number of electrons in grid  = density x area
real_electrons_in_grid = electron_density * grid_length * grid_length

#  how many real electrons each simulated electron represent
electron_weight = real_electrons_in_grid / number_of_electrons

# area each grid point responsible for; edge points cover area of half a square
edge_weight = numpy.ones(grid_points_side)
edge_weight[0] = 0.5
edge_weight[-1] = 0.5

# shape (rows = y, columns = x), grid point area responsible, same as x_grid and y_grid from meshgrid
grid_point_area = grid_spacing * grid_spacing * numpy.outer(edge_weight, edge_weight)

# electric potential at every grid point, (volts), from 0 to applied_voltage, (q/l) * V, uniform e-field in neg e-direction
potential = (numpy.copy(x_grid) / grid_length )* applied_voltage

# field at every grid point, (volts per metre), init, empty arrays
grid_field_x = numpy.zeros((grid_points_side, grid_points_side))
grid_field_y = numpy.zeros((grid_points_side, grid_points_side))

# solve Poisson per these many steps, charge small change in 10 femto second
steps_between_poisson_solves = 10
# Jacobi passes per solve; each iteration closer to solved potential
jacobi_iterations = 200

#+++++++++++++++++++++++++++++++++++++++++++
#Plotting initial
#+++++++++++++++++++++++++++++++++++++++++++


# Draw a dot at each grid point
pyplot.scatter(x_grid/1e-9, y_grid/1e-9, color="lightgrey", label="grid points")

pyplot.scatter(electron_x_positions/1e-9, electron_y_positions/1e-9, color="blue", s=10, label="electrons")




#+++++++++++++++++++++++++++++++++++++++++++
#Main loop Monte Carlo
#+++++++++++++++++++++++++++++++++++++++++++

for step in range(number_of_steps):


    # clip boundaries, round the num of distance from left edge, finds nearest grid point in column
    # relative pos then grid index kept within boundary of grid
    
    nearest_column_x_point = numpy.clip(numpy.round((electron_x_positions - left_edge) / grid_spacing).astype(int), 0, grid_points_side - 1)
  
    nearest_row_y_point = numpy.clip(numpy.round((electron_y_positions - bottom_edge) / grid_spacing).astype(int), 0, grid_points_side - 1)

    #  charge only calculates on steps where Poisson solved
    if step % steps_between_poisson_solves == 0:
      
        # init electron count array of zeros
        electron_count = numpy.zeros((grid_points_side, grid_points_side))

        # update electron_count 2D array with electrons nearest to each grid point on (x,y), add 1
        numpy.add.at(electron_count, (nearest_row_y_point, nearest_column_x_point), 1)

      
        # find real electrons per square metre at each grid point, using weight found before
        electron_sheet_density = electron_count * electron_weight / grid_point_area
      
        # net sheet charge density eqn (= electron charge *(e_density - e_sheet_density))
        sheet_charge = electron_charge * (electron_density - electron_sheet_density)
      
        # charge per volume, for Poisson's equation
        volume_charge = sheet_charge / two_deg_thickness

        # solve Poisson's equation using Jacobi iteration
        for iteration in range(jacobi_iterations):
          
            # inside each point = average of 4 neighbours + charge term
            # eqn used: phi = 1/4 (left + right + up + down + spacing^2 x charge / permittivity)
            # finite difference approximation of Poisson's equation
            potential[1:-1, 1:-1] = 0.25 * (potential[1:-1, :-2] + potential[1:-1, 2:] + potential[:-2, 1:-1] + potential[2:, 1:-1] 
                                            + grid_spacing * grid_spacing * volume_charge[1:-1, 1:-1] / gaas_permittivity)
          
            # left edge (source) at 0 V
            potential[:, 0] = 0.0
          
            # right edge (drain) fixed applied voltage; electrons pulled to higher potential, so left to right
            potential[:, -1] = applied_voltage
          
            # Added bottom wall, no field crosses it, so bottom row = row just inside
            potential[0, 1:-1] = potential[1, 1:-1]
          
            # Added top wall, no e-field crossing out, top row efield =  efleid row just inside
            potential[-1, 1:-1] = potential[-2, 1:-1]

        # ADDED: ---------- Step 3: electric field from the potential, E = -slope of potential ----------
        # ADDED: numpy.gradient gives the slope along rows (y) first, then along columns (x)
        potential_slope_y, potential_slope_x = numpy.gradient(potential, grid_spacing)
        grid_field_x = -potential_slope_x
        grid_field_y = -potential_slope_y

    # ADDED: ---------- Step 4: push each electron with the field at its grid point ----------
  
    # ADDED: field felt by each electron
    electron_field_x = grid_field_x[nearest_row, nearest_column]
    electron_field_y = grid_field_y[nearest_row, nearest_column]
    # ADDED: a = -eE/m*, minus sign because the electron's charge is negative
    electron_x_v = electron_x_v - electron_charge * electron_field_x / effective_mass * time_step
    # ADDED: the field now also has a y part, so electrons get pushed up/down too
    electron_y_v = electron_y_v - electron_charge * electron_field_y / effective_mass * time_step



    electron_x_positions = electron_x_positions + electron_x_v * time_step
    electron_y_positions = electron_y_positions + electron_y_v * time_step

    # electrons going over the top of the 2D grid
    over_top = electron_y_positions > top_edge
    # mirror by the distance overshot, update positions
    electron_y_positions = numpy.where(over_top, 2.0 * top_edge - electron_y_positions, electron_y_positions)
    # flip volcity
    electron_y_velocities = numpy.where(over_top, -electron_y_v, electron_y_v)

    # same for going beow
    under_bottom = electron_y_positions < bottom_edge
    electron_y_positions = numpy.where(under_bottom, 2.0 * bottom_edge - electron_y_positions, electron_y_positions)
    electron_y_velocities = numpy.where(under_bottom, -electron_y_v, electron_y_v)

    # electrons going out the right end, has "left" the grid or 2DEG
    out_right = electron_x_positions > right_edge
    # Count the electron that has "left" the grid
    electrons_out_right = out_right + numpy.sum(out_right)
    
    # electrons that scattered backwards past left edge
    out_left = electron_x_positions < left_edge
    
    # Every electron that left through either left or right side
    outside_grid = out_right | out_left
    
    # number of electrons to replace
    number_to_replace = numpy.sum(outside_grid)

    # replacement electron enters from left, same as initialisation
    electron_x_positions[outside_grid] = left_edge
    
    # random y-axis coordinate
    electron_y_positions[outside_grid] = random_number.uniform(bottom_edge, top_edge, number_to_replace)
    
    # random direction for replacement
    replace_angles = random_number.uniform(-numpy.pi / 2.0, numpy.pi / 2.0, number_to_replace)

    electron_x_v[outside_grid] = electron_speed * numpy.cos(replace_angles)
    electron_y_v[outside_grid] = electron_speed * numpy.sin(replace_angles)

    # pick a random numbers, for probability 
    random_numbers = random_number.uniform(0.0, 1.0, number_of_electrons)
    # check if its below the scattering probability
    scattered = random_numbers < scattering_probability
    
   
    scattered_angles = random_number.uniform(0.0, 2.0 * numpy.pi, number_of_electrons)
  
    #speed unchanged only direction changes
    electron_x_v = numpy.where(scattered, electron_speed * numpy.cos(scattered_angles), electron_x_v)
    electron_y_v = numpy.where(scattered, electron_speed * numpy.sin(scattered_angles), electron_y_v)

#+++++++++++++++++++++++++++++++++++++++++++
#Plotting final
#+++++++++++++++++++++++++++++++++++++++++++

# Every 50 steps, draw electron trail
if step % 50 == 0:
    pyplot.scatter(electron_x_positions /10e-9, electron_y_positions /10e-9, color="blue", s=2, alpha=0.3)

# final positions of the electrons in red
pyplot.scatter(electron_x_positions /10e-9, electron_y_positions /10e-9, color="red", s=10, label="end")


# Make one unit in x equal one unit in y
pyplot.gca().set_aspect("equal")

# Label the horizontal and vertical axes
pyplot.xlabel("x (nm)")
pyplot.ylabel("y (nm)")

# Save the plot as an image file, trimmed empty space
pyplot.savefig("electron.png", bbox_inches="tight")
