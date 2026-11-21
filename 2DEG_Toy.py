import numpy
import matplotlib
matplotlib.use("Agg") #remove window disp, saves only final png, for speed-up
import matplotlib.pyplot as pyplot

import scipy.sparse
import scipy.sparse.linalg

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
# reduced Planck constant, (joule*second)
reduced_planck_constant = 1.054571817e-34


# effective mass in GaAs, mass in GaAs
effective_mass = 0.067 * free_electron_mass
# how easily electrons move in GaAs 300K, in square metres per volt seconds,
electron_mobility = 0.8
# num of electrons per square metre of the 2DEG typically
electron_density = 3e15
# The voltage applied on top of left edge and right edge of AlGaAs/GaAs sandwhich
# value large enough to not be lost as noise (if too small, it does not contribute much to thermal speed)
applied_voltage = 10e-3


#+++++++++++++++++++++++++++++++++++++++++++
#Initialisation of Grid and Grid points
#+++++++++++++++++++++++++++++++++++++++++++

# grid too close repel, electrons "smooths out" ~8 nm, too close repel
# at least 2 grid points across it [2]
# 3 nm spacing gives ~2.6;
# 31 points makes 30 gaps, so 31 points keeps the device 90 nm long
grid_points_side = 31

grid_spacing = 3e-9


x_positions = numpy.arange(grid_points_side) * grid_spacing
y_positions = numpy.arange(grid_points_side) * grid_spacing

#creates grid of points of x and y dimension,
#and coordinates for each point in grid
x_grid, y_grid = numpy.meshgrid(x_positions, y_positions)

#+++++++++++++++++++++++++++++++++++++++++++
#Electron Initialisation
#+++++++++++++++++++++++++++++++++++++++++++

# more electrons so each carries less charge, and random jitter is smaller
# published runs use tens of thousands [1], [2]
number_of_electrons = 20000
left_edge = x_positions[0]
right_edge = x_positions[-1]
bottom_edge = y_positions[0]
top_edge = y_positions[-1]

grid_length = right_edge-left_edge

random_number = numpy.random.default_rng()

# start electrons spread randomly [1] if all at x = 0, first Poisson solve see no electrons affect calculations.
electron_x_positions = random_number.uniform(left_edge, right_edge, number_of_electrons)

electron_y_positions = random_number.uniform(bottom_edge, top_edge, number_of_electrons)

# femto second scale, average scattering is 0.3 pico so need smaller to "see" whats happening
# but not too small                              
time_step = 1e-15


# let simulation settle, then average over the long settled stretch [1]
# first 3 ps thrown away; needs > 5 collisions ~1.5 ps to settle, 
number_of_steps = 15000
settling_steps = 3000

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

# Boltzmann constant 
boltzmann_constant = 1.380649e-23
# room temp
temperature = 300.0

# thermal_speed sqrt(kT/m), using boltzman distribution for more realism
thermal_speed = numpy.sqrt(boltzmann_constant * temperature / effective_mass)


#flux-weighted Maxwellian (Rayleigh)
random_values = random_number.uniform(0.0, 1.0, number_of_electrons)

# electrons already in device: random directions, Gaussian for velocity [1], [3]
# Rayleigh is only for electrons entering through a contact [3] only for injection, faster ones cross into device more often
electron_x_v = thermal_speed * random_number.normal(0.0, 1.0, number_of_electrons)

electron_y_v = thermal_speed * random_number.normal(0.0, 1.0, number_of_electrons)

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


# electrons move back-forth fast (190 fs/cycle, plasma freq ~3.3e13 rad/s)
# update the voltage often enough to track this movement [2]
# plasma frequency*interval < 0.2; every 5 fs gives ~0.17 [2]
steps_between_poisson_solves = 5



# build Poisson matrix

# Field crosses left and right edges, but removed cus potential fixed (Dirichlet)
poisson_array_unknown_sizing = grid_points_side * (grid_points_side - 2)


# indexing the points of array, indexes of points to solve for poisson
def indexing(i,j):
    return i * (grid_points_side - 2) + (j - 1)



poisson_rows = [] #eqn
poisson_cols = [] #coordinate
poisson_values = [] #coefficient

# boundary contribution to the right-hand side (from left: 0, right:10mV voltages)
boundary_rhs = numpy.zeros(poisson_array_unknown_sizing)







# finite-difference Poisson for each point: 4 x (point) - left - right - up - down = spacing^2 x charge / permittivity
for i in range(grid_points_side):
    for j in range(1, grid_points_side - 1):
      
        k = indexing(i, j)
      
        # left and right neighbours
        # if it's a contact, known voltage goes to the right-hand side
        for jj in (j - 1, j + 1):
            if jj == 0:
                boundary_rhs[k] += 0.0
            elif jj == grid_points_side - 1:
                boundary_rhs[k] += applied_voltage
            else:
                poisson_rows.append(k)
                poisson_cols.append(indexing(i, jj))
                poisson_values.append(-1.0)
              
        # up and down neighbours;
        # if pass wall, use the mirror-image point instead
        # mirror, assumes voltage flat across the wall and no field through it
        # the mirror point is added twice, the matrix adds them up to -2
        # since mirrored now 2 above, in upper row 
      
        for ii in (i - 1, i + 1):
            if ii < 0:
                ii = 1
            elif ii > grid_points_side - 1:
                ii = grid_points_side - 2
              
            poisson_rows.append(k)
            poisson_cols.append(indexing(ii, j))
            poisson_values.append(-1.0)
    
        poisson_rows.append(k)
        poisson_cols.append(k)
        poisson_values.append(4.0)

####
matrix_rows = poisson_rows
matrix_cols = poisson_cols
matrix_values = poisson_values

#builds matrix with the recorded values
poisson_matrix = scipy.sparse.csc_matrix((matrix_values, (matrix_rows, matrix_cols)),shape=(poisson_array_unknown_sizing, poisson_array_unknown_sizing))
#factorise matrix to be later solved for voltage at every point inside grid in loop
poisson_solver = scipy.sparse.linalg.splu(poisson_matrix)


#+++++++++++++++++++++++++++++++++++++++++++
#Plotting initial
#+++++++++++++++++++++++++++++++++++++++++++


# Draw a dot at each grid point
pyplot.scatter(x_grid/1e-9, y_grid/1e-9, color="lightgrey", label="grid points")

pyplot.scatter(electron_x_positions/1e-9, electron_y_positions/1e-9, color="blue", s=10, label="electrons")

# only draw 200, not full electrons
plot_count = 200
pyplot.scatter(electron_x_positions[:plot_count]/1e-9, electron_y_positions[:plot_count]/1e-9, color="blue", s=10, label="electrons")


drift_velocity_sum = 0.0
samples = 0

# stores electron path to output once at end as saved png, speed-up
path_x = []
path_y = []



#+++++++++++++++++++++++++++++++++++++++++++
#Main loop 
#+++++++++++++++++++++++++++++++++++++++++++

for step in range(number_of_steps):




    # share individual electron's charge between 4 corners of grid-square [2]
    # closer corner has bigger share, this is smoother than all-to-nearest-point
    # position measured in grid-squares
  
    grid_x = (electron_x_positions - left_edge) / grid_spacing
    grid_y = (electron_y_positions - bottom_edge) / grid_spacing
  
    # electron coordinates in grid square, ans kept within the device, 2D-mesh grid
    grid_square_x = numpy.clip(numpy.floor(grid_x).astype(int), 0, grid_points_side - 2)
    grid_square_y = numpy.clip(numpy.floor(grid_y).astype(int), 0, grid_points_side - 2)
  
    # length across square for electron
    fractions_x = grid_x - grid_square_x
    fractions_y = grid_y - grid_square_y
  
    # share of each corner where 1 is 100%, 0.7 70% so on
    q00 = (1 - fractions_x) * (1 - fractions_y)
    q10 = fractions_x * (1 - fractions_y)
    q01 = (1 - fractions_x) * fractions_y
    q11 = fractions_x * fractions_y

    #  charge only calculates on steps where Poisson solved
    if step % steps_between_poisson_solves == 0:
      
        # init electron count array of zeros
        electron_count = numpy.zeros((grid_points_side, grid_points_side))

        # update electron_count 2D array with electrons nearest to each grid point on (x,y), add 1
        numpy.add.at(electron_count, (nearest_rows_y_point, nearest_columns_x_point), 1)

        # add based on corner "percentage"
        numpy.add.at(electron_count, (grid_square_y, grid_square_x), q00)
        numpy.add.at(electron_count, (grid_square_y, grid_square_x + 1), q10)
        numpy.add.at(electron_count, (grid_square_y + 1, grid_square_x), q01)
        numpy.add.at(electron_count, (grid_square_y + 1, grid_square_x + 1), q11)
      
        # find real electrons per square metre at each grid point, using weight found before
        electron_sheet_density = electron_count * electron_weight / grid_point_area
      
        # net sheet charge density eqn (= electron charge *(e_density - e_sheet_density))
        sheet_charge = electron_charge * (electron_density - electron_sheet_density)
      
        # charge per volume, for Poisson's equation
        volume_charge = sheet_charge / two_deg_thickness

        # matrix_built with coefficients * voltage inside solved for = right hand side
        # solve Poisson's equation directly (4 phi - neighbours = spacing^2 x charge / permittivity)

        # volume charge does not consider walls and flatten from 2D to 1D, to get RHS voltage, voltage from e and voltage from contacts
        right_hand_side = boundary_rhs + (grid_spacing * grid_spacing / gaas_permittivity) * volume_charge[:, 1:-1].ravel()
        
        # solve for volage and put back into the grid for other calculation
        potential[:, 1:-1] = poisson_solver.solve(right_hand_side).reshape(grid_points_side, grid_points_side - 2)


          
        # left edge (source) at 0 V
        potential[:, 0] = 0.0
          
        # right edge (drain) fixed applied voltage; electrons pulled to higher potential, so left to right
        potential[:, -1] = applied_voltage
          


        
        # calculate electric field–potential relationship
        potential_slope_y, potential_slope_x = numpy.gradient(potential, grid_spacing)
        grid_field_x = -potential_slope_x
        grid_field_y = -potential_slope_y
      
        grid_field_y[0, :] = 0.0
        grid_field_y[-1, :] = 0.0

    # push by field, update v and p of electron
  
    # field felt by each electron
    electron_field_x = (q00 * grid_field_x[grid_square_y, grid_square_x] + q10 * grid_field_x[grid_square_y, grid_square_x + 1]
                        + q01 * grid_field_x[grid_square_y + 1, grid_square_x] + q11 * grid_field_x[grid_square_y + 1, grid_square_x + 1])
    electron_field_y = (q00 * grid_field_y[grid_square_y, grid_square_x] + q10 * grid_field_y[grid_square_y, grid_square_x + 1] 
                        + q01 * grid_field_y[grid_square_y + 1, grid_square_x] + q11 * grid_field_y[grid_square_y + 1, grid_square_x + 1])
  
    # a = -eE/m*, F=ma, F= -qE
    electron_x_v = electron_x_v - electron_charge * electron_field_x / effective_mass * time_step
    electron_y_v = electron_y_v - electron_charge * electron_field_y / effective_mass * time_step



    electron_x_positions = electron_x_positions + electron_x_v * time_step
    electron_y_positions = electron_y_positions + electron_y_v * time_step

    # electrons going over the top of the 2D grid
    over_top = electron_y_positions > top_edge
    # mirror by the distance overshot, update positions
    electron_y_positions = numpy.where(over_top, 2.0 * top_edge - electron_y_positions, electron_y_positions)
    # flip volcity
    electron_y_v = numpy.where(over_top, -electron_y_v, electron_y_v)

    # same for going beow
    under_bottom = electron_y_positions < bottom_edge
    electron_y_positions = numpy.where(under_bottom, 2.0 * bottom_edge - electron_y_positions, electron_y_positions)
    electron_y_v = numpy.where(under_bottom, -electron_y_v, electron_y_v)

    # electrons going out the right end, has "left" the grid or 2DEG
    out_right = electron_x_positions > right_edge
    # Count the electron that has "left" the grid
    electrons_out_right = electrons_out_right + numpy.sum(out_right)
    
    # electrons that scattered backwards past left edge
    out_left = electron_x_positions < left_edge
    
    # Every electron that left through either left or right side
    outside_grid = out_right | out_left
    
    # number of electrons to replace
    number_to_replace = numpy.sum(outside_grid)



    left_injection_fraction = 0.5
   
    
    # True = enters from left, False = enters from right, electron injection either side
    enter_from_left = random_number.uniform(0.0, 1.0, number_to_replace) < left_injection_fraction

    # left edge or right edge, injected
    electron_x_positions[outside_grid] = numpy.where(enter_from_left, left_edge, right_edge)

    # random y-axis coordinate assigned
    electron_y_positions[outside_grid] = random_number.uniform(bottom_edge, top_edge, number_to_replace)

    # random new speed, boltzman maxwell dist
  
    random_values = random_number.uniform(0.0, 1.0, number_to_replace)
    new_speed = thermal_speed * numpy.sqrt(-2.0 * numpy.log(1.0 - random_values))

    # left contact moves right, right contact moves left (sign: -ve)
    electron_x_v[outside_grid] = numpy.where(enter_from_left, new_speed, -new_speed)

    # sidewards speed, random using boltzman dist
    electron_y_v[outside_grid] = thermal_speed * random_number.normal(0.0, 1.0, number_to_replace)







  
    # pick a random numbers, for probability 
    random_numbers = random_number.uniform(0.0, 1.0, number_of_electrons)
    # check if its below the scattering probability
    scattered = random_numbers < scattering_probability
    
   
    # phonon scattering, new rand velocity using maxwell-boltzman distributio 
    # speed, direction reset,
    # selectrons heated by the field relax to lattice temperature, relaxation
    electron_x_v = numpy.where(scattered, thermal_speed * random_number.normal(0.0, 1.0, number_of_electrons), electron_x_v)
    electron_y_v = numpy.where(scattered, thermal_speed * random_number.normal(0.0, 1.0, number_of_electrons), electron_y_v)

#+++++++++++++++++++++++++++++++++++++++++++
#Plotting final
#+++++++++++++++++++++++++++++++++++++++++++

    # save electrons' positions for the trail plot
    if step % 50 == 0:
        path_x.append(electron_x_positions[:plot_count].copy())
        path_y.append(electron_y_positions[:plot_count].copy())

    # once settled, add (average * speed) each step [1]
    if step >= settling_steps:
        drift_velocity_sum += electron_x_v.mean()
        samples += 1


# draw electron path once
pyplot.scatter(numpy.concatenate(path_x) / 1e-9, numpy.concatenate(path_y) / 1e-9, color="blue", s=2, alpha=0.3)

# final positions of the electrons in red
pyplot.scatter(electron_x_positions[:plot_count] /1e-9, electron_y_positions[:plot_count] /1e-9, color="red", s=10, label="end")


# Make one unit in x equal one unit in y
pyplot.gca().set_aspect("equal")

# Label the horizontal and vertical axes
pyplot.xlabel("x (nm)")
pyplot.ylabel("y (nm)")

# Save the plot as an image file, trimmed empty space
pyplot.savefig("electron.png", bbox_inches="tight")

drift_velocity = drift_velocity_sum / samples
print(f"drift velocity = {drift_velocity:.3g} m/s")
print(f"apparent mobility = {drift_velocity / (applied_voltage / grid_length):.3g} m^2/Vs (input {electron_mobility})")


