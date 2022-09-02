import numpy as np

# Inputs:
RESTART = False
Start_on_iter = 0 # If RESTART = False , this is ignored

# Unit cell size:
UC_Lx = 200.
UC_Ly = 30.
UC_H = 60.
dz = 2.

# Grid size
Nx = 20 # Number of cells in X direction
Ny = 3 # Number of cells in Y direction
Nz = int(round(UC_H/dz))
Number_of_elements_per_edge = 5
DIM = [ UC_Lx , UC_Ly , UC_H ]
Ns = [ Nx , Ny , Nz , Number_of_elements_per_edge ]

# Exclude external boundary edges during optimization?
EXCLUDE_BOUNDARY = True
N_connect_layers = 0 # Specify how many layers should be denoted as connection layers
if EXCLUDE_BOUNDARY:
	# Specify boundary thickness:
	Boundary_thickness = 0.75

# If EXCLUDE_BOUNDARY == True, please specify how the boundaries should be marked
def GetExclusionSet( Ns ):
	Nx , Ny , Nz , N_ele_per_edge = Ns
	N_edges = ( Nx + 1 ) * Ny + Nx * ( Ny + 1 )
	ExcSet = np.zeros( N_edges , dtype=bool )

	# # All bottom horizontal edges
	# ExcSet[ : Nx ] = True

	# All top horizontal edges
	ExcSet[ Nx * Ny : Nx * ( Ny + 1 ) ] = True

	# All left and right vertical edges
	for i in range( Ny ):
		offset = Nx * ( Ny + 1 )
		# ExcSet[ offset + i * ( Nx + 1 ) ] = True
		ExcSet[ offset + (i+1) * ( Nx + 1 ) - 1 ] = True
	return ExcSet



# Periodic lattice design?
PERIODIC = False
N_cell_x = 1
N_cell_y = 1

# Whether to tie the lattice to the top and bottom face plates or not
WRITE_TIE = False

# Optimization settings
filter_radius = 40. # Radius of the mesh independence filter
max_top_opt_iter = 25 # Max number of iteration for topology optimization
max_inner_solver_iter = 4000 # Max number of iteration for thickness update

# Shell thickness parameters
ini_shell_thickness = 0.75
max_shell_thickness = 2.
threshold_thickness = 0.4
max_dt = 0.08 # Max thickness change per iteration
# Mass of final design / mass of initial design
target_ratio_of_initial_design = 1.

# Choose to update thickness
Objective_type = 2
# If 1: Thickness update is done by min std( E / t )
# If 2: Thickness update is done by min std( t / V - E_e / E_tot )


base_file_name = 'case2_obj' + str(Objective_type) + '_max' + str(max_top_opt_iter)

# Abaqus settings
ncpu = 8 # Number of cpus to use when running sim