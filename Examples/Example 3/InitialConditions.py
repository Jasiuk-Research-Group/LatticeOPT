# Inputs:
RESTART = False
Start_on_iter = 0 # If RESTART = False , this is ignored

# Unit cell size:
UC_Lx = 100.
UC_Ly = 100.
UC_H = 8.5
dz = 0.5

# Grid size
Nx = 14 # Number of cells in X direction
Ny = 14 # Number of cells in Y direction
Nz = int(round(UC_H/dz))
Number_of_elements_per_edge = 5
DIM = [ UC_Lx , UC_Ly , UC_H ]
Ns = [ Nx , Ny , Nz , Number_of_elements_per_edge ]

# Exclude external boundary edges during optimization?
EXCLUDE_BOUNDARY = True
N_connect_layers = 1 # Specify how many layers should be denoted as connection layers
if EXCLUDE_BOUNDARY:
	# Specify boundary thickness:
	Boundary_thickness = 0.25

# Periodic lattice design?
PERIODIC = False
N_cell_x = 1
N_cell_y = 1

# Whether to tie the lattice to the top and bottom face plates or not
WRITE_TIE = True

# Optimization settings
filter_radius = 10. # Radius of the mesh independence filter
max_top_opt_iter = 25 # Max number of iteration for topology optimization
max_inner_solver_iter = 500 # Max number of iteration for thickness update

# Shell thickness parameters
ini_shell_thickness = 0.25
max_shell_thickness = 1.25
threshold_thickness = 0.05
max_dt = 0.05 # Max thickness change per iteration
# Mass of final design / mass of initial design
target_ratio_of_initial_design = 1.

# Choose to update thickness
Objective_type = 1
# If 1: Thickness update is done by min std( E / t )
# If 2: Thickness update is done by min std( t / V - E_e / E_tot )


base_file_name = 'case3_obj' + str(Objective_type) + '_max' + str(max_top_opt_iter)

# Abaqus settings
ncpu = 8 # Number of cpus to use when running sim