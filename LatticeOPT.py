from LatticeOPT_lib import *
from InitialConditions import *
import cProfile

################################################################################################
if PERIODIC and not EXCLUDE_BOUNDARY:
	print('Currently, periodic option needs to exclude outer boundaries')
	exit()

if EXCLUDE_BOUNDARY:
	ExcSet = GetExclusionSet( Ns )
	ConnectSet = GetConnectionSet( Ns , N_connect_layers )
else:
	ExcSet = np.zeros( ( Nx + 1 ) * Ny + Nx * ( Ny + 1 ) , dtype = bool )
	ConnectSet = ExcSet.copy()

if PERIODIC:
	ComputePeriodicMtx( Nx , Ny , N_cell_x , N_cell_y , ExcSet )

# Get lattice wall lengths
wall_length = np.zeros_like( ExcSet )
offset = Nx * ( Ny + 1 )
wall_length[ :offset ] = ( UC_Lx / N_cell_x ) / Nx # Horizontal edges
wall_length[ offset: ] = ( UC_Ly / N_cell_y ) / Ny # Vertical edges

if not RESTART:
	# Initial conditions
	# Edge thickness array, uniform thickness
	thickness = np.ones( ( Nx + 1 ) * Ny + Nx * ( Ny + 1 ) ) * ini_shell_thickness
		
	# Get initial sum of thickness
	initial_sum = np.dot( np.multiply( thickness , wall_length ) , 1. - ExcSet ) # Only sum those are included in optimization
	thickness *= target_ratio_of_initial_design # Linearly scale the thickness to satisfy the mass constraint directly on step 1

	# Specify distribution
	thickness_dist = np.ones( [ len(thickness) , 1 ] ) * thickness[0]

	# Set exclusion set
	if EXCLUDE_BOUNDARY:
		thickness[ ExcSet ] = Boundary_thickness
		thickness_dist[ ExcSet , : ] = Boundary_thickness

	# Save initial design
	f = open('DesignIteration-1.npy' , 'wb')
	np.save( f , np.array([None,thickness,[thickness_dist]],dtype=object) )
	f.close()
	print('Saved initial material distribution!')

	f=open('OptLog_' + base_file_name + '.txt','w')
	f.close()
	Start_on_iter = 0
else:
	f = open('DesignIteration-1.npy' , 'rb')
	_,thickness,_ = np.load( f , allow_pickle=True )
	f.close()
Target_sum_A = np.dot( np.multiply( thickness , wall_length ) , 1. - ExcSet ) * target_ratio_of_initial_design



# Write out base inp
WriteBaseINP( base_file_name , DIM , Ns , WRITE_TIE )

# Write mesh filter
# ComputeMtx( DIM , Ns , filter_radius )


################################################################################################################
# Begin optimization!
for itr in range( Start_on_iter , Start_on_iter + max_top_opt_iter ):
	st = time.time()
	print('\n\n                                                         >>> Starting iteration ' + str(itr) + ' / ' + str(max_top_opt_iter) )

	# Write new inp with section assigment
	jn = SectionAssignment( base_file_name , itr , Ns )	

	# Run simulation and gather data from Abaqus
	f=open('currItr','wb')
	itr_dict = { 'jn' : jn ,
				 'ncpu' : ncpu ,
				 'Ns' : Ns ,
				 'max_inner_solver_iter' : max_inner_solver_iter ,
				 'm_new' : Target_sum_A ,
				 'base_file_name' : base_file_name ,
				 'ExcSet' : ExcSet ,
				 'ConnectSet' : ConnectSet ,
				 'PERIODIC' : PERIODIC ,
				 'wall_length' : wall_length }
	pickle.dump( itr_dict , f , protocol=2 )
	f.close()

	# Run and check status of the simulation
	if itr > 0:
		os.system('abaqus cae noGUI=GetFEData.py')
		CheckStatus( jn , itr )

	# Update thickness
	no_improve = DesignUpdate( itr , max_shell_thickness , threshold_thickness , max_dt , Objective_type )

	# Check if the design changes
	if no_improve:
		print('No improvements from optimization! Ending process!')
		break

	# Stats
	iter_time = time.time() - st
	print( 'Optimization iteration took ' , iter_time , 's' )

# Final cleanup
for i in range( 1 , 200 ):
	try:
		os.remove( 'abaqus.rpy.' + str(i) )
	except:
		pass