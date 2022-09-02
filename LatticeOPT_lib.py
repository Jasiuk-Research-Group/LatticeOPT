import numpy as np
import pickle
import os
try:
	import matplotlib.pyplot as plt
	from matplotlib import cm
	from sklearn.preprocessing import normalize
	from scipy.sparse import *
	import tensorflow as tf
	from scipy.optimize import minimize
	from scipy.stats import linregress
	import time
	from scipy.optimize import NonlinearConstraint
except:
	pass


# Return the indices of the two nodes of the edge
def EdgeConnectivity( Type , i , j ):
	# [ P1,x   P1,y ]
	# [ P2,x   P2,y ]
	if Type == 0:
		return np.array([i,j,i,j+1]).reshape([2,2])
	elif Type == 1:
		return np.array([i,j,i+1,j]).reshape([2,2])
	else:
		return None


# Plot current design
def PlotDesign( MatDist , itr ):
	fig = plt.figure()
	s = MatDist.shape
	cx = np.linspace( 0. , 1. , s[2] )
	cy = np.linspace( 0. , 1. , s[1] )
	row_bdry = [  0 , 1 , 1 , 1 ]
	col_bdry = [  1 , 0 , 1 , 1 ]
	for Type , rb , cb in zip( range(2) , row_bdry , col_bdry ):
		for i in range( s[1] - rb ):
			for j in range( s[2] - cb ):
				if MatDist[Type,i,j]:
					edge_pts = EdgeConnectivity( Type , i , j )
					plt.plot( cx[ edge_pts[:,1] ] , cy[ edge_pts[:,0] ] , 'r' , linewidth='2' )
	ax = fig.gca()
	# ax.set_xticks(cx)
	# ax.set_yticks(cy)
	ax.set_aspect('equal')
	# plt.grid()
	# plt.show()
	if itr == -1:
		plt.savefig( 'InitialCondition.png' , dpi = 600 )
	else:
		plt.savefig( 'DesignIteration'+str(itr)+'.pdf' , dpi = 600 )
	plt.close()	


# Plot design edges, colored by an array
def PlotEdgeColor( s , cArry , threshold , itr , ExcSet ):
	NP = 30

	fig = plt.figure()
	cx = np.linspace( 0. , 1. , s[1] )
	cy = np.linspace( 0. , 1. , s[0] )
	row_bdry = [  0 , 1 , 1 , 1 ]
	col_bdry = [  1 , 0 , 1 , 1 ]
	eid = 0

	PX , PY , CC = [] , [] , []
	for Type , rb , cb in zip( range(2) , row_bdry , col_bdry ):
		for i in range( s[0] - rb ):
			for j in range( s[1] - cb ):
				if ( threshold == None ) or ( ( threshold != None ) and ( cArry[eid] > threshold ) ) :
					edge_pts = EdgeConnectivity( Type , i , j )
					if ExcSet[ eid ]:
						plt.plot( [cx[ edge_pts[0,1] ] , cx[ edge_pts[1,1] ]] , [cy[ edge_pts[0,0] ] , cy[ edge_pts[1,0] ]] , 'k' )
					else:
						PX += list( np.linspace( cx[ edge_pts[0,1] ] , cx[ edge_pts[1,1] ] , NP ) )
						PY += list( np.linspace( cy[ edge_pts[0,0] ] , cy[ edge_pts[1,0] ] , NP ) )
						CC += list( cArry[eid] * np.ones(NP) )
				eid += 1

	plt.scatter( PX , PY , s=30 , c=CC , cmap = cm.jet )
	fig.gca().set_aspect('equal')
	m = cm.ScalarMappable(cmap=cm.jet)
	m.set_array(CC)
	plt.colorbar(m)

	# Save figure
	if itr == -1:
		plt.savefig( 'InitialCondition.png' , dpi = 600 )
	else:
		plt.savefig( 'ThicknessDesignIteration'+str(itr)+'.png' , dpi = 600 )
	plt.close()


# Convert 3D indices to 1D node ID
def Three2One( i , j  , k , Nx , Ny , Nz ):
	return k * ( Nx+1 ) * ( Ny + 1 ) + j * ( Nx + 1 ) + i + 1 # Node IDs are 1-based


# Write out base inp
def WriteBaseINP( fn , DIM , Ns , WRITE_TIE ):
	UC_Lx , UC_Ly , UC_H = DIM
	Nx , Ny , Nz , N_ele_per_edge = Ns

	# Step 1: Build grid points
	# Note that extra nodes are put in to retain the tensor product structure. These extra nodes (with no elements associated with them) are ignored by Abaqus
	x = np.linspace( 0. , UC_Lx , Nx*N_ele_per_edge+1 )
	y = np.linspace( 0. , UC_Ly , Ny*N_ele_per_edge+1 )
	z = np.linspace( 0. , UC_H  , Nz+1 )
	Z , Y , X = np.meshgrid( z , y , x , indexing='ij'  )

	# Step 2: Write an inp with BCs but no element sets for section assigment
	f = open(fn+'_base.inp','w')
	# Write header
	f.write('*Heading\n*Preprint, echo=NO, model=NO, history=NO, contact=NO\n*Part, name=PART-1\n*Node\n')
	# Write nodes
	for i , x , y , z in zip( range(len(X.flatten())) , X.flatten() , Y.flatten() , Z.flatten() ):
		f.write( str(i+1) + ' , ' + str(x) + ' , ' + str(y) + ' , ' + str(z) + '\n' )

	# Write elements
	ElementBreakPoints = [ 0 ]
	f.write('*Element, type=S4R\n')
	eid = 0
	# All horizontal elements
	for k in range( Nz ):
		for j in range( Ny + 1 ):
			j *= N_ele_per_edge
			for i in range( Nx * N_ele_per_edge ):
				n1 = Three2One( i , j  , k , Nx*N_ele_per_edge , Ny*N_ele_per_edge , Nz )
				n2 = Three2One( i + 1 , j  , k , Nx*N_ele_per_edge , Ny*N_ele_per_edge , Nz )
				n3 = Three2One( i + 1 , j  , k + 1, Nx*N_ele_per_edge , Ny*N_ele_per_edge , Nz )
				n4 = Three2One( i , j  , k + 1, Nx*N_ele_per_edge , Ny*N_ele_per_edge , Nz )
				eid += 1
				f.write( str(eid) + ' , ' + str(n1) + ' , ' + str(n2) + ' , ' + str(n3) + ' , ' + str(n4) + '\n' )
	ElementBreakPoints.append( eid )

	# All vertical elements
	for k in range( Nz ):
		for j in range( Ny * N_ele_per_edge ):	
			for i in range( Nx + 1 ):
				i *= N_ele_per_edge
				n1 = Three2One( i , j  , k , Nx*N_ele_per_edge , Ny*N_ele_per_edge , Nz )
				n2 = Three2One( i , j + 1 , k , Nx*N_ele_per_edge , Ny*N_ele_per_edge , Nz )
				n3 = Three2One( i , j + 1 , k + 1, Nx*N_ele_per_edge , Ny*N_ele_per_edge , Nz )
				n4 = Three2One( i , j  , k + 1, Nx*N_ele_per_edge , Ny*N_ele_per_edge , Nz )
				eid += 1
				f.write( str(eid) + ' , ' + str(n1) + ' , ' + str(n2) + ' , ' + str(n3) + ' , ' + str(n4) + '\n' )
	ElementBreakPoints.append( eid ) # This gives total element count
	print( 'Mesh has a total of ' , eid , ' elements')  
	print( 'Design has a total of ' , ( Nx + 1 ) * Ny + Nx * ( Ny + 1 ) , ' edges' )

	# End this part
	f.write('*End Part\n')

	# Now, fetch BCs from the other inp file
	f2 = open('BoundaryCond.inp','r')
	L = f2.readlines()
	f2.close()
	for l in L:
		f.write(l)
		if '*Assembly, name=Assembly' in l:
			f.write('*Instance, name=PART-1-1, part=PART-1\n*End Instance\n')

			# Write tie constraints
			if WRITE_TIE:
				f.write('*Surface, type=NODE, name=PART-1-1_TOP_TIE_CNS_, internal\nPART-1-1.TOP_TIE, 1.\n')
				f.write('*Surface, type=NODE, name=PART-1-1_BOT_TIE_CNS_, internal\nPART-1-1.BOT_TIE, 1.\n')
				f.write('*Tie, name=Constraint-1, adjust=no\nPART-1-1_TOP_TIE_CNS_, PLATE-1.NEG\n')
				f.write('*Tie, name=Constraint-2, adjust=no\nPART-1-1_BOT_TIE_CNS_, PLATE-2.POS\n')
		if '*End Assembly' in l:
			f.write('*Section Controls, name=EC-1, ELEMENT DELETION=YES\n1., 1., 1.\n')
	f.close()

	# Save break points
	f=open('ElemBreakPt','wb')
	pickle.dump( ElementBreakPoints , f , protocol=2 )
	f.close()


# Write element sets based on MatDist
def ElemSetFormDesign( MatDist , Ns , thickness_dist ):
	Nx , Ny , Nz , N_ele_per_edge = Ns
	# Read saved break points
	f=open('ElemBreakPt','rb')
	ElementBreakPoints = pickle.load(f)
	f.close()
	N_elements = ElementBreakPoints[-1]

	row_bdry = [  0 , 1 , 1 , 1 ]
	col_bdry = [  1 , 0 , 1 , 1 ]
	Edge_sets , Shell_thickness = {} , {}
	for idx , EdgeSet , StartPt in zip( range(2) , MatDist , ElementBreakPoints ): 
		s = EdgeSet.shape
		for i in range( s[0] - row_bdry[idx] ):
			for j in range( s[1] - col_bdry[idx] ):
				my_edge_set , my_edge_thickness = [] , {}

				if idx == 0:
					edge_id = i * Nx + j
					for l in range( N_ele_per_edge ):
						base_ID = i * Nx * N_ele_per_edge + j * N_ele_per_edge + l + 1
						for k in range( Nz ):
							my_element_id = ElementBreakPoints[0] + base_ID + k * ( Ny + 1 ) * Nx * N_ele_per_edge
							my_normalized_z_coord = ( k + 0.5 ) / float( Nz ) # At element center, from 0 to 1
							my_edge_set += [ my_element_id ]
							my_edge_thickness[ my_element_id ] = GetShellThickness( thickness_dist , edge_id , my_normalized_z_coord ) # Store element thickness

				elif idx == 1:
					edge_id = Nx * ( Ny + 1 ) + i * ( Nx + 1 ) + j
					for l in range( N_ele_per_edge ):
						base_ID = ( i * N_ele_per_edge + l ) * ( Nx + 1 ) + j + 1
						for k in range( Nz ):
							my_element_id = ElementBreakPoints[1] + base_ID + k * ( Nx + 1 ) * Ny * N_ele_per_edge
							my_normalized_z_coord = ( k + 0.5 ) / float( Nz ) 
							my_edge_set += [ my_element_id ]
							my_edge_thickness[ my_element_id ] = GetShellThickness( thickness_dist , edge_id , my_normalized_z_coord ) # Store element thickness

				Edge_sets[ edge_id ] = my_edge_set
				Shell_thickness[ edge_id ] = my_edge_thickness
	return Edge_sets , Shell_thickness


# Write set to Abaqus inp
def WriteSET( f , SET ):
	lc = 0; string = ''
	for my_eid in SET:
		string += str( my_eid ) + ', '
		if lc == 7:
			string += '\n'
			lc = 0
		lc += 1
	if string[-1] == ' ':
		string += '\n'
	f.write( string )


# Write element set for section assigment
def SectionAssignment( fn , itr , Ns ):
	Nx , Ny , Nz , N_ele_per_edge = Ns
	MatDist = np.ones( [ 2 , Ny+1 , Nx + 1 ] , dtype=int )

	# Obtain updated design from last iteration
	f = open('DesignIteration' + str( itr - 1 ) + '.npy' , 'rb')
	ph1 , ph2 , _ = np.load( f , allow_pickle=True )
	f.close()
	thickness_dist = _[0]

	# Get edge element sets
	Edge_sets , Shell_thickness = ElemSetFormDesign( MatDist , Ns , thickness_dist )

	# Read template
	f=open(fn+'_base.inp','r')
	Lines=f.readlines()
	f.close()

	# Write new inp with section assignment
	written = False
	f = open(fn+'_itr'+str(itr)+'.inp' ,'w')
	for i in range(len(Lines)):
		l = Lines[i]
		if ('*End Part' in l) and not written:
			written = True
			# Write element sets
			solid_set = []
			for edge_id , edge_set in Edge_sets.items():
				header = '*Elset, elset=Edge-' + str(edge_id) + '\n'
				f.write(header)
				WriteSET( f , edge_set )

				# Write shell thickness distribution
				dist_name = 'EdgeDist-' + str(edge_id)
				f.write('*Distribution, name=' + dist_name + ', location=ELEMENT, Table=' + dist_name + '_Table\n,     1.\n')
				for eid , e_thick in Shell_thickness[ edge_id ].items():
					f.write( str(eid) + ' , ' + str(e_thick) + '\n' )

					# Form a solid set
					if e_thick > 1e-4:
						solid_set += [ eid ]

				# Write shell section with thickness distribution
				f.write( '*Shell Section, elset=Edge-' + str(edge_id) + ', material=MAT, controls=EC-1, shell thickness=' + dist_name + '\n1., 5\n' )
				
			# Write solid element set
			f.write('*Elset, elset=SOLID\n')
			WriteSET( f , solid_set )

			# Write top and bottom tie node set
			f.write('*Nset, nset=BOT_tie\n')
			N_nodes_per_layer = (Nx*N_ele_per_edge+1)*(Ny*N_ele_per_edge+1)
			WriteSET( f , np.arange( 1 , N_nodes_per_layer+1 ) )
			f.write('*Nset, nset=TOP_tie\n')
			WriteSET( f , np.arange( Nz*N_nodes_per_layer+1 , (Nz+1)*N_nodes_per_layer+1 ) )

		f.write(l)
		if ('*End Assembly' in l): # Some book-keeping
			for edge_id in Edge_sets.keys():
				dist_name = 'EdgeDist-' + str(edge_id)
				f.write( '*Distribution Table, name=' + dist_name + '_Table\nlength\n' )
	f.close()
	return fn+'_itr'+str(itr)


# Check simulation status
def CheckStatus( jn , itr ):
	try:
		f = open(jn+'.log' , 'r')
		L = f.readlines()
		f.close()
		if 'COMPLETED' not in L[-1]:
			print('                                                         WARNING: Itr ' + str(itr) + ' failed, continue w/ existing data' )
		else:
			print('                                                         Itr ' + str(itr) + ' completed successfully' )
			# Clean up
			List_delete = [ '.prt' , '.sim' , '.sta' , '.abq' , '.sel' , '.mdl' , '.stt' , '.res' , '.dat' , '.pac' , '.msg' , '.com' , '.log' ]
			for D in List_delete:
				try:
				    os.remove( jn + D )
				except:
				    pass
	except:
		print('                                                         WARNING: Itr ' + str(itr) + ' failed, continue w/ existing data' )


# Compute filter
def ComputeMtx( DIM , Ns , filter_radius ):
	# Return 3 csr matrices:
	# P: Projection mtx, projects edge variables to nodes
	# W: Weighted sum mtx, performs weighted sum of nodal variables
	# C: Connectivity mtx, projects nodal variables to edges
	# Final structure:
	# edge_var_processed = C @ ( W @ ( P @ edge_var ) )
	# And edge_var = [ H edges , V edges ].T

	# Example edge numbering for a 2x2 cell case
	#  
	#      5      6
	#    -----------
	# 10 |   11|   | 12
	#    |  3  | 4 |
	#    -----------
	#  7 |   8 |   | 9
	#    |     |   |
	#    -----------
	#       1     2

	UC_Lx , UC_Ly , UC_H = DIM
	Nx , Ny , Nz , N_ele_per_edge = Ns

	N_edges = ( Nx + 1 ) * Ny + Nx * ( Ny + 1 ) #+ 2 * Nx * Ny
	N_nodes = ( Nx + 1 ) * ( Ny + 1 )

	# P and C matrix
	row_bdry = [  0 , 1 , 1 , 1 ]
	col_bdry = [  1 , 0 , 1 , 1 ]
	edge_id = 0
	pi , pj = [] , []
	for Type , rb , cb in zip( range(2) , row_bdry , col_bdry ):
		for i in range( Ny + 1 - rb ):
			for j in range( Nx + 1 - cb ):
				edge_pts = EdgeConnectivity( Type , i , j )
				node1_idx = edge_pts[0,0] * ( Nx + 1 ) + edge_pts[0,1]
				node2_idx = edge_pts[1,0] * ( Nx + 1 ) + edge_pts[1,1]

				pi.append( node1_idx )
				pj.append( edge_id )
				pi.append( node2_idx )
				pj.append( edge_id )
				edge_id += 1
	pv = np.ones_like( pi )
	Connectivity = coo_matrix( (pv, (pi, pj)), shape=(N_nodes, N_edges))
	P = normalize( Connectivity , norm='l1', axis=1).tocsr() # Normalize row-wise
	C = normalize( Connectivity , norm='l1', axis=0).T.tocsr() # Normalize column-wise
	
	# W matrix
	de , delta = [0,0] , [0,0]
	for i in range(2):
		de[i] = DIM[i] / Ns[i] # Uniform grid size
		# Number of row/cols to go up and down
		delta[i] = int(np.ceil( filter_radius  / de[i] ))

	wi , wj , wv = [] , [] , []
	for i in range( Ny+1 ):
		for j in range( Nx+1 ):
			my_idx_1D = i * ( Nx + 1 ) + j # Node index, 1D

			# Now, search in a square grid for all candidate NODES
			for di in range( max( 0 , i - delta[1] ) , min( i + delta[1] + 2, Ny + 1 ) ):
				for dj in range( max( 0 , j - delta[0] ) , min( j + delta[0] + 2, Nx + 1 ) ):
					curr_node_idx_1D = di * (Nx+1) + dj # Each row has Nx+1 nodes

					r_ij = np.sqrt( np.power( ( di - i )*de[1] , 2. ) + np.power( ( dj - j )*de[0] , 2. ) )
					if r_ij < filter_radius:
						node_wts = filter_radius - r_ij

						wi.append( my_idx_1D )
						wj.append( curr_node_idx_1D )
						wv.append( node_wts )
	W = normalize( coo_matrix( (wv, (wi, wj)), shape=(N_nodes, N_nodes)) , norm='l1', axis=1).tocsr() # Normalize row-wise
	
	f = open('Filter.npy' , 'wb')
	np.save( f , np.array([P,C,W],dtype=object) )
	f.close()


# Compute periodic wrappers
def ComputePeriodicMtx( Nx , Ny , N_cell_x , N_cell_y , ExcSet ):
	# Test if divisible
	if ( Nx % N_cell_x != 0 ) or ( Ny % N_cell_y != 0 ):
		print('Grid size and unit cell size incompatible!')
		exit()

	x_edges_per_cell = Nx // N_cell_x
	y_edges_per_cell = Ny // N_cell_y

	# Modify exclude set ExcSet
	# Horizontal edges
	for row in range( N_cell_y - 1 ):
		offset = ( row + 1 ) * y_edges_per_cell * Nx
		ExcSet[ offset : offset + Nx ] = True
	# Vertical edges
	offset = Nx * ( Ny + 1 )
	for col in range( N_cell_x - 1 ):
		col_offset = offset + ( col + 1 ) * x_edges_per_cell
		for j in range( Ny ):
			ExcSet[ col_offset + j * ( Nx + 1 ) ] = True

	# Compute the periodic mask and its reverse
	free_x_edges_per_cell = x_edges_per_cell - 1

	# Define masks
	forward_mask = np.zeros( len(ExcSet) )

	# Horizontal edges
	offset = Nx * ( Ny + 1 )
	for idx in range( offset ):
		reduced_idx = idx % ( Nx * y_edges_per_cell )
		reduced_idx2 = reduced_idx % x_edges_per_cell
		reduced_row = reduced_idx // Nx
		forward_mask[ idx ] = reduced_idx2 + ( reduced_row - 1 ) * x_edges_per_cell
	
	# Vertical edges
	reduced_offset = int(np.max( forward_mask))
	for row in range( Ny ):
		for col in range( Nx + 1 ):
			reduced_row = row % y_edges_per_cell
			reduced_col = col % x_edges_per_cell
			forward_mask[ offset + row * (Nx+1) + col ] = reduced_offset + reduced_row * free_x_edges_per_cell + reduced_col
	forward_mask[ ExcSet ] = -1

	N_free_edges_in_UC = int(np.max(forward_mask)) + 1
	print( 'Unit cell has ' , N_free_edges_in_UC , ' designable edges' )


	# Build inverse map as a matrix
	inverse_mask = np.zeros( [ len(ExcSet) , N_free_edges_in_UC ] )
	for i , j in enumerate( forward_mask ):
		j = int(j)
		if j == -1:
			continue
		inverse_mask[ i , j ] = 1.

	# Forward map as a matrix
	forward_mask = normalize( inverse_mask.T , norm='l1', axis=1)

	f = open('PeriodicFilter.npy' , 'wb')
	np.save( f , np.array([forward_mask,inverse_mask],dtype=object) )
	f.close()


# Get exclusion set for boundary edges
def GetExclusionSet( Ns ):
	Nx , Ny , Nz , N_ele_per_edge = Ns
	N_edges = ( Nx + 1 ) * Ny + Nx * ( Ny + 1 )
	ExcSet = np.zeros( N_edges , dtype=bool )

	# All bottom horizontal edges
	ExcSet[ : Nx ] = True

	# All top horizontal edges
	ExcSet[ Nx * Ny : Nx * ( Ny + 1 ) ] = True

	# All left and right vertical edges
	for i in range( Ny ):
		offset = Nx * ( Ny + 1 )
		ExcSet[ offset + i * ( Nx + 1 ) ] = True
		ExcSet[ offset + (i+1) * ( Nx + 1 ) - 1 ] = True
	return ExcSet


# Get set for boundary-edge connectors
def GetConnectionSet( Ns , N_layers ):
	Nx , Ny , Nz , N_ele_per_edge = Ns
	N_edges = ( Nx + 1 ) * Ny + Nx * ( Ny + 1 )
	offset = Nx * ( Ny + 1 )
	ConnectSet = np.zeros( N_edges )

	# All N_layers of bottom edges, excluding the boundary
	ConnectSet[ Nx : Nx * (N_layers+1) ] = 1.
	for L in range( N_layers ):
		ConnectSet[ offset + 1 + L * (Nx+1) : offset + 1 + L * (Nx+1) + (Nx-1) ] = 1.

	# All N_layers of top horizontal edges, excluding the boundary
	ConnectSet[ Nx * Ny - Nx * N_layers : Nx * Ny ] = 1.
	for L in range( Ny - N_layers , Ny ):
		ConnectSet[ offset + 1 + L * (Nx+1) : offset + 1 + L * (Nx+1) + (Nx-1) ] = 1.

	# All N_layers of left and right vertical edges, excluding the boundary
	for i in range( 1 + N_layers , Ny - N_layers ):
		ConnectSet[ i * Nx : i * Nx + N_layers ] = 1.
		ConnectSet[ (i+1) * Nx - N_layers : (i+1) * Nx  ] = 1.

	for i in range( N_layers , Ny - N_layers ):
		ConnectSet[ 1 + offset + i * (Nx+1) : 1 + offset + i * (Nx+1) + N_layers ] = 1.
		ConnectSet[ -1 + offset + (i+1) * (Nx+1) - N_layers : offset + (i+1) * (Nx+1) - 1  ] = 1.
	return ConnectSet


# Get shell thickness as a function of normalized Z coord
def GetShellThickness( thickness_dist , edge_id , my_normalized_z_coord ):
	control_pts = thickness_dist[ edge_id , : ]
	if len( control_pts ) == 1: # If constant thickness along z direction
		return control_pts[0]
	else:
		x = np.linspace(0.,1.,len(control_pts))
		return np.interp( my_normalized_z_coord , x , control_pts )


# Thickness updating scheme
def UpdateThickness( edge_energy , t_old , thickness_weight , wl , m_new , t_max , t_trheshold , max_dt , Objective_type , max_inner_solver_iter , itr ):
	N = len( t_old ) # Get total number of edges
	curr_energy = tf.convert_to_tensor(edge_energy,dtype=np.float64)
	curr_thickness = tf.convert_to_tensor(t_old,dtype=np.float64)
	t_weight = tf.convert_to_tensor(thickness_weight,dtype=np.float64)
	curr_length = tf.convert_to_tensor(wl,dtype=np.float64)

	##############################  Objective function  ################################################
	global obj_list # Use this to monitor objective function progress
	obj_list = []
	if Objective_type == 1:
		def tfOBJ( inp ):
			return tf.math.reduce_std( tf.math.divide( curr_energy , tf.math.multiply( inp , curr_length ) ) )

	elif Objective_type == 2:
		def tfOBJ( inp ):
			edge_energy_sum = np.sum( edge_energy )
			curr_area = tf.math.multiply( inp , curr_length )
			area_sum = tf.math.reduce_sum( curr_area )
			# BESO sensitivity = V_e / V_tot - E_e / E_tot
			sensitivity = (curr_area / area_sum - edge_energy / edge_energy_sum)
			return tf.math.reduce_std( sensitivity ) * 1000.

	def OBJ_jac( X ):
		inp = tf.Variable( X , dtype=tf.float64 )
		with tf.GradientTape() as tape:
			std = tfOBJ( inp )
			obj_list.append( std.numpy() ) # Store objective value every time gradient is evaluated 
		# Get gradient 
		return tape.gradient( std , inp ).numpy()

	def OBJ_hess( X ):
		inp = tf.Variable( X , dtype=tf.float64 )
		with tf.GradientTape() as tape:
			with tf.GradientTape() as tape2:
				std = tfOBJ( inp )
			# Get gradient 
			grad = tape2.gradient( std , inp )
		# Get hess
		return tape.gradient( grad , inp ).numpy()

	def OBJ( X ):
		inp = tf.convert_to_tensor(X,dtype=np.float64)
		return tfOBJ( inp ).numpy()

	# We use this call back to do early stopping if objective is increasing 
	def OBJ_watch_dog( arg1 , arg2 ):
		# n_check = 5
		# if len( obj_list ) > 2 * n_check:
		# 	if np.mean( obj_list[-n_check:] ) > np.mean( obj_list[-2*n_check:-n_check] ):
		# 		print('Objective function increasing, early stopping...')
		# 		return True
		return False

	##############################  Mass constraint  ################################################
	def tfMass( inp ):
		return tf.math.reduce_sum( tf.math.multiply( inp , curr_length ) )

	def Mass_jac( X ):
		inp = tf.Variable( X , dtype=tf.float64 )
		with tf.GradientTape() as tape:
			mass = tfMass( inp )
		# Get gradient 
		return tape.gradient( mass , inp ).numpy()

	def Mass_hess( X , v ):
		nc = len(X)
		return np.zeros([ nc , nc ])

	def Mass( X ):
		inp = tf.convert_to_tensor(X,dtype=np.float64)
		return tfMass( inp ).numpy()

	# Mass constraint is enforced as an equality constraint
	mass_cons = NonlinearConstraint( Mass , m_new , m_new , jac=Mass_jac , hess=Mass_hess)


	##############################  Connectedness constraint  ################################################
	connect_threshold = 0.95
	def tfConnect( inp ):
		normalized_weights = tf.math.divide( t_weight , curr_thickness )
		return tf.math.reduce_mean( tf.math.multiply( normalized_weights , inp ) )

	def Connect_jac( X ):
		inp = tf.Variable( X , dtype=tf.float64 )
		with tf.GradientTape() as tape:
			connect = tfConnect( inp )
		# Get gradient 
		return tape.gradient( connect , inp ).numpy()

	def Connect_hess( X , v ):
		# inp = tf.Variable( X , dtype=tf.float64 )
		# with tf.GradientTape() as tape:
		# 	with tf.GradientTape() as tape2:
		# 		connect = tfConnect( inp )
		# 	# Get gradient 
		# 	grad = tape2.gradient( connect , inp )
		# # Get hess
		# return v[0] * tape.gradient( grad , inp ).numpy()
		nc = len(X)
		return np.zeros([ nc , nc ])

	def Connect( X ):
		inp = tf.convert_to_tensor(X,dtype=np.float64)
		return tfConnect( inp ).numpy()

	connect_cons = NonlinearConstraint( Connect , connect_threshold , 1.25 , jac=Connect_jac, hess=Connect_hess)

	# # Define constraints
	# if np.sum( thickness_weight ) > 1e-6:
	# 	cons = [ mass_cons , connect_cons ]
	# else:
	# 	cons = [ mass_cons ]
	cons = [ mass_cons ]


	# Define bounds
	bnds = []
	for i in range( N ):
		lb = max( 0. , t_old[i] - max_dt )
		ub = min( t_max , t_old[i] + max_dt )
		bnds.append( ( lb , ub ) )

	# Optimize
	res = minimize(OBJ, t_old, method='trust-constr', jac=OBJ_jac, hess = OBJ_hess,
		constraints=cons, options={'gtol': 1e-5, 'disp': True,'maxiter': max_inner_solver_iter},bounds=bnds , callback = OBJ_watch_dog )
	print('Final function value is : ' , '{:1.4f}'.format(res.fun)  )
	print('Disconnectedness metric is : ' , '{:1.4f}'.format( Connect( res.x ) )  )

	# Threshold
	t_new = (res.x).copy()
	remove_flag = ( t_new < t_trheshold )
	t_new[ remove_flag ] = 1.e-6
	return t_new


# Update design
def DesignUpdate( itr , t_max , t_trheshold , max_dt , Objective_type ):
	# Obtain updated design from last iteration
	f = open('DesignIteration' + str( itr - 1 ) + '.npy' , 'rb')
	edge_energy_old , t_old , _ = np.load( f , allow_pickle=True )
	f.close()
	thickness_dist = _[0]

	# Get current iteration info
	f=open('currItr','rb')
	itr_dict = pickle.load(f)
	itr_dict['jn'] = str( itr_dict['jn'] ); itr_dict['base_file_name'] = str(itr_dict['base_file_name'])
	f.close()
	Nx , Ny , Nz , N_ele_per_edge = itr_dict['Ns']
	IncSet = np.logical_not( itr_dict['ExcSet'] ) # All edges that need to be included in optimization
	wall_length = itr_dict['wall_length'] # Length of all lattice walls

	# Read in raw segment energies
	f=open('Energies_'+itr_dict['jn'],'rb')
	curr_segment_energy ,curr_segment_ese,curr_segment_epd,curr_segment_edmd , curr_edge_status = np.load( f )
	f.close()
	# Convert to edge energy, which is the sum of all segments in that edge
	curr_energy = np.sum( curr_segment_energy , axis = -1 ) # Segments are in the last axis
	curr_ese = np.sum( curr_segment_ese , axis = -1 )
	curr_pd = np.sum( curr_segment_epd , axis = -1 )
	curr_dmd = np.sum( curr_segment_edmd , axis = -1 )

	# Get edge status at last increment
	last_edge_status = 1. - np.sum( curr_edge_status , axis = -1 )[-1,:] # Percent disconnectedness
	# last_edge_status = np.sum( curr_edge_status , axis = -1 )[-1,:] # Percent connectedness


	# Total energy at last step
	tot_e = np.sum( curr_energy[-1,:] )
	# objFun and mass ratio 
	curr_sum_t = np.sum( np.multiply( t_old[IncSet] , wall_length[IncSet] ) ) # When calculating current sum of thickness, only include the edges in the optimization
	objFun = tot_e / np.sum( np.multiply( t_old , wall_length ) ) # When calculating SEA, include all edges
	print( 'Current specific absorption = ' , '{:1.3f}'.format(objFun) )


	# Specific absorption per edge, computed using edge end energy
	end_energy = DefineEndEnergy( curr_energy , curr_ese , curr_pd , curr_dmd )

	# Average with old data
	if itr > 0:
		f = open('DesignIteration' + str( itr - 2 ) + '.npy' , 'rb')
		ph1 , thickness_prev , _ = np.load( f , allow_pickle=True )
		f.close()

		for eid in range( len(t_old) ):
			# Check if an edge is removed or added
			#            Last iteration                        Current iteration
			if ( thickness_prev[eid] > t_trheshold ) == ( t_old[eid] > t_trheshold ):
				# Only average if the statuses of the edge are same in both steps 
				end_energy[ eid ] = 0.5 * ( end_energy[ eid ] + edge_energy_old[ eid ] )

	# Prepare connectedness constraint
	# For the constraint, true edge thickness is weighted by its percent connectedness
	thickness_weight = itr_dict['ConnectSet'] * last_edge_status

	# Thickness update through optimization
	thickness_new = t_old.copy()
	# Note that we only pass the data for all included edges to the optimization

	# Data arrays for optimization
	if itr_dict['PERIODIC']:
		# Load periodic wrapper
		f = open('PeriodicFilter.npy' , 'rb')
		forward_mask,inverse_mask = np.load( f , allow_pickle=True )
		f.close()

		# Get periodic image
		energy = forward_mask @ end_energy
		t_guess = forward_mask @ t_old
		wl = forward_mask @ wall_length
		t_wt = forward_mask @ thickness_weight
		new_mass = itr_dict['m_new'] / np.sum( inverse_mask[:,0] ) # Mass of each unit cell

	else:
		energy = end_energy[ IncSet ] # Energy measure
		t_guess = t_old[ IncSet ] # Initial guess for thickness
		t_wt = thickness_weight[ IncSet ] # Current thickness weight
		wl = wall_length[ IncSet ]
		new_mass = itr_dict['m_new']

	# Perform optimization
	updated_values = UpdateThickness( energy , t_guess , t_wt , wl , \
									new_mass , t_max , t_trheshold , max_dt , Objective_type , \
									itr_dict['max_inner_solver_iter'] , itr )

	# Put updated values back
	if itr_dict['PERIODIC']:
		t_periodic = inverse_mask @ updated_values
		thickness_new[ IncSet ] = t_periodic[ IncSet ]
	else:
		thickness_new[ IncSet ] = updated_values


	# Apply linear mass scaling
	ratio = np.sum( np.multiply( thickness_new[IncSet] , wall_length[IncSet] ) ) / itr_dict['m_new']
	thickness_new[IncSet] /= ratio

	# Observed max thickness change
	observed_max = np.max( np.abs( thickness_new[IncSet] - t_old[IncSet] ) )
	print( '\nSanity check:\nMin wall thickness = ' , '{:1.3f}'.format(np.min( thickness_new[IncSet])) , ', Max wall thickness = ' , '{:1.3f}'.format(np.max(thickness_new[IncSet])) )
	print( 'Max thickness change = ' , '{:1.4f}'.format(observed_max) , '\nActual mass / expected mass = ' , '{:1.2f}'.format( ratio * 100. ) )

	# Update thickness distribution
	thickness_dist_new = np.ones( [ len(thickness_new) , 1 ] )
	thickness_dist_new[ : , 0 ] = thickness_new

	
	# Save the updated design
	f = open('DesignIteration'+str(itr)+'.npy' , 'wb')
	np.save( f , np.array([end_energy,thickness_new,[thickness_dist_new]],dtype=object) )
	f.close()

	# Write log
	f=open('OptLog_' + itr_dict['base_file_name'] + '.txt','a')
	line = str(itr) + ',' + str(objFun) + ',' + str(curr_sum_t) + ',' + str(observed_max) + '\n'
	f.write(line)
	f.close()

	# Save design figure
	PlotEdgeColor( [Ny+1,Nx+1] , thickness_new , t_trheshold , itr , itr_dict['ExcSet'] )
	return np.allclose( t_old , thickness_new , atol = 1e-2 )


# Define end energy
def DefineEndEnergy( curr_energy , curr_ese , curr_pd , curr_dmd ):
	# Just use the last total energy
	# end_energy = curr_energy[-1]

	# Use the mean of the product of PD and DMD
	# end_energy = curr_pd[-1] * curr_dmd[-1]

	# Use mean total energy
	# end_energy = np.mean( curr_energy , axis = 0 )	

	# Weighted average of total energy
	ss = curr_energy.shape
	weights = np.linspace( 0.25 , 1. , ss[0] )
	weights /= np.sum( weights ) # Normalize
	end_energy = np.zeros( ss[1] )
	for eid in range( ss[1] ):
		end_energy[ eid ] = np.dot( curr_energy[ : , eid ] , weights )
	return end_energy


# Get indices for weighted edge graph
def GetEdgeIndices( s ):
	N_nodes = np.prod(s)
	# EdgeGraph = np.zeros( [ N_nodes , N_nodes ] )
	row_bdry = [  0 , 1 , 1 , 1 ]
	col_bdry = [  1 , 0 , 1 , 1 ]
	idx1 , idx2 = [],[]
	for Type , rb , cb in zip( range(2) , row_bdry , col_bdry ):
		for i in range( s[0] - rb ):
			for j in range( s[1] - cb ):
				edge_pts = EdgeConnectivity( Type , i , j )
				p1_1d_idx = edge_pts[0,0] * s[1] + edge_pts[0,1]
				p2_1d_idx = edge_pts[1,0] * s[1] + edge_pts[1,1]

				idx1.append( [p1_1d_idx , p2_1d_idx] )
				idx2.append( [p2_1d_idx , p1_1d_idx] )
	return idx1 + idx2