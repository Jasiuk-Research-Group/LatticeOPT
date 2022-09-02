# Import all libraries
from abaqus import *
from abaqusConstants import *
from caeModules import *
from odbAccess import *
import visualization
import mesh
import regionToolset
import numpy as np
import pickle
from os.path import exists

f=open('currItr','rb')
itr_dict = pickle.load(f)
itr_dict['jn'] = str( itr_dict['jn'] )
f.close()

# Read saved break points
f=open('ElemBreakPt','rb')
ElementBreakPoints = pickle.load(f)
f.close()

# Compute constants
itr_dict['Ns'] = [ int(_) for _ in itr_dict['Ns'] ]
Nx , Ny , Nz , N_ele_per_edge = itr_dict['Ns']
N_ele_per_Z_segment = Nz
N_edges = ( Nx + 1 ) * Ny + Nx * ( Ny + 1 )
type_edges_per_layer = [ Nx * ( Ny + 1 ) , ( Nx + 1 ) * Ny ]
offsets = [ 0 , Nx * ( Ny + 1 ) ]

INST_NAME = 'PART-1-1'

# Run simulation in Abaqus
mdb.JobFromInputFile(name=itr_dict['jn'], 
    inputFileName=itr_dict['jn']+'.inp', 
    type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, 
    memory=90, memoryUnits=PERCENTAGE, explicitPrecision=SINGLE, 
    nodalOutputPrecision=SINGLE, userSubroutine='', scratch='', 
    resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, numDomains=itr_dict['ncpu'], 
    activateLoadBalancing=False, multiprocessingMode=DEFAULT, numCpus=itr_dict['ncpu'])

# Submit
mdb.jobs[itr_dict['jn']].submit(consistencyChecking=OFF)
mdb.jobs[itr_dict['jn']].waitForCompletion()


# Open target odb
Odb = session.openOdb(name= itr_dict['jn'] + '.odb' )

# Go through all time steps
N_frames = len(Odb.steps['Step-1'].frames)

# Store energy quantities per edge
curr_energy = np.zeros( [ N_frames - 1 , N_edges , 1 ] )
curr_ese = np.zeros_like( curr_energy )
curr_epd = np.zeros_like( curr_energy )
curr_edmd = np.zeros_like( curr_energy )
curr_edge_status_tmp = np.ones( [ N_frames - 1 , N_edges , Nz ] )

# print('\n\n\n\nNEW OUTPUT:')
for f in range( 1 , N_frames ): # Let's skip the first time step, which is just all zeros
	# Frame to get result from
	FRAME = Odb.steps['Step-1'].frames[f]

	ESE = FRAME.fieldOutputs['ELSE'].values
	EPD = FRAME.fieldOutputs['ELPD'].values
	EDMD = FRAME.fieldOutputs['ELDMD'].values
	STATUS = FRAME.fieldOutputs['STATUS'].values

	for e in range( len( ESE ) ):
		instance_name = ESE[e].instance.name
		if instance_name == INST_NAME:
			# Current element label, 1 based
			elemLabel = ESE[e].elementLabel

			# Find out where this entry should go
			for g in range( 2 ):
				if elemLabel <= ElementBreakPoints[g+1]: # Determine edge type
					shifted_index = ( elemLabel - ElementBreakPoints[g] - 1 ) % ( type_edges_per_layer[g] * N_ele_per_edge )
					layer_index = ( elemLabel - ElementBreakPoints[g] - 1 ) // ( type_edges_per_layer[g] * N_ele_per_edge )
					segment_id = layer_index // N_ele_per_Z_segment

					if g == 0:
						edge_idx_1d = offsets[g] + shifted_index // N_ele_per_edge
					else:
						row = shifted_index // ( ( Nx + 1 ) * N_ele_per_edge )
						col = shifted_index % ( Nx + 1 )
						edge_idx_1d = offsets[g] + row * ( Nx + 1 ) + col

					# All energy components
					ese = ESE[e].data
					pd = EPD[e].data
					dmd = EDMD[e].data
					stat = STATUS[e].data

					curr_energy[ f - 1 , edge_idx_1d , segment_id ] += ( ese + pd + dmd ) # Sum elastic and plastic and damage strain energies

					# Store individual components
					curr_ese[ f - 1 , edge_idx_1d , segment_id ] += ese
					curr_epd[ f - 1 , edge_idx_1d , segment_id ] += pd
					curr_edmd[ f - 1 , edge_idx_1d , segment_id ] += dmd

					# Check if the element has failed or not
					if stat < 0.5:
						# If an element has failed, mark the corresponding layer as disconnected
						curr_edge_status_tmp[ f - 1 , edge_idx_1d , layer_index ] = 0. 
					break

	# Symmetrize segment data
	curr_energy[ f - 1 , : , : ] = 0.5 * ( curr_energy[ f - 1 , : , : ] + np.flip( curr_energy[ f - 1 , : , : ] , axis = -1 ) )

# Average over all layers
curr_edge_status = np.expand_dims( np.mean( curr_edge_status_tmp , axis = -1 ) , axis = -1 )

f=open('Energies_'+itr_dict['jn'],'wb')
np.save( f , [curr_energy,curr_ese,curr_epd,curr_edmd,curr_edge_status] )
f.close()