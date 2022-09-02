from LatticeOPT_lib import *

Plot_Secondary = True
base = 'case3_obj1_max25'


n_itr = int( base[-2:] )
itr = np.arange(n_itr)

# Read core energy data
tot_core_e , tot_core_pd , tot_core_dmd , mean_edge_status = np.zeros([n_itr,20]), np.zeros([n_itr,20]), np.zeros([n_itr,20]), np.zeros([n_itr,20])
for i in range( n_itr ):
	jn = base + '_itr' + str(i)
	f=open('Energies_'+jn,'rb')
	curr_segment_energy ,curr_segment_ese,curr_segment_epd,curr_segment_edmd , curr_edge_status = np.load( f )
	f.close()

	LL = curr_segment_energy.shape[0]

	tot_core_e[i,-LL:] = np.sum( curr_segment_energy[:,:,0] , axis = 1 ) / 1000.
	tot_core_pd[i,-LL:] = np.sum( curr_segment_epd[:,:,0] , axis = 1 ) / 1000.
	tot_core_dmd[i,-LL:] = np.sum( curr_segment_edmd[:,:,0] , axis = 1 ) / 1000.
	mean_edge_status[i,-LL:] = np.mean( curr_edge_status[:,:,0] , axis = 1 )

# # Read total energy data
# PD_ , DMD_ , SE_ , KE_ , WK_ = np.load( 'HistoryOutput.npy' )
# tot_e , tot_pd , tot_dmd , tot_ke , tot_wk = np.zeros([n_itr,20]), np.zeros([n_itr,20]), np.zeros([n_itr,20]), np.zeros([n_itr,20]), np.zeros([n_itr,20])

# # Interp to same shape
# x_in = np.linspace( 0. , 1. , 201 )
# x_out = np.linspace( 0. , 1. , 20 )

# for i in range( 25 ):
# 	tot_e[i,:] = np.interp( x_out , x_in , PD_[i,:] + DMD_[i,:] + SE_[i,:] )
# 	tot_pd[i,:] = np.interp( x_out , x_in , PD_[i,:] )
# 	tot_dmd[i,:] = np.interp( x_out , x_in , DMD_[i,:] )
# 	tot_ke[i,:] = np.interp( x_out , x_in , KE_[i,:] )
# 	tot_wk[i,:] = np.interp( x_out , x_in , WK_[i,:] )

# # Get face sheet energy
# tot_face_e = tot_e - tot_core_e
# tot_face_pd = tot_pd - tot_core_pd
# tot_face_dmd = tot_dmd - tot_core_dmd


# Figure 1: Core energy absorption stats
plt.figure()
y = tot_core_e[:,-1] / tot_core_e[0,-1]
plt.plot( itr , y , 'tab:orange' )
# mid = np.argmin(y)
# plt.plot( itr[mid] , y[mid] , 'rp' )
y = tot_core_dmd[:,-1] / tot_core_dmd[0,-1]
plt.plot( itr , y , 'tab:purple' )
# plt.plot( itr[mid] , y[mid] , 'kd' )
y = tot_core_pd[:,-1] / tot_core_pd[0,-1]
plt.plot( itr , y , 'tab:gray' )
plt.xlabel('Design iterations', fontsize=16)
plt.ylabel('Normalized energies [-]', fontsize=16)
plt.legend(['Core energy absorption' , 'Core damage dissipation' , 'Core plastic dissipation' ], fontsize=13)
plt.savefig( 'CoreStats' + '.pdf', format = 'pdf', dpi=600 , bbox_inches="tight" )

# Figure 2: Core edge connectedness
plt.figure()
plt.plot( itr , mean_edge_status[:,4] , 'tab:orange' )
plt.plot( itr , mean_edge_status[:,9] , 'tab:purple' )
plt.plot( itr , mean_edge_status[:,14] , 'tab:gray' )
plt.plot( itr , mean_edge_status[:,-1] , 'tab:olive' )
plt.xlabel('Design iterations', fontsize=16)
plt.ylabel('Mean edge connectedness measure [-]', fontsize=14)
plt.legend(['t=0.025 ms' , 't=0.05 ms' , 't=0.075 ms' , 't=0.1 ms' ], fontsize=13)
plt.savefig( 'CoreConnect' + '.pdf', format = 'pdf', dpi=600 , bbox_inches="tight" )


# Figure 3: Whole system energy
plt.figure()
y = tot_wk[:,-1] / tot_wk[0,-1]
plt.plot( itr , y , 'tab:orange' )
mid = np.argmin(y)
plt.plot( itr[mid] , y[mid] , 'rp' )

y = tot_dmd[:,-1] / tot_dmd[0,-1]
plt.plot( itr , y , 'tab:purple' )
plt.plot( itr[mid] , y[mid] , 'kd' )

y = tot_pd[:,-1] / tot_pd[0,-1]
plt.plot( itr , y , 'tab:gray' )

plt.xlabel('Design iterations', fontsize=16)
plt.ylabel('Normalized energies [-]', fontsize=16)
plt.legend(['External work, system' , 'Minimum compliance design' , 'Damage dissipation, system' , 'Minimum damage design' , 'Plastic dissipation, system' ], fontsize=13)
plt.savefig( 'SystemStats' + '.pdf', format = 'pdf', dpi=600 , bbox_inches="tight" )


# fig, axs = plt.subplots(2, 2, sharex='row')
# s = 0 # Number of points to skip
# axs[0,0].plot(itr[s:],tot_e[s:],'k-o')
# axs[0,0].set_title('Total energy absorption [J]')
# axs[0,1].plot(itr[s:],tot_pd[s:],'k-s')
# axs[0,1].set_title('Total plastic energy [J]')
# axs[1,0].plot(itr[s:],tot_dmd[s:],'k-s')
# axs[1,0].set_title('Total damage energy [J]')
# axs[1,1].plot(itr[s:],mean_edge_status[s:],'k-s')
# axs[1,1].set_title('Mean edge connectedness')

skip = 10
special_set = [ 0 , np.argmin(tot_core_e[:,-1]) , skip+np.argmax(tot_core_e[:,-1][skip:]) ]
print( special_set )


if Plot_Secondary:
	# Read energy data
	tot_e , tot_pd , tot_dmd = np.zeros([3,5]) , np.zeros([3,5]) , np.zeros([3,5])
	for i in range( 3 ):
		for j in range( 5 ):
			jn = base + '_itr' + str( special_set[i] ) + '-mul-' + str(j)
			f=open('Energies_'+jn,'rb')
			curr_segment_energy ,curr_segment_ese,curr_segment_epd,curr_segment_edmd , curr_edge_status = np.load( f )
			f.close()

			tot_e[i,j] = ( np.sum(curr_segment_energy[-1,:,:]) / 1000. )
			tot_pd[i,j] = ( np.sum(curr_segment_epd[-1,:,:]) / 1000. )
			tot_dmd[i,j] = ( np.sum(curr_segment_edmd[-1,:,:]) / 1000. )

	cc = [ 'k-s' , 'g-x' , 'm-d' ]
	x = np.linspace( 1. , 2. , 5 )
	fig, axs = plt.subplots(1, 3, sharex='row')
	for i in range( 3 ):
		axs[0].plot( x ,tot_e[i,:] , cc[i] )
		axs[1].plot( x ,tot_pd[i,:] , cc[i] )
		axs[2].plot( x ,tot_dmd[i,:] , cc[i] )

	axs[0].set_title('Total energy absorption [J]')
	axs[1].set_title('Total plastic energy [J]')
	axs[2].set_title('Total damage energy [J]')

	for ax in axs:
		ax.legend(['Initial','Lowest SEA','Optimized'])

plt.show()