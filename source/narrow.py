import os
import sys
import numpy as np
from pdb import set_trace
from matplotlib import rc as matplotlibrc
from matplotlib.patches import FancyArrowPatch
import matplotlib.pyplot as plt
# plt.switch_backend('agg')
from itertools import product, combinations
from scipy import signal
import pickle
from scipy.interpolate import interp1d

figheight      = 6
figwidth      = 6
lineWidth      = 1.5
textFontSize   = 18
gcafontSize    = 15
legendFontSize = 20

TRACE 	 = 1


if TRACE:


	nfile=50
	xa_list=np.zeros(nfile)
	ya_list=np.zeros(nfile)	
	t_list =np.zeros(nfile)	

	xa_list2=np.zeros(nfile)
	ya_list2=np.zeros(nfile)	
	t_list2 =np.zeros(nfile)

	for i in range (0,nfile):
		pkl_name	= './data/data_'+str(i+1)+'.pkl'
		data 		= pickle.load( open( pkl_name, "rb" ),encoding="latin1" )

		N 		= data['N']	
		t 		= data['time']				
		z 		= data['z']		

		x=np.real(z)
		y=np.imag(z)
		xh=x[0:np.int(N/2)]
		yh=y[0:np.int(N/2)]		

		interpolate_y = interp1d(xh,yh, kind='cubic')
		ya_list[i]	= interpolate_y(0.)*2
		xa_list[i]	= np.max(x)*2
		t_list[i] 	= t

		#######################################################################

		pkl_name	= './data/data2_'+str(i+1)+'.pkl'
		data 		= pickle.load( open( pkl_name, "rb" ),encoding="latin1" )

		N 		= data['N']	
		t 		= data['time']				
		z 		= data['z']		

		x=np.real(z)
		y=np.imag(z)
		xh=x[0:np.int(N/2)]
		yh=y[0:np.int(N/2)]		

		interpolate_y = interp1d(xh,yh, kind='cubic')
		ya_list2[i]	= interpolate_y(0.)*2
		xa_list2[i]	= np.max(x)*2
		t_list2[i] 	= t

	fig = plt.figure(0, figsize=(figwidth*1.6,figheight*0.6))
	ax  = fig.add_subplot(121,alpha=1)#	

	plt.plot(t_list,ya_list,'k--')
	plt.plot(t_list2,ya_list2,'k-')
	ax.set_xlabel('t',fontsize=1.2*gcafontSize)
	ax.set_ylabel('h',fontsize=1.2*gcafontSize)

	ax  = fig.add_subplot(122,alpha=1)#	
	plt.plot(t_list,xa_list,'k--')
	plt.plot(t_list2,xa_list2,'k-')
	ax.set_xlabel('t',fontsize=1.2*gcafontSize)
	ax.set_ylabel('w',fontsize=1.2*gcafontSize)

	figname='test.png'
	plt.tight_layout()		
	plt.savefig(figname)
	# plt.show()
	plt.close()
	print (figname," saved!")

	# Qc=2e-3
	# T=(5e-3)**2*(2e-3)*(2*np.pi)/Qc
	# tf=t*T
	# tf=t*1.0
	# set_trace()
	# time_string='tf = '+str(round(tf,3))+' s'
	# N_string='N = '+str(len(z))
	# ax.text(-9,-9,time_string,fontsize=gcafontSize)
	# ax.text(0,-9,N_string,fontsize=gcafontSize)	
	# ax.set_xlabel('x',fontsize=1.2*gcafontSize)
	# ax.set_ylabel('y',fontsize=1.2*gcafontSize)
