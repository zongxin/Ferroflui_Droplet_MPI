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

import module as mod
import parameter as para

# matplotlibrc('text.latex', preamble='\usepackage{color}')
# matplotlibrc('text',usetex=True)
# matplotlibrc('font', family='serif')

figheight      = 6
figwidth       = 6
lineWidth      = 1.5
textFontSize   = 14
gcafontSize    = 10
legendFontSize = 10

def write_figure(t):

	z 		 	=	mod.zg
	z0 		 	=	mod.z0
	kappa		=	mod.kappa
	alpha		=	para.alpha
	R 			=	para.R

	fig = plt.figure(0, figsize=(figwidth*1.3,figheight*1.3))
	ax  = plt.subplot(221,alpha=1)#

	ax.grid(False)
	ax.set_xlim([-1.5,1.5])
	ax.set_ylim([-1.5,1.5])	
	# plt.axis('equal')	
	plt.plot(np.real(z0),np.imag(z0),'r')		
	plt.plot(np.real(z),np.imag(z),'k')
	plt.scatter(np.real(z[::2]),np.imag(z[::2]),c='k',marker='.',s=8.)	
		
	ax.set_xlabel('x',fontsize=1.2*gcafontSize)
	ax.set_ylabel('y',fontsize=1.2*gcafontSize)

	time_string='t='+str(round(t,5))	
	ax.text(0.5,-1.2,time_string, fontsize=1.*gcafontSize)

	ax.set_aspect('equal')



	ax  = plt.subplot(222,alpha=1)#
	plt.plot(alpha,kappa,'k',lw=lineWidth)


	ax.set_xlabel('alpha',fontsize=1.2*gcafontSize)
	ax.set_ylabel('kappa',fontsize=1.2*gcafontSize)



	ax  = plt.subplot(223,alpha=1)#
	m=40
	z_hat=np.abs(np.fft.fft(np.abs(z)))/len(z)
	z_hat0=np.abs(np.fft.fft(np.abs(z0)))/len(z)

	plt.plot(np.array(range(1,m)),z_hat[1:m],'k',marker='.',linewidth=lineWidth)
	plt.plot(np.array(range(1,m)),z_hat0[1:m],'r--',marker='.')	


	ax  = plt.subplot(224,alpha=1)#
	# plt.axes(ax)
	plt.plot(para.alpha,np.abs(z)-R,'k',linewidth=lineWidth)
	plt.plot(para.alpha,np.abs(mod.z0)-R,'r--')	
	ax.set_xlabel('alpha',fontsize=1.2*gcafontSize)
	ax.set_ylabel('r',fontsize=1.2*gcafontSize)
	figname="frame_"+str(mod.plotj)+'.png'

	plt.tight_layout()		
	plt.savefig('./movie/'+figname)
	# plt.show()
	plt.close()
	print ("\n*****************************************************"	)
	print (figname," saved!")
	mod.plotj=mod.plotj+1

def write_data(t):

	pkl_name='./data/data2_'+str(mod.dataj)+'.pkl'
	data={}
	data['time']	= t	
	data['N']		= para.Ng	
	data['z0']		= mod.z0		
	data['z']		= mod.zg
	data['gamma']	= mod.gamma	
	data['kappa']	= mod.kappa	
	# data['tracer']	= tracer				
	output = open(pkl_name, 'wb')
	pickle.dump(data, output)
	output.close()
	print (pkl_name," saved!")
	print ("*****************************************************"	)
	mod.dataj=mod.dataj+1

def print_screen(i,t):
	print ("\n====================================================\n")
	print ("TIME ADVANCE	:	",i)
	print ('TIME 		:	',t)
	print ('MAX 	VEL 	:	',np.max(mod.umax))#*para.dt/mod.ds)			
	print ('CFL  	Conv 	:	',np.max(mod.umax)*para.dt/mod.ds)	
	print ('CFL  	Capl 	:	',np.max(mod.umax)*para.dt/(mod.ds)**3)		
	if (para.CN==1):
		print ('MAX 	k_CN	:	',mod.k_CN		)				
