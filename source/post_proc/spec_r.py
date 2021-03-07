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
from scipy.signal import argrelextrema
import functions as func
import matplotlib.collections as collections

figheight      = 6
figwidth      = 6
lineWidth      = 1.8
textFontSize   = 18
gcafontSize    = 15
legendFontSize = 20


Nf= 1339
tracer=np.zeros(Nf,dtype='complex')
k4_mag=np.zeros(Nf)	
t_list=np.zeros(Nf)
plotj=0

for i in range (0,Nf):
	if (i%10==0):
		plotj=plotj+1
		pkl_name	= './../data/data2_'+str(i+1)+'.pkl'
		data 		= pickle.load( open( pkl_name, "rb" ) )

		N 		= data['N']	
		t 		= data['time']
		z0 		= data['z0']				
		z 		= data['z']		

		alpha = np.linspace(0,2*np.pi,N, endpoint=False)

		#=======================================================================

		fig = plt.figure(0, figsize=(figwidth*2*0.8,figheight*0.8))


		ax  = fig.add_subplot(121,alpha=1)#
		plt.axes(ax)
		plt.plot(alpha,np.abs(z),'k',linewidth=lineWidth)
		plt.plot(alpha,np.abs(z0),'r--')	
		ax.set_xlabel('alpha',fontsize=1.2*gcafontSize)
		ax.set_ylabel('r',fontsize=1.2*gcafontSize)
		ax.set_ylim([0.97,1.03])
		time_string='t='+str(t)
		ax.text(4,0.975,time_string, fontsize=1.*gcafontSize)

		ax  = fig.add_subplot(122,alpha=1)#
		plt.axes(ax)
		ax.grid(False)
		# a=1.5
		# ax.set_xlim([-a,a])
		# ax.set_ylim([-a,a])	
		m=40
		z_hat=(np.abs(np.fft.fft(np.abs(z)))/256)**2 *2
		if (i==0):
			z_hat0=(np.abs(np.fft.fft(np.abs(z0)))/256)**2 *2
		plt.plot(np.array(range(1,m)),z_hat[1:m],'k',marker='.',linewidth=lineWidth)
		plt.plot(np.array(range(1,m)),z_hat0[1:m],'r--',marker='.')	
		ax.set_xlabel('k',fontsize=1.2*gcafontSize)
		ax.set_ylabel('E(k,t)',fontsize=1.2*gcafontSize)		

		ax.set_ylim([0,10e-7])


		figname="ttttest_"+str(plotj)+'.png'
		plt.tight_layout()		
		plt.savefig('./../tttest/'+figname)
		# plt.show()
		plt.close()
		print figname," saved!"
	# set_trace()