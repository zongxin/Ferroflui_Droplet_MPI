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

Data_output = 2

Nf= 400
m=40

area=np.zeros(Nf)
dataj=0
for i in range (0,Nf):
	# if (i%Data_output==0):

	dataj=dataj+1

	rpkl_name	= './../data/data2_'+str(i+1)+'.pkl'
	data 		= pickle.load( open( rpkl_name, "rb" ) )

	N 		= data['N']	
	t 		= data['time']
	z0 		= data['z0']				
	z 		= data['z']		

	zlong=np.hstack((z,z[0]))
	dx=np.real(zlong[1:]-zlong[:-1])
	dy=np.imag(zlong[1:]-zlong[:-1])
	ds=(dx**2+dy**2)**0.5

	z_s 	= func.cal_der(z,ds,1)

	xs=np.real(z_s)
	ys=np.imag(z_s)	

	x=np.real(z)
	y=np.imag(z)	

	area[i]=np.sum((x*ys-y*xs)*ds)


set_trace()
plt.plot(area/area[0])
plt.show()



