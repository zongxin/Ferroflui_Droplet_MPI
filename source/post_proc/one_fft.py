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
from scipy.interpolate import interp1d


Nf= 179
m=40
dataj=0

fftdata=np.zeros((Nf,m),dtype='complex')
time   = np.zeros((Nf))

for i in range (0,Nf):


	dataj=dataj+1

	rpkl_name	= './../data/data2_'+str(i+1)+'.pkl'
	data 		= pickle.load( open( rpkl_name, "rb" ) )#,encoding="latin1"

	t 	= data['time']				
	z 	= data['z']	


	zlong=np.concatenate((z,z[0:1]))

	ang=np.angle(z)
	ang2=np.hstack((ang[-1],ang[0:-1]))
	a=np.where((ang-ang2)<0)[0]
	a=np.int(a[0]*1.0)


	ang3=ang*1.0
	ang3[a:]=ang3[a:]+2*np.pi

	ang_int=np.linspace(0.,2*np.pi,512,endpoint=False)
	ang3long=np.concatenate((ang3[-400:]-2*np.pi,ang3,ang3[:80]+2*np.pi))
	zlong=np.concatenate((z[-400:],z,z[:80]))	

	interpolate_z = interp1d(ang3long,zlong, kind='cubic')		
	new_z	 	 =interpolate_z(ang_int)
	z_hat=np.fft.fft(np.abs(new_z))/len(new_z)
	time[i]=t
	fftdata[i,:]=z_hat[0:m]

wpkl_name='fft_one.pkl'
data={}
data['time']	= time
data['z_hat']	= fftdata			
output = open(wpkl_name, 'wb')
pickle.dump(data, output)
output.close()
print (wpkl_name," saved!")



