import os
import sys
import numpy as np
from pdb import set_trace
from matplotlib import rc as matplotlibrc
from matplotlib.patches import FancyArrowPatch
import matplotlib.pyplot as plt
from itertools import product, combinations
from scipy import signal
# from scipy.integrate import solve_bvp
# from scipy.integrate import solve_ivp
import pickle
from scipy.interpolate import interp1d
import functions as func
from scipy.signal import hilbert

# matplotlibrc('text.latex', preamble=r'\usepackage{color}')
matplotlibrc('text',usetex=True)
matplotlibrc('font', family='serif')

figheight      = 6
figwidth      = 6
lineWidth      = 1.5
textFontSize   = 18
gcafontSize    = 15
legendFontSize = 20


def cal_mean(f,da):
	return np.sum(f)*da/(2*np.pi)


Nf= 245
nname='k7_NB1_1'

As=np.zeros(Nf)
Sk=np.zeros(Nf)
time=np.zeros(Nf)

for i in range (0,Nf):

	rpkl_name	= './../data/data2_'+str(i+1)+'.pkl'
	data 		= pickle.load( open( rpkl_name, "rb" ) )#,encoding="latin1"

	time[i] 	= data['time']				
	z 		= data['z']	


	zlong=np.concatenate((z,z[0:1]))

	ang=np.angle(z)
	ang2=np.hstack((ang[-1],ang[0:-1]))
	a=np.where((ang-ang2)<0)[0]
	a=np.int(a[0]*1.0)

	# ang_end=np.abs(ang[0]-ang[1])

	ang3=ang*1.0
	ang3[a:]=ang3[a:]+2*np.pi
	ang3=np.hstack((ang3,ang3[0]+2*np.pi))
	ang_int=np.linspace(ang3[0],ang3[0]+2*np.pi,512,endpoint=False)
	# set_trace()
	interpolate_z = interp1d(ang3,zlong, kind='cubic')		
	new_z	 	 =interpolate_z(ang_int)

	da=ang_int[1]-ang_int[0]
	fss=np.abs(new_z)
	fm=cal_mean(fss,da)
	f_p=fss-fm
	Sk[i]=cal_mean(f_p**3,da)/cal_mean(f_p**2,da)**1.5
	hilb=np.imag(hilbert(f_p))
	As[i]=cal_mean(hilb**3,da)/cal_mean(f_p**2,da)**1.5

	# set_trace()
fig = plt.figure(0,figsize=(figwidth*1.0,figheight*0.8))	
ax  = fig.add_subplot(111,alpha=1)#
plt.plot(time,As)
plt.plot(time,Sk)
time_string=r'$A_s=$'+str(As[-2])
ax.text(0.01,0.7,time_string)
time_string=r'$S_k=$'+str(Sk[-2])
ax.text(0.01,0.8,time_string)
ax.set_xlabel(r'$t$',fontsize=1.2*gcafontSize)
ax.set_ylabel(r'$A_s,S_k$',fontsize=1.2*gcafontSize)
figname=nname+'.png'
plt.tight_layout()		
plt.savefig(figname)
plt.show()
plt.close()
print (figname," saved!")


pkl_name=nname+'.pkl'
data={}

data['time']	= time
data['As']		= As		
data['Sk']		= Sk			
output = open(pkl_name, 'wb')
pickle.dump(data, output)
output.close()
print (pkl_name," saved!")

