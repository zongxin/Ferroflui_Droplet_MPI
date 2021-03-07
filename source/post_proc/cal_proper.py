import os
import sys
import numpy as np
from pdb import set_trace
from matplotlib import rc as matplotlibrc
from matplotlib.patches import FancyArrowPatch
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from itertools import product, combinations
from scipy import signal
import pickle
import functions as func

matplotlibrc('text.latex', preamble=r'\usepackage{color}')
matplotlibrc('text',usetex=True)
matplotlibrc('font', family='serif')

figheight      = 8
figwidth       = 8
lineWidth      = 2
gcafontSize    = 18
legendFontSize = 20





Nf=184

I=np.zeros(Nf)
C=np.zeros(Nf)
Time=np.zeros(Nf)
for i in range (0,Nf):
	pkl_name	= './data/data2_'+str(i+1)+'.pkl'
	data 		= pickle.load( open( pkl_name, "rb" ) )

	N 		= data['N']	
	t 		= data['time']
	z0 		= data['z0']				
	z 		= data['z']		
	gamma	= data['gamma']	
	kappa	= data['kappa']	

	ds=np.abs(z[2]-z[1])
	z[2]-z[1]
	L = N*ds

	A=func.cal_area(z,ds)
	I[i]=L**2/(4*np.pi*A)
	C[i]=np.max(np.abs(z))/np.min(np.abs(z))
	Time[i]=t

plt.plot(Time,I)
plt.plot(Time,C)
plt.savefig('test.png')
plt.show()



