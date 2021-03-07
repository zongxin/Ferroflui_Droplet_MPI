import os
import sys
import numpy as np
from pdb import set_trace
from matplotlib import rc as matplotlibrc
from matplotlib.patches import FancyArrowPatch
import matplotlib.pyplot as plt
from itertools import product, combinations
from scipy import signal
import pickle

figheight      = 10
figwidth       = 10
lineWidth      = 1.5
textFontSize   = 18
gcafontSize    = 15
legendFontSize = 20

Nfile 	= 500
area 	= np.zeros((Nfile))
area2 	= np.zeros((Nfile))
for i in range (Nfile):
	pkl_name	= './data/data_'+str(i+1)+'.pkl'
	data 		= pickle.load( open( pkl_name, "rb" ) )

	pkl_name2	= '../CN_elliptical_256/data/data_'+str(i+1)+'.pkl'
	data2 		= pickle.load( open( pkl_name, "rb" ) )

	z 	= data['z']	
	z2 	= data['z']		
	if (i==0):	
		z0		= 1.0*z
		N 		= len(z)
		alpha 	= np.linspace(0,2*np.pi,N, endpoint=False)
	if (i==Nfile-1):
		t 		= data['time']	

	r 		= np.abs(z)
	area[i]	= np.sum(r**2*(np.pi/N))

	r2 		= np.abs(z2)
	area2[i]	= np.sum(r2**2*(np.pi/N))	
tlist	= np.linspace(0,t,Nfile)		


Nfile3 	= 144
area3 	= np.zeros((Nfile3))
for i in range (Nfile3):
	pkl_name	= '../RK_elliptical_512/data/data_'+str(i+1)+'.pkl'
	data 		= pickle.load( open( pkl_name, "rb" ) )

	z 	= data['z']	
	
	if (i==0):	
		z0		= 1.0*z
		N 		= len(z)
		alpha 	= np.linspace(0,2*np.pi,N, endpoint=False)
	if (i==Nfile3-1):
		t 		= data['time']	

	r 		= np.abs(z)
	area3[i]	= np.sum(r**2*(np.pi/N))
tlist3	= np.linspace(0,t,Nfile3)



Nfile4 	= 90
area4 	= np.zeros((Nfile4))
for i in range (Nfile4):
	pkl_name	= '../CN_elliptical_512/data/data_'+str(i+1)+'.pkl'
	data 		= pickle.load( open( pkl_name, "rb" ) )

	z 	= data['z']	
	
	if (i==0):	
		z0		= 1.0*z
		N 		= len(z)
		alpha 	= np.linspace(0,2*np.pi,N, endpoint=False)
	if (i==Nfile4-1):
		t 		= data['time']	

	r 		= np.abs(z)
	area4[i]	= np.sum(r**2*(np.pi/N))
tlist4	= np.linspace(0,t,Nfile4)

jj=4
fig = plt.figure(0, figsize=(figwidth*0.8,figheight*0.6))
ax  = fig.add_subplot(111,alpha=1)#
plt.axes(ax)
ax.plot(tlist,area/area[0],'k--',lw=lineWidth,label="RK2 N256")
ax.plot(tlist,area2/area2[0],'g--',lw=lineWidth,label="CN  N256")
ax.plot(tlist3[::jj],area3[::jj]/area3[0],'k-',marker='o',mfc="None",lw=lineWidth,label="RK2 N512")
ax.plot(tlist4[::jj],area4[::jj]/area4[0],'g-',marker='o',mec='g',mfc="None",lw=lineWidth,label="CN  N512")

ax.axhline(1.0,c='r',linestyle='--')
ax.set_xlabel('t',fontsize=1.2*gcafontSize)
ax.set_ylabel('ratio V/V0',fontsize=1.2*gcafontSize)

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels,fontsize=1.0*gcafontSize,loc=0,numpoints=1)
#ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#        ncol=2, mode="expand", borderaxespad=0.)
#ax.set_xticks()


figname= "val_elliptical.png"
plt.tight_layout()		
plt.savefig('./movie/'+figname)
plt.show()
plt.close()
print figname," saved!"
