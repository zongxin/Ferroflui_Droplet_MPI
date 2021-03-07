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
from matplotlib import rcParams, cycler

figheight      = 5
figwidth       = 4*2
lineWidth      = 1.5
textFontSize   = 18
gcafontSize    = 15
legendFontSize = 20

Nf= 1
m =40

Energy=np.zeros((m,Nf))
tlist =np.zeros(Nf)

for i in range (0,Nf):

	rpkl_name	= './fft_data/fft_'+str(i+1)+'.pkl'
	data 		= pickle.load( open( rpkl_name, "rb" ),encoding='latin1' )
	tlist[i] 		= data['time']			
	z_hat 	= np.abs(data['z_hat'])
	z_hat 	= z_hat	/512
	Energy[:,i] = z_hat
	Energy[0,i] = z_hat[0]

set_trace()
cmap = plt.cm.viridis_r
rcParams['axes.prop_cycle'] = cycler(color=cmap(np.linspace(0, 1, 15)))

fig = plt.figure(0, figsize=(figwidth,figheight))

ax  = fig.add_subplot(111,alpha=1)#


for i in range (1,16):
	plt.plot(tlist,Energy[i,:],lw=lineWidth,label='k='+str(i))
	# set_trace()
ax.set_xlabel('time',fontsize=1.2*gcafontSize)
ax.set_ylabel('E(k,t)',fontsize=1.2*gcafontSize)

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels,fontsize=0.7*gcafontSize,bbox_to_anchor=(1.04,1), loc="upper left",numpoints=1)

# ax  = fig.add_subplot(212,alpha=1)#
# for i in range (1,10):
# 	plt.plot(Energy[i,:],lw=lineWidth)
# ax.set_ylim([0,3e-5])
# ax.set_xlabel('time',fontsize=1.2*gcafontSize)
# ax.set_ylabel('E(k,t)',fontsize=1.2*gcafontSize)

figname='Energy.png'
plt.tight_layout()		
plt.savefig(figname)
plt.show()
plt.close()
print (figname," saved!")


