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

figheight      = 4
figwidth       = 5
lineWidth      = 1.5
textFontSize   = 18
gcafontSize    = 15
legendFontSize = 20

nname='vf_NB1_1'
Nf= 179
Ne=40
phase=np.zeros((5,Nf))
w_list=np.zeros((5,Nf))
vf_list=np.zeros((5,Nf))
Sk=np.zeros(Nf)
time=np.zeros(Nf)
m=np.array(range(5))


rpkl_name	= 'fft_one.pkl'
data 		= pickle.load( open( rpkl_name, "rb" ) )#,encoding="latin1"

time     	= data['time']				
z_hat 		= data['z_hat']	

for i in range (0,Nf):
	phase[m,i]=np.angle(z_hat[i,(m+1)*7])

for im in range(5):
	phase1=phase[im,0:-1]
	phase2=phase[im,1:]
	a=np.where((phase1-phase2)<0)[0]
	# set_trace()
	for ia,aa in enumerate(a):
		phase[im,aa+1:]=phase[im,aa+1:]-2*np.pi
	w_list[im,:]=-np.gradient(phase[im,:],time)
	vf_list[im,:]=w_list[im,:]/((im+1)*7) 


fig = plt.figure(0,figsize=(figwidth*2.0,figheight*0.8))

ax  = fig.add_subplot(121,alpha=1)#
for im in range(5):
	plt.plot(time,w_list[im,:])

ax.set_ylim([0,5000])
ax.set_xlabel(r'$t$',fontsize=1.2*gcafontSize)
ax.set_ylabel(r'$d\omega/dt$',fontsize=1.2*gcafontSize)
plt.setp(ax.get_xticklabels(),fontsize=1.2*gcafontSize)
plt.setp(ax.get_yticklabels(), fontsize=1.2*gcafontSize)
ax.set_xticks([0,0.01,0.02,0.03])

ax  = fig.add_subplot(122,alpha=1)#
for im in range(5):
	plt.plot(time,vf_list[im,:])

wpkl_name=nname+'.pkl'
data={}
data['time']	= time
data['z_hat']	= vf_list			
output = open(wpkl_name, 'wb')
pickle.dump(data, output)
output.close()
print (wpkl_name," saved!")


time_string=r'$v_f=$'+str(np.mean(vf_list[:,-3]))
ax.text(0.01,90,time_string)
ax.set_ylim([60,100])
ax.set_xlabel(r'$t$',fontsize=1.2*gcafontSize)
ax.set_ylabel(r'$v_f$',fontsize=1.2*gcafontSize)
plt.setp(ax.get_xticklabels(),fontsize=1.2*gcafontSize)
plt.setp(ax.get_yticklabels(), fontsize=1.2*gcafontSize)
ax.set_xticks([0,0.01,0.02,0.03])
figname=nname+'.png'
plt.tight_layout()		
plt.savefig(figname)
plt.show()
plt.close()
print (figname," saved!")





# plt.plot(time,As)
# plt.plot(time,Sk)

# plt.show()

# pkl_name='./post_proc/H_NB_1_37.pkl'
# data={}

# data['time']	= time
# data['As']		= As		
# data['Sk']		= Sk			
# output = open(pkl_name, 'wb')
# pickle.dump(data, output)
# output.close()
# print pkl_name," saved!"

