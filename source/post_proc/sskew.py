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

Nf= 50

Nss=250
Ns=Nss*7
zss_list=np.zeros((7,Nss+1),dtype=complex)
As=np.zeros(Nf)
St=np.zeros(Nf)
time=np.zeros(Nf)
for i in range (0,Nf):

	rpkl_name	= './data/data2_'+str(i+1)+'.pkl'
	data 		= pickle.load( open( rpkl_name, "rb" ) )#,encoding="latin1"

	time[i] 	= data['time']				
	z 		= data['z']	

	zma=np.argmax(np.abs(z))

	zm=np.concatenate((z[zma:],z[:zma]))
	gamma_old=np.zeros(np.shape(z))

	z,_,s,ds=func.get_mesh(zm,gamma_old,Ns)

	for k in range(7):
		if (k<6):
			zss_list[k,:]=z[k*Nss:(k+1)*Nss+1]
		if (k==6)	:
			zss_list[k,:-1]=z[k*Nss:(k+1)*Nss]
			zss_list[k,-1]=z[0]
	

	# zss=z[0:Nss+1]
	fss=np.mean(np.abs(zss_list),axis=0)
	# set_trace()
	angle=np.angle(zss_list[0,:])

	angler=angle[1:]

	if any(angler<angle[:-1]):
		aa=np.where(angler<angle[:-1])
		a1=aa[0][0]+1
		angle[a1:]=angle[a1:]+np.pi*2

	
	angle=angle-angle[0]

	im=np.argmin(np.abs(fss))

	interpolate_a2 = interp1d(np.abs(fss[im:]),angle[im:], kind='cubic')
	interpolate_a1 = interp1d(np.abs(fss[:im])[::-1],angle[:im][::-1], kind='cubic')

	# set_trace()
	x1	  =interpolate_a1([1.0])
	x2	  =interpolate_a2([1.0])
	b1=x1
	b2=np.pi*2/7-x2

	As[i]=b1/b2-1
	a1=np.abs(np.max(fss)-1)
	a2=np.abs(np.min(fss)-1)	
	St[i]=a1/a2-1	
	print i,b1,b2,b1/b2-1
	# plt.scatter(time[i],b1,c='k')
	# plt.scatter(time[i],b2,c='r')	
	# fig = plt.figure(0)	
	# ax  = fig.add_subplot(111,alpha=1)#
	# plt.axes(ax)
	# plt.plot(angle,np.abs(fss))
	# # # set_trace()
	# plt.axvline(x1,c='k')
	# plt.axvline(x2,c='r')
	# plt.axhline(1.0,c='r')	
	# ax.set_xlim([0,0.9])
	# plt.show()
	# plt.close()
plt.plot(time,As)
plt.plot(time,St)

plt.show()

pkl_name='./post_proc/NB_1_37.pkl'
data={}

data['time']	= time
data['As']		= As		
data['St']		= St
data['angle']	= angle
data['fss'] 	= fss			
output = open(pkl_name, 'wb')
pickle.dump(data, output)
output.close()
print pkl_name," saved!"

