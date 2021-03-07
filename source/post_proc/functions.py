import os
import sys
import numpy as np
from pdb import set_trace
from matplotlib import rc as matplotlibrc
from matplotlib.patches import FancyArrowPatch
import matplotlib.pyplot as plt
from itertools import product, combinations
from scipy import signal
from scipy.interpolate import interp1d
from scipy.interpolate import Rbf
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# from mpl_toolkits.mplot3d import proj3d
# from scipy.interpolate import RegularGridInterpolator
# from mpl_toolkits.mplot3d import Axes3D


Plot_1D_vel 	=False
Plot_quiver_vel =False

def cal_PV(z_face,z_alpha,gamma,ds):
	# current PV contain both 1/j2pi, and z_a term
	N 	   = len(z_face)
	dalpha = 2*np.pi/N

	even_list = z_face[::2]
	odd_list  = z_face[1::2]

	even_gamma = gamma[::2]
	odd_gamma  = gamma[1::2]

	odd_arr,even_arr 	= np.meshgrid(odd_list,even_list)
	even_gamma_arrA,even_gamma_arrB    = np.meshgrid(even_gamma,even_gamma)	
	odd_gamma_arrA, odd_gamma_arrB     = np.meshgrid(odd_gamma,odd_gamma)


	#even_distance
	distance 	= even_arr-odd_arr	
	even_matrix = 1./distance*odd_gamma_arrA
	even_int 	= np.sum(even_matrix,axis=1)*2.*ds/(2.*np.pi*1j)


	#odd_distance
	distance 	= odd_arr-even_arr

	odd_matrix  = 1./distance*even_gamma_arrB
	odd_int 	= np.sum(odd_matrix,axis=0)*2.*ds/(2.*np.pi*1j)	


	PV 		 = np.empty(z_face.shape,dtype=complex)	
	PV[::2]  = even_int
	PV[1::2] = odd_int

	return PV

def get_mesh(zl,gamma_l,Ns):
	



	zlong		= np.concatenate((mod.zg,mod.zg[0:1]))
	gamma_long	= np.concatenate((mod.gamma_g,mod.gamma_g[0:1]))	
	dz			= zlong[1:]-zlong[:-1]
	ds 			= np.abs(dz)

	s_old		= np.cumsum(ds)
	s_old		= np.concatenate(([0.0],s_old))
	so 			= s_old[:-1]
	L 			= np.sum(ds)
	so_l 		= so[mod.extent_int]

	interpolate_z 		= interp1d(so_l,zl, kind='cubic')
	interpolate_gamma 	= interp1d(so_l,gamma1, kind='cubic')

	ds 			= L/Ns
	s 			= np.array(range(Ns))*ds
	sn 	 		= s[mod.extent]	

	zl  				= interpolate_z(sn)
	gamma_l				= interpolate_gamma(sn)


	# if (size>1):

	# else:
	# 	zlong=np.concatenate((z,z[0:1]))
	# 	gamma_long=np.concatenate((gamma,gamma[0:1]))	
	# 	dz=zlong[1:]-zlong[:-1]
	# 	ds=np.abs(dz)
	# 	s_old=np.cumsum(ds)
	# 	s_old=np.concatenate(([0.0],s_old))

	# 	interpolate_z = interp1d(s_old,zlong, kind='cubic')
	# 	interpolate_gamma = interp1d(s_old,gamma_long, kind='cubic')

	# 	L=np.sum(ds)
	# 	ds=L/Ns
	# 	s=np.array(range(Ns))*ds

	# 	mesh_z	  =interpolate_z(s)
	# 	mesh_gamma=interpolate_gamma(s)

	return zl,gamma_l,s,ds

def cal_area(z_face,ds):

	#A=0.5 int_c (x.ys-y.xs)ds

	z_s = cal_der(z_face,ds,1) #dz/d alpha=dxda+idyda
	x_s = np.real(z_s)
	y_s = np.imag(z_s)	

	x   = np.real(z_face)
	y   = np.imag(z_face)
	
	S=0.5*np.sum((x*y_s-y*x_s))*ds

	return S


def cal_der(y,ds,n_dir):

	if   (n_dir == 1):
		ylong 	= np.concatenate((y,y[0:2]))
		dy 		= (ylong[2:]-ylong[0:-2])/(2*ds)
		yy 		= np.empty(y.shape,dtype=complex)
		yy[1:]	= dy[:-1]
		yy[0] 	= dy[-1] 

	elif (n_dir == 2):		
		ylong 	= np.concatenate((y,y[0:2]))
		dy 		= (ylong[2:]+ylong[0:-2]-2*ylong[1:-1])/(ds*ds)
		yy 		= np.empty(y.shape,dtype=complex)
		yy[1:] 	= dy[:-1]
		yy[0] 	= dy[-1]				
	else:
		yy 		= "None"

	return yy		

def cal_kappa(z_face,ds):
	z_a = cal_der(z_face,ds,1) #dz/d alpha=dxda+idyda
	x_a = np.real(z_a)
	y_a = np.imag(z_a)	
	s_a = np.abs(z_a) # ds/ d alpha= (x_a^2+y_a^2)^0.5

	z_aa = cal_der(z_face,ds,2) #dz/d alpha=dxda+idyda
	x_aa = np.real(z_aa)
	y_aa = np.imag(z_aa)	

	kappa=(x_a*y_aa-x_aa*y_a)/s_a**3

	return kappa	


# def filter_a(zold):
# 	z_hat=np.fft.fft(zold)
# 	mask=np.where(np.abs(z_hat)<=1e-3)
# 	print "MASK!!!!=",np.shape(mask)	
# 	print "MASK!!!!=",mask

# 	z_hat[mask]=0+0j
# 	znew=np.fft.ifft(z_hat)
# 	set_trace()
# 	# znew2=np.fft.ifft(np.fft.fft(zold))
# 	return znew		

def filter_b(zold):
	N = len(zold)
	x_hat=np.fft.fft(np.real(zold))
	y_hat=np.fft.fft(np.imag(zold))

	x_hat[N/2-10:N/2+10+1]=0+0j
	y_hat[N/2-10:N/2+10+1]=0+0j	

	xnew=np.real(np.fft.ifft(x_hat))
	ynew=np.real(np.fft.ifft(y_hat))	
	znew = xnew+1j*ynew
	# set_trace()
	return znew	

def get_tracer(z,z2,tracer0):

	N=len(z)
	NN=np.int(N/4)
	# if (time==0):
	xt0=np.real(tracer0)
	yt0=np.imag(tracer0)


	x 	= np.real(z)
	y 	= np.imag(z)
	xn 	= np.real(z2)
	yn 	= np.imag(z2)

	rbfx = Rbf(x, y, xn)#, epsilon=2
	rbfy = Rbf(x, y, yn)

	xt  = rbfx(xt0, yt0)
	yt  = rbfy(xt0, yt0)

	alpha=yn/xn
	alpha=alpha[:NN]
	theta=yt/xt

	interpolate_x = interp1d(alpha,xn[:NN], kind='cubic')
	xt	  =interpolate_x(theta)
	interpolate_y = interp1d(alpha,yn[:NN], kind='cubic')
	yt	  =interpolate_y(theta)	
	# plt.close()
	# plt.plot(x,y,'r')
	# plt.plot(xn,yn,'k')
	# plt.scatter(xt0,yt0,c='r',marker='+')
	# plt.scatter(xt,yt,c='k',marker='+')	
	# plt.savefig('test.png')
	return xt+1j*yt


def cal_phase(z):


	z_hat=np.fft.fft(np.abs(z))
	z_ang=np.arctan(np.imag(z_hat)/np.real(z_hat))

	x = np.real(z_hat)
	y = np.imag(z_hat)
	
	m1=np.where((x<0)&(y>0))
	z_ang[m1]=z_ang[m1]+np.pi

	m2=np.where((x<0)&(y<0))
	z_ang[m2]=z_ang[m2]-np.pi
	# set_trace()

	mask = np.where(np.abs(z_hat)<1e-6)
	z_ang[mask]=0
	return z_ang

# def filter_b(zold):
# 	N = len(zold)
# 	z_hat=np.fft.fft(zold)

# 	z_hat[N/2-10:N/2+10+1]=0+0j
# 	znew=np.fft.ifft(z_hat)
# 	# set_trace()
# 	return znew	