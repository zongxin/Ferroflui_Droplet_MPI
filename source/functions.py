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

import parameter as para
import module as mod


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()



Plot_1D_vel 	=False
Plot_quiver_vel =False

def cal_PV(zl,zg,gamma_l,gamma_g,ds):

	# zl=np.array(range(256))
	# zg=np.array(range(256))	
	# gamma_l =np.array(range(256))*0.1
	# gamma_g =np.array(range(256))*0.1

	# current PV contain both 1/j2pi, and z_a term

	Ng 	   		= para.Ng
	dalpha 		= 2*np.pi/Ng
	Nl 			= mod.Nl

	z_e_l 		= zl[::2]
	z_o_l  		= zl[1::2]

	z_e_g 		= zg[::2]
	z_o_g  		= zg[1::2]	

	gamma_e_l 	= gamma_l[::2]
	gamma_o_l  	= gamma_l[1::2]
	
	gamma_e_g 	= gamma_g[::2]
	gamma_o_g  	= gamma_g[1::2]	



	#even_distance
	A,B 				=	np.meshgrid(gamma_o_l,gamma_o_g)
	even_arr,odd_arr 	=	np.meshgrid(z_e_l,z_o_g)
	distance 			= 	even_arr-odd_arr	
	even_matrix  		= 	1./distance*B
	even_int 			= 	np.sum(even_matrix,axis=0)*2.*ds/(2.*np.pi*1j)


	#odd_distance
	A,B 				=	np.meshgrid(gamma_e_l,gamma_e_g)
	odd_arr,even_arr 	=	np.meshgrid(z_o_l,z_e_g)
	distance 			= 	odd_arr-even_arr	
	odd_matrix  		= 	1./distance*B
	odd_int 			= 	np.sum(odd_matrix,axis=0)*2.*ds/(2.*np.pi*1j)	

	



	PV 		 = np.empty(Nl,dtype=complex)	
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
	so_l 		= so[mod.extent_in]
	if (rank==0):
		so_l[:para.nghost_in]=so_l[:para.nghost_in]-L
	if (rank==size-1):
		so_l[-para.nghost_in:]=so_l[-para.nghost_in:]+L	
	# set_trace()
	interpolate_z 		= interp1d(so_l,zl, kind='cubic')
	interpolate_gamma 	= interp1d(so_l,gamma_l, kind='cubic')


	ds 			= L/Ns
	s 			= np.array(range(Ns))*ds
	sn 	 		= s[mod.extent]	
	if (rank==0):
		sn[:para.nghost]=sn[:para.nghost]-L
	if (rank==size-1):
		sn[-para.nghost:]=sn[-para.nghost:]+L	

	# set_trace()
	zl  				= interpolate_z(sn)
	gamma_l				= interpolate_gamma(sn)
	gamma_l 			= gamma_l[para.nghost:-para.nghost]

	recv 				= comm.allgather(zl[para.nghost:-para.nghost])
	mod.zg 				= np.array(recv).flatten()

	recv 				= comm.allgather(gamma_l)
	mod.gamma_g 		= np.array(recv).flatten()

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

	z_s = cal_der(z_face,ds,1,mod.PBC) #dz/d alpha=dxda+idyda
	x_s = np.real(z_s)
	y_s = np.imag(z_s)	

	x   = np.real(z_face)
	y   = np.imag(z_face)
	
	S=0.5*np.sum((x*y_s-y*x_s))*ds

	return S


def cal_der(y,ds,n_dir):

	if   (n_dir == 1):

		dy 			= (y[2:]-y[0:-2])/(2*ds)
		yy 			= np.empty(y.shape,dtype=complex)
		yy[1:-1]	= dy
		yy[0] 		= dy[0] 
		yy[-1] 		= dy[-1]

	elif (n_dir == 2):		

		dy 			= (y[2:]+y[0:-2]-2*y[1:-1])/(ds*ds)
		yy 			= np.empty(y.shape,dtype=complex)
		yy[1:-1]	= dy
		yy[0] 		= dy[0] 
		yy[-1] 		= dy[-1]	

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


def filter_a(zold):
	z_hat=np.fft.fft(zold)
	mask=np.where(np.abs(z_hat)<=1e-3)
	print ("MASK!!!!=",np.shape(mask)	)
	print ("MASK!!!!=",mask)

	z_hat[mask]=0+0j
	znew=np.fft.ifft(z_hat)
	set_trace()
	# znew2=np.fft.ifft(np.fft.fft(zold))
	return znew		

def filter_b(zold):
	N = len(zold)
	x_hat=np.fft.fft(np.real(zold))
	y_hat=np.fft.fft(np.imag(zold))

	x_hat[int(N/2-10):int(N/2+10+1)]=0+0j
	y_hat[int(N/2-10):int(N/2+10+1)]=0+0j	

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