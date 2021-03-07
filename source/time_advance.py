import os
import sys
import numpy as np
from pdb import set_trace
from matplotlib import rc as matplotlibrc
from matplotlib.patches import FancyArrowPatch
import matplotlib.pyplot as plt
# plt.switch_backend('agg')
# plt.switch_backend('TkAgg')
from itertools import product, combinations
from scipy import signal
import functions as func
import flowvar as flow
from mpi4py import MPI

import parameter as para
import module as mod

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


def Crank_Nicolson(zl_l,ds,dt,gamma_old):
#	BRO!  BE CAREFUL OF CONGUGATE!

	dt2		= dt*0.5
	alpha 	= 0.#dt2*b_dot/(2*b)
	a_new 	= 0.#dt2*b_dot1/(2*b1)
	zk0_c	= np.conj(zl_l)
	gamma0	= 1.0 * gamma_old
	z_s0 	= func.cal_der(zl_l,ds,1)

	ngost		= para.nghost
	zl_sh		= zl_l[ngost:-ngost]
	zk0_c		= zk0_c[ngost:-ngost]	
	z_s0		= z_s0[ngost:-ngost]

	PV0			= func.cal_PV(zl_sh,mod.zg,gamma0,mod.gamma_g,ds)	

	k 		= 0
	it_err	= 1.0
	z_sk 	= 1.0 * z_s0
	gamma_k = 1.0 * gamma0
	PV_k	= 1.0 * PV0
	zknew 	= 1.0 * zk0_c

	while ((it_err>1e-6)&(k<=100)):

		zk_old		= 1.0*zknew
		gamma_store	= 1.0*gamma_k

		#zknew		= (1-alpha)/(1+a_new)*zk0_c\
		zknew		= zk0_c- (dt2/(1+alpha))*(gamma_k/(2*z_sk)+gamma0/(2*z_s0))\
					  + dt2/(1+alpha)*(PV_k+PV0)

		recv 		= comm.allgather(zknew)
		mod.zg 		= np.array(recv).flatten()

		zl_k_lo		= mod.zg[mod.extent] #long
		zl_k_sh 	= zl_k_lo[ngost:-ngost]
		zknew 		= zl_k_sh*1.0

		z_sk 		= func.cal_der(np.conj(zl_k_lo),ds,1)	
		z_sk 		= z_sk[ngost:-ngost]

		mod.zg		=	np.conj(mod.zg)
		uu,gamma_k,kappa 	= flow.cal_vel(np.conj(zl_k_lo),gamma_store)
		PV_k				= func.cal_PV (np.conj(zl_k_sh),mod.zg,gamma_k,mod.gamma_g,ds)

		k 	   	= k+1
		it_err 	= np.mean(np.abs(zk_old-zknew)**2)**0.5
		# print ('stop for CN=='	)
		# set_trace()
	# if (rank==0):
	# 	plt.plot(1/(1+alpha)*(PV_k+PV0))
	# 	# set_trace()
	# 	# plt.plot(mod.gamma_g)
	# 	plt.grid('on')	
	# 	plt.savefig('test_cn.png')
	# exit()


	mod.umax	=	comm.reduce(np.max(np.abs(uu)), op=MPI.MAX, root=0)
	mod.k_CN	=	k

	return kappa




def Runge_Kutta(zl_l,ds,dt,gamma_old):
	ngost		= para.nghost
	if (para.RK==1):

		u1,gamma,kappa	= flow.cal_vel(zl_l,gamma_old)



		gamma_old		= gamma*1.0
		z1				= zl_l+u1*dt	
		z = 1.0*z1
		# zout=z1*1.0	
		# z1				= func.filter_b(z1)

	if (para.RK==2):

		zl_sh		= zl_l[ngost:-ngost]

		# Stage 1
		# set_trace()
		u1,gamma,kappa	= flow.cal_vel(zl_l,gamma_old)
		gamma_old		= gamma*1.0
		# set_trace()
		z1				= zl_sh+u1*dt	

		recv 		= comm.allgather(z1)
		mod.zg 		= np.array(recv).flatten()	
		mod.zg		= func.filter_b(mod.zg)
		z1			= mod.zg[mod.extent] #long


		u2,gamma,kappa	= flow.cal_vel(z1,gamma_old)
		gamma_old		= gamma*1.0
		z2	= zl_sh+(u1+u2)*dt/2	

		recv 		= comm.allgather(z2)
		mod.zg 		= np.array(recv).flatten()	
		mod.zg		= func.filter_b(mod.zg)
		z3			= mod.zg[mod.extent] #long

		uu=u1+u2

		mod.umax=comm.reduce(np.max(np.abs(uu)), op=MPI.MAX, root=0)

		_,gamma,kappa	= flow.cal_vel(z3,gamma_old)	


	if (para.RK==4):

		# Stage 1
		u1,gamma1,kappa	= flow.cal_vel(z,gamma_old)
		z1				= z+u1*dt*0.5

		# Stage 2
		u2,gamma2,kappa	= flow.cal_vel(z1,gamma_old)
		z2				= z+u2*dt*0.5

		# Stage 3
		u3,gamma3,kappa	= flow.cal_vel(z2,gamma_old)
		z3				= z+u3*dt

		# Stage 4
		u4, gamma4,kappa= flow.cal_vel(z3,gamma_old)


		z 				= z+(u1/6+u2/3+u3/3+u4/6)*dt
		_, gamma,kappa  = flow.cal_vel(z,gamma_old)	
		# z				= func.filter_b(zz)	

	return kappa




def Adams_Bashforth(z,z_old,ds,dt,gamma_old):
	u1,gamma,kappa		= flow.cal_vel(z,gamma_old,)
	u_old,gamma,kappa	= flow.cal_vel(z_old,gamma_old)	
	z					= z+0.5*dt*(3*u1-u_old)		
	_,gamma,kappa		= flow.cal_vel(z,gamma_old)	

	return z,gamma,kappa




def advance(zl_lo,gamma_old,t,dt,ds):

	if (para.CN == 1):
		if (t==0):

			_,gamma_old,_	= flow.cal_vel(zl_lo,gamma_old)

		kappa				= Crank_Nicolson(zl_lo,ds,dt,gamma_old)


		recv 				= comm.gather(kappa, root=0)
		mod.kappa 			= np.array(recv).flatten()


	if (para.RK > 0):

		kappa	= Runge_Kutta(zl_lo,ds,dt,gamma_old)

		recv 				= comm.gather(kappa, root=0)
		mod.kappa 			= np.array(recv).flatten()		
		# tracer 				= func.get_tracer(z,z2,tracer_old)

		# tracer 				= tracer_old 
		# tracer_old 			= tracer*1.0

	if (para.AB == 1):
		# parameter for AB
		# b 	 	  =	1.0
		# b1 	      =	1.0	
		# b_old 	  =	1.0		
		# b_dot 	  =	0.0
		# b_dot1    =	0.0
		# b_dot_old =	0.0

		if (t==0):	
			z_n 				= z*1.0
			z_new,gamma,mod.kappa	= Runge_Kutta(z_n,ds,dt)	
			z 					= z_new*1.0	
		else:
			z_n_1				= z_n  *1.0
			z_n 				= z_new*1.0
			z_new,gamma,mod.kappa	= Adams_Bashforth(z_n,z_n_1,ds,dt,gamma_old)
			z 					= z_new*1.0		
