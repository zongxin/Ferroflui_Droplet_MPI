import os
import sys
import numpy as np
from pdb import set_trace
from matplotlib import rc as matplotlibrc
from matplotlib.patches import FancyArrowPatch
import matplotlib.pyplot as plt
from itertools import product, combinations
from scipy import signal
import functions as func
from mpi4py import MPI

import parameter as para
import module as mod




def cal_vel(zl_lo,gamma_0):


	comm = MPI.COMM_WORLD

	# zl_lo 		:local z+ghost cell
	# gamma_0 		:local z without ghost cell
	# This function include the gather of gamma
	# Outputs has no ghost cell

	tau			= para.tau
	NB1			= para.NB1
	NB2			= para.NB2
	chi			= para.chi
	NORM_STRESS	= para.NORM_STRESS
	ngost		= para.nghost

	ds			= mod.ds
	b 			= mod.b
	b_dot		= mod.b_dot
	Err_out 	= mod.Err_out



	# gamma_0 : initial guess
	z_s 		= func.cal_der(zl_lo,ds,1)

	kappa 		= func.cal_kappa(zl_lo,ds)
	kappa_s 	= np.real(func.cal_der(kappa,ds,1))

	phi_s,psi_s,phi_norm_s,psi_norm_s,phi_psi_s=cal_mag(zl_lo,z_s,ds)

	R_const 	= b*b*(tau*kappa_s - NB1/2*phi_s- NB2/2*psi_s)\
					-b_dot/(2.*b)*np.real(z_s*np.conj(zl_lo))

#	REMOVE GHOST CELL	
	zl_sh		= zl_lo[ngost:-ngost]	
	z_s			= z_s[ngost:-ngost]		
	R_const		= R_const[ngost:-ngost]				
	phi_norm_s	= phi_norm_s[ngost:-ngost]
	psi_norm_s	= psi_norm_s[ngost:-ngost]
	phi_psi_s	= phi_psi_s[ngost:-ngost]
	kappa 		= kappa [ngost:-ngost]

	if (NORM_STRESS):
		# set_trace()
		R_const 	= R_const-b*b*chi/2*(NB1*phi_norm_s\
										+ NB2*psi_norm_s\
										-(np.sqrt(NB1*NB2))*phi_psi_s)

	gamma_new 	  = gamma_0*1.0

	PV_new 	  = func.cal_PV(zl_sh,mod.zg,gamma_new,mod.gamma_g,ds)

	subiter	  = 0
	it_err    = 1.0
	mod.pv_err= 1.0
	while ((mod.pv_err>1e-6)&(subiter<100)):
		

		gamma 	  		= gamma_new*1.0
		recv 			= comm.allgather(gamma)
		mod.gamma_g 	= np.array(recv).flatten()


		PV 		  = PV_new*1.0
		PV_new 	  = func.cal_PV(zl_sh,mod.zg,gamma,mod.gamma_g,ds)
		


		gamma_new  = (R_const+np.real(PV_new*z_s))*2

		it_err	  = np.mean((np.abs(PV_new-PV)/np.abs(PV))**2)**0.5

		if (subiter == 0):
			it_err	=  1.0
		subiter     =  subiter+1

		mod.pv_err=comm.allreduce(it_err, op=MPI.MAX)

	gamma  = gamma - np.mean(mod.gamma_g)
	mod.gamma_g 	= mod.gamma_g - np.mean(mod.gamma_g)

	u_conj = -b_dot/(2*b)*np.conj(zl_sh)-gamma/(2*z_s)+func.cal_PV(zl_sh,mod.zg,gamma,mod.gamma_g,ds)

	# print ('stop for var=='	)
	# set_trace()
	if (Err_out==1):
		print ("Subiterations: ",subiter,"	Error:",it_err)#,"	N=",len(z)

	return np.conj(u_conj),gamma,kappa


def cal_mag(z,z_s,ds):

	phi 		= 1./np.abs(z)**2
	phi_s 		= np.real(func.cal_der(phi,ds,1))
	phi_norm 	= (np.real(z_s*np.conj(z))/np.abs(z)**2)**2
	phi_norm_s 	= np.real(func.cal_der(phi_norm,ds,1))	

	psi 		= np.abs(z)**2
	psi_s 		= np.real(func.cal_der(psi,ds,1))
	psi_norm 	= (np.imag(z_s*np.conj(z)))**2
	psi_norm_s 	= np.real(func.cal_der(psi_norm,ds,1))	

	phi_psi 	= (np.imag(z**2)*(-np.real(z_s**2))+np.imag(z_s**2)*np.real(z**2))/np.abs(z)**2
	phi_psi_s	= np.real(func.cal_der(phi_psi,ds,1))	

	return phi_s,psi_s,phi_norm_s,psi_norm_s,phi_psi_s


