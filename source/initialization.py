import os
import sys
import numpy as np
from pdb import set_trace
import module as mod
import parameter as para
from mpi4py import MPI
comm = MPI.COMM_WORLD
import pickle

rank = comm.Get_rank()
size = comm.Get_size()


def ini_shape():

	if (para.Restart==1):
		data 			= pickle.load( open( para.Resname, "rb" ) )

		mod.t 			= data['time']			
		mod.zg 			= data['z']		
		mod.z0 			= data['z0']			
		mod.gamma		= data['gamma']	

		mod.dataj		= para.Res_i
		mod.plotj		= para.Res_i
		mod.Ng 			= len(mod.zg)
	else:	
		alpha 	=	para.alpha
		R 		=	para.R
		n 		=	para.n

		# xn 		=	(-0.052548103890057304-0.004128940307543485j)
		# xn2 	=	(0.017920232977232567+0.002658398993006049j)
		# xn3 	=	(-0.008535189439188302-0.0016199392387838797j)
		# xn4 	=	(0.004667277501091407+0.0009935965353182817j)
		# xn5 	=	(-0.002735911613585555-0.0006034506962295917j)

		xn 				=	1./500
		xn2,xn3,xn4,xn5 =	0,0,0,0

		x0 		=	-(np.abs(xn)**2+np.abs(xn2)**2+np.abs(xn3)**2+np.abs(xn4)**2+np.abs(xn5)**2)/R
			
		xi		=	x0+ xn*np.exp(1j*n*alpha)+np.conj(xn)*np.exp(-1j*n*alpha)\
					+xn2*np.exp(1j*n*2*alpha)+np.conj(xn2)*np.exp(-1j*n*2*alpha)\
					+xn3*np.exp(1j*n*3*alpha)+np.conj(xn3)*np.exp(-1j*n*3*alpha)\
					+xn4*np.exp(1j*n*4*alpha)+np.conj(xn4)*np.exp(-1j*n*4*alpha)\
					+xn5*np.exp(1j*n*5*alpha)+np.conj(xn5)*np.exp(-1j*n*5*alpha)
	
		RR 		=	R+xi
		mod.z0 	=	RR*np.cos(alpha)+1j*RR*np.sin(alpha)
		mod.zg 	=	mod.z0*1.0

		#initialize gamma
		gamma_0    = 0.1
		mod.gamma  = gamma_0*np.ones(para.Ng)


	if (rank==0):
		print(' ')
		print('--->','PURTURBATION INITIALIZATION FINISHED!')

	return mod.extent_in


def ini_mpi(size,rank):	

	mod.Nl 			= int(para.Ng/size)	
	mod.extent 		= get_MPI_extent(para.Ng,rank,size)
	mod.extent_in 	= get_MPI_extent_int(para.Ng,rank,size)

	if (rank==0):
		print(' ')
		print('--->','MPI INITIALIZATION FINISHED!')

	return mod.extent_in

def get_MPI_extent(Nglobal, myrank, totalsize ):

	Nl = int(Nglobal/totalsize)
	ngost=para.nghost
	G_list=np.arange(0, Nglobal, dtype=int)


	if (totalsize==1):
		extent = np.concatenate([G_list[-ngost:],G_list,G_list[:ngost]])

	else:	
		if myrank == 0:
			extent = np.concatenate([G_list[-ngost:],G_list[:Nl*(myrank+1)+ngost]])

		elif myrank == (totalsize-1):
			extent = np.concatenate([G_list[Nl*myrank-ngost:],G_list[:ngost]])

		else:
			extent = G_list[Nl*(myrank)-ngost:Nl*(myrank+1)+ngost]

	return extent

def get_MPI_extent_int(Nglobal, myrank, totalsize ):



	Nl = int(Nglobal/totalsize)
	ngost=para.nghost_in
	G_list=np.arange(0, Nglobal, dtype=int)

	if (totalsize==1):
		extent = np.concatenate([G_list[-ngost:],G_list,G_list[:ngost]])

	else:
		if myrank == 0:
			extent = np.concatenate([G_list[-ngost:],G_list[:Nl*(myrank+1)+ngost]])

		elif myrank == (totalsize-1):
			extent = np.concatenate([G_list[Nl*myrank-ngost:],G_list[:ngost]])

		else:
			extent = G_list[Nl*(myrank)-ngost:Nl*(myrank+1)+ngost]

	return extent	


