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
from scipy.interpolate import interp1d
from mpi4py import MPI

import parameter as para
import module as mod
import initialization as init
import out_put as wrt

import functions as func
import time_advance as adv


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()



Ng	=	para.Ng
dt 	=	para.dt

init.ini_shape()

extent_int=init.ini_mpi(size,rank)




# print('myrank:',rank,'myextent:',mod.extent)


t 	=	mod.t

#	Start time advance 
if (rank==0):
	print(' ')
	print('--->','TIME ADVANCE STARTED!')

for i in range (5000000):


	zl_lo 		=	mod.zg[extent_int]
	gamma_l_lo	=	mod.gamma_g[extent_int]


	zl_lo,gamma_l_sh,s,mod.ds=func.get_mesh(zl_lo,gamma_l_lo,Ng)

	adv.advance(zl_lo,gamma_l_sh,t,dt,mod.ds)



	if (rank==0):
		if (i%para.Plot_output==0):

			wrt.write_figure(t)

		if (i%para.Data_output==0):
		
			wrt.write_data(t)

		if (i%para.Text_output==0):		

			wrt.print_screen(i,t)

	i=i+1	
	t=t+dt	
