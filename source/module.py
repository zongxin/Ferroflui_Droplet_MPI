import os
import sys
import numpy as np
from pdb import set_trace
import parameter as para


Ng			=	para.Ng
Nl			=	Ng*1

t 			=	0.
b 	    	=	np.exp(0.)
b_dot   	=	0.
ds 			=	0.


z0			=	np.zeros(Ng,dtype=complex)
zg			=	np.zeros(Ng,dtype=complex)
zkg			=	np.zeros(Ng,dtype=complex)

gamma0		=	np.zeros(Ng)
gamma_g		=	np.zeros(Ng)
gammak_g	=	np.zeros(Ng)
kappa		=	np.zeros(Ng)

#	peroidic BC
PBC 		=	1

#	MPI extent
extent		=	np.zeros(Ng)
extent_in	=	np.zeros(Ng)

#	Iteration info
k_CN 		= 0
umax 		= 0.
it_err		= 0.
pv_err		= 0.

#	Out put control
plotj		= 1
dataj		= 1
Err_out 	= 0

#		Tracer
tracer0 	=	z0[50]
tracer  	=	tracer0*1.0