import os
import sys
import numpy as np
from pdb import set_trace

# Restart key

Restart   	= 1
Res_i		= 57
Resname		= "./data/data2_"+str(Res_i)+".pkl"


# TIME SCHEME
#-------------------

RK 	  		= 0
CN    		= 1
AB    		= 0


Plot_output = 2000#125
Text_output = 20
Data_output = 2000#125

R 	  		= 1
Ng 	  		= 1024
dt 	  		= 1e-7

tau	  		= 1.
chi   		= 1.0
NB1	  		= 1.0
km 	  		= 7
NB2	  		= (3*km**2-1+2*NB1/R)/(2*(1+chi)*R**3)

print ('====================================')
print ('NB1,NB2,NB1*NB2')
print (NB1,NB2,NB1*NB2)
print ('====================================')

NB1	  		= NB1*2
NB2	  		= NB2*2

n 			= 7
alpha 		= np.linspace(0,2*np.pi,Ng, endpoint=False)

nghost 		= 10
nghost_in	= 20

NORM_STRESS = 1