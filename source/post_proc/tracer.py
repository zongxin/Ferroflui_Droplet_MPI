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

def get_tracer(z,z2,tracer0)

	if (time==0):
		xt0=np.real(tracer0)
		yt0=np.imag(tracer0)


	x 	= np.real(z)
	y 	= np.imag(z)
	xn 	= np.real(z2)
	yn 	= np.imag(z2)

	rbfx = Rbf(x, y, xn, epsilon=2)
	rbfy = Rbf(x, y, yn, epsilon=2)

	xt  = rbfx(xt0, yt0)
	yt  = rbfy(xt0, yt0)

	return xt+1j*yt



