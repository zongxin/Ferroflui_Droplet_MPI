import os
import sys
import numpy as np
from pdb import set_trace
from matplotlib import rc as matplotlibrc
from matplotlib.patches import FancyArrowPatch
import matplotlib.pyplot as plt
# plt.switch_backend('agg')
from itertools import product, combinations
from scipy import signal
import pickle
from scipy.signal import argrelextrema
import functions as func
import matplotlib.collections as collections

from scipy.interpolate import interp1d

matplotlibrc('text',usetex=True)
matplotlibrc('font', family='serif')

figheight      = 6
figwidth      = 6
lineWidth      = 1.5
textFontSize   = 18
gcafontSize    = 20
legendFontSize = 20

EVO_CURV = 1
EVO 	 = 0
TRACE 	 = 0
FFT 	 = 0
TRACER 	 = 0
WAVE 	 = 0
if EVO_CURV:
	for i in range (0,231):
		pkl_name	= './data/data2_'+str(i+1)+'.pkl'
		data 		= pickle.load( open( pkl_name, "rb" ) )

		N 		= data['N']	
		t 		= data['time']
		z0 		= data['z0']				
		z 		= data['z']		
		gamma	= data['gamma']	
		kappa	= data['kappa']	

		angle=np.angle(z)
		angle[np.where(angle<0)]=angle[np.where(angle<0)]+2*np.pi
		xi0=np.argmin(angle)
		z=np.concatenate((z[xi0:],z[:xi0]))

		alpha = np.linspace(0,2*np.pi,N, endpoint=False)

		fig = plt.figure(0, figsize=(figwidth*1.5,figheight*0.7))
		ax  = fig.add_subplot(121,alpha=1)#
		plt.axes(ax)
		ax.grid(False)
		a=1.5
		ax.set_xlim([-a,a])
		ax.set_ylim([-a,a])	
		# plt.axis('equal')	


		plt.plot(np.real(z0),np.imag(z0),'r')		
		plt.plot(np.real(z),np.imag(z),'k')
		plt.scatter(np.real(z[::2]),np.imag(z[::2]),c='k',marker='.',s=8.)		
		ax.set_xlabel(r'$x$',fontsize=1.2*gcafontSize)
		ax.set_ylabel(r'$y$',fontsize=1.2*gcafontSize)
		plt.setp(ax.get_xticklabels(),fontsize=1.*gcafontSize)
		plt.setp(ax.get_yticklabels(), fontsize=1.*gcafontSize)

		ax  = fig.add_subplot(122,alpha=1)#

		plt.axes(ax)
		# plt.plot(alpha,np.abs(z),'k',lw=lineWidth)
		angle=np.angle(z)
		angle[np.where(angle<0)]=angle[np.where(angle<0)]+2*np.pi
		xi0=np.argmin(angle)
		z=np.concatenate((z[xi0:],z[:xi0]))
		angle=np.concatenate((angle[xi0:],angle[:xi0]))
		
		rr=np.abs(z)
		plt.plot(angle,rr,'k-')
		plt.plot(angle,np.abs(z0),'r-')
		ax.set_ylim([0.8,1.3])
		time_string=r'$t=$'+str(round(t,4))
		ax.text(3,0.82,time_string,fontsize=0.8*gcafontSize)					

		ax.set_xlabel(r'$\alpha$',fontsize=1.2*gcafontSize)
		ax.set_ylabel(r'$h$',fontsize=1.2*gcafontSize)
		plt.setp(ax.get_xticklabels(),fontsize=1.*gcafontSize)
		plt.setp(ax.get_yticklabels(), fontsize=1.*gcafontSize)
		figname="frame_"+str(i)+'.png'
		plt.tight_layout()		
		plt.savefig('./tttest/'+figname)
		# plt.show()
		plt.close()
		print (figname," saved!")

if EVO:
	for i in range (500-1,500):
		pkl_name	= './data/data_'+str(i+1)+'.pkl'
		data 		= pickle.load( open( pkl_name, "rb" ) )

		N 		= data['N']	
		t 		= data['time']
		z0 		= data['z0']				
		z 		= data['z']		
		gamma	= data['gamma']	
		kappa	= data['kappa']	

		alpha = np.linspace(0,2*np.pi,N, endpoint=False)

		fig = plt.figure(0, figsize=(figwidth*0.8,figheight*0.75))
		ax  = fig.add_subplot(111,alpha=1)#
		plt.axes(ax)
		ax.grid(False)
		ax.set_xlim([-1.5,1.5])
		ax.set_ylim([-1.5,1.5])	
		# plt.axis('equal')	
		plt.plot(np.real(z0),np.imag(z0),'r')	
		plt.scatter(np.real(z0[::2]),np.imag(z0[::2]),color='r',marker='.',s=8.)		
		plt.plot(np.real(z),np.imag(z),'k')
		plt.scatter(np.real(z[::2]),np.imag(z[::2]),c='k',marker='.',s=8.)		
		ax.set_xlabel('x',fontsize=1.2*gcafontSize)
		ax.set_ylabel('y',fontsize=1.2*gcafontSize)


		figname="path_"+str(i)+'.png'
		plt.tight_layout()		
		plt.savefig('./movie/'+figname)
		# plt.show()
		plt.close()
		print (figname," saved!")


if TRACE:
	fig = plt.figure(0, figsize=(figwidth*0.8,figheight*0.75))
	ax  = fig.add_subplot(111,alpha=1)#	

	for i in range (0,578):
		pkl_name	= './data/data2_'+str(i+1)+'.pkl'
		data 		= pickle.load( open( pkl_name, "rb" ) )

		N 		= data['N']	
		t 		= data['time']
		z0 		= data['z0']				
		z 		= data['z']		
		gamma	= data['gamma']	
		kappa	= data['kappa']	

		alpha = np.linspace(0,2*np.pi,N, endpoint=False)

		plt.axes(ax)
		ax.grid(False)
		a=3
		ax.set_xlim([-a,a])
		ax.set_ylim([-a,a])	
		# plt.axis('equal')	
		if (i==0):
			plt.plot(np.real(z0),np.imag(z0),'r')
		elif ((i%50==0)):		
			plt.plot(np.real(z),np.imag(z),'k',alpha=0.7)
			# print i 
		# elif ((i>30)&(i%30==0)):		
		# 	plt.plot(np.real(z),np.imag(z),'k',alpha=0.7)	
		# 	print i 					
			# plt.scatter(np.real(z[::8]),np.imag(z[::8]),c='k',marker='.',s=8.)	

	# Qc=2e-3
	# T=(5e-3)**2*(2e-3)*(2*np.pi)/Qc
	# tf=t*T
	tf=t*1.0
	set_trace()
	time_string='tf = '+str(round(tf,3))+' s'
	N_string='N = '+str(len(z))
	ax.text(-9,-9,time_string,fontsize=gcafontSize)
	ax.text(0,-9,N_string,fontsize=gcafontSize)	
	ax.set_xlabel('x',fontsize=1.2*gcafontSize)
	ax.set_ylabel('y',fontsize=1.2*gcafontSize)
	figname='trace2.png'
	plt.tight_layout()		
	plt.savefig(figname)
	# plt.show()
	plt.close()
	print (figname," saved!")



if FFT:

	Nf= 59
	tracer=np.zeros(Nf,dtype='complex')
	k4_mag=np.zeros(Nf)	
	t_list=np.zeros(Nf)
	plotj=0
	for i in range (0,Nf):
		pkl_name	= './data/data2_'+str(i+1)+'.pkl'
		data 		= pickle.load( open( pkl_name, "rb" ) )
		tracer[i]  = data['tracer']

	for i in range (0,Nf):
		if (i%2==0):
			plotj=plotj+1
			pkl_name	= './data/data2_'+str(i+1)+'.pkl'
			data 		= pickle.load( open( pkl_name, "rb" ) )

			N 		= data['N']	
			t 		= data['time']
			z0 		= data['z0']				
			z 		= data['z']		
			gamma	= data['gamma']	
			kappa	= data['kappa']	
			traceri  = data['tracer']

			alpha = np.linspace(0,2*np.pi,N, endpoint=False)

			#=======================================================================

			fig = plt.figure(0, figsize=(figwidth*1.5*0.8,figheight*1.4*0.8))
			ax  = fig.add_subplot(221,alpha=1)#
			plt.axes(ax)
			ax.grid(False)
			a=1.5
			ax.set_xlim([-a,a])
			ax.set_ylim([-a,a])	
			# plt.axis('equal')	
			plt.plot(np.real(z0),np.imag(z0),'r')	
			# plt.plot(np.cos(alpha),np.sin(alpha),'b--')			
			plt.plot(np.real(z),np.imag(z),'k')
			# plt.plot(np.real(tracer),np.imag(tracer),'g')		
			plt.scatter(np.real(z[::4]),np.imag(z[::4]),c='k',marker='.',s=8.)		
			plt.scatter(np.real(z[0]),np.imag(z[0]),c='r',marker='o',s=10.)				
			# plt.scatter(np.real(traceri),np.imag(traceri),c='g',edgecolors='g',marker='o',s=8.)			
			ax.set_xlabel('x',fontsize=1.2*gcafontSize)
			ax.set_ylabel('y',fontsize=1.2*gcafontSize)

			#=======================================================================

			ax  = fig.add_subplot(222,alpha=1)#
			plt.axes(ax)
			plt.plot(alpha,kappa,'k',lw=lineWidth)
			time_string='t='+str(t)

			b =4000
			ax.set_ylim([-2,8])
			ax.text(3,-1.5,time_string,fontsize=1.*gcafontSize)					

			ax.set_xlabel('alpha',fontsize=1.2*gcafontSize)
			ax.set_ylabel('kappa',fontsize=1.2*gcafontSize)

			#=======================================================================

			ax  = fig.add_subplot(223,alpha=1)#
			y1=0
			y2=2

			m=40
			mlist=np.array(range(0,m))
			plt.axes(ax)		
			# collection = collections.BrokenBarHCollection.span_where(
			#     np.array(range(0,m)), ymin=y1, ymax=y2,where=mlist<=12, facecolor='grey', alpha=0.5)
			# ax.add_collection(collection)



			z_hat=np.abs(np.fft.fft(np.abs(z)))
			if (i==0):
				z_hat0=np.abs(np.fft.fft(np.abs(z0)))

			ax.set_ylim([y1,y2])		
			ax.set_xlim([0,m])	
			plt.plot(np.array(range(1,m)),z_hat[1:m],'k',marker='.',linewidth=lineWidth)
			plt.plot(np.array(range(1,m)),z_hat0[1:m],'r--',marker='.')	
			# ax.axvline(12,c='k',ls='--')
			# ax.axvline(7,c='k',ls='--')	
			ax.set_xlabel('k',fontsize=1.2*gcafontSize)
			ax.set_ylabel('fft',fontsize=1.2*gcafontSize)
			k4_mag[i]  = z_hat[4]
			t_list[i]  = t

			#=======================================================================

			ax  = fig.add_subplot(224,alpha=1)#
			plt.axes(ax)
			# aa=np.where(k4_mag>0)

			# bb=aa[-1]
			# plt.plot(t_list[aa],k4_mag[aa],'k-')
			# plt.scatter(t_list[aa][-1],k4_mag[aa][-1],c='k',marker='o',s=8.)		
			# ax.set_ylim([0,6])		
			# ax.set_xlim([0,800])
			plt.plot(alpha,np.abs(z0),'r-',linewidth=lineWidth)	
			plt.plot(alpha,np.abs(z),'k',linewidth=lineWidth)
			
			aa=0.1
			ax.set_ylim([1-aa,1.05+aa])			
			ax.set_xlabel('alpha',fontsize=1.2*gcafontSize)
			ax.set_ylabel('r',fontsize=1.2*gcafontSize)



			# z_ang=func.cal_phase(z)
			# # set_trace()
			# if (i==0):
			# 	z_ang0=func.cal_phase(z0)
			# plt.plot(np.array(range(1,m)),z_ang[1:m],'k',marker='.',linewidth=lineWidth)
			# plt.plot(np.array(range(1,m)),z_ang0[1:m],'r--',marker='.')	


			figname="ttttest_"+str(plotj)+'.png'
			plt.tight_layout()		
			plt.savefig('./tttest/'+figname)
			# plt.show()
			plt.close()
			print (figname," saved!")
		# set_trace()
if TRACER:


	Nf= 83
	tracer=np.zeros(Nf,dtype='complex')
	t_list=np.zeros(Nf)
	I=np.zeros(Nf)
	C=np.zeros(Nf)

	for i in range (0,Nf):
		pkl_name	= './data/data2_'+str(i+1)+'.pkl'
		data 		= pickle.load( open( pkl_name, "rb" ) )
		tracer[i]   = data['tracer']
		t_list[i] 	= data['time']
		z 			= data['z']	
		z0 			= data['z0']		
		N 			= data['N']	

		ds=np.abs(z[2]-z[1])
		z[2]-z[1]
		L = N*ds

		A=func.cal_area(z,ds)
		I[i]=L**2/(4*np.pi*A)
		C[i]=np.max(np.abs(z))/np.min(np.abs(z))

	ymag=np.abs(tracer)
	Lmax=argrelextrema(ymag, np.greater)



	fig = plt.figure(0, figsize=(figwidth*1.6,figheight*1.5))

	#-------------------------------------------------------------

	ax  = fig.add_subplot(221,alpha=1)#	

	alpha = np.linspace(0,2*np.pi,N, endpoint=False)
	plt.plot(np.real(z0),np.imag(z0),'r-')		
	plt.plot(np.real(tracer),np.imag(tracer),c='g')
	plt.plot(np.cos(alpha),np.sin(alpha),'b--')	
	# ax.set_ylim([-1.2,1.2])		
	# ax.set_xlim([-1.2,1.2])
	ax.set_ylim([0.8,1.2])		
	ax.set_xlim([-0.2,0.2])
	ax.set_xlabel('x',fontsize=1.2*gcafontSize)
	ax.set_ylabel('y',fontsize=1.2*gcafontSize)

	#-------------------------------------------------------------

	ax  = fig.add_subplot(222,alpha=1)#	

	plt.plot(t_list,np.abs(tracer),c='g')
	plt.plot(t_list,np.ones(Nf),'k--')
	ax.set_xlabel('x',fontsize=1.2*gcafontSize)
	ax.set_ylabel('y',fontsize=1.2*gcafontSize)
	ax.set_ylim([0.8,1.15])		
	ax.set_xlim([0,800])
	#-------------------------------------------------------------

	ax  = fig.add_subplot(223,alpha=1)#	

	plt.plot(t_list,I,'k-')
	plt.plot(t_list,C,'k--')
	ax.set_xlabel('x',fontsize=1.2*gcafontSize)
	ax.set_ylabel('y',fontsize=1.2*gcafontSize)
	ax.set_xlim([0,800])
	ax.set_ylim([1.0,1.2])	

	figname='tracer.png'
	plt.tight_layout()		
	plt.savefig(figname)
	plt.show()
	plt.close()
	print (figname," saved!")
	print (t_list[Lmax])

if WAVE:

	# LaTeX setup
	matplotlibrc('text.latex', preamble=r'\usepackage{color}')
	matplotlibrc('text',usetex=True)
	matplotlibrc('font', family='serif')

	Nf=50
	# t=np.array(range(Nf))*100*0.0000005
	# a=np.loadtxt("./data/Displacement_0.txt")
	Np=512
	# Npindex=np.array(range(Np))
	x=np.linspace(0,1,Np)*np.pi*2


	Yd=np.zeros((Nf,Np))

	fig = plt.figure(0, figsize=(figwidth*1,figheight*1.0))
	ax=fig.add_subplot(1,1,1)

	for i in range (0,Nf):
		pkl_name	= './data/data2_'+str(i+1)+'.pkl'
		data 		= pickle.load( open( pkl_name, "rb" ) )
		# tracer[i]   = data['tracer']
		# t_list[i] 	= data['time']
		z 			= data['z']	
		z0 			= data['z0']		
		N 			= data['N']	




		angle=np.angle(z)
		angle[np.where(angle<0)]=angle[np.where(angle<0)]+2*np.pi
		xi0=np.argmin(angle)

		# interpolate_x = interp1d(yy,xx, kind='cubic')
		# print i
		# set_trace()	
		# x0	= interpolate_x(0.)
		# dis=np.abs(xx-x0)

		z=np.concatenate((z[xi0:],z[:xi0]))
		angle=np.concatenate((angle[xi0:],angle[:xi0]))
		
		Yd[i,:]=(np.abs(z)-1)*20
		
		# if (i>=55):
		# 	plt.plot(angle,Yd[i,:]+np.int(i/1)*0.8,'g-',)
		# else:	
		plt.plot(angle,Yd[i,:]+np.int(i/1)*0.8,'k-')
	# ax.set_xlabel('x',fontsize=1.2*gcafontSize)
	# ax.set_ylabel('y',fontsize=1.2*gcafontSize)
	# ax.set_ylim([0.8,1.15])		
	ax.set_xlim([0,np.pi*2])
	ax.set_xticks([0,np.pi,np.pi*2])
	# ax.set_yticks([0.,150.*0.8,200.*0.8,300.*0.8,400.*0.8])	
	# # ax.set_yticks([])
	# ax.set_yticklabels([0,150,200,300,400])
	ax.set_xticklabels(['0', r'$\pi$', r'2$\pi$'],fontsize=1.2*gcafontSize)
	# ax.spines['left'].set_visible(False)	
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.xaxis.set_ticks_position('bottom')
	# ax.axis('off')				
	figname = "wave.png"
	plt.tight_layout()
	plt.savefig('./post_proc/'+figname)
	plt.show()
	plt.close()	