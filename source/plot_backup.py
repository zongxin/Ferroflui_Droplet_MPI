	if (Plot_quiver&(i==0)):
		fig = plt.figure(0, figsize=(figwidth*1.5,figheight*0.7))
		ax  = fig.add_subplot(121,alpha=1)#
		plt.axes(ax)
		ax.grid(False)

		plt.quiver(np.real(z),np.imag(z),np.real(u1),np.imag(u1))
		ax.set_xlabel('x',fontsize=1.2*gcafontSize)
		ax.set_ylabel('y',fontsize=1.2*gcafontSize)


		ax  = fig.add_subplot(122,alpha=1)#
		plt.axes(ax)
		plt.plot(alpha,kappa,'k',lw=lineWidth)
		ax.set_ylim([0,2])
		ax.set_xlabel('alpha',fontsize=1.2*gcafontSize)
		ax.set_ylabel('kappa',fontsize=1.2*gcafontSize)
		time_string='t='+str(t)
		ax.text(5,1.2,time_string,fontsize=1.*gcafontSize)
		figname="quiver.png"
		plt.tight_layout()		
		plt.savefig('./figure/'+figname)
		# plt.show()
		plt.close()
		print figname," saved!"
		plotj=plotj+1