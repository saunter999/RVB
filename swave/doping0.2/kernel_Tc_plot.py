#!/usr/bin/env python
from scipy import *
from pylab import *



def swave(Nk):
	gk=[]
	ks=linspace(-pi,pi,Nk)	
	for kx in ks[:-1]:
	   for ky in ks[:-1]:
	     gk.append(cos(kx)+cos(ky))
	gk=array(gk)
	return gk

def dwave(Nk):
	dk=[]
	ks=linspace(-pi,pi,Nk)	
	for kx in ks[:-1]:
	   for ky in ks[:-1]:
	     dk.append(cos(kx)-cos(ky))
	dk=array(dk)
	return dk

def tanhevl(x):                                               
    if x>15:return 1.                                         
    if x<-15:return -1.                                       
    return tanh(x)      


if __name__=="__main__":

	mu=-0.00209931545345;dop=0.05
	Nk=40;t=0.5;J=0.1
	Ns=(Nk-1.)**2
	gk=swave(Nk)
	e0k=-2.0*t*dop*gk
	dk=dwave(Nk)
	Tm=linspace(1e-6,1.0,10)
	kernel=[]
	for T in Tm:
	    sumk=0.0
	    for idx,ek in enumerate(e0k):
		Ek=abs(ek-mu)
		#sumk+=tanhevl(0.5*Ek/T)*gk[idx]**2/Ek
		sumk+=tanhevl(0.5*Ek/T)*dk[idx]**2/Ek
#		sumk+=tanhevl(0.5*Ek/T)/Ek
	    sumk=J/2.0*sumk/Ns
	    kernel.append(sumk-1.0)
	plot(Tm,kernel,'o-')
	show()
