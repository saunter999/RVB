#!/usr/bin/env python
from scipy import *
from pylab import *
from scipy import optimize




def swave(Nk):
	gk=[]
	ks=linspace(-pi,pi,Nk)	
	for kx in ks[:-1]:
	   for ky in ks[:-1]:
	     gk.append(cos(kx)+cos(ky))
	gk=array(gk)
	return gk

def tanhevl(x):                                               
    if x>15:return 1.                                         
    if x<-15:return -1.                                       
    return tanh(x)      

def sceqs(mu,delta,temp,doping):
	deltap=0.0;mup=0.0
	gapksum=0.0
	for idx,ek in enumerate(e0k):
	    #Ek=( (ek-mu)**2+(J*gk[idx]*delta)**2 )**0.5
	    Ek=( (ek-mu)**2+(J*gk[idx]*delta)**2 )**0.5
	    #gapksum+=tanhevl(0.5*Ek/temp)*delta*gk[idx]**2/Ek
	    gapksum+=tanhevl(0.5*Ek/temp)*delta*gk[idx]**2/Ek
	gapksum=gapksum/Ns
	deltap=gapksum*J/2.0
	mup=optimize.brentq(funcmu,-20.,20.,args=(delta,temp,doping))
	return (mup,deltap)

def funcmu(mu,delta,temp,doping):
	sumk=0.0
	for idx,ek in enumerate(e0k):
	    Ek=( (ek-mu)**2+(J*gk[idx]*delta)**2 )**0.5
	    #Ek=( (ek-mu)**2+(J*gk[idx]*delta)**2 )**0.5
	    sumk+=tanhevl(0.5*Ek/temp)*(ek-mu)/Ek
	sumk=sumk/Ns
	return sumk-doping
	
def scTceqs(mu,Tc,doping):
	mup=optimize.brentq(funcmu,-20.,20.,args=(0.0,Tc,doping))
	#Tcp=optimize.brentq(funcTc,1e-9,1.0,args=(mu,))
	sumk=0.0
	for idx,ek in enumerate(e0k):
	    Ek=abs(ek-mu)
	    sumk+=tanhevl(0.5*Ek/Tc)*gk[idx]**2/Ek
	    #sumk+=tanhevl(0.5*Ek/Tc)*gk[idx]**2/Ek
	sumk=J/2.0*sumk/Ns
	Tcp=Tc*sumk
	return(mup,Tcp)

def funcTc(Tc,mu):
	sumk=0.0
	for idx,ek in enumerate(e0k):
	    Ek=abs(ek-mu)
	    sumk+=tanhevl(0.5*Ek/Tc)*gk[idx]**2/Ek
	    #sumk+=tanhevl(0.5*Ek/Tc)*gk[idx]**2/Ek
	sumk=J/2.0*sumk/Ns
	return sumk-1.0
	
	
if __name__=="__main__":
	"""
	Mean field soltion of RVB theory of t-j model.
	1)We solve (Delta,mu) as a function of doping d at a given temperature T.
	2)We determine Tc of RVB singlet(not superconducting Tc...)as a function of doping d
	Reference:The resonating valence bond state and high Tc superconductivity-A mean field theory. 1987 Anderson
	"""
	Nk=100;t=0.5;
	Ns=(Nk-1.0)**2
	J=0.1
	ev_to_k=11600
	gk=swave(Nk)
	
	dopls=linspace(0.05,0.2,1)  ## Numerically not stable/accurate for low doping dope<0.05..Physically,Anderson argues that the theory fails in low doping
	##solving self-consistent equaitons by iteration
	fgap=open("Gap_doping.dat",'w')
	print>>fgap,"#","temperature","RVB Gap,","Tc_RVB,","chemical potential mu","#num of iterations"
	#1)We solve (Delta,mu) as a function of doping d for several temperatures T.
	Tm=linspace(1e-5,0.01,30)
	for dop in dopls:
	  print>>fgap,"#","doping=",dop
	  e0k=-2.0*t*dop*gk
	  mu=-0.50;Delta=1.0;   ## Be careful to check that converged solution is not depending on inital values of parameters.
	  Nit=1000;err=1e-7
	  for Temp in Tm:
	    print "____________________________"
	    print "T=",Temp
	    it=1;
	    while True:
	       (mup,Deltap)=sceqs(mu,Delta,Temp,dop)	
	       it+=1
	       diffmu=mup-mu
	       diffDelta=Deltap-Delta
	       if abs(diffmu)<err or abs(diffDelta)<err or it==Nit:
		    print>>fgap,Temp,Delta,Delta*ev_to_k,mu,it
		    break;
	       mu=(mup+mu)/2.
	       Delta=(Deltap+Delta)/2.
	       #print "mu","Delta","iteration","diffmu",'diffDelta'
	       #print mu,Delta,it,diffmu,diffDelta
	  
	#2)We determine Tc of RVB singlet(not superconducting Tc...)as a function of doping d.
	### not good to determine Tc for s wave pairing....(because Tc is usually very small ???)
	fTc=open("Tc_doping.dat",'w')
	print>>fTc,"#doping","Tc(ev)","Tc(K)","mu","#num of iterations"
	for dop in dopls:
	  print "____________________________"
	  mu=1.0;Tc=1.0;
	  it=1;Nit=100;err=1e-7
	  while True:
		 (mup,Tcp)=scTceqs(mu,Tc,dop)
		 it+=1
		 diffmu=mup/mu-1.0
		 diffTc=Tcp/Tc-1.0
		 if abs(diffmu)<err or abs(diffTc)<err or it==Nit:
		      print>>fTc,dop,Tc,Tc*ev_to_k,mu,it
		      break;
		 mu=(mup+mu)/2.
		 Tc=(Tcp+Tc)/2.
