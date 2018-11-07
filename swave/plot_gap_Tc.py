#!/usr/bin/env python
from scipy import *
from pylab import *



if __name__=="__main__":
	ev_to_k=11604
	data1=loadtxt("Gap_doping.dat").transpose()
	data2=loadtxt("Tc_doping.dat").transpose()
	plot(data1[0]*ev_to_k,data1[1]*ev_to_k,"o-")
	xlabel("T[K]",size=19)
	ylabel("Gap[K]",size=19)
	title("Doping=0.05",size=19)
	axvline(x=data2[1])	
	show()
