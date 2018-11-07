#!/usr/bin/env python
from scipy import *
from pylab import *



if __name__=="__main__":
	data1=loadtxt("Gap_doping.dat").transpose()
	data2=loadtxt("Tc_doping.dat").transpose()
	plot(data1[0],data1[1],"o-")
	xlabel("doping")
	ylabel("Gap")
	axvline(x=data2[1])	
	show()
