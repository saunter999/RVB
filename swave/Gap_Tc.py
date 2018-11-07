#!/usr/bin/env python
from scipy import *
from pylab import *




if __name__=="__main__":
	data=loadtxt("doping_Tc.dat").transpose()
	plot(data[0],data[1],'o-')
	xlabel("doping",size=19)
	ylabel("Tc[K]",size=19)
	ylim([0,20])
	show()
