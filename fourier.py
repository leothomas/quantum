#!/usr/bin/env python2

from qutip import *
import numpy as np
import scipy as scp
import matplotlib.pyplot as plt

step_rec = 0.001
def fourrier_transform_rec(func, t):
	def ft(E,crt):
		if E <= 4:
			return ft(E+step_rec, step_rec*(crt+(func(E)*np.exp(-1j*t*E))))
		else:
			return crt
	return ft 

def fourrier_transform_loop(func,t):
	Emin = -5
	Emax = 5
	crt = 0
	E = Emin
	step = 0.001
	while E<=Emax:
		try:
			crt = crt + (func(E)*np.exp(-1j*t*E))
		except:
			print "singular matrix at E = ", E
		E = E+step
		#if np.abs(crt) >= 1000:
		#	print "E --> ", E
	return step*crt

def fourier_transform(func, t):
	return scp.integrate.quad(lambda x: func(x)*np.exp(-1j*t*x), -scp.inf, scp.inf)
