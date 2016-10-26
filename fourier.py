#!/usr/bin/env python2

from qutip import *
import numpy as np
import scipy as scp
import matplotlib.pyplot as plt

# matrix can be any hamiltonian (in matrix form, NOT qubit object). This can be 
# changed if necessary
def GE(matrix, E):
		threshold = 0.00001
		n = len(matrix)
		Greens = (E*qeye(n) - matrix).full()
		# cond(Greens) increases as the matrix approaches singularity. 
		# These conditionals catch near singular matrix and subtract 
		# i*Epsilon such that the determinant is never = 0
		# THIS IS THE MOST COMPUTATIONALLY INTENSIVE ASPECT OF THIS PROGRAM
		if np.linalg.cond(Greens) < 1/threshold:
			return np.linalg.inv(Greens)
		else:
			print "Bad condition matrix at E = ", E
			#return np.array([[0,0],[0,0]])
			return np.linalg.inv((1j*0.0001*qeye(n) - matrix).full())

# tail recursive implementation of a fourrier transform.
# Python limits number of recursions, so this algorithm is not 
# particularly helpful. I have left if for reference
step_rec = 0.001
def fourier_transform_rec(func, t):
	def ft(E,crt):
		if E <= 4:
			return ft(E+step_rec, step_rec*(crt+(func(E)*np.exp(-1j*t*E))))
		else:
			return crt
	return ft 

# iterative implementation of fourrier transform using left hand rectangle
# rule for integration
def fourier_transform_loop(func, Emax, dx, t):
	# initialize E to lower integral limit, and current sum to 0
	E = -Emax
	crt = 0
	while E<=Emax:
		try:
			# func(E) => whichever element of G(E) we are interested in: G11,
			# G12, etc, evalutes it at thay particular E value, and adds it
			# to the current sum. Catches an exception if the matrix is singular
			# or if G(E) returns an otherwise invalid value. 
			crt = crt + (func(E)*np.exp(-1j*t*E))
		except:
			print "singular matrix at E = ", E
		E = E+dx
		# multiplies total sum by dx
	return dx*crt

# implementation using built in integration function to integrate more efficiently
# and using greater limits. Issues with return type of func(x)
def fourier_transform(func, t):
	return scp.integrate.quad(lambda x: func(x)*np.exp(-1j*t*x), -scp.inf, scp.inf)

# implementation using trapz, numerical integration using trapezoids. Run-time
# issues present. Havent had time to debug. 
def fourier_trans(func, t):
	dx = 0.001
	x = np.arange(-15, 15, dx)
	fun = np.zeros(len(x))
	for i in range(0, len(x)):
		fun[i] = func(x[i])
		i=i+1
	return scp.trapz(fun * np.exp(-1j*t*x))