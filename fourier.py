#!/usr/bin/env python2

from qutip import *
import numpy as np
import scipy as scp
import matplotlib.pyplot as plt

def inverse_2x2(m):
	if (m[0][0]*m[1][1] != m[0][1]*m[1][0]):
		return 1/(m[0][0]*m[1][1]- m[0][1]*m[1][0]) * \
					np.array([[m[1][1], -m[0][1]], [-m[1][0], m[0][0]]] ) 
	else:
		print "Error: Singular matrix. E: ", E
# matrix can be any hamiltonian (in matrix form, NOT qubit object). This can be 
# changed if necessary
def GE(matrix, E):
		n = len(matrix)
		Greens = ( (E-10**(-10)*1j)*qeye(n) - matrix).full()
		#Greens = (E*qeye(n) - matrix).full()
		#if np.linalg.cond(Greens) < 1/1000:
		#		Greens = ( (E-10**(-10)*1j)*qeye(n) - matrix).full()
		return np.linalg.pinv(Greens)

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
