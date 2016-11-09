#!/bin/bash/env python
import numpy as np
import scipy as scp
import matplotlib.pyplot as plt
import fourier as fourier
from qutip import *

class qubit:

	def __init__(self, Eup, Edown, T):
		self.Eup = Eup
		self.Edown = Edown
		self.T = T

	def twolevel(self):
		return Qobj([[self.Eup, self.T], [np.conj(self.T), self.Edown]])

	# curried functions: define fun_down takes a qubit as input and returns a 
	# function of E. This way when we definte func = fun_down(qubit) func is a 
	# function itself, and we can then call func(E) for any E and get the expectation
	def GE11(self):
		def fu(E):
			GEqubit = fourier.GE(self.twolevel().full(), E)
			return GEqubit[0][0]
		return fu
	def GE12(self):
		def fd(E):
			GEqubit = fourier.GE(self.twolevel().full(), E)
			return GEqubit[0][1]
		return fd 

class channel:

	def __init__(self, n):
		self.n = n
		self.matrix = qeye(n).full()

	# input is either an array of energies or a 2x2 matrix. If later, Ematrix 
	# becomes the energies, if not, the array values are inputed along the main 
	# diagonal of matrix. 
	def setEnergies(self, energies):
		if len(energies) != self.n: 
			print "Given energy array/matrix not same size as channel"
			print "\n self.n: ", self.n
			print "\n len(energies): ", len(energies)
		
		elif energies.ndim == 2:
			self.matrix = Ematrix
		
		elif energies.ndim == 1: 
			for i in range(0,len(energies)):
				self.matrix[i][i] = energies[i]
		else: 
			print "Invalid input matrix dimensions"

	def asQobj(self):
		return Qobj(self.matrix)
