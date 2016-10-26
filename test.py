#!/usr/bin/env python2

from qutip import *
import numpy as np
import scipy as scp
import matplotlib.pyplot as plt
# fourier contains necesary methods for Green's function analysis of 
# the qubit dynamics
import fourier as fourier
# qubit contains necesary method fo qubit representation and operations, such
# as 2x2 representation, and object initialization. It also contains hard-
# coded functions for modeling the time evolution of the qubit, as were found
# symbolically (by hand) in qubit.pdf
import qubit as qubit

# initializing 2 qubit objects (Eup, Edown, coupling)
qubitR = qubit.qubit(2.0, 2.0, 0.5)
qubitL = qubit.qubit(0.5, 1.0, 1.0)

qubitEnt = tensor(qubitR.twolevel(), qubitL.twolevel())
print qubitR.twolevel()
print qubitL.twolevel()
print qubitEnt
print "ptrace(0): ", qubitEnt.ptrace(0)
print "ptrace(1): ", qubitEnt.ptrace(1)

# arranging an axis along which to output the time evolution
t = np.arange(0,10*np.pi,0.01)

# calculating the fourier transform using the hardcoded functions from qubit.py
expectsUpR = qubitR.spinUpFrac(t)
expectsDownR = qubitR.spinDownFrac(t)

fig, hardcoded = plt.subplots()
hardcoded.plot( t, expectsUpR, label='Expect |0> - |G(t)_11|^2/2pi')
hardcoded.plot( t, expectsDownR, label='Expect |1> - |G(t)_12|^2/2pi')
plt.ylabel('Expectation')
plt.xlabel('t')
plt.title('Time evolution of the expectation value of a qubit')
legend = hardcoded.legend(loc = 'upper center', bbox_to_anchor=(0.5,-0.1)\
	,ncol=2, shadow = True)
fig.savefig('Expectation values', bbox_extra_artists=(legend,), bbox_inches='tight')
fig.subplots_adjust(bottom=0.2)
frame = legend.get_frame()
frame.set_facecolor('0.90')
hardcoded.grid()
# creating curried functions for G(E)[1][1] and G(E)[1][2]
# up and down are now functions that can be evaluated as: up(E) and down(E)
# for any E
up = qubitR.fun_up()
down = qubitR.fun_down()

# calculating fourier transform using an iterative method. This may take up
# to several minutes, since calculating the matrix inverse for every single 
# value of E (as we integrate over it) is very computationally intensive, 
# especially as the matrix becomes near-singular
transUp = fourier.fourier_transform_loop(up, 15, 0.001, t)
transDown = fourier.fourier_transform_loop(down, 15, 0.001, t)

fig, fourier0 = plt.subplots()
fourier0.plot(t, np.power(np.abs(transUp),2)/(4*np.pi**2))
fourier0.plot(t, np.power(np.abs(transDown),2)/(4*np.pi**2))
fourier0.set_title("-15 to 15, dx = 0.001")
fourier0.grid()
"""
transUp = fourier.fourier_transform_loop(up, 5, 0.001, t)
transDown = fourier.fourier_transform_loop(down, 5, 0.001, t)

fig, fourier2 = plt.subplots()
fourier2.plot(t, np.power(np.abs(transUp),2)/(4*np.pi**2))
fourier2.plot(t, np.power(np.abs(transDown),2)/(4*np.pi**2))
fourier2.set_title("-5 to 5, dx = 0.001")
fourier2.grid()

transUp = fourier.fourier_transform_loop(up, 50, 0.01, t)
transDown = fourier.fourier_transform_loop(down, 50, 0.01, t)

fig, fourier3 = plt.subplots()
fourier3.plot(t, np.power(np.abs(transUp),2)/(4*np.pi**2))
fourier3.plot(t, np.power(np.abs(transDown),2)/(4*np.pi**2))
fourier3.set_title("-50 to 50, dx = 0.01")
fourier3.grid()
"""
plt.show()
# this took 8 min to run...