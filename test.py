#!/usr/bin/env python2

from qutip import *
import numpy as np
import scipy as scp
import matplotlib.pyplot as plt

import fourier as fr
import qubit as qubit

qubitR = qubit.qubit(1.0, 1.0, 1.0)
qubitL = qubit.qubit(0.5, 1.0, 1.0)

iters = 400
tlist = np.linspace(0,10*np.pi,iters)
#expectsDownR = expect(qubitR.state(tlist), basis(2,0))
expectsUpR = qubitR.spinUpFrac(tlist)
#expectsDownR = qubitR.spinDownFrac(tlist)
expectsDownR = qubitR.spinDownFrac(tlist)

fig, hardcoded = plt.subplots()
hardcoded.plot( tlist, expectsUpR, label='Expect |0> - |G(t)_11|^2/2pi')
hardcoded.plot( tlist, expectsDownR, label='Expect |1> - |G(t)_12|^2/2pi')
plt.ylabel('Expectation')
plt.xlabel('t')
plt.title('Time evolution of the expectation value of a qubit')
legend = hardcoded.legend(loc = 'upper center', bbox_to_anchor=(0.5,-0.1)\
	,ncol=2, shadow = True)
fig.savefig('Expectation values', bbox_extra_artists=(legend,), bbox_inches='tight')
fig.subplots_adjust(bottom=0.2)
frame = legend.get_frame()
frame.set_facecolor('0.90')

print "\n green's function of the fucking qubit"

print qubitR.fun_up (0)
print qubitR.fun_down(0)

t = np.arange(0,15*np.pi,0.01)
transUp = fr.fourrier_transform_loop(qubitR.fun_up,t)
transDown = fr.fourrier_transform_loop(qubitR.fun_down, t)

fig, fourier = plt.subplots()
fourier.plot(t, np.power(np.abs(transUp),2)/(4*np.pi**2))
fourier.plot(t, np.power(np.abs(transDown),2)/(4*np.pi**2))

fourier.grid()
plt.show()