#!/bin/bash/env python

from qutip import *
import numpy as np
import scipy as scp
import matplotlib.pyplot as plt
import qubit as qubit

# |0>
up_basis = basis(2,0)
# |1>
down_basis = basis(2,1)

epsilon_up = np.sqrt(.5)
epsilon_down = np.sqrt(.5)


# testing normalization and Qobj functions
vec = (up_basis + down_basis).unit()
qubit = Qobj([[epsilon_up], [epsilon_down]])
qubit = Qobj(vec)

# creating a qubit with equal expecation values of being found in |0> or |1>
# note: qu0 = qu1
qu0 = 0.5 * ket2dm(basis(2,0)) + 0.5*ket2dm(basis(2,1))
qu1 = (ket2dm(basis(2,0)) + ket2dm(basis(2,1))).unit()
# producing  a complex valued qubit
qu2 = (np.complex(3,4)*ket2dm(basis(2,0)) + np.complex(3,4)*ket2dm(basis(2,1))).unit()

print "\n basis up: " , up_basis
print "\n basis down: ",down_basis
#print "\n Superposition: ", vec
#print "\n qubit: ", qubit.full()
print "\n qu1: ", qu1
print "\n qu2: ", qu2

# trace distance/fidelity between 2 impure states (qu1 and qu2) and 
# 2 pure states (|0>, |1>) for comparison 
# For pure state, T = sqrt(1-(F^2))
print "\n trace distance: ", tracedist(qu1,qu2)
print "\n fidelity: ", fidelity(qu1, qu2) 
print" \n purestate(up_basis/down_basis)?  ", tracedist(up_basis, down_basis) == \
			np.sqrt(1- (np.power(fidelity(up_basis, down_basis), 2)))
print" \n purestate (qu1/qu2)? ", tracedist(qu1, qu2) == \
			np.sqrt(1- (np.power(fidelity(qu1, qu2), 2)))

# expectation values of qubits in each pure states
print"\n expect qu1 in up_basis: ", expect(qu1, up_basis)
print"\n expect qu1 in down_basis: ", expect(qu1, down_basis)
print"\n expect q2 in upbasis: ", expect(qu2, up_basis)
print"\n expect q2 in down_basis: ", expect(qu2, down_basis)

# Followed an online example for the following. Some questions...
two_spin = tensor(up_basis, down_basis)
qu1_tensor = tensor(qu1,qu1)
qu2_tensor = tensor(qu2,qu2)

# This should be the engtangled quantity H = tensor(qu1,qu2)
qu_test = tensor(qu1, qu2)

# expecation value in basis_up TENSOR basis_down (why?)
print"\n expect qu1 in 2basis: ", expect(qu1_tensor, two_spin)
print"\n expect qu2 in 2basis: ", expect(qu2_tensor, two_spin)
print"\n entangled qu1/qu2 (?): in 2basis: ", expect(qu_test,two_spin)


qubitR = qubit.Qubit(1.0, 0.5, 1.0)
qubitL = qubit.Qubit(0.5, 1.0, 1.0)
qubitEnt = tensor(qubitR.twolevel(), qubitL.twolevel())

print "\n QubitR: ", qubitR.twolevel()
print "\n QubitL: ", qubitL.twolevel()

print "\n Qubit entangled: ", qubitEnt 
print "\n Entangled qubit expect in up_basis: ", expect(qubitEnt, tensor(qeye(2), basis(2,1)))
print "\n Entangled qubit expect in tensor(up_basis,down_basis): ", expect(qubitEnt, tensor(basis(2,0), basis(2,1)))
print "\n Qubit ptrace: ", qubitEnt.ptrace(0)

channel = tensor (qubitR.twolevel(), qeye(2), qubitL.twolevel())
print "\n channel: ", channel
print "\n channel trace: ", channel.ptrace(0)

print "\n trace dist: ", tracedist(qubitEnt.ptrace(0), channel.ptrace(0))
#print "\n trace dist: ", tracedist(qubitEnt.ptrace(0), qubitEnt.ptrace(1))


#H = tensor(sigmaz(), qubitR.twolevel()) + tensor(qubitL.twolevel(), sigmaz()) \
#		+ 0.5 * tensor(sigmax(), sigmax())
#print "H : ", H
#print "H expect: ", expect(H, tensor(basis(2,0), basis(2,1)))
#print "H ptrace: ", H.ptrace(0)


#print "G11", qubit.Gt_11(np.pi)
#print "G22", qubit.Gt_22(np.pi)
#print "G12", qubit.Gt_12(np.pi)

iters = 400
tlist = np.linspace(0,10*np.pi,iters)
#expectsDownR = expect(qubitR.state(tlist), basis(2,0))
expectsUpR = qubitR.spinUpFrac(tlist)
#expectsDownR = qubitR.spinDownFrac(tlist)
expectsDownR = qubitR.spinDownFrac(tlist)

	# trace of tensor and tracedistance between L and R remains constant

	#print "tracedistance (L/R): ", tracedist(L,R)
	#print "tensor: ", tensor(L, R)
	#print "trace: ", tensor(L,R).trace()

fig, axes = plt.subplots()
axes.plot( tlist, expectsUpR, label='Expect |0> - |G(t)_11|^2/2pi')
axes.plot( tlist, expectsDownR, label='Expect |1> - |G(t)_12|^2/2pi')
plt.ylabel('Expectation')
plt.xlabel('t')
plt.title('Time evolution of the expectation value of a qubit')
legend = axes.legend(loc = 'upper center', bbox_to_anchor=(0.5,-0.1)\
	,ncol=2, shadow = True)
fig.savefig('Expectation values', bbox_extra_artists=(legend,), bbox_inches='tight')
fig.subplots_adjust(bottom=0.2)
frame = legend.get_frame()
frame.set_facecolor('0.90')

plt.show()

