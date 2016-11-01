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
qubitR = qubit.qubit(-2.0, 1.0,  0.75)
qubitL = qubit.qubit(0.5, -0.1, 0.2)

qubitEnt = tensor(qubitR.twolevel(), qubitL.twolevel())
print qubitR.twolevel()
print qubitL.twolevel()
print qubitEnt
print "ptrace(0): ", qubitEnt.ptrace(0)
print "ptrace(1): ", qubitEnt.ptrace(1)
print "tracedist (0, L) : ", tracedist(qubitEnt.ptrace(0), qubitL.twolevel())
print "tracedist (1, R) : ", tracedist(qubitEnt.ptrace(1), qubitR.twolevel())
print "tracedist (1, L) : ", tracedist(qubitEnt.ptrace(1), qubitL.twolevel())
print "tracedist (0, R) : ", tracedist(qubitEnt.ptrace(0), qubitR.twolevel())

channel_length = 5
rand_ch = qubit.channel(channel_length)
energies = np.random.uniform(low=(-1.0), high=1.0, size =(channel_length))
rand_ch.setEnergies(energies)
system_rand = tensor(qubitL.twolevel(), rand_ch.asQobj(), qubitR.twolevel())
ptrace_rand = system_rand.ptrace(2)
result_rand = qubit.qubit(ptrace_rand.full()[0][0],\
	ptrace_rand.full()[1][1],ptrace_rand.full()[0][1])

id_ch = qubit.channel(channel_length)
e = 1
energies0 = np.zeros(channel_length)
for i in range(0, len(energies0)):
	energies0[i] = e
#energies0 = np.array([e,e,e,e,e,e,e,e,e,e])

id_ch.setEnergies(energies0)
system_id = tensor(qubitL.twolevel(), id_ch.asQobj(), qubitR.twolevel())
ptrace_id = system_id.ptrace(2)
result_id = qubit.qubit(ptrace_id.full()[0][0], \
	ptrace_id.full()[1][1], ptrace_id.full()[0][1])

controlR = qubitEnt.ptrace(1)
ctrlR = qubit.qubit(controlR.full()[0][0], \
	controlR.full()[1][1], controlR.full()[0][1])

#print result_rand.twolevel()
#print result_id.twolevel()
#print system_rand.ptrace(2)


# arranging an axis along which to output the time evolution
t = np.arange(0,50,0.01)
"""
up = qubitR.GE11()
down = qubitR.GE12()

transUp = fourier.fourier_transform_loop(up, 5, 0.001, t)
transDown = fourier.fourier_transform_loop(down, 5, 0.001, t)

fig, temp = plt.subplots()
temp.plot(t, np.power(np.abs(transUp),2)/(4*np.pi**2))
temp.plot(t, np.power(np.abs(transDown),2)/(4*np.pi**2))
"""
# creating curried functions for G(E)[1][1] and G(E)[1][2]
# up and down are now functions that can be evaluated as: G11(E) and G12(E)

#GE11 = qubitR.GE11()
#GE12 = qubitR.GE12()

GE11_rand = result_rand.GE11()
GE11_id = result_id.GE11()
GE11R = qubitR.GE11()
GE11_ctrl = ctrlR.GE11()

transGE11_rand = fourier.fourier_transform_loop(GE11_rand, 5, 0.001, t)
transGE11_id = fourier.fourier_transform_loop(GE11_id, 5, 0.001, t)
transGE11R = fourier.fourier_transform_loop(GE11R, 5, 0.001, t)
transGE11_ctrl = fourier.fourier_transform_loop(GE11_ctrl, 5, 0.001, t)

fig, fourier1 = plt.subplots()
#fourier3.plot(t, np.power(np.abs(transUp),2)/(4*np.pi**2))
#fourier3.plot(t, np.power(np.abs(transDown),2)/(4*np.pi**2))
rand11 = fourier1.plot(t, np.power(np.abs(transGE11_rand),2)/(4*np.pi**2), label = "Random energies channel")
unif11 = fourier1.plot(t, np.power(np.abs(transGE11_id),2)/(4*np.pi**2), label = "Uniform energies channel")
fourier1.set_title("G_11(t)")
#fourier1.set_legend(handles =[rand11, unif11, orig11])
fourier1.grid()


GE12_rand = result_rand.GE12()
GE12_id = result_id.GE12()
GE12R = qubitR.GE12()
GE12_ctrl = ctrlR.GE12()

transGE12_rand = fourier.fourier_transform_loop(GE12_rand, 5, 0.001, t)
transGE12_id = fourier.fourier_transform_loop(GE12_id, 5, 0.001, t)
transGE12R = fourier.fourier_transform_loop(GE12R, 5, 0.001, t)
transGE12_ctrl = fourier.fourier_transform_loop(GE12_ctrl, 5, 0.001, t)

fig, fourier2 = plt.subplots()
rand12 = fourier2.plot(t, np.power(np.abs(transGE12_rand),2)/(4*np.pi**2), label = "random energies channel")
unif12 = fourier2.plot(t, np.power(np.abs(transGE12_id),2)/(4*np.pi**2), label = "uniform energies channel")
fourier2.set_title("G_12(t)")
#fourier2.legend(handles =[rand12, unif12, orig12])
fourier2.grid()

fig, fourier3 = plt.subplots()
orig12 = fourier3.plot(t, np.power(np.abs(transGE12R),2)/(4*np.pi**2), label= "original qubit")
ctrl12 = fourier3.plot(t, np.power(np.abs(transGE12_ctrl),2)/(4*np.pi**2), label= "ideally entangled")
fourier3.grid()

fig, fourier4 = plt.subplots()
orig11 = fourier4.plot(t, np.power(np.abs(transGE11R),2)/(4*np.pi**2), label= "Original qubit")
ctrl11 = fourier4.plot(t, np.power(np.abs(transGE11_ctrl),2)/(4*np.pi**2), label= "ideally entangled")
fourier4.grid()

# this took 8 min to run...
plt.show()

"""
E = np.arange(-5, 5, 0.001)
GE_11 = map(GE11, E)
GE_12 = map(GE12, E)


fig, GE11 = plt.subplots()
GE11.plot(E, np.real(GE_11))
GE11.plot(E, np.imag(GE_11))
GE11.set_title("G(E)_11")
GE11.grid()


fig, GE12 = plt.subplots()
GE12.plot(E, np.real(GE_12))
GE12.plot(E, np.imag(GE_12))
GE12.set_title("G(E)_12")
GE12.grid()
"""