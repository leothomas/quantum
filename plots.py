import numpy as np
import scipy as scp
import matplotlib.pyplot as plt
from qutip import *
# fourier contains necesary methods for Green's function analysis of 
# the qubit dynamics
import fourier as fourier
# qubit contains necesary method fo qubit representation and operations, such
# as 2x2 representation, and object initialization. It also contains hard-
# coded functions for modeling the time evolution of the qubit, as were found
# symbolically (by hand) in qubit.pdf
import qubit as qubit

qt0 = qubit.qubit(np.sqrt(.2), np.sqrt(.8), 0.75)
qt1 = qubit.qubit(0.4, -0.6, 0.75)

entangled = tensor(qt0.twolevel(), qt1.twolevel()) 

print "\n qubit0: ", qt0.twolevel()
print "\n qubit1: ", qt1.twolevel()
print "\n entangled (tensored) qubit: ", entangled

entangled0 = entangled.ptrace(0)
entangled1 = entangled.ptrace(1)

print "\n ptrace(0): ", entangled0
print "\n ptrace(1): ", entangled1

qt0_ent = qubit.qubit(entangled0.full()[0][0], \
		entangled0.full()[1][1], entangled0.full()[0][1])
qt1_ent = qubit.qubit(entangled1.full()[0][0], \
		entangled1.full()[1][1], entangled1.full()[0][1])

channel_length = 5
rand_ch = qubit.channel(channel_length)
energies = np.random.uniform(low=(-1.0), high=1.0, size =(channel_length))
#energies = np.ones(channel_length)
rand_ch.setEnergies(energies)
for i in range(0, len(rand_ch.matrix)):
	if i < len(rand_ch.matrix)-1:
		rand_ch.matrix[i+1][i] = 1.0
		rand_ch.matrix[i][i+1] = 1.0

print rand_ch.matrix

system_rand = tensor(qt0.twolevel(), rand_ch.asQobj(), qt1.twolevel())
ptrace_rand = system_rand.ptrace(2)
#qt1_ent = qubit.qubit(ptrace_rand.full()[0][0],\
#	ptrace_rand.full()[1][1],ptrace_rand.full()[0][1])


print qt0_ent.twolevel()
print qt1_ent.twolevel()

t = np.arange(0, 10*np.pi, 0.01)
"""
GE11 = qt1_ent.GE11()
GE12 = qt1_ent.GE12()
#GE112 = qt1.GE12()

transGE11 = fourier.fourier_transform_loop(GE11, 5, 0.001, t)
transGE12 = fourier.fourier_transform_loop(GE12, 5, 0.001, t)
#transGE112 = fourier.fourier_transform_loop(GE112, 5, 0.001, t)

fig, transforms = plt.subplots()
transforms.plot(t, np.power(np.abs(transGE11), 2)/(4*np.pi**2), label="GE11")
transforms.plot(t, np.power(np.abs(transGE12), 2)/(4*np.pi**2), label="GE12 (0)") 
#transforms.plot(t, np.power(np.abs(transGE112),2 )/(4*np.pi**2), label="GE12 (1)")
#transforms.plot(t, np.abs(transGE11), label="GE11")
#transforms.plot(t, np.abs(transGE12), label="GE12") 
transforms.grid()
"""
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
plt.show()