from qutip import *
import numpy as np
import scipy as scp
import matplotlib.pyplot as plt
# fourier contains necesary methods for Green's function analysis of 
# the qubit dynamics
#import fourier as fourier
# qubit contains necesary method fo qubit representation and operations, such
# as 2x2 representation, and object initialization. It also contains hard-
# coded functions for modeling the time evolution of the qubit, as were found
# symbolically (by hand) in qubit.pdf
import qubit as qubit

qt = qubit.qubit(np.sqrt(.5), -np.sqrt(.5), 1.0)

E = np.arange(-4, 4, 0.001)
t = np.arange(0, 10*np.pi, 0.01)

up = qt.GE11()
down = qt.GE12()

#ft0 = ftup(0)
#ft1 = ftup(1)
ftup = map(up, E)
ftdown = map(down, E)

def ft(ft, t, E):
	return 0.001 * np.dot(ft, np.exp(-1j*t*E))

print ft(ftup, 0, E)

#transup = scp.fftpack.fft(ftup)
#transdown = scp.fftpack.fft(ftdown)
transup = np.zeros(len(t))
transdown = np.zeros(len(t))
for i in range(0, len(t)):
	transup[i] = ft(ftup,t[i], E)
	transdown[i] = ft(ftdown,t[i], E)

fig, fourier = plt.subplots()
fourier.plot(t, np.abs(transup))
fourier.plot(t, np.abs(transdown))
fourier.grid()
plt.show()

"""
def fbase(x): return 0
def ftsum(ft):
	def ftsumh(t):
		#return lambda t: sum((ft)(t))
		return lambda t: np.sum(map(ft, t))
		#crt = fbase(t)
		#for i in range(0, len(ft)):
		#		crt = crt + ft[i](t)
		#return crt
	return ftsumh

#ftsum = lambda t: ft0(t) + ft1(t)

fsumup = ftsum(ftup)
fsumdown = ftsum(ftdown)
fig, fourier = plt.subplots()
fourier.plot(t, np.power(np.abs(0.001*fsumup(t)),2)/(4*np.pi**2))
fourier.plot(t, np.power(np.abs(0.001*fsumdown(t)),2)/(4*np.pi**2))
fourier.grid()
plt.show()

f = lambda func,E,t : ft_help( func, E, t)


def ft (func):
	def ft_help (t):
		return np.dot(map(func, E), np.exp(-1j*t*E))
	return ft_help

def help(E, t):
	return np.exp(-1j*t*E)
GE11t = ft(up)
GE12t = ft(down)
"""
"""
GE11 = map(up, E)
GE12 = map(down,E)
fig, energy1 = plt.subplots()
energy1.plot(E, GE11)
energy1.grid()
fig, energy2 = plt.subplots()
energy2.plot(E, GE12)
energy2.grid()
"""
"""
#ft_up = scp.fftpack.fft(GE11)
#ft_down = scp.fftpack.fft(GE12)
trans_up = map(GE11t, t)
trans_down = map(GE12t, t)

trans_up = trans_up*0.001
trans_down = trans_down*0.001

fig, fourier = plt.subplots()
fourier.plot(t, np.abs(trans_up))
fourier.plot(t, np.abs(trans_down))
fourier.grid()
plt.show
"""