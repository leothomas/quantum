import numpy as np
import scipy as scp 
import matplotlib.pyplot as plt
from functions import *

# Define Hamiltonian operators for two qubits, left and right
qubitL = qubit(1, 1, 0.25)
qubitR = qubit(0.7, 0.7, 0.75)

# integration step 
dE = 0.001

# energy range to integrate fourier transformes over
E = np.arange(-5, 5, dE)

# time scale to output fourier transforms over
t = np.arange(0, 50, 0.1)

# GE (Green's energy function) maped over E(energy range). 
#Makes a 3d array, with "pages" corresponding to the evolution of GE. 
# GEL(0) = GE(E = -5), GEL(1) = GE(E = -5+dE), etc...
GEL = np.array(list(map(GE(qubitL), E)))
GER = np.array(list(map(GE(qubitR), E)))

# plot energy functions
fig, energy1 = plt.subplots()
energy1.plot(E, np.array(list(map(getByIndex(0,0), GEL))))
energy1.grid()

fig, energy2 = plt.subplots()
energy2.plot(E, np.array(list(map(getByIndex(0,0), GER))))
energy2.grid()

# fourrier transform (time evolution) of the Green's energy 
# function, evaluated over the range specified by t (time space)
# G(E) -> G(t)
transformsL = fourier(GEL, E, dE, t)
transformsR = fourier(GER, E, dE, t)

# takes first column from G(t) as state vector Psi(t) and returns
# time evolution of the expectation values, rho(t)=| Psi(t) >< Psi(t) |
expect_valsL = np.array(list(map(density_matrix, transformsL)))
expect_valsR = np.array(list(map(density_matrix, transformsR)))

# N = number of links in chain (simulation time increases exponentially
# with this)
N = 6
# u = energy dependence matrix 
u = np.zeros(shape=(4*N, 4*N), dtype = np.complex)
energyDep =0.1

# the following locations always determine the interaction points of 
# the spin_down energy of the left qubit and first element of the chain, 
# and the spin_up energy of the right qubit and last element of the chain
u[0][0] = energyDep
u[1][1] = energyDep
u[(2*N) -2][(2*N)-2] = energyDep
u[(4*N)-2][(4*N)-2] = energyDep

# defines a disorder ranges for the channel. Will be plotted as: 
# difference between entangled and un-entangled qubits' time 
# evolutions vs. disorder ranges. 
# The following defines 10 disorder ranges between 0 and 1, and 5 
# between 1 and 2.25
d1= np.arange(0, 1.0, 0.1)
d2 = np.arange(1.0, 2.25, 0.25)

disorder_range = np.append(d1,d2)

# to test a single disorder range, comment out the above line, and 
# uncomment the line below, changing 0 to range dersired 

#disorder_range = [0]

# --- IMPT: you will also have to comment out the functions
# that write the outputted data to .txt files (anything with "open", 
# "read", "write", "close", etc). These expect a specific 
# data format and writting difference ranges to those files will cause
# the graphing functions to crash

diffsL = open("diffsL.txt", 'a')
diffsR = open("diffsR.txt", 'a')

# arrays will hold the difference between the unentangled and entangled 
# expectation values' time evolutions
diff_left = np.zeros(len(disorder_range))
diff_right = np.zeros(len(disorder_range))


for j in np.arange(0, len(disorder_range)):
	
	# channel(length, type, energyBound): 
	# types: 0 --> 0 energies with unit coupling constants
	# 		 1 --> unit energies with random coupling constatns
	# anything else --> random energies with 1 coupling constants

	channel = channel(N, -1, disorder_range[j])
	
	#qubitChannel = combined Hamiltonains of both qubits and channel
	qubitChannel = couple([qubitL, channel, qubitR], u)

	# Green's energy for the entangled system
	GECh = np.array(list(map(GE(qubitChannel), E)))
	
	# Fourrier transforms of each element of the Greens' function:
	# G(E) -> G(t)
	transformsCh = fourier(GECh, E, dE, t)

	# time evolution of the density matrix of the entangled system
	rho_ch = np.array(list(map(density_matrix, transformsCh)))

	"""
	# Concurrence not yet useful/ accurate 

	concur_ch = np.array(list(map(concurence, rho_ch)))
	fig, concur1 = plt.subplots()
	concur1.plot(t, concur_ch)
	concur1.set_title("concurence of density matrix -> 0 for pure states")
	concur1.grid()
	"""

	# will hold exepctation values of recoverd qubits
	recover_left_ch = np.zeros(shape=(len(rho_ch), 2, 2), dtype= np.complex)
	recover_right_ch = np.zeros(shape=(len(rho_ch), 2, 2), dtype= np.complex)

	# partial trace operations 
	for i in range(len(recover_left_ch)):
		recover_left_ch[i] = ptrace(rho_ch[i], [1,2], [2,N,2])
		recover_right_ch[i] = ptrace(rho_ch[i], [2,3], [2,N,2])

	"""

	# Von Neuman entropy functions, not sure if works correctly 

	entropyL = np.array(list(map(entropy, recover_left_ch)))
	entropyR = np.array(list(map(entropy, recover_right_ch)))

	
	fig, entropy1 = plt.subplots()
	entropy1.plot(t, entropyL, label = 'left qubit')
	entropy1.plot(t, entropyR, label = 'right qubit')
	entropy1.set_title("Von neuman entropy of density matrix -> 0 for pure state")
	entropy1.grid()
	"""
		
	# Un-comment the following functions to see the time evolutions
	# of the entangled qubits, for each disorder specified by the 
	# disorder range

	"""
	titleLeft = "Traced out channel + right (recover left) - e range = " + str(disorder_range[j])
	titleRight = "Traced out channel + left (recover right) - e range = " + str(disorder_range[j])

	fig, fourier6 = plt.subplots()
	fourier6.plot(t, np.abs(np.array(list(map(getByIndex(0,0), recover_left_ch)))))
	fourier6.set_ylim([0,1])
	fourier6.set_title(titleLeft)
	fourier6.grid()

	fig, fourier7 = plt.subplots()
	fourier7.plot(t, np.abs(np.array(list(map(getByIndex(0,0), recover_right_ch)))))
	fourier7.plot(t, np.abs(np.array(list(map(getByIndex(1,1), recover_right_ch)))))
	fourier7.set_ylim([0,1])
	fourier7.set_title(titleRight)
	fourier7.grid()
	"""
	
	# returns the average difference between the entangled and unentangled 
	# qubits, since expectation values are a percentage values, their 
	# average difference is also a percentage value

	diff_left[j] = avg_diff(np.abs(np.array(list(map(getByIndex(0,0), recover_left_ch)))), \
		np.abs(np.array(list(map(getByIndex(0,0), expect_valsL)))))
	
	diff_right[j] = avg_diff(np.abs(np.array(list(map(getByIndex(0,0), recover_right_ch)))), \
		np.abs(np.array(list(map(getByIndex(0,0), expect_valsR)))))

	
	# writes these differences to the files
	diffsL.write('{:.12f}'.format(diff_left[j]) + " ")
	diffsR.write('{:.12f}'.format(diff_right[j]) + " ")

	# print the average differences
	print("\n Average difference (%) between expectation values of single qubit and \
		entangled qubit at disorder = ", disorder_range[j] )
	print ("\n left qubit (up): ", diff_left[j])
	print ("\n right qubit (up): ", diff_right[j])

print("\n dlu: \n",diff_left)
print("\n dru: \n",diff_right)

diffsL.write("\n")
diffsR.write("\n")

diffsL.close()
diffsR.close()

