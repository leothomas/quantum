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

# plotting time evolution of qubit L probabilites
fig, fourier11 = plt.subplots()
fourier11.plot(t, np.abs(np.array(list(map(getByIndex(0,0), expect_valsL)))))
fourier11.plot(t, np.abs(np.array(list(map(getByIndex(1,1), expect_valsL)))))
fourier11.set_title("Expectation Values of Left Qubit (1,1) and (2,2)")
fourier11.grid()

fig, fourier22 = plt.subplots()
fourier22.plot(t, np.abs(np.array(list(map(getByIndex(0,0), expect_valsR)))))
fourier22.plot(t, np.abs(np.array(list(map(getByIndex(1,1), expect_valsR)))))
fourier22.set_title("Expectation Values of Right Qubit (1,1) and (2,2)")
fourier22.grid()

"""
The following code couples two qubits with an energy dependence,
finds the time evolution of the entangled system and uses the partial
trace operation to find dynamics specific to either sub-system. 
The accuracy of these operations can be veified by setting the energy 
dependence betweent the two qubits and verified that the dynamics 
returned by the partial trace operation correctly match those 
of the original, unentangled qubits in the above code.
"""

# u = 4x4 matrix of zeros, with energy dependence at location, 
# depending on system 
u = np.zeros(shape = (4,4), dtype = np.complex)
energyDep = 0.75
u[2][2] = energyDep
u[1][1] = -energyDep

# Hamiltonians of the coupled qubit system
qubitCoupled = couple([qubitR, qubitL], u)

# Green's energy function of coupled hamiltonian
GEC = np.array(list(map(GE(qubitCoupled), E)))

"""
# Debugging methods to print various steps of the operations, all return 
# the first "page" of their operation:

print ("\n Semi-tensored qubit: \n")
print (np.kron(np.eye(2), qubitR))
print ("\n")
print (np.kron(qubitL, np.eye(2)))


print ("\n qubitL: \n", qubitL)
print ("\n qubitR: \n", qubitR)
print ("\n tensored qubit: \n", qubitCoupled)
print ("\n Green's function: \n", GE(qubitCoupled)(0))
print ("\n tr_q2(q1 X q2) = qubitL: \n", ptrace(qubitCoupled, [2],[2,2]))
print ("\n tr_q1(q2 X q1) = qubitR:\n", ptrace(qubitCoupled, [1], [2,2]))
"""

#temporarily holds function being transformed over
tempGE= np.zeros(len(GEC), dtype = np.complex)

# fourrier transform of elements of the coupled system's Green's energy
# G(E) -> G(t)
transforms = fourier(GEC, E, dE, t)

# Column 1 of G(t) -> state vector psi(t), and density matrix rho(t) =
# |psi(t)><psi(t)|
rho = np.array(list(map(density_matrix, transforms)))

# arrays which will hold traced out qubit dynamics
recover_left = np.zeros(shape=(len(rho), 2, 2), dtype= np.complex)
recover_right = np.zeros(shape=(len(rho), 2, 2), dtype= np.complex)

# partial trace operation: for each "page" in rho(t)
# ptrace( entangled_system, [index_of_element(s)_ to_trace out], 
#			[dimensions_of_sub_systems])
for i in range(len(recover_left)):
	recover_left[i] = ptrace(rho[i], [2], [2,2])
	recover_right[i] = ptrace(rho[i], [1], [2,2])

# plot traced/recovered qubit dynamics
fig, fourier2 = plt.subplots()
fourier2.plot(t, np.abs(np.array(list(map(getByIndex(0,0), recover_left)))))
fourier2.plot(t, np.abs(np.array(list(map(getByIndex(1,1), recover_left)))))
fourier2.set_ylim([0,1])
fourier2.set_title("Traced out right (recover left)- (1,1) and (2,2)")
fourier2.grid()

fig, fourier4 = plt.subplots()
fourier4.plot(t, np.abs(np.array(list(map(getByIndex(0,0), recover_right)))))
fourier4.plot(t, np.abs(np.array(list(map(getByIndex(1,1), recover_right)))))
fourier4.set_ylim([0,1])
fourier4.set_title("Traced out left (recover right) - (1,1) and (1,2)")
fourier4.grid()

plt.show()
