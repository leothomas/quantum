
import numpy as np
import scipy as scp
import matplotlib.pyplot as plt

# creates qubit Hamiltonian
def qubit(Eup, Edown, Tau):
	return np.array([[Eup, Tau],[np.conj(Tau), Edown]])

# hard coded inverse of 2x2 matrix (faster) 
def inverse_2x2(m):
	if (m[0][0]*m[1][1] != m[0][1]*m[1][0]):
		return 1/(m[0][0]*m[1][1]- m[0][1]*m[1][0]) * \
					np.array([[m[1][1], -m[0][1]], [-m[1][0], m[0][0]]] ) 
	else:
		print ("Error: Singular matrix. E: ", E)

# GE= inverse(E-H)
# subtract (i*10^-10) to avoid singularity
# pinv => pseudo-inverse. Uses least squares regression to solve 
# inv(A) * I = A for any Hamiltonian greater than 2x2
def GE(qubit):
	def GEhelp(E):
		n = len(qubit)
		Greens = ( (E-(10**(-7)*1j))*np.eye(n) - qubit)
		if len(Greens) == 2:
			return inverse_2x2(Greens)
		else: 
			return np.linalg.pinv(Greens)
	return GEhelp

# combines 2 or 3 system Hamiltonians, with energy dependences 
# specified by U, which should be of same size as the combined 
# Hamiltonian
def couple(Hlist, U):
	if (len(Hlist) ==2):
		return (np.kron(Hlist[0], np.eye(len(Hlist[1])) ) + np.kron(np.eye(len(Hlist[0])), Hlist[1]) + U)
	elif (len(Hlist) == 3):
		return (np.kron(Hlist[0], np.kron(np.eye(len(Hlist[1])), np.eye(len(Hlist[2])))) + \
				np.kron(np.eye(len(Hlist[0])), np.kron(Hlist[1], np.eye(len(Hlist[2])))) + \
				np.kron(np.kron(np.eye(len(Hlist[0])), np.eye(len(Hlist[1]))), Hlist[2]) + U)

# Returns the fourrier transform of each element of G(E), evaluated over
# interval specified by t.
def fourier(GE, E, dE, t):
	transforms = np.zeros(shape=(len(t), len(GE[0]), len(GE[0][0])))
	func = np.zeros(len(GE))
	for i in range(len(GE[0])):
		for j in range(len(GE[0][0])):
			func = np.array(list(map(getByIndex(i,j), GE)))
			transform = fourier_transform(func, E, dE)
			transforms = insert(transforms, np.array(list(map(transform, t))), i, j)
	return transforms


# fourier transform takes an array of GE evaluated over range E
# returns function of t, can be evaluated at later time 
def fourier_transform(func, E, dE):
	def ft_help(t):
		return dE * np.sum(np.exp(-1j*E*t) * np.transpose(func))
	return ft_help

def normalize(x):
	#norm = np.sqrt(sum((i**2) for i in x))
	norm = np.linalg.norm(x)
	return x/norm

# returns normalized density matrix, given time evolution of Green's
# function 
def density_matrix(n):
	# get first column of n
	psi = n[0]
	norm = np.linalg.norm(psi)
	psi = psi/norm
	return np.multiply(np.transpose(np.conj(psi)), psi.reshape(len(psi),1))

# hard coded partial trace operation. Need help here for taking the partial
# trace of qubits tensored with channel -> confirmed to function with GE's and 
# fourier transformes at time t = 0
# x denotes qubit to trace out (0 = Left, 1 = Right)
def ptrace2q(x):
	def ptraceHelp(n):
		if (x == 1 ):
			return np.array([[n[0][0]+n[1][1], n[0][2]+n[1][3]], \
				[n[2][0]+n[3][1], n[2][2]+n[3][3]]])
		elif (x == 0 ):
			return np.array([[n[0][0]+n[2][2],n[0][1]+n[2][3]],\
				[n[1][0]+n[3][2],n[1][1]+n[3][3]]])
	return ptraceHelp

#takes partial trace of matrix p
# traceout = subsystems to trace out
# dimensions of tensor space
# p = matrix 
# traceout = index of syb-system to traceout
# dims = array of dimensions 
# traceout the first of 2 qubits coupled to a length five channel:
# ptrace(rho, 1, [2,5,2])

# This is converted from matlab code:
# https://github.com/CoryGroup/quantum-utils-matlab/blob/master/src/ptrace.m
def ptrace(p, traceout, dims):

	"""
	NOTE: ptrace(2) of q1 X (channel X q2) ==
		 ptrace(2) of ptrace(2) of (q1 X q2) X channel

	"""
	dims1 = [i+1 for i, x in enumerate(dims) if x ==1]
	
	traceout = np.setdiff1d(traceout, dims1)
	dims = np.array([x for x in dims if x !=1])
	n = len(dims)
	rdims = dims[::-1]
	keep = [x+1 for x in range(n)]
	keep = np.setdiff1d(keep, traceout)
	dimtrace = np.prod([dims[x-1] for x in traceout])
	dimkeep = len(p)/dimtrace

	"""
	# TODO: incorporate this check if p is a state vectors
	# 		and not a density matrix. For now assuming p is 
	#		always the result of a tensor product, we should 
	# 		have no problems. I, on the other hand, have a lot of problems

	if (len(p) == 1 or len(p[0]) ==1):
	
		perm = n+1 -(np.append((keep[::-1]),traceout)) 
		# print 
		x = np.reshape(np.transpose(np.reshape(p,rdims),perm), np.append[dimkeep,dimtrace])
		x = x*np.transpose(x)
	
		print (x)
	
	else:
		"""
	
	perm  = n + 1 - np.append(keep[::-1], np.append(keep[::-1]-n, np.append(traceout,traceout-n)))
	#t2 = np.transpose(np.transpose(t1, perm-1), [0,1,3,2])
	reshapedims = np.append(dimkeep,np.append(dimkeep,dimtrace**2))
	temp = np.reshape(np.transpose(np.reshape(np.transpose(p),np.append(rdims,rdims)), perm-1), reshapedims)
	# create this array: [1:(dimtrace+1):dimtrace**2]
	arr= np.arange(1, (dimtrace**2)+1, dimtrace+1)
	x = np.sum(temp[:,:,arr-1],2);
	
	return x


# creates a channel following the specifications of e: 
# e == 0, energies = 0, and coupling is unit
# e == 1, energies = 1, and coupling is random
# e == anyhting else, energies are uniformly, randomly distributed 
#	and coupling constants are unit
def channel(length, e, energyBound):
	
	if e ==0: 
		energies = np.zeros(length)
		couplings = 0.25*np.ones(length-1)
	elif e == 1:
		energies = np.ones(length)
		#couplings = [1,2,3]
		couplings = np.random.uniform(low=0, high=(energyBound), size=length-1)
	else:
		energies= np.random.uniform(low=(-energyBound), high=(energyBound), size=length)
		couplings = 0.25* np.ones(length-1)
		#couplings = np.random.uniform(low=0, high=(energyBound), size = length-1)

	matrix = np.eye(length) * np.transpose(energies)
	matrix  = matrix + 0j
	for i in range(0, length-1):
		matrix[i][i+1] = couplings[i]
		matrix[i+1][i] = couplings[i]

	#print ("\n channel: \n", matrix)
	return matrix

# obtain element at position i,j, from "page" n of a 3D matrix: 
# ie: for i in transforms: getByIndex(1,1)(i)
# would return the elements at 1,1 of each "page" in transforms
def getByIndex(i, j):
	def biHelp(n):
		return n[i][j]
	return biHelp

# returns ith column of "page" n of a 3D matrix. See above.
def getByColumn(i):
	def bcHelp(n):
		a = np.zeros(len(n), dtype = np.complex)
		for k in range(len(n)):
			a[k] = n[k][i]
		return a
	return bcHelp


# TODO: find out how to do this using mapping
#def insertByIndex()

# inserts array a in 3D array at position i,j
# ie: a[k] is now the element at position i,j of the kth page of n
def insert(n, a, i, j):
	n = n+0j
	for k in range(len(a)):
		n[k][i][j] = a[k]
	return n

# average difference between two arrays, element-wise
def avg_diff(a, b):
	return (np.mean(np.abs(np.subtract(a,b))))

# concurrence  and von neuman entropy --> not sure if very usefull 
# measures of entanglement dynamics as of right now
def concurence(a):
	eigvals = np.sort(np.linalg.eigvals(a))
	first = eigvals[0]
	tempsum = 0
	for i in range(1, len(eigvals)-1):
		tempsum = tempsum + eigvals[i]
	if (first- tempsum >0):
		return first-tempsum
	else:
		return 0

def entropy(rho):
	return -np.trace(rho * np.log2(rho))