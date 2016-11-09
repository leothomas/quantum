import numpy as np
import scipy as scp 
import matplotlib.pyplot as plt

# creates qubit
def qubit(Eup, Edown, Tau):
	return np.array([[Eup, Tau],[np.conj(Tau), Edown]])

# hard coded inverse of 2x2 matrix is faster to do this way 
def inverse_2x2(m):
	if (m[0][0]*m[1][1] != m[0][1]*m[1][0]):
		return 1/(m[0][0]*m[1][1]- m[0][1]*m[1][0]) * \
					np.array([[m[1][1], -m[0][1]], [-m[1][0], m[0][0]]] ) 
	else:
		print ("Error: Singular matrix. E: ", E)
# GE= inverse(E-H)
# subtract (i*10^-10) to avoid singularity
# pinv => pseudo-inverse. Uses least squares regression to solve 
# inv(A) * I = A
def GE(qubit):
	def GEhelp(E):
		n = len(qubit)
		Greens = ( (E-10**(-10)*1j)*np.eye(n) - qubit)
		if len(Greens) == 2:
			return inverse_2x2(Greens)
		else: 
			return np.linalg.pinv(Greens)
	return GEhelp

# obtaining element at ij position of each "page" of fourier transforms 
def getByIndex(i, j):
	def biHelp(n):
		return n[i][j]
	return biHelp
# fourier transform takes an array of GE evaluated over range E
# returns function of t, can be evaluated at later time 
def fourier_transform(func, E):
	def ft_help(t):
		return dE * np.sum(np.exp(-1j*E*t) * np.transpose(func))
	return ft_help

# creates a channel of uniformly distributed random values between -1.0 and 1.0
# along the main diagonal and coupling constants of 1.0 everywhere
def channel(length):
	energies= np.random.uniform(low=(-1.0), high=(1.0), size=length)
	#energies = np.ones(length)
	matrix = np.eye(length) * np.transpose(energies)
	for i in range(0, length-1):
		matrix[i][i+1] = 1.0 	
		matrix[i+1][i] = 1.0
	return matrix

# hard coded partial trace operation. Need help here for taking the partial
# trace of qubits tensored with channel -> confirmed to function with GE's and 
# fourier transformes at time t = 0
# x denotes qubit to trace out (0 = Left, 1 = Right)
def ptrace2q(x):
	def ptraceHelp(n):
		if (x == 1 ):
			return np.array([[n[0][0]+n[1][1],n[0][2]+n[1][3]],\
				[n[2][0]+n[3][1],n[2][2]+n[3][3]]])
		elif (x == 0 ):
			return np.array([[n[0][0]+n[2][2],n[0][1]+n[2][3]],\
				[n[1][0]+n[3][2],n[1][1]+n[3][3]]])
	return ptraceHelp

# testing begins here:

# define two qubits, if Eup + Edown != 1, the partial trace returned is 
# scaled by a factor of Eup + Edown. 
qubitR = qubit(0.6, -0.4, 0.19)
qubitL = qubit(-0.2, 0.8, 1.0)

# integration step 
dE = 0.001
# energy range to integrate fourier transformes over
E = np.arange(-5, 5, dE)
# time scale to output fourier transforms over
t = np.arange(0, 50, 0.01)

# GE(qubitL) maped over E. Make a 3d array, with "pages" corresponding 
# to the E-evolution of GE. page(0) = GE(E = -5), page(1) = GE(E = -5+dE)
GEL = np.array(list(map(GE(qubitL), E)))
GER = np.array(list(map(GE(qubitR), E)))

# fourier transforms of left qubit
# returns a function of t
GE11t = fourier_transform(np.array(list(map(getByIndex(0,0), GEL))), E)
GE12t = fourier_transform(np.array(list(map(getByIndex(0,1), GEL))), E)

# evaluating fourier transforms over the desired values of t
trans_11 = np.array(list(map(GE11t,t)))
trans_12 = np.array(list(map(GE12t,t)))

# plotting time evolution of qubit L probabilites
fig, fourier = plt.subplots()
fourier.plot(t, np.abs(trans_11)/(4*np.pi**2))
fourier.plot(t, np.abs(trans_12)/(4*np.pi**2))
fourier.set_title("Expectation Values of Right Qubit")
fourier.grid()

# tensor product of qubitL and qubitR
qubitCoupled = np.kron(qubitL, qubitR)
print ("\n qubitL: \n", qubitL)
print ("\n qubitR: \n", qubitR)
print ("\n tensored qubit: \n", qubitCoupled)
print ("\n tr_R(rho) = qubitL: \n", ptrace2q(1)(qubitCoupled))
print ("\n tr_L(rho) = qubitR:\n", ptrace2q(0)(qubitCoupled))

# fourier transform of GE(Left) and GE(right)
# to be compared with ptraces of fourier transform (GE (tensor(L,R))
fourierL = np.zeros(shape=(len(GEL[0]), len(GEL[0][0])))
for i in range(len(GEL[0])):
	for j in range(len(GEL[0][0])):
		func = np.array(list(map(getByIndex(i,j), GEL)))
		fourierL[i][j] = fourier_transform(func, E)(0)

print ("\n fourier transform of qubitL at t =0: \n", fourierL)

fourierR = np.zeros(shape=(len(GER[0]), len(GER[0][0])))
for i in range(len(GER[0])):
	for j in range(len(GER[0][0])):
		func = np.array(list(map(getByIndex(i,j), GER)))
		fourierR[i][j] = fourier_transform(func, E)(0)

print ("\n fourier transform of qubitR at t=0: \n", fourierR)
 
# coupling factor between qubits - untested as of yet
u = np.zeros(shape=(len(qubitCoupled), len(qubitCoupled[0])))
couplingConst = 0
u[2][2] = couplingConst

qubitCoupled = qubitCoupled + u

# GEC(HR tensor HL)
GEC = np.array(list(map(GE(qubitCoupled), E)))

# fourier transforms of GE(coupled) at t =0
fourierCoupled = np.zeros(shape=(len(GEC[0]),len(GEC[0][0])))
for i in range(len(GEC[0])):
	for j in range(len(GEC[0][0])):
		func = np.array(list(map(getByIndex(i,j),GEC)))
		fourierCoupled[i][j]= fourier_transform(func, E)(0)

print ("\n fourier transforms of coupledQubit at t=0: \n", fourierCoupled)
print ("\n traced out right of fourier transforms of coupledQubit: ")
print ("\n ",ptrace2q(1)(fourierCoupled))
print ("\n traced out left of fourier transforms of Coupled Qubit: ")
print ("\n ", ptrace2q(0)(fourierCoupled))

# transforms "book" same as GE "pages"
transforms = np.zeros(shape=(len(t), len(GEC[0]), len(GEC[0][0])))

# TODO: find out how to do this using mapping
#def insertByIndex()

def insert(n, a, i, j):
	for k in range(len(a)):
		n[k][i][j] = a[k]
	return n

#temporarily holds function being transformed over
tempGE= np.zeros(len(GEC))

for i in range(len(GEC[0])):
	for j in range(len(GEC[0][0])):
		# gets GE to transform over
		tempGE = np.array(list(map(getByIndex(i,j), GEC)))
		
		# apply suspended evaluation of fourier transforms
		tempGEt= fourier_transform(tempGE, E)

		#evaluate transform at all times in t		
		trans= np.array(list(map(tempGEt,t)))
		

		#adding transform to "book"
		transforms = insert(transforms,trans, i,j)

		#figure plotting methods: 

		#fig, temp = plt.subplots()
		#temp.plot(t, np.abs(trans))
		#title = ("transforms( " + str(i) + ", " + str(j) + " )")
		#temp.set_title(title)

# trace out right qubit. This should match the fourier transform of
# the original left qubit. It doesn't. Suprised? I'm not.
trace_left = np.array(list(map(ptrace2q(1), transforms)))

fig, fourier2 = plt.subplots()
fourier2.plot(t, np.abs(np.array(list(map(getByIndex(0,0), trace_left)))))
fourier2.plot(t, np.abs(np.array(list(map(getByIndex(0,1), trace_left)))))
fourier2.set_title("Traced out left qubit")
fourier2.grid()

fig, fourier3 = plt.subplots()
fourier3.plot(t, np.abs(np.array(list(map(getByIndex(1,0), trace_left)))))
fourier3.plot(t, np.abs(np.array(list(map(getByIndex(1,1), trace_left)))))
fourier3.set_title("Traced out left qubit")
fourier3.grid()

"""
for i in range(len(trace_left[0])):
	for j in range(len(trace_left[0][0])):
		ij = getByIndex(i,j)
		trans = np.array(list(map(ij, trace_left)))
		fig, temp = plt.subplots()
		temp.plot(t, np.abs(trans))
		temp.grid()
"""
"""
for i in range(len(GEC[0])):
	for j in range(len(GEC[0][0])):
		tempGE = getByIndex(GEC,i,j)
		#tempGE = GEC[:][i][j]
		tempGEt= fourier_transform(tempGE, E)		
		trans= np.array(list(map(tempGEt,t)))
		fig, temp = plt.subplots()
		transforms = insert(transforms,trans, i,j)
		#transforms[:][i][j]=np.array(list(map(tempGEt,t)))
		
"""
for i in range(len(transforms[0])):
	for j in range(len(transforms[0][0])):
		fig, temp = plt.subplots()
		temp.plot(t, np.abs(np.array(list(map(getByIndex(i,j), transforms)))))
		title = ("transforms( " + str(i) + ", " + str(j) + " )")
		temp.set_title(title)
		temp.grid()

plt.show()