import numpy as np
import scipy as scp 
import matplotlib.pyplot as plt
import qutip as qt

def qubit(Eup, Edown, Tau):
	return np.array([[Eup, Tau],[np.conj(Tau), Edown]])

qubitR = qubit(0.4, -0.6, 0.25)
qubitL = qubit(np.sqrt(.2), -np.sqrt(.8), 0.75)

def inverse_2x2(m):
	if (m[0][0]*m[1][1] != m[0][1]*m[1][0]):
		return 1/(m[0][0]*m[1][1]- m[0][1]*m[1][0]) * \
					np.array([[m[1][1], -m[0][1]], [-m[1][0], m[0][0]]] ) 
	else:
		print ("Error: Singular matrix. E: ", E)

def GE(qubit):
	def GEhelp(E):
		n = len(qubit)
		Greens = ( (E-10**(-10)*1j)*np.eye(n) - qubit)
		if len(Greens) == 2:
			return inverse_2x2(Greens)
		else: 
			return np.linalg.pinv(Greens)
	return GEhelp

dE = 0.001
E = np.arange(-5, 5, dE)
t = np.arange(0, 50, 0.01)

def byIndex(i, j):
	def bihelp(GE):
		return GE[i][j]
	return bihelp

GER = np.array(list(map(GE(qubitR), E)))

GE11 = np.array(list(map(byIndex(0,0), GER)))
GE12 = np.array(list(map(byIndex(0,1), GER)))

# 5 seconds for all mapings

def fourier_transform(func, E):
	def ft_help(t):
		return dE * np.sum(np.exp(-1j*E*t) * np.transpose(func))
	return ft_help
"""
GE11t = fourier_transform(GE11, E)
GE12t = fourier_transform(GE12, E)
trans_11 = np.array(list(map(GE11t,t)))
trans_12 = np.array(list(map(GE12t,t)))

fig, fourier = plt.subplots()
fourier.plot(t, np.abs(trans_11))
fourier.plot(t, np.abs(trans_12))
fourier.grid()
plt.show()
"""

def channel(length):
	energies= np.random.uniform(low=(-1.0), high=(1.0), size=length)
	#energies = np.ones(length)
	matrix = np.eye(length) * np.transpose(energies)
	for i in range(0, length-1):
		matrix[i][i+1] = 1.0 	
		matrix[i+1][i] = 1.0
	return matrix

def numElem(n): return len(n) * len(n[0])

qL = qt.Qobj(qubitL)
qR = qt.Qobj(qubitR)

qubitCoupled = qt.tensor(qL, qR)
print "\n", qL
print "\n", qR
print "\n", qubitCoupled
print "\n", qubitCoupled.ptrace(0)
print "\n", qubitCoupled.ptrace(1)



"""
qubitR = qubit(0.4, -0.6, 0.25)
qubitL = qubit(np.sqrt(.2), -np.sqrt(.8), 0.75)

***PARTIAL TRACES RETURNED ZERO MATRIX WHEN EUP = -EDOWN***

Quantum object: dims = [[2], [2]], shape = [2, 2], type = oper, isherm = True
Qobj data =
[[ 0.4472136   0.75      ]
 [ 0.75       -0.89442719]]

Quantum object: dims = [[2], [2]], shape = [2, 2], type = oper, isherm = True
Qobj data =
[[ 0.4   0.25]
 [ 0.25 -0.6 ]]

Quantum object: dims = [[2, 2], [2, 2]], shape = [4, 4], type = oper, isherm = True
Qobj data =
[[ 0.17888544  0.1118034   0.3         0.1875    ]
 [ 0.1118034  -0.26832816  0.1875     -0.45      ]
 [ 0.3         0.1875     -0.35777088 -0.2236068 ]
 [ 0.1875     -0.45       -0.2236068   0.53665631]]

Quantum object: dims = [[2], [2]], shape = [2, 2], type = oper, isherm = True
Qobj data =
[[-0.08944272 -0.15      ]
 [-0.15        0.17888544]]

Quantum object: dims = [[2], [2]], shape = [2, 2], type = oper, isherm = True
Qobj data =
[[-0.17888544 -0.1118034 ]
 [-0.1118034   0.26832816]]
 
"""



"""
qubitCoupled = np.kron(qubitL, qubitR)

print (qubitL)
print (qubitR)
print (qubitCoupled)
"""


GEC = np.array(list(map(GE(qubitCoupled.full()), E)))

tempGE= np.zeros(len(GEC))

#for i in range(len(GEC)):
	#GEC[i].ptrace(0)
	#GEC[i].ptrace(1)


for j in range(len(GEC[0])):
	for k in range(len(GEC[0][0])):
		for i in range(len(GEC)):
			tempGE[i] = GEC[i][j][k]

		fig, temp = plt.subplots()
		tempGEt = fourier_transform(tempGE, E)
		transform = np.array(list(map(tempGEt,t)))
		temp.plot(t, np.abs(transform))
		temp.grid()

plt.show()