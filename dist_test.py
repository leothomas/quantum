import numpy as np
import scipy as scp
import matplotlib.pyplot as plt

from functions import *

def val_test():

	crt_min = 5000
	crt_vals = [-2,-2,-2]
	distances = np.array([])
	vals = np.array([[]])

	for i in np.arange(-1.0, 1.0, 0.25):
		for j in np.arange(-1.0, 1.0, 0.25):
			for k in np.arange(0.0, 1.0, 0.25):
				
				qubitL = qubit(i, j, k)
				qubitR = qubit(i, j, k)

				# integration step 
				dE = 0.001

				# energy range to integrate fourier transformes over
				E = np.arange(-5, 5, dE)

				# time scale to output fourier transforms over
				t = np.arange(0, 50, 0.1)

				t = np.arange(0, 50, 0.1)

				# create a "fouriers" function that takes the fourier transform
				# of each element of a matrix, and returns a matrix of suspended eval
				# functions

				# GE (Green's energy function) maped over E. Make a 3d array, with "pages" 
				# corresponding to the E-evolution of GE. 
				# page(0) = GE(E = -5), page(1) = GE(E = -5+dE), etc...
				GEL = np.array(list(map(GE(qubitL), E)))
				#GER = np.array(list(map(GE(qubitR), E)))

				transforms = fourier(GEL, E, dE, t)

				expect_vals = np.array(list(map(density_matrix, transforms)))

				# plotting time evolution of qubit L probabilites

				expect_og_up = np.abs(np.array(list(map(getByIndex(0,0), expect_vals))))
				expect_og_down = np.abs(np.array(list(map(getByIndex(1,1), expect_vals))))

				# tensor product of qubitL and qubitR
				qubitCoupled = np.kron(qubitR, qubitL)
				# coupling factor between qubits - untested as of yet
				u = np.zeros(shape=(len(qubitCoupled), len(qubitCoupled[0])))
				couplingConst = 0.0
				u[2][2] = couplingConst

				qubitCoupled = qubitCoupled + u
				# GEC(HR tensor HL)
				GEC = np.array(list(map(GE(qubitCoupled), E)))

				tempGE= np.zeros(len(GEC))
				# transforms "book" same as GE "pages"
				#transforms = np.zeros(shape=(len(t), len(GEC[0]), len(GEC[0][0])))
				transforms = fourier(GEC, E, dE, t)
				# psi colummn vector will create a 4x4 (hopefully) matrix
				rho = np.array(list(map(density_matrix, transforms)))

				# trace out right qubit. This should match the fourier transform of
				# the original left qubit. It doesn't. Suprised? I'm not.
				#ptrace_vec = np.vectorize(ptrace, excluded=['traceout', 'dims'])
				recover_left = np.zeros(shape=(len(rho), 2, 2), dtype= np.complex)
				recover_right = np.zeros(shape=(len(rho), 2, 2), dtype= np.complex)

				for ind in range(len(recover_left)):
					recover_left[ind] = ptrace(rho[ind], [2], [2,2])
					recover_right[ind] = ptrace(rho[ind], [1], [2,2])

				expect_rec_up = np.abs(np.array(list(map(getByIndex(0,0), recover_left))))
				expect_rec_down = np.abs(np.array(list(map(getByIndex(0,0), recover_left))))

				dist_up = measure_dist(expect_og_up, expect_rec_up)
				dist_down = measure_dist(expect_og_down, expect_rec_down)

				np.append(distances, (dist_up+dist_down))
				np.append(vals, [i,j,k])

				if (dist_down + dist_up <= crt_min):
					crt_min = dist_down + dist_up
					crt_vals = [i,j,k]
					print ("\n crt min dist: ", crt_min )
					print ("\n crt minv vals: ", crt_vals)


	print ("min values: [Eup, Edown, T]: ", vals[distances.index(min(distances))])


val_test()