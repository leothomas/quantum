import numpy as np
import scipy as scp
import matplotlib.pyplot as plt
from qutip import *

b = Bloch()

point = [1/np.sqrt(3), 1/np.sqrt(3), 1/np.sqrt(3)]
b.add_vectors(point)
b.add_states(basis(2,0))
b.show()

"""
b3d = Bloch3d()

b3d.add_states(basis(2,0), basis(2,1))
b3d.show()
"""