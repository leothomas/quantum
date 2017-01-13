import numpy as np
import scipy as scp
from scipy import optimize
import matplotlib.pyplot as plt
from functions import *

"""

# Following code defines a range over which to graph the best fit
# line. However, its not yet implemented because the function to 
# which the data should be fitted is not yet clearly defined. 


d1= np.arange(0, 1.0, 0.1)
d2 = np.arange(1.0, 2.25, 0.25)
#disorder_range = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.25]
disorder_range = np.append(d1,d2)

e = disorder_range
eFit = np.append([-0.1], e)
eFit = np.append(e, [2.35])
"""

# load values from text fiels
dr = np.loadtxt('diffsR.txt')
dl = np.loadtxt('diffsL.txt')

# transpose to have rows corresponding to each disorder value, 
# as opposed to rows corresponding to each trial, with one 
# value for each disorder.
dr = np.transpose(dr)

# mean value and standard deviation for each disorder values
# (right qubit)
yr = [np.mean(i) for i in dr]
#yre =[np.std(i)/np.sqrt(len(i)) for i in dr]
yre = [np.std(i) for i in dr]
# (left qubit)
dl = np.transpose(dl)
yl = [np.mean(i) for i in dl]
#yle = [np.std(i)/np.sqrt(len(i)) for i in dl]
yle = [np.std(i) for i in dl]


fig, diffR = plt.subplots()

diffR.errorbar(e, yr, yre, linestyle='None', marker = '.')

diffR.set_xlim([-0.10, 2.35])
diffR.set_ylim([0.0, 0.4])
diffR.set_title("Avg. difference in right qubit expectation values")



# Function used to fit the average variance values as a function of 
# disorder

def func(x,a,b):
	return a*np.exp(b*x)

"""

# Fitting and plotting data. This is not implemented since the function 
# to which the data should be fitted is not clear.

poptR, pcov = scp.optimize.curve_fit(func, e, yr, maxfev=1200 )
fitHandle = ("f(x) = " + str('{:.3f}'.format(poptR[0])) + " e ^(" + str('{:.3f}'.format(poptR[1])) + "x)")
drFit = func(eFit, poptR[0], poptR[1])

# Plot fit line
fitLineR = diffR.plot(eFit, drFit, '-', label =fitHandle)

"""

# Notice variance decreasing as a function of energy distribution
# since as the disorder range allows for energies greater than the 
# channel's bandwidth, the probability of having a "closed" channel 
# increases exponentially and channel fidelity decreases. 

fig, errors = plt.subplots()
errors.scatter(e, yre)
errors.set_title("Variance as a function of energy distribution")
errors.grid()

#diffR.axvline(x=0.5, color='red', linestyle='dashed')
diffR.plot([0.5, 0.5], [0.0, 0.27], color = 'red', linestyle = 'dashed')
diffR.plot([0.5, 0.5], [0.32, 0.4], color = 'red', linestyle = 'dashed')
diffR.grid()

# Fitting and plotting the fit line for the variance as a function of
# energy distribution

poptL, pcov = scp.optimize.curve_fit(func, e, yl)
dlFit = func(eFit, poptL[0], poptL[1])
fitHandle = ("f(x) = " + str('{:.3f}'.format(poptL[0])) + "/ x^(" + str('{:.3f}'.format(poptL[1])) + ")")

# plotting differences for the left qubit. 

fig, diffL = plt.subplots()

diffL.errorbar(e, yl, yle, linestyle='None', marker = '.')
fitLineL = diffL.plot(eFit, dlFit, '-', label = fitHandle)
diffL.grid()
diffL.set_xlim([-0.1,2.35 ])
diffL.set_title("Avg. difference in left qubit expectation values")
plt.legend(handles = [fitLineL[0]])
plt.show()
 	