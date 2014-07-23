"""
Read timeseries from DATA_PATH folder and retrieves the
eigenvalues of the transition matrix. Needs to be tested
and optimized.
"""

import numpy as np
import matplotlib.pyplot as p
import os

# this folder should exist (and timeseries be populated)
# note: there should be more timesteps in the timeseries
# than the number of molecules
DATA_PATH = 'timeseries/'
FIG_PATH = 'fig/'

# both should not be set true at the same time
AVERAGED = False
SENSITIVITY = True

for i, fname in enumerate(os.listdir(DATA_PATH)):
    data = np.loadtxt(DATA_PATH+fname, skiprows=1)
    nmols = data.shape[1] - 1
    eigenvalues = np.linspace(0, 0, nmols)

    # data = time, molA, molB...
    time = data[:, 0]
    mols = data[:, 1:]
    mat = np.eye(nmols)

    print "Filename: %s, nmols = %d, max_time = %d" % (fname, nmols, time[-1])
    assert time[-1] > nmols, 'max_time should be larger or equal to nmols'

    # testing purposes
    if AVERAGED:
        for k in range(16):
            ti = k*200
            tf = nmols + k*200
            A = mols[ti:tf, :]

            for j in range(nmols):
                b = data[ti:tf, j]
                x = np.linalg.solve(A, b)
                mat[j] = x

            eigenvalues += np.absolute(np.linalg.eig(mat)[0])
            avg_eigenvalues = eigenvalues/(k+1)
            eigen = len(np.where(avg_eigenvalues > 1.)[0])
            print "%02d eigenvalues > 1. after %d iterations" % (eigen, k)

    if SENSITIVITY:
        for k in range(16):
            ti = k*200
            tf = nmols + k*200

            # we just need as many timesteps as there are molecules
            A = mols[ti:tf, :]

            for j in range(nmols):
                b = data[ti:tf, j]
                x = np.linalg.solve(A, b)
                mat[j] = x

            eigenvalues = np.linalg.eig(mat)[0]
            eigen = len(np.where(np.absolute(eigenvalues) > 1.)[0])
            print "%02d eigenvalues > 1. " % (eigen)

            p.figure()
            p.plot(np.real(eigenvalues), np.imag(eigenvalues), 'bo')
            p.xlim(-3, 3)
            p.ylim(-3, 3)
            p.savefig(FIG_PATH+fname+'_'+str(ti)+'.png', bbox_inches='tight')


    # we just want to do the full analysis just yet
    if i == 0:
        break
