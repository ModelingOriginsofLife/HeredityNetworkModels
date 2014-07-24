import numpy as np
import matplotlib.pyplot as p

data = np.loadtxt('transition.mat')
eig = np.linalg.eig(data)[0]

p.plot(np.real(eig), np.imag(eig), 'bo')
p.savefig('transition.png', bbox_inches='tight')
