import mol
import matplotlib.pyplot as plt
import numpy as np


maxdistance = 10

m = mol.generate_matrix(10)

# Calculate the mutual information between a step and n steps away
t = np.arange(1.0,maxdistance+1,1.0)
plt.plot(t,mol.range_results_optimized(m))




plt.show()
