import math
import sympy as S
from sympy.core.symbol import Symbol
from sympy import Matrix, Float, solve_linear_system, solve
import numpy as NP
import numpy.random as random
import matplotlib.pyplot as plt

### This files implements useful functions to compute the mutual information (Shanon and Weaver 1949)
### of a system at two different moments. It is also possible to get the mutual information 
### as a function of the separtion between two states.

### However, we are working with "stationary" distribution. If need be, we normalize by the
### total population at a given step.

### Other functions could be implemented to try to track the amount of information transmitted
### from one time step to the next, such as the transfer entropy (Schreiber 2000)

### One reason to consider other metrics is that it seems the "amount" of mutual information is
### (roughly) the same no mater how many cliques (or attractors in a more general system) we have.
### The rational would be that we only measure how independent the system is with itself at a later
### time. The ways in which it stays correlated is not taken into account.

# Determine if the matrix is normalized
def isNormalized(m):
	for line in m:
		if NP.round(sum(line),10) != 1.0:
			return False
	return True

# Compute the probability to be in a given state at equilibrium
# In a non-normalized case, this is equivalent to the probability
# of picking a member of a specific species once the linear regime
# is reached.
def compute_proba(matrix):
	eqs = []
	syms = []
	lastline = []
	for i in range(len(matrix)-1): #we want to avoid overdetermination
		line = []
		for j in range(len(matrix)):
			#print "value for i: %i" % i
			if (i==j):
				line.append(matrix[j][i]-1.0)
			else:
				line.append(matrix[j][i])
		# We are summing in the form a1 p1 + ... + an pn = 0
		line.append(0)
		syms.append(Symbol("p%i" % i))
		eqs.append(line)
		lastline.append(1)
	syms.append(Symbol("p%i" % (len(matrix)-1)))
	lastline.append(1)
	lastline.append(1)
	eqs.append(lastline)
	#print syms
	#print eqs	
	return solve_linear_system(Matrix(eqs),*syms)

# In the case the matrix was not normalized. Requires sympy ver 0.7.4 or later
# Note that we are returning only the first possible stationary distribution
# There might be more than one (eg if we have disconnected parts)
def compute_stationary(matrix):
	eqs = []
	syms = []
	totalpop = 0
	for i in range(len(matrix)):
		syms.append(Symbol("p%i" % i))
	for i in range(len(matrix) - 1): #we don't overdetermine either
		val = 0
		for j in range(len(matrix)):
			val += matrix[j][i] * syms[j]
		eqs.append(val)
		totalpop += val
	#we need the last line only for the totalpop
	val = 0	
	for i in range(len(matrix)):
		val += matrix[i][len(matrix) - 1] * syms[i]
	totalpop += val
	for i in range(len(matrix)-1):
		eqs[i] -= totalpop*syms[i]
	eqs.append(sum(syms)-1.0)
	res = solve(eqs,*syms,dict=True)
	return removeNegative(res)[0]

# Non-normalized matrixes tend to have solutions with negative concentrations
# We don't take those into consideration
def removeNegative(res):
	newres = res
	for sol in res:
		for value in sol:
			if (sol[value] < 0.0):
				newres.remove(sol)
				break
	return newres

# Probability to pick a molecule of type i in the system and of type j steps timesteps
# later. If the system has no structure, this should be equal to the probability
# of type i time type j
def pij(i,j,mat,statics=[],steps=1):
	value = 0.0
	if(statics==[]):
		proba = compute_proba(mat) if isNormalized(mat) else compute_stationary(mat)
		value = proba[Symbol("p%i" % i)]
	else:
		value = statics[Symbol("p%i" % i)]
	size = len(mat)
	path = NP.matrix(mat)
	for index in range(steps -1):
		path *= NP.matrix(mat)
	return value*path.getA()[i][j]
			
# Compute the mutual information between state(n) and state(n+steps)
def mutual_information(matrix,steps=1,statics=[]):
	value = 0.0
	if statics==[]:
		statics = compute_proba(matrix) if isNormalized(matrix) else compute_stationary(matrix)
	popincrease = newpopsize(statics,matrix)
	for i in range(len(matrix)):
		for j in range(len(matrix[0])):
			p = pij(i,j,matrix,statics,steps)/(popincrease**steps) # This normalization is only valid
									       # because we are at the stationary
									       # distribution
			if p>0:
				value += p * math.log(p/(statics[Symbol("p%i" % i)]*statics[Symbol("p%i" % j)]))
	return value 

#For tests, generate random normalized transition matrix of size n
def generate_matrix(size=10):
	mat = []
	for i in range(size):
		line = random.ranf(size)
		line = line/sum(line)
		mat.append(line)
	return mat

#Generate a range of results
def range_results(m,start=1,stop=10):
	l = []
	statics = compute_proba(m) if isNormalized(m) else compute_stationary(m)
	for i in range(stop-start+1):
		l.append(mutual_information(m,i+1,statics))
	return l

# Same as before, but we reimplement the algo to be able to reuse intermediate computations
def range_results_optimized(m,start=1,stop=10):
	l = []
	statics = compute_proba(m) if isNormalized(m) else compute_stationary(m)
	path = NP.matrix(m)
	popincrease = newpopsize(statics,matrix)
	for i in range(stop-start+1):
		value = 0.0
		for a in range(len(m)):
			for b in range(len(m[0])):
				p = statics[Symbol("p%i" % a)]*path.getA()[a][b]/popincrease
				if p>0:
					value += p * math.log(p/(statics[Symbol("p%i" % a)]*statics[Symbol("p%i" % b)]))
		path *= NP.matrix(m)
		l.append(value)
	return l

# Convenience function to show results. If any specific representation is agreed upon in the future,
# it should go in there
def plot_results(ar):
	plt.plot(ar)
	plt.show()

# For tests, mainly. While "normalizing" (that is making the sum of the weights leaving a specific node
# equal to 1) change the static distribution, it should keep properties like the presence of cliques.
# There is obviously a problem if a node has no transition starting from it.
def normalize(mat):
	res = []
	for i in range(len(mat)):
		res.append(mat[i]/sum(mat[i]))
	return res

# Used to compute the increase of population over time
def newpopsize(pop,mat):
	if(len(pop) != len(mat)):
		raise(Exception("Population vector and transition matrix should be the same size"))
	popsize = 0.0	
	for i in range(len(pop)):
		for j in range(len(mat)):
			popsize += pop[j]*mat[j][i]
	return popsize
