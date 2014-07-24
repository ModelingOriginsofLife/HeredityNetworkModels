import math
import sympy as S
from sympy.core.symbol import Symbol
from sympy import Matrix, solve_linear_system
import numpy as NP
import numpy.random as random
import matplotlib.pyplot as plt

# We have a given amount of files representing various realisations of the system



### Reading a file that contains time series that should look like t mol1 mol2 ...
##for index in range(len(sys.argv) - 2):
##	file = sys.argv[i+1]  
##	data = np.loadtxt(file)
##
##	time = data[:,0]

###we are trying to build the transition matrix

##	for i in range(len(data[0])-1):
	

# Let say we already have the transition matrix instead, that is in a file as a n*n matrix
# separated by blanks

#data = NP.loadtxt(sys.argv[1])
# data[i][j] represents the probability to go from state i to state j
# Considering the amount of molecules we are dealing with, this might hurt.

# Compute the probability to be in a given state at equilibrium
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

def pij(i,j,mat,statics=[],steps=1):
	value = 0.0
	if(statics==[]):
		value = compute_proba(mat)[Symbol("p%i" % i)]
	else:
		value = statics[Symbol("p%i" % i)]
	size = len(mat)
	path = NP.matrix(mat)
	for index in range(steps -1):
		path *= NP.matrix(mat)
	return value*path.getA()[i][j]
			
# Compute the mutual information between state(n) and state(n+1)
def mutual_information(matrix,steps=1,statics=[]):
	value = 0.0
	if statics==[]:
		statics = compute_proba(matrix)
	for i in range(len(matrix)):
		for j in range(len(matrix[0])):
			p = pij(i,j,matrix,statics,steps)
			if p>0:
				value += p * math.log(p/(statics[Symbol("p%i" % i)]*statics[Symbol("p%i" % j)]))
	return value

# Other possible metric: compute 

#For tests, generate random transition matrix of size n
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
	statics = compute_proba(m)
	for i in range(stop-start+1):
		l.append(mutual_information(m,i+1,statics))
	return l

def range_results_optimized(m,start=1,stop=10):
	l = []
	statics = compute_proba(m)
	path = NP.matrix(m)
	for i in range(stop-start+1):
		value = 0.0
		for a in range(len(m)):
			for b in range(len(m[0])):
				p = statics[Symbol("p%i" % a)]*path.getA()[a][b]
				if p>0:
					value += p * math.log(p/(statics[Symbol("p%i" % a)]*statics[Symbol("p%i" % b)]))
		path *= NP.matrix(m)
		l.append(value)
	return l

#Show results fast
def plot_results(ar):
	plt.plot(ar)
	plt.show()

def normalize(mat):
	res = []
	for i in range(len(mat)):
		res.append(mat[i]/sum(mat[i]))
	return res
