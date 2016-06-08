import numpy as np
import sympy as sp
from functools import singledispatch

###### types to string #######
@singledispatch
def t2s(arg):
	return 'error'
@t2s.register(int)
def ri(arg):
	return 'i'
@t2s.register(float)
def rf(arg):
	return 'f'
@t2s.register(np.ndarray)
def ra(arg):
	return 'a'

###### initialize ############
n = 5
D = np.random.randn(n**2).reshape(n,n)
for i in np.arange(n):
	for j in np.arange(i):
		D[i][j] = D[j][i]

def culA():
	a = np.zeros(n**2).reshape(n,n)
	dr2 = dc2 = np.zeros(n)
	for i in np.arange(n):
		for j in np.arange(n):
			dr2[i] += D[j][i]**2
			dc2[i] += D[i][j]**2
	da2 = 0
	for i in np.arange(n):
		da2 += dr2[i]
	da2 = da2 / n**2
	dr2 = dr2 / 2
	dc2 = dc2 / 2

	for i in np.arange(n):
		for j in np.arange(n):
			a[i][j] = ((dc2[i]+dr2[j]) - da2 - D[i][j]*2)/2

	return a

A = culA()

L,X = np.linalg.eigh(A)
L.sort()
L = L[::-1]
L_pos = np.array([L[i] if L[i]>0 else 0 for i in range(n)])
d = np.count_nonzero(L_pos)	# d ... the number of positive values
Ln = np.sqrt(np.diag(L_pos))

P = X.dot(Ln)

f2 = np.array(L[0:d])
f2[::2] = 0
f1 = L[0:d] - f2
e1 = f1 / np.linalg.norm(f1)
e2 = f2 / np.linalg.norm(f2)

p0 = P[0,0:d]
x = np.dot(p0,e1)
y = np.dot(p0,e2)

######## upDating #############

# *args = (x,y,x2,y2,n,p,e1,e2)
def agi(*args):
	# input length matching
	if(len(args) != 7):
		raise Exception("input length error")

	# input type matching
	s = ""
	for i in np.arange(7):
		s = s + t2s(args[i])

	if(s != "ffffaaa"):
		raise Exception("input type error")

	# assign
	x,y,x2,y2,p,e1,e2 = args

	# Define return values
	a1,b1,c1,a2,b2,c2,t,s = sp.symbols('a1 b1 c1 a2 b2 c2 t s')
	E1 = a1*e1 + b1*e2 + c1*p
	E2 = a2*e1 + b2*e2 + c2*p
	r  = s*e1 + t*e2

	# satisfy construction
	ans = sp.solve([
		np.dot(E1,E1) - 1,
		np.dot(E2,E2) - 1,
		np.dot(E1,E2),
		np.dot(r,r) - 1,
		np.dot(E1,r) - np.dot(e1,r),
		np.dot(E2,r) - np.dot(e2,r),
		np.dot(p,E1) - x2,
		np.dot(p,E2) - y2
		],
		[a1,b1,c1,a2,b2,c2,t,s]) 

	# return the answer
	print(ans)
	return 0

answer = agi(x,y,1.0,1.0,p0,e1,e2)