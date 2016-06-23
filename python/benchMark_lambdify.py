import numpy as np
import sympy as sp
from sympy.utilities.lambdify import lambdify
from sympy import Matrix
from scipy import optimize as opt

edge = np.genfromtxt('csv/adjacency.csv', delimiter=",").astype(np.int16)
eN = len(edge)
P = np.genfromtxt('csv/mdSpace.csv', delimiter=",")
n = len(P)
L = np.genfromtxt('csv/eigVals.csv', delimiter=",")
L_pos = np.array([L[i] if L[i] > 0 else 0 for i in range(n)])
d = np.count_nonzero(L_pos)
Ln = np.sqrt(L_pos)
f2 = np.array(Ln[0:d])
f2[::2] = 0
f1 = Ln[0:d] - f2
e1 = (f1 / np.linalg.norm(f1)).reshape(d,1)
e2 = (f2 / np.linalg.norm(f2)).reshape(d,1)

print("init: ready")

# sympy
a1,b1,c1,a2,b2,c2,t,s = sp.symbols('a1 b1 c1 a2 b2 c2 t s')   # variables
x2_s,y2_s = sp.symbols('x2_s y2_s')  # values
P_i = sp.MatrixSymbol('P_i', d, 1)
E1 = sp.MatrixSymbol('E1', d, 1)
E2 = sp.MatrixSymbol('E2', d, 1)
var = (x2_s,y2_s,P_i,E1,E2,a1,b1,c1,a2,b2,c2,t,s)

_E1 = a1*sp.Matrix(E1) + b1*sp.Matrix(E2) + c1*sp.Matrix(P_i)
_E2 = a2*sp.Matrix(E1) + b2*sp.Matrix(E2) + c2*sp.Matrix(P_i)
R = s*sp.Matrix(E1) + t*sp.Matrix(E2)

f = Matrix([
		_E1.dot(_E1) - 1,
		_E2.dot(_E2) - 1,
		_E1.dot(_E2),
		R.dot(R) - 1,
		_E1.dot(R) - sp.Matrix(E1).dot(R),
		_E2.dot(R) - sp.Matrix(E2).dot(R),
		sp.Matrix(P_i).dot(_E1) - x2_s,
		sp.Matrix(P_i).dot(_E2) - y2_s
		])

func = sp.simplify(sp.Matrix.norm(f))

lam_f = lambdify(var, func, 'numpy')

def lam(x2, y2, p, e_1, e_2):
    return lambda a1,b1,c1,a2,b2,c2,t,s: \
        lam_f(x2, y2, p, e_1, e_2, a1, b1, c1, a2, b2, c2, t, s)

arr = np.array([1, 1, 1, 1, 1, 1, 1, 1])
print("lambda: ready")

X_sample = 3 * np.random.random_sample((100, 1)) - 1.5
Y_sample = 3 * np.random.random_sample((100, 1)) - 1.5

print("ready")

import time
start = time.time()

for i in range(100):
    global e1, e2
    temp1 = e1.reshape(d,1)
    temp2 = e2.reshape(d,1)
    f2 = lam(X_sample[i], Y_sample[i], P[14].reshape(d, 1), temp1, temp2)
    def g(args): return f2(*args)
    res = opt.minimize(g, arr, method='L-BFGS-B')
    if (res.success):
        e1 = res.x[0] * temp1 + res.x[1] * temp2 + res.x[2] * P[14].reshape(d, 1)
        e2 = res.x[3] * temp1 + res.x[4] * temp2 + res.x[5] * P[14].reshape(d, 1)

elapsed_time = time.time() - start
print("elapsed_time: " + str(elapsed_time) + "[sec]")