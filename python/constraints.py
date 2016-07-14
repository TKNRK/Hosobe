import numpy as np
import sympy as sp
from scipy import optimize as opt
from matplotlib import pyplot as plt
from sympy.utilities.lambdify import lambdify
from sympy import *

#  initialize
edge = np.genfromtxt('csv/adjacency.csv', delimiter=",").astype(np.int64)
P = np.genfromtxt('csv/mdSpace.csv', delimiter=",")
n = len(P)
L = np.genfromtxt('csv/eigVals.csv', delimiter=",")
L_pos = np.array([L[i] if L[i] > 0 else 0 for i in range(n)])
d = np.count_nonzero(L_pos)
Ln = np.sqrt(L_pos)

f1 = np.array(Ln[0:d])
f2 = np.array(Ln[0:d])
f3 = np.array(Ln[0:d])
f1[0::3] = 0 ; f1[1::3] = 0
f2[1::3] = 0 ; f2[2::3] = 0
f3[2::3] = 0 ; f3[0::3] = 0
e1 = (f1 / np.linalg.norm(f1)).reshape(d,1)
e2 = (f2 / np.linalg.norm(f2)).reshape(d,1)
e3 = (f3 / np.linalg.norm(f3)).reshape(d,1)
temp1 = e1
temp2 = e2
temp3 = e3

Xs = np.zeros(n)
Ys = np.zeros(n)
Zs = np.zeros(n)

def update_points():
    for i in np.arange(n):
        global Xs, Ys, Zs
        p0 = P[i, 0:d]
        Xs[i] = np.dot(p0, e1)
        Ys[i] = np.dot(p0, e2)
        Zs[i] = np.dot(p0, e3)

update_points()

print("init: ready")
print(Xs[14],Ys[14],Zs[14])

identifier = "14"

# sympy
a1, b1, c1, d1, a2, b2, c2, d2, a3, b3, c3, d3, t1, s1, u1, t2, s2, u2 = sp.symbols(
	'a1 b1 c1 d1 a2 b2 c2 d2 a3 b3 c3 d3 t1 s1 u1 t2 s2 u2')  # variables
x2_s, y2_s, z2_s = sp.symbols('x2_s y2_s z2_s')  # values
P_i = sp.MatrixSymbol('P_i', d, 1)
E1 = sp.MatrixSymbol('E1', d, 1)
E2 = sp.MatrixSymbol('E2', d, 1)
E3 = sp.MatrixSymbol('E3', d, 1)
_var = (x2_s, y2_s, z2_s, P_i, E1, E2, E3, a1, b1, c1, d1, a2, b2, c2, d2, a3, b3, c3, d3, t1, s1, u1, t2, s2, u2)

E0 = sp.Matrix(P_i - x2_s * E1 - y2_s * E2 - z2_s * E3)
E0 = sp.simplify(E0 / sp.Matrix.norm(E0))

_E1 = sp.Matrix(a1 * E1 + b1 * E2 + c1 * E3 + d1 * E0)
_E2 = sp.Matrix(a2 * E1 + b2 * E2 + c2 * E3 + d2 * E0)
_E3 = sp.Matrix(a3 * E1 + b3 * E2 + c3 * E3 + d3 * E0)
R1 = sp.Matrix(t1 * E1 + s1 * E2 + u1 * E3)
R2 = sp.Matrix(t2 * E1 + s2 * E2 + u2 * E3)

_f = Matrix([
	_E1.dot(_E1) - 1,
	_E2.dot(_E2) - 1,
	_E3.dot(_E3) - 1,
	_E1.dot(_E2),
	_E2.dot(_E3),
	_E3.dot(_E1),
	R1.dot(R1) - 1,
	R2.dot(R2) - 1,
	R1.dot(R2),
	_E1.dot(R1) - sp.Matrix(E1).dot(R1),
	_E2.dot(R1) - sp.Matrix(E2).dot(R1),
	_E3.dot(R1) - sp.Matrix(E3).dot(R1),
	_E1.dot(R2) - sp.Matrix(E1).dot(R2),
	_E2.dot(R2) - sp.Matrix(E2).dot(R2),
	_E3.dot(R2) - sp.Matrix(E3).dot(R2),
	sp.Matrix(P_i).dot(_E1) - x2_s,
	sp.Matrix(P_i).dot(_E2) - y2_s,
	sp.Matrix(P_i).dot(_E3) - z2_s
])

#sp.refine(_f,sp.Q.zero(Matrix(E1).dot(Matrix(E1))))

_func = sp.Matrix.norm(_f)

_lam_f = lambdify(_var, _func, 'numpy')


def _lam(x2, y2, z2, p, e_1, e_2, e_3):
	return lambda a1, b1, c1, d1, a2, b2, c2, d2, a3, b3, c3, d3, t1, s1, u1, t2, s2, u2: \
		_lam_f(x2, y2, z2, p, e_1, e_2, e_3, a1, b1, c1, d1, a2, b2, c2, d2, a3, b3, c3, d3, t1, s1, u1, t2, s2, u2)

print("end")

ons = np.ones(18)
_arr = np.array(ons)

_f2 = _lam(0, 0, 0, P[int(identifier)].reshape(d, 1), e1, e2, e3)

def _g(args): return _f2(*args)

res = opt.minimize(_g, _arr, method='L-BFGS-B',options={'ftol':1e-3})
print(res)
e1 = res.x[0] * temp1 + res.x[1] * temp2 + res.x[2] * temp3 + res.x[3] * P[int(identifier)].reshape(d, 1)
e2 = res.x[4] * temp1 + res.x[5] * temp2 + res.x[6] * temp3 + res.x[7] * P[int(identifier)].reshape(d, 1)
e3 = res.x[8] * temp1 + res.x[9] * temp2 + res.x[10] * temp3 + res.x[11] * P[int(identifier)].reshape(d, 1)
update_points()
print(Xs[14], Ys[14], Zs[14])