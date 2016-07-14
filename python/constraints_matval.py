import numpy as np
import sympy as sp
from scipy import optimize as opt
from sympy.utilities.lambdify import lambdify
from sympy import *
from tkinter import *
import time
from functools import *
from sympy.matrices.expressions.matmul import remove_ids

#  initialize
_wid = _hei = 700  # window's width and height
wid = hei = 500  # canvas's width and height

# load adjacency and multi-dimensional space
edge = np.genfromtxt('csv/adjacency.csv', delimiter=",").astype(np.int64)
edge_num = len(edge)
P = np.genfromtxt('csv/mdSpace.csv', delimiter=",")
node_num, high_dim = P.shape

dim = 2

def genE():
    L = np.sqrt(np.genfromtxt('csv/eigVals.csv', delimiter=",")[0:high_dim])
    base = np.zeros(high_dim * dim).reshape(dim, high_dim)
    for i in range(high_dim): base[i % dim][i] = 1
    E = base*L
    return E

Es = genE()
print(Es.shape)
temp1 = Es[0]
temp2 = Es[1]

Xs = np.zeros(node_num)
Ys = np.zeros(node_num)
Xs_scaled = np.zeros(node_num)
Ys_scaled = np.zeros(node_num)
boundingX = 0
boundingY = 0

def scale(pnt,bool):
	if(bool): return wid*(pnt + boundingX/2)/boundingX + (_wid - wid)/2
	else: return (hei-100)*(boundingY/2 - pnt)/boundingY + (_hei - hei)/2
def unscale(pnt,bool):
	if (bool): return boundingX * ((pnt - (_wid - wid)/2) - wid / 2) / wid
	else: return boundingY * ((pnt - (_hei - hei)/2) - (hei - 100) / 2) / (100 - hei)
def update_points():
	for i in range(node_num):
		global Xs, Ys, boundingX, boundingY
		p0 = P[i, 0:high_dim]
		Xs[i] = np.dot(p0, Es[0]) ; Ys[i] = np.dot(p0, Es[1])
	boundingX = max([np.amax(Xs), np.amin(Xs)]) * 2
	boundingY = max([np.amax(Ys), np.amin(Ys)]) * 2
	for i in range(node_num):
		Xs_scaled[i] = scale(Xs[i], True);Ys_scaled[i] = scale(Ys[i], False)

update_points()

print("init: ready")

# sympy
A = MatrixSymbol('A',high_dim,high_dim)
q = MatrixSymbol('q', 1, dim)  # values
E = sp.MatrixSymbol('E', high_dim, high_dim) # = (E[0], E[1], ..., E[dim-1] P[thisID] 0 ...).T
var = (q,E,A)

_E = E * A  # = (e1' e2' e3' ... R1 R2 ...)

constraints1 = remove_ids(refine(_E.T * _E, Q.orthogonal(E)))
constraints2 = remove_ids(refine(E * _E, Q.orthogonal(E)))

bases_e = Matrix(constraints1[0:dim,0:dim] - Identity(dim))
_bases_e = Matrix.norm(bases_e)

bases_r = Matrix(constraints1[dim:2*dim-1,dim:2*dim-1] - Identity(dim-1))
_bases_r = Matrix.norm(bases_r)

# _Ei.dot(Rj) - sp.Matrix(Ei).dot(Rj),
eMulR =  Matrix(constraints1[0:dim,dim:2*dim-1] - constraints2[0:dim,dim:2*dim-1])
_eMulR = Matrix.norm(eMulR)


# sp.Matrix(P_i).dot(_Ej) - wj
pew = Matrix(constraints2[dim,0:dim] - q)
_pew = Matrix.norm(pew)

f = Matrix([_bases_e,_bases_r,_eMulR,_pew])

func = sum(f)
lam_f = lambdify(var, func, 'numpy')

def lam(q_,E, P_i):
	E_ = np.zeros(high_dim*high_dim).reshape(high_dim,high_dim)
	partOfE_ = np.r_[E,P_i]
	E_[0:dim+2,0:high_dim] = partOfE_
	return lambda A_: lam_f(q_,E_,A_)

arr_init = np.zeros(high_dim * high_dim).reshape(high_dim, high_dim)

def arr_initializer(a,b):
    print(type(a), a.shape, a)
    arr_init[0:dim+1,0:dim] = a		# E' variables ( E'[0] = this[0:dim+1,0] * (E:E0) )
    arr_init[0:dim,dim:2*dim] = b	# R  variables ( ignore )
    return f2(arr_init)

"""
init_E = np.r_[np.diag([1 for i in range(dim)]),np.zeros(dim).reshape(1,dim)].reshape(dim*(dim+1),)
init_R = np.ones(dim*(dim-1)).reshape(dim,dim-1)
init = np.array([init_E, init_R])
"""
init_E = np.r_[np.diag([1 for i in range(dim)]),np.zeros(dim).reshape(1,dim)].reshape(dim*(dim+1),)
init_R = np.ones(dim*(dim-1))
init = np.r_[init_E, init_R]

f2 = lam(np.array([[0,0]]), np.r_[Es,P[14].reshape(1,high_dim)], P[14].reshape(1,high_dim))

def g(args):
	arg1 = args[0:2*dim+dim].reshape(dim+1,dim)
	arg2 = args[2*dim+dim:4*dim].reshape(dim,dim-1)
	return arr_initializer(arg1,arg2)

res = opt.minimize(g, init, method='L-BFGS-B',options={'ftol':1e-3})
print(res)
p_i = P[14] / np.linalg.norm(P[14])
if(res.success):
    Es[0] = res.x[0] * temp1 + res.x[2] * temp2 + res.x[4] * p_i
    Es[1] = res.x[1] * temp1 + res.x[3] * temp2 + res.x[5] * p_i
    update_points()
    print(Xs[14])
    print(Ys[14])