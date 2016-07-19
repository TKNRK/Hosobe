import numpy as np
import sympy as sp
from scipy import optimize as opt
from sympy.utilities.lambdify import lambdify
from sympy import *
from tkinter import *
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
temp1 = Es[0]
temp2 = Es[1]

Xs = np.zeros(node_num)
Ys = np.zeros(node_num)

def update_points():
	for i in range(node_num):
		global Xs, Ys
		p0 = P[i, 0:high_dim]
		Xs[i] = np.dot(p0, Es[0]) ; Ys[i] = np.dot(p0, Es[1])

update_points()

print("init: ready")

# sympy
A = MatrixSymbol('A',12,12)	# coefficient matrix : Aのi列目が ei' の係数
xs = MatrixSymbol('xs', 1, 12)  # updated point
q = MatrixSymbol('q', 1, 2)  # updated point
E = MatrixSymbol('E', 12,12) # = (e1 e2 e0 P[thisID] 0 ...) : 縦ベクトルの列
var = (q,xs,E,A)

_E = E * A  # = (e1' e2' R1 0 ...) : 縦ベクトルの列

constraints1 = remove_ids(refine(_E.T * _E, Q.orthogonal(E)))
constraints2 = remove_ids(refine(E.T * _E, Q.orthogonal(E)))
constraints3 = E.T * _E

# ei' * ej' = δij (クロネッカーのデルタ)
bases_e = Matrix(constraints1[0:2,0:2] - Identity(2))
_bases_e = bases_e[0,0]**2 + bases_e[0,1]**2 + bases_e[1,1]**2

# Ri * Rj = δij (クロネッカーのデルタ)
bases_r = Matrix(constraints1[2:3,2:3] - Identity(1))
_bases_r = bases_r[0,0]**2

# _Ei.dot(Rj) - sp.Matrix(Ei).dot(Rj),
eMulR =  Matrix(constraints1[0:2,2:3] - constraints2[0:2,2:3])
_eMulR = eMulR[0,0]**2 + eMulR[1,0]**2

# sp.Matrix(P_i).dot(_Ej) - wj
pew = Matrix(constraints3[3,0:2] - q)
#_pew = pew[0,1]**2 + pew[0,0]**2
_pew = Matrix((xs * A[:,0] - Identity(1)) + (xs * A[:,1] - Identity(1)))[0,0]

f = Matrix([_bases_e,_bases_r,_eMulR,_pew])

func = sum(f)
lam_f = lambdify(var, func, 'numpy')
#print(lambdastr(lam_f))

def lam(q_,Esub, P_i):
    """
    q_    :(1,low_dim) ドラッグされた行先の点の座標
    Esub  :(low_dim,high_dim) 画面が更新される前の基底
    P_i   :(1,high_dim) ドラッグされた点に対応する高次元座標
    """
    # a[:,0:3]
    E0 = P_i - q_[0,0]*Esub[0].reshape(1,12) - q_[0,1]*Esub[0].reshape(1,12)
    xs_ = np.array([[np.linalg.norm(E0),Xs[14],Ys[14],0,0,0,0,0,0,0,0,0]])  # ?
    E0 = E0 / np.linalg.norm(E0)
    E_ = np.zeros(12*12).reshape(12,12)
    partOfE_ = np.r_[np.r_[Esub,E0],P_i]
    E_[0:4,0:12] = partOfE_
    return lambda A_: lam_f(q_, xs_, E_.T, A_)

arr_init = np.zeros(high_dim * high_dim).reshape(high_dim, high_dim)

def arr_initializer(a,b):
    arr_init[0:3,0:2] = a		# E' variables ( E'[0] = this[0:dim+1,0] * (E:E0) )
    arr_init[0:2,2:3] = b	# R  variables ( ignore )
    return f2(arr_init)

init = np.array([1,0,0,0,1,0,1,1])

f2 = lam(np.array([[0,0]]), Es, P[14].reshape(1,12))

def g(args):
	arg1 = args[0:6].reshape(3,2)
	arg2 = args[6:8].reshape(2,1)
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

point = np.array([[5,10]]).reshape(1,2)
p00 = np.random.randn(12).reshape(1,12)
e00 = np.random.randn(12).reshape(1,12)
e000 = np.random.randn(12).reshape(1,12)
e0 = p00 - 5*e00 - 10*e000
e0 = e0 / np.linalg.norm(e0)
sampleA = np.zeros(144).reshape(12,12)
sampleA[0,0] = 1
sampleA[0,1] = 1
sampleA[0,2] = 1
sampleA[1,0] = 1
sampleA[1,1] = 1
sampleA[1,2] = 1
sampleA[2,0] = 1
sampleA[2,1] = 1
sampleE = np.zeros(144).reshape(12,12)
sampleE[:,0] = e00
sampleE[:,1] = e000
sampleE[:,2] = e0
print(lam_f(point,e0,sampleE,sampleA))
print(g(np.array([1,1,1,1,1,1,1,1])))