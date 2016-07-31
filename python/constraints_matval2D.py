import numpy as np
from scipy import optimize as opt
from sympy.utilities.lambdify import lambdify, lambdastr
from sympy import *
from sympy.matrices.expressions.matmul import remove_ids

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
q = MatrixSymbol('q',1,2)  # updated position
p_i = MatrixSymbol('p_i', 1, 12)  # selected point (high_dim) (P[thisID])
E = MatrixSymbol('E', 12,12) # = (e1 e2 e0 0 ...) : 縦ベクトルの列
A_inputEs = MatrixSymbol('A_inputEs', 3, 2)  # ei' の係数行列
A_inputRs = MatrixSymbol('A_inputRs', 2, 1)  # Ri の係数行列
var = (q, p_i, E, A_inputEs, A_inputRs)
A = Matrix(np.zeros(144).reshape(12,12))
A[0:3,0:2] = Matrix(A_inputEs)
A[0:2,2:3] = Matrix(A_inputRs)

_E = E * A  # = (e1' e2' R1 0 ...) : 縦ベクトルの列

constraints1 = (refine(_E.T * _E, Q.orthogonal(E))).doit()
constraints2 = (refine(E.T * _E, Q.orthogonal(E))).doit()
constraints3 = p_i * _E

# ei' * ej' = δij (クロネッカーのデルタ)
bases_e = Matrix(constraints1[0:2,0:2] - Matrix(Identity(2)))
_bases_e = bases_e[0,0]**2 + bases_e[0,1]**2 + bases_e[1,1]**2

# Ri * Rj = δij (クロネッカーのデルタ)
bases_r = Matrix(constraints1[2:3,2:3] - Matrix(Identity(1)))
_bases_r = bases_r[0,0]**2

# _Ei.dot(Rj) - sp.Matrix(Ei).dot(Rj),
eMulR =  Matrix(constraints1[0:2,2:3] - constraints2[0:2,2:3])
_eMulR = eMulR[0,0]**2 + eMulR[1,0]**2

# sp.Matrix(P_i).dot(_Ej) - wj
pew = Matrix(constraints3[0,0:2] - q)
_pew = pew[0,0]**2 + pew[0,1]**2

f = Matrix([_bases_e,_bases_r,_eMulR,_pew])

func = _bases_e + _bases_r + _eMulR + _pew
lam_f2 = lambdify(var, func, 'numpy')
def lam_f(*args):
    qq, pp, ee, aae, aar = args
    return lam_f2(*args)


def lam(q_, P_i, Esub):
    """
    q_    :(1,low_dim) ドラッグされた行先の点の座標
    Esub  :(low_dim,high_dim) 画面が更新される前の基底
    P_i   :(1,high_dim) ドラッグされた点に対応する高次元座標
    """
    # a[:,0:3]
    E0_temp = P_i - (q_[0,0]*Esub[0,:].reshape(1,12) + q_[0,1]*Esub[1,:].reshape(1,12))
    E0 = E0_temp / np.linalg.norm(E0_temp)
    E_ = np.zeros(12*12).reshape(12,12)
    partOfE_ = np.r_[Esub,E0]
    E_[0:3,0:12] = partOfE_
    return lambda A_1, A_2: lam_f(q_, P_i, E_.T, A_1, A_2)

zahyou = np.array([[-6,3]])
f2 = lam(zahyou, P[14].reshape(1,12), Es)

def g(args):
    """
    最適化関数
    :param args: 行列を行ベクトルに崩したもの
    :return: f2 に args を叩き込んだもの
    """
    arr1 = args[0:6].reshape(3,2)  # E' variables ( E'[0] = this[0:dim+1,0] * (E:E0) )
    arr2 = args[6:8].reshape(2,1)  # R  variables ( ignore )
    return f2(arr1,arr2)

init = np.array([1,0,0,1,0,0,0,0])
res = opt.minimize(g, init, method='L-BFGS-B')
print(res)

print("original x: "+str(Xs[14]))
print("original y: "+str(Ys[14]))
temp1 = Es[0].copy() ; temp2 = Es[1].copy()
p_iii = P[14].reshape(1,12) - (zahyou[0,0]*Es[0,:].reshape(1,12) + zahyou[0,1]*Es[1,:].reshape(1,12))
p_ii = p_iii / np.linalg.norm(p_iii)

Es[0] = res.x[0] * temp1 + res.x[2] * temp2 + res.x[4] * p_ii
Es[1] = res.x[1] * temp1 + res.x[3] * temp2 + res.x[5] * p_ii

update_points()
print("updated x: "+str(Xs[14]))
print("updated y: "+str(Ys[14]))
