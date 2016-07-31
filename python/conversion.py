import numpy as np
from sympy import *
from sympy.matrices.expressions.matmul import remove_ids

point = np.array([5,10]).reshape(1,2)
p00 = np.random.randn(12).reshape(1,12)
e00 = np.random.randn(12).reshape(1,12)
e000 = np.random.randn(12).reshape(1,12)
e0_temp = p00 - (point[0,0]*e00 + point[0,1]*e000)
e0 = e0_temp / np.linalg.norm(e0_temp)
sampleA = np.zeros(144).reshape(12,12)
sampleA[0,0] = 1;sampleA[0,1] = 1;sampleA[0,2] = 1;sampleA[1,0] = 1
sampleA[1,1] = 1;sampleA[1,2] = 1;sampleA[2,0] = 1;sampleA[2,1] = 1
sampleE = np.zeros(144).reshape(12,12)
sampleE[:,0] = e00
sampleE[:,1] = e000
sampleE[:,2] = e0
EE = sampleE[:,0:2].T

# sympy
#A = MatrixSymbol('A',12,12)	# coefficient matrix : Aのi列目が ei' の係数
A_inputEs = MatrixSymbol('A_inputEs', 3, 2)
A_inputRs = MatrixSymbol('A_inputRs', 2, 1)
q = MatrixSymbol('q',1,2)
p_i = MatrixSymbol('p_i', 1, 12)  # updated point
E = MatrixSymbol('E', 12,12) # = (e1 e2 e0 P[thisID] 0 ...) : 縦ベクトルの列
var = (q,p_i,E,A_inputEs,A_inputRs)
A = Matrix(np.zeros(144).reshape(12,12))
A[0:3,0:2] = Matrix(A_inputEs)
A[0:2,0:1] = Matrix(A_inputRs)

_E = E * A  # = (e1' e2' R1 0 ...) : 縦ベクトルの列

constraints1 = remove_ids(refine(_E.T * _E, Q.orthogonal(E)))
constraints2 = remove_ids(refine(E.T * _E, Q.orthogonal(E)))
constraints3 = p_i * _E

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
pew = Matrix(constraints3[0,0:2] - q)
_pew = pew[0,1]**2 + pew[0,0]**2
#_pew = Matrix((xs * A[:,0] - Identity(1)) + (xs * A[:,1] - Identity(1)))[0,0]

f = Matrix([_bases_e,_bases_r,_eMulR,_pew])

func = sum(f)
lam_f2 = lambdify(var, func, 'numpy')
def lam_f(*args):
    return lam_f2(*args)
#print(lambdastr(lam_f))

def lam(q_, Esub, P_i):
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
    return lambda A_E, A_R: lam_f(q_, P_i, E_.T, A_E, A_R)

f2 = lam(point, EE, p00)

init = np.array([1,0,0,0,1,0,1,1])

def g(args):
	arg1 = args[0:6].reshape(3,2)
	arg2 = args[6:8].reshape(2,1)
	return f2(arg1,arg2)

print("lambda: ready")


print(lam_f(point, p00, sampleE, sampleA[0:3,0:2], sampleA[0:2,0:1]))
print(g(np.array([1,1,1,1,1,1,1,1])))