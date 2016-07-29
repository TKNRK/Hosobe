import numpy as np
import sympy as sp
from scipy import optimize as opt
from sympy.utilities.lambdify import lambdify
from sympy import Matrix, MatrixSymbol, refine, Identity, Q
from tkinter import *
import time

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

func = _bases_e + _bases_r + _eMulR + _pew
lam_f2 = lambdify(var, func, 'numpy')
def lam_f(*args):
    qq, pp, ee, aae, aar = args
    return lam_f2(*args)


def lam(q_, P_i, Esub, E_0):
    """
    q_    :(1,low_dim) ドラッグされた行先の点の座標
    Esub  :(low_dim,high_dim) 画面が更新される前の基底
    P_i   :(1,high_dim) ドラッグされた点に対応する高次元座標
    E_0   :(1,high_dim) P_i - q_*Esub を正規化したもの
    """
    # a[:,0:3]
    E_ = np.zeros(high_dim*high_dim).reshape(high_dim,high_dim)
    partOfE_ = np.r_[Esub,E_0]
    E_[0:dim+1,0:high_dim] = partOfE_
    return lambda A_1, A_2: lam_f(q_, P_i, E_.T, A_1, A_2)

arr_init = np.array([1,0,0,1,0,0,0,0])
print("lambda: ready")

######## Graph Drawing ########
root = Tk()
w = Canvas(root, width=_wid, height=_hei, bg='White')
w.pack()
circles = []
lines = []
r = 10

# 初期描画
for e in edge:
	lines.append(w.create_line(Xs_scaled[e[0]-1], Ys_scaled[e[0]-1],Xs_scaled[e[1]-1], Ys_scaled[e[1]-1], fill='Black', tags='edge'))
for i in range(node_num):
	circles.append(w.create_oval(Xs_scaled[i] - r, Ys_scaled[i] - r, Xs_scaled[i] + r, Ys_scaled[i] + r, fill="White", tags='node'))

# 移動
def move_node(event):
    global Es
    temp1 = Es[0]
    temp2 = Es[1]
    x2 = unscale(event.x,True)
    y2 = unscale(event.y,False)
    thisID = event.widget.find_withtag(CURRENT)[0] - (edge_num+1)
    E_0 = P[thisID] - (x2*Es[0] + y2*Es[1])
    E_0 = E_0 / np.linalg.norm(E_0)
    position = np.array([[x2,y2]])
    f2 = lam(position, P[thisID].reshape(1, high_dim), Es, E_0.reshape(1, high_dim))
    def g(args):
        """
        最適化関数
        :param args: 行列を行ベクトルに崩したもの
        :return: f2 に args を叩き込んだもの
        """
        arr1 = args[0:6].reshape(3, 2)  # E' variables ( E'[0] = this[0:dim+1,0] * (E:E0) )
        arr2 = args[6:8].reshape(2, 1)  # R  variables ( ignore )
        return f2(arr1, arr2)
    res = opt.minimize(g, arr_init, method='L-BFGS-B',options={'ftol':1e-3})
    if(res.success):
        Es[0] = res.x[0] * temp1 + res.x[2] * temp2 + res.x[4] * E_0
        Es[1] = res.x[1] * temp1 + res.x[3] * temp2 + res.x[5] * E_0
        update_points()
        for i in range(node_num):
            w.coords(circles[i], Xs_scaled[i] - r, Ys_scaled[i] - r, Xs_scaled[i] + r, Ys_scaled[i] + r)
        for i in range(edge_num):
            w.coords(lines[i], Xs_scaled[edge[i][0]-1], Ys_scaled[edge[i][0]-1],Xs_scaled[edge[i][1]-1], Ys_scaled[edge[i][1]-1])

# バインディング
w.tag_bind('node', '<Button1-Motion>', move_node)

root.mainloop()
