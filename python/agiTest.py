import numpy as np
import sympy as sp
from scipy import optimize as opt
from sympy.utilities.lambdify import lambdify
from sympy import Matrix
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

dim = 3

def genE():
    L = np.sqrt(np.genfromtxt('csv/eigVals.csv', delimiter=",")[0:high_dim])
    base = np.zeros(high_dim * dim).reshape(dim, high_dim)
    for i in range(high_dim): base[i % dim][i] = 1
    E = base*L
    return E

E = genE()

Xs = np.zeros(node_num)
Ys = np.zeros(node_num)
Zs = np.zeros(node_num)
Xs_scaled = np.zeros(node_num)
Ys_scaled = np.zeros(node_num)
Zs_scaled = np.zeros(node_num)
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
		global Xs, Ys, Zs, boundingX, boundingY
		p0 = P[i, 0:high_dim]
		Xs[i] = np.dot(p0, E[0]) ; Ys[i] = np.dot(p0, E[1]) ; Zs[i] = np.dot(p0, E[2])
	boundingX = max([np.amax(Xs), np.amin(Xs)]) * 2
	boundingY = max([np.amax(Ys), np.amin(Ys)]) * 2
	for i in range(node_num):
		Xs_scaled[i] = scale(Xs[i], True);Ys_scaled[i] = scale(Ys[i], False)

update_points()

print("init: ready")

# sympy
a1, b1, c1, d1, a2, b2, c2, d2, a3, b3, c3, d3, t1, s1, u1, t2, s2, u2 = sp.symbols(
	'a1 b1 c1 d1 a2 b2 c2 d2 a3 b3 c3 d3 t1 s1 u1 t2 s2 u2')  # variables
x2_s, y2_s, z2_s = sp.symbols('x2_s y2_s z2_s')  # values
P_i = sp.MatrixSymbol('P_i', high_dim, 1)
E1 = sp.MatrixSymbol('E1', high_dim, 1)
E2 = sp.MatrixSymbol('E2', high_dim, 1)
E3 = sp.MatrixSymbol('E3', high_dim, 1)
_var = (x2_s, y2_s, z2_s, P_i, E1, E2, E3, a1, b1, c1, d1, a2, b2, c2, d2, a3, b3, c3, d3, t1, s1, u1, t2, s2, u2)

E0 = sp.Matrix(P_i - x2_s * E1 - y2_s * E2 - z2_s * E3)
E0 = sp.simplify(E0 / sp.Matrix.norm(E0))

_E1 = sp.Matrix(a1 * E1 + b1 * E2 + c1 * E3 + d1 * E0)
_E2 = sp.Matrix(a2 * E1 + b2 * E2 + c2 * E3 + d2 * E0)
_E3 = sp.Matrix(a3 * E1 + b3 * E2 + c3 * E3 + d3 * E0)
R1 = sp.Matrix(t1 * E1 + s1 * E2 + u1 * E3)
R2 = sp.Matrix(t2 * E1 + s2 * E2 + u2 * E3)

print("Variable declaration finished.")

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

_func = sp.Matrix.norm(_f)

print("culcurated the norm")

_lam_f = lambdify(_var, _func, 'numpy')

print("lambdify ends.")

def _lam(x2, y2, z2, p, e_1, e_2, e_3):
	return lambda a1, b1, c1, d1, a2, b2, c2, d2, a3, b3, c3, d3, t1, s1, u1, t2, s2, u2: \
		_lam_f(x2, y2, z2, p, e_1, e_2, e_3, a1, b1, c1, d1, a2, b2, c2, d2, a3, b3, c3, d3, t1, s1, u1, t2, s2, u2)

ons = np.ones(18)
_arr = np.array(ons)

print("lambda : ready")

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
    global E
    x2 = unscale(event.x,True)
    y2 = unscale(event.y,False)
    thisID = event.widget.find_withtag(CURRENT)[0] - (edge_num+1)
    _f2 = _lam(x2, y2, Zs[thisID], P[int(thisID)].reshape(high_dim, 1), E[0].reshape(high_dim, 1),E[1].reshape(high_dim, 1),E[2].reshape(high_dim, 1))
    def _g(args): return _f2(*args)
    res = opt.minimize(_g, _arr, method='L-BFGS-B', options={'ftol': 1e-2})
    if(res.success):
        temp1 = E[0]
        temp2 = E[1]
        temp3 = E[2]
        e0 = P[thisID] - x2 * temp1 - y2 * temp2 - Zs[thisID] * temp3
        e0 = e0 / np.linalg.norm(e0)
        E[0] = res.x[0] * temp1 + res.x[1] * temp2 + res.x[2] * temp3 + res.x[3] * e0
        E[1] = res.x[4] * temp1 + res.x[5] * temp2 + res.x[6] * temp3 + res.x[7] * e0
        E[2] = res.x[8] * temp1 + res.x[9] * temp2 + res.x[10] * temp3 + res.x[11] * e0
        update_points()
        for i in range(node_num):
            w.coords(circles[i], Xs_scaled[i] - r, Ys_scaled[i] - r, Xs_scaled[i] + r, Ys_scaled[i] + r)
        for i in range(edge_num):
            w.coords(lines[i], Xs_scaled[edge[i][0]-1], Ys_scaled[edge[i][0]-1],Xs_scaled[edge[i][1]-1], Ys_scaled[edge[i][1]-1])

# バインディング
w.tag_bind('node', '<Button1-Motion>', move_node)

root.mainloop()
