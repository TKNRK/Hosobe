import numpy as np
import sympy as sp
from scipy import optimize as opt
from sympy.utilities.lambdify import lambdify
from sympy import *
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
E = genE()
Xs = np.zeros(node_num)
Ys = np.zeros(node_num)
Xs_scaled = np.zeros(node_num)
Ys_scaled = np.zeros(node_num)
boundingH = 0
boundingV = 0
def scale(pnt,bool):
	if(bool): return wid*(pnt + boundingH / 2) / boundingH + (_wid - wid) / 2
	else: return (hei-100)*(boundingV / 2 - pnt) / boundingV + (_hei - hei) / 2
def unscale(pnt,bool):
	if (bool): return boundingH * ((pnt - (_wid - wid) / 2) - wid / 2) / wid
	else: return boundingV * ((pnt - (_hei - hei) / 2) - (hei - 100) / 2) / (100 - hei)
def update_points():
	for i in range(node_num):
		global Xs, Ys, boundingH, boundingV
		p0 = P[i, 0:high_dim]
		Xs[i] = np.dot(p0, E[0]) ; Ys[i] = np.dot(p0, E[1])
	boundingX = max([np.amax(Xs), np.amin(Xs)]) * 2
	boundingY = max([np.amax(Ys), np.amin(Ys)]) * 2
	for i in range(node_num):
		Xs_scaled[i] = scale(Xs[i], True);Ys_scaled[i] = scale(Ys[i], False)
update_points()
print("init: ready")

# sympy
a1,b1,c1,a2,b2,c2,t,s = symbols('a1 b1 c1 a2 b2 c2 t s')   # variables
x2_s,y2_s = symbols('x2_s y2_s')  # values
P_i = MatrixSymbol('P_i', high_dim, 1)
_E = MatrixSymbol('_E', high_dim, high_dim)
var = (x2_s,y2_s,P_i,_E,a1,b1,c1,a2,b2,c2,t,s)

_E1 = a1*_E + b1*_E + c1*P_i
_E2 = a2*_E + b2*_E + c2*P_i
R = s*_E + t*_E

constraints1 = refine(E,Q.orthogonal(E))

f = Matrix([
		_E1.dot(_E1) - 1,
		_E2.dot(_E2) - 1,
		_E1.dot(_E2),
		R.dot(R) - 1,
		_E1.dot(R) - sp.Matrix(_E).dot(R),
		_E2.dot(R) - sp.Matrix(_E).dot(R),
		sp.Matrix(P_i).dot(_E1) - x2_s,
		sp.Matrix(P_i).dot(_E2) - y2_s
		])

func = sp.Matrix.norm(f)
lam_f = lambdify(var, func, 'numpy')

def lam(x2, y2, p, e_1, e_2):
    return lambda a1,b1,c1,a2,b2,c2,t,s: \
        lam_f(x2, y2, p, e_1, e_2, a1, b1, c1, a2, b2, c2, t, s)

arr_init = np.array([1, 0, 0, 0, 1, 0, 1, 1])
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
    global E
    temp1 = E[0]
    temp2 = E[1]
    x2 = unscale(event.x,True)
    y2 = unscale(event.y,False)
    thisID = event.widget.find_withtag(CURRENT)[0] - (edge_num+1)
    f2 = lam(x2, y2, P[thisID].reshape(high_dim, 1), E[0].reshape(high_dim, 1),E[1].reshape(high_dim, 1))
    def g(args): return f2(*args)
    res = opt.minimize(g, arr_init, method='L-BFGS-B',options={'ftol':1e-3})
    if(res.success):
        E[0] = res.x[0] * temp1 + res.x[1] * temp2 + res.x[2] * P[thisID].reshape(high_dim, 1)
        E[1] = res.x[3] * temp1 + res.x[4] * temp2 + res.x[5] * P[thisID].reshape(high_dim, 1)
        update_points()
        for i in range(node_num):
            w.coords(circles[i], Xs_scaled[i] - r, Ys_scaled[i] - r, Xs_scaled[i] + r, Ys_scaled[i] + r)
        for i in range(edge_num):
            w.coords(lines[i], Xs_scaled[edge[i][0]-1], Ys_scaled[edge[i][0]-1],Xs_scaled[edge[i][1]-1], Ys_scaled[edge[i][1]-1])

# バインディング
w.tag_bind('node', '<Button1-Motion>', move_node)

root.mainloop()
