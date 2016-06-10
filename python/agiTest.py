import numpy as np
import sympy as sp
from scipy import optimize as opt
from sympy.utilities.lambdify import lambdify
from sympy import Matrix
from tkinter import *


#  initialize
wid = 500  # view's width
hei = 500  # view's height
edge = np.genfromtxt('csv/adjacency.csv', delimiter=",").astype(np.int64)
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
temp1 = e1
temp2 = e2

print(P[0].shape)
print(type(P[0]))

Xs = np.zeros(n)
Ys = np.zeros(n)
Xs_scaled = np.zeros(n)
Ys_scaled = np.zeros(n)
boundingX = 0
boundingY = 0

def scale(pnt,bool):
	if(bool):
		return wid*(pnt + boundingX/2)/boundingX
	else:
		return (hei-100)*(boundingY/2 - pnt)/boundingY

def unscale(pnt,bool):
	if (bool):
	    return boundingX * (pnt - wid / 2) / wid
	else:
	    return boundingY * (pnt - (hei - 100) / 2) / (100 - hei)

def update_points():
	for i in range(n):
		global Xs, Ys, boundingX, boundingY
		p0 = P[i, 0:d]
		Xs[i] = np.dot(p0, e1) ; Ys[i] = np.dot(p0, e2)
	boundingX = max([np.amax(Xs), np.amin(Xs)]) * 2
	boundingY = max([np.amax(Ys), np.amin(Ys)]) * 2
	for i in range(n):
		Xs_scaled[i] = scale(Xs[i], True);Ys_scaled[i] = scale(Ys[i], False)

update_points()

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

func = sp.Matrix.norm(f)

lam_f = lambdify(var, func, 'numpy')

def lam(x2, y2, p, e_1, e_2):
    return lambda a1,b1,c1,a2,b2,c2,t,s: \
        lam_f(x2, y2, p, e_1, e_2, a1, b1, c1, a2, b2, c2, t, s)

arr = np.array([1, 1, 1, 1, 1, 1, 1, 1])
print("lambda: ready")

######## Graph Drawing ########

root = Tk()

w = Canvas(root, width=wid, height=hei, bg='White')
w.pack()

circles = []
lines = []
r = 10

# 初期描画

for e in edge:
	lines.append(w.create_line(Xs_scaled[e[0]-1], Ys_scaled[e[0]-1],Xs_scaled[e[1]-1], Ys_scaled[e[1]-1], fill='Black', tags='edge'))

for i in range(n):
	circles.append(w.create_oval(Xs_scaled[i] - r, Ys_scaled[i] - r, Xs_scaled[i] + r, Ys_scaled[i] + r, fill="White", tags='node'))

# 移動
def move_node(event):
    global e1, e2, temp1, temp2
    x = event.x
    y = event.y
    x2 = unscale(int(x),True)
    y2 = unscale(y,False)
    thisID = event.widget.find_withtag(CURRENT)[0] - (eN+1)
    f2 = lam(x2, y2, P[thisID].reshape(d, 1), e1, e2)
    def g(args): return f2(*args)
    res = opt.minimize(g, arr, method='L-BFGS-B')
    print(res.success)
    if(res.success):
        e1 = res.x[0] * temp1 + res.x[1] * temp2 + res.x[2] * P[thisID].reshape(d, 1)
        e2 = res.x[3] * temp1 + res.x[4] * temp2 + res.x[5] * P[thisID].reshape(d, 1)
        temp1 = e1
        temp2 = e2
        update_points()
        for i in range(n):
            w.coords(circles[i], Xs_scaled[i] - r, Ys_scaled[i] - r, Xs_scaled[i] + r, Ys_scaled[i] + r)
        for i in range(len(edge)):
            w.coords(lines[i], Xs_scaled[edge[i][0]-1], Ys_scaled[edge[i][0]-1],Xs_scaled[edge[i][1]-1], Ys_scaled[edge[i][1]-1])

# バインディング
w.tag_bind('node', '<Button1-Motion>', move_node)

root.mainloop()
