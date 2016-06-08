import numpy as np
import sympy as sp
from scipy import optimize as opt
from matplotlib import pyplot as plt
from sympy.utilities.lambdify import lambdify
from sympy import Matrix

#  initialize
edge = np.genfromtxt('csv/adjacency.csv', delimiter=",").astype(np.int64)
P = np.genfromtxt('csv/mdSpace.csv', delimiter=",")
n = len(P)
L = np.genfromtxt('csv/eigVals.csv', delimiter=",")
L_pos = np.array([L[i] if L[i] > 0 else 0 for i in range(n)])
d = np.count_nonzero(L_pos)
Ln = np.sqrt(L_pos)
f2 = np.array(Ln[0:d])
f2[::2] = 0
f1 = Ln[0:d] - f2
e1 = f1 / np.linalg.norm(f1)
e2 = f2 / np.linalg.norm(f2)
temp1 = e1
temp2 = e2

Xs = np.array([])
Ys = np.array([])
for i in np.arange(n):
    p0 = P[i,0:d]
    Xs = np.append(Xs,np.dot(p0,e1))
    Ys = np.append(Ys,np.dot(p0,e2))

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

lam_f = lambdify(var, f, 'numpy')

def lam(x2, y2, p, e_1, e_2):
    return lambda a1,b1,c1,a2,b2,c2,t,s: \
        np.linalg.norm(lam_f(x2, y2, sp.Matrix(p), sp.Matrix(e_1), sp.Matrix(e_2), a1, b1, c1, a2, b2, c2, t, s))

print("lambda: ready")

######## Graph Drawing ########
identifier = ""
arr = np.array([1, 1, 1, 1, 1, 1, 1, 1])

class DraggableCircle:
	def __init__(self, circle):
		self.circle = circle
		self.press = None

	def connect(self):
		'connect to all the events we need'
		self.cidpress = self.circle.figure.canvas.mpl_connect(
			'button_press_event', self.on_press)
		self.cidrelease = self.circle.figure.canvas.mpl_connect(
			'button_release_event', self.on_release)
		self.cidmotion = self.circle.figure.canvas.mpl_connect(
			'motion_notify_event', self.on_motion)

	def on_press(self, event):
		'on button press we will see if the mouse is over us and store some data'
		if event.inaxes != self.circle.axes: return

		contains, attrd = self.circle.contains(event)
		if not contains: return
		print('event contains', self.circle.center)
		x0, y0 = self.circle.center
		self.press = x0, y0, event.xdata, event.ydata
		global identifier
		identifier = self.circle.get_label()

	def on_motion(self, event):
		'on motion we will move the rect if the mouse is over us'
		if self.press is None: return
		if event.inaxes != self.circle.axes: return
		if (self.circle.get_label() == identifier):
			global e1, e2, Xs, Ys
			f2 = lam(event.xdata, event.ydata, P[int(identifier)], e1, e2)
			def g(args): return f2(*args)
			res = opt.minimize(g, arr, method='L-BFGS-B')
			print(res)
			e1 = res.x[0] * temp1 + res.x[1] * temp2 + res.x[2] * P[0]
			e2 = res.x[3] * temp1 + res.x[4] * temp2 + res.x[5] * P[0]
			# print('x0=%f, xpress=%f, event.xdata=%f, dx=%f, x0+dx=%f' %
			#      (x0, xpress, event.xdata, dx, x0+dx))
			"""
			self.circle.center = (x0 + dx, y0 + dy)
			self.circle.figure.canvas.draw()
			for i in range(len(edge)):
				if (edge[i, 0] - 1 == int(identifier)):
					edges[i].set_xdata((x0 + dx, edges[i].get_xdata()[1]))
					edges[i].set_ydata((y0 + dy, edges[i].get_ydata()[1]))
				if (edge[i, 1] - 1 == int(identifier)):
					edges[i].set_xdata((edges[i].get_xdata()[0], x0 + dx))
					edges[i].set_ydata((edges[i].get_ydata()[0], y0 + dy))
			"""

	def on_release(self, event):
		'on release we reset the press data'
		self.press = None
		self.circle.figure.canvas.draw()

	def disconnect(self):
		'disconnect all the stored connection ids'
		self.circle.figure.canvas.mpl_disconnect(self.cidpress)
		self.circle.figure.canvas.mpl_disconnect(self.cidrelease)
		self.circle.figure.canvas.mpl_disconnect(self.cidmotion)


gca = plt.gca()
nodes = np.array([])
edges = np.array([])

for i in np.arange(n):
	circle = plt.Circle((Xs[i], Ys[i]), radius=0.2, fc='y', label=str(i))
	gca.add_patch(circle)
	nodes = np.append(nodes, circle)

for e in edge:
	line = plt.Line2D((Xs[e[0] - 1], Xs[e[1] - 1]), (Ys[e[0] - 1], Ys[e[1] - 1]), lw=1)
	gca.add_line(line)
	edges = np.append(edges, line)

dcs = []
for node in nodes:
	dc = DraggableCircle(node)
	dc.connect()
	dcs.append(dc)

plt.axis('scaled')
plt.show()
