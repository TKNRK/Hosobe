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
e1 = (f1 / np.linalg.norm(f1)).reshape(d,1)
e2 = (f2 / np.linalg.norm(f2)).reshape(d,1)
temp1 = e1
temp2 = e2

print(P[0].shape)
print(type(P[0]))

Xs = np.zeros(n)
Ys = np.zeros(n)

def update_points():
	for i in np.arange(n):
		global Xs, Ys
		p0 = P[i, 0:d]
		Xs[i] = np.dot(p0, e1)
		Ys[i] = np.dot(p0, e2)

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
			global e1, e2, temp1, temp2
			f2 = lam(event.xdata, event.ydata, P[int(identifier)].reshape(d,1), e1, e2)
			def g(args): return f2(*args)
			res = opt.minimize(g, arr, method='L-BFGS-B')
			print(res)
			e1 = res.x[0] * temp1 + res.x[1] * temp2 + res.x[2] * P[int(identifier)].reshape(d,1)
			e2 = res.x[3] * temp1 + res.x[4] * temp2 + res.x[5] * P[int(identifier)].reshape(d,1)
			print(type(e1))
			print(e1.shape)
			temp1 = e1
			temp2 = e2
			update_points()
			for i in np.arange(n):
				nodes[i].center = (Xs[i], Ys[i])
				nodes[i].figure.canvas.draw()
			for i in range(len(edge)):
				edges[i].set_xdata((Xs[edge[i,0] - 1], Xs[edge[i,1] - 1]))
				edges[i].set_ydata((Ys[edge[i,0] - 1], Ys[edge[i,1] - 1]))

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
