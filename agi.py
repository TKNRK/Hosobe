import numpy as np
import sympy as sp
import scipy as cp
from matplotlib import pyplot as plt
import DraggableCircle.py

###### initialize ############

edge = np.genfromtxt('csv/adjacency.csv', delimiter=",").astype(np.int64) 
P = np.genfromtxt('csv/mdSpace.csv', delimiter=",") 
n = len(P)
L =np.genfromtxt('csv/eigVals.csv', delimiter=",") 
L_pos = np.array([L[i] if L[i]>0 else 0 for i in range(n)])
d = np.count_nonzero(L_pos)	# d ... the number of positive values
Ln = np.sqrt(L_pos)

f2 = np.array(Ln[0:d])
f2[::2] = 0
f1 = Ln[0:d] - f2
e1 = f1 / np.linalg.norm(f1)
e2 = f2 / np.linalg.norm(f2)

Xs = np.array([])
Ys = np.array([])

for i in np.arange(n):
    p0 = P[i,0:d]
    Xs = np.append(Xs,np.dot(p0,e1))
    Ys = np.append(Ys,np.dot(p0,e2))

######## Graph Drawing ########

identifier = ""

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
        if(self.circle.get_label() == identifier):
            x0, y0, xpress, ypress = self.press
            dx = event.xdata - xpress
            dy = event.ydata - ypress
            #print('x0=%f, xpress=%f, event.xdata=%f, dx=%f, x0+dx=%f' %
            #      (x0, xpress, event.xdata, dx, x0+dx))
            self.circle.center = (x0+dx,y0+dy)
            self.circle.figure.canvas.draw()
            for i in range(len(edge)):
                if(edge[i,0]-1 == int(identifier)): 
                    edges[i].set_xdata((x0+dx,edges[i].get_xdata()[1]))
                    edges[i].set_ydata((y0+dy,edges[i].get_ydata()[1]))
                if(edge[i,1]-1 == int(identifier)): 
                    edges[i].set_xdata((edges[i].get_xdata()[0],x0+dx))
                    edges[i].set_ydata((edges[i].get_ydata()[0],y0+dy))
          
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
    circle = plt.Circle((Xs[i], Ys[i]), radius=0.2, fc='y',label=str(i))
    gca.add_patch(circle)
    nodes = np.append(nodes,circle)

for e in edge:
    line = plt.Line2D((Xs[e[0]-1], Xs[e[1]-1]), (Ys[e[0]-1],Ys[e[1]-1]), lw=1)
    gca.add_line(line)
    edges = np.append(edges,line)
    
dcs = []
for node in nodes:
    dc = DraggableCircle(node)
    dc.connect()
    dcs.append(dc)

plt.axis('scaled')
plt.show()



######## upDating #############

# newton
# *args = (x,y,x2,y2,n,p,e1,e2)
def agiN(*args):

	# assign
	x,y,x2,y2,p,e1,e2 = args

	# Define return values
	a1,b1,c1,a2,b2,c2,t,s = sp.symbols('a1 b1 c1 a2 b2 c2 t s')
	var = np.array([a1,b1,c1,a2,b2,c2,t,s])
	E1 = a1*e1 + b1*e2 + c1*p
	E2 = a2*e1 + b2*e2 + c2*p
	r  = s*e1 + t*e2

	# functions
	f = np.array([
		np.dot(E1,E1) - 1,
		np.dot(E2,E2) - 1,
		np.dot(E1,E2),
		np.dot(r,r) - 1,
		np.dot(E1,r) - np.dot(e1,r),
		np.dot(E2,r) - np.dot(e2,r),
		np.dot(p,E1) - x2,
		np.dot(p,E2) - y2
	])
	# Hessian
	df = np.array([])
	for i in range(len(var)):
		for j in range(len(f)):
			df = np.append(df, sp.diff(f[j],var[i]))
	df = df.reshape(len(var),len(f))

	# return the answer
	return 0

# BFGS
def agiB(*args):
	# assign
	x,y,x2,y2,p,e1,e2 = args

	# Define return values
	a1,b1,c1,a2,b2,c2,t,s = sp.symbols('a1 b1 c1 a2 b2 c2 t s')
	var = [a1,b1,c1,a2,b2,c2,t,s]
	E1 = a1*e1 + b1*e2 + c1*p
	E2 = a2*e1 + b2*e2 + c2*p
	r  = s*e1 + t*e2

	# functions
	f = np.array([
		np.dot(E1,E1) - 1,
		np.dot(E2,E2) - 1,
		np.dot(E1,E2),
		np.dot(r,r) - 1,
		np.dot(E1,r) - np.dot(e1,r),
		np.dot(E2,r) - np.dot(e2,r),
		np.dot(p,E1) - x2,
		np.dot(p,E2) - y2
	])

	func = np.dot(f,f)
	cp.minimize(func,var(0),method='L-BFGS-B',maxiter=5)