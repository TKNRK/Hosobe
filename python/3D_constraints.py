d = 12
e1 = np.array([1,0,0,1,0,0,1,0,0,1,0,0])
e2 = np.array([0,1,0,0,1,0,0,1,0,0,1,0])
e3 = np.array([0,0,1,0,0,1,0,0,1,0,0,1])
temp1 = e1; temp2 = e2; temp3 = e3



# sympy
a1,b1,c1,d1,a2,b2,c2,d2,a3,b3,c3,d3,t,s,u = sp.symbols('a1 b1 c1 d1 a2 b2 c2 d2 a3 b3 c3 d3 t s u')   # variables
x2_s,y2_s,z2_s = sp.symbols('x2_s y2_s z2_s')  # values
P_i = sp.MatrixSymbol('P_i', d, 1)
E1 = sp.MatrixSymbol('E1', d, 1)
E2 = sp.MatrixSymbol('E2', d, 1)
E3 = sp.MatrixSymbol('E3', d, 1)
var = (x2_s,y2_s,z2_s,P_i,E1,E2,a1,b1,c1,d1,a2,b2,c2,d2,a3,b3,c3,d3,t,s,u)

_E1 = sp.Matrix(a1*E1 + b1*E2 + c1*E3 + d1*P_i)
_E2 = sp.Matrix(a2*E1 + b2*E2 + c2*E3 + d2*P_i)
_E3 = sp.Matrix(a3*E1 + b3*E2 + c3*E3 + d3*P_i)
R = sp.Matrix(t*E1 + s*E2 + u*E3)

f = Matrix([
		_E1.dot(_E1) - 1,
		_E2.dot(_E2) - 1,
		_E3.dot(_E3) - 1,
		_E1.dot(_E2),
		_E2.dot(_E3),
		_E3.dot(_E1),
		R.dot(R) - 1,
		_E1.dot(R) - sp.Matrix(E1).dot(R),
		_E2.dot(R) - sp.Matrix(E2).dot(R),
		_E3.dot(R) - sp.Matrix(E3).dot(R),
		sp.Matrix(P_i).dot(_E1) - x2_s,
		sp.Matrix(P_i).dot(_E2) - y2_s,
		sp.Matrix(P_i).dot(_E3) - z2_s
		])

func = sp.Matrix.norm(f)

lam_f = lambdify(var, func, 'numpy')

def lam(x2, y2, z3, p, e_1, e_2, e_3):
    return lambda a1,b1,c1,a2,b2,c2,t,s: \
        lam_f(x2, y2, p, e_1, e_2, a1, b1, c1, a2, b2, c2, t, s)

arr = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

f2 = lam(event.xdata, event.ydata, event.zdata, P[int(identifier)].reshape(d,1), e1, e2, e3)
def g(args): return f2(*args)
res = opt.minimize(g, arr, method='L-BFGS-B')
print(res)
e1 = res.x[0] * temp1 + res.x[1] * temp2 + res.x[2] * temp3 + res.x[3] * P[int(identifier)].reshape(d,1)
e2 = res.x[4] * temp1 + res.x[5] * temp2 + res.x[6] * temp3 + res.x[7] * P[int(identifier)].reshape(d,1)
e3 = res.x[8] * temp1 + res.x[9] * temp2 + res.x[10] * temp3 + res.x[11] * P[int(identifier)].reshape(d,1)