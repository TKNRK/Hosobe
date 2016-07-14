"""if __name__ == "__main__":

    import Rot_3D
    import numpy as np

    a = np.array([1,1,1])
    b = np.array([2,2,2])
    P = np.arange(36).reshape(12,3)
    r = Rot_3D.rot_3D(a,b)
    #print(r.rotation(P))
   """
from sympy import *
from scipy import optimize as opt
import numpy as np

B = MatrixSymbol('B',2,2)
C = MatrixSymbol('C',2,2)
_B = Matrix(B)
_C = Matrix(C)

fa = _B[0,1]**2 + _B[0,0]**2 + _C[0,0]**2 + _C[0,1]**2

function = lambdify((B,C),fa)

initial1 = np.array([[1,1],[1,1]])
initial2 = np.array([[1,1],[1,1]])
initial = np.array([initial1, initial2])

def function2(args):
    arg1 = args[0:4].reshape(2,2)
    arg2 = args[4:8].reshape(2,2)
    return function(arg1,arg2)

res = opt.minimize(function2, initial, method='L-BFGS-B',options={'ftol':1e-3})
print(res)