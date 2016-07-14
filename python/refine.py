from sympy import *
from functools import *
from sympy.matrices.expressions.matmul import remove_ids

#l = Symbol('l')
l = 5
T = MatrixSymbol('T', l, l)
M = MatrixSymbol('M', l, l) # M = [[e1],[e2],[e3],0]
tM = M * T

facts = Q.orthogonal(M)

t2M2 = remove_ids(refine(tM.T * tM, facts))

print(t2M2)
"""
a = Symbol('a',real=True)
b = Symbol('b',real=True)
e1 = a*M[:,0] + b*M[:,1]

print(simplify(refine(e1.T * e1, Q.orthogonal(M))))
"""

#print(reduce(lambda a, x: a + x, [ a2M2[i, j] for i in range(0, 3) for j in range(0, 3)]))