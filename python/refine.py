from sympy import *
from functools import *

l = Symbol('l')
T = MatrixSymbol('T', l, l)
M = MatrixSymbol('M', l, l)
tM = T * M

facts = Q.orthogonal(M) and Q.real_elements(T)

t2M2 = refine(tM * tM.T, facts)
print(t2M2)
"""
a = Symbol('a',real=True)
b = Symbol('b',real=True)
e1 = a*M[:,0] + b*M[:,1]

print(simplify(refine(e1.T * e1, Q.orthogonal(M))))
"""

#print(reduce(lambda a, x: a + x, [ a2M2[i, j] for i in range(0, 3) for j in range(0, 3)]))