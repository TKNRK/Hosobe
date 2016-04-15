import numpy as np

n = 10
d = np.random.randn(n**2).reshape(n,n)
for i in np.arange(n):
	for j in np.arange(i):
		d[i][j] = d[j][i]

def culA():
	a = np.zeros(n**2).reshape(n,n)
	dr2 = dc2 = np.zeros(n)
	for i in np.arange(n):
		for j in np.arange(n):
			dr2[i] += d[j][i]**2
			dc2[i] += d[i][j]**2
	da2 = 0
	for i in np.arange(n):
		da2 += dr2[i]
	da2 = da2 / n**2
	dr2 = dr2 / 2
	dc2 = dc2 / 2

	for i in np.arange(n):
		for j in np.arange(n):
			a[i][j] = ((dc2[i]+dr2[j]) - da2 - d[i][j]*2)/2

	return a

a = culA()

L,X = np.linalg.eigh(a)
L.sort()
L = L[::-1]
Ln = np.sqrt(np.diag(L))

P = X * Ln

f2 = np.array(L)
f2[::2] = 0
f1 = L - f2
e1 = f1 / np.linalg.norm(f1)
e2 = f2 / np.linalg.norm(f2)
