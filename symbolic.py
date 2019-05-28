from sympy import *
from numpy import *
from sympy.abc import a, b, c, d, x

# A matrix in general for t = 3
t = 3
e = 4
k = 3
M = full((13, 10), S(0))
p = Poly(a*x**3 + b*x**2 + c*x + d, x)
p = p ** e
count = 0

for i in range(k+1):
	for j in range(i+1):
		q = Poly(x ** j, x)
		q *= p
		refSize = full((e*t + 1), S(0))
		coeffs = array(q.all_coeffs())
		print("COEFFS FOR i = " + str(i) + " j = " + str(j))
		print(coeffs)
		refSize[-coeffs.shape[0]:] = coeffs
		M[:, count] = refSize
		count += 1
	p = p.diff()

M = Matrix(M)
print(M.rank())



