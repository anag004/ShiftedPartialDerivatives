from sympy import *
from sympy.abc import a, b, c
from numpy import *

M = array([
	[1, 0, 0, 0, 0, 0, 0, 0, 0],
	[0, 1, 0, 0, 0, 0, 0, 0, 0],
	[0, 0, 1, 0, 0, 0, 0, 0, 0],
	[0, 0, 1, 6, 15, 20, 15, 6, 1],
	[0, 1, 7, 21, 35, 35, 7, 1, 0], 
	[1, 8, 28, 56, 80, 56, 28, 8, 1]
	])

print(linalg.matrix_rank(M))
print(Matrix(M).rank())