from numpy.linalg import matrix_rank
import numpy as np
from sympy import *
from scipy.special import binom
from flint import *
from progress.bar import Bar
import random
from sympy.abc import x, y
from math import ceil

P = 1000000000 + 7

def get_coeffs(p, deg):
    ans = []
    for ax in range(0, deg + 1):
        ans.append(p.coeff_monomial((x ** ax) * (y ** (deg - ax))))
    return ans

# # Returns the relevant matrix 
# def get_mat(k, f):
#     d = f.total_degree()
#     count = 0
#     # Create a matrix with (d+1)*(d+2)/2 entries
#     M = np.zeros(((k+1)**2, (d+1)))
#     xmon = poly(x, domain=FF(P))
#     ymon = poly(y, domain=FF(P))
    
#     for b in range(0, k+1):
#         for a in range(0, k+1):
#             ftmp = f
#             for ax in range(b):
#                 ftmp = ftmp.diff(x)
#             for bx in range(k-b):
#                 ftmp = ftmp.diff(y)
#             M[count] = get_coeffs((xmon ** a) * (ymon ** (k - a)) * ftmp, d)
#             count += 1
#     return M

# Returns the relevant matrix -- differentiate by k2, mult by k1 
def get_mat(k1, k2, f):
    d = f.total_degree()
    count = 0
    # Create a matrix with (d+1)*(d+2)/2 entries
    M = np.zeros(((k1+1)*(k2 + 1), (d-k2+k1+1)))
    # print("SHAPE EXPECTED: " + str(M.shape))
    xmon = poly(x, domain=FF(P))
    ymon = poly(y, domain=FF(P))

    for b in range(0, k2+1):
        for a in range(0, k1+1):
            ftmp = f
            for ax in range(b):
                ftmp = ftmp.diff(x)
            for bx in range(k2-b):
                ftmp = ftmp.diff(y)
            M[count] = get_coeffs((xmon ** a) * (ymon ** (k1 - a)) * ftmp, d-k2+k1)
            count += 1
    return M

def D(k1, k2, f):
    X = get_mat(k1, k2, f).astype(int).tolist()
    M = nmod_mat(X, P)
    # print(M)
    return M.rank()


# Generate a random polynomial of degree t
def randq(t):
    arr = None
    while True:
        arr = np.random.randint(low=0, high=P, size=(t+1))
        if arr[0] != 0:
            break 
    ans = poly(arr[0] * x**t, domain=FF(P))
    for i in range(1, t+1):
        if (arr[i] % P) != 0:
            ans += poly((arr[i] % P) * (x ** (t - i)) * (y ** i), domain=FF(P))
    return ans

def genq(t):
    if t == 0: 
        return 1
    else:
        ans = poly(x**t, domain=FF(P))
        for i in range(1, t+1):
            ans += poly(((1) % P) * (x ** (t - i)) * (y ** i), domain=FF(P))
        return ans
        # return poly(x ** t + x**(t-1)*y , domain = FF(P))

def randf(s, e, t):
    res = (randq(t) ** e) 
    for j in range(s - 1):
        # Generate a randomQ
        qj = randq(t)
        # Generate a random aj
        res += (qj ** e)
    return res

def randr(roots, lead):
    ans = poly(lead*x - lead*roots[0]*y, domain=FF(P))
    for i in range(1, len(roots)):
        ftmp = poly(x - roots[i]*y, domain=FF(P))
        ans *= ftmp 
    return ans

def bcirc(arr):
    n = len(arr) - 1
    for i in range(n+1):
        arr[i] = float(arr[i])/binom(n, i)
    M = np.zeros((n+1, n+1))
    for i in range(n+1):
        M[i] = arr 
        tmp = arr[0]
        arr.append(arr.pop(0))
    return np.linalg.det(M)

# klim = 10
# for k in range(1, klim):
#     ans = (k + 1)*(k+2)*(2*k + 3)/6
#     # q = randq(t)
#     mypow = (ans + 10)/2
#     mydeg = ans + 10
#     # print(q)
#     print("Power: " + str(mypow))
#     f = randq(mydeg)
#     print(f)
#     print("Evaluating for k = " + str(k))
#     print("GOT: " + str(D(k, f)) + " EXPECTED: " + str(ans))

# e = 10
# t = 3
# s = 5
# k = 3
# n = 100
# bar = Bar('Processing', max=n)
# for trial in range(n):
#     # Generate s random polynomials, get values
#     M = np.zeros((s, e*t + 1))
#     for count in range(s):
#         q = randq(t) ** e 
#         i = np.random.randint(low=0, high=k)
#         j = np.random.randint(low=0, high=k)
#         for ax in range(i):
#             q = q.diff(x)
#         for ay in range(k-i):
#             q = q.diff(y)
#         xmon = poly(x, domain=FF(P))
#         ymon = poly(y, domain=FF(P))
#         M[count] = get_coeffs((xmon ** j) * (ymon ** (k - j)) * q, e*t)

#     M = nmod_mat(M.astype(int).tolist(), P)
#     ans = M.rank()
#     if (ans != s):
#         print("Fails => expected " + str(s) + " got " + str(ans))
#     bar.next()
# bar.finish()

# # Test that dimension is full for random polynomials
# d = 100
# f = randq(d)
# for k in range((int)(sqrt(d+1)-1)):
#     ans = D(k, f)
#     if (ans != (k+1)**2):
#         print("Wrong expected " + str((k+1)**2) + " got " + str(ans))
#     else:
#         print("Success for k = "  + str(k))

# Test that dimension is full for random polynomials q^e
# t = 10
# e = 14
# print("Generating random polynomial...")
# f = genq(t) ** e
# print("Beginning computation....")
# k = t-2
# ans = D(k, f)
# if (ans != (k+1)**2):
#     print("Wrong expected " + str((k+1)**2) + " got " + str(ans))
# else:
#     print("Success for k = "  + str(k))


# Test that Rij is independent for fixed j
# xmon = poly(x, domain=FF(P))
# ymon = poly(y, domain=FF(P))
# s = 10
# t = 3
# k = 40
# i = k
# j = k
# M = np.zeros((s, k*t + 1))
# bar = Bar('Processing', max=s)
# for l in range(s):
#     q = randq(t) ** k
#     for ax in range(i):
#         q = q.diff(x)
#     for bx in range(j-i):
#         q = q.diff(y)
#     M[l] = get_coeffs((xmon ** j) * (ymon ** (k - j)) * q, k*t)
#     bar.next()
# bar.finish()

# # Test that some subset of polynomials is independent
# e = 5
# t = 5
# k = 4
# M = np.zeros((2*k+1, e*t + 1))
# q = randq(t) ** e
# count = 0
# xmon = poly(x, domain=FF(P))
# ymon = poly(y, domain=FF(P))
# for j in range(k):
#     qtmp = q
#     for l in range(k):
#         qtmp = qtmp.diff(x)
#     M[count] = get_coeffs((xmon ** j) * (ymon ** (k-j)) * qtmp, e*t)
#     count += 1
# for j in range(k+1):
#     qtmp = q
#     for l in range(j):
#         qtmp = qtmp.diff(x)
#     for l in range(k-j):
#         qtmp = qtmp.diff(y)
#     M[count] = get_coeffs((xmon ** k) * qtmp, e*t)
#     count += 1
# M = nmod_mat(M.astype(int).tolist(), P)
# print(M.rank())


# M = M.astype(int).tolist()
# print(nmod_mat(M, P).rank())

# Test the dimension conjecture -- t <= e + 2
# trials = 10
# e = 10
# t = 8
# print("Generating for e = " + str(e) + " and t = " + str(t))
# f = [(randq(t) ** e) for i in range(trials)]
# print("Starting testing...")
# for k in range(e):
#     pa = 0
#     fa = 0
#     expected = 0
#     if (k <= t-2):
#         expected = (k+1)**2
#     else:
#         expected = k*t + 1
#     bar=Bar('Processing', max=trials)
#     for fxy in f:
#         ans = D(k, fxy)
#         if (ans == expected):
#             pa+=1
#         else:
#             fa+=1
#         bar.next()
#     bar.finish()
#     print("For k = " + str(k) + " passed : " + str(pa) + " failed : " + str(fa))



# # Testing the conjecture that f1, f2, ..., fn have distinct roots => shifts of these are LI
# n = 10
# k = 5
# d = (n-1)*(k+1)   
# xmon = poly(x, domain=FF(P))
# ymon = poly(y, domain=FF(P))
# M = np.zeros((n*(k+1), d + k + 1))
# count = 0 
# root = 1
# for i in range(n):
#     # rootarr = []
#     # for l in range(d):
#     #     rootarr.append(root)
#     #     root += 1
#     f = randr(np.random.randint(low=0, high=P, size=d), random.randint(1, P))
#     print(f)
#     for j in range(k+1):
#         print("Generating count = " + str(count))
#         M[count] = get_coeffs((xmon ** j) * (ymon ** (k - j)) * f , d + k)
#         count += 1
# M = nmod_mat(M.astype(int).tolist(), P)
# print("Expected: " + str(n*(k+1)) + " Got: " + str(M.rank()))



# while (True):
#     roots = [int(t) for t in raw_input().split(" ")]
#     f = randr(roots, 1)
#     print(f)
#     print(get_coeffs(f, len(roots)))

# Test the conjecture that k2 = 1, k1 = t - 2, 2*(t-1) < e*t + 1
# e = 10
# tlim = 3
# k2 = 10
# for t in range(3, tlim+1):
#     k1 = t - 2
#     f = randq(t) ** e
#     print("Expected " + str((k1 + 1)*(k2 + 1)) + " Got " + str(D(k1, k2, f)))   

# Code to test PtQkt + P2t\del=tQkt
def getDim(pt, p2t, k, t):
    xmon = poly(x, domain=FF(P))
    ymon = poly(y, domain=FF(P))
    count = 0
    M = np.zeros(((t+1)*(k*t+1), (k+1)*t + 1))
    for j in range(k*t + 1):
        q = (xmon ** j) * (ymon ** (k*t - j))
        q1 = q
        for i in range(t+1):
            for ax in range(i):
                q1 = q.diff(x)
            for bx in range(t - i):
                q1 = q.diff(y)
            f = pt*q + p2t*q1
            # print(f)

            coeffs = get_coeffs(f, (k+1)*t)
            # print(coeffs)
            M[count] = coeffs
            count += 1
    M = nmod_mat(M.astype(int).tolist(), P)
    # print(M)
    return M.rank()

def krank(f):
    k = 0
    while(True):
        tmp = D(k, k, f)
        if (tmp != (k+1)**2):
            break
        k += 1
    return k-1


# t = 1
# klim = 100
# p = 0
# f = 0
# p1 = randr([5], 1)
# p2 = randr([10, 10], 1)
# for k in range(1, klim+1):
#     print("Testing for k = " + str(k) + " t = " + str(t))
#     ans  = getDim(p1, p2, k, t)
#     print("Expected = " + str((k+t)*t + 1) + " got " + str(ans))
#     if (ans == (k+t)*t + 1):
#         p += 1
#     else:
#         f += 1
# print("Passed " + str(p) + " failed " + str(f))

# Find a range of values
# t = 2
# e = 7
# k1 = 2
# k2 = 4
# s = 2
# M = 100

# assert k1 >= t - 2
# assert s*(k2*t - k2 + k1 + 1) <= (k1 + 1)*(k2 + 1)
# assert s*(k2*t - k2 + k1 + 1) <= e*t + 1
# assert k2 <= e

# # a = random.randint(1, M) 
# a = 1
# tmp = poly(x ** t + (a *  y) ** t, domain=FF(P))
# print(a)
# q = tmp ** e

# for j in range(2, s + 1):
#     # a = random.randint(1, M)
#     a = -1
#     tmp = poly(x ** t + (a * y) ** t, domain=FF(P))
#     print(a)
#     q += tmp ** e

# ans = D(k1, k2, q)
# print("ANSWER : " + str(D(k1, k2, q)))
# assert ans == s*(k2*t - k2 + k1 + 1)

# # Test that for random degree t polynomials this holds
# t = 3
# e = 150
# k1 = 30
# k2 = 17
# s = 3
# M = 100

# assert k1 >= t - 2
# assert 2*s*(k2*t - k2 + k1 + 1) <= (k1 + 1)*(k2 + 1)
# assert 2*s*(k2*t - k2 + k1 + 1) <= e*t + 1 + k1 - k2 
# assert k2 < e

# print("Testing for e = " + str(e) + " k1 = " + str(k1) + " k2 = " + str(k2) + " s = " + str(s))

# # q = randq(t) ** e
# a = random.randint(1, M) 
# print(a)
# tmp = poly(x ** t + a*(y ** t), domain=FF(P))
# q = tmp ** e

# for j in range(2, s + 1):
#     a = random.randint(1, M)
#     print(a) 
#     tmp = poly(x ** t +  a*(y ** t), domain=FF(P))
#     q += tmp ** e

# ans = D(k1, k2, q)
# print("ANSWER : " + str(D(k1, k2, q)))
# print("EXPECTED : " + str(s*(k2*t - k2 + k1 + 1)))
# assert ans == s*(k2*t - k2 + k1 + 1)

# t = 3
# k1lim = 20
# k2lim = 20
# slim = 6
# M = P
# p = 0
# f = 0
# N = 3
# # assert k1 >= t - 2
# # assert s*(k2*t - k2 + k1 + 1) <= (k1 + 1)*(k2 + 1)
# # assert s*(k2*t - k2 + k1 + 1) <= e*t + 1

# for k1 in range(t-2, k1lim + 1):
#     for k2 in range(0, k2lim+1):
#         for s in range(3, min(slim, int((k1 + 1)*(k2 + 1)/(2*(k2*t - k2 + k1 + 1)))) + 1):
#             e = 2*max(int((s*(k2*t - k2 + k1 + 1) - 1 + k2 - k1)/t) + 1, k2)
#             print("Testing for e = " + str(e) + " k1 = " + str(k1) + " k2 = " + str(k2) + " s = " + str(s))
#             ans=0
#             for i in range(N):
#                 p1 = 0
#                 a = random.randint(1, M) 
#                 b = random.randint(1, M)
#                 tmp = poly(b * x ** 2 + y ** 2  + b * a * x * y, domain=FF(P))
#                 q = (tmp ** e)

#                 # q = randq(t) ** e
            

#                 for j in range(2, s + 1):
#                     a = random.randint(1, M)
#                     b = random.randint(1, M)
#                     tmp = poly(b * x ** 2 + y ** 2  + b * a * x * y, domain=FF(P))
#                     q += (tmp ** e)

#                     # q += randq(t) ** e
                    
#                 ans = D(k1, k2, q)
#                 if (ans == s*(k2*t - k2 + k1 + 1)):
#                     p1 += 1
#             if p1 == 0:
#                 f += 1
#                 print("EXPECTED: " + str(s*(k2*t - k2 + k1 + 1)) + " GOT: " + str(ans))
#             else:
#                 p += 1

# print("PASSED: " + str(p))
# print("FAILED: " + str(f))

# print(D(2, 4, poly(2*x**14 + 42*x**10*y**4 + 70*x**6*y**8 + 14*x**2*y**12, domain=FF(P))))

# Function to test conjecture with a gap of n
# def testDimAdd(s, t, n):
#     M = P
#     k1 = int(ceil(2*n*s*(t-1)-1))
#     k2 = int(ceil(2*n*s-1))
#     e = int(max(ceil((float)(n*s*(k2*t - k2 + k1 + 1)-k1+k2-1)/t), k2+1))
#     # print("Choosing k1 = " + str(k1) + " k2 = " + str(k2) + " e = " + str(e))
#     assert n*s*(k2*t - k2 + k1 + 1) <= (k1 + 1)*(k2 + 1)
#     assert n*s*(k2*t - k2 + k1 + 1) <= e*t + 1 + k1 - k2
#     assert k2 < e 
#     a = random.randint(1, M) 
#     tmp = poly(x ** t + a*(y ** t), domain=FF(P))
#     q = tmp ** e
#     for j in range(2, s+1):
#         a = random.randint(1, M) 
#         tmp = poly(x ** t + a*(y ** t), domain=FF(P))
#         q += tmp ** e
#     return [D(k1, k2, q), s*(k2*t - k2 + k1 + 1)]

# slim = 10
# beg = 1.0
# incr = 0.1
# t = 2

# for s in range(slim + 1):
#     curr = beg
#     while True:
#         tmp = testDimAdd(s, t, curr)
#         if (tmp[0] == tmp[1]):
#             print("For s = " + str(s) + " ratio = " + str(curr))
#             break
#         else:
#             curr += incr

# Get the null space of a particular polynomial
def getNullSpace(f, k1, k2):
    # The matrix of coefficients
    X = get_mat(k1, k2, f)
    # The dimension of the operator space
    N = X.shape[0]
    # The coefficient type
    K = X.shape[1]
    print("SHAPE: N = " + str(N) + " K = " + str(K))
    # The symbols we will use
    eq_vars = []
    for b in range(0, k2+1):
        for a in range(0, k1+1):
            eq_vars.append(symbols('x_' + str(b) + '_' + str(a)))
    assert N == len(eq_vars)
    # The list of equations to solve
    eq_list = []
    for i in range(K):
        eq = S(0)
        # Iterate over the operators
        for j in range(N):
            eq += eq_vars[j] * X[j][i]
        eq_list.append(Eq(eq, 0))
    return solve(eq_list, eq_vars)

# Test the conjecture D_k+t-1, 2 = D_k, 1 (comp) x^t-1 dely , y^t-1 delx
def testNullSpace2Conjecture(k1, k2, t, p, A, B):
    xmon = poly(x, domain=FF(P))
    ymon = poly(y, domain=FF(P))
    p1 = (xmon ** (t-1)) * (p.diff(y)) - poly(A * (y ** (t - 1)), domain=FF(P)) * (p.diff(x))
    p2 = (xmon ** (t-1)) *(p.diff(y)) - poly(B * (y ** (t - 1)), domain=FF(P)) * (p.diff(x))
    # print(p1)
    # print(p2)
    # p1 = x ** (t - 1) * p.diff(y)
    # p2 = y ** (t - 1) * p.diff(x)
    d = p1.total_degree() 
    M = np.zeros((2*(k1 - t + 2)*k2, d + k1 - t - k2 + 3))
    print("SHAPE: " + str(M.shape))
    M[:(k1 - t + 2)*k2][:] = get_mat(k1 - t  + 1, k2 - 1, p1)
    M[(k1 - t + 2)*k2:][:] = get_mat(k1 - t + 1, k2 - 1, p2)
    X = nmod_mat(M.astype(int).tolist(), P)
    # print("SINGLE: " + str(D(k1 - t + 1, k2 - 1, p1)))
    return X.rank()

e = 10
k1 = 10
k2 = 3
t = 3
A = 1
B = 2
p = randq(e * t)
assert((k1 + 1)*(k2 + 1) <= 2*(k1 - t + 2)*k2)
q1 = poly(x ** t + A * (y ** t), domain=FF(P)) ** e
q2 = poly(x ** t + B * (y ** t), domain=FF(P)) ** e
# print("dim(D(q1^e + q2^e)) = " + str(D(k1, k2, q1 + q2)))
# print("dim(D(q1^e)) + dim(D(q2^e)) = " + str(D(k1, k2, q1) + D(k1, k2, q2)))
print("GOT: "+ str(testNullSpace2Conjecture(k1, k2, t, p, A, B)))
print("EXPECTED: " + str(D(k1, k2, p)))

# assert(D(k1, k2, q1 + q2) == D(k1, k2, q1) + D(k1, k2, q2))
