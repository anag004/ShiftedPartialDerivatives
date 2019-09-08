from scipy import *
from numpy.linalg import matrix_rank
from flint import *
from progress.bar import Bar

P = 1000000000 + 7

# Takes the derivative of polynomial f wrt x
def getderiv(f):
    return nmod_poly(poly1d(f.coeffs()[::-1]).deriv().c[::-1].tolist(), P)

# Takes the derivative of the polynomial wrt y
def getderivy(f):
    ansCoeffs = []
    coeffArr = f.coeffs()
    for i in range(f.degree() + 1):
        ansCoeffs.append((f.degree() - i) * coeffArr[i])
    return nmod_poly(ansCoeffs, P)

# Returns the relevant matrix
def get_mat(k, f):
    # print("Random polynomial: ") 
    # print(f)
    d = f.degree()
    count = 0
    # Create a matrix with (d+1) entries
    M = zeros((d+1, (k+1)**2))
    for i in range(k+1):
        ftmp = f 
        for l in range(i):
            ftmp = getderiv(ftmp)
        for l in range(k-i):
            ftmp = getderivy(ftmp)

        for j in range(k+1):
            # Multiply by x**j
            xj = nmod_poly([0, 1], P)
            xj = xj ** j
            res = xj * ftmp
            # print("For i: " + str(i) + " and j: " + str(j) + " the polynomial is: ")
            # print(res)
            coeffs = array(res.coeffs()[::-1])
            refSize = array([0 for r in range(d + 1)])
            if (coeffs.shape[0] == 0):
                coeffs = array([0])
            refSize[-coeffs.shape[0]:] = coeffs
            M[:, count] = refSize
            count += 1
    # print(M)
    return M

def D(k, f):
    X = get_mat(k, f).astype(int).tolist()
    M = nmod_mat(X, P)
    return M.rank()


# Generate a random polynomial of degree t
def randq(t):
    while True:
        arr = random.randint(low=0, high=100, size=(t+1)).tolist()
        ans = nmod_poly(arr, P)
        if ans.degree() == t:
            return ans


def randf(s, e, t):
    res = nmod_poly([0], P)
    for j in range(s):
        # Generate a randomQ
        qj = randq(t)
        # Generate a random aj
        aj = random.randint(low=1, high=5)
        res = res + aj * (qj ** e)
    return res

def conj1(s, e, t, k, b):
    # width = (k+1)*(l+1)+(k+1)*k/2
    # M = zeros(((l + e*t + 1), s*width))
    # ans = min(s*width, s*(l+k*t+1), l + e*t + 1)
    # arr = []
    # for i in range(s):
    #     q = randq(t) 
    #     arr.append(q)
    #     q = q ** e
    #     Mi = get_mat(k, l, q)
    #     M[:, i*width:(i+1)*width] = Mi 
    # res = nmod_mat(M.astype(int).tolist(), P).rank()
    # if b:
    #     print("EXPECTED: " + str(ans) + " GOT: " + str(res))   
    # if res != ans:
    #     print("FAILED: " + str(res))
    #     print("EXPECTED: " + str(ans) + " GOT: " + str(res))
    #     print(arr)
    f = randf(s, e, t)
    # print("Random polynomial: ")
    # print(f)
    ans = min((k+1)**2, s*(k*t+1), e*t + 1)
    res = D(k, f)
    if b and res != ans:
        print("EXPECTED: " + str(ans) + " GOT: " + str(res))   
    return res == ans

def conj2(s, e, t, k, b, n):
    # Generate n random Qjs and average their t values
    r = 0
    for i in range(n):
        qj = randq(t)
        r = max(r, D(k, qj ** e))
    

    if b:
        print("EXPECTED: " + str(min((s*r), (k+1)**2, e*t+1)))

    # Generate random fs and find the d-value
    sr = 0
    for i in range(n):
        f = randf(s, e, t)
        sr = max(sr, D(k, f))
    
    if b:
        print("GOT: " + str(sr))

# conj1(2, 8, 2, 3, True)

# p = 0
# f = 0
# for i in range(100):
#     if conj1(2, 33, 3, 16, True):
#         p += 1
#     else:
#         f += 1
# print("PASSED: " + str(p) + " FAILED: " + str(f))
    
# elim = 100
# klim = 50  
# slim = 10
# t = 5
# p = 0
# f = 0

# bar = Bar('Processing', max=slim-1)

# for s in range(1, slim+1):
#     for k in range(1, klim + 1):
#         if ((k+1)**2 > s*(k*t + 1)):
#             for e in range(1, elim + 1):
#                 if ((e*t + 1) > s*(k*t + 1)):
#                     if (conj1(s, e, t, k,  True)):
#                         p += 1
#                         # print("TRYING " + "e = " + str(e) + " s = " + str(s) + " k = " + str(k) + " l = " + str(l) + str(" PASSED"))
#                     else:
#                         print("TRYING " + "e = " + str(e) + " s = " + str(s) + " k = " + str(k) + str(" FAILED"))
#                         f += 1
#                     break
#     bar.next()
# bar.finish()
# print("PASSED: " + str(p))
# print("FAILED: " + str(f))

# for e in range(1, elim):
#     for k in range(e*t/2):
#         slim = min((k+1)**2, e*t + 1)/((k*t + 1)) - 1
#         # print("S range is " + str(slim))
#         for s in range(1, slim):
#             print("TRYING " + "e = " + str(e) + " s = " + str(s) + " k = " + str(k))
#             if (conj1(s, e, t, k,  True)):
#                 p += 1
#                 # print("TRYING " + "e = " + str(e) + " s = " + str(s) + " k = " + str(k) + " l = " + str(l) + str(" PASSED"))
#             else:
#                 print("FAILED")
#                 f += 1
#     bar.next()
# bar.finish()
# print("PASSED: " + str(p))
# print("FAILED: " + str(f))

bar = Bar('Processing', max=99)
# Test that for a simple polynomial dimension is kt + 1
n = 10
for t in range(3, 100):
    k = t - 2
    e = 2*k
    for i in range(n):
        q = (randq(t) ** e)
        ans = D(k, q)
        if (k*t + 1 != ans):
            print("Failed for t = " + str(t) + " expected " + str(k*t + 1) + " got " + str(ans))
    bar.next()
bar.finish()
