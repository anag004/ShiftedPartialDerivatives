from scipy import *
from numpy.linalg import matrix_rank
from flint import *
from progress.bar import Bar

P = 1000000000 + 7

# Takes the derivative of polynomial f
def getderiv(f):
	return nmod_poly(poly1d(f.coeffs()[::-1]).deriv().c[::-1].tolist(), P)

def get_mat(k, l, f):
    d = f.degree()
    count = 0
    M = zeros((l+d+1, (k+1)*(l+1)+(k+1)*k/2))
    for i in range(k+1):
        for j in range(i+l+1):
			
            # Get the polynomial x^j
			xj = nmod_poly([0, 1], P)
			xj = xj ** j
			res = xj * f
			coeffs = array(res.coeffs()[::-1])
			refSize = array([0 for r in range(l + d + 1)])
			refSize[-coeffs.shape[0]:] = coeffs
			M[:, count] = refSize
			count += 1
        # Differentiate f
        f = getderiv(f)
    return M

def D(k, l, f):
	X = get_mat(k, l, f).astype(int).tolist()
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

def conj1(s, e, t, k, l, b):
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
    ans = min((k+1)*(l+1)+(k+1)*k/2, s*(l+k*t+1), l + e*t + 1)
    res = D(k, l, f)
    if b and res != ans:
        print("EXPECTED: " + str(ans) + " GOT: " + str(res))   
    return res == ans

def conj2(s, e, t, k, l, b, n):
    # Generate n random Qjs and average their t values
    r = 0
    for i in range(n):
        qj = randq(t)
        r = max(r, D(k, l, qj ** e))
    

    if b:
        print("EXPECTED: " + str(min((s*r), k*(k+1)/2 + (k+1)*(l+1), l+e*t+1)))

    # Generate random fs and find the d-value
    sr = 0
    for i in range(n):
        f = randf(s, e, t)
        sr = max(sr, D(k, l, f))
    
    if b:
        print("GOT: " + str(sr))

elim = 250
llim = 200
klim = 200  
t = 2
p = 0
f = 0
bar = Bar('Processing', max=elim-1)
for e in range(200, elim):
    for l in range(llim):
        for k in range(e):
            slim = min((k*(k+1)+(l+1)*(k+1)), l + e*t + 1)/(4*(l + k*t + 1))
            for s in range(1, slim):

                if (conj1(s, e, t, k, l, True)):
                    p += 1
                    # print("TRYING " + "e = " + str(e) + " s = " + str(s) + " k = " + str(k) + " l = " + str(l) + str(" PASSED"))
                else:
                    print("TRYING " + "e = " + str(e) + " s = " + str(s) + " k = " + str(k) + " l = " + str(l) + str(" FAILED"))
                    f += 1
    bar.next()
bar.finish()
print("PASSED: " + str(p))
print("FAILED: " + str(f))