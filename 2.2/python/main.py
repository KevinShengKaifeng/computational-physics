from CGsolve import vector_add, inner_product
from Newton_solve import newton_solve, newton_bisection_solve, D
from math import sqrt, e, pi, factorial
from bisection_solve import bisection_solve
from plot import funcs_plot
from QR_method import pragmatic_QR_method


def polynomial_value(clist, x):
    value = 0
    for i in clist:
        value = value*x+i
    return value


# H(n+1)=x*H(n)-n/2*H(n-1)
def Hermite_coefficient(n):
    if n == 0:
        return [1]
    elif n == 1:
        return [1, 0]
    list0 = [1]
    list1 = [1, 0]
    for i in range(n-1):
        nlist = vector_add(list1+[0], inner_product(-(i+1)/2, [0, 0]+list0))
        list0 = list1
        list1 = nlist
    return list1


def recursive_Hermite(n, x0):
    if n == 0:
        return 1
    elif n == 1:
        return x0
    x1, x2 = 1, x0
    for i in range(1, n):
        x = x0*x2-i/2*x1
        x1 = x2
        x2 = x
    return x2


def remove_root(clist, x0):
    nclist = [0]
    for i in clist:
        nclist.append(nclist[-1]*x0+i)
    if abs(nclist[-1]) > 1e-5:
        print("polynomial don't have the given root!")
    return nclist[1:-1]


def range_finder(func, origin_range):
    a, b = origin_range
    f = func(b)
    if f*func(a) <= 0:
        return origin_range
    N = 15
    for i in range(1, N):
        for j in range(1, 2**i):
            c = a+(b-a)/2**i*j
            fc = func(c)
            if fc*f < 0:
                return a+(b-a)/2**i*(j-1), c
    print("can't find root range!")


# 已知xn<sqrt(2n)
def polynomial_find_roots(clist, root_range=()):
    roots = []
    if not root_range:
        root_range = (0, sqrt(2 * len(clist) - 2))
    search_time = 200*len(clist)
    while clist[-1] == 0:
        roots.append(0)
        clist = clist[:-1]
    ifeven = True
    for i in range((len(clist)-1)//2):
        if clist[1+2*i] != 0:
            ifeven = False
            break
    while len(clist) > 1:

        def poly(x):
            return polynomial_value(clist, x)

#        funcs_plot(poly, xrange=(root_range[0],root_range[1]-2), N=300)
        x0 = newton_solve(poly, 0, root_range, search_time)
        roots.append(x0)
        clist = remove_root(clist, x0)
        if ifeven:
            roots.append(-x0)
            clist = remove_root(clist, -x0)
    roots.sort()
    return roots


def pick_find_roots(func, n, root_range=()):
    roots = []
    if not root_range:
        root_range = (0.1, sqrt(2 * n))
    search_time = 200*n
    while True:
        def nfunc(x):
            y = func(x)
            for i in roots:
                y /= (x-i)
            return y
        x0 = newton_solve(nfunc, 0.1, root_range, search_time)
        if x0 is None:
            print(len(roots))
            roots.sort()
            return roots
        roots.append(x0)


def ergodic_find_roots(func, root_range, deltax=1e-1, accuracy=1e-8):
    a, b = root_range
    if func(a) >= 0:
        sign = 1
    else:
        sign = -1
    ranges = []
    x = a
    while x <= b:
        if func(x)*sign <= 0:
            ranges.append((x-deltax, x))
            sign *= -1
        x += deltax
    print(len(ranges))
    roots = []
    for i in ranges:
        root = newton_bisection_solve(func, i, accuracy)
        roots.append(root)
    return roots


def energy_and_wave_function():
    def phi(n, x):
        return recursive_Hermite(n, x)*e**(-x**2/2)*2**(n/2)/pi**(1/4)/sqrt(factorial(n))

    def base_func(i, x):
        return (-1)**i*phi(n1, x)/(x-result[i])/D(result[i], lambda y: phi(n1, y))

    def Phi(k, z):
        sum = 0
        for i in range(n1):
            sum += base_func(i, z/s)*re[1][k][i]
        return sum

    s = 0.6
    for n1 in [12, 24, 48, 96]:
        for n2 in [1, n1//4, n1//3]:
            ran = (-sqrt(2 * n1), sqrt(2 * n1))
            result = ergodic_find_roots(lambda x: recursive_Hermite(n1, x), ran)
            matrix = []
            for i in range(n1):
                row = []
                for j in range(n1):
                    if i == j:
                        row.append((4*n1-1-2*result[i]**2)/6/s**2+(s*result[i])**2+(s*result[i])**4)
                    else:
                        row.append((-1)**(i-j)*(2/(result[i]-result[j])**2-1/2)/s**2)
                matrix.append(row)
            re = pragmatic_QR_method(matrix)
            print("N=%s,n=%s,En=%s" % (n1, n2, re[0][n2-1]))
            funcs_plot(lambda x: Phi(n2, x), xrange=ran, xname="z", yname="Phi", N=401)


def hermite_roots():
    for N in [24,48,96,244]:
        ran = (-sqrt(2 * N), sqrt(2 * N))
        result = ergodic_find_roots(lambda x: recursive_Hermite(N, x), ran)
        print(N, result[N//2], result[-1])


energy_and_wave_function()
