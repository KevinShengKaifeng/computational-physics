from Gaussian_elimination import gelimination_solve


def h(i, x):
    return x[i+1]-x[i]


def splinei(x0, i, x, y, M):
    return (M[i]*(x[i+1]-x0)**3+M[i+1]*(x0-x[i])**3)/6/(h(i, x))\
        + (y[i]-M[i]*h(i, x)**2/6)*(x[i+1]-x0)/h(i, x)\
        + (y[i+1]-M[i+1]*h(i, x)**2/6)*(x0-x[i])/h(i, x)


def spline(x0, x, y, M):
    n = len(x) - 1
    if x0 < x[0]:
        return splinei(x0, 0, x, y, M)
    for i in range(n):
        if x0 < x[i+1]:
            return splinei(x0, i, x, y, M)
    return splinei(x0, n-1, x, y, M)


def f(x, y, *n):
    if len(n) == 1:
        return y[n[0]]
    else:
        return (f(x, y, *n[:-2], n[-1])-f(x, y, *n[:-1]))/(x[n[-1]]-x[n[-2]])


def empty_matrix(size):
    a = []
    for i in range(size):
        ai = []
        for j in range(size):
            ai.append(0)
        a.append(ai)
    return a


def spline_interpolation(x, y):
    n = len(x) - 1
    miu = [0]
    lamb = [0]
    d = [0]
    for i in range(n - 1):
        miu.append((h(i, x)) / (x[i + 2] - x[i]))
        lamb.append((h(i + 1, x)) / (x[i + 2] - x[i]))
        d.append(6 * f(x, y, i, i + 1, i + 2))
    miu.append(0)
    d.append(0)
    A = empty_matrix(n + 1)
    for i in range(n):
        A[i][i + 1] = lamb[i]
        A[i][i] = 2
        A[i + 1][i] = miu[i + 1]
    A[n][n] = 2
    M = gelimination_solve(A, d)
    return lambda x0: spline(x0, x, y, M)
