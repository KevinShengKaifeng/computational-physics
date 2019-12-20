from math import e, sqrt
from subsection_Simpson_integrate import Simpson_integrate


def A(r):
    return e**(-r/3)*r*(6-r)


def D(f, deltar, n=1):
    return lambda ri: (f(ri+deltar)-f(ri-deltar))/2/deltar


def func(x, deltar=0.005):
    k = 2*sqrt(6)/3**5
    return k**2*(A(x)**2*(1-x)-A(x)/2*D(lambda r: r**2*D(A, deltar)(r), deltar)(x))


print(Simpson_integrate(func, 0, 60, 6000))
