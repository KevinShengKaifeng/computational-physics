import numpy as np
import matplotlib.pyplot as plt


def complex_func_plot(func, a, b, N=100):
    x = [a]
    for i in range(N):
        x.append(a+(b-a)/N*(i+1))
    y = []
    z = []
    iscomplex = False
    for i in x:
        v = func(i)
        y.append(v.real)
        z.append(v.imag)
        if v.imag != 0:
            iscomplex = True
    plt.plot(x, y)
    if iscomplex:
        plt.plot(x, z)
    plt.show()
#    plt.savefig("fig.eps")


def funcs_plot(*func, xrange, xname='', yname='', N=100):
    a, b = xrange
    x = [a]
    for i in range(N):
        x.append(a + (b - a) / N * (i + 1))
    ys = []
    for f in func:
        y = []
        for i in x:
            y.append(f(i))
        ys.append(y)
        if f.__name__ != "<lambda>":
            plt.plot(x, y, label=f.__name__)
        else:
            plt.plot(x, y)
    plt.xlabel(xname)
    plt.ylabel(yname)
    plt.legend()
    plt.show()


if __name__ == '__main__':
    complex_func_plot(lambda x: np.e**(1j*x), 0, 2*np.pi)
