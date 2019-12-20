from math import sqrt, sinh, cosh, cos
import math
from plot import funcs_plot
from subsection_Simpson_integrate import Simpson_integrate
from spline_interpolation import spline_interpolation
from Newton_solve import newton_solve, D


def sin(x):
    if type(x) == complex:
        return sin(x.real)*cosh(x.imag)+cos(x.real)*sinh(x.imag)*1j
    elif type(x) == int or type(x) == float:
        return math.sin(x)


def A(t):
    return A0*sin(omega*t)*sin(omega*t/4)**2


# 鞍点方程(pz+A0sin(omega*t)sin(omega*t/4)^2)^2/2+Ip=0
def saddle_f(t, pz=1):
    return (A(t)+pz)**2/2+Ip


def bind_pz(saddle_f, pz):
    def saddle_f_pz(x):
        return saddle_f(x, pz)
    return saddle_f_pz


def find_roots(func, T, *, n=50):
    roots = []
    for i in range(n):
        r = newton_solve(func, 5j+i/n*T)
        if not r:
            continue
        if r.imag < 0:
            continue
        while r.real < 0:
            r += T
        while r.real > T:
            r -= T
        isincluded = False
        for ri in roots:
            if abs(ri-r) < 1e-6:
                isincluded = True
        if not isincluded:
            roots.append(r)
    return roots


def S(t, pz):
    saddle_f_pz = bind_pz(saddle_f, pz)
    return Simpson_integrate(saddle_f_pz, 0, t, 200)


def find_saddle_points(pz):
    n0 = 50
    roots = []
    saddle_f_pz = bind_pz(saddle_f, pz)
    while len(roots) != 6:
        roots = find_roots(saddle_f_pz, 4 * math.pi / omega, n=n0)
        n0 += 5
        if n0 > 100:
            print("roots finding error!")
            break
    return roots


def M_saddle_point_method(Ek):
    pz = sqrt(2*Ek)
    saddle_f_pz = bind_pz(saddle_f, pz)
    roots = find_saddle_points(pz)
    suma = 0
    for t in roots:
        suma += math.e**(1j*S(t, pz))/D(t, saddle_f_pz)
    return -(2*Ip)**(5/4)/sqrt(2)*suma


def M_square_saddle_point(Ek):
    return abs(M_saddle_point_method(Ek)) ** 2


def M_square_saddle_point_p(pz):
    return M_square_saddle_point(pz**2/2)


def M_direct_integrate(Ek):
    pz = sqrt(2 * Ek)
    x = []
    y = []
    for i in range(100):
        x.append(i/100*4*math.pi/omega)
        y.append(S(x[-1], pz))
    St = spline_interpolation(x, y)
#    funcs_plot(St, xrange=(0, 4*math.pi/omega), N=100)

    def integrate_func(t):
        r = -(pz+A(t))*D(t, A)/((pz+A(t))**2+2*Ip)**3*math.e**(1j*St(t))
        return r
    return 2**(7/2)*(2*Ip)**(5/4)/math.pi*Simpson_integrate(integrate_func, 0, 4*math.pi/omega, 2000)


def M_square_direct_integrate(Ek):
    return abs(M_direct_integrate(Ek)) ** 2


def M_square_direct_integrate_p(pz):
    return M_square_direct_integrate(pz**2/2)


def ratio(Ek):
    return M_square_direct_integrate(Ek)/M_square_saddle_point(Ek)


lam = 1600
I0 = 2e14
E0 = sqrt(I0/3.5094448314e16)
Ip = 13.6/27.2
omega = 45.5633525316/lam
A0 = E0/omega
roots = find_saddle_points(1)
with open("output.txt", 'w') as output:
    for i in roots:
        output.write(str(i))
        output.write('\n')
funcs_plot(ratio,
            xrange=(0, 2), xname='Ek(a.u.)', yname='ratio', N=400)
#funcs_plot(M_square_saddle_point, M_square_direct_integrate,
#            xrange=(0, 2), xname='Ek(a.u.)', yname='|M0p(tf,ti)|^2', N=400)
#funcs_plot(M_square_saddle_point_p, M_square_direct_integrate_p,
#           xrange=(0, 1), xname='pz(a.u.)', yname='|M0p(tf,ti)|^2', N=400)
