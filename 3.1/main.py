from ran import Ran
from numpy import log, cos, sin, sqrt, pi
import numpy
from Newton_solve import ergodic_find_roots
from pylab import hist2d, hist
import pylab


def gaussian(ran_object, sd=1, mean=0):
    u = ran_object.float_rand()
    v = ran_object.float_rand()
    return sqrt(-2*log(u))*cos(2*pi*v)*sd+mean


def test():
    ran = Ran()
    x = []
    y = []
    for i in range(10000):
        point = w_dist(ran, 0.08)
        x.append(point[0])
        y.append(point[1])
    hist2d(x, y, bins=100)
    pylab.colorbar()
    pylab.show()


def e_dist(ran_object, e_max):
    alpha = 1-1/sqrt(1+e_max)
    while True:
        u = ran_object.float_rand()
        e = 1/(1-alpha*u)**2-1

        def g_func(e):
            return numpy.e**(-2/3/e)*((e+1)/e)**(3/2)

        if e_max > 0.8:
            g = g_func(e)/g_func(0.8)
        else:
            g = g_func(e)/g_func(e_max)
        if ran_object.float_rand() <= g:
            return e


def w_dist(ran_object, e_max=10):
    e = e_dist(ran_object, e_max)
    v = gaussian(ran_object, sqrt(e/2))
    return e, v


def A_line(t):
    return numpy.array((cos(pi/2*t/2/T)**2*A0*sin(omega*t), 0))


def A_ellipse(t):
    return cos(pi/2*t/2/T)**2*A0/sqrt(5)*numpy.array((sin(omega * t)*2, cos(omega * t)))


def D(func, t, delta_t=0.1):
    return (func(t+delta_t/2)-func(t-delta_t/2))/delta_t


def final_state(start_time, v0):
    step = int((2*T-start_time)//2+1)
    delta_t = (2*110.23-start_time)/step
    A = A_ellipse
    E0 = -D(A, start_time)
    r = -0.5/sum(E0**2)*E0
    v0d = numpy.array((-r[1], r[0]))
    v = v0*v0d/sqrt(sum(v0d**2))
    for i in range(step):
        t = start_time+i*delta_t
        next_v = v+A(t+delta_t)-A(t)-r/sqrt(sum(r**2)+0.04)**3*delta_t
        next_r = r+(v+next_v)/2*delta_t
        for j in range(1):
            next_v = v+A(t+delta_t)-A(t)-(r/sqrt(sum(r**2)+0.04)**3+next_r/sqrt(sum(next_r**2)+0.04)**3)/2*delta_t
            next_r = r+(v+next_v)/2*delta_t
        r = next_r
        v = next_v
    return r, v


def p_dist():
    ran = Ran()
    e_max_l = 0.07553
    e_max_e = 0.03410
    e_max = e_max_e
    A = A_ellipse
    x = []
    y = []
    for i in range(100):
        e, v0 = w_dist(ran, e_max)

        def func(t):
            return sqrt(sum(D(A, t)**2))-e

        roots = ergodic_find_roots(func, (-2*T, 2*T), 0.1, 0.01)
        if not roots:
            continue
        start_time = roots[int(ran.float_rand()*len(roots)//len(roots))]
        rf, pf = final_state(start_time, v0)
        L = rf[0]*pf[1]-rf[1]*pf[0]
        if sum(pf**2)-2/sqrt(sum(rf**2)) < 0:
            continue
        p_inf = sqrt(sum(pf**2)-2/sqrt(sum(rf**2)))
        a = numpy.array((pf[1]*L-rf[0]/sqrt(sum(rf**2)), -pf[0]*L-rf[1]/sqrt(sum(rf**2))))
        p_inf = p_inf/(1+p_inf**2*L**2)*(numpy.array((-a[1], a[0]))*L*p_inf-a)
        x.append(p_inf[0])
        y.append(p_inf[1])
    hist2d(x, y, bins=150, range=[[-1.5, 1.5], [-1.5, 1.5]])
    pylab.colorbar()
    pylab.show()
    hist(x, bins=150, range=[-1.5, 1.5])
    pylab.show()
    hist(y, bins=150, range=[-1.5, 1.5])
    pylab.show()


A0 = 1.325
omega = 0.057
T = 110.23
p_dist()