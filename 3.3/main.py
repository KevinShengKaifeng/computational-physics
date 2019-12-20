import numpy
from numpy import array, sqrt, zeros, ones, sin, cos, pi, arange, log
from power_method import tridiagonal_inverse_power_solve, tridiagonal_dot
from tridiagonal_bisection import bisection_solve
from chase_solve import chase_solve
from matplotlib import pyplot as plt


def stationary_state(k, x_max=100):
    delta_x = 0.1
    n = 2*int(x_max/delta_x)
    a = zeros(n)
    b = zeros(n)
    x = []
    for i in range(n):
        a[i] = 1/delta_x**2-1/sqrt(1+((-n//2+i)*delta_x)**2)
        b[i] = -1/2/delta_x**2
        x.append((-n//2+i)*delta_x)
    b[0] = 0
    energy = bisection_solve(a, b, k)
    a -= ones(n)*energy
    wave_func = tridiagonal_inverse_power_solve(a, b)[1]
    wave_func /= sqrt(wave_integrate(wave_func, wave_func))
#    plt.plot(x, array(wave_func)**2)
#    plt.show()
    return energy, array(wave_func)


def wave_integrate(w1, w2, delta_x=0.1):
    s = 0
    for i in range(len(w1)):
        s += w1[i].conjugate()*w2[i]
    return s*delta_x


def attosecond_laser():
    delta_x = 0.1
    x_max = 100
    n = 2*int(x_max/delta_x)
    omega = 1
    E0 = sqrt(1e13/3.5094448314e16)
    N = 20
    t_max = N*2*pi/omega
    delta_t = 0.05
    t_step = int(t_max/delta_t)
    ground_wave_func = stationary_state(1, x_max)[1]
    wave_func = ground_wave_func.copy()
    p = []
    t = []
    x = []
    for i in range(n):
        x0 = (-n//2+i)*delta_x
        x.append(x0)
    x = array(x)
    for j in range(t_step):
        time = j*delta_t
        a = ones(n)
        b = ones(n)
        a *= 1/delta_x**2
        a += -1/sqrt(1+x**2)+x*E0*cos(omega*time/2/N)**2*sin(omega*time)
        b *= -1/2/delta_x**2
        b[0] = 0
        wave_func = tridiagonal_dot(ones(n)-1j*delta_t/2*a, -1j*delta_t/2*b, wave_func)
        wave_func = chase_solve(1j*delta_t/2*b, ones(n)+1j*delta_t/2*a, 1j*delta_t/2*b, wave_func)
#        if j%100 == 0:
#            plt.plot(x, abs(array(wave_func)) ** 2)
#            plt.show()
        t.append(time)
        p.append(abs(wave_integrate(ground_wave_func, wave_func)))
    plt.plot(t, p)
    plt.show()
    ionize_wave_func = wave_func-abs(wave_integrate(ground_wave_func, wave_func))*ground_wave_func
    k = arange(-1.5, 1.5, 0.01)
    pk = []
    for i in k:
        wk = numpy.e**(1j*i*x)
        pk.append(abs(wave_integrate(wk, ionize_wave_func)))
    plt.plot(k, pk)
    plt.show()


def infrared_laser():
    delta_x = 0.1
    x_max = 200
    n = 2*int(x_max/delta_x)
    omega = 45.5633525316/300
    E0 = sqrt(2e14/3.5094448314e16)
    N = 50
    t_max = N*2*pi/omega/2
    delta_t = 0.05
    t_step = int(t_max/delta_t)
    ground_wave_func = stationary_state(1, x_max)[1]
    wave_func = ground_wave_func.copy()
    t = []
    x = []
    for i in range(n):
        x0 = (-n//2+i)*delta_x
        x.append(x0)
    x = array(x)
    d = []
    acc = []
    absorb_func = numpy.e ** (-(abs(x) - 0.75 * x_max) ** 2 / 0.2 ** 2)
    for i in range(int(0.25 * n / 2), int(1.75 * n / 2) + 1):
        absorb_func[i] = 1
    for j in range(t_step):
        time = j*delta_t
        a = ones(n)
        b = ones(n)
        a *= 1/delta_x**2
        E = E0*cos(omega*time/2/N)**2*sin(omega*time)
        a += -1/sqrt(1+x**2)+x*E
        b *= -1/2/delta_x**2
        b[0] = 0
        wave_func = tridiagonal_dot(ones(n)-1j*delta_t/2*a, -1j*delta_t/2*b, wave_func)
        wave_func = chase_solve(1j*delta_t/2*b, ones(n)+1j*delta_t/2*a, 1j*delta_t/2*b, wave_func)
        wave_func *= absorb_func
#        if j%100 == 0:
#            plt.plot(x, abs(array(wave_func)) ** 2)
#            plt.show()
        d.append(wave_integrate(wave_func, x*wave_func))
        acc.append(wave_integrate(wave_func, (-x/(x**2+1)**(3/2)+E)*wave_func))
        t.append(time)
    t = array(t)
    w = arange(1, 15, 0.01)
    w *= omega
    plt.plot(t, d)
    plt.show()
    plt.plot(t, acc)
    plt.show()
    Awd = []
    Awa = []
    for i in w:
        wt = numpy.e**(-1j*i*t)
        Awd.append(abs(wave_integrate(wt, d, delta_t)))
        Awa.append(abs(wave_integrate(wt, acc, delta_t)))
    Awd *= w**2/sqrt(2*pi)
    Awa /= sqrt(2*pi)
    plt.plot(arange(1, 15, 0.01), 2*log(Awd)/log(10), label="dipole moment")
    plt.plot(arange(1, 15, 0.01), 2*log(Awa)/log(10), label="accelerate")
    plt.legend()
    plt.show()


infrared_laser()
