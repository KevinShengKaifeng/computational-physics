from numpy import sin, pi, sqrt, log, cos, cosh, tanh
import numpy
import matplotlib.pyplot as plt


def Ek(k, q, p):
    Q = 0
    for j in range(1, len(q)-1):
        Q += sin(pi*k*j/(len(q)-1))*q[j]
    Q *= sqrt(2/(len(q)-2))
    P = 0
    for j in range(1, len(q)-1):
        P += sin(pi*k*j/(len(q)-1))*p[j]
    P *= sqrt(2/(len(q)-2))
    omegak = 2*sin(pi*k/2/(len(q)-1))
    return (omegak**2*Q**2+P**2)/2


def set_qp(Q0_list, n):
    q = numpy.zeros(n + 2)
    p = numpy.zeros(n + 2)
    for k in range(len(Q0_list)):
        for j in range(1, 1 + n):
            q[j] += sqrt(2 / (n)) * sin(pi * (k + 1) * j / (n + 1)) * Q0_list[k]
    return q, p


def set_orphan_qp(n):
    qorphan = numpy.zeros(n+2)
    porphan = numpy.zeros(n+2)
    B, k0 = 0.5, 11
    omegak0 = 2*sin(pi*k0/2/(n+1))
    for i in range(1, 1+n):
        qorphan[i] = B*cos(pi*k0*(i-n/2)/(n+1))/cosh(sqrt(3/2)*B*omegak0*(i-n/2))
        porphan[i] = B/cosh(sqrt(3/2)*B*omegak0*(i-n/2))*(omegak0*(1+3/16*omegak0**2*B**2)*sin(pi*k0*(i-n/2)/(n+1))
            +sqrt(3/2)*B*cos(pi*k0*(i-n/2)/(n+1))*sin(pi*k0/(n+1))*tanh(sqrt(3/2)*B*omegak0*(i-n/2)))
    return qorphan, porphan


def alpha_FPU(alpha, n, q, p):
    omega1 = 2*sin(pi/2/(n+1))
    N = 400
    tm = N*2*pi/omega1
    step = N*500
    delta_t = tm/step
    next_p = p.copy()
    next_q = q.copy()
    E = [[],[],[],[]]
    tlist = []
    for i in range(step):
        for j in range(1, n+1):
            next_p[j] = p[j]+delta_t*((q[j+1]-2*q[j]+q[j-1])-alpha*(q[j]-q[j+1])**2+alpha*(q[j-1]-q[j])**2)
        next_q = q+delta_t/2*(p+next_p)
        for m in range(1):
            for j in range(1, n+1):
                next_p[j] = p[j]+delta_t/2*((q[j+1]-2*q[j]+q[j-1])-alpha*(q[j]-q[j+1])**2+alpha*(q[j-1]-q[j])**2+
                    (next_q[j+1]-2*next_q[j]+next_q[j-1])-alpha*(next_q[j]-next_q[j+1])**2+alpha*(next_q[j-1]-next_q[j])**2)
            next_q = q+delta_t/2*(p+next_p)
        p = next_p.copy()
        q = next_q.copy()
        for m in range(4):
            E[m].append(Ek(m+1, q, p))
        tlist.append(i*N/step)
#        if i%10==0:
#            plt.plot(q)
#            plt.show()
    for i in range(4):
        plt.plot(tlist, E[i], label="k=%s" % i+1)
    plt.legend()
    plt.show()


def beta_FPU(beta, n, q, p):
    omega1 = 2*sin(pi/2/(n+1))
    N = 500
    tm = N*2*pi/omega1
    step = N*500
    delta_t = tm/step
    next_p = p.copy()
    next_q = q.copy()
    E = [[],[],[],[],[]]
    tlist = []
    for i in range(step):
        for j in range(1, n+1):
            next_p[j] = p[j]+delta_t*((q[j+1]-2*q[j]+q[j-1])-beta*(q[j]-q[j+1])**3+beta*(q[j-1]-q[j])**3)
        next_q = q+delta_t/2*(p+next_p)
        for m in range(1):
            for j in range(1, n+1):
                next_p[j] = p[j]+delta_t/2*((q[j+1]-2*q[j]+q[j-1])-beta*(q[j]-q[j+1])**3+beta*(q[j-1]-q[j])**3+
                    (next_q[j+1]-2*next_q[j]+next_q[j-1])-beta*(next_q[j]-next_q[j+1])**3+beta*(next_q[j-1]-next_q[j])**3)
            next_q = q+delta_t/2*(p+next_p)
        p = next_p.copy()
        q = next_q.copy()
        for m in range(5):
            E[m].append(log(Ek(m+9, q, p)))
        tlist.append((i+1)*N/step)
        if i%500==0:
            plt.plot(q)
            plt.show()
#    for i in range(3):
#        plt.plot(tlist, E[i], label="k=%s" % (2*i+1))
    for i in range(5):
        plt.plot(tlist[500:], E[i][500:], label="k=%s" % (i+9))
    plt.legend()
    plt.show()


Q0_list = numpy.zeros(11)
Q0_list[-1] = 1
#alpha_FPU(0.25, 32, *set_qp([20], 32))
beta_FPU(1, 128, *set_orphan_qp(128))
