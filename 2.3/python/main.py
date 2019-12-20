from numpy import sqrt, sin, cos, pi
import numpy
from ran import Ran


def statistical_test():
    ran = Ran()
    x = []
    y = []
    for i in range(1000):
        x.append(ran.float_rand())
        y.append(ran.float_rand())
    plt.plot(x, y, 'r.')
    plt.show()


def discrete_wander(step):
    ran = Ran()
    x, y = 0, 0
    xl = [0]
    yl = [0]
    print("n=0:\tr=(0,0),R=0")
    for i in range(step):
        r = ran.float_rand()
        if r < 1/4:
            x += 1
        elif r < 1/2:
            x -= 1
        elif r < 3/4:
            y += 1
        else:
            y -= 1
        print("n=%s:\tr=(%s,%s),R=%s" % (i+1, x, y, sqrt(x**2+y**2)))
        xl.append(x)
        yl.append(y)
    plt.plot(xl, yl, marker='o')
    plt.annotate('final position', xy=(xl[-1], yl[-1]), xytext=(xl[-1]+1/2, yl[-1]+1/2),
                 arrowprops=dict(arrowstyle='->'))
    plt.annotate('start position', xy=(0, 0), xytext=(1/2, 1/2),
                 arrowprops=dict(arrowstyle='->'))
    plt.show()


def continuous_wander_with_plot(step):
    ran = Ran()
    x, y = 0, 0
    xl = [0]
    yl = [0]
    for i in range(step):
        theta = ran.float_rand()*2*pi
        x += cos(theta)
        y += sin(theta)
        xl.append(x)
        yl.append(y)
    plt.plot(xl, yl, marker='o')
    plt.annotate('final position', xy=(xl[-1], yl[-1]), xytext=(xl[-1] + 1 / 2, yl[-1] + 1 / 2),
                 arrowprops=dict(arrowstyle='->'))
    plt.annotate('start position', xy=(0, 0), xytext=(1 / 2, 1 / 2),
                 arrowprops=dict(arrowstyle='->'))
    plt.show()


def continuous_wander(step, seed):
    ran = Ran(seed)
    x, y = 0, 0
    rl = [0]
    for i in range(step):
        theta = ran.float_rand() * 2 * pi
        x += cos(theta)
        y += sin(theta)
        rl.append(sqrt(x ** 2 + y ** 2)**2)
    return rl


def expectation(n, step):
    expectation_of_r = numpy.zeros(step+1)
    ran = Ran()
    for i in range(n):
        expectation_of_r += numpy.array(continuous_wander(step, ran.int_rand()))
    expectation_of_r /= n
    plt.plot(expectation_of_r, marker='.')
    plt.xlabel('N')
    plt.ylabel('E(R)')
    nl = numpy.arange(step+1)**2
    a = sum(expectation_of_r*sqrt(nl))/sum(nl)
    fyl = a*sqrt(nl)
    plt.plot(fyl, label="a=%s" % a)
    plt.legend()
    plt.show()
    return expectation_of_r


pass
