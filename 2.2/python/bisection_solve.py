from math import log
import time


def bisection_solve(func, root_range, accuracy=1e-8):
    a, b = root_range
    if func(a) == 0:
        return a
    elif func(b) == 0:
        return b
    assert func(a)*func(b) < 0, "root range error!"
    c = (a+b)/2
    while b - c > accuracy:
        if func(c) == 0:
            return c
        elif func(c)*func(b) < 0:
            a = c
        else:
            b = c
        c = (a + b) / 2
    return c


if __name__ == '__main__':
    start = time.clock()
    print(bisection_solve(log, (0.5, 2)))
    end = time.clock()
    print("solve time:", end - start, 's')
