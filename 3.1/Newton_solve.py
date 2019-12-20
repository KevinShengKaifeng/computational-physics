from math import sin
import time


def D(x, func, delta=1e-3):
    return (func(x+delta/2)-func(x-delta/2))/delta


def newton_solve_fixed_x0(func, x0):
    max_step = 30
    deltamin = 1e-8
    deltaxmax = abs(x0)+100
    rndeltax = abs(x0)+1
    derivative_min = 1e-3
    rnderivative = 2e-5
    i = 0
    while i < max_step:
        d = D(x0, func)
        if abs(d) < derivative_min:
            d = rnderivative
        deltax = -func(x0) / d
        if abs(deltax) < deltamin:
            return x0
        if abs(deltax) > deltaxmax:
            deltax = rndeltax*deltax/abs(deltax)
        x0 += deltax
        i += 1
    return


def newton_solve(func, x0=0, root_range=(), find_time=10):
    solve = newton_solve_fixed_x0(func, x0)
    if solve is not None:
        return solve
    elif not root_range:
        print("root not found!")
        return
    for i in range(find_time):
        solve = newton_solve_fixed_x0(func, root_range[0]+(root_range[1]-root_range[0])/find_time*(i+1))
        if solve is not None:
            return solve
    print("root not found!")
    return


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
    roots = []
    for i in ranges:
        root = newton_bisection_solve(func, i, accuracy)
        roots.append(root)
    return roots


def newton_bisection_solve(func, root_range, accuracy=1e-8):
    a, b = root_range
    if func(a) == 0:
        return a
    elif func(b) == 0:
        return b
    assert func(a) * func(b) < 0, "root range error!"
    c = (a+b)/2
    while b-c > accuracy:
        if func(c) == 0:
            return c
        elif func(c)*func(b) < 0:
            if c-a < accuracy:
                return c
            a = c
        else:
            b = c
        d = D(c, func)
        if d == 0:
            c = (a + b) / 2
        c = a - func(a) / d
        if c < a or c > b:
            c = (a + b) / 2
    return c


if __name__ == '__main__':
    start = time.clock()
    print(newton_bisection_solve(sin, (1, 4)))
    end = time.clock()
    print("solve time:", end - start, 's')
