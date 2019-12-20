def sign(a, b, x):
    s = 0
    q = a[0] - x
    for i in range(len(a)):
        if q < 0:
            s += 1
        if i < len(a)-1:
            q = a[i+1]-x-b[i+1]**2/q
    return s


def bisection_solve(a, b, m):
    u = 0
    n = len(a)
    for i in range(n):
        if i == n-1:
            s = abs(a[-1])+abs(b[-1])
        else:
            s = abs(a[i]) + abs(b[i]) + abs(b[i + 1])
        if u < s:
            u = s
    l = -u
    step = 100
    accuracy = 0.001
    for i in range(step):
        r = (l+u)/2
        if sign(a, b, r) >= m:
            u = r
        else:
            l = r
        if (u-l) < accuracy:
            break
    return (l+u)/2
