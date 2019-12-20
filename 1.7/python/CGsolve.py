from math import sqrt


def typeof(a):
    if isinstance(a, list):
        if isinstance(a[0], list):
            # a is a matrix
            return 2
        elif isinstance(a[0], dict):
            # a is a scarce matrix
            return 2.5
        else:
            # a is a vector
            return 1
    elif isinstance(a, dict):
        # a is a scarce vector
        return 1.5
    else:
        # a is a scalar
        return 0


def iter_product(a, b):
    # 返回c，c为b与a的每项内积的结果
    assert typeof(a) > typeof(b), 'iter_product error!'
    c = []
    for i in a:
        c.append(inner_product(i, b))
    return c


def inner_product(a, b):
    if typeof(a) < typeof(b):
        return inner_product(b, a)
    if typeof(a) == 0:
        return a*b
    elif typeof(a) == 1:
        if typeof(b) == 0:
            return iter_product(a, b)
        else:
            assert len(a) == len(b), 'vector size error!'
            c = 0
            for i in range(len(a)):
                c += a[i]*b[i]
            return c
    elif typeof(a) == 1.5:
        if typeof(b) == 0:
            c = a.copy()
            for i in c:
                c[i] *= b
            return c
        elif typeof(b) == 1:
            c = 0
            for i in a:
                c += a[i]*b[i]
            return c
        elif typeof(b) ==1.5:
            c = 0
            for i in a:
                if i in b:
                    c += a[i]*b[i]
            return c
    elif typeof(a) == 2:
        if typeof(b) < typeof(a):
            return iter_product(a, b)
        elif typeof(b) == 2:
            assert len(a[0]) == len(b), 'matrix size error!'
            c = []
            for i in range(len(a)):
                line = []
                for j in range(len(a)):
                    line.append(0)
                    for k in range(len(a)):
                        line[j] += a[i][k] * b[k][j]
                c.append(line)
            return c
    elif typeof(a) == 2.5:
        c = a.copy()
        if typeof(b) <= 1.5:
            return iter_product(a, b)
        elif typeof(b) == 2:
            c = []
            for i in range(len(a)):
                line = []
                for j in range(len(a)):
                    line.append(0)
                    for k in a[i]:
                        line[j] += a[i][k] * b[k][j]
                c.append(line)
            return c
        elif typeof(b) == 2.5:
            pass


def vector_add(a, b):
    if typeof(a) < typeof(b):
        return vector_add(b, a)
    assert (typeof(a)-typeof(b)) < 1, 'vector_add type error!'
    if typeof(a) == 1:
        c = []
        for i in range(len(a)):
            c.append(a[i]+b[i])
        return c
    elif typeof(a) == 1.5:
        c = b.copy()
        if typeof(b) == 1.5:
            for i in a:
                if i in c:
                    c[i] += a[i]
                else:
                    c[i] = a[i]
        elif typeof(b) == 1:
            for i in a:
                c[i] += a[i]
        return c
    elif (typeof(a) == 2)or(typeof(a) == 2.5):
        c = []
        for i in range(len(a)):
            c.append(vector_add(a[i], b[i]))
        return c


def norm(a):
    if typeof(a) <= 1.5:
        return sqrt(inner_product(a, a))


def CGsolve(A, b, limit=0):
    x = b.copy()
    for i in range(len(x)):
        x[i] = 0
    r = b.copy()
    p = r.copy()
    for i in range(len(b)):
        alpha = inner_product(r, r)/inner_product(p, inner_product(A, p))
        x = vector_add(x, inner_product(p, alpha))
        r1 = vector_add(r, inner_product(inner_product(A, p), -alpha))
        beta = inner_product(r1, r1)/inner_product(r, r)
        r = r1
        p = vector_add(r1, inner_product(p, beta))
        if norm(inner_product(p, alpha)) <= limit:
            print("end after %s iteration" % (i+1))
            break
    return x
