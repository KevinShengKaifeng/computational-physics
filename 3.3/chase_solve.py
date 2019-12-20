def chase_solve(a, b, c, f):
    for i in range(len(a)-1):
        c[i+1] = c[i+1]/b[i]
        b[i+1] = b[i+1]-c[i+1]*a[i+1]
        f[i+1] = f[i+1]-c[i+1]*f[i]
    a[-1] = f[-1]/b[-1]
    for i in range(len(a)-1):
        a[-2-i] = f[-2-i]/b[-2-i]-c[-1-i]*a[-1-i]
    return a


def test1():
    b = [-2,-2,-2,-2]
    a = [1,1,1,1]
    c = [1,1,1,1]
    f = [0,-1,-1,0]
    print(chase_solve(a,b,c,f))


if __name__ == "__main__":
    test1()
