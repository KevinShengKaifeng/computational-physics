def Simpson_int(f, a, b):
    return (b-a)/6*(f(a)+4*f((a+b)/2)+f(b))


def Simpson_integrate(f, a, b, N=200):
    if b == a:
        return 0
    integral = 0
    delta = (b-a)/N
    for i in range(N//2):
        integral += Simpson_int(f, a+2*i*delta, a+2*(i+1)*delta)
    if N % 2 == 1:
        integral += Simpson_int(f, b-delta, b)
    return integral
