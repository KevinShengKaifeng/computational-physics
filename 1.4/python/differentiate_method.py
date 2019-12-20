from CGsolve import vector_add, inner_product
from math import sqrt, pi, cos, sin


def legendre_coefficient(l):
    if l == 0:
        return [1]
    elif l == 1:
        return [1, 0]
    else:
        last_two = [[1], [1, 0]]
        for i in range(1, l):
            last = last_two[1]
            last_two[1].append(0)
            last_two[1] = vector_add(inner_product(last_two[1],(2*i+1)/(i+1)),inner_product([0, 0]+last_two[0], -i/(i+1)))
            last_two[0] = last
        return last_two[1]
        
        
def associated_legendre_coefficient(l, m):
    if m == 0:
        return legendre_coefficient(l)
    coefficient = legendre_coefficient(l)
    for i in range(m):
        newc = []
        order = len(coefficient)-1
        for j in range(order):
            newc.append(coefficient[j]*(order-j))
        coefficient = newc
    return coefficient


def associated_legendre_value(l, m, theta):
    x = cos(theta)
    if m == 0:
        return polynomial_value(legendre_coefficient(l), x)*sqrt((2*l+1)/(4*pi))
    coefficient = legendre_coefficient(l)
    for i in range(m):
        newc = []
        order = len(coefficient)-1
        for j in range(order):
            newc.append(coefficient[j]*(order-j))
        coefficient = newc
        coefficient = inner_product(coefficient, abs(sin(theta))/sqrt((l+i+1)*(l-i)))
    return polynomial_value(coefficient, x)*sqrt((2*l+1)/(4*pi))
        
        
def polynomial_value(clist, x):
    value = 0
    for i in clist:
        value = value*x+i
    return value
    
    
def SphericalHarmonics(l, m, theta, phi):
    Y = associated_legendre_value(l, m, theta)
    return Y*cos(m*phi), Y*sin(m*phi)

    
if __name__ == '__main__':
    print(SphericalHarmonics(100,90,3*pi/10,pi/6))