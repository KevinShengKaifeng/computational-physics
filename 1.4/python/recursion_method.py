from math import sqrt, cos, sin, pi


def legendre_function(l, x):
    assert l >= 0, 'l value error!'
    if l == 0:
        return 1
    elif l == 1:
        return x
    else:
        p_values = [1, x]
        for i in range(1, l):
            last_p = p_values[1]
            p_values[1] = ((2*i+1)/(i+1)*x*p_values[1]-i/(i+1)*p_values[0])
            p_values[0] = last_p
        return p_values[1]


def normalized_associated_legendre_function(l, m, theta):
    # 为了控制数值的大小，将球谐函数系数加入连带勒让德函数计算
    assert m >= 0 and m <= l, 'm value error!'
    x = cos(theta)
    if m == 0:
        return legendre_function(l, x)*sqrt((2*l+1)/(4*pi))
    elif x == 1 or x == -1:
        return 0
    else:
        values = []
        for i in range(2*m+1):
            values.append(legendre_function(l-m+i, x))
        m0 = 0
        while m0 < m:
            m0 += 1
            values = maddone(values, l-m+m0, m0, theta)
            for i in range(len(values)):
                values[i] /= sqrt((l+m0)*(l-m0+1))
#            print(values)
        return values[0]*sqrt((2*l+1)/(4*pi))
           
           
def maddone(values, firstl, m0, theta):
    values0 = []
    for i in range(len(values)-2):
        values0.append(((firstl+m0+i)/(2*(firstl+i)+1)*(firstl+m0+i-1)*values[i]-(firstl-m0+2+i)/(2*(firstl+i)+1)*(firstl-m0+1+i)*values[i+2])/sin(theta))
    return values0
    
    
def SphericalHarmonics(l, m, theta, phi):
    Y = normalized_associated_legendre_function(l, m, theta)
    return Y*cos(m*phi), Y*sin(m*phi)
    

if __name__ == '__main__':
    print(SphericalHarmonics(100,1,pi/1000,pi/6))
