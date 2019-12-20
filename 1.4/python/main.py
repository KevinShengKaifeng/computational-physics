import differentiate_method, recursion_method, expansion_method
from math import pi


def SphericalHarmonics(l, m, theta, phi):
    assert theta > 0 and theta < pi, 'theta error!'
    if m > l-10:
        f = differentiate_method.SphericalHarmonics
    elif m < 6:
        f = recursion_method.SphericalHarmonics
    elif theta < (pi/500) or theta > (499*pi/500):
        f = expansion_method.SphericalHarmonics
    elif m > 0.7*l:
        f = differentiate_method.SphericalHarmonics
    else:
        f = recursion_method.SphericalHarmonics
    re, im = f(l, m, theta, phi)
    if re!=0 and im!=0:
        if abs(re/im) < 1e-10:
            re = 0
        elif abs(im/re) < 1e-10:
            im = 0
    return re, im
        
        
thetal = [pi/1000, 3*pi/10, 501*pi/1000]
ll = [100, 500, 1000]
with open("output.txt", 'w') as output:
    for theta in thetal:
        for l in ll:
    #        output.write("l:%s, m:%s, theta:%spi, phi:pi/6" % (l,1,theta/pi))
            output.write(str(SphericalHarmonics(l,1,theta,pi/6))+'\n')
    #        output.write("l:%s, m:%s, theta:%spi, phi:pi/6" % (l,l/100,theta/pi))
            output.write(str(SphericalHarmonics(l,l//100,theta,pi/6))+'\n')
    #        output.write("l:%s, m:%s, theta:%spi, phi:pi/6" % (l,l/10,theta/pi))
            output.write(str(SphericalHarmonics(l,l//10,theta,pi/6))+'\n')
    #        output.write("l:%s, m:%s, theta:%spi, phi:pi/6" % (l,l-1,theta/pi))
            output.write(str(SphericalHarmonics(l,l-1,theta,pi/6))+'\n')
print("finished")
