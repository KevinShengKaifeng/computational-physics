from spline_interpolation import spline_interpolation


x = [0, 3, 5, 7, 9, 11, 12, 13, 14, 15]
y = [0, 1.2, 1.7, 2.0, 2.1, 2.0, 1.8, 1.2, 1.0, 1.6]
interpolation_func = spline_interpolation(x, y)
with open("output.txt", 'w') as output:
    deltax = 0.1
    xi = x[0]-1
    while xi <= x[-1]+1:
        output.write("%s\n%s\n" % (xi, interpolation_func(xi)))
        xi += deltax
