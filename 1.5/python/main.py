import math


def ax(vx, vy):
    global B
    return -B*vx*(vx**2+vy**2)**(1/2)


def ay(vx, vy):
    global B
    g = 9.8
    return -B*vy*(vx**2+vy**2)**(1/2)-g


def shoot1(output, theta0=30/180*math.pi, *, deltat=0.05, v0=700):
    t, x, y, vx, vy = 0, 0, 0, v0*math.cos(theta0), v0*math.sin(theta0)
    Nmax = 4000
    while y >= 0:
        output.write("%s\n%s\n%s\n%s\n%s\n%s\n%s\n" % (x, y, vx, vy, ax(vx, vy), ay(vx, vy), t))
        x += vx * deltat
        y += vy * deltat
        vx += ax(vx, vy) * deltat
        vy += ay(vx, vy) * deltat
        t += deltat
        if t > Nmax*deltat:
            print("calculation end after %s steps, t=%ss" % (Nmax, t))
            return


def shoot2(output, theta0=30/180*math.pi, *, deltat=0.05, v0=700):
    t, x, y, vx, vy = 0, 0, 0, v0*math.cos(theta0), v0*math.sin(theta0)
    Nmax = 4000
    while y >= 0:
        output.write("%s\n%s\n%s\n%s\n%s\n%s\n%s\n" % (x, y, vx, vy, ax(vx, vy), ay(vx, vy), t))
        x += (vx+ax(vx, vy) * deltat/2) * deltat
        y += (vy+ay(vx, vy) * deltat/2) * deltat
        vx += (ax(vx, vy)+ax(vx+ax(vx, vy)*deltat, vy+ay(vx, vy)*deltat))/2 * deltat
        vy += (ay(vx, vy)+ay(vx+ax(vx, vy)*deltat, vy+ay(vx, vy)*deltat))/2 * deltat
        t += deltat
        if t > Nmax*deltat:
            print("calculation end after %s steps, t=%ss" % (Nmax, t))
            return


shoot = shoot2
with open("output30.txt", 'w') as output:
    B = 4*10**-5
    shoot(output)
with open("output30_no_resistance.txt", 'w') as output:
    B = 0
    shoot(output)
with open("output40.txt", 'w') as output:
    B = 4*10**-5
    shoot(output, 40/180*math.pi)
with open("output40_no_resistance.txt", 'w') as output:
    B = 0
    shoot(output, 40/180*math.pi)
with open("output50.txt", 'w') as output:
    B = 4 * 10 ** -5
    shoot(output, 50/180*math.pi)
with open("output50_no_resistance.txt", 'w') as output:
    B = 0
    shoot(output, 50/180*math.pi)
