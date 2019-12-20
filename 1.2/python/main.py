def decay(deltat, output):
    global tao, N0
    n, t = N0, 0
    while t <= 5:
        output.write("%s\n%s\n" % (t, n))
        t += deltat
        n += -n / tao * deltat


tao, N0 = 1.0, 100
with open("output0.4.txt", 'w') as output:
    decay(0.4, output)
with open("output0.2.txt", 'w') as output:
    decay(0.2, output)
with open("output0.1.txt", 'w') as output:
    decay(0.1, output)
with open("output0.05.txt", 'w') as output:
    decay(0.05, output)
