import time
from ISRsolve import ISRsolve
from CGsolve import CGsolve


def generator(size):
    assert size >= 4 and size % 2 == 0, 'size error!'
    A = [[5], [4, 6], [1, 4, 6]]
    for i in range(3, size):
        line = []
        for j in range(i-2):
            line.append(0)
        line += [1, 4, 6]
        A.append(line)
    A[size-1][size-1] = 5
    b = [60, ]
    for i in range(size-2):
        b.append(120)
    b.append(60)
    return A, b


def generator2(size):
    A = [{0: 5}]
    for i in range(size-2):
        A.append({i+1: 6})
    A.append({size-1: 5})
    for i in range(1, size):
        A[i][i-1] = 4
        A[i-1][i] = 4
    for i in range(2, size):
        A[i][i-2] = 1
        A[i-2][i] = 1
    b = [60, ]
    for i in range(size - 2):
        b.append(120)
    b.append(60)
    return A, b


start = time.clock()
solution = ISRsolve(*generator(100))
end = time.clock()
print("solve time:", end - start, 's'+"\nsolution:")
with open("output.txt", 'w') as output:
    for i in solution:
        output.write(str(i)+'\n')