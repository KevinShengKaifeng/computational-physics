import time
from CGsolve import CGsolve, vector_add


def generator(size):
    if size%2:
        print("size error!")
        raise ValueError("size error!")
    A = []
    b = []
    for i in range(size):
        A.append({size-1-i: 1/2})
        b.append(1.5)
    for i in range(size-1):
        A[i][i] = 3
        A[i][i+1] = -1
        A[i+1][i] = -1
    A[size-1][size-1] = 3
    b[0], b[size//2-1], b[size//2], b[-1] = 2.5, 1, 1, 2.5
    A = vector_add(A, empty_matrix(size))
    with open("data.txt", 'w') as output:
        output.write(str(size)+'\n')
        for i in range(size):
            for j in range(size):
                if j in A[i]:
                    output.write(str(A[i][j])+' ')
                else:
                    output.write(str(0)+' ')
            output.write('\n')
        for i in range(size):
            output.write(str(b[i])+' ')
    return A, b


def empty_matrix(size):
    a = []
    for i in range(size):
        ai = []
        for j in range(size):
            ai.append(0)
        a.append(ai)
    return a


start = time.clock()
solution = CGsolve(*generator(100), 1e-6)
end = time.clock()
print("solve time:", end - start, 's'+"\nsolution:")
with open("output.txt", 'w') as output:
    for i in solution:
        output.write(str(i)+'\n')
