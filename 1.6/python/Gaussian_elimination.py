import time


def find_max(column, A):
    n = len(A)
    rcount = column
    maxrow = column
    while rcount < n:
        if abs(A[rcount][column]) > abs(A[maxrow][column]):
            maxrow = rcount
        rcount += 1
    return maxrow


def exchange(l1, l2, A, b):
    r1, b1 = A[l1], b[l1]
    A[l1], b[l1] = A[l2], b[l2]
    A[l2], b[l2] = r1, b1


def lineadd(target, source, multiple, A, b):
    n = len(A)
    i = source
    while i < n:
        A[target][i] += A[source][i]*multiple
        i += 1
    b[target] += b[source] * multiple


def gelimination_solve(A, b):
    n = len(A)
    ccount = 0
    # 选主元消元法将系数矩阵化为上三角矩阵
    while ccount < n-1:
        exchange(find_max(ccount, A), ccount, A, b)
        if A[ccount][ccount] == 0:
            raise "strange coefficient matrix"
        i = ccount+1
        while i < n:
            lineadd(i, ccount, -A[i][ccount]/A[ccount][ccount], A, b)
            i += 1
        ccount += 1
    if A[n-1][n-1] == 0:
        raise "strange coefficient matrix"
    # 回代
    ccount = n-1
    while ccount > 0:
        i = ccount-1
        while i >= 0:
            b[i] -= b[ccount]*A[i][ccount]/A[ccount][ccount]
            A[i][ccount] = 0
            i -= 1
        ccount -= 1
    # 求解，对角矩阵归一化
    ccount = 0
    while ccount < n:
        b[ccount] /= A[ccount][ccount]
        ccount += 1
    return b


def load_line(linestr):
    line = linestr.split()
    i = 0
    while i < len(line):
        line[i] = int(line[i])
        i += 1
    return line


def load_matrix(datafile):
    A, b = [], []
    with open(datafile, 'r') as data:
        n = int(data.readline())
        i = 0
        while i < n:
            A.append(load_line(data.readline()))
            i += 1
        b = load_line(data.readline())
    return A, b


if __name__ == '__main__':
    start = time.clock()
    solution = gelimination_solve(*load_matrix("Gaussian_elimination_data.txt"))
    end = time.clock()
    print("solve time: ", end - start, 's'+"\nsolution:")
    for xi in solution:
        print(xi)
