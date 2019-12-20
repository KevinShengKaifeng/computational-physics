def ISRsolve(A, b):
    # Improved Square Root solve, A should be symmetric positive definite
    # solve lij and di, store in A
    n = len(A)
    for i in range(1, n):
        t = []
        for j in range(i):
            t.append(A[i][j])
            if j >= 1:
                for k in range(j):
                    t[j] -= t[k]*A[j][k]
            A[i][j] = t[j]/A[j][j]
        for k in range(i):
            A[i][i] -= t[k]*A[i][k]
    # solve xi, store in b
    for i in range(1, n):
        for k in range(i):
            b[i] -= A[i][k]*b[k]
    b[n-1] /= A[n-1][n-1]
    for i in range(n-1)[::-1]:
        b[i] /= A[i][i]
        for k in range(i+1, n):
            b[i] -= A[k][i]*b[k]
    return b
