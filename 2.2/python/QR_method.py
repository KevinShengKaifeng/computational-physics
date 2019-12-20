import numpy as np
from numpy import sqrt


def exchange_row_line(mat, i, j):
    rowi = np.array(mat[i])
    mat[i] = mat[j]
    mat[j] = rowi
    columni = np.array(mat[:, i])
    mat[:, i] = mat[:, j]
    mat[:, j] = columni


def Householder_transform(mat):
    size = len(mat)
    transform_matrix = np.eye(size)
    for i in range(size-2):
        all_zero = False
        for j in range(i+1, size):
            if mat[:, i][j] != 0:
                if j != i+1:
                    exchange_row_line(mat, i+1, j)
                    exchange_row_line(transform_matrix, i+1, j)
                break
            if j == size-1:
                all_zero = True
                break
        if all_zero:
            continue
        c = mat[:, i][i+1:]
        sigma = c[0]/abs(c[0])*sqrt(np.inner(c, c))
        e1 = np.zeros(len(c))
        e1[0] = 1
        u = c+sigma*e1
        beta = sigma*(sigma+c[0])
        uni = np.eye(len(c))
        R = uni-np.outer(u, u)/beta
        U = np.eye(size)
        for j in range(len(c)):
            U[i+j+1, i+1:] = R[j]
        mat = np.dot(U, np.dot(mat, U))
        transform_matrix = np.dot(U, transform_matrix)
    return mat, transform_matrix


def Householder_matrix_find_QR(mat):
    min_value = 1e-18
    size = len(mat)
    Q = np.eye(size)
    for i in range(size-1):
        xi = mat[i, i]
        xk = mat[i+1, i]
        if abs(xk) < min_value:
            continue
        c = xi/sqrt(xi**2+xk**2)
        s = xk/sqrt(xi**2+xk**2)
        G = np.eye(size)
        G[i, i], G[i+1, i+1] = c, c
        G[i, i+1], G[i+1, i] = s, -s
        Q = np.dot(G, Q)
        mat = np.dot(G, mat)
    return Q.T, mat


# 仅当特征值均为实数时收敛
def pragmatic_QR_method(mat):
    if type(mat) is list:
        mat = np.array(mat)
    mat_save = mat.copy()
    mat, transform_mat = Householder_transform(mat)
    size = len(mat)
    eigenvalues = []
    min_value = 1e-15
    max_step = 30
    for i in range(1, size):
        step = 0
        while abs(mat[-1, -2]) > min_value:
            s = mat[-1, -1]
            for j in range(len(mat)):
                mat[j, j] -= s
            Q0, R0 = Householder_matrix_find_QR(mat)
            mat = np.dot(R0, Q0)
            for j in range(len(mat)):
                mat[j, j] += s
            step += 1
            if step > max_step:
                print("iteration can't converge!")
                return
        eigenvalues.append(mat[-1, -1])
        mat = mat[:-1, :-1]
    eigenvalues.append(mat[0, 0])
    eigenvalues.sort()
    eigenvectors = find_eigenvectors(mat_save, eigenvalues)
    return eigenvalues, eigenvectors


def find_eigenvectors(mat, eigenvalues):
    size = len(mat)
    eigenvectors = []
    for i in eigenvalues:
        nmat = np.linalg.inv(mat-(i-1e-4)*np.eye(size))
        max_step = 50
        vector = np.ones(size)/sqrt(size)
        last_vector = vector
        for j in range(max_step):
            vector = np.dot(nmat, vector)
            vector /= sqrt(np.inner(vector, vector))*vector[0]/abs(vector[0])
            delta = vector-last_vector
            if sqrt(np.inner(delta, delta)) < 1e-10:
                eigenvectors.append(vector)
                break
            last_vector = vector
    if len(eigenvectors) != size:
        print("can't find all eigenvectors!")
    return eigenvectors


def make_symmetry(mat):
    return (mat+mat.T)/2


if __name__ == "__main__":
    x = make_symmetry(np.random.rand(100, 100))
    n = pragmatic_QR_method(x)
    print(n[0])
    print(len(n[1]))
