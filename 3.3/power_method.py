from numpy import ones, array, zeros, eye
from chase_solve import chase_solve


def power_method_solve(A, x=[]):
    if not x:
        x = ones(len(A))
    else:
        x /= max(x)
    eigenvalue = 0
    step = 20
    accuracy = 0.001
    for i in range(step):
        x = A.dot(x)
        if abs(max(x)-eigenvalue) < accuracy:
            eigenvalue = max(x)
            x /= eigenvalue
            break
        eigenvalue = max(x)
        x /= eigenvalue
    return eigenvalue, x


def tridiagonal_dot(a, b, x):
    n = len(a)
    nx = a.copy()
    for j in range(n):
        if j == 0:
            nx[j] = x[j]*a[j]+x[j+1]*b[j+1]
        elif j == n-1:
            nx[j] = x[j-1]*b[j]+x[j]*a[j]
        else:
            nx[j] = x[j-1]*b[j]+x[j]*a[j]+x[j+1]*b[j+1]
    return nx


def tridiagonal_power_solve(a, b, x=[]):
    n = len(a)
    if not x:
        x = ones(n)
    else:
        x /= max(x)
    eigenvalue = 0
    step = 20
    accuracy = 0.001
    for i in range(step):
        x = tridiagonal_dot(a, b, x)
        if abs(max(x)-eigenvalue) < accuracy:
            eigenvalue = max(x)
            x = x/eigenvalue
            break
        eigenvalue = max(x)
        x = x/eigenvalue
    return 1/eigenvalue, x


def tridiagonal_inverse_power_solve(a, b, x=[]):
    n = len(a)
    if not x:
        x = ones(n)
    else:
        x /= max(x)
    eigenvalue = 0
    step = 20
    accuracy = 0.001
    A = eye(n)
    for i in range(n):
        A[i][i] = a[i]
        if i == 0:
            continue
        A[i-1][i] = b[i]
        A[i][i-1] = b[i]
    for i in range(step):
        x = chase_solve(b.copy(), a.copy(), b.copy(), x)
        if abs(max(x)-eigenvalue) < accuracy:
            eigenvalue = max(x)
            x /= eigenvalue
            break
        eigenvalue = max(x)
        x /= eigenvalue
    return eigenvalue, x


def test1():
    A = array([[3,1],[1,3]])
    return power_method_solve(A)


if __name__ == "__main__":
    print(test1())
