import numpy as np
from datetime import datetime
from enum import Enum
from math import floor, log

# minimum value of a matrix element
MIN = 1
# maximum value of a matrix element
MAX = 20
# matrix size
SIZE = 128
# previous size potency of 2
POTENCY_SIZE = 1 if SIZE == 1 else (lambda n: 2 ** int(floor(log(SIZE - 1, 2))))(SIZE)

class MatrixMethod(Enum):
    STANDARD = 1
    STRASSEN = 2

def print_matrix(matrix):
    for i in matrix:
        print('\t'.join(map(str, i)))

def gen_random_matrix(size_a, size_b):
    return np.random.randint(MIN, MAX, (size_a, size_b)).tolist()

def standard_method(m1, m2):
    res = [[0 for x in range(len(m1))] for y in range(len(m2))]

    for i in range(len(m1)):
        for j in range(len(m2[0])):
            for k in range(len(m2)):
                res[i][j] += m1[i][k] * m2[k][j]

    return res

def split(matrix):
    first = matrix
    second = matrix
    third = matrix
    forth = matrix

    while(len(first) > len(matrix)/2):
        first = first[:len(first)/2]
        second = second[:len(second)/2]
        third = third[len(third)/2:]
        forth = forth[len(forth)/2:]

    while(len(first[0]) > len(matrix[0])/2):
        for i in range(len(first[0])/2):
            first[i] = first[i][:len(first[i])/2]
            second[i] = second[i][len(second[i])/2:]
            third[i] = third[i][:len(third[i])/2]
            forth[i] = forth[i][len(forth[i])/2:]

    return first, second, third, forth
            
def add_matrix(a, b):
    if type(a) == int:
        d = a + b
    else:
        d = []
        for i in range(len(a)):
            c = []
            for j in range(len(a[0])):
                c.append(a[i][j] + b[i][j])
            d.append(c)
    return d

def sub_matrix(a, b):
    if type(a) == int:
        d = a - b
    else:
        d = []
        for i in range(len(a)):
            c = []
            for j in range(len(a[0])):
                c.append(a[i][j] - b[i][j])
            d.append(c)
    return d

def strassen(N, M, size):
    if size == 1:
        XY = [[0]]
        XY[0][0] = N[0][0] * M[0][0]
        return XY
    elif size == POTENCY_SIZE:
        return standard_method(N, M)
    else:
        A, B, C, D = split(N)
        E, F, G, H = split(M)

        # P1 = A(F - H) | AF - AH
        P1 = strassen(A, sub_matrix(F, H), size/2)

        # P2 = (A + B)H
        P2 = strassen(add_matrix(A, B), H, size/2)

        # P3 = (C + D)E
        P3 = strassen(add_matrix(C, D), E, size/2)

        # P4 = D(G - E)
        P4 = strassen(D, sub_matrix(G, E), size/2)

        # P5 = (A + D)(E + H)
        P5 = strassen(add_matrix(A, D), add_matrix(E, H), size/2)

        # P6 = (B - D)(G + H)
        P6 = strassen(sub_matrix(B, D), add_matrix(G, H), size/2)

        # P7 = (A - C)(E + F)
        P7 = strassen(sub_matrix(A, C), add_matrix(E, F), size/2)

        # XY_1 = P5 + P4 - P2 + P6
        XY_1 = add_matrix(add_matrix(P5, sub_matrix(P4, P2)), P6)

        # XY_2 = P1 + P2
        XY_2 = add_matrix(P1, P2)

        # XY_3 = P3 + P4
        XY_3 = add_matrix(P3, P4)

        # XY_4 = P1 + P5 - P3 - P7
        XY_4 = sub_matrix(add_matrix(P1, sub_matrix(P5, P3)), P7)

        XY = [[0 for row in range(len(XY_1) * 2)] for col in range(len(XY_1) * 2)]

        for i in range(len(XY_1)):
            for j in range(len(XY_1)):
                XY[i][j] = XY_1[i][j]
                XY[i][j+len(XY_1)] = XY_2[i][j]
                XY[i+len(XY_1)][j] = XY_3[i][j]
                XY[i+len(XY_1)][j+len(XY_1)] = XY_4[i][j]
    
        return XY

def strassen_method(A, B):
    return strassen(A, B, len(A))

def multiply_matrix(m1, m2, method):
    if method == MatrixMethod.STANDARD:
        return standard_method(m1, m2)
    elif method == MatrixMethod.STRASSEN:
        return strassen_method(m1, m2)
    else:
        print("O metodo [{}] nao eh permitido ou nao esta configurado".format(method))

    return None

if __name__ == '__main__':
    m1 = gen_random_matrix(SIZE, SIZE)
    m2 = gen_random_matrix(SIZE, SIZE)
    
    full_strassen = []
    full_standard = []

    for i in range(30):
        now = datetime.now()
        multiply_matrix(m1, m2, MatrixMethod.STRASSEN)
        full_strassen.append((datetime.now() - now).total_seconds())

        now = datetime.now()
        multiply_matrix(m1, m2, MatrixMethod.STANDARD)
        full_standard.append((datetime.now() - now).total_seconds())

    average_strassen = sum(full_strassen) / len(full_strassen)
    average_standard = sum(full_standard) / len(full_standard)

    print("Strassen levou: {0:.10f} segundos".format(average_strassen))
    print("Standard levou: {0:.10f} segundos".format(average_standard))
    