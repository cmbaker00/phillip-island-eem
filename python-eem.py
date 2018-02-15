from openpyxl import load_workbook
from string import ascii_lowercase
import itertools
import numpy as np
import numpy.matlib

def read_matrix(fname):
    wb = load_workbook(filename = fname)
    sheets = wb.sheetnames
    ws = wb[sheets[0]]
    n_row = num_rows(ws)
    n_col = num_cols(ws)
    if n_row != n_col:
        raise Exception('Number of rows and columns is not consistant')
    m = np.array(read_all_rows(ws, 2, n_row, 2, n_col))
    return m

def read_vector(fname):
    wb = load_workbook(filename = fname)
    sheets = wb.sheetnames
    ws = wb[sheets[0]]
    n_row = num_rows(ws, st = 1)
    r = np.array(read_all_rows(ws, 1, n_row, 2, 2))
    return r


def read_row(ws, row, st, ed):
    counter = 1
    output = []
    for s in iter_all_strings():
        if counter >= st:
            output.append(ws[s + str(row)].value)
        if counter == ed:
            return output
        counter += 1

def read_all_rows(ws, rows, rowe, st, ed):
    arr = []
    for row in range(rows, rowe + 1):
        c_row = read_row(ws, row, st, ed)
        c_row = none_2_zero(c_row)
        arr.append(c_row)
    return arr

def none_2_zero(arr):
    for i in range(len(arr)):
        if arr[i] == None:
            arr[i] = 0
    return arr

def num_rows(ws, st = 2):
    flag = True
    while flag:
        if ws['A'+str(st)].value == None:
            return st - 1
        else:
            st += 1

def num_cols(ws, st = 2):
    count = 1
    for s in iter_all_strings():
        if count < st:
            count += 1
            continue
        if ws[s + str(1)].value == None:
            return count - 1
        else:
            count += 1

def iter_all_strings():
    size = 1
    while True:
        for s in itertools.product(ascii_lowercase, repeat=size):
            yield "".join((s.upper() for s in s))
        size +=1

def equilibrium_state(A,r):
    iA = np.linalg.inv(A)
    n = -np.matmul(iA,r)
    return n

def calc_jacobian(A,r,n):
    i_matrix = np.eye(len(n))

    J = i_matrix*r + A*np.matlib.repmat(n,1,len(n)) + i_matrix*np.matmul(A,n)

    return J

def calc_stability(A, r, n):
    J = calc_jacobian(A, r, n)
    ev = np.real(np.linalg.eig(J)[0])
    max_eig = np.max(ev)
    if max_eig < 0:
        return True
    else:
        return False

A = read_matrix('Phillip_islands_community.xlsx')
r = read_vector('Phillip_islands_r.xlsx')
n = equilibrium_state(A,r)
J = calc_jacobian(A,r,n)
st = calc_stability(A, r, n)
print(st)