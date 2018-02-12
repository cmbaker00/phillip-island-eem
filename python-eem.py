from openpyxl import load_workbook
from string import ascii_lowercase
import itertools
import numpy as np

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
    for row in range(rows, rowe):
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

print(read_matrix('Phillip_islands_community.xlsx'))