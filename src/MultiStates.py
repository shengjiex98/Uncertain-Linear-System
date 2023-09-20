import os,sys
PROJECT_ROOT = os.environ['ULS_ROOT_DIR']

sys.path.append(PROJECT_ROOT)

import numpy as np
import control as ctrl
from sympy import symbols, Matrix
from lib.OrderUncertainties import *
from Benchmarks import Bench, sys_variables

def get_sensitivity(A, B):
    mat = OrdUnc(A)
    if mat.determineCase() == 1:
        unc = mat.distinctPos(B)
    else:
        unc = mat.multSig(B)
    return unc

def all_sensitivity(bench, K):
    # A_bar_i = A + BKE_i where E_i = I(n) with the (i, i) entry swapped to lambda
    for i in range(bench.nx):
        E = Matrix(np.eye(bench.nx))
        E[i, i] = symbols('l')
        A = bench.A + bench.B * K * bench.E

def find_AB(sys, x):
    nx = sys.A.shape[0]
    l = symbols('l')
    E = np.eye(nx)
    E[x,x] = l
    print(nx)

def main():
    all_sensitivity(sys_variables['F1'], 1)

main()
