import os,sys
PROJECT_ROOT = os.environ['ULS_ROOT_DIR']

sys.path.append(PROJECT_ROOT)

import numpy as np
import control
from sympy import symbols, Matrix
from lib.OrderUncertainties import *
from Benchmarks import sys_variables

class Bench:
    def __init__(self, A, B, C=None, D=None) -> None:
        self.A = A
        self.B = B
        self.nx = A.shape[1]
        self.nu = B.shape[1]
        self.C = C if C else np.eye

def get_sensitivity(A, B):
    mat = OrdUnc(A)
    if mat.determineCase() == 1:
        unc = mat.distinctPos(B)
    else:
        unc = mat.multSig(B)
    return unc

def all_sensitivity(sys, K):
    # A_bar_i = A + BKE_i where E_i = I(n) with 
    for i in range(sys.nx):
        E_i = Matrix(np.eye(sys.nx))
        E_i[i, i] = symbols('l')

def find_AB(sys, x):
    nx = sys.A.shape[0]
    l = symbols('l')
    E = np.eye(nx)
    E[x,x] = l
    print(nx)

def main():
    find_AB("")

main()
