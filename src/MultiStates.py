import os,sys
PROJECT_ROOT = os.environ['ULS_ROOT_DIR']
sys.path.append(PROJECT_ROOT)

from lib.OrderUncertainties import *
import numpy as np
import control as ct
import sympy
from sympy import symbols
from Benchmarks import sys_variables

def get_sensitivity(A, B):
    mat = OrdUnc(A)
    if mat.determineCase() == 1:
        unc = mat.distinctPos(B)
    else:
        unc = mat.multSig(B)
    return unc

def all_sensitivity(A):
    pass

def find_AB(sys, x):
    nx = sys.A.shape[0]
    l = symbols('l')
    E = np.eye(nx)
    E[x,x] = l
    print(nx)

def main():
    find_AB("")

main()
