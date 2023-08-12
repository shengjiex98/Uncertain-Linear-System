import os,sys
PROJECT_ROOT = os.environ['ULS_ROOT_DIR']
sys.path.append(PROJECT_ROOT)

from lib.OrderUncertainties import *

def get_sensitivity(A, B):
    mat = OrdUnc(A)
    if mat.determineCase() == 1:
        unc = mat.distinctPos(B)
    else:
        unc = mat.multSig(B)
    return unc

def all_sensitivity(A):
    pass
