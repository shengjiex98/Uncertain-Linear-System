'''
Author: Bineet Ghosh
Email: ghosh.bineet22@gmail.com
- Generates random half-spaces to be used as a template for
computation of U
'''
import os,sys
PROJECT_ROOT = os.environ['ULS_ROOT_DIR']
sys.path.append(PROJECT_ROOT)
from Parameters import *

import numpy as np
import random
import math



class PredGen:
    '''
    Generates random half-spaces to be used as a template for
    computation of U
    '''

    def __init__(self,A,n,samples):
        self.A=A # Generate Random matrix close to this matrix
        self.n=n # Create random matrix of dimension n
        self.nSamples=samples # Number of random predicates

    def getRandPred(self):
        '''
        Returns random samples of the uncertain linear
        system.
        '''
        sam=[]

        for i in range(self.nSamples-1):
            sam.append(self.getRandMat())
        sam.append(self.A)

        return sam

    def getRandMat(self):
        '''
        Returns a random matrix representing
        predicates
        '''

        rMat=np.zeros((self.n,self.n),dtype=float)

        for i in range(self.n):
            for j in range(self.n):
                a=self.A[i][j]*(1-CLOSE)
                b=self.A[i][j]*(1+CLOSE)
                lb=min(a,b)
                ub=max(a,b)
                if lb==0 and ub==0:
                    lb=-EPSILON
                    ub=EPSILON
                    rd=random.randrange(math.ceil(lb*PRECISION),math.floor(ub*PRECISION))
                    r=float(rd/PRECISION)
                    #r=0
                else:
                    rd=random.randrange(math.ceil(lb*PRECISION),math.floor(ub*PRECISION))
                    r=float(rd/PRECISION)
                    if (r>ub or r<lb):
                        print("Gotcha!!")
                        print("r :",r)
                        print("[lb,ub]: ",lb,ub)
                        exit(0)
                rMat[i][j]=r

        return rMat
