import os,sys
PROJECT_ROOT = os.environ['ULS_ROOT_DIR']
sys.path.append(PROJECT_ROOT)

import pickle
import numpy as np
import math
from scipy.io import loadmat
from scipy.sparse import csr_matrix, csc_matrix


from Parameters import *
from lib.SplitMet import *
from lib.OrderUncertainties import *
from lib.VisualizationComp import *
from lib.StarOperations import *

class AnasthesiaPKPD:

    @staticmethod
    def createMatrix(A,B,mode,h):
        ''' Creates a single matrix based on
        . or +.
        In case of . a rough approximation is
        done'''

        n1=np.size(A,0)
        if (np.size(B)>0):
            n2=np.size(B,1)
        else:
            n2=0
        n=n1+n2

        C=np.zeros((n,n),dtype=np.float)
        if mode=='+':
            for i in range(n1):
                for j in range(n1):
                    C[i][j]=A[i][j]
            for i in range(n1):
                j2=0
                for j in range(n1,n1+n2):
                    C[i][j]=B[i][j2]
                    j2=j2+1
            for i in range(n1,n1+n2):
                C[i][i]=1
        elif mode=='.':
            I=np.zeros((n1,n1),dtype=np.float)
            for i in range(n1):
                I[i][i]=1
            A2=h*A
            A2=np.add(I,A2)
            B2=h*B
            for i in range(n1):
                for j in range(n1):
                    C[i][j]=A2[i][j]
            for i in range(n1):
                j2=0
                for j in range(n1,n1+n2):
                    C[i][j]=B2[i][j2]
                    j2=j2+1
            for i in range(n1,n1+n2):
                C[i][i]=1

        return C

    def getDynamics(weight=25):
        #With weight of 25 Kg
        v1=458.4*weight
        k10=0.1527*pow(weight,-0.3)
        k12=0.114
        k13=0.0419
        k21=0.055
        k31=0.0033
        kd=40 #with Td=20s
        A=np.array([
        [-(k10+k12+k13),k12,k13,0],
        [k21,-k21,0,0],
        [k31,0,-k31,0],
        [kd,0,0,-kd]
        ])
        B=np.array([
        [1/v1],
        [0],
        [0],
        [0]
        ])

        return (A,B)

    def getReachSetC1C2():
        '''
        Returns the reachable set for Comparment 1 and 2, i.e., C1 and C2
        '''
        C=[2,2,2,2,2]
        V=np.array([
        [1,0,0,0,0],
        [0,1,0,0,0],
        [0,0,1,0,0],
        [0,0,0,1,0],
        [0,0,0,0,1],
        ])
        #P=[(0,6),(0,10),(0,10),(1,8),(2,2)]
        #P=[(2,4),(4,6),(4,6),(3,5),(1,3)]
        P=[(-1,1),(-1,1),(-1,1),(-1,1),(-1,1)]
        initialSet=(C,V,P)

        p=1.8
        Er={
        (1,0): [1-(p/100),1+(p/100)],
        (1,1): [1-(p/100),1+(p/100)],
        (0,2): [1-(p/100),1+(p/100)],
        (2,2): [1-(p/100),1+(p/100)]
        }
        T=3

        (dynA,dynB)=AnasthesiaPKPD.getDynamics()
        A=AnasthesiaPKPD.createMatrix(dynA,dynB,'.',0.01)
        #np.set_printoptions(precision=3)

        print(">> STATUS: Computing Reachable Sets . . .")
        time_taken=time.time()
        rs=Split(A,Er,initialSet,T)
        (reachORS,reachRS)=rs.getReachableSetAllList()
        time_taken=time.time()-time_taken
        print("\tTime Taken: ",time_taken)
        print(">> STATUS: Reachable Sets Computed!")
        #exit(0)
        rs_cora_tool = loadmat(PROJECT_ROOT+'/cora_comp/pkpd_rs_cora.mat')
        print((csr_matrix(rs_cora_tool['Z_c']).toarray().reshape(-1)))
        #exit()
        coraZonoGSize=csr_matrix(rs_cora_tool['Z_G']).toarray().shape[1]

        Visualize.vizRSComp([reachORS[-1]],[(csr_matrix(rs_cora_tool['Z_c']).toarray().reshape(-1),csr_matrix(rs_cora_tool['Z_G']).toarray(),[(-1,1)]*coraZonoGSize)],0,1,fname="pkpd_1-2_cora_comp")








if True:
    AnasthesiaPKPD.getReachSetC1C2()
