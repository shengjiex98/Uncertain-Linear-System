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

class CD:

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

    def getDynamics():
        dynamics = loadmat(DATA_PATH+'/large_benchmarks/CDplayer.mat')
        A = dynamics['A'].toarray()
        B = csr_matrix(dynamics['B']).toarray()


        return (A,B)

        return (A,B)

    def getReachSet():

        C=[0]*122
        V=np.identity(122)
        P=[(-1,1)]*122
        initialSet=(C,V,P)

        p=1.8
        Er={
        (61,58): [1-(p/100),1+(p/100)],
        (62,57): [1-(p/100),1+(p/100)],
        (63,56): [1-(p/100),1+(p/100)],
        (64,55): [1-(p/100),1+(p/100)],
        (65,54): [1-(p/100),1+(p/100)],
        (66,53): [1-(p/100),1+(p/100)]
        }
        T=150
        initialSet=(C,V,P)


        (dynA,dynB)=CD.getDynamics()
        A=CD.createMatrix(dynA,dynB,'.',0.01)
        #np.set_printoptions(precision=3)

        print(">> STATUS: Computing Reachable Sets . . .")
        time_taken=time.time()
        rs=Split(A,Er,initialSet,T)
        #(reachORS,reachRS)=rs.getReachableSetAllList()
        (reachORS,reachRS)=rs.getReachableSetRed()
        time_taken=time.time()-time_taken
        print("\tTime Taken: ",time_taken)
        print(">> STATUS: Reachable Sets Computed!")
        #exit(0)
        rs_cora_tool = loadmat(PROJECT_ROOT+'/cora_comp/cd_rs_cora.mat')
        print((csr_matrix(rs_cora_tool['Z_c']).toarray().reshape(-1)))
        #exit()
        coraZonoGSize=csr_matrix(rs_cora_tool['Z_G']).toarray().shape[1]

        Visualize.vizRSComp([reachORS[-1]],[(csr_matrix(rs_cora_tool['Z_c']).toarray().reshape(-1),csr_matrix(rs_cora_tool['Z_G']).toarray(),[(-1,1)]*coraZonoGSize)],0,1,fname="pkpd_1-2_cora_comp")
        #Visualize.vizRSComp([reachORS[-1]],[],0,1,fname="pkpd_1-2_cora_comp")








if True:
    CD.getReachSet()
