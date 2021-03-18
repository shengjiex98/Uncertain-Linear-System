'''
This file provides APIs to drive the Robust Metric
computation.
'''
from RobustMetric import *
from Benchmarks import *

mode='.'
h=0.01

class DriverRM:
    '''
    APIs to compute the robustness metric for various Benchmarks
    '''
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

    def dcConv():
        dynA=Benchmarks.DCConv.A
        dynB=Benchmarks.DCConv.B
        A=DriverRM.createMatrix(dynA,dynB,mode,h)
        cells=[(0,2),(1,1)]
        #cells=[(0,0),(0,1),(0,2),(1,0),(1,1),(1,2),(2,0),(2,1),(2,2)]
        C=[0,0,0]
        V=np.array([
        [1,0,0],
        [0,1,0],
        [0,0,1]
        ])
        P=[(-1,1),(-1,1),(1,1)]
        rs=(C,V,P)
        T=4
        th1=0
        th2=2
        name="DcConv"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        P_unsafe=[(2,3),(-1,1),(1,1)]
        unsafe=(C,V,P_unsafe)
        rm=RobustMetric(A,cells)
        #metric=rm.getRobustMetric(rs,T,unsafe)
        rm.compareHeus(100,rs,T)
        #print(metric)

    def holesCXc():
        dynA=Benchmarks.HolesCXc.A
        dynB=Benchmarks.HolesCXc.B
        A=DriverRM.createMatrix(dynA,dynB,mode,h)
        cells=[(0,3),(1,2),(3,2),(4,3),(5,5),(9,9)]
        cells=[]
        for i in range(10):
            for j in range(10):
                cells.append((i,j))
        C=[0,0,0,0,0,0,0,0,0,0]
        V=np.array([
        [1,0,0,0,0,0,0,0,0,0],
        [0,1,0,0,0,0,0,0,0,0],
        [0,0,1,0,0,0,0,0,0,0],
        [0,0,0,1,0,0,0,0,0,0],
        [0,0,0,0,1,0,0,0,0,0],
        [0,0,0,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,1,0,0,0],
        [0,0,0,0,0,0,0,1,0,0],
        [0,0,0,0,0,0,0,0,1,0],
        [0,0,0,0,0,0,0,0,0,1],
        ])
        P=[(-1,1),(-1,1),(1,1),(-1,1),(-1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]
        rs=(C,V,P)
        T=10
        th1=0
        th2=1
        name="HolesCXc"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        P_unsafe=[(-1,1),(-1,1),(1,1),(-1,1),(-1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]
        unsafe=(C,V,P_unsafe)
        rm=RobustMetric(A,cells)
        #metric=rm.getRobustMetric(rs,T,unsafe)
        #print(metric)
        rm.compareHeus(1000000,rs,T)
        #sp2.printReachableSetAll(th1,th2,name+"_Org")

    def acc():
        dynA=Benchmarks.ACC.A
        dynB=Benchmarks.ACC.B
        A=DriverRM.createMatrix(dynA,dynB,mode,h)
        P=10
        cells=[(0,3),(2,3)]
        C=[0,0,0,0]
        V=np.array([
        [1,0,0,0],
        [0,1,0,0],
        [0,0,1,0],
        [0,0,0,1],
        ])
        P=[(0,35),(5,50),(0,35),(1,1)]
        rs=(C,V,P)
        T=50
        th1=0
        th2=1
        name="ACC"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        P_unsafe=[(40,45),(5,50),(0,35),(1,1)]
        unsafe=(C,V,P_unsafe)
        rm=RobustMetric(A,cells)
        #metric=rm.getRobustMetric(rs,T,unsafe)
        rm.compareHeus(20000,rs,T)
        #print(metric)

    def laneChange():
        dynA=Benchmarks.LaneChange.A
        dynB=Benchmarks.LaneChange.B
        A=DriverRM.createMatrix(dynA,dynB,mode,h)
        cells=[(3,6),(4,6),(5,6)]
        C=[0,0,0,0,0,0,0]
        V=np.array([
        [1,0,0,0,0,0,0],
        [0,1,0,0,0,0,0],
        [0,0,1,0,0,0,0],
        [0,0,0,1,0,0,0],
        [0,0,0,0,1,0,0],
        [0,0,0,0,0,1,0],
        [0,0,0,0,0,0,1],
        ])
        P=[(0,50),(3,3.5),(0,0),(20,30),(0,0),(0,0),(1,1)]
        rs=(C,V,P)
        T=50
        th1=0
        th2=1
        start=0
        n=3
        step=0.01
        name="LaneChange"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        P_unsafe=[(52,60),(3,3.5),(0,0),(20,30),(0,0),(0,0),(1,1)]
        unsafe=(C,V,P_unsafe)
        rm=RobustMetric(A,cells)
        #metric=rm.getRobustMetric(rs,T,unsafe)
        rm.compareHeus(20,rs,T)
        print(metric)

    def pkpd2():
        dynA=Benchmarks.PKPD2.A
        dynB=Benchmarks.PKPD2.B
        A=DriverRM.createMatrix(dynA,dynB,mode,h)
        cells=[(0,4),(0,0)]
        C=[5,5,5,5,5]
        V=np.array([
        [1,0,0,0,0],
        [0,1,0,0,0],
        [0,0,1,0,0],
        [0,0,0,1,0],
        [0,0,0,0,1],
        ])
        P=[(1,6),(0,10),(0,10),(1,8),(0,200)]
        rs=(C,V,P)
        T=30
        name="PKPD2"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        P_unsafe=[(7,9),(0,10),(0,10),(1,8),(0,200)]
        unsafe=(C,V,P_unsafe)
        rm=RobustMetric(A,cells)
        #metric=rm.getRobustMetric(rs,T,unsafe)
        rm.compareHeus(20,rs,T)
        print(metric)
        #sp2=Split(A,E,rs,T)
        #sp2.printReachableSetAll(th1,th2,name+"_Org")

    def giradI():
        dynA=Benchmarks.GiradI.A
        dynB=Benchmarks.GiradI.B
        A=DriverRM.createMatrix(dynA,dynB,mode,h)
        cells=[(0,0),(1,0)]
        C=[0,0]
        V=np.array([
        [1,0],
        [0,1]
        ])
        P=[(0.9,1.1),(-0.1,0.1)]
        rs=(C,V,P)
        T=10
        th1=0
        th2=1
        start=0
        n=200
        step=0.01
        name="GiradI-Full"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        P_unsafe=[(1.2,1.7),(-0.1,0.1)]
        unsafe=(C,V,P_unsafe)
        rm=RobustMetric(A,cells)
        #metric=rm.getRobustMetric(rs,T,unsafe)
        rm.compareHeus(20,rs,T)
        print(metric)



DriverRM.pkpd2()
