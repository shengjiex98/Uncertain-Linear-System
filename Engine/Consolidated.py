'''
Author: Bineet Ghosh, under supervision of Dr. Sridhar Duggirala
Email: ghosh.bineet22@gmail.com
- Given a Linear System with uncertainties, compute the reachable set of the
uncertain system using the Splitting Method.
Documentation: Not yet available. (TODO)
'''

from SplitMet import *
from BloatAPI import *

class SplitBloat:
    '''
    APIs for both Splitting and Bloating
    '''

    def __init__(self,A,E,theta,T,n):
        self.A=A # Linear system (without perturbation)
        self.Er=E
        '''
        Following dictionary data-structure has been used to represent
        a error matrix:
        {
            (i,j): [a,b]
        }
        Indicating array[i][j] has a perturbation of -(1-a)% to (b-1)%
        For eg: [0.9,1,1] means a perturbation of -10% and +10%
        '''
        self.Theta=theta # The initial set
        self.T=T #Max number of steps
        self.methodName="SplitBloat"
        self.name=n

    def getStats(self,th1,th2,start,n,step,methodList,p='slow'):
        sp=Split(self.A,self.Er,self.Theta,self.T)
        sp.printReachableSetAll(th1,th2,self.name)
        self.plotTimeCompare(start,n,step,methodList,p)

    def matrixify(self):
        '''
        Converts the dictionary of error with relative value
        to another dictionary of error with absolute value
        '''
        E_new={}
        for key in self.Er:
            lb=-(1-self.Er[key][0])*self.A[key[0]][key[1]]
            ub=(1-self.Er[key][1])*self.A[key[0]][key[1]]
            E_new[key]=[lb,ub]

        return E_new


    def plotTimeCompare(self,start,n,step,methodList,p='slow'):
        '''
        Plots the bloating factor with time, upto time tBound.
        Comparing all the techniques mentioned in methodList
        in the same plot.
        '''

        plotX=[]
        plotY=[]
        i=1
        time_taken_kagstrom1=0
        time_taken_kagstrom2=0
        time_taken_loan=0
        t_kagstrom1=0
        t_kagstrom2=0
        t_loan=0
        matE=self.matrixify()
        for m in methodList:
            if m.lower()=='kagstrom1':
                t_kagstrom1=time.time()
                bloat=BloatKagstrom(self.A,matE)
                (plotX,plotY)=bloat.computeBloatingFactorWithTime(start,n,step,p)
                time_taken_kagstrom1=time_taken_kagstrom1+(time.time()-t_kagstrom1)
            elif m.lower()=='kagstrom2':
                t_kagstrom2=time.time()
                bloat=BloatKagstrom(self.A,matE)
                (plotX,plotY)=bloat.computeBloatingFactor2WithTime(start,n,step,p)
                time_taken_kagstrom2=time_taken_kagstrom2+(time.time()-t_kagstrom2)
            elif m.lower()=='loan':
                t_loan=time.time()
                bloat=BloatLoan(self.A,matE)
                (plotX,plotY)=bloat.computeBloatingFactorWithTime(start,n,step,p)
                time_taken_loan=time_taken_loan+(time.time()-t_loan)
            plotY=list(filter(lambda x: math.inf!=x,plotY))
            plotX=plotX[:len(plotY)]
            plt.plot(plotX,plotY,label=m)
            i=i+1

        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("Time")
        plt.ylabel("Bloating Factor")
        plt.legend(loc='best')
        plt.savefig("BloatPlot/"+self.name+"Bloat")
        #plt.show()
        plt.close()
        print("-")
        print("Kagstrom 1: ",time_taken_kagstrom1)
        print("Kagstrom 2: ",time_taken_kagstrom2)
        print("Loan: ",time_taken_loan)
