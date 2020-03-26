'''
Author: Bineet Ghosh, under supervision of Dr. Sridhar Duggirala
Email: ghosh.bineet22@gmail.com
- Given a Linear System with uncertainties, compute the reachable set of the
uncertain system using the Splitting Method.
Documentation: Not yet available. (TODO)
'''

import math
import numpy as np
import copy

class Grid:
    '''
    Grids the given inital set into pieces
    '''

    def __init__(self,Theta,st1,st2):

        self.Theta=Theta
        self.state1=st1
        self.state2=st2

    def splitStar(self):
        '''
        Performs gridding just by dividing the predicate.
        -6<=alpha_i<=6 => -6<=alpha_i<=0, 0<=alpha_i<=-6
        '''
        C=self.Theta[0]
        V=self.Theta[1]
        P=self.Theta[2]

        l1=abs(P[self.state1][0]-P[self.state1][1])
        mid1=math.ceil(l1/2)
        l2=abs(P[self.state2][0]-P[self.state2][1])
        mid2=math.ceil(l2/2)

        p1=(P[self.state1][0],P[self.state1][0]+mid1)
        p2=(P[self.state2][0],P[self.state2][0]+mid2)
        p3=(P[self.state1][0]+mid1,P[self.state1][1])
        p4=(P[self.state2][0]+mid2,P[self.state2][1])

        P1=copy.copy(P)
        P2=copy.copy(P)
        P3=copy.copy(P)
        P4=copy.copy(P)

        P1[self.state1]=p1
        P1[self.state2]=p2
        P2[self.state1]=p1
        P2[self.state2]=p4
        P3[self.state1]=p3
        P3[self.state2]=p2
        P4[self.state1]=p3
        P4[self.state2]=p4

        listStars=[(C,V,P1),(C,V,P2),(C,V,P3),(C,V,P4)]

        return listStars

    def printStars(listStars):

        for star in listStars:
            print("---Star---")
            print("Center: ",star[0])
            print("Vector: \n",star[1])
            print("Predicate: ",star[2])

if False:
    C=[0,0,0]
    V=np.array([
    [1,-1,1],
    [-1,1,0],
    [1,0,-1],
    ])
    P=[(-6,5),(-7,9),(-2,9)]
    rs=(C,V,P)
    Grid.printStars(Grid(rs,0,2).splitStar())
