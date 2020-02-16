'''
Author: Bineet Ghosh, under supervision of Dr. Sridhar Duggirala
Email: ghosh.bineet22@gmail.com

- Given a Linear System with uncertainties, compute the reachable set of the
uncertain system using the Splitting Method.

Note that this is not the correct version, it uses interval arithmetic.

Documentation: Not yet available. (TODO)
'''

import matplotlib.pyplot as plt
import numpy as np
import math
from gurobipy import *

class Visualization:
    '''
    Given a generaized star, and two states; this class provides APIs to
    visualize the star.
    '''

    def __init__(self,th1,th2,st1,st2):
        self.theta1=th1 # X Axis
        self.theta2=th2 # Y Axis
        self.star1=st1 # Perturbation Free Star
        self.star2=st2 # Perturbed Star
        self.n=st1[1].shape[0]

    def getPlots(self):
        print("Under Construction!!")

        C=self.star1[0]
        V=self.star1[1]
        P=self.star1[2]
        X_list=[]
        Y_list=[]

        C2=self.star2[0]
        V2=self.star2[1]
        P2=self.star2[2]
        X_list2=[]
        Y_list2=[]

        semiDefFlag=False
        model = Model("qp")
        model.setParam( 'OutputFlag', False )

        # Create Predicate Variables
        predVars=[]
        for i in range(self.n):
            name="Pred"+str(i)
            predVars.append(model.addVar(-GRB.INFINITY,GRB.INFINITY,name=name,vtype='C'))
        #-----------------------

        # Axes Variables
        X=model.addVar(-GRB.INFINITY,GRB.INFINITY,name="X",vtype='C')
        Y=model.addVar(-GRB.INFINITY,GRB.INFINITY,name="Y",vtype='C')
        #-------------------------

        # Create the Star Constraints
        objX=0
        for i in range(self.n):
            objX=objX+(predVars[i]*V[self.theta1][i])
        objX=C[self.theta1]+objX

        objY=0
        for i in range(self.n):
            objY=objY+(predVars[i]*V[self.theta2][i])
        objY=C[self.theta2]+objY

        model.addConstr(X==objX,"X Axis")
        model.addConstr(Y==objY,"Y Axis")
        #-----------------------------------

        # Predicate Constraints
        for i in range(self.n):
            a=P[i][0][0]
            b=P[i][0][1]
            if a==b:
                model.addConstr(predVars[i]==a,name)
            else:
                model.addConstr(predVars[i]>=min(a,b),name+".1")
                model.addConstr(predVars[i]<=max(a,b),name+".2")
        #-----------------------------------


        # Quadrant Specific Constraints

        # 1st Quadrant

        obj=X+Y
        model.setObjective(obj,GRB.MAXIMIZE)

        for an in range(90):
            if an==0:
                model.addConstr(X==0,"Angle")
            else:
                m=math.tan(math.radians(an))
                model.addConstr(Y==m*X,"Angle")
            try:
                model.optimize()
                #model.write("dump.bas")
                status = model.Status
                if status==GRB.Status.UNBOUNDED:
                    print("UNBOUNDED ")
                else:
                    if status == GRB.Status.INF_OR_UNBD or \
                       status == GRB.Status.INFEASIBLE  or \
                       status == GRB.Status.UNBOUNDED:
                        print('**The model cannot be solved because it is infeasible or unbounded**')
                    else:
                        xVal=model.getVarByName("X").x
                        yVal=model.getVarByName("Y").x
                        X_list.append(xVal)
                        Y_list.append(yVal)
            except:
                semiDefFlag=True

            if semiDefFlag==True:
                print("Shoot!!")

            semiDefFlag=False
            model.remove(model.getConstrByName("Angle"))
        #-----------------------------

        # 2nd Quadrant

        obj=X-Y
        model.setObjective(obj,GRB.MAXIMIZE)

        for an in range(0,-90,-1):
            if an==-90:
                model.addConstr(X==0,"Angle")
            else:
                m=math.tan(math.radians(an))
                model.addConstr(Y==m*X,"Angle")
            try:
                model.optimize()
                #model.write("dump.bas")
                status = model.Status
                if status==GRB.Status.UNBOUNDED:
                    print("UNBOUNDED ")
                else:
                    if status == GRB.Status.INF_OR_UNBD or \
                       status == GRB.Status.INFEASIBLE  or \
                       status == GRB.Status.UNBOUNDED:
                        print('**The model cannot be solved because it is infeasible or unbounded**')
                    else:
                        xVal=model.getVarByName("X").x
                        yVal=model.getVarByName("Y").x
                        X_list.append(xVal)
                        Y_list.append(yVal)
            except:
                semiDefFlag=True

            if semiDefFlag==True:
                print("Shoot!!")

            semiDefFlag=False
            model.remove(model.getConstrByName("Angle"))
        #-----------------------------

        # 3rd Quadrant

        obj=-X-Y
        model.setObjective(obj,GRB.MAXIMIZE)

        for an in range(-90,-180,-1):
            if an==-90:
                model.addConstr(X==0,"Angle")
            else:
                m=math.tan(math.radians(an))
                model.addConstr(Y==m*X,"Angle")
            try:
                model.optimize()
                #model.write("dump.bas")
                status = model.Status
                if status==GRB.Status.UNBOUNDED:
                    print("UNBOUNDED ")
                else:
                    if status == GRB.Status.INF_OR_UNBD or \
                       status == GRB.Status.INFEASIBLE  or \
                       status == GRB.Status.UNBOUNDED:
                        print('**The model cannot be solved because it is infeasible or unbounded**')
                    else:
                        xVal=model.getVarByName("X").x
                        yVal=model.getVarByName("Y").x
                        X_list.append(xVal)
                        Y_list.append(yVal)
            except:
                semiDefFlag=True

            if semiDefFlag==True:
                print("Shoot!!")

            semiDefFlag=False
            model.remove(model.getConstrByName("Angle"))
        #-----------------------------

        # 4th Quadrant

        obj=-X+Y
        model.setObjective(obj,GRB.MAXIMIZE)

        for an in range(90,180):
            if an==90:
                model.addConstr(X==0,"Angle")
            else:
                m=math.tan(math.radians(an))
                model.addConstr(Y==m*X,"Angle")
            try:
                model.optimize()
                #model.write("dump.bas")
                status = model.Status
                if status==GRB.Status.UNBOUNDED:
                    print("UNBOUNDED ")
                else:
                    if status == GRB.Status.INF_OR_UNBD or \
                       status == GRB.Status.INFEASIBLE  or \
                       status == GRB.Status.UNBOUNDED:
                        print('**The model cannot be solved because it is infeasible or unbounded**')
                    else:
                        xVal=model.getVarByName("X").x
                        yVal=model.getVarByName("Y").x
                        X_list.append(xVal)
                        Y_list.append(yVal)
            except:
                semiDefFlag=True

            if semiDefFlag==True:
                print("Shoot!!")

            semiDefFlag=False
            model.remove(model.getConstrByName("Angle"))
        #-----------------------------


        model=0
        semiDefFlag=False
        model2 = Model("qp")
        model2.setParam( 'OutputFlag', False )

        # Create Predicate Variables
        predVars=[]
        for i in range(self.n):
            name="Pred"+str(i)
            predVars.append(model2.addVar(-GRB.INFINITY,GRB.INFINITY,name=name,vtype='C'))
        #-----------------------

        # Axes Variables
        X=model2.addVar(-GRB.INFINITY,GRB.INFINITY,name="X",vtype='C')
        Y=model2.addVar(-GRB.INFINITY,GRB.INFINITY,name="Y",vtype='C')
        #-------------------------

        # Create the Star Constraints
        objX=0
        for i in range(self.n):
            objX=objX+(predVars[i]*V2[self.theta1][i])
        objX=C2[self.theta1]+objX

        objY=0
        for i in range(self.n):
            objY=objY+(predVars[i]*V2[self.theta2][i])
        objY=C2[self.theta2]+objY

        model2.addConstr(X==objX,"X Axis")
        model2.addConstr(Y==objY,"Y Axis")
        #-----------------------------------

        # Predicate Constraints
        for i in range(self.n):
            a=P2[i][0][0]
            b=P2[i][0][1]
            if a==b:
                model2.addConstr(predVars[i]==a,name)
            else:
                model2.addConstr(predVars[i]>=min(a,b),name+".1")
                model2.addConstr(predVars[i]<=max(a,b),name+".2")
        #-----------------------------------

        obj=X+Y
        model2.setObjective(obj,GRB.MAXIMIZE)

        for an in range(90):
            if an==0:
                model2.addConstr(X==0,"Angle")
            else:
                m=math.tan(math.radians(an))
                model2.addConstr(Y==m*X,"Angle")
            try:
                model2.optimize()
                #model.write("dump.bas")
                status = model2.Status
                if status==GRB.Status.UNBOUNDED:
                    print("UNBOUNDED ")
                else:
                    if status == GRB.Status.INF_OR_UNBD or \
                       status == GRB.Status.INFEASIBLE  or \
                       status == GRB.Status.UNBOUNDED:
                        print('**The model cannot be solved because it is infeasible or unbounded**')
                    else:
                        xVal=model2.getVarByName("X").x
                        yVal=model2.getVarByName("Y").x
                        X_list2.append(xVal)
                        Y_list2.append(yVal)
            except:
                semiDefFlag=True

            if semiDefFlag==True:
                print("Shoot!!")

            semiDefFlag=False
            model2.remove(model2.getConstrByName("Angle"))
        #-----------------------------

        # 2nd Quadrant

        obj=X-Y
        model2.setObjective(obj,GRB.MAXIMIZE)

        for an in range(0,-90,-1):
            if an==-90:
                model2.addConstr(X==0,"Angle")
            else:
                m=math.tan(math.radians(an))
                model2.addConstr(Y==m*X,"Angle")
            try:
                model2.optimize()
                #model.write("dump.bas")
                status = model2.Status
                if status==GRB.Status.UNBOUNDED:
                    print("UNBOUNDED ")
                else:
                    if status == GRB.Status.INF_OR_UNBD or \
                       status == GRB.Status.INFEASIBLE  or \
                       status == GRB.Status.UNBOUNDED:
                        print('**The model cannot be solved because it is infeasible or unbounded**')
                    else:
                        xVal=model2.getVarByName("X").x
                        yVal=model2.getVarByName("Y").x
                        X_list2.append(xVal)
                        Y_list2.append(yVal)
            except:
                semiDefFlag=True

            if semiDefFlag==True:
                print("Shoot!!")

            semiDefFlag=False
            model2.remove(model2.getConstrByName("Angle"))
        #-----------------------------

        # 3rd Quadrant

        obj=-X-Y
        model2.setObjective(obj,GRB.MAXIMIZE)

        for an in range(-90,-180,-1):
            if an==-90:
                model2.addConstr(X==0,"Angle")
            else:
                m=math.tan(math.radians(an))
                model2.addConstr(Y==m*X,"Angle")
            try:
                model2.optimize()
                #model.write("dump.bas")
                status = model2.Status
                if status==GRB.Status.UNBOUNDED:
                    print("UNBOUNDED ")
                else:
                    if status == GRB.Status.INF_OR_UNBD or \
                       status == GRB.Status.INFEASIBLE  or \
                       status == GRB.Status.UNBOUNDED:
                        print('**The model cannot be solved because it is infeasible or unbounded**')
                    else:
                        xVal=model2.getVarByName("X").x
                        yVal=model2.getVarByName("Y").x
                        X_list2.append(xVal)
                        Y_list2.append(yVal)
            except:
                semiDefFlag=True

            if semiDefFlag==True:
                print("Shoot!!")

            semiDefFlag=False
            model2.remove(model2.getConstrByName("Angle"))
        #-----------------------------

        # 4th Quadrant

        obj=-X+Y
        model2.setObjective(obj,GRB.MAXIMIZE)

        for an in range(90,180):
            if an==90:
                model2.addConstr(X==0,"Angle")
            else:
                m=math.tan(math.radians(an))
                model2.addConstr(Y==m*X,"Angle")
            try:
                model2.optimize()
                #model.write("dump.bas")
                status = model2.Status
                if status==GRB.Status.UNBOUNDED:
                    print("UNBOUNDED ")
                else:
                    if status == GRB.Status.INF_OR_UNBD or \
                       status == GRB.Status.INFEASIBLE  or \
                       status == GRB.Status.UNBOUNDED:
                        print('**The model cannot be solved because it is infeasible or unbounded**')
                    else:
                        xVal=model2.getVarByName("X").x
                        yVal=model2.getVarByName("Y").x
                        X_list2.append(xVal)
                        Y_list2.append(yVal)
            except:
                semiDefFlag=True

            if semiDefFlag==True:
                print("Shoot!!")

            semiDefFlag=False
            model2.remove(model2.getConstrByName("Angle"))
        #-----------------------------


        #------------------------------

        return (X_list,Y_list,X_list2,Y_list2)



    def displayPlot(self):
        (X1,Y1,X2,Y2)=self.getPlots()
        plt.axes()
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("State "+str(self.theta1))
        plt.ylabel("State "+str(self.theta2))
        plt.plot(X1,Y1,'bo',label="Unperturbed")
        plt.plot(X2,Y2,'r+',label="Perturbed")
        plt.axis('scaled')
        plt.legend()
        plt.show()


if False:
    C=[0,0]
    V=np.array([
    [1,0],
    [0,1]
    ])
    P=np.array([
    [(-1e-7,1e-7)],
    [(-1e-7,1e-7)]
    ])
    rs=(C,V,P)
    C2=[0,0]
    V2=np.array([
    [0,1],
    [-1,1]
    ])
    P2=np.array([
    [(-1e-6,1e-6)],
    [(-1e-6,1e-6)]
    ])
    rs2=(C2,V2,P2)
    v=Visualization(0,1,rs,rs2)
    v.displayPlot()
