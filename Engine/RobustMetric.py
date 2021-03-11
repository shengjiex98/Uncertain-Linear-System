'''
This file provides APIs to compute the robustness metric of
a linear uncertain system, given an unsafe set
'''

from OrderUncertainties import *
from SplitMet import *
from VisualizationReachSet import *

class RobustMetric:
    '''
    Compute the robustness metric of a given linear uncertain
    system and a corresponding unsafe set
    '''

    def __init__(self,A,cells):
        self.A=A # Linear system (without perturbation)
        self.cells=cells

    def distributePerturbation(self,p):
        '''
        Given a perturbation budget `p`, disribute it in the cells
        mentioned in `self.cells`
        cells=[...,(i,j),...], where (i,j) is a cell in `self.A`

        # Algorithm

        - for a cell c_i \in cells, plug p, B[c_i]=p; and call the ordering algorithm
        - obtain the \sigma_{c_i} from the ordering algoritm
        - order the cells based on the obtained sigma_ values from the above step
            - highest sigma, gets the lowest propotion of p, and so.
        '''
        Er={}
        sensitivityMat=OrdUnc(self.A).getSVSentivity(p)
        sm=0
        for i in self.cells:
            sm+=sensitivityMat[i]
        for i in range(len(self.cells)):
            tmp=(p/sm)*sensitivityMat[self.cells[-1-i]]
            #print(self.cells[i],": ",tmp)
            Er[self.cells[i]]=[(100-tmp)/100,(100+tmp)/100]
        #print(Er)

        return Er

    def checkIntersection(star,unsafe):
        '''
        Returns:
            - False, if `star` and `unsafe` doesn't intersect
            - True, otherwise

        Algorithm:
            - Encode `star` and `unsafe` as Gurobi variables
            - Find a point which lies on both `star` and `unsafe`
                - If such a point exists, i.e. the optimization is feasible: return True
                - Otherwise: return False
        '''
        intersectFlag=True

        model = Model("qp")

        # Encode `star` . . .

        C_star=star[0]
        V_star=star[1]
        P_star=star[2]
        n_star=len(C_star)
        vecSize_star=V_star.shape[1]

        ## Encode predicate variables
        predVars_star=[]
        for i in range(vecSize_star):
            name="pred_"+str(i)
            predVars_star.append(model.addVar(-GRB.INFINITY,GRB.INFINITY,name=name,vtype='C'))
            model.addConstr(predVars_star[i]>=P_star[i][0],name+".1")
            model.addConstr(predVars_star[i]<=P_star[i][1],name+".2")

        ## Encode the states
        stateVars_star=[]
        for i in range(n_star):
            #stateVars_star.append(model.addVar(-GRB.INFINITY,GRB.INFINITY,name="State_RS_"+str(i),vtype='C'))
            obj=C_star[i]
            for j in range(vecSize_star):
                obj=obj+(V_star[i][j]*predVars_star[j])
            #model.addConstr(stateVars_star[i]==obj)
            stateVars_star.append(obj)

        # . . . Encode `star` end

        # Encode `unsafe` . . .

        C_unsafe=unsafe[0]
        V_unsafe=unsafe[1]
        P_unsafe=unsafe[2]
        n_unsafe=len(C_unsafe)
        vecSize_unsafe=V_unsafe.shape[1]

        ## Encode predicate variables
        predVars_unsafe=[]
        for i in range(vecSize_unsafe):
            name="predU_"+str(i)
            predVars_unsafe.append(model.addVar(-GRB.INFINITY,GRB.INFINITY,name=name,vtype='C'))
            model.addConstr(predVars_unsafe[i]>=P_unsafe[i][0],name+".1")
            model.addConstr(predVars_unsafe[i]<=P_unsafe[i][1],name+".2")

        ## Encode the states
        stateVars_unsafe=[]
        for i in range(n_unsafe):
            #stateVars_unsafe.append(model.addVar(-GRB.INFINITY,GRB.INFINITY,name="State_Unsafe_"+str(i),vtype='C'))
            obj=C_unsafe[i]
            for j in range(vecSize_unsafe):
                obj=obj+(V_unsafe[i][j]*predVars_unsafe[j])
            #model.addConstr(stateVars_unsafe[i]==obj)
            stateVars_unsafe.append(obj)

        # . . . Encode `unsafe` end

        # Add constraints for finding a point that lies on both `star` and `unsafe`
        for i in range(n_star):
            model.addConstr(stateVars_star[i]==stateVars_unsafe[i])
            #model.addConstr(stateVars_star[i]>=stateVars_unsafe[i])
            #model.addConstr(stateVars_star[i]<=stateVars_unsafe[i])

        # Set Objective
        model.setObjective(stateVars_star[0]) # Can be anything --- we just care about feasibility

        # Check feasibility
        model.optimize()

        print()

        status = model.Status
        if status==GRB.Status.UNBOUNDED:
            print("UNBOUNDED ")
        else:
            if status == GRB.Status.INF_OR_UNBD or \
               status == GRB.Status.INFEASIBLE  or \
               status == GRB.Status.UNBOUNDED:
                intersectFlag=False
            else:
                intersectFlag=True
                '''for v in model.getVars():
                    print('%s %g' % (v.varName, v.x))'''

        return intersectFlag

    def getRobustMetric(self,IS,t,unsafe):
        '''
        Given initial set, unsafe set: return the robustness metric

        # Algorithm:
        - Keep increasing p:
            - \Lambda=self.distributePerturbation(p)
            - if unsafe:
                break
        '''
        p=100
        delta=100

        it=1

        t_taken=time.time()
        while True:
            #print(p)
            Er=self.distributePerturbation(p)
            #Er={}
            #print(Er)
            sp=Split(self.A,Er,IS,t)
            RS=sp.getReachableSetAll()
            fg_intersect=RobustMetric.checkIntersection(RS,unsafe)
            if fg_intersect==True:
                t_taken=time.time()-t_taken
                if True:
                    print("* Number of iterations: ",it)
                    print("* Total Perturbation: ",p-delta)
                    print("* Total Time Taken: ",t_taken)
                Visualization.displayRSnUnsafe(RS,unsafe,th1=0,th2=1,name="RS")
                return ((p-delta),Er)
            #Visualization.displayRSnUnsafe(RS,unsafe,th1=0,th2=1,name="RS")
            p+=delta
            it+=1

        return -404


# Test Case
if False:
    A=np.array([
    [20,3],
    [2,51.2]
    ])
    cells=[(0,1),(1,0)]
    C=[0,0]
    V=np.array([
    [1,0],
    [0,1]
    ])
    P=[(-1,1),(-1,1)]
    rs=(C,V,P)
    t=2
    P_unsafe=[(600,800),(600,800)]
    unsafe=(C,V,P_unsafe)
    rm=RobustMetric(A,cells)
    #rm.distributePerturbation(10)
    metric=rm.getRobustMetric(rs,t,unsafe)
    print(metric)
