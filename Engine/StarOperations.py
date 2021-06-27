import numpy as np
import numpy.linalg as LA
from gurobipy import *

class StarOp:
    '''
    This class provides APIs for performing various operations on
    generalized stars (representing a polytope)
    '''

    def bloatFact(star,times):
        '''
        bloat `star` by `times`

        Algorithm:
         - Increase the length of ALL the predicates on both the sides by:
            (((n-th root of `times`)*(length of the predicate))-(length of the predicate))/2
         - The Anchor and Basis Vector remains the same.
        '''

        C=star[0] # Anchor of the star
        V=star[1] # Basis vectors of the star
        P=star[2] # predicate of the star
        n=len(C) # dimension of the system

        n_root=(times**(1/n))
        P_bloat=[]
        for pred in P:
            l_p=abs(pred[1]-pred[0])
            incFact=((n_root*l_p)-l_p)/2
            #print(incFact)
            p_min=min(pred[0],pred[1])-incFact
            p_max=max(pred[0],pred[1])+incFact
            P_bloat.append((p_min,p_max))

        bloatedStar=(C,V,P_bloat)

        #print(P_bloat)

        #VisualizeRS.vizBloating(star,bloatedStar)

        return bloatedStar

    def volume(star):
        '''
        Computes the so-called volume of `star`

        Algorithm:
            - detV = Compute  abs(det(V)), V: set of basis vectors of the star
            - vol=detV * (length of each predicates)
        '''

        C=star[0] # Anchor of the star
        V=star[1] # Basis vectors of the star
        P=star[2] # predicate of the star
        n=len(C) # dimension of the system

        detV=abs(LA.det(V))

        vol=1
        for pred in P:
            l_p=abs(pred[1]-pred[0]) # length of each predicate
            vol=vol*l_p
        vol=vol*detV

        return vol

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
            obj=C_star[i]
            for j in range(vecSize_star):
                obj=obj+(V_star[i][j]*predVars_star[j])
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
            obj=C_unsafe[i]
            for j in range(vecSize_unsafe):
                obj=obj+(V_unsafe[i][j]*predVars_unsafe[j])
            stateVars_unsafe.append(obj)

        # . . . Encode `unsafe` end

        # Add constraints for finding a point that lies on both `star` and `unsafe`
        for i in range(n_star-1):
            model.addConstr(stateVars_star[i]==stateVars_unsafe[i])
            #model.addConstr(stateVars_star[i]>=stateVars_unsafe[i])
            #model.addConstr(stateVars_star[i]<=stateVars_unsafe[i])

        # Set Objective
        model.setObjective(stateVars_star[0]) # Can be anything --- we just care about feasibility

        # Check feasibility
        model.optimize()

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

        return intersectFlag

    def checkIntersectionPoints(star,unsafePoints):
        '''
        Returns:
            - False, if `star` and no points in `unsafePoints` intersect
            - True, otherwise

        Algorithm:
            - Encode `star` and as Gurobi variables
            - Find a point which lies on both `star` and `unsafe`[i]
                - If such a point exists, i.e. the optimization is feasible: return True
                - Otherwise: return False
        '''
        for unsafePoint in unsafePoints:
            X_unsafe=unsafePoint[0]
            Y_unsafe=unsafePoint[1]
            rs=StarOp.checkIntersectionPoint(star,(X_unsafe,Y_unsafe))
            if rs==True:
                return True

        return False

    def checkIntersectionPoint(star,unsafePoint):
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

        X_unsafe=unsafePoint[0]
        Y_unsafe=unsafePoint[1]

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
            obj=C_star[i]
            for j in range(vecSize_star):
                obj=obj+(V_star[i][j]*predVars_star[j])
            stateVars_star.append(obj)

        # . . . Encode `star` end

        # Add constraints for finding a point that lies on both `star` and `unsafe`
        model.addConstr(stateVars_star[0]==X_unsafe)
        model.addConstr(stateVars_star[1]==Y_unsafe)

        # Set Objective
        model.setObjective(stateVars_star[0]) # Can be anything --- we just care about feasibility

        # Check feasibility
        model.optimize()

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

        return intersectFlag


if False:
    C=[0,0]
    V=np.array([
    [0,1],
    [1,1]
    ])
    #P=[(6,11),(5,15),(5,15),(6,13)]
    P=[(-1,1),(-1,1)]
    rs=(C,V,P)
    C_u=[3,3]
    V_u=np.array([
    [0,1,0],
    [1,1,2]
    ])
    P_u=[(-1,1),(-1,1),[-2,1]]
    unsafe=(C_u,V_u,P_u)
    VisualizeRS.vizBloating(rs,unsafe)
    print("\n\nIntersection: ",StarOp.checkIntersection(rs,unsafe))
