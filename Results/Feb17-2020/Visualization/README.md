## Description

Given a `<C,V,P>` star representation as _Reachable Set_ and two state variables, this tool plots the reachable set in 2-D space. It uses 360 optimization formulation for 360 angles.

* The blue circles are used to represent the un-perturbed set
* The red plus is used to represent the perturbed set.

----------------



## Example 1

C=[0,0]
V=np.array([
[1,0],
[0,1]
])
P=[(-1e-7,1e-7), (-1e-7,1e-7)]
rs=(C,V,P)
C2=[0,0]
V2=np.array([
[0,1],
[-1,1]
])
P2=[(-1e-6,1e-6), (-1e-6,1e-6)]
rs2=(C2,V2,P2)
v=Visualization(0,1,rs,rs2)
v.displayPlot()



## Example 2

C=[0,0]
V=np.array([
[0,-1],
[-1,0]
])
P=[(-1,1), (-2,2)]
rs=(C,V,P)
C2=[0,0]
V2=np.array([
[0,1],
[1,-1]
])
P2=[(-4,4), (-7,7)]
rs2=(C2,V2,P2)
v=Visualization(0,1,rs,rs2)
v.displayPlot()