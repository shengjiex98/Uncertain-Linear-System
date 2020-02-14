## Description

Computes U based on the splitting based method. _Generalized Star_ representation is used for the computation. 

Given a `<C,V,P>` star representation as _Over-approximated Reachable Set_. U is returned as `<A_c C, A_cV, P'>` , where `A_c` is the center (point) matrix and `P'` is computed using Gurobi _s.t_ it over-approximates the effect of perturbation given `<C,V,P>`

