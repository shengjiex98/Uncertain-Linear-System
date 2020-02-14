\ Model qp
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - IS1 + [ 2 Pert(0,_1) * IS1 ] / 2 
Subject To
 Pert-C(0,_1).1: Pert(0,_1) >= 0.9
 Pert-C(0,_1).2: Pert(0,_1) <= 1.1
 ReachSet-C-0.1: IS0 <= 2
 ReachSet-C-0.1: IS0 >= 1
 ReachSet-C-1.1: IS1 <= 2
 ReachSet-C-1.1: IS1 >= 1
Bounds
 Pert(0,_1) free
 IS0 free
 IS1 free
End
