\ Model qp
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 0.03 IS0 + [ 0.06 Pert(0,_1) * IS0 ] / 2 
Subject To
 Pert-C(0,_1).1: Pert(0,_1) >= 0.9
 Pert-C(0,_1).2: Pert(0,_1) <= 1.1
 ReachSet-C-0.1: IS0 >= 0.108325959739223
 ReachSet-C-0.2: IS0 <= 0.2334886199644169
 ReachSet-C-1.1: IS1 >= 0.0172045439111268
 ReachSet-C-1.2: IS1 <= 0.3246100357925154
Bounds
 Pert(0,_1) free
 IS0 free
 IS1 free
End
