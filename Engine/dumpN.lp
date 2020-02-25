\ Model qp
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - Pred1 + [ 2 Pert(0,_1) * Pred1 ] / 2 
Subject To
 Pert-C(0,_1).1: Pert(0,_1) >= 0.98
 Pert-C(0,_1).2: Pert(0,_1) <= 1.02
 Pred-C-0.1: Pred0 <= 5
 Pred-C-0.2: Pred0 >= -5
 Pred-C-1.1: Pred1 <= 5
 Pred-C-1.2: Pred1 >= -5
Bounds
 Pert(0,_1) free
 Pred0 free
 Pred1 free
End
