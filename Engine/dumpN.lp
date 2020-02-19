\ Model qp
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 0.032 IS0 - 0.032 IS1 + [ 0.064 Pert(4,_0) * IS0
   + 0.064 Pert(4,_1) * IS1 ] / 2 
Subject To
 Pert-C(4,_0).1: Pert(4,_0) >= -9
 Pert-C(4,_0).2: Pert(4,_0) <= 11
 Pert-C(4,_1).1: Pert(4,_1) >= -9
 Pert-C(4,_1).2: Pert(4,_1) <= 11
 ReachSet-C-0.1: IS0 <= 1
 ReachSet-C-0.1: IS0 >= -1
 ReachSet-C-1.1: IS1 <= 1
 ReachSet-C-1.1: IS1 >= 1
 ReachSet-C-2.1: IS2 <= 1
 ReachSet-C-2.1: IS2 >= 1
 ReachSet-C-3.1: IS3 <= 1
 ReachSet-C-3.1: IS3 >= 1
 ReachSet-C-4.1: IS4 <= 1
 ReachSet-C-4.1: IS4 >= 1
Bounds
 Pert(4,_0) free
 Pert(4,_1) free
 IS0 free
 IS1 free
 IS2 free
 IS3 free
 IS4 free
End
