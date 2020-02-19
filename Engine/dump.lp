\ Model qp
\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
  IS0 + IS1 + IS2 + IS3 + IS4
Subject To
 Predicate-C-0.1: IS0 <= 0
 Predicate-C-1.1: IS1 <= 0
 Predicate-C-2.1: 1.01 IS2 <= 0
 Predicate-C-3.1: 0.01 IS1 + 1.01 IS3 <= 0
 Predicate-C-4.1: 0.032 IS0 + 0.032 IS1 + 1.01 IS4 <= -0.64
Bounds
 IS0 free
 IS1 free
 IS2 free
 IS3 free
 IS4 free
End
