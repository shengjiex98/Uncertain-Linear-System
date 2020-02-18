\ Model qp
\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
  IS0 + IS1 + IS2 + IS3 + IS4 + IS5
Subject To
 Predicate-C-0.1: IS0 + 0.01 IS2 <= 0
 Predicate-C-1.1: IS1 + 0.01 IS3 <= 0
 Predicate-C-2.1: IS2 + 2e-05 IS4 <= 0
 Predicate-C-3.1: IS3 + 2e-05 IS5 <= 0
 Predicate-C-4.1: IS4 <= 0
 Predicate-C-5.1: IS5 <= 0
Bounds
 IS0 free
 IS1 free
 IS2 free
 IS3 free
 IS4 free
 IS5 free
End
