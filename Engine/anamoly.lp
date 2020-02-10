\ Model qp
\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
 
Subject To
 Pert-C(2,_0).1: Pert(2,_0) >= 0.9
 Pert-C(2,_0).2: Pert(2,_0) <= 1.1
 Pert-C(2,_3).1: Pert(2,_3) >= 0.9
 Pert-C(2,_3).2: Pert(2,_3) <= 1.1
 Pert-C(3,_2).1: Pert(3,_2) >= 0.9
 Pert-C(3,_2).2: Pert(3,_2) <= 1.1
 ReachSet-C-0.1: IS0 >= 2.745427655134477
 ReachSet-C-0.2: IS0 <= 4.61139707606489
 ReachSet-C-1.1: IS1 >= -1.483672152727975
 ReachSet-C-1.2: IS1 <= -0.0468428101955415
 ReachSet-C-2.1: IS2 >= 0.8672287172015348
 ReachSet-C-2.2: IS2 <= 16.47744811710425
 ReachSet-C-3.1: IS3 >= -17.32019778368522
 ReachSet-C-3.2: IS3 <= -5.383286579920052
 ReachSet-C-4: IS4 = 1
 ReachSet-C-5: IS5 = 1
Bounds
 Pert(2,_0) free
 Pert(2,_3) free
 Pert(3,_2) free
 IS0 free
 IS1 free
 IS2 free
 IS3 free
 IS4 free
 IS5 free
End
