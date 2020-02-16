\ Model qp
\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
  X + Y
Subject To
 X_Axis: - Pred0 + X = 0
 Y_Axis: - Pred1 + Y = 0
 Pred1.1: Pred0 >= -5
 Pred1.2: Pred0 <= 5
 Pred1.1: Pred1 >= -5
 Pred1.2: Pred1 <= 5
 Angle: - X + Y = 0
Bounds
 Pred0 free
 Pred1 free
 X free
 Y free
End
