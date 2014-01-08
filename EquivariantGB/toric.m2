restart
needsPackage "FourTiTwo"
needsPackage "EquivariantGB"

R = buildERing({symbol y},{2},QQ,3)
S = buildERing({symbol x},{1},QQ,3)

F = {x_0^2*x_1}
M = buildEMonomialMap(S,R,F)

egbToric M

