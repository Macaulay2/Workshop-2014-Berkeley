restart
needsPackage "EquivariantGB"

R = buildERing({symbol y},{2},QQ,2)
S = buildERing({symbol x},{1},QQ,2)

F = {x_0*x_1}
F = {x_0^2*x_1}
M = buildEMonomialMap(S,R,F)

G = egbToric M






