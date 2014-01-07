restart
needsPackage "EquivariantGB"
needsPackage "FourTiTwo"

R = buildERing({symbol y},{2},QQ,3)
S = buildERing({symbol x},{1},QQ,3)

F = {x_0^2*x_1}
M = buildEMonomialMap(S,R,F)

M(y_(0,1)*y_(1,0))
A = transpose matrix (((flatten entries M(vars R)) / exponents) / flatten)
toricGroebner(A, R)

