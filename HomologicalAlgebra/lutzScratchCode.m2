A = QQ[x,y,z]
M = A^1
rawDual M
presentation M

--source of non-free reflexives: high syzygy of residue field of non-regular Gor rings.

kk = ZZ/ ideal (101)
R = kk[x]/ideal(x^2)
S = oo
K = coker map(S^1, , {gens S})
resolution(K, LengthLimit =>6)
C = oo
K1 = dual K
K2 = dual K1
presentation K1
presentation K2


-- testing extend
restart
kk = ZZ/ ideal (101)
R = kk[x]
M = coker matrix {{x}}
C = resolution( M, LengthLimit =>6)
C.dd
A = matrix{{x}}
extend(C,C,A)
