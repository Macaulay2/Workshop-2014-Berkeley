A = QQ[x,y,z]
M = A^1
rawDual M
presentation M

--source of non-free reflexives: high syzygy of residue field of non-regular Gor rings.

kk = ZZ/ ideal (101);
R = kk[x]/ideal(x^2)
S = oo
K = coker map(S^1, , {gens S})
resolution(K, LengthLimit =>6)
C = oo
