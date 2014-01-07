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

--rough sketch of construction 3.6

restart
omega = method()
omega(ZZ,ChainComplex) := (n,C) -> (
     coker C.dd_(n + 1)
     )

--
-- building check for double-dual map
--

restart
kerDD = method()
kerDD(Module) := M -> (
     DM = coker dual presentation M;
     Ext^1(DM, ring M)
     )

cokerDD = method()
cokerDD(Module) := M -> (
     DM = coker dual presentation M;
     Ext^2(DM, ring M)
     )

doubleDualMap = method()
doubleDualMap(Module):= M -> (
--     map(dual (dual M), M, id_(ambient M))
     map(dual (dual M), M, id_(cover M))
     )

kk = ZZ/ ideal (101)
R = kk[x]/ideal (x^2)
M = coker matrix {{x}}


--the following breaks the code for both
-- map(dual (dual M), M, id_(cover M))  and
-- map(dual (dual M), M, id_(ambient M))
R = ZZ
L = ZZ^1
N = coker matrix {{2}}
M = N++L
A = doubleDualMap(M)
K = kerDD(M)
K = prune kerDD(M)
C = cokerDD(M)
C = prune cokerDD(M)
prune ker A == K
prune coker A == C
