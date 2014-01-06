restart

A = QQ[x,y,z]

-- example 1

M = A^1

N = dual dual M

isIsomorphism(map(N,M,id_(A^1)))

-- example 2

M = coker vars A

ambient M

N = Hom(Hom(M,A^1), A^1)

ambient N

ambient M


map(Hom(Hom(M,A^1), A^1), M, id_(A^1) )

dual M
dual dual M


--source of non-free reflexives:  high syzygy
-- of residue field of non-regular Gor ring.

k = ZZ / 101  

S = k[x] / ideal(x^2)

-- this is our non-regular Gor ring.

m = coker vars S
res m

res(m, LengthLimit => 10)

ambient dual dual m

isIsomorphism map(dual dual m, m, id_(ambient m))

-- returns true because m is supposed
-- to be reflexive (because is maximal CM over 
-- Gor ring)

ambient m
dual dual m

-- example of extend

help extend

C = res(m, LengthLimit => 10)

extend(C,C, id_(C_0))


-- make some trivial methods trying to understand
-- 3.6 of the paper

omega = method()

omega(ZZ,ChainComplex) := (n,C) -> (
    coker C.dd_(n + 1)
    )

C
omega(1,C)

omega(2,C)

needsPackage"SpectralSequences"

truncate(C,-1)
truncate(C,1)

g = 2

G = omega(g,C)

L = res(dual G, LengthLimit => 10)

truncate(C,-8)


sigmaGminus1 = (dual truncate(C,-(length C - g)))[g-1]

dual G

omega(1-g, dual C)

ambient dual G

ambient omega(1-g, dual C)

-- the canonical map we want to lift ?!
map(dual G, omega(1-g, dual C), id_(
	ambient omega(1-g, dual C))
)
