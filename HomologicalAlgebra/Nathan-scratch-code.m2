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

