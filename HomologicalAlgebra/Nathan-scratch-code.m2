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

isIsomorphism map(dual dual m, m, id_(cover m))

-- returns true because m is supposed
-- to be reflexive (because is maximal CM over 
-- Gor ring)

ambient m
dual dual m

-- example of extend
help extend
C = res(m, LengthLimit => 10)
extend(C,C, id_(C_0))

-- we can truncate the complex
-- this is the notation P < g in the paper.
-- truncate a chain complex
-- appears in spectral sequence package

needsPackage"SpectralSequences"

truncate(C,-1)
truncate(C,1)

restart
needsPackage"SpectralSequences"
-- load this so we can truncate
k = ZZ / 101  
S = k[x] / ideal(x^2)
kk = coker vars S
g = 2
C = res(kk, LengthLimit => 10)

-- make some trivial methods trying to understand
-- 3.6 of the paper

omega = method()

omega(ZZ,ChainComplex) := (n,C) -> (
    coker C.dd_(n + 1)
    )

C
omega(1,C)
omega(2,C)
G = omega(g,C)
L = res(dual G, LengthLimit => 10)
truncate(C,-8)
sigmaGminus1 = (dual truncate(C,-(length C - g)))[g-1]
dual G
omega(1-g, dual C)
ambient dual G
ambient omega(1-g, dual C)

-- trying to compute the canonical map
-- in 3.6
-- the canonical map we want to lift ?!
map(dual G, omega(1-g, dual C), id_(
	cover omega(1-g, dual C)))

ambient (image dual C.dd_g)
ambient dual G

dual C.dd_g
-- or is the canonical map one of the maps
-- a or b below??
a = map(dual G, omega(1-g, dual C), dual (C.dd_g))

b = inducedMap(dual G, omega(1-g, dual C), dual (C.dd_g))

dual (C.dd_g)

omega(1-g, dual C)

dual b

source a
source b
target a
target b

(dual C).dd_(-2)
coker (dual C).dd_(-1)

Cdual = dual C

Cdual_(-3)

image Cdual.dd_(-2)
coker Cdual.dd_(-1)

-- question:
-- why is one of these the zero map and one
-- of them not the zero map??
inducedMap(image Cdual.dd_(-2), 
    coker Cdual.dd_(-1), Cdual.dd_(-2)) 
map(image Cdual.dd_(-2), coker Cdual.dd_(-1), Cdual.dd_(-2)) 

gens coker Cdual.dd_(-1)

gens image Cdual.dd_(-2)

