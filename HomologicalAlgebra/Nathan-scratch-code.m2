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
cover m
dual m
cover dual dual m

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

ambient kk


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


restart
R = ZZ

matrix(R,{{2}})

M = R^1 ++ coker matrix(R,{{2}})
cover M

presentation dual M

ambient dual dual M

Hom(Hom(M,R^1),R^1)

N = ZZ^2
dual N
dual dual N


-- to do:  try to 
-- understand Franks biduality code.

-- is there an easy way to
-- augment a resolution
-- to a chain complex which includes the 
-- module we are resolving??
A = QQ[x,y,z]
I = coker vars A
C = res I
C_0
inducedMap(I, C_0, id_(C_0))

-- trying to understand pshFwd method;
-- then want to push forward a chain complex

kk = QQ
R = kk[a,b]
S = kk[z,t]


--- trying to understand Ext^1(M,N) where M and N are modules over a PID
-- step 1:
-- resolve M and N we get
-- 0 -->F_1 --> F_0 --> M --> 0
-- 0 -->G_1 --> G_0 --> N --> 0

-- step 2:
-- want to convert elements \psi \in Hom(F_0,G_1) to extensions M_\psi

-- step 3:  determine when \psi and \psi' determine equivalent extensions
-- enough to check if \psi - \psi' map to zero under the canonical map
-- Hom(F_0,G_1) --> Im(Hom(F_0,G_0) --> Hom(F_0,G_1)) where this map 
-- comes from the Hom complex.

R = ZZ/11[x]

makeMandN := (m,n) -> {coker(matrix(R,{{x^m}})), coker(matrix(R,{{x^n}}))}
makeMandN(3,4)
makeMandN(11,10)

L = makeMandN(3,4)
needsPackage"SpectralSequences"
L#1
L#0

C = Hom(complete res L#0, complete res L#1)

C_1 .cache.components

C_1 .cache.indices

image C.dd_0

HH_0 C
prune HH_(-1) C

C
image C.dd_0
prune C_(-1) / (image C.dd_0)


HH C
prune HH C
Ext^1(L#0,L#1)
prune Ext(L#0,L#1)

prune Ext^0(L#0,L#1)
C

help Ext
