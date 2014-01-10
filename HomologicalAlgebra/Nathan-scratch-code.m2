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

help genericMatrix
M = coker(matrix(R,x^3))
N = coker(matrix(R,x^4))
FM = complete res M
FN = complet res N
coker
C_(-1)

---
---
-- Trying to constuct extensions --


-- try to understand constructing extensions
restart
R = QQ[x]
A = coker matrix({{x^4}})
B = coker matrix({{x^3}})

FA = complete res A
MA = ker FA.dd_0
Hom(MA,B)
myMatrix = cover Hom(MA,B)

cover MA
cover B
R
random(R^1,R^{1:-4})

-- we want:
f = (- random(R^1,R^{1:-2}))
E = coker inducedMap(FA_0 ++ B, MA, inducedMap(FA_0,MA,id_(FA_0)) || f)

-- this is push out

map(E, B, 0*id_(R^1) || id_B)
map(A,E, inducedMap(A, FA_0,id_(FA_0)) | 0*id_(R^1))

chainComplex{map(A,E, inducedMap(A, FA_0,id_(FA_0)) | 0*id_(R^1)), map(E, B, 0*id_(R^1) || id_B)
 }
prune HH oo
-- so above seems to be an extension.



restart
R = QQ[x]
constructRandomExtensions = method()
constructRandomExtensions(ZZ,ZZ,ZZ) := (a,b,d) -> (
    A := coker matrix({{x^a}});
    B := coker matrix({{x^b}});
    f := (- random(R^1,R^{1:-d}));
    FA := complete res A;
    MA := ker FA.dd_0;
    E := coker inducedMap(FA_0 ++ B, MA, inducedMap(FA_0,MA,id_(FA_0)) || f);
    chainComplex{map(A,E, inducedMap(A, FA_0,id_(FA_0)) | 0*id_(R^1)), map(E, B, 0*id_(R^1) || id_B)
    })

Extension = constructRandomExtensions(6,5,3)

prune HH Extension

-- scratch extension code --

constructExtensions = method()



---- COR scratch
restart
needsPackage "SpectralSequences"
Q = QQ[a,b,c]
K = koszul vars Q
K' = koszul matrix {{a^3,b^3,c^3}}
Rmod = coker matrix {{a^3,b^3,c^3}}
R = Q/ideal{a^3,b^3,c^3}
KR = koszul vars R
M = coker matrix {{a,b,0},{0,b,c}}
pfM = pushForward(phi,M)
res pfM
Mres = res(M,LengthLimit=>7)
phi = map(R,Q)
pfMres = chainComplex apply(drop(spots Mres,1), i -> map(pushForward(phi,target Mres.dd_i),
	                   pushForward(phi,source Mres.dd_i),
			   matrix applyTable(entries Mres.dd_i, f -> sub(f,Q))));
filtCOR = (filteredComplex pfMres) ** (Rmod ** K);
filtCOR = (filteredComplex pfMres) ** KR;
ssCOR = prune spectralSequence filtCOR
ssCOR^0
ssCOR2 = ssCOR^2;
support ssCOR^2

help "spots"
keys ssCOR^2
select(spots ssCOR^3 .dd, l -> ssCOR^3 .dd_l != 0) 
X = ssCOR^3 .dd_{5,0}
ssCOR^3 .dd.degree
source X
target X
select(spots ssCOR^5.dd, l -> true) 

---
---

A = QQ[x,y,z]

M = coker vars A

F = res M

N = coker matrix(A,{{x^2,y,z}})

G = res N

H = Hom(F_1,G_0)

ourMap = Hom(F_1, inducedMap(N, G_0, id_(G_0)))

source ourMap == H

L = target ourMap

(ourMap) * (H_0)
(ourMap) * (H_1)
(ourMap) * (H_2)

mPrime = (ourMap) * (H_0) + (ourMap) * (H_1) + (ourMap) * (H_2)

ourMap*(H_0 + H_1 + H_2)

H_0 + H_1 + H_2
-- mPrime is an element of L
-- How to check if mPrime is zero??
ourMap

mPrime 
mPrime == 0_L

H
help"genericMatrix"
vars A

ring L

lre = flatten(entries mm)
lre 
numgens H

lre#1

ourMap*(sum apply(numgens H , i-> (lre#i) * H_i)) 

-- something like this will work.
-- assume given m \in Hom(F_1,G_0) which is a cycle
-- are we assuming that m is a matrix of an appropriate size??
-- and defines an extension 0 --> N --> ? --> M --> 0
isExtTrivial = (m, M, N) -> (
    ourMap := Hom( )
    )


flatten mm
flatten entries mm
# (gens L)

L

l = sum apply(flatten entries mm, i -> (ourMap) * i)
l == 0_L

ker relations M
