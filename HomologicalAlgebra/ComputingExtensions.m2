restart
-- here are some experimental scripts related to computting extensions.
-- assume for the moment that F ( --> A ) and G ( --> ) B are chain complexes with homology in
-- degree zero only and H_0 F = A and H_0 G = B  

needsPackage"SpectralSequences" --  might need the Hom code in SS pacakge

constructRandomExtensions = method()
constructRandomExtensions(ZZ,ChainComplex,ChainComplex) := (d,F,G) -> (
    B := coker G.dd_1;
    z := Hom(ker F.dd_0,B);
    S := z.ring;
    numCycleGens := numgens z ;
    randomRingElements := flatten entries random(S^1,S^{numCycleGens:-d});
    randomElementOfZ := sum apply(numCycleGens, i -> randomRingElements#i * z_{i});
    E := coker map (F_0 ++ B, ker F.dd_0, id_(F_0) || - (super randomElementOfZ)) ;
    e := chainComplex(map(coker F.dd_1,E, id_(F_0) | 0*id_(G_0)), inducedMap(E,B, 0*id_(F_0) || id_(G_0)))
)

-- try example
S = QQ[a,b,c]
A = coker matrix{{a^3,b^4,c^2}}
B = coker matrix{{a^1,b^6,c^3} }
F = complete res A
G = complete res B
coker G.dd_1
E = constructRandomExtensions(3,F,G)
prune HH E

-- try example over a PID

R = QQ[x]
A = coker matrix{{x^5}}
B = coker matrix{{x^4}}
F = complete res A
G = complete res B
E = constructRandomExtensions(3,F,G)
prune HH E

-- so constructRandomExtensions seems to work more generally.
restart
R = QQ[z,y]
F = res ideal(z^2,z*y)
G = res ideal(z*y,y^2)
constructRandomExtensions(3,F,G)

constructExtensions = method()
-- assume for the moment that F ( --> A ) and G ( --> ) B are chain complexes with homology in
-- degree zero only and H_0 F = A and H_0 G = B  
-- want a matrix M of the appropriate size m x n where m = rank F_1, n = rank G_0
-- for the moment we assume that M is an element of Z below
-- might be better to do this with a list of ring elements,
-- assume that the 
constructExtensions(Matrix,ChainComplex,ChainComplex) := (M,F,G) -> (   
    B := coker G.dd_1;
    Z := Hom(image F.dd_1,coker G.dd_1);
    R := Z.ring;
    E := coker map (F_0 ++ B, ker F.dd_0, id_(F_0) || - M) ;
    e := chainComplex(map(A,E, id_(F_0) | 0*id_(G_0)), inducedMap(E,B, 0*id_(F_0) || id_(G_0)))
    )

-- try some examples 
R = QQ[x]
A = coker matrix{{x^5,x^3,x^2}}
B = coker matrix{{x^4,x}}
F = complete res A
G = complete res B
M = matrix{{x^4 + x^3}}
Z = Hom(ker F.dd_0,coker G.dd_1)
E = constructExtensions(matrix{{x^2 + x^3}}, F, G)
prune HH E

-- try a more general example
S = QQ[a,b,c]
A = coker matrix{{a^3,b^4,c^2}}
B = coker matrix{{a^1,b^6,c^3}}
F = complete res A
G = complete res B
Z = Hom(ker F.dd_0,coker G.dd_1)
gens Z
(a^2+ b^3 ) * Z_{0} -- this seems to be an element of Z
E = constructExtensions((a^2+ b^3 ) * Z_{0},F,G)
prune HH E


-- now want to decide if a given extension is zero
-- input element of  Hom(ker F.dd_0,coker G.dd_1)
isExtensionZero = method()
isExtensionZero(List,ChainComplex,ChainComplex) := (L,F,G) -> (
    Z := Hom(image F.dd_1,coker G.dd_1);
    B := coker G.dd_1;
    A := coker F.dd_1;
    M := image F.dd_1;
    j := inducedMap(F_0, M ,id_(F_0));
    myMod := coker Hom(j, B);
    ppi := inducedMap(myMod, Z, id_(ambient Z)) ; 
    (ppi * (sum apply(numgens Z, i -> (L#i) * Z_{i}))) == 0
)

-- try some examples

R = QQ[x]
A = coker matrix{{x^5}}
B = coker matrix{{x^4}}
F = complete res A
G = complete res B
Z = Hom(image F.dd_1,coker G.dd_1);
M = image F.dd_1

isExtensionZero({x^4 + x^3},F,G)
isExtensionZero({x^4 + x^3}, F, G)
isExtensionZero({x^4}, F, G)
isExtensionZero({0_R}, F, G)
isExtensionZero({x^6}, F, G)
isExtensionZero({x^1}, F, G)
isExtensionZero({x^2}, F, G)
isExtensionZero({x^3}, F, G)
isExtensionZero({x^1}, F, G)
isExtensionZero({0}, F, G)



-- so maybe this is OK.

-- over a PID we could compute the type of an extension.  Do this later.

-- An alternative way to do the above is to use the Hom complex and
-- start with elements of Hom(F_1,G_0)