--natural map from a module to its double dual courtesy of Frank Moore
bidualityMap = method()
bidualityMap(Module) := M -> (
   R := ring M;
   Md := Hom(M,R^1);
   Mdd := Hom(Md,R^1);
   W := coker relations Mdd;
   gensMdd := map(W,source gens Mdd,gens Mdd);
   tempMap := matrix apply(numgens M, i -> 
       {map(R^1,Md,matrix {apply(numgens Md, j ->
            (homomorphism Md_{j})*(M_{i}))},Degree=>(degrees M)_i)});
   retVal := map(Mdd,M,(transpose matrix tempMap) // gensMdd);
   retVal
   )

-- find the nth syzygy of a chain complex; the cokernel of the (n+1)'st
-- differential
omega = method()
omega(ZZ,ChainComplex) := (n,C) -> (
     coker C.dd_(n + 1)
)

-- code by Kat to truncate a chain complex; all components of the complex
-- in degrees >= g will be made to 0.
truncateComplex = method()
truncateComplex(ZZ,ChainComplex) := (g,P) -> (
     if g > max P then return P
--     if g <= min P then return image(0*id_P)
     else (
	  K:=new ChainComplex;
	  K.ring=P.ring;
	  for i from min P+1 to max P do (
	       if i < g then K.dd_i=P.dd_i;
	       );
	  return K
	  )
)

--take a resolution of a module, append the module in degree -1
augmentResolution = method(TypicalValue => ChainComplex);
augmentResolution (Module,ChainComplex) := (M,C) -> (
    if not (prune M)== (prune (coker C.dd_1)) then error "expected module
    to be the homology in degree 0";
    if not min C==0 then error "expected complex concentrated 
    in non-negative degrees";
    mapsList = ();
    mapsList = append(mapsList, map(M, C_0,id_(C_0)));
    for j from 1 to max C do (
        f_j = C.dd_j;
        mapsList = append(mapsList, f_j); 
    );
    augC := chainComplex(mapsList);
    augC
    )

--latest version of consruction 3.6
buildCR = 
method(TypicalValue => ChainComplex
    , Options => {LengthLimit => 2}
    )
buildCR (ZZ,Module):= opts -> (g,M) -> (
    if g <= 0
    then error "expected g>0 for a syzygy of positive degree";
    n:= opts.LengthLimit;
--======================     
--Step 1: construct kappa
--Step 1.1: construct Q
    P := resolution(M, LengthLimit=>max(n,g+2));
    Pd := dual P;
    Pt := truncateComplex(g, P);
    Ptd := dual Pt;
    Q := Ptd[-(g-1)];
--Step 1.2: construct L - is this necessary? I don't think so.
    G := omega(g,P);
    Gd = dual G;
    L := resolution(Gd, LengthLimit => max(g+2, n));
--Step 1.3: construct the natural surjection; first factor
    toLiftFirstFactor := map(image(Pd.dd_(-g+1)), omega(1-g, Pd), id_(Pd_(1-g)));
--Step 1.4: construct the natural injection; second factor
    K := kernel(Pd.dd_(-g));
    I := image(Pd.dd_(-g+1));
    h := gens K;
    phi := Pd.dd_(-g+1)//h;
    toLiftSecondFactor := map (K, I, phi);
--and construct kappa.     
    kappa := toLiftSecondFactor * toLiftFirstFactor;
--=========================
--Step 2: lift kappa
--Step 2.1: build source of lifting
    D = truncateComplex(g,resolution(source kappa, LengthLimit => g));
--Step 2.1.1: augment with source kappa
    D' = augmentResolution(source kappa, D);
--Step 2.2: build target of lifting
    E = resolution(target kappa, LengthLimit => max(g+2, n));
--Step 2.2.1: augment with target kappa
    E' = augmentResolution(target kappa, E);
--and lift kappa     
    kappaLifted = extend(E', D', kappa)[1];
--==========================
--Step 3: build the degree g differential
--Step 3.1: build the factors
    w := map(G, P_g, id_(P_g));
    d := bidualityMap(G);
    lambda := map(Gd, E_0, id_(E_0));
--Step 3.2: combine
    alpha := (dual lambda)*d*w;
--Step 3.3: adjust source and target, might not be necessary
    alpha' := map((dual E_0),P_g,alpha);
--============================
--Step 4: build the source of the chain complex map
    S := new ChainComplex;
    S.ring = M.ring;
--Step 4.1: build the modules    
    for i from (g-1-max(g+2,n)) to g-1 do (
	S_i = dual E_(g-i-1); 
    );
    for i from g to max(g+2,n) do (
	S_i = P_i;
    );
--Step 4.2: build the differentials
    for i from (g-1-max(g+2,n)) to g-1 do (
	S.dd_i = dual E.dd_(g-i);
    );   
    for i from g+1 to max(g+2,n) do (
    	S.dd_i = P.dd_i;	
    );
    S.dd_g = alpha';
--=============================
--Step 5: build the target of the chain complex map
    T := new ChainComplex;
    T.ring = M.ring;
--Step 5.1: build the modules        
    for i from (g-1-max(g+2,n)) to g-1 do (
    	T_i = (dual D_(g-i-1));
    );
    for i from g to max(g+2,n) do (
	T_i = P_i;
    );
--Step 5.2: build the differentials
    for i from (g-1-max(g+2,n)) to g-1 do (
	T.dd_i = dual D.dd_(g-i);
    );   
    for i from g+1 to max(g+2,n) do (
    	T.dd_i = P.dd_i;	
    );
    T.dd_g = map(T_(g-1),P_g,P.dd_g);
--=============================
--Step 6: build the maps between the source and target
--    for i from (g-1-max(g+2,n)) to -1 do(
--	f_i = map(T_i,S_i, 0*(dual kappaLifted_(-g+1+i)));
--    );
    for i from (g-1-max(g+2,n)) to g-1 do(
	f_i = dual kappaLifted_(g-i-1);
    );	
    for i from g to max(g+2,n) do (
	f_i = id_(P_i);
    );
--==============================    
--Step 7: output T, S, and f_i as a hash table
output := new MutableHashTable;
output.ff = new MutableHashTable;
output.source = S;
output.target = T;
for i from (g-1-max(g+2,n)) to max(g+2,n) do (
    output.ff#i = f_i
    );
--==============================
--Step 8: make the chain complex map    
    cRes := map (T,S,i-> null); 
for i from (g-1-max(g+2,n)) to max(g+2,n) do (
    cRes_i = f_i
    );    
    cRes
--output
    )

buildMaps = method()
buildMaps(ZZ) := j -> (
    mapsList := ();
    if (j > max(g+2,n) 	or j< (g-1)-max(g+2,n)) 
    then error "integer out of bounds";
    for i from (g-1-max(g+2,n)) to max(g+2,n) do (
    if i == j then mapsList = append(mapsList, f_i)
    else mapsList = append(mapsList, null);
    );
    mapsList
    )
buildComplex = method()
buildComplex(ZZ) := t -> (
    mapsList := buildMaps(t);
    C = map(P,S,i -> mapsList_(i+(g-1+1)));
    C
    )
sign = method()
sign(ZZ) := j -> (
    if j >= 0 then return 1 else return (-1))

end
--------Test Code----------
restart
load "construction3-6v3.m2"
R = QQ[x,y,z]/ideal(x*y*z)
M = coker vars R
g = 3
n = 5
C = buildCR(g,M) --the output, as a chain complex map
--it works, that is it runs without errors, but it remains to check
--that the lifting kappaLifted in degrees 0, 1, 2 is correct.

--a test for strict equality of sources/targets
for i from (g-1-max(g+2,n)) to max(g+2,n) do (
      print (i, source C_i === C.source_i, target C_i === C.target_i)
      )
--a less strict test;
for i from (g-1-max(g+2,n)) to max(g+2,n) do (
      print (i, source C_i == C.source_i, target C_i == C.target_i)
      )
