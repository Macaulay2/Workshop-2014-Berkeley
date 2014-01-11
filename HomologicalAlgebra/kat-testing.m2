--natural map from a module to its double dual courtesy of Frank Moore
bidualityMap = M -> (
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

--version 4 is below, it incorporates some of the previous methods
--and also (hopefully) fixes the indexing issue.
constructionV4 = 
method(TypicalValue => ChainComplex
     , Options => {LengthLimit => 2}
     )
constructionV4 (ZZ,Module):= 
     opts -> 
     (g,M) -> (
     if g <= 0
     then error "expected g>0 for a syzygy of positive degree";
     n:= opts.LengthLimit;
     P := resolution(M, LengthLimit=>max(n,g+2));
--     P := resolution(M, LengthLimit=>g+2);
     G := omega(g,P);
     Pd := dual P;
     toLiftFirstFactor := map(image(Pd.dd_(-g+1)), omega(1-g, Pd), id_(Pd_(1-g)));
     K := kernel(Pd.dd_(-g));
     I := image(Pd.dd_(-g+1));
     h := gens K;
     phi := Pd.dd_(-g+1)//h;
     toLiftSecondFactor := map (K, I, phi);
     kappa := toLiftSecondFactor * toLiftFirstFactor;
     Pt := truncateComplex(g, P);
     Ptd := dual Pt;
     Q := Ptd[-(g-1)];
--     Qd = dual Q;
--     augP := P;
--    augP.dd_(0) := map(M, P_0,id_(P_0)); for some reason,
--    this gives an error if there's a colon, but not if it isn't there
--     augP.dd_(0) = map(M,P_0,id_(P_0));
--     augPs = augP[-1];
     Gd := dual G;
     L := resolution(Gd, LengthLimit => max(g+2, n));
     augL := L;
     augL.dd_0 = map(Gd, L_0,id_(L_0));
     augLs = augL[-1];
     kappaLifted := (extend(augLs,Q,kappa))[1];
     --gives error maps not composable
     w := map(G, P_g, id_(P_g));
     d := bidualityMap(G);
     Ld := dual L;
--     L := resolution(Gd, LengthLimit => g+2);
    lambda := map(Gd, L_0, id_(L_0));
    lambdaDual := dual lambda; -- now implemented in M2 1.6
-- here are the significant changes to version 3;
-- build the source of the chain complex map
    S := new ChainComplex;
-- put in the ring
    S.ring = M.ring;
-- put in the modules
    for i from (g-1-max(g+2,n)) to g-1 do (
	S_i = Ld_(-g+1+i); 
	);
    for i from g to max(g+2,n) do (
	S_i = P_i;
	);
-- S_g = P_g;
--put in the differentials
    for i from (g-1-max(g+2,n)) to g-1 do (
	S.dd_i = Ld.dd_(-g+1+i);
	);  
    for i from g+1 to max(g+2,n) do (
    	S.dd_i = P.dd_i;	
    	);
    S.dd_g = lambdaDual*d*w;
--build the target of the chain complex map
--actually this is already built as it is P    
--build the maps between the source and target
--and name the maps consistently
    for i from (g-1-max(g+2,n)) to g-1 do(
	f_i = dual kappaLifted_(-g+1+i);
	);
    for i from g to max(g+2,n) do (
	f_i = id_(P_i);
	);
--make the chain complex map
    cRes := map (P,S,i-> f_i); 
    cRes
    )

end

--------Test Code----------
restart
load "kat-testing.m2"

-- example for testing constructionV2
R = QQ[x,y,z]/ideal(x*y*z)
M = coker map(R^1,,{gens R})
g=3
n=5
constructionV4(g,M)


--This code here checks if the source and target of the f_i maps are what they should be
for i from (g-1-max(g+2,n)) to max(g+2,n) do (
      print (i, source f_i === S_i, target f_i === P_i)
      )


     



