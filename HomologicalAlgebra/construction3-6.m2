--natural map from a module to its double dual courtesy of Frank Moore
restart
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

CompleteResolution = new Type of MutableHashTable
globalAssignment CompleteResolution


-- what follow is version 0.1 of construction 3.6 from Avramov & Martsinkovsky
-- later versions will turn the construction into an object of the type
-- CompleteResolution
construction = method()
construction (ZZ, Module) := (g, M) -> (
     G = omega(g, M);
     Gd = dual G;
     L = resolution(Gd, LengthLimit => g+2); -- adjust later to suit user input
     P = resolution(M, LengthLimit => g+2);  -- ditto
     Pd = dual P;
     liftFirstFactor = map     
     truncateDualShift(g, P); -- via Kat's contribution
     
