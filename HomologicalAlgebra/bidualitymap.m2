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

completeRes = method()
completeRes (Module,ZZ) := (M,n) -> (
   --- R is a Gorenstein ring, and M a module over it.  This function
   --- returns the complete resolution to degree n
   R := ring M;
   Mres := res(M, LengthLimit => n);
   C := res(coker transpose Mres.dd_n, LengthLimit=>2*n);
   D := Hom(C,R^1);
   D[n - max C]
)

