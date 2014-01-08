newPackage("NCAlgebraV2",
     Headline => "Additional functions for NCAlgebra",
     Version => "0.1",
     Date => "January 6th, 2014",
     Authors => {
	  {Name => "Frank Moore",
	   HomePage => "http://www.math.wfu.edu/Faculty/Moore.html",
	   Email => "moorewf@wfu.edu"},
	  {Name => "Andrew Conner",
	   HomePage => "http://www.math.wfu.edu/Faculty/Conner.html",
	   Email => "connerab@wfu.edu"},
          {Name => "Courtney Gibbons",
	   HomePage => "http://people.hamilton.edu/cgibbons",
	   Email => "crgibbon@hamilton.edu"}},
     AuxiliaryFiles => true,
     DebuggingMode => true
     )

export {subQuotientAsCokernel,homologyAsCokernel,identityMap,NCChainComplex}

debug needsPackage "NCAlgebra"

subQuotientAsCokernel = method()
subQuotientAsCokernel (NCMatrix, NCMatrix) := (M,N) -> (
   --- following Algorithm 6.3.1 in Boehm
   L := M | N;
   kerL := rightKernelBergman(L);
   rowsMN := #(M.source);
   tempKer := rightMingens (kerL^(toList(0..(rowsMN-1))));
   tempKer
)

homologyAsCokernel = method()
homologyAsCokernel(NCMatrix,NCMatrix) := (M,N) -> (
    if M*N != 0 then return "Error: maps do not compose to zero"
    else (
    B := N.ring;
    Z := Z = zeroMap((N.target),(N.source),B);
    kerM := rightKernelBergman(M);
    subQuotientAsCokernel(kerM,N)
    )
)

TEST ///
restart
needsPackage "NCAlgebraV2"
debug needsPackage "NCAlgebra"
B = threeDimSklyanin(QQ,{1,1,-1},{x,y,z})
M = ncMatrix {{x^2,y^2,z^3}}
Msyz = rightKernelBergman M
test = subQuotientAsCokernel(Msyz,Msyz)
///

--NCMatrix ** Matrix := 
--Matrix ** NCMatrix := 
NCMatrix ** NCMatrix := (M,N) -> (
   entriesM := entries M;
   MtensN := ncMatrix applyTable(entriesM, e -> e*N);
   --- now we must assignDegrees to make make them compatible
   --- with the maps M and N
   newSource := flatten apply(#(M.source), i ->
         apply(#(N.source), j -> ((M.source)#i)+((N.source)#j)));
   newTarget := flatten apply(#(M.target), i ->
         apply(#(N.target), j -> ((M.target)#i)+((N.target)#j)));
   assignDegrees(MtensN,newTarget,newSource)
)

NCMatrix ++ NCMatrix := (M,N) -> (
   B := ring M;
   urZero := zeroMap(M.target,N.source,B);
   lrZero := zeroMap(N.target,M.source,B);
   ds := ncMatrix {{M,urZero},{lrZero,N}};
   assignDegrees(ds,M.target | N.target, M.source | N.source)
)

-------------------------------------------
--- NCChainComplex Methods ----------------
-------------------------------------------
NCChainComplex = new Type of HashTable

resolution NCMatrix := opts -> M -> (
   i := 0;
   numSyz := if opts#LengthLimit === infinity then numgens ring M else opts#LengthLimit;
   currentM := M;
   syzList := {M} | while (i < numSyz and currentM != 0) list (
      newM := rightKernelBergman currentM;
      currentM = newM;
      currentM
   ) do i = i+1;
   new NCChainComplex from apply(#syzList, i -> (i,syzList#i))
)

betti NCChainComplex := opts -> C -> (
    firstbettis := flatten apply(
    	keys (tally (C#0).target), 
    	i -> {(0,(C#0).target,i) => (tally (C#0).target)_i}
    );
    lastbettis := flatten flatten apply(#C-1, j -> 
	apply(
    	    keys (tally (C#j).source), 
    	    i -> {(j+1,(C#j).source,i) => (tally (C#j).source)_i}
	    )
	);
    L := firstbettis | lastbettis;
    B := new BettiTally from L
)


TEST ///
restart
needsPackage "NCAlgebraV2"
needsPackage "NCAlgebra"
B = threeDimSklyanin(QQ,{1,1,-1},{x,y,z})
M = ncMatrix {{x,y,z}}
Mres = res M
betti Mres
N = ncMatrix {{x,y,z^2}}
Nres = res N
betti Nres
--- The resolution of the following NCMatrix has a 0 kernel 
--- in homological degree 2
L = ncMatrix {{x^2,x*z,y}};
Lres = res L
betti Lres
rightKernelBergman(Lres#2)
///

-- Twist --     
NCMatrix Array := (M,n) -> (
    if #n != 1 then return "Error: Please enter a single integer" else
    M**(assignDegrees(ncMatrix {{promote(1,M.ring)}},{-1*n#0},{-1*n#0}))
    )

TEST ///
restart
needsPackage "NCAlgebraV2"
needsPackage "NCAlgebra"
B = threeDimSklyanin(QQ,{1,1,-1},{x,y,z})
M = ncMatrix {{x,y,z}}
M[1]
(M[1]).source
///

identityMap = method()
identityMap (List, NCRing) := (L,R) -> (
   n := #L;
   B := coefficientRing R;
   I := ncMatrix applyTable(entries id_(B^n), e -> promote(e,R));
   assignDegrees(I,L,L)
)

identityMap (ZZ,NCRing) := (n,R) -> identityMap(toList(n:0),R)

zeroMap = method()
zeroMap (List, List, NCRing) := (tar,src,B) -> (
   R := coefficientRing B;
   myZero := ncMatrix applyTable(entries map(R^#tar,R^#src,0), e -> promote(e,B));
   assignDegrees(myZero,tar,src);
   myZero
)

Hom (NCMatrix,NCMatrix,ZZ) := (M,N,d) -> (
   B := ring M;
   Nsyz := rightKernelBergman N;  -- be careful if Nsyz is zero!
   L1 := identityMap(N.target,B);
   K1 := L1 ** (transpose M);
   L2 := identityMap(M.source,B);
   L3 := identityMap(M.target,B);
   L4 := identityMap(N.source,B);
   K2 := N ** (transpose L2);
   K3 := N ** (transpose L3);
   K4 := L4 ** (transpose M);
   K5 := Nsyz ** (transpose L2);
   myZeroMap := zeroMap(K3.target,K5.source,B);
   K1ent := entries K1;
   K2ent := entries K2;
   K3ent := entries K3;
   K4ent := entries K4;
   K5ent := entries K5;
   myZeroMapEnt := entries myZeroMap;
   --K := K1|K2;
   --H := (K3 | myZeroMap) || (K4 | K5);
   --H = K3 || K4   -- do this if Nsyz == 0
   K1' := matrix apply(#(K1.target), i -> apply(#(K1.source), j -> 
	leftMultiplicationMap(K1ent#i#j, d - (K1.source)#j, d - (K1.target)#i)));
   K2' := matrix apply(#(K2.target), i -> apply(#(K2.source), j -> 
	rightMultiplicationMap(-K2ent#i#j, d - (K2.source)#j, d - (K2.target)#i)));
   K3' := matrix apply(#(K3.target), i -> apply(#(K3.source), j -> 
	rightMultiplicationMap(-K3ent#i#j, d - (K3.source)#j, d - (K3.target)#i)));
   K4' := matrix apply(#(K4.target), i -> apply(#(K4.source), j -> 
	leftMultiplicationMap(-K4ent#i#j, d - (K4.source)#j, d - (K4.target)#i)));
   K5' := matrix apply(#(K5.target), i -> apply(#(K5.source), j -> 
	rightMultiplicationMap(-K5ent#i#j, d - (K5.source)#j, d - (K5.target)#i)));
   myZeroMap' := matrix apply(#(myZeroMap.target), i -> apply(#(myZeroMap.source), j -> 
   	   rightMultiplicationMap(-myZeroMapEnt#i#j, d - (myZeroMap.source)#j, d - (myZeroMap.target)#i)));
   K' := K1'|K2';
   H' := matrix {{K3',myZeroMap'},{K4',K5'}};
   --H' = matrix {{K3'},{K4'}}  -- do this if Nsyz == 0
   myHom := prune ((ker K') / (image H'));
   homGens := mingens image(gens image myHom.cache.pruningMap)^(toList(0..(numgens source K1' - 1)));
   basisMatr := fold(apply(#(K1.source), i -> basis(d-(K1.source)#i,B)), (a,b) -> a ++ b);
   flattenedMatrs := basisMatr * homGens;
   retVal := apply(apply(#(flattenedMatrs.source), i -> flatten entries flattenedMatrs_{i}), L -> ncMatrix pack(#(N.target),L));
   retVal
)

TEST ///
restart
needsPackage "NCAlgebraV2"
needsPackage "NCAlgebra"
B = threeDimSklyanin(QQ,{1,1,-1},{x,y,z})
M = ncMatrix {{x,y}}
N = ncMatrix {{x^2,y^2}}
subQuotientAsCokernel(M,N)
///

TEST ///
restart
needsPackage "NCAlgebraV2"
needsPackage "NCAlgebra"
B = threeDimSklyanin(QQ,{1,1,-1},{x,y,z})
R = coefficientRing B
M = ncMatrix {{x,y,0},{0,y,z}}
N = ncMatrix {{x,y},{x,y}}
Hom(M,N,2)
///

end

--- bug fix/performance/interface improvements
------------------------------------
--- Testing!

--- additions in the near future
------------------------------------
--- skewPolynomialRing with NCRing bases
--- Finish right mingens
--- Implement left kernels and mingens etc (opposite ring now done)
--- Make sure that trivial ideals are handled correctly
--- Make sure constructions over the tensor algebra are handled correctly
--- isFiniteDimensional?
--- 'basis' for f.d. algebras?
--- Generating set for algebras not over a field
--- fix coefficients to return a pair again.

--- other things to add or work on in due time
-----------------------------------
--- notation to refer to a certain graded piece of an algebra e.g. A_3
--- a dimension function for graded pieces dim(A,17)
--- Dare I say it, Diamond Lemma?
--- Make Quotients of Quotients work.
--- NCRingMap kernels (to a certain degree)  -- Not sure I can do this with
---   Bergman, can't use block orders in bergman.
--- anick resolution
--- NCModules (?) (including module gb (via simple), hilbert series, modulebettinumbers)
--- Work on reduction code a bit?
--- Hilbert series for modules (and sided ideals)
--- enveloping algebras, tensor products of algebra, Hochschild (co)homology?

--- symbols in hash tables of exported types
----
restart
uninstallPackage "NCAlgebraV2"
installPackage "NCAlgebraV2"
needsPackage "NCAlgebraV2"
viewHelp "NCAlgebraV2"

restart
uninstallPackage "NCAlgebra"
installPackage "NCAlgebra"
needsPackage "NCAlgebra"
viewHelp "NCAlgebra"

