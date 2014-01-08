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

export {subQuotientAsCokernel,identityMap,NCChainComplex}


needsPackage "NCAlgebra"


subQuotientAsCokernel = method()
subQuotientAsCokernel (NCMatrix, NCMatrix) := (M,N) -> (
   --- following Algorithm 6.3.1 in Boehm
   L := M | N;
   kerL := rightKernelBergman(L);
   rowsMN := #(M.source);
   kerL^(toList(0..(rowsMN-1)))
)

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
   assignDegrees(MtensN,newSource,newTarget)
)

Hom (NCMatrix,NCMatrix,ZZ) := (M,N,d) -> (
   R := coefficientRing ring M;
   sourceM := M.source;
   sourceN := N.source;
   targetM := M.target;
   targetN := N.target;
   error "err";
   R
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
///


identityMap = method()
identityMap (List, NCRing) := (L,R) -> (
   n := #L;
   B := coefficientRing R;
   I := ncMatrix applyTable(entries id_(B^n), e -> promote(e,R));
   assignDegrees(I,toList(n:0),L)
)

TEST ///
restart
needsPackage "NCAlgebraV2"
needsPackage "NCAlgebra"
B = threeDimSklyanin(QQ,{1,1,-1},{x,y,z})
M = ncMatrix {{x,y}}
N = ncMatrix {{x^2,y^2}}
subQuotientAsCokernel(M,N)

restart
needsPackage "NCAlgebraV2"
needsPackage "NCAlgebra"
R = QQ[w]/ideal(w^2+w+1)
B = threeDimSklyanin(R,{1,1,-1},{x,y,z})
M = ncMatrix {{x,y,0},{0,y,z}}
N = ncMatrix {{x,y}}
Hom(M,N,1)
L1 = identityMap({0},B)
K1 = L1 ** (transpose M)
L2 = identityMap({0,0,0},B)
K2 = L2 ** N
K = K1 | -K2
kerK = rightKernelBergman K
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

