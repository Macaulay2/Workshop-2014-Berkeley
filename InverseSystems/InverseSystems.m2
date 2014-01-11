newPackage(    
    "InverseSystems",
    Version => "0.1", 
    Date => "January 10, 2014",
    Authors => {{Name => "", 
	    Email => "", 
	    HomePage => ""}},
    Headline => "computes inverse systems",
    DebuggingMode => false
    )
     	       
	    
export {newFromDual,modFromDual,newToDual,newToDualTrunc,intersectInverseSystems}
        	
        	
newFromDual = method()
modFromDual = method()
newToDual = method()    
newToDualTrunc = method()
intersectInverseSystems = method(Dispatch => Thing)

newFromDual Matrix := Matrix => (f) -> (
    S := ring f;
    g := lcm first entries compress flatten monomials f;
    f1 := contract(transpose f, map(S^{degree g},S^1,matrix{{g}}));
    e := first exponents g;
    I := ideal for i to -1 + numgens S list S_i^(e_i+1);
    R := S / I;
    presentation prune trim coker (sub(gens ker sub(f1,R),S) | ((target f) ** gens I))
    )

modFromDual Matrix := Matrix => (F) -> (
    S := ring F;
    d := max for m in first entries compress flatten monomials F list sum first exponents m;
    L := matrix {{1_S}}; 
    for i from 1 to d+1 do L = L | gens power(ideal vars S,i); 
    pM := (L ** id_(target F)) * syz sub(contract((transpose L) * L,transpose F),vars S-vars S);
    presentation prune trim coker pM
    )

newToDual (ZZ,Matrix) := Matrix => (d,f) -> (
    S := ring f;
    g := product apply(generators S, v -> v^d);
    I := ideal for i to -1 + numgens S list S_i^(d+1);
    R := S / I;
    transpose contract(
	transpose mingens image (sub(gens ker sub(transpose f,R),S)),
	map(S^{degree g},S^1,matrix{{g}})))

newToDualTrunc (ZZ,Matrix) := Matrix => (d,f) -> (
    newToDual (d,f | ((target f) ** gens power(ideal vars ring f,d+1)))
    )

intersectInverseSystems Sequence :=  Matrix  => L -> intersectInverseSystems toList L

intersectInverseSystems List := Matrix => L -> (
    -- first check that all modules have the same target
    -- and the same base ring
    if #L === 0 then error "expected at least one argument";
    M := newFromDual(first L);
    for F in L do M = M | newFromDual(F);
    d := sum first exponents lcm first entries monomials flatten M;
    newToDual(d,M)
    )


beginDocumentation()


doc ///
   Key
      InverseSystems
   Headline
      for computations with Macaulay inverse systems	
   Description
      Text
         This package makes it possible to do computations with Macaulay inverse systems 
	 modules and ideals in a polynomial ring. 
      	  
      Example
    
   Caveat
      The package uses the contraction action of the polynomial ring on itself. 
   SeeAlso
      fromDual
      toDual
///

doc ///
   Key 
      newFromDual
   Headline 
      computes the submodule annihilating an inverse system
   Usage 
      N = newFromDual(F)
   Inputs
      F:Matrix 
         F is a matrix where the columns are the generators of an inverse system in   
	 the inverse system of a free module. 
   Outputs
      N:Matrix
         N is a presentation matrix for the module that has F as its inverse system.
   Description
      Text
      
      Example 
         F = matrix(QQ[x,y,z],{{x^4,y^4},{z^4,x^2*y^2}})
         N = newFromDual(F)
      Text 
         The command {\tt newFromDual} extends the command {\tt fromDual} so that it also 
         works for modules.
///

end
TEST ///
-- test code and assertions here
-- may have as many TEST sections as needed
///
