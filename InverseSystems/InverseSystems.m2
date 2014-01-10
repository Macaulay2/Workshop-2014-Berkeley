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
    newToDual (d,f | ((target f) ** gens power(ideal vars ring f,d)))
    )

intersectInverseSystems Sequence :=  Matrix  => L -> intersectInverseSystems toList L

intersectInverseSystems List :=  Matrix => L -> (
    -- first check that all modules have the same target
    -- and the same base ring
    if #L === 0 then error "expected at least one argument";
    M := matrix {{}};
    for F in L do M = M | newFromDual(F);
    d := sum first exponents lcm first entries monomials M;
    newToDual(d,M)
    )


beginDocumentation()

document {
    Key => InverseSystems,
    Headline => "for computations with Macaulay inverse systems",
    Caveat => "",
    SeeAlso => ""
    }

document {
       Key => newFromDual,
       Headline => "computes the submodule annihilating an inverse system",
       Usage => "N = newFromDual(F)",
       Inputs => {
	  "F" => {"F ", TO Matrix, " with entries in a polynomial ring"},
	  },
       Outputs => { 
	   "N" => {"N ", TO Matrix, " with entries in a polynomial ring"},
	  }, 
       PARA{}, "The input matrix ", TT "F", " represents the generators of 
       an inverse system and the output, ", TT "N", " is a presentation matrix 
       for the submodule annihilating ", TT " F", 
       
       EXAMPLE {
	   "F = matrix(QQ[x,y,z],{{x^4,y^4},{z^4,x^2*y^2}}",
	   "N = newFromDual(F)"
	   }
}

end
TEST ///
-- test code and assertions here
-- may have as many TEST sections as needed
///
