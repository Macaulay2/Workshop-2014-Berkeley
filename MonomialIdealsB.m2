-- Copyright 2014:  Sonja Mapes
-- You may redistribute this file under the terms of the GNU General Public
-- License as published by the Free Software Foundation, either version 2
-- of the License, or any later version.

------------------------------------------
------------------------------------------
-- Header
------------------------------------------
------------------------------------------

--if version#"VERSION" <= "1.4" then (
    needsPackage "SimplicialComplexes";
    needsPackage "Graphs";
    needsPackage "Posets";
--    )

newPackage(
    "MonomialIdealsB",
        Version => "1.0",
        Date => "9. January 2014",
        Authors => {
            {Name => "Sonja Mapes", Email => "smapes1@nd.edu", HomePage => "http://www.nd.edu/~smapes1/"}
        },
        DebuggingMode => true, Reload => true,
        Headline => "Package for processing monomial idea routines",
   --     Configuration => {
   --         "DefaultPDFViewer" => "open", -- "open" for Macs and "evince" for Linux
   --         "DefaultPrecompute" => true,
   --         "DefaultSuppressLabels" => true
   --         },
        DebuggingMode => false
    --    if version#"VERSION" > "1.0" then PackageExports => {
    --        "SimplicialComplexes",
    --        "Graphs",
    --        "Posets"
    --        }
        )

--if version#"VERSION" <= "1.4" then (
    needsPackage "SimplicialComplexes";
    needsPackage "Graphs";
    needsPackage "Posets";
  --  )

-- Load configurations
--posets'PDFViewer = if instance((options Posets).Configuration#"DefaultPDFViewer", String) then (options Posets).Configuration#"DefaultPDFViewer" else "open";
--posets'Precompute = if instance((options Posets).Configuration#"DefaultPrecompute", Boolean) then (options Posets).Configuration#"DefaultPrecompute" else true;
posets'SuppressLabels = if instance((options Posets).Configuration#"DefaultSuppressLabels", Boolean) then (options Posets).Configuration#"DefaultSuppressLabels" else true;

export {
   --symbols
   "monomialLabels",
   
   -- coordinatizations
   "minimalSQFIdeal",
   "labeling", 
   "possCoord",
   "minimalNSIdeal",  
   -- polarizations
   "polarRing",
   "polarization"
    }

 ------------------------------------------
------------------------------------------
-- Methods
------------------------------------------
------------------------------------------

------------------------------------------
--Coordinatizations:
------------------------------------------

-- takes a monomial ideal and outputs the minimal sqfree monomial ideal
-- input: monomial ideal
-- output: monomial ideal
minimalSQFIdeal = method()
minimalSQFIdeal(Ideal) := Ideal => (I) -> (
     lcmLat := lcmLattice I;
     meetIrred := meetIrreducibles lcmLat;
     ats := atoms lcmLat;     
     notAbove := apply(ats, atom-> select(#meetIrred, mis -> (meetIrred_mis)%atom != 0)); 
     s:= getSymbol "s";
     minSQFring := QQ(monoid[apply(#meetIrred, i-> s_i)]);
     ideal apply (notAbove, notAboves -> product apply (notAboves, i-> minSQFring_i))
     )

-- need option so that user can input own var  


-- function that assigns a labeling to a lattice
-- input lattice and labels and it just adds it to the cache so that the lattice is now "labeled" i.e. the data is now
-- saved as part of the data structure
-- labels need to be defined as monomials in some ring
labeling = method()
labeling(Poset,List) := (Poset) => (L, l) -> (
     if isLattice L then L.cache.monomialLabels = l else error "L is not a lattice";
     L
     ) 
 
 
 -- function that takes a labeling and produces the corresponding monomial ideal
-- inputs:  labeled lattice
-- outputs:  monomial ideal 
possCoord = method()
possCoord(Poset) := Ideal => (L) -> ( 
     if L.cache.monomialLabels == {} then error "L is not labeled" else
     ats := atoms L;     
     notAbove := apply(ats, atom-> select(#L.GroundSet, i -> i != 0 and compare(L,atom,L.GroundSet#(i)) == false )); 
     ideal apply (notAbove, notAboves -> product apply (notAboves, i-> L.cache.monomialLabels_i))
     )

-- function returns minimal non squarefree monomial ideal associated to a finite atomic lattice
minimalNSIdeal = method()
minimalNSIdeal(Poset) := Ideal => P -> (
     	Q := indexLabeling P;  
	Q' := dropElements(subposet(Q, meetIrreducibles Q), maximalElements Q); -- induced subposet on mi
	--n := #maximalChains Q';
	maxChains := symbol maxChains;
	OneHolder := symbol OneHolder;
	j := 0;
	D := apply(maximalChains Q', i -> length i);
	maxChains = select(maximalChains Q', i -> length i == max D); -- maxChains+D helps identify chains of max length
	L := new MutableList from apply(#Q.GroundSet, i -> OneHolder); 
	while maxChains!= {} do (
		apply(first maxChains, i -> L#i = j);
		j = j+1;
		Q' = subposet(Q', toList ( set(Q'.GroundSet) - set(first maxChains)));
		maxChains = maximalChains Q';
		); 
	n:= # unique toList select(L, placeholder -> placeholder =!= OneHolder);
	x := getSymbol "x";
	R := QQ(monoid[x_0..x_(n-1)]); 
  	L' := apply(L, placeholder -> if placeholder === OneHolder then 1 else R_(placeholder));
	-- convert L to monomials    
	labeling(Q, toList L');
	possCoord Q	
	)

------------------------------------------
--Polarizations:
------------------------------------------

polarRing = method ()
polarRing (List, List) := PolynomialRing => (L,M) -> (
    	if #L != #M then error "Inputs must be of same length";
     	A := apply(M, i -> i_0);
     	while L != {} do(
     		l := first L;
     		m := first M;
     		A = A | {m_0..m_(l-1)};
     		L = take(L,-(#L-1));
     		M = take(M,-(#M-1));
     		);
     	A = unique splice A;
     	R := QQ new Array from A
     	)

polarization = method ()
polarization (MonomialIdeal) := MonomialIdeal => M -> (
	--In this first part we are just constructing the ring
	L := matrix flatten apply(flatten entries gens M, i -> exponents i);
	L = apply(toList(0..numrows L - 1), i -> max flatten entries L _ {i}); 
	R := polarRing(L, toList vars(0..#L - 1));
	--This gives us each entry of the ideal in terms of exponents
	E := flatten apply(flatten entries gens M, i -> exponents i);
	B := {};
	C := {};
	while E != {} do (
		e := first E;
		N := hashTable  apply(toList(0..#L-1), i -> (vars(i), e_i));
		N = apply(pairs N, i -> toList i); 
		while N != {} do (
			n := first N;
			l := n_0;
			t := n_1;		
			if t != 0 then B = B | toList{l_0..l_(t-1)};
			N = take(N, -(#N-1));
			); 	
		C = C | {B};
		B = {};
		E = take(E,-(#E-1));
		); 
	C = apply(C, i -> splice i) 
	)


 
------------------------------------------
------------------------------------------
-- Documentation
------------------------------------------
------------------------------------------

beginDocumentation()

-- Front Page
doc ///
    Key
        MonomialIdealsB
    Headline
        a package for working with Monomial Ideals
    Description
        Text
          This package has functions relavant to monomial ideals.
        Text
            @SUBSECTION "Contributors"@
            --
          The following people have generously contributed code to the package:
          --@HREF("?","Jack Burkart")@,
	  --@HREF("http://www.nd.edu/~smapes1/","Sonja Mapes")@,
          --@HREF("?","Lindsay Piechnik")@
///



------------------------------------------
--Coordinatizations:
------------------------------------------
 
doc ///
    Key
       minimalSQFIdeal
        (minimalSQFIdeal,Ideal)
    Headline
        takes a monomial ideal and outputs the minimal sqfree monomial ideal
    Usage
        J = minimalSQRIdeal I
    Inputs
       I:Ideal
    Outputs
       J:Ideal
    Description
        Text
            This produces the minimal square free monomial ideal with the 
	    same lcm lattice as the given ideal. For the definitions, see
    -- math/0511032, \"Minimal monomial ideals and linear resolutions\", by "Jeffery Phan."
       Example
          R = QQ[x,y,z];
          I = ideal (x^4, x*y, y^7);
	  minimalSQFIdeal I
 
    SeeAlso
       minimalNSIdeal
       possCoord
       labeling
	
///

doc ///
    Key
      labeling
        (labeling,Poset,List)
    Headline
        assigns a labeling to a lattice
    Usage
        K = labeling(L,l)
    Inputs
        L:Poset
	l:List 
    Outputs
       K:Poset
    Description
        Text
            This takes a lattice and assigns a label to each position in lattice 
	    (note the order you list the labels in has to match the list of the ground set).
	    For the definitions, see
      -- math/1009.1430, \"Finite atomic lattices and resolutions of monomial ideals\", by "Sonja Mapes."
        Example
	    L = booleanLattice 3;
--	    L.l;
            l = {1,a,b,1,c,1,1,1};  -- note here the order you list the labels in has to match the list of the ground set
            labeling (L, l)
            
    SeeAlso
       minimalSQFIdeal
       minimalNSIdeal
       possCoord
       
	
///


doc ///
    Key
      possCoord
        (possCoord,Poset)
    Headline
        takes a labeling and produces the corresponding monomial ideal
    Usage
        I = possCoord L
    Inputs
        L:Poset
    Outputs
      I:Ideal
    Description
        Text
            This takes a labeled lattice and produces the corresponding monomial idea.  
	  as described in 
     	 -- math/1009.1430, \"Finite atomic lattices and resolutions of monomial ideals\", by "Sonja Mapes."
        Example
	   S = QQ[a,b,c];
	   L = booleanLattice 3
	   L.GroundSet;
	   groundSetMonomials = {1,a,b,1,c,1,1,1};  -- note here the order you list the labels in has to match the list of the ground set
	   labeling (L, groundSetMonomials);
	-- to see the result of using the labeling function type the following
	   L.cache.monomialLabels;
	   possCoord (L)
            
    SeeAlso
       minimalSQFIdeal
       minimalNSIdeal
       labeling
       
	
///

doc ///
    Key
      minimalNSIdeal
        (minimalNSIdeal,Poset)
    Headline
        takes a lattice and produces the most non square free minimal monomial ideal for that lattice
    Usage
        I = minimalNSIdeal P
    Inputs
        P:Poset
    Outputs
      I:Ideal
    Description
        Text
            This takes a lattice and produces the corresponding monomial idea with highest possible powers of each variable.  
	  For definitions and context see
    -- math/1009.1430, \"Finite atomic lattices and resolutions of monomial ideals\", by "Sonja Mapes."
        Example
	     P = poset {{a,b},{a,c},{a,d},{b,e},{c,e},{d,e}};
	     minimalNSIdeal P
            
    SeeAlso
       minimalNSIdeal
       labeling
       possCoord
	
///


------------------------------------------
--Polarizations:
------------------------------------------


--polarRing
doc ///
	Key
		polarRing
		(polarRing, List, List)
	Headline
		Returns a ring suitable for polarization of a monomial
	Usage
		R = polarRing (L,M)
	Inputs
		L:List
			A list of exponents
		M:List
			A list of variable names for the ring being created
	Outputs
		R:PolynomialRing
			A polynomial ring corresponding to the polarazation
	Description
		Text
			This method is used inside of the polarization method to make a suitable
			polynomial ring, and the user may find it useful to use this on small examples.
			The first list should be a list of non-negative integers, and then the next list 
			should be a list of the same size as the first with a list of desired non-indexed variables
		Example
			R = polarRing({2,4,5},{x,y,z})
			R = polarRing({1,1,1},{a,b,c})
		Text
			By convention, we use the rule that a 0 in the exponent list will still give out just one variable
			so that we don't lose information about the existence of some variable in the polynomial ring
		Example
			R = polarRing({0,0,3},{x,y,z})
	SeeAlso
		polarization
///

--polarization
---NOTE: THIS FUNCTION DOES NOT CURRENTLY WORK (HENCE THE EXAMPLES BEING COMMENTED OUT)---
doc ///
	Key
		polarization
		(polarization, MonomialIdeal)
	Headline
		Returns the polarization ideal of a given ideal
	Usage
		I = polarization M
		I = polarization m
	Inputs
		M:Ideal
			This ideal should be in a polynomial ring of less than 53 variables
	Outputs
		I:Ideal
			The new polarization ideal in the appropriate ring
	Description
		Text
			This method takes a monomial ideal and outputs its corresponding polarization ideal. This
			is the ideal generated by taking each generator of the input ideal and indexing the monomial 
			by its exponents, as in the following example.
	--	Example
	--		R = QQ[x,y,z]
	--		polarization(ideal(x^2, x*y, y^3))
 		Text
			It is important to note that we keep z in the new ring so that we do not lose any information about the
			original ring that the input ideal was in. The user may also just input any monomial if they are interested
			in the polarization of one element
	--	Example
	--		R = QQ[x,y,z]
	--		polarization(x^2*y*z^5)
	SeeAlso
		polarRing
			--For the code behind the construction of the ring
///



------------------------------------------
------------------------------------------
-- Tests
------------------------------------------
------------------------------------------

------------------------------------------
--Coordinatizations:
------------------------------------------

---minimalSQFIdeal test
TEST ///

	R = QQ[x,y,z]
        I = ideal (x^4, x*y, y^7)
	J = minimalSQFIdeal I
	J' =  ideal (s_3*s_4,s_0*s_1,s_1*s_3)
	assert(J === J')
	///

--- labeling test
-TEST ///
	L = booleanLattice 3
	L.l
        l = {1,a,b,1,c,1,1,1}  
        K = labeling (L, l)
	K'= --what I know the answer to be
	assert(K === K')
       	///
----------I DON'T KNOW HOW TO MAKE A LABELED POSET TO TEST THIS AGAINST-----------

--- possCoord test
TEST ///
        S = QQ[a,b,c]
	L = booleanLattice 3
	L.GroundSet
	groundSetMonomials = {1,a,b,1,c,1,1,1} 
	labeling (L, groundSetMonomials)
	I = possCoord L
     	J = ideal (b*c, a*c, a*b)
	assert(I === J)
       ///

---  minimalNSIdeal test
TEST ///
    P = poset{{0,1},{0,2},{0,3},{0,4},{0,5},{1,6},{2,6},{2,7},{2,8},{3,6},{4,8},{4,9},{5,9},{6,10},{7,10},{7,11},{8,11},{8,12},{9,12},{10,13},{11,13},{12,13}}
    I = minimalNSIdeal P
    J = ideal (x_1*x_2*(x_3)^3,x_0*x_1*(x_3)^2,x_0*x_2*(x_3)^3,(x_0)^3*x_1*x_3,(x_0)^3*x_1*x_2)
    assert(I === J)
///

end
------------------------------------------
--Polarizations:
------------------------------------------

--polarRing test
TEST ///
L = {1,2,3}
L' = {0,2,2}
M = {x,y,z}
assert(polarRing(L,M) === QQ[x_0, y_0, y_1, z_0, z_1, z_2])
assert(polarRing(L',M) === QQ[x_0, y_0, y_1, z_0, z_2])
///

-- polarization test
TEST ///
R = QQ[x,y,z]
I = monomialIdeal(x^2, y^2*x, z^4*x^3*y)
J = monomialIdeal(x^2, y^2)
m = x*y*z^4
p = x*y^3
assert(polarization I === monomialIdeal(x_0*x_1, x_0*y_0*y_1, x_0*x_1*x_2*y_0*z_0*z_1*z_2*z_3))
assert(polarization J === monomialIdeal(x_0*x_1, y_0*y_1))
assert(polarization m === monomialIdeal(x_0*y_0*z_0*z_1*z_2*z_3))
assert(polarization p === monomialIdeal(x_0*y_0*y_1*y_2))
///





end