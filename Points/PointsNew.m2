-- -*- coding: utf-8 -*-
newPackage(
	"PointsNew",
    	Version => "1.9", 
    	Date => "9 January 2014",
    	Authors => {
	     {Name => "Samuel Lundqvist", Email => "samuel@math.su.se", HomePage => "http://www.math.su.se/~samuel"}
	     },
    	Headline => "Computing with sets of affine points and functionals (i.e. FGLM conversion)",
    	DebuggingMode => true
    	)
-- Current developers Samuel Lundqvist
-- Past developers Gregory G. Smith, Mike Stillman, Stein A. Strømme
-- Contributors 
-- Acknowledgements Based upon the old Points package by the three past developers above.
-- Especially, the code for reduceColumn and part of "Points" is taken from there. 
export {
     points, --default std, option for GB.
     pointsByIntersection,
     makeRingMaps,
     nfPoints,
     separators,
     FGLM,
     standardMonomials, --standardMonomials
     GBFromMatrices
     --bordermonomials
     }


debug Core

makeRingMaps = method (TypicalValue => List)
makeRingMaps (Matrix, Ring) := List => (M,R) -> (
     K := coefficientRing R;
     pts := entries transpose M;
     apply(pts, p -> map(K, R, p))
     )

reduceColumn := (M,Mchange,H,c) -> (
     -- M is a mutable matrix
     -- Mchange is either null, or a matrix with same number of columns as M
     -- H is a hash table: H#r == c if column c has pivot for row r 
     -- returns true if the element reduces to 0
     M = raw M;
     if Mchange =!= null then Mchange = raw Mchange;
     r := rawNumberOfRows M - 1;
     while r >= 0 do (
	  a := M_(r,c);
	  if a != 0 then (
	       -- is there a pivot?
	       if not H#?r then (
		    b := 1//a;
		    rawMatrixColumnScale(M, b, c, false);
		    if Mchange =!= null then rawMatrixColumnScale(Mchange, b, c, false);		    
		    H#r = c;
		    return false;
		    )
	       else (
	       	    pivotc := H#r;
	       	    rawMatrixColumnChange(M, c, -a, pivotc, false);
		    if Mchange =!= null then rawMatrixColumnChange(Mchange, c, -a, pivotc, false);
	       ));
     	  r = r-1;
	  );
     true
     )



multMatrices = method (TypicalValue => List)
multMatrices (List, GroebnerBasis, Ring) := (std, Gb, R) -> (
     --Gb := forceGB matrix {gblist};
     K := coefficientRing R;
     s := #std;
     --print s;
     --m := 0;
     l := 0;
     RtoK := map(K,R);
     bb := {}; --The list of matrices to return 
     for i from 0 to #(gens R)-1 do (
	  m :=  mutableMatrix map(K^s, K^s, 0);
	  for j from 0 to s-1 do (
	       l = flatten entries (coefficients ((R_i * std#j) % Gb, Monomials => std))_1;
	       scan(#l, k -> m_(j,k) = RtoK(l#k));
	  );	
     	  bb = append (bb, matrix m); --We don't want it to be mutable anymore.
     );	   
     return bb;
)
	       

getCoefficientVector := (varIndex, cv, K, multMatrix, s, connectionvector) -> (
     -- varIndex is the index of the variable we are multiplying with monomial with.
     -- cv is coefficient vector of the monomial.
     -- K is the field that we are performing computations over. 
     -- multMatrix is the multiplication matrix with respect to the variable.
     -- s is the k-dimension of the ring.
     -- returns the coefficientvector of variable*monomial from multMatrix.
     -- varIndex = - 1 represents the unit
     if (varIndex == -1) then (
	 --l := matrix {flatten ({1_K,toList((s-1):0_K)})};
	 -- return l	
	  return matrix {connectionvector};
	  );
       return cv * multMatrix;
     )
	       
	       
addNewMonomial := (M,col,monom,maps) -> (
     -- M is an s by s+1 matrix, s=#points
     -- monom is a monomial
     -- maps is a list of s ring maps, which will give the values
     --  of the monom at the points
     -- replaces the 'col' column of M with the values of monom
     --    at the s points.
     scan(#maps, i -> M_(i,col) = maps#i monom)
     )

pointsByIntersection = method(TypicalValue => List)
pointsByIntersection (Matrix,Ring) := (M,R) -> (
     flatten entries gens gb intersect apply (
       entries transpose M, p -> ideal apply(#p, i -> R_i - p#i)))


 
splitGens = method()
splitGens (Matrix, Ring) := (M,R) -> (
       K := coefficientRing R;
     s := numgens source M;
     Fs := makeRingMaps(M,R);
     P := mutableMatrix map(K^s, K^(s+1), 0);
     PC := mutableMatrix map(K^(s+1), K^(s+1), 0);
     for i from 0 to s-1 do PC_(i,i) = 1_K;
     H := new MutableHashTable; -- used in the column reduction step
     thiscol := 0;
     esslist := {};   --The essential variables listed sorted ascending.
     nonesspoly := 0_K; -- The nonessential variables. Not necessarily sorted.
     l := sum(apply(toList(0 .. numgens R - 1),i-> R_i));
     m := 0;
	 while l != 0 do (
	     m = someTerms(l,-1,1);
    	     l = l - m;
             addNewMonomial(P,thiscol,m,Fs);
	     rawMatrixColumnScale(raw PC, raw(0_K), thiscol, false);
	     PC_(thiscol,thiscol) = 1_K;
	     isLT := reduceColumn(P,PC,H,thiscol);
	     if isLT then (
		 nonesspoly = nonesspoly + m;
	     --do nothing;
	     ) else (
	       -- we modify L, Lhash, thiscol, and also PC
	       esslist = append(esslist, m);
	      if (#esslist == s-1) then (
		  return (esslist, flatten entries (monomials (nonesspoly + l)));
		  );
	        thiscol = thiscol + 1;
	       );
	);
  if (nonesspoly == 0) then (
      return (esslist, {});
      );
  return (esslist, flatten entries (monomials (nonesspoly)));
  )
  


points = method(Options =>{groebnerBasis => false})
points (Matrix,Ring) := Sequence => o -> (M,R) -> (
     -- The columns of M form the points.  M should be a matrix of size
     -- n by s, where n is the number of variables of R
     K := coefficientRing R;
     s := numgens source M;
     -- The local data structures:
     -- (P,PC) is the matrix which contains the elements to be reduced
     -- Fs is used to evaluate monomials at the points
     -- H is a hash table used in Gaussian elimination: it contains the
     --    pivot columns for each row
     -- L is the sum of monomials which is still to be done
     -- Lhash is a hashtable: Lhash#monom = i means that only 
     --    essgens_i*monom, ..., essgens_n*monom should be considered
     -- G is a list of GB elements
     -- inG is the ideal of initial monomials for the GB
     essgens := {};
     nonessgens := {};
     splitGensCall := false;
     -- Rhough estimate whether or not we gain speed using the call to splitGens.
     if (s < numgens R and numgens R > 50) then (
     	 (essgens,nonessgens) = splitGens(M,R);
	 splitGensCall = true;
	 ) else (
	 essgens = sort gens R; --sorts in ascending wrt the monomial order of R
	 splitGensCall = false;
	 );	 
     Fs := makeRingMaps(M,R);
     P := mutableMatrix map(K^s, K^(s+1), 0);
     PC := mutableMatrix map(K^(s+1), K^(s+1), 0);
     for i from 0 to s-1 do PC_(i,i) = 1_K;
     H := new MutableHashTable; -- used in the column reduction step
     Lhash := new MutableHashTable; -- used to determine which monomials come next
     L := 1_R;
     Lhash#L = 0; -- start with multiplication by essgen_0.
     thiscol := 0;
     G := {};
     inG := trim ideal(0_R);
     inGB := forceGB gens inG;
     Q := {}; -- the list of standard monomials
     nL := 1;
     monom :=0;
     g := 0;
     isLT := 0;
     f:=0;
     while L != 0 do (
	  -- First step: get the monomial to consider
	  L = L % inGB;
	  monom = someTerms(L,-1,1);
	  L = L - monom;
	  -- Now fix up the matrices P, PC
          addNewMonomial(P,thiscol,monom,Fs);
	  if (o#groebnerBasis == true) then ( --if GB, we need to remember the ops
	      rawMatrixColumnScale(raw PC, raw(0_K), thiscol, false);
	      PC_(thiscol,thiscol) = 1_K;
          );
	  isLT = reduceColumn(P,PC,H,thiscol);
	  if isLT then (
	       -- we add to G, inG	  
	       inG = inG + ideal(monom);
	       inGB = forceGB gens inG;
	       if (o#groebnerBasis == true) then (
		   g = sum apply(toList(0..thiscol-1), i -> PC_(i,thiscol) * Q_i);
	       	   G = append(G, PC_(thiscol,thiscol) * monom + g);
	     	 );
	       )
	  else (
	       -- we modify L, Lhash, thiscol, and also PC
	       Q = append(Q, monom);
	       if (length Q == s and o#groebnerBasis == false) then ( --stds are done	     
     		   return (Q, inverse transpose matrix{apply(
			   Fs, f -> f(transpose matrix{Q}))});    		   
		   );
	       f = sum apply(toList(Lhash#monom .. #essgens - 1), i -> (
			 newmon := monom * essgens_i;
			 Lhash#newmon = i;
			 newmon));
	       nL = nL + size(f);
	       L = L + f;
	       thiscol = thiscol + 1;
	       )
	  );
      --We end up by computing the Groebner basis elements with leadterm
      --in nonessgen. For efficiency purposes, they were not added to L in the
      --first step. One could also use that we already have reduced them in
      --the getEssGen-procedure.
      if (splitGensCall == true and o#groebnerBasis == true) then (
    	  L = sum apply(toList(0..#nonessgens -1), i-> (nonessgens_i));	      
          while (L != 0 ) do (
	      monom = someTerms(L,-1,1);
	      L = L - monom;
	      addNewMonomial(P,thiscol,monom,Fs);
	      rawMatrixColumnScale(raw PC, raw(0_K), thiscol, false);
	      PC_(thiscol,thiscol) = 1_K;
              reduceColumn(P,PC,H,thiscol);	  
	      inG = inG + ideal(monom);
	    -- inGB = forceGB gens inG;
	    g = sum apply(toList(0..thiscol-1), i -> PC_(i,thiscol) * Q_i);
	    G = append(G, PC_(thiscol,thiscol) * monom + g);
      	    );
    	);  
    --     print("number of monomials considered = "|nL);
    if (o#groebnerBasis == false) then (
     	stds := transpose matrix{Q};
     	A := inverse transpose matrix{apply(Fs, f -> f stds)};
    	return (Q,A);
    	) else (
    	return (Q,inG,G);
     	)
    )


-- The separators of the points as linear combinations of the standard monomials
-- stds are the standard monomials returned by pointsMat
-- A is the inverse of the matrix B where b_(i,j) = std_i(p_j
separators = method()
separators (Matrix, Matrix) := (stds, A) -> (
     transpose (A) * (matrix entries stds)
     )


-- The normal form of a polynomial using Ainv and linear algebra
-- p is the polynomial of which we want to compute nf
-- phi is the list of ring maps returned from makeRingMaps 
-- stds is the list of standard monomials.
-- A is the inverse of the matrix B where b_(i,j) = std_i(p_j).
nfPoints = method();
nfPoints (RingElement, Matrix, List, Matrix) := (f, pointsMatrix, std, A) ->
   nfPoints(f, makeRingMaps(pointsMatrix, ring f), std, A);
nfPoints (RingElement, List, List, Matrix) := (p, phi, stds, A) -> (
     --Evaluate the vector on the points
     v := transpose matrix {apply (phi, r -> r p)};
     w := A * v;
     --Fix the stds
     stdsniceform := matrix {stds};
     --return the normal form
     first (first entries (stdsniceform*w))
     )


-- Return the stdmons as a list
standardMonomials = method()
standardMonomials(PolynomialRing, GroebnerBasis) := (S,Gb) -> (
     I := monomialIdeal(leadTerm (Gb));
     basisSmodI := flatten (entries (basis (S/I)));
     -- we want the monomials, not the residues, so we map from S/I to S. 
     SmodItoS := map(S,S/I);
     apply(basisSmodI, i -> SmodItoS(i))
     )
-- FGLM is computed by calling GBFromMatrices. This is done by 
-- setting up the multiplication matrices and then
-- connecting them by setting 1 = (1,0,...,0) * (e_1, ..., e_m)^t which is ok
-- since we assume e_1 < ... < e_m from which it follows that e_1 = 1.
FGLM = method()
FGLM (GroebnerBasis, PolynomialRing, Option) := (GS,S,monOrd) -> (   
     basisS := standardMonomials (S,GS);
     s := length basisS;
      K := coefficientRing S;
      connectionvector := flatten toList (1_K,toList(s-1:0_K));  --ok even if s-1<0
    return GBFromMatrices(multMatrices(basisS, GS, S), 
	connectionvector, S, monOrd);
    )

--Removes elements at the end of the list for which supp(m) > #copies.
--These elements are multiplies of ini(I) and should not be considered, cf the
--trick in the FGLM-paper.
removeElements := (L) -> (
    Lleast := 0;
     while L != {} do (
	  Lleast = first(L);
       if (#support(Lleast#0) > (Lleast#1)#2) then (
       L = drop(L,1);
       ) else (
       return L;
       );
  );
  return L;   
)

-- Computes a Groebner basis from the multiplication matrices
-- row j of the matrix i is  c_1,...,c_m, where 
-- x_i * e_j = c_1 e_1 + ... + c_m e_m. We do not need to know the basis
-- e_1,...,e_m, but we need to know how 1 is represented in the basis in order
-- to start the computation. The vector connectionvector has this info, i.e. 
-- connectionvector * (e_1,...,e_m)^t = 1. 
-- Examples: In FGLM, the connection vector is (1,0...,0). With respect to a 
-- separatorbasis, the connection vector is (1,...,1).
GBFromMatrices = method()
GBFromMatrices (List, List, PolynomialRing, Option) := 
(mm,connectionvector,S,monOrd) -> (      
     --from now on, we will compute over the ring R;
     R := newRing(S, monOrd);
     K := coefficientRing R;
     StoK := map(K,S);
     s := length connectionvector;   
     -- The local data structures:
     -- (P,PC) is the matrix which contains the elements to be reduced
     -- Fs is used to evaluate monomials at the points
     -- H is a hash table used in Gaussian elimination: it contains the
     --    pivot columns for each row
     -- L is the monomials to be done. An element in L is of the form
     -- (monom, (variable, coefflist, i)), where monom and variable is in R and
     -- the coefficientlist is the coefficient of (monom/variable) in basisS and
     -- i is the number of copies of the element in L.
     -- monom is only used for keeping L sorted and to be able to use MergePairs.   
     -- G is a list of GB elements
     -- inG is the ideal of initial monomials for the GB
     --Fs := makeRingMaps(M,R);
     P := mutableMatrix map(K^s, K^(s+1), 0);
     PC := mutableMatrix map(K^(s+1), K^(s+1), 0);
     for i from 0 to s-1 do PC_(i,i) = 1_K;
     H := new MutableHashTable; -- used in the column reduction step     
     -- essgens := getEss(gens R);
     essgens := gens R; --When removing this, the list below must also be changed 
     -- i.e. All noness should be inserted into L.
     L := {(1_R, (-1, {1},1))}; -- the list of potential elements. The list {1} is dummy.  
     Q := {}; -- the list of standard monomials
     thiscol := 0;
     G := {}; -- the list of Groebner basis elements to return
     numList := 1;
    Lleast := 0;
    cv := 0;
    monom := 0;
    while L != {} do (
	  -- First step: get the monomial to consider
	  Lleast = L#0;
	  L = drop(L,1);
	   --Pick the minimal elementet in L.
	  -- Now fix up the matrices P, PC
	  -- We are multiplying with (Lleast#1)#0;
          cv = getCoefficientVector(
	       (Lleast#1)#0, (Lleast#1)#1, K, mm#((Lleast#1)#0),s,connectionvector);	  
	  scan(s, i -> P_(i,thiscol) = cv_(0,i)); --add the column to P
	  rawMatrixColumnScale(raw PC, raw(0_K), thiscol, false);
	  PC_(thiscol,thiscol) = 1_K;
          isLT := reduceColumn(P,PC,H,thiscol);
	  monom = (Lleast#0);
	  if isLT then (
	       -- we add to G
	       g := sum apply(toList(0..thiscol-1), i -> PC_(i,thiscol) * Q_i);
	       G = append(G, PC_(thiscol,thiscol) * monom + g);
	       )
	  else (
	       -- we modify L and thiscol
	       Q = append(Q, monom);	  
	       newList := {};
	       for i from 0 to #essgens - 1 do (
	       	  numList = numList + 1;
	       	  newList = prepend((monom * essgens_i, (i, cv,1)),newList);
		  );
	       L = mergePairs(L,newList, (v,w)->(v#0,v#1,v#2 + w#2));

	       thiscol = thiscol + 1;
	       );
	  L = removeElements(L);
	  );
     --print("number of monomials considered = " numList);
     (Q,G,R)
     )



     
     
beginDocumentation()

document {
     Key => PointsNew,
     "A package to compute with points in affine spaces",
     }
document {
     Key => {nfPoints, (nfPoints,RingElement,List,List,Matrix)},
     Headline => "Normal form wrt standard monomials using linear algebra",
     Usage => "nfPoints(p,phi,std,A)",
     Inputs => {
     	  "p" => RingElement => "The polynomial for which we want to compute the normal form",
	  "phi" => List => "The ring maps.",
	  "std" => Matrix => "The standard monomials as a list",
	  "A" => Matrix => "The matrix returned by pointsstd, i.e. the inverse of the matrix defined by evaluating the standard monomials on the input points",
	  },
     Outputs => {RingElement => "The normal form of f wrt the standard monomials"},
     "Computing normal forms with respect to a vanishing ideal of points should be done by linear algebra and not by means of a Gröbner basis. 
     The timing below indicates that the speedup is drastic, even for toy examples. ",
     EXAMPLE lines ///
     
    loadPackage "PointsNew"
     pointsMatrix = random(ZZ^20, ZZ^30)      
      R = QQ[x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20]
      time (std,A) = points(pointsMatrix,R);
      g =     x1^3*x2^3;
      timing (nf1g = nfPoints(f, pointsMatrix, stdnew, Anew);)
      --If we could also call nfPoints with the ringmaps.
      --This is the recommended if numerous calls are to be done
      phi = makeRingMaps(pointsMatrix,R);
      timing (nf2g = nfPoints(f, pointsMatrix, stdnew, Anew);)
      nf1g == nf2g
      --  True 
      -- Now compute the normal form by means of the Groebner basis.
      timing ((Q,inG,G) = points(pointsMatrix,R, options=>groebnerBasis);)
      Gb = forceGB (matrix {G});
      timing(nf3f = f % Gb;)
      nf3f == nf2f
      --True
     ///
     }

document {
     Key => {makeRingMaps, (makeRingMaps,Matrix,Ring)},
     Headline => "evaluation on points",
     Usage => "makeRingMaps(M,R)",
     Inputs => {
     	  "M" => Matrix => "in which each column consists of the coordinates of a point",
	  "R" => PolynomialRing => "coordinate ring of the affine space containing the points",
	  },
     Outputs => {List => "of ring maps corresponding to evaluations at each point"},
     "Giving the coordinates of a point in affine space is equivalent to giving a
     ring map from the polynomial ring to the ground field: evaluation at the point.  Given a
     finite collection of points encoded as the columns of a matrix,
     this function returns a corresponding list of ring maps.",
     EXAMPLE lines ///
     M = random(ZZ^3, ZZ^5)
     R = QQ[x,y,z]
     phi = makeRingMaps(M,R)
     phi#2
     ///
     }



document {
     Key => {points, (points,Matrix,Ring)},
     Headline => "produces the ideal and initial ideal from the coordinates
     of a finite set of points",
     Usage => "(Q,inG,G) = points(M,R)",
     Inputs => {
     	  "M" => Matrix => "in which each column consists of the coordinates of a point",
	  "R" => PolynomialRing => "coordinate ring of the affine space containing the points",
	  },
     Outputs => {
          "Q" => List => "list of standard monomials",
 	  "inG" => Ideal => "initial ideal of the set of points",
 	  "G" => List => "list of generators for Grobner basis for ideal of points"
 	  },
     "This function uses the Buchberger-Moeller algorithm to compute a grobner basis
     for the ideal of a finite number of points in affine space.  Here is a simple
     example.",
           -- If we want the Groebner basis, we give that as an option. Now the 
      -- output is changed and the matrix A is replaced by generators for in(I)
      -- and the Groebner basis as a list.
     EXAMPLE lines ///
      M = random(ZZ^3, ZZ^5)
     R = QQ[x,y,z]
     (Q,inG,G) = points(M,R)
     monomialIdeal G == inG
     ///,
     PARA{},
     "Next a larger example that shows that the Buchberger-Moeller algorithm in ",
     TT "points", " may be faster than the alternative method using the intersection
     of the ideals for each point.",
     EXAMPLE lines ///
     R = ZZ/32003[vars(0..4), MonomialOrder=>Lex]
     M = random(ZZ^5, ZZ^150)
     time J = pointsByIntersection(M,R);
     time C = points(M,R);
     J == C_2  
     ///,
     SeeAlso => {pointsByIntersection},
       Caveat => "Program does not check that the points are distinct."
     }


document {
     Key => {pointsByIntersection, (pointsByIntersection,Matrix,Ring)},
     Headline => "Computes ideal of point set by intersecting maximal ideals",
     Usage => "pointsByIntersection(M,R)",
     Inputs => {
     	  "M" => Matrix => "in which each column consists of the coordinates of a point",
	  "R" => PolynomialRing => "coordinate ring of the affine space containing the points",
	  },
     Outputs => {
 	  List => "Grobner basis for ideal of a finite set of points",
 	  },
     "This function computes the ideal of a finite set of points by intersecting
     the ideals for each point.  The coordinates of the points are the columns in
     the input matrix ", TT "M", ".",
     EXAMPLE lines ///
     M = random(ZZ^3, ZZ^5)
     R = QQ[x,y,z]
     pointsByIntersection(M,R)
     ///,
     SeeAlso => {points},
     }

document {
     Key => {FGLM,  (FGLM,  GroebnerBasis, PolynomialRing, Option)},
       Headline => "Uses the FGLM algorithm to change a Groebner basis for a zero-dimensional ideal wrt to a monomial ordering mo1 to another 
     Groebner basis with respect to a monomial ordering mo2.",
     Usage => "G2 = points(std,G1,R,mo2)",
     Inputs => {
	 "G1" => GroebnerBasis => "A Groebner basis in R1",	
	  "R1" => PolynomialRing => "A polynomial ring",	  
	  "mo2" => Option =>"The output monomial ordering"
	  },
     Outputs => {  "Q2" => List => "A list of the standard monomials wrt to mo2",
	 "G2" => List => "A list of the Groebner basis wrt to mo2",
	 "R2" => PolynomialRing => "The polynomial ring where G2 lives"	   	 
		  },
     EXAMPLE lines ///
     M = random(ZZ^20, ZZ^40);
          
     -- 1. Compute a Gröbner basis of I(M) with respect 
     -- to DegRevLex using the BM-algorithm.
     R1 = QQ[x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,
	 x12,x13,x14,x15,x16,x17,x18,x19,x20];
     timing ((Q1,inG1,G1) = points(M,R1,groebnerBasis=>true);)
     Gb1 = forceGB matrix {G1};
     I1 = ideal gens Gb1;
     
     -- 2. Convert the basis to a Lex-basis using FGLM.
     timing( (Q2,G2,R2) = FGLM(Gb1,R1, MonomialOrder => Lex);)
     
     -- 3. Compute a Gröbner basis of I(M) with respect 
     -- to Lex using the BM-algorithm.
     R3 = newRing(R1, MonomialOrder => Lex)  
     timing((Q3,inG3,G3) = points(M,R3,groebnerBasis=>true);)
     
     --Map the result from FGLM (which is in R2) to R3
     R2toR3 = map(R3,R2);  
     G2mappedtoR3 = apply(G2, p -> R2toR3(p));
    --Check that they are equal
     gens forceGB matrix  {sort G2mappedtoR3} == gens forceGB matrix {sort G3}
     -- true
     
     -- CONTINUE HERE. Q: Determine the order in which the functions return
     -- elements.
     -- Now, to show that we gain speed, 
     -- compute a Lex Gröbner basis from the generators of IR.
     -- First map I from R1 to R2
     R1toR2 = map(S2,R);
     IS2 = ideal RtoS2 gens IR;
     timing (BuchbergerLex = gens gb IS2)
    -- 30.9826 seconds
     --Check that the result agrees with the FGLM result
     sort(BuchbergerLex) == gens forceGB matrix  {sort FGLMLexGb}
     -- true
     ///,
     SeeAlso => {points},
     }


TEST ///
R = ZZ/32003[vars(0..4), MonomialOrder=>Lex]
M = matrix(ZZ/32003,  {{0, -9, 4, -2, -4, -9, -10, 6, -8, 0}, 
            {1, 0, -10, 9, 3, -4, 1, 1, -10, -3}, 
	    {5, 7, -4, -5, -7, 7, 4, 6, -3, 2}, 
	    {2, 8, 6, -6, 4, 3, 8, -10, 7, 8}, 
	    {-9, -9, 0, 4, -3, 9, 4, 4, -4, -4}})
phi = makeRingMaps(M,R)
apply (gens(R),r->phi#2 r)
assert ( {4, -10, -4, 6, 0} == apply (gens(R),r->phi#2 r) )

J = pointsByIntersection(M,R);
C = points(M,R);
assert ( J == C_2 )
assert ( C_1 == ideal(e^6,d*e^3,d^2*e,d^3,c,b,a) )
assert ( C_0 == sort apply (standardPairs monomialIdeal C_2, p -> p#0) )
assert (
     (pointsMat(M,R))#0 == 
      matrix(ZZ/32003, {{1, -9, 81, -729, 6561, 4957, 2, -18, 162, 4}, {1, -9, 81, -729, 6561,
      4957, 8, -72, 648, 64}, {1, 0, 0, 0, 0, 0, 6, 0, 0, 36}, {1, 4, 16, 64, 256, 1024,
      -6, -24, -96, 36}, {1, -3, 9, -27, 81, -243, 4, -12, 36, 16}, {1, 9, 81, 729, 6561,
      -4957, 3, 27, 243, 9}, {1, 4, 16, 64, 256, 1024, 8, 32, 128, 64}, {1, 4, 16, 64,
      256, 1024, -10, -40, -160, 100}, {1, -4, 16, -64, 256, -1024, 7, -28, 112, 49}, {1,
      -4, 16, -64, 256, -1024, 8, -32, 128, 64}})
)
assert ( first entries transpose (pointsMat(M,R))#1 == C_0 )
///
end
toString C_1
restart
errorDepth = 0


uninstallPackage "Points"
installPackage "Points"
R = ZZ/32003[vars(0..4), MonomialOrder=>Lex]
M = matrix(ZZ/32003,  {{0, -9, 4, -2, -4, -9, -10, 6, -8, 0}, 
            {1, 0, -10, 9, 3, -4, 1, 1, -10, -3}, 
	    {5, 7, -4, -5, -7, 7, 4, 6, -3, 2}, 
	    {2, 8, 6, -6, 4, 3, 8, -10, 7, 8}, 
	    {-9, -9, 0, 4, -3, 9, 4, 4, -4, -4}})

phi = makeRingMaps(M,R)
apply (gens(R),r->phi#2 r)
assert ( {4, -10, -4, 6, 0} == apply (gens(R),r->phi#2 r) )


phi#2
time J = pointsByIntersection(M,R)
transpose matrix{oo}



time C = points(M,R)
transpose gens ideal C_2

M = random(ZZ^3, ZZ^5)
R = QQ[x,y,z]
phi = makeRingMaps(M,R)
apply (gens(R),r->phi#2 r)
phi#2

R = ZZ/32003[vars(0..4), MonomialOrder=>Lex]
M = random(ZZ^5, ZZ^150)

time J = pointsByIntersection(M,R);
transpose matrix{oo}

time C = points(M,R);
transpose gens ideal C_2
assert(J == C_2)

R = ZZ/32003[vars(0..4)]

K = ZZ/32003
R = K[vars(0..7), MonomialOrder=>Lex]
R = K[vars(0..7)]
M = random(K^8, K^500)
time C = points(M,R);
time J = pointsByIntersection(M,R);
assert(C_2 == J)

K = ZZ/32003
R = K[x_0 .. x_39]
M = random(K^40, K^80)
time C = points(M,R);


getColumnChange oo_0
apply(Fs, f -> f(a*b*c*d))
B = sort basis(0,2,R)
B = sum(flatten entries basis(0,2,R))
B = matrix{reverse terms B}
P = transpose matrix {apply(Fs, f -> f (transpose B))}
B * syz 
transpose oo
 -- column reduction:

P = mutableMatrix P 
H = new MutableHashTable
reduceColumn(P,null,H,0)
reduceColumn(P,null,H,1)
P
reduceColumn(P,null,H,2)
reduceColumn(P,null,H,3)
reduceColumn(P,null,H,4)
reduceColumn(P,null,H,5)
reduceColumn(P,null,H,6)
reduceColumn(P,null,H,7)
reduceColumn(P,null,H,8)
reduceColumn(P,null,H,9)
P
reduceColumn(P,null,H,10)
reduceColumn(P,null,H,11)
reduceColumn(P,null,H,12)
P

M = matrix{{1,2,3,4}}

K = ZZ/32003
M ** K

-- Local Variables:
-- compile-command: "make -C $M2BUILDDIR/Macaulay2/packages PACKAGES=Points pre-install"
-- End:
