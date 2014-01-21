newPackage(
	"IdealClosure",
	Version => "0.5",
	Date => "1-20-14",
	Authors => {
		{  Name => "Douglas Leonard", Email => "leonada@auburn.edu" },
		{  Name => "Nicholas Armenoff", Email => "nicholas.armenoff@uky.edu" }
	},
	Headline => "Local view of integral closures of ideals",
	DebuggingMode => true,
	Reload => true
)

export {
	qthIdealClosure, --written before the January 2014 workshop
-- to be written	rationalIdealClosure
	qthReesAlgebra, -- started at the workshop, finished a week after
-- to be written        rationalreesalgebra
        reesring,    --setup for qthReesAlgebra   
	Exponent,
	Idealsize,
	Modulesize,
	Units,
	Verb
};

-- PURPOSE : To implement Leonard's qth-power algorithm 
--           to compute the integral closure C(I,Q(A)) 
--           of an ideal I of the integral closure Q(A)
--           of a domain A=R/J,
--           an integral extension of a polynomial ring P
--	     over ZZ/q as qthIdealClosure.
--           (rationalIdealClosure, when written will be a lifting to QQ)
--           (reesalgebra, when written, will be the non-homogeneous Rees-algebra-type
--            ideal that can be used to determine in which C(I^n,Q(A)) an element of C(A) lives)
--   INPUT : multivariate polynomial ring, inputring, over ZZ/q
--           an ideal, quot, defining an affine domain,
--           a list, idealgens, of an ideal in the domain, 
--           wt, a list corresponding to the local monomial order weights 
--  OUTPUT : list of generators of the closure for qthIdealClosure or rationalIdealClosure
--           ideal of relations for reesalgebra
--           
--	     
--           
-- COMMENTS :   
--            (1)  lift to QQ not written. 
--            (2) The conjecture is that Q=q^Exponent 
--                at least the degree of the extension of A over P is sufficient.
--            (3) Macaulay2 isn't really set up for local monomial orderings.

------------------------------------------------------------------------
-----------------------  Code  -----------------------------------------
------------------------------------------------------------------------

------------------------------------------------------------------------
--qthIdealClosure code for Macaulay2 as of August 27, 2012--------------
--with the example inputs used to test it added since-------------------
------------------------------------------------------------------------
qthIdealClosure=method(
     Options=>{
	  Exponent=>1,Idealsize=>30,Modulesize=>40,Units=>0,Verb=>0
	  });
-------------------------------------------------------------
qthIdealClosure(Ring,List,List,List):=
   o->(inputring,quot,idealgens,wt)->(
   ----------------------------------------------------------
   --new ring setup------------------------------------------
   ----------------------------------------------------------
   q:=char inputring;
   Q:=q^(o.Exponent);
   idealsize:=o.Idealsize;
   modulesize:=o.Modulesize;
   verbosity := o.Verb;
   extwt:=for i to #wt-1 list ((for j to idealsize+modulesize-1 list 0)|wt#i);
   varno:=numgens inputring;
   T := local T;
   G := local G;
   B := local B;
   X := local X;
   U := local U;
   R:=ZZ/q[T_0..T_1,
        G_0..G_(idealsize-1),
        B_0..B_(modulesize-1),
	X_0..X_(varno-1-o.Units),
	U_0..U_(o.Units-1),
        Degrees=>(for i to 1 list{1,0})|
	         (for i to idealsize+modulesize-1 list {1,1})|
		 (for i to varno-1 list {1,0}),
        MonomialOrder=>{
           (2),
           (Weights=> (for i to idealsize-1 list -1)|(for i to modulesize+varno-1-o.Units list 0),
            Weights=> (for i to idealsize-1 list 0)|(for i to modulesize-1 list -1)|(for i to varno-1-o.Units list 0),
               for i to #wt-1 list Weights=> extwt#i
           ),
           (o.Units),
           Position=>Up},
	   Global=>false];
   ------import ideal----------------------------------------
   phi:=map(R,inputring,matrix{for i to varno-1-o.Units list X_i}|matrix{for i to o.Units-1 list U_i});
   newidealgens:=idealgens/phi;
   phiquot:=quot/phi;
   phiideal:=ideal(phiquot);
   phiinv:=map(inputring,R,
      matrix{{0,0}}|
      matrix{for i to idealsize-1 list 0}|
      matrix{for i to modulesize-1 list 0}|
      matrix{for i to numgens(inputring)-1 list inputring_i}
      );
   ----------------------------------------------------------
   --initialization------------------------------------------
   ----------------------------------------------------------
   List1:=List2:=List3:=List4:=List5:=List6:={};
   modulegens:=newmodulgens:=degreelist:=oldidealgens:=LMold:=LMnew:={};
   i:=j:=k:=0;
   deg1:=0;
   poly0:=poly1:=poly2:=poly3:=0;
   coef:=0;
   IG:=IB:=JB:=ideal(); 
   GB:=matrix{{0}};
   ----------------------------------------------------------
   --loop through increasing sequence of ideals--------------
   ----------------------------------------------------------
   while newidealgens!=oldidealgens do(
      oldidealgens=newidealgens;
      degreelist=for i to varno-1-o.Units list 
         max( for j to #newidealgens-1 list (degree(X_i,newidealgens#j)));
      if quot!={0} then(
         degreelist=for i to varno-1-o.Units list(
            deg1=degreelist#i;
            for j to #phiquot-1 do(
               deg1=max{deg1,(degree(X_i,phiquot#j))}
            );
            deg1
         );
      );
      List1={1_R}; 
      for i to varno-1-o.Units do(
         List1=flatten for j to degreelist#i list 
            for k to #List1-1 list(
	       poly0=(List1#k*X_i^j);
	       poly1=poly0%phiideal;
	       if leadMonomial(poly0)==leadMonomial(poly1) then poly1 else continue
         );	   
      );
      List1=reverse sort unique List1;
      --------------------------------------------------------------
      --reduction modulo the input ideal, I-------------------------
      --------------------------------------------------------------
      List2=for i to #List1-1 list {List1#i,((List1#i)^Q)%phiideal};
      IG=ideal(for i to #newidealgens-1 list T_0*newidealgens#i-T_1*G_i); 
      List3=for i to #List2-1 list 
         {(List2#i)#0,(List2#i)#1,(T_0*(List2#i)#1)%IG} ;
      --------------------------------------------------- 
      --loop through descending sequence of modules------
      ---------------------------------------------------
      modulegens={0_R};
      newmodulegens:={1_R};
      LMold={0};
      LMnew=for i to #newmodulegens-1 list leadMonomial(newmodulegens#i);
      while LMold!=LMnew do(
         modulegens=newmodulegens;
          LMold=LMnew;
         if (verbosity >= 1) then print(apply(modulegens, v->phiinv(v)));
         --reduction by module elements---------------------------
         IB=ideal(for i to #modulegens-1 list T_0*T_1*modulegens#i-T_1^2*B_i);
         JB=IB+T_0^(Q-1)*IG+ideal(T_1^Q)+phiideal+(ideal(T_0,T_1))^(Q+1);
         GB=gens gb JB;
         List4=reverse sort unique for i to #List3-1 list(
           {(List3#i)#0,(List3#i)#1,(T_0^(Q-1)*(List3#i)#2)%JB}
         );
         --row-reduction---------------------------
	 k=0;
	 while k<#List4 do(
            poly0=(List4#k)#0;
            poly1=(List4#k)#1;
            poly2=(List4#k)#2;
            j=0;
            while j<k do(
               if poly2!=0 
               and (List4#j)#2!=0 
               and leadMonomial(poly2)==leadMonomial((List4#j)#2) 
               then(
                  coef=leadTerm(poly2)//leadTerm((List4#j)#2);
                  poly0=poly0-coef*((List4#j)#0);
                  poly1=poly1-coef*((List4#j)#1);
                  poly2=(poly2-coef*((List4#j)#2))%JB;
                  j=0;
               )
               else(
                  j=j+1;
	       );--end if/else	  
            );--end while 
            List4=for i to # List4-1 list if i!=k then  List4#i else {poly0,poly1,poly2};
	    if poly2==0 then(
	       List4=for i to #List4-1 list(
		  if i>k and leadMonomial((List4#i)#0)%ideal(leadMonomial(poly0))==0 then continue else List4#i
	       );--end for/list
	    );--end if poly2
            k=k+1;
         );--end while k
         --extracting new module generators-------------------
         List6=for i to #(List4)-1 list 
            if (List4#i)#2==0 then ((List4#i)#0)%phiideal else continue;
         newmodulegens=reverse sort flatten entries gens gb (
            ideal(List6)+ideal(newidealgens)
         );
         newmodulegens=for i to #newmodulegens-1 list(
            poly0=newmodulegens#(#newmodulegens-1-i);
            for j to #newmodulegens-i-2 do( 
               if poly0!=0 then (
	          if (leadMonomial(poly0))%ideal(leadMonomial(newmodulegens#j))==0 then poly0=0;
	       );
	    );
	    poly0
	 );
         newmodulegens=reverse sort unique for i to #newmodulegens-1 list(
	       if newmodulegens#i!=0 then newmodulegens#i else continue
	       );
         LMnew=for i to #newmodulegens-1 list leadMonomial(newmodulegens#i);
      );
       newidealgens=newmodulegens
   );
   --------------------------------------------------------------
   --map back to input ring and output---------------------------
   --------------------------------------------------------------
   for i to #newidealgens-1 list phiinv(newidealgens#i)
);
-----------------------------------------------------------------------------
--unexported functions needed for reesalgebra--------------------------------
-----------------------------------------------------------------------------
 setcomplement = (m,l)-> (
    for i to #l-1 do(
       m=for j to #m-1 list 
 	 if leadMonomial(l#i)==leadMonomial(m#j) then continue else m#j
    );
    m
 );
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
 weightsFromMatrix = (m)->(
    apply(entries m , v->(Weights=>v))
 );
------------------------------------------------------------------------------------
 zeroMatrix = (a,b)->(
    matrix(for i to a-1 list(
	      for j to b-1 list 0)
	   )
 );
------------------------------------------------------------------------------------
 grevLex = (m)->( 
    matrix(for i to m-1 list(
	      for j to m-1 list(
		 if i+j<m then 1 else 0)
	      )
	  )
 );
------------------------------------------------------------------------------------
wtGrevlex = (matr)-> (
   grev := for i to numColumns(matr)-numRows(matr)-1 list (
              for j to numColumns(matr)-1 list(
	         if i+j < numColumns(matr)-numRows(matr) then 1 else 0
                 )
              );
              matr||matrix(grev)
          );
------------------------------------------------------------------------------------
--ring from weights and char--------------------------------------------------------
------------------------------------------------------------------------------------

reesring = method();
reesring(List,Ring) := (wt,kk)->(
------------------------------------------------------------------------------------
-- setup variable names and local weights from wt-----------------------------------
------------------------------------------------------------------------------------
symboldj := 0;
     mat := matrix(wt);    
     Mat := if numRows(mat) == numColumns(mat) then mat else wtGrevlex(mat);    
 symbold := (for i to numColumns(mat)-numRows(mat)-1 list "y")|
           (for i to numRows(mat)-1 list "x");
 symbold = for j to numColumns(mat)-1 list(
              symboldj=symbold#j;
              for i to numRows(mat)-1 do(
 	        symboldj=symboldj|toString(-mat_(i,j))
              );
              symboldj
           );
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
      r := kk[symbold,
             MonomialOrder=>{
                apply(entries Mat, v->(Weights=>v))},
             Global=>false];
      (mat,Mat,r)	 
 );
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
--extend the matrix of local weights to weights for the ideal-----------------------
------------------------------------------------------------------------------------
------------------------------------------------------------------------
qthReesAlgebra = method(TypicalValue => Matrix,
     Options=>{
	  Exponent => 1,Idealsize => 30,Modulesize => 40,Units => 0,Verb => 0
	  });
-------------------------------------------------------------
qthReesAlgebra(List,PolynomialRing,List,List):=
             o->(wt,R,J,I)->(
         rees := reesring(wt,coefficientRing(R));
         mat := rees#0;
         Mat := rees#1;
	 rng := rees#2;
    symbnewj := 0;
          kk := coefficientRing(rng);
------------------------------------------------------------------------------------    
 matnew := mat*matrix(for i to #gens(R)-1 list(
                       for j to #I-1 list(
 	                 degree(R_i,leadMonomial(I#j))
 	              )
                    )
             );
------------------------------------------------------------------------------------	 
 Matnew := if numRows(mat) == numColumns(mat) then(
     (Mat|matnew)||
     (zeroMatrix(#I,numColumns(Mat))|grevLex(#I))
     ) else(
     (Mat|(matnew||zeroMatrix(numRows(Mat)-numRows(matnew),#I)))||
     (zeroMatrix(#I,numColumns(Mat))|grevLex(#I)));
------------------------------------------------------------------------------------
--rees0-----------------------------------------------------------------------------
------------------------------------------------------------------------------------
symbnew := for i to #I-1 list "g1";
symbnew = gens(R)|for j to numColumns(matnew)-1 list(
              symbnewj = symbnew#j;
              for i to numRows(matnew)-1 do(
	         symbnewj = symbnewj|toString(-matnew_(i,j))
              );
              symbnewj
           );   
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
rees0ring := kk[symbnew,
                  MonomialOrder=>{
                       apply(entries Matnew, v->(Weights=>v))},
                  Global=>false];
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
-- set up local generator function variable t for rees(I)---------------------------
------------------------------------------------------------------------------------
--t:=getSymbol "t";
     rt := kk{t};
     use rt;
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
rees0quotientring := tensor(rt,
                      rees0ring/sub(ideal(J),rees0ring)
                    );           
use rees0quotientring;		
newreesideal := ideal(for i to #I-1 list 
                       sub(I#i,rees0quotientring)
                       -rees0quotientring_(#gens(R)+i+1)*t
                    );
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
time     ic1 := qthIdealClosure(R,J,I,wt);
       icold := ic1;
        Inew := I;
           k := 0;
 newreesring := rees0ring;
 Iold:=Inew;
 Matold:=Matnew; 
 symbold:=symbnew; 
 oldreesring:=newreesring; 
 oldreesideal:=newreesideal;
 newreesquotientring:=newreesring;
 Gnew:=matrix{{}};
 icnew:=icold;	   
------------------------------------------------------------------------------------
--loop from here to compute C(I^k) recursively until rees algebra is computed-------
------------------------------------------------------------------------------------
        Inew = setcomplement(icold,Inew); 
        while Inew != {} do( print(Inew);
                   k = k+1;
               Iold  = Inew;
             Matold  = Matnew; 
            symbold  = symbnew; 
        oldreesring  = newreesring;
       oldreesideal  = newreesideal; 
-----------------------------------------------------------------------------------
-----------------------------------------------------------------------------------
              matnew = mat*matrix(for i to #gens(R)-1 list(
                                     for j to #Inew-1 list(
	                                degree(R_i,Inew#j)
	                             )
                                  )
                       );
              Matnew = 
(Matold|(matnew||zeroMatrix(numRows(Matold)-numRows(matnew),numColumns(matnew))))||
	  (zeroMatrix(#Inew,numColumns(Matold))|grevLex(#Inew));
------------------------------------------------------------------------------------  
             symbnew = for i to #Inew-1 list "g"|k;
             symbnew = symbold|for j to numColumns(matnew)-1 list(
                          symbnewj = symbnew#j;
                          for i to numRows(matnew)-1 do(
	                     symbnewj = symbnewj|toString(-matnew_(i,j))
                          );
                          symbnewj
                       );   
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
       newreesring  = kk[symbnew,
                           MonomialOrder=>{
                              apply(entries Matnew, v->(Weights=>v))},
                           Global=>false];
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
newreesquotientring = tensor(rt,
                         newreesring/sub(ideal(J),newreesring)
                      );           
       newreesideal = sub(oldreesideal,newreesquotientring)+
                      ideal(for i to #Inew-1 list 
                         sub(Inew#i,newreesquotientring)
                         -newreesquotientring_(#gens(oldreesring)+i+1)*(sub(t,newreesquotientring))^k
                      );
               Gnew = gens gb newreesideal;
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
                Inew = flatten entries gens gb ((ideal(icold)*ideal(ic1)));
print(k);
time           icnew = qthIdealClosure(R,J,Inew,wt);
               icnew = for i to #icnew-1 list icnew#i%sub(ideal(J),R);
               icnew = flatten entries gens gb ideal icnew;
               icold = icnew;
                Inew = setcomplement(icold,Inew);		
         );
         Gnew
);
-------------------------------------------------------------------------------
----------------------  Documentation -----------------------------------------
-------------------------------------------------------------------------------
beginDocumentation()

document {
     Key => {
	  qthIdealClosure, (qthIdealClosure, Ring, List, List, List)
	  },
     Headline => "computes integral closures of ideals in positive characteristic, q",
     Usage => "(polynomials) = qthIdealClosure(inputring,quot,idealgens,wt)",
     Inputs => {
	  "inputring" => "inputring over ZZ/q with local ordering",
	  "quot" => "list of generators defining the ideal of relations of the affine domain",
	  "idealgens" => "set of generators for the ideal",
	  "wt"=> "the local weight list used above"
	  },
     Outputs => {
	  "polynomials" => List => "generating set for the closure"
	  },
     EXAMPLE {
		"q := 2;",
		"inputring1=(ZZ/q)[x,y,MonomialOrder=>{
			Weights=>{-1,-1},
			Weights=>{ 0,-1}},
				Global=>false];",
		"idealgens1={x^2+x^5,y^2+y^5};",
		"wt={{-1,-1},{0,-1}};",
		"qthIdealClosure(inputring1,{0},idealgens1,wt)"
	 },
     SeeAlso => {integralClosure},
	 Subnodes => {
        "Options",
		TO "Exponent",
		TO "Idealsize",
		TO "Modulesize",
		TO "Units",
		TO "Verb"
     },
    PARA { "The QthPower approach is to compute an increasing sequence of ideals
	   culminating with C(I,Q(A)), the closure of the input I 
	   as an ideal of the integral closure Q(A) of the input ring A.
	   Each ideal is last of a descending sequence of modules defined roughly by
	   M_{i+1}=< f in M_i : f^Q in M_i^{Q-1}I>,
	   which can be computed in positive char q by linear methods when Q is a power of q.
	   reesalgebra is an ideal which allows for determining n such that r in C(I^n,Q(A))
	   for any input r in Q(A).
	   See the paper by Leonard about this.
	   rationalIntegral Closure is a lifting of results over ZZ/q fro various q to
	   results over QQ."
    }
}
-------------------------------------------------------------------------------
document {
	Key => { Exponent },
	Headline => "runs the algorithm with Q = q^Exponent; default is 1",
	Usage => "qthIdealClosure(inputring, quot, idealgens, wt, Exponent => n)"
}
-------------------------------------------------------------------------------
document {
	Key => { Idealsize },
	Headline => "upper bound on number of ideal generators; default is 30",
	Usage => "qthIdealClosure(inputring, quot, idealgens, wt, Idealsize => n)"
}
-------------------------------------------------------------------------------
document {
	Key => { Modulesize },
	Headline => "upper bound on number of module generators; default is 40",
	Usage => "qthIdealClosure(inputring, quot, idealgens, wt, Modulesize => n)"
}
-------------------------------------------------------------------------------
document {
	Key => { Units },
	Headline => "declares the last n variables as units; default is 0",
	Usage => "qthIdealClosure(inputring, quot, idealgens, wt, Units => n)"
}
-------------------------------------------------------------------------------
document {
	Key => { Verb },
	Headline => "prints bases for intermediate modules if greater than 0; default is 0",
	Usage => "qthIdealClosure(inputring, quot, idealgens, wt, Verb => n)"
}
-------------------------------------------------------------------------------
document {
     Key => {
	  qthReesAlgebra, (qthReesAlgebra, List, PolynomialRing, List, List)
	  },
     Headline => "computes the non-homogeneous rees algebra in positive characteristic",
     Usage => "(polynomials) = reesalgebra(wt,R,J,I)",
     Inputs => {
	  "wt" => "list of weights used for local monomial ordering",
	  "R" => "ring with variables having local weight names from reesring",
	  "J" => "quotient for integrally closed A=R/J, {0} for polynomial ring",
	  "I"=> "ideal of A=R/J"
	  },
     Outputs => {
	  "polynomials" => Matrix => "gens gb rees algebra"
	  },
     EXAMPLE {
	       "wt := {{-9,-4,-4},{-3,-4,0}}",
               "kk := ZZ/5",
	       "R := (reesring(wt,kk))#2",
               "J := {y93^4+x44^3*x40^6}",
               "I := {x44,x40}",
               "rr := qthReesAlgebra(wt,R,J,I)",
               "(y93^2*x44*x40)%rr"
	 },
     SeeAlso => {reesAlgebra},
	 Subnodes => {
        "Options",
		TO "Exponent",
		TO "Idealsize",
		TO "Modulesize",
		TO "Units",
		TO "Verb"
    },
    PARA { "The integral closure of ideals paper contains 
	   a non-homogeneous version of the classical Rees algebra.
	   This has the map defining the ideal internal to it,
	   if {\tt r in R} is in {\tt C(I^d,A)}, 
	   it will have 1, not d images of r,
	   and the output can be used easily to determine the max d
	   such that {\tt r in R} is in {\tt C(I^d)}.
	   It calls qthIdealClosure once for each d, as long as
	   {\tt C(I^d)\neq C(I^(d-1)C(I)}."
    }
}
-------------------------------------------------------------------------------
--document {
--	Key => { Exponent },
--	Headline => "runs the algorithm with Q = q^Exponent; default is 1",
--	Usage => "reesalgebra(wt,R,J,I, Exponent => n)"
--}
-------------------------------------------------------------------------------
--document {
--	Key => { Idealsize },
--	Headline => "upper bound on number of ideal generators; default is 30",
--	Usage => "reesalgebra(wt,R,J,I, Idealsize => n)"
--}
-------------------------------------------------------------------------------
--document {
--	Key => { Modulesize },
--	Headline => "upper bound on number of module generators; default is 40",
--	Usage => "reesalgebra(wt,R,J,I, Modulesize => n)"
--}
-------------------------------------------------------------------------------
--document {
--	Key => { Units },
--	Headline => "declares the last n variables as units; default is 0",
--	Usage => "reesalgebra(wt,R,J,I, Units => n)"
--}
-------------------------------------------------------------------------------
--document {
--	Key => { Verb },
--	Headline => "prints bases for intermediate modules if greater than 0; default is 0",
--	Usage => "reesalgebra(wt,R,J,I, Verb => n)"
--}
----------------------------------------------------------------------------
document {
     Key => {
	  reesring, (reesring, List, Ring)
	  },
     Headline => "computes a polynomial ring with weights as given,
                  with variable names reflecting thos weights,
		  over the coeffecient ring given, with Global false",
     Usage => "(mat,Mat,R) = reesring(wt,kk)",
     Inputs => {
	  "wt" => "list of weights used for local monomial ordering",
	  "kk" => "ZZ/q or QQ"
	  },
     Outputs => {
--	  "mat" => Matrix => "wt as a matrix",
--	  "Mat" => "mat extended to give induced weights to the ideal generator names as well",
	  "R" => "with proper variable names and induced weights"
	  },
     EXAMPLE {
	       "wt := {{-9,-4,-4},{-3,-4,0}}",
               "kk := ZZ/5",
	       "(reesring(wt,kk))"
	 },
    PARA { "This sets up the initial rees algebra_0
	   for the reesalgebra method above."
         },
}
-----------------------------------------------------------------------------------
-----------------------------------------------------------------
--examples tested------------------------------------------------
-----------------------------------------------------------------  
TEST///
-----------------------------------------------------------------
--example 1 from paper by Leonard--------------------------------
--  this was the initial motivation for the local approach-------
-----------------------------------------------------------------
q=2;
inputring1=ZZ/q[x,y,MonomialOrder=>{
	Weights=>{-1,-1},
	Weights=>{ 0,-1}},
        Global=>false];
idealgens1={x^2+x^5,y^2+y^5};
wt={{-1,-1},{0,-1}};
IC1=qthIdealClosure(inputring1,{0},idealgens1,wt);
assert(IC1 == {x^2,x*y,y^2});
---------------------------------------------------------------
--example 2 ---------------------------------------------------
--from the same Integral closures of ideals paper--------------
--and this for ideals in quotient rings------------------------
---------------------------------------------------------------
q=2;
inputring2=ZZ/q[y90,x31,x20,MonomialOrder=>{
               Weights=>{-9,-3,-2},
               Weights=>{ 0,-1, 0},
               Weights=>{ 1, 0, 0}},
	       Global=>false];
quot2={y90^2+y90*x31^3+x20^9};
idealgens2={x31^2,x31*x20,x20^3};
wt2={{-9,-3,-2},
    {-0,-1, 0},
    { 1, 0, 0}};
IC2a=qthIdealClosure(inputring2,quot2,idealgens2,wt2,Idealsize=>60,Modulesize=>60);
IC2b=idealClosure(inputring2,quot2,flatten entries gens gb (ideal(IC2a))^2,wt2);
IC2c=idealClosure(inputring2,quot2,flatten entries gens gb (ideal(IC2b)*ideal(IC2a)),wt2);
assert(IC2a == {x31*x20,x20^3,x31^2,y90});
-----------------------------------------------------------------
--example 3------------------------------------------------------
--need for homogenization ---------------------------------------
--since 1-x^2,1-y^2 are local units------------------------------
--h^2-x^2,h^2-y^2 are not----------------------------------------
-----------------------------------------------------------------
q=3;
inputring3=ZZ/q[x,y,h,MonomialOrder=>{
	Weights=>{-1,-1,-1},
	Weights=>{ 0,-1,-1},
	Weights=>{ 0, 0,-1}},
        Global=>false];
idealgens3={y^2-h^2,x^2-h^2}
wt3={{-1,-1,-1},
     { 0,-1,-1},
     { 0, 0,-1}};
IC3=qthIdealClosure(inputring3,{0},idealgens3,wt3);
assert(IC3 == {x^2-y^2,y^2-h^2});
----------------------------------------------------------------
--example 4 with units------------------------------------------
----------------------------------------------------------------
--q=2;
--inputring4=ZZ/q[x,y,a,b,MonomialOrder=>{
--	  Weights=>{-1,-1,-1,-1},
--	   Weights=>{ 0,-1,-1,-1},
--	   Weights=>{ 0,0,-1,-1},
--	   Weights=>{0,0,0,-1}},
--	  Global=>false];     
--idealgens4={x^2+x*y^3*a,y^2+x^3*y*b};
--idealgens4={x^2*a,y^2*b}
--wt4={{-1,-1,-1,-1},
--    { 0,-1,-1,-1},{0,0,-1,-1},{0,0,0,-1}}
--IC4=qthIdealClosure(inputring4,{0},idealgens4,wt4,Units=>2);	       
----------------------------------------------------------------
--example 5----------------------------------------------------- 
----------------------------------------------------------------
q=5;
inputring5=ZZ/q[x,y,z,MonomialOrder=>{
	  Weights=>{-1,-1,-1},
	  Weights=>{ 0,-1,-1},
	  Weights=>{ 0, 0,-1}},
          Global=>false];
idealgens5={x^3,
           y^4,
           z^3};
wt5={{-1,-1,-1},
     { 0,-1,-1},
     { 0, 0,-1}};
IC5=qthIdealClosure(inputring5,{0},idealgens5,wt5);
assert(IC5 =={x^3,x^2*z,x*z^2,z^3,x^2*y^2,x*y^3,x*y^2*z,y^4,y^3*z,y^2*z^2});
----------------------------------------------------------------
--example 6, from ????------------------------------------------
----------------------------------------------------------------
q=2;
inputring6=ZZ/q[r,x,y,v,MonomialOrder=>{
	  Weights=>{-2,-2,-3,-1},
	  Weights=>{-1,-2,-3, 0},
	  Weights=>{ 1, 1, 0, 0},
	  Weights=>{ 1, 0, 0, 0}},
	  Global=>false];
quot6={r*x-y*v,x^5+x^2*v-r*y,r^2-x^4*v-x*v^2};
idealgens6={y^2,x*v^2+x^4*v};
wt6={{-2,-2,-3,-1},
     {-1,-2,-3, 0},
     { 1, 1, 0, 0},
     { 1, 0, 0, 0}};
IC6a=qthIdealClosure(inputring6,quot6,idealgens6,wt6);
IC6b=qthIdealClosure(inputring6,quot6,flatten entries gens gb (ideal(IC6a))^2,wt6);
assert(IC6a == {x*v^2+x^4*v, r*y,x*y*v + x^4*y,y^2});
-----------------------------------------------------------------
--example 7, from ???--------------------------------------------
-----------------------------------------------------------------
q=11;
inputring7=ZZ/q[x,y,MonomialOrder=>{
	Weights=>{-1,-1},
	Weights=>{0,-1}},
    Global=>false];
quot7={0};
idealgens7={x^5-x^3*y,
           x^4*y+y^3,
           x^3*y^2+x*y^3,
           x^2*y^3-y^4,
           x*y^4,
           y^5};
wt7={{-1,-1},
     { 0,-1}};
IC7=qthIdealClosure(inputring7,quot7,idealgens7,wt7,
     Idealsize=>6,Modulesize=>6);
assert(IC7 == {y^3, x^3*y-x^5, x^2*y^2, x^6});
-------------------------------------------------------
--example 8, used in my paper, but really from????-----
-------------------------------------------------------
q=2;
inputring8=ZZ/q[x,y,z];
idealgens8={x^5*y^4,
            y^5*z^4,
            z^5*x^4};
wt8={{-1,-1,-1},
     { 0,-1,-1},
     { 0, 0,-1}};
IC8=qthIdealClosure(inputring8,{0},idealgens8,wt8,Exponent=>2);
assert(IC8 == {x^5*y^4, x^4*y^4*z, x^4*y^3*z^2, x^4*y^2*z^3, x^4*y*z^4, x^4*z^5,
       x^3*y^4*z^2, x^3*y^3*z^3, x^3*y^2*z^4, x^2*y^4*z^3, x^2*y^3*z^4,
       x*y^4*z^4, y^5*z^4});
--------------------------------------------------------
--example 9, from ????----------------------------------
--------------------------------------------------------
q=3;
inputring9=ZZ/q[x,y,z,w];
idealgens9={x^2-y^2+x*z,
           x*y-y*z+x*w,
	   x*z-z^2+y*w,
	   x*w,
	   y^2-2*x*z, 
	   z^2*w-y*w,
	   y*z*w};
wt9={{-1,-1,-1,-1},
     { 0,-1,-1,-1},
     { 0, 0,-1,-1},
     { 0, 0, 0,-1}}
IC9=qthIdealClosure(inputring9,{0},idealgens9,wt9);
assert({IC9 == {x^2-x*z, x*y-y*z, x*z-z^2+y*w, x*w, y^2+z^2-y*w, y*w, z^2*w});
-------------------------------------------------------
--example 10  integrally closed, returns a local GB----
-------------------------------------------------------
q=2;
inputring10=ZZ/q[z,y,x,MonomialOrder=>{
	  Weights=>{-1,-1,-1},
	  Weights=>{ 0,-1,-1},
	  Weights=>{ 0, 0,-1}},
          Global=>false];
idealgens10= {x*z+x*y*z+y^2*z+x^3,x*y+x*y*z+x*z^2+y^3};
wt10={{-1,-1,-1},
      { 0,-1,-1},
      { 0, 0,-1}};
IC10=qthIdealClosure(inputring10,{0},idealgens10,wt10);
assert(IC10 == {z*x+z*y^2+z*y*x+x^3, 
	        y*x+z^2*x+z*y*x+y^3,
                z^3*y^2+z^2*y^3+z*y^4+y^3*x^2});
---------------------------------------------------------
--example 11 returns a local GB--------------------------
---------------------------------------------------------
q=5;
inputring11=ZZ/q[y,x,MonomialOrder=>{
	  Weights=>{-1,-1},
	  Weights=>{0,-1}},
          Global=>false];
idealgens11= {x^2+y^3+x^5,y^2*x+x^3+y^5};
wt11={{-1,-1},
      { 0,-1}};
IC11=qthIdealClosure(inputring11,{0},idealgens11,wt11);
assert(IC11 == {x^2+y^3, y^2*x, y^4});
---------------------------------------------------------
--example 12---------------------------------------------
---------------------------------------------------------
q=2;
inputring12=ZZ/q[x,y,MonomialOrder=>{
	Weights=>{-1,-1},
	Weights=>{0,-1}},
        Global=>false];
idealgens12={x^3+x^2*y^2,y^3+x^2*y^2}
wt12={{-1,-1},
      { 0,-1}}
IC12=qthIdealClosure(inputring12,{0},idealgens12,wt12);
assert(IC12 == {x^3, x^2*y, x*y^2, y^3});
----------------------------------------------------
--example 13 ---------------------------------------
----------------------------------------------------
q=2;
inputring13=ZZ/q[z5,z4,z3,z2,y,x2,x1,MonomialOrder=>{
	Weights=>{1,1,1,1,1,0,0},
	Weights=>{1,1,1,1,0,0,0},
	Weights=>{1,1,1,0,0,0,0},
	Weights=>{1,1,0,0,0,0,0},
	Weights=>{1,0,0,0,0,0,0},
	Weights=>{-1,-1,-1,-1,-1,-1,-1},
	Weights=>{ 0, 0, 0, 0, 0, 0,-1}},
    Global=>false];
quot13={y^6+y^3*x2^2*x1+x2^6,y^2+z2*x2,
     y*z2+z3*x2,
     y*(z3+x1)+z4*x2,
     y*z4+z5*x2,
     y*z5+x2^2}
idealgens13={y*x1,x2^2}
wt13={{ 1, 1, 1, 1, 1, 0, 0},
      { 1, 1, 1, 1, 0, 0, 0},
      { 1, 1, 1, 0, 0, 0, 0},
      { 1, 1, 0, 0, 0, 0, 0},
      { 1, 0, 0, 0, 0, 0, 0},
      {-1,-1,-1,-1,-1,-1,-1},
      { 0, 0, 0, 0, 0, 0,-1}}
IC13=qthIdealClosure(inputring13,quot13,idealgens13,wt13);
assert(IC13 == {z5*x2,z4*x2,z3*x2,z2*x2,y*x2,y*x1,x2^2});
----------------------------------------------------
--example 14----------------------------------------
----------------------------------------------------
q=2;
inputring14=ZZ/q[x,y,MonomialOrder=>{
	  Weights=>{-1,-1},
	  Weights=>{ 0,-1}},
	  Global=>false];     
idealgens14={x^2+y^3,y^5*(1+x*y)};
wt14={{-1,-1},
      { 0,-1}};
IC14=qthIdealClosure(inputring14,{0},idealgens14,wt14);	       
assert(IC14 == {x^2+y^3, x*y^4, y^5});
-----------------------------------------------------------
--example 15-----------------------------------------------
-----------------------------------------------------------
q=2;
inputring15=ZZ/q[x,y,a,b,MonomialOrder=>{
	  (Weights=>{-1,-1},
	   Weights=>{ 0,-1}),
	   (2)},
	  Global=>false];   
idealgens15={x^3+y^5,x^5*y+y^4};
wt15={{-1,-1},
      { 0,-1}};
IC15=qthIdealClosure(inputring15,{0},idealgens15,wt15,Units=>2);	       
assert(IC15 == {x^3, x^2*y^2, x*y^3, y^4});
---------------------------------------------------------
--example 16, integral closure of conductor--------------
---------------------------------------------------------
q=2;
inputring16=ZZ/q[y,x,MonomialOrder=>{
	  Weights=>{-1,-2},
	   Weights=>{0,-1}
	   },
	  Global=>false];
quot16={y^8+y^2*x^3+x^9};     
idealgens16={x^13};
wt16={{-1,-2},
      { 0,-1}};
IC16=qthIdealClosure(inputring16,quot16,idealgens16,wt16);	       
assert(IC16 == {y^7+y*x^3+y^6*x^3+x^6+y^5*x^6+y^4*x^9, 
	        y^6*x^10, 
		y^4*x^11, 
		y^2*x^12,
                    x^13});
------------------------------------------------------------
--example 17------------------------------------------------
------------------------------------------------------------
q=2;
inputring17=ZZ/q[y,x,MonomialOrder=>{
	  Weights=>{-1,-3},
	   Weights=>{0,-1}
	   },
	  Global=>false];
quot17={y^5+y^2*(x^4+x)+y*x^2+x^12};     
idealgens17={x^2};
wt17={{-1,-3},
      { 0,-1}};
IC17=qthIdealClosure(inputring17,quot17,idealgens17,wt17);	       
assert(IC17 == {y^4+y*x, y^3*x, x^2});
///
