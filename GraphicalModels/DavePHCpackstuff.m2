--Code
{* 
restart;
loadPackage("Graphs");
loadPackage("GraphicalModels");
loadPackage("GraphicalModelsMLE");
loadPackage "PHCpack"
G = mixedGraph(digraph {{1,{2}},{2,{3,4}}},bigraph {{1,3},{2,4}});
R=gaussianRing(G);
I=gaussianVanishingIdeal(R);
U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}};
M = sampleCovarianceMatrix(U);
J = scoreEquationsFromCovarianceMatrix(G,U);	
paramR = ring J;	
-- scoreEquations must return an ideal of R
-- also, the input must be R instead of G and it should not create again R
Inc = map (R,paramR);
use R;
J = Inc J;
Mparam = entries gaussianParametrization R;
L0 = Mparam#0;
Lu = flatten toList apply (1..#Mparam-1,i->drop(Mparam#i,{0,i-1})); 
Lparam = L0|Lu;
S = gaussianRing(4);
F = map(R,S,Lparam);
SJ = preimage (F,J);
NN=numgens S;
tempR = CC[zz_1 .. zz_NN];
tempF = map(tempR,S,vars tempR);
tempL = flatten entries gens tempF SJ;
solveSystem(tempL)

restart
loadPackage("Graphs");
loadPackage("GraphicalModels");
loadPackage("GraphicalModelsMLE");
loadPackage "PHCpack";
G = mixedGraph(digraph {{1,{2}},{2,{3,4}}},bigraph {{1,3},{2,4}});
U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}};
solveScoreEquationsInPHCpack = (G,U) -> (
R:=gaussianRing(G);    
J:=scoreEquationsFromCovarianceMatrix(G,U);    
paramR := ring J;	
-- scoreEquations must return an ideal of R
-- also, the input must be R instead of G and it should not create again R
Inc := map (R,paramR);
use R;
J = Inc J;
Mparam := entries gaussianParametrization R;
L0 := Mparam#0;
Lu := flatten toList apply (1..#Mparam-1,i->drop(Mparam#i,{0,i-1})); 
Lparam := L0|Lu;
S := gaussianRing(4);
F := map(R,S,Lparam);
SJ := preimage (F,J);
NN:=numgens S;
tempR := CC[local zz_1 .. local zz_NN];
tempF := map(tempR,S,vars tempR);
tempL := flatten entries gens tempF SJ;
solveSystem(tempL)
);
solveScoreEquationsInPHCpack(G,U)
*}

--Transcripts
Fordham-David-Swinarski:~ davids$ M2
Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases,
               PrimaryDecomposition, ReesAlgebra, TangentCone

i1 : restart;
Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases,
               PrimaryDecomposition, ReesAlgebra, TangentCone

i1 : loadPackage("Graphs");

i2 : loadPackage("GraphicalModels");

i3 : loadPackage("GraphicalModelsMLE");
--loading configuration for package "Bertini" from file /Users/davids/Library/Application Support/Macaulay2/init-Bertini.m2

i4 : loadPackage "PHCpack"
--loading configuration for package "PHCpack" from file /Users/davids/Library/Application Support/Macaulay2/init-PHCpack.m2
--warning: symbol "Verbose" in Core.Dictionary is shadowed by a symbol in PHCpack.Dictionary
--  use the synonym Core$Verbose

o4 = PHCpack

o4 : Package

i5 : G = mixedGraph(digraph {{1,{2}},{2,{3,4}}},bigraph {{1,3},{2,4}});

i6 : R=gaussianRing(G);

i7 : I=gaussianVanishingIdeal(R);

o7 : Ideal of R

i8 : U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}};

i9 : M = sampleCovarianceMatrix(U);

              4        4
o9 : Matrix QQ  <--- QQ

i10 : J = scoreEquationsFromCovarianceMatrix(G,U);

o10 : Ideal of QQ[l   , l   , l   , p   , p   , p   , p   , p   , p   ]
                   1,2   2,4   2,3   1,1   2,2   3,3   4,4   1,3   2,4

i11 : paramR = ring J;

i12 : -- scoreEquations must return an ideal of R
      -- also, the input must be R instead of G and it should not create again R
      Inc = map (R,paramR);

o12 : RingMap R <--- paramR

i13 : use R;

i14 : J = Inc J;

o14 : Ideal of R

i15 : Mparam = entries gaussianParametrization R;

i16 : L0 = Mparam#0;

i17 : Lu = flatten toList apply (1..#Mparam-1,i->drop(Mparam#i,{0,i-1})); 

i18 : Lparam = L0|Lu;

i19 : S = gaussianRing(4);

i20 : F = map(R,S,Lparam);

o20 : RingMap R <--- S

i21 : SJ = preimage (F,J);

o21 : Ideal of S

i22 : NN=numgens S;

i23 : tempR = CC[zz_1 .. zz_NN];

i24 : tempF = map(tempR,S,vars tempR);

o24 : RingMap tempR <--- S

i25 : tempL = flatten entries gens tempF SJ;

i26 : solveSystem(tempL)
using temporary files /var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-26102-0/0PHCinput and /var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-26102-0/0PHCoutput

o26 = {{7.1875, -1.625, -1.8125, 5.875, 1.25, .875, -2.75, 1.6875, -2.23091,
      --------------------------------------------------------------------------
      7.25}, {7.1875, 9.85141e-16, -.955973, -5.54708e-15, 1.25, .658867, -2.75,
      --------------------------------------------------------------------------
      1.45966, -1.44951, 7.25}}

o26 : List

i27 : restart
Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases,
               PrimaryDecomposition, ReesAlgebra, TangentCone

i1 : restart
Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases,
               PrimaryDecomposition, ReesAlgebra, TangentCone

i1 : loadPackage("Graphs");

i2 : loadPackage("GraphicalModels");

i3 : loadPackage("GraphicalModelsMLE");
--loading configuration for package "Bertini" from file /Users/davids/Library/Application Support/Macaulay2/init-Bertini.m2

i4 : loadPackage "PHCpack";
--loading configuration for package "PHCpack" from file /Users/davids/Library/Application Support/Macaulay2/init-PHCpack.m2
--warning: symbol "Verbose" in Core.Dictionary is shadowed by a symbol in PHCpack.Dictionary
--  use the synonym Core$Verbose

i5 : G = mixedGraph(digraph {{1,{2}},{2,{3,4}}},bigraph {{1,3},{2,4}});

i6 : U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}};

i7 : solveScoreEquationsInPHCpack = (G,U) -> (
     R:=gaussianRing(G);    
     J:=scoreEquationsFromCovarianceMatrix(G,U);    
     paramR := ring J;
     -- scoreEquations must return an ideal of R
     -- also, the input must be R instead of G and it should not create again R
     Inc := map (R,paramR);
     use R;
     J = Inc J;
     Mparam := entries gaussianParametrization R;
     L0 := Mparam#0;
     Lu := flatten toList apply (1..#Mparam-1,i->drop(Mparam#i,{0,i-1})); 
     Lparam := L0|Lu;
     S := gaussianRing(4);
     F := map(R,S,Lparam);
     SJ := preimage (F,J);
     NN:=numgens S;
     zz := local zz;
     tempR := CC[zz_1..zz_NN];
     tempF := map(tempR,S,vars tempR);
     tempL := flatten entries gens tempF SJ;
     solveSystem(tempL)
     );

i8 : solveScoreEquationsInPHCpack(G,U)
warning: clearing value of symbol zz to allow access to subscripted variables based on it
       : debug with expression   debug 6274   or with command line option   --debug 6274
using temporary files /var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-26102-0/0PHCinput and /var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-26102-0/0PHCoutput
stdio:27:1:(3):[1]: error: key not found in hash table