restart;
loadPackage("Graphs");
loadPackage("GraphicalModels");
loadPackage("Bertini");
loadPackage("GraphicalModelsMLE");
loadPackage "PHCpack"

G = mixedGraph(digraph {{1,{2}},{2,{3,4}}},bigraph {{1,3},{2,4}});
R=gaussianRing(G);
I=gaussianVanishingIdeal(R)
U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}};
M = sampleCovarianceMatrix(U)	
J = scoreEquationsFromCovarianceMatrix(G,U);	
paramR = ring J	
-- scoreEquations must return an ideal of R
-- also, the input must be R instead of G and it should not create again R
Inc = map (R,paramR)
use R
J = Inc J
Mparam = entries gaussianParametrization R
L0 = Mparam#0
Lu = flatten toList apply (1..#Mparam-1,i->drop(Mparam#i,{0,i-1})) 
Lparam = L0|Lu
S = gaussianRing(4);
F = map(R,S,Lparam)
SJ = preimage (F,J)

NN=numgens S;
tempR = CC[zz_1 .. zz_NN];
tempF = map(tempR,S,vars tempR);
tempL = flatten entries gens tempF SJ

solveSystem(tempL)

bertiniZeroDimSolve(tempL)

NN:=numgens S;
tempR := QQ[local zz_1 .. local zz_NN];
tempF := map(tempR,S,vars tempR);
tempL := flatten entries gens tempF SJ
bertiniZeroDimSolve(tempL)

----------------------------------------

CS = CC monoid S
CF = map (CS,S, vars CS)
SJ = CF SJ
CL = flatten entries gens SJ
solveSystem CL
