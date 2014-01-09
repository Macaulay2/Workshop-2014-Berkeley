restart;
loadPackage("Graphs");
loadPackage("GraphicalModels");
loadPackage("Bertini");
loadPackage("GraphicalModelsMLE");
G = mixedGraph(digraph {{1,{2}},{2,{3,4}}},bigraph {{1,3},{2,4}});
R=gaussianRing(G);
I=gaussianVanishingIdeal(R)
U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}};
J = scoreEquationsFromCovarianceMatrix(G,U);	
	
-- scoreEquations must return an ideal of R
-- also, the input must be R instead of G and it should not create again R

use R
J =  ideal(12992*p_(1,3)+3105,1495*p_(4,4)-1102*p_(2,4)-897,126063*p_(3,3)+146687*p_(1,3),359950*p_(2,2)+50531*p_(2,4)-101774,27*p_(1,1)+203*p_(1,3
      ),203*l_(2,3)-107,5*l_(2,4)+16*p_(2,4)+11,179975*l_(1,2)+62192*p_(2,4)+13182,9568*p_(2,4)^2-2204*p_(2,4)-897,12992*p_(1,3)*p_(2,4)+3105*p_(2,4),
      168792064*p_(1,3)^2-9641025)

Mparam = entries gaussianParametrization R
L0 = Mparam#0
Lu = flatten toList apply (1..#Mparam-1,i->drop(Mparam#i,{0,i-1})) 
Lparam = L0|Lu
S = gaussianRing(4);
F = map(R,S,Lparam)
SJ = preimage (F,J)

NN=numgens S;
tempR = QQ[zz_1 .. zz_NN];
tempF = map(tempR,S,vars tempR);
tempL = flatten entries gens tempF SJ
bertiniZeroDimSolve(tempL)

NN:=numgens S;
tempR := QQ[local zz_1 .. local zz_NN];
tempF := map(tempR,S,vars tempR);
tempL := flatten entries gens tempF SJ
bertiniZeroDimSolve(tempL)
