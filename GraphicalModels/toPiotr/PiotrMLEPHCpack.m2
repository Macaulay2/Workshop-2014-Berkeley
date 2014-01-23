loadPackage("Graphs");
loadPackage("GraphicalModels");
--loadPackage("Bertini");
loadPackage("PHCpack");
loadPackage("GraphicalModelsMLE");
debug GraphicalModelsMLE;

scoreEquationsFromCovarianceMatrixVerbose = (G, U) -> (
    R:=gaussianRing(G);
    S := sampleCovarianceMatrix(U);    
    print concatenate("S=",toString S) << endl;
    use R;   --This shouldn't be necessary once we fix a bug in GraphicalModels 
    Lambda := directedEdgesMatrix R;   
    print concatenate("Lambda=",toString Lambda) << endl;
    -- d is equal to the number of vertices in G
    d := numRows Lambda;
    Omega := bidirectedEdgesMatrix R;
    print concatenate("Omega=",toString Omega) << endl;    
    -- move to a new ring, lpR, which does not have the s variables
    numSvars:=lift(d*(d+1)/2,ZZ);
    --lp rings is the ring without the s variables
    lpRvarlist:=apply(numgens(R)-numSvars,i->(gens(R))_i);
    KK:=coefficientRing(R);
    lpR:=KK[lpRvarlist];
    lpRTarget:=apply(numgens(R),i-> if i<= numgens(R)-numSvars-1 then (gens(lpR))_i else 0);
    F:=map(lpR,R,lpRTarget);
    Lambda = matRtolpR(Lambda,F);
    Omega = matRtolpR(Omega,F);
    FR := frac(lpR);
    A := inverse (id_(lpR^d)-Lambda);  --What is the ring of A?
    Sigma := (transpose A) * Omega * A;  
    Dpoly:=det Sigma;
    Sigma=substitute(Sigma,FR);
    print concatenate("Sigma=",toString Sigma) << endl;      
    D:=det Sigma;
    print concatenate("D=",toString D) << endl;      
    n:=#U;
    Sigmainv := Sigma^-1;  
    print concatenate("Sigmainv=",toString Sigmainv) << endl;  
    logLikelihoodFirstTermDerivative := (transpose JacobianMatrixOfRationalFunction(D))*matrix{{-n/(2*D)}};
    print concatenate("logLikelihoodFirstTermDerivative=",toString logLikelihoodFirstTermDerivative) << endl;  		
    logLikelihoodSecondTerm:=(n/2)*trace(S*Sigmainv);
    print concatenate("Trace term =",toString(logLikelihoodSecondTerm))<< endl;
    logLikelihoodSecondTermDerivative := transpose JacobianMatrixOfRationalFunction(logLikelihoodSecondTerm);   
    print concatenate("logLikelihoodSecondTermDerivative=",toString logLikelihoodSecondTermDerivative) << endl;    
    logLikelihoodDerivative := logLikelihoodFirstTermDerivative - logLikelihoodSecondTermDerivative;  
    LLD:=flatten entries(logLikelihoodDerivative);
    J:=ideal apply(#LLD, i -> lift(numerator(LLD_i),lpR));   
    J = saturate(J, Dpoly);
    return J
);


{*
solveScoreEquationsInBertini = (G,U) -> (
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
bertiniZeroDimSolve(tempL)
)
*}

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
)

