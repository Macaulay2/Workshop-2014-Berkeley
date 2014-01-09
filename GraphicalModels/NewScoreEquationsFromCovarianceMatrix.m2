NewScoreEquationsFromCovarianceMatrix = (G, U) -> (
    S := sampleCovarianceMatrix(U);   
    R := gaussianRing(G); 
    use R;   
    -- Lambda
    Lambda := directedEdgesMatrix R;   
    -- d is equal to the number of vertices in G
    d := numRows L;
    -- Omega
    Omega := bidirectedEdgesMatrix R;
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
    A := inverse (id_(lpR^d)-Lambda);
    Sigma := (transpose A) * Omega * A;
    Dpoly:=det Sigma;
    Sigma=substitute(Sigma,FR);
    print toString(S) << endl;
    D:=det Sigma;
    n:=#U;
    Sigmainv := Sigma^-1;  
    logLikelihoodFirstTermDerivative := transpose JacobianMatrixOfRationalFunction(D)*matrix{{-n/(2*D)}};	
    logLikelihoodSecondTerm:=(n/2)*trace(S*Sigmainv);
    logLikelihoodSecondTermDerivative := transpose JacobianMatrixOfRationalFunction(logLikelihoodSecondTerm);   
    logLikelihoodDerivative := logLikelihoodFirstTermDerivative - logLikelihoodSecondTermDerivative;  
    LLD:=flatten entries(logLikelihoodDerivative);
    J:=ideal apply(#LLD, i -> lift(numerator(LLD_i),lpR));   
    J = saturate(J, Dpoly);
    return J
);