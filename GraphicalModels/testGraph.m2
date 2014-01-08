restart
needsPackage "Bertini"
needsPackage "GraphicalModels"
load "sharedFunctions.m2"

--G = mixedGraph(digraph {{1,{2,3}},{2,{3}},{3,{4}}},bigraph {{1,2},{2,4}});
G = mixedGraph(digraph {{1,{2,3}},{2,{3}},{3,{4}}},bigraph {{3,4}});
--G = mixedGraph(digraph {{1, {2}}, {2, {3}}}, bigraph {{2, 3}});
-- number of vertices in G
--d = 4;
-- list whose rows are the sampels
--U = {matrix{{1, 2, -10}}, matrix{{-1, -20, 5}}, matrix{{3, -50, 2}}, matrix{{-1, -4, -30}}};
U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}

matRtolpR = (M) -> (
E:=entries(M);    
return matrix apply(#E, i -> apply(#(E_i), j -> F(E_i_j)))    
);

matZZtoQQ = (M) -> (
E:=entries(M);    
return matrix apply(#E, i -> apply(#(E_i), j -> (1/1)*E_i_j))    
);

sampleCovarianceMatrix = (U) -> (
    n = #U;
    U=apply(#U, i -> if ring(U_i)===ZZ then matZZtoQQ(U_i) else U_i);
    Ubar = matrix{{(1/n)}} * sum(U);
    return ((1/n)*(sum apply(n, i -> (transpose (U#i-Ubar))*(U#i-Ubar))));        
);

JacobianMatrixOfRationalFunction = (F) -> (
f:=numerator(F);
g:=denominator(F);
R:=ring(f);
answer:=diff(vars(R), f) * g - diff(vars(R), g)*f;
answer=substitute(answer, FR);
matrix({{(1/g)^2}})*answer
);

MLEmixedGraph = (G, U) -> (
    V = sampleCovarianceMatrix(U);
    R = gaussianRing G;
    I = gaussianVanishingIdeal R;
    -- Lambda
    L = directedEdgesMatrix R;
    -- d is equal to the number of vertices in G
    d = numRows L;
    -- Omega
    W = bidirectedEdgesMatrix R;
    -- move to a new ring, lpR, which does not have the s variables
numSvars=lift(d*(d+1)/2,ZZ);
--lp rings is the ring without the s variables
lpRvarlist=apply(numgens(R)-numSvars,i->(gens(R))_i);
KK=coefficientRing(R);
lpR=KK[lpRvarlist];
lpRTarget=apply(numgens(R),i-> if i<= numgens(R)-numSvars-1 then (gens(lpR))_i else 0);
F=map(lpR,R,lpRTarget);
I=F(I);    
    L = matRtolpR L;
    W = matRtolpR W;
    K = inverse (id_(lpR^d)-L);
    S = (transpose K) * W * K;
    FR = frac(lpR);
    Sinv = inverse substitute(S, FR);
    C1 = trace(Sinv * V)/2;
    f = numerator(C1);
    g = denominator(C1);
    CrazyExp = JacobianMatrixOfRationalFunction(trace(Sinv * V)/2);
    LL = (substitute((jacobian(matrix{{det(S)}})), FR))*matrix{{(-(#U)/(2*det(S)))}} - (transpose CrazyExp);
--    I1 = ideal(0);
--    for i to (numgens target LL)-1 do (I1 = I1 + ideal(numerator(LL_(i, 0))));
--    --Sigma = matrix{{s_(1,1), s_(1,2), s_(1,3)}, {s_(1,2), s_(2,2), s_(2,3)}, {s_(1,3), s_(2,3), s_(3,3)}}
--    use R;
--    Sigma = matrix{{s_(1,1), s_(1,2), s_(1,3), s_(1,4)}, {s_(1,2), s_(2,2), s_(2,3), s_(2,4)}, {s_(1,3), s_(2,3), s_(3,3), s_(3,4)}, 
--        {s_(1,4), s_(2,4), s_(3,4), s_(4,4)}};
--    J = ideal(Sigma - S);
ELL=flatten entries(LL);
J=ideal apply(#ELL, i -> lift(numerator(ELL_i),lpR));
J = saturate(J, det(S));
--    print I1;
    print toString(J) << endl;
--    print ring(I1);
    print toString(ring(J)) << endl;
--    I1 = I1 + J;
--    --I = I:det(S)
--    I1 = saturate(I1, det(S));
    print toString(dim J) << endl;
    degree J;
    if (dim J == 0) then (
        print concatenate("The dimension of the ideal of likelihood equations is 0 and its degree is ", toString (degree J));
    ) else (
        print concatenate("The dimension of the ideal of likelihood equations is", toString (dim J),
	    "and its degree is", toString (degree J), ". The model is not identifiable.");
    );
    print toString(J == radical J) << endl;
    print toString(mingens J) << endl;
);

MLEmixedGraph(G, U)
