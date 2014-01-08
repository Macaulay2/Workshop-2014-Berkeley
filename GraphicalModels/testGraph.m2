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

sampleCovarianceMatrix = (U) -> (
    n = #U;
    Ubar = matrix{{(1/n)}} * sum(U);
    return ((1/n)*(sum apply(n, i -> (transpose (U#i-Ubar))*(U#i-Ubar))));        
);

MLEmixedGraph = (G, U) -> (
    V = sampleCovarianceMatrix(U);
    R = gaussianRing G;
    I = gaussianVanishingIdeal R;
    -- Lambda
    L = directedEdgesMatrix R;
    -- d is equal to the number of vertices in G
    d = numgens source L;
    -- Omega
    W = bidirectedEdgesMatrix R;
    K = inverse (id_(R^d)-L);
    S = (transpose K) * W * K;
    FR = frac(R);
    Sinv = inverse substitute(S, FR);
    C1 = trace(Sinv * V)/2;
    f = numerator(C1);
    g = denominator(C1);
    CrazyExp = JacobianMatrixOfRationalFunction(trace(Sinv * V)/2);
    LL = (substitute((jacobian(matrix{{det(S)}})), FR))*matrix{{(-(#U)/(2*det(S)))}} - (transpose CrazyExp);
    I1 = ideal(0);
    for i to (numgens target LL)-1 do (I1 = I1 + ideal(numerator(LL_(i, 0))));
    --Sigma = matrix{{s_(1,1), s_(1,2), s_(1,3)}, {s_(1,2), s_(2,2), s_(2,3)}, {s_(1,3), s_(2,3), s_(3,3)}}
    use R;
    Sigma = matrix{{s_(1,1), s_(1,2), s_(1,3), s_(1,4)}, {s_(1,2), s_(2,2), s_(2,3), s_(2,4)}, {s_(1,3), s_(2,3), s_(3,3), s_(3,4)}, 
        {s_(1,4), s_(2,4), s_(3,4), s_(4,4)}};
    J = ideal(Sigma - S);
    print I1;
    print J;
    print ring(I1);
    print ring(J);
    I1 = I1 + J;
    --I = I:det(S)
    I1 = saturate(I1, det(S));
    print dim I1;
    degree I1;
    if (dim I1 == 0) then (
        print concatenate("The dimension of the ideal of likelihood equations is 0 and its degree is ", toString (degree I1));
    ) else (
        print concatenate("The dimension of the ideal of likelihood equations is", toString (dim I1),
	    "and its degree is", toString (degree I1), ". The model is not identifiable.");
    );
    I1 == radical I1;
    mingens I1;
);

MLEmixedGraph(G, U)
