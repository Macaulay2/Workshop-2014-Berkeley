
p. 38

Prop 2.1.12 p. 39 has the version of the log likelihood function



NewScoreEquationsFromCovarianceMatrix = (R, U) -> (
    S := sampleCovarianceMatrix(U);    
    use R;   --This shouldn't be necessary once we fix a bug in GraphicalModels 
    Lambda := directedEdgesMatrix R;   
    -- d is equal to the number of vertices in G
    d := numRows Lambda;
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
    A := inverse (id_(lpR^d)-Lambda);  --What is the ring of A?
    Sigma := (transpose A) * Omega * A;
    Dpoly:=det Sigma;
    Sigma=substitute(Sigma,FR);
    D:=det Sigma;
    n:=#U;
    Sigmainv := Sigma^-1;  
    logLikelihoodFirstTermDerivative := (transpose JacobianMatrixOfRationalFunction(D))*matrix{{-n/(2*D)}};	
    logLikelihoodSecondTerm:=(n/2)*trace(S*Sigmainv);
    logLikelihoodSecondTermDerivative := transpose JacobianMatrixOfRationalFunction(logLikelihoodSecondTerm);   
    logLikelihoodDerivative := logLikelihoodFirstTermDerivative - logLikelihoodSecondTermDerivative;  
    LLD:=flatten entries(logLikelihoodDerivative);
    J:=ideal apply(#LLD, i -> lift(numerator(LLD_i),lpR));   
    J = saturate(J, Dpoly);
    return J
);

--------------------------------------------
-- turn the entries of an integer matrix into rational numbers
--------------------------------------------

matZZtoQQ = (M) -> (
    E:=entries(M);    
    return matrix apply(#E, i -> apply(#(E_i), j -> (1/1)*E_i_j))    
);

--------------------------------------------
-- change the entries of matrix M from ring R to ring lpR 
-- via the map F:R-->lpR
--------------------------------------------

matRtolpR = (M,F) -> (
    E:=entries(M);    
    return matrix apply(#E, i -> apply(#(E_i), j -> F(E_i_j)))    
);

sampleCovarianceMatrix = method();
sampleCovarianceMatrix(List) := (U) -> (
    n := #U;
    U = apply(#U, i -> if ring(U_i)===ZZ then matZZtoQQ(U_i) else U_i);
    Ubar := matrix{{(1/n)}} * sum(U);
    return ((1/n)*(sum apply(n, i -> (transpose (U#i-Ubar))*(U#i-Ubar))));        
);


JacobianMatrixOfRationalFunction = method();
JacobianMatrixOfRationalFunction(RingElement) := (F) -> (
    f:=numerator(F);
    g:=denominator(F);
    R:=ring(f);
    answer:=diff(vars(R), f) * g - diff(vars(R), g)*f;
    answer=substitute(answer, ring(F));
    matrix({{(1/g)^2}})*answer
);

-----------------------------------
-- This function takes a mixed graph and a permutation and
-- creates the same graph but with permuted vertices.
-----------------------------------
permuteVerticesOfGraph = method();
permuteVerticesOfGraph(MixedGraph, List) := (G, perm) -> (
    perm1 = apply(#perm, i -> (perm#i +1));
    edgesUn = edges(Un);
    edgesBi = edges(Bi);
    edgesDi = edges(Di);
    newEdgesUn = {};
    for i to (#edgesUn-1) do (
	ed = toList(edgesUn#i);
	ed1 = apply(#ed, i -> (perm1#(ed#i-1)));
	newEdgesUn = append(newEdgesUn, ed1);
    );
    newEdgesBi = {};
    for i to (#edgesUn-1) do (
        ed = toList(edgesBi#i);
        ed1 = apply(#ed, i -> (perm1#(ed#i-1)));
        newEdgesBi = append(newEdgesBi, ed1);
    );
    newEdgesDi = {};
    for i to (#edgesDi-1) do (
       ed = edgesDi#i;
       ed1 = apply(#ed, i -> (perm1#(ed#i-1)));
       newEdgesDi = append(newEdgesDi, ed1);
    );
    newG = mixedGraph(graph(newEdgesUn), digraph(newEdgesDi), bigraph(newEdgesBi));
    return newG;
);

------------------------------------------
-- This function separates the vertices into two groups: the undirected
-- and the bidirected subgraphs. If a vertex does not have any adjacent
-- undirected or bidirected edge, then, it is distributed to one of the
-- two sets in a way to not create collisions.
-----------------------------------------
orderingOfMixedGraphUndirectedBidirectedSeparation = method();
orderingOfMixedGraphUndirectedBidirectedSeparation(MixedGraph) := (G) -> (
    Un = G#graph#Graph;	   
    Bi = G#graph#Bigraph;    
    Di = G#graph#Digraph;
    vertUn = vertices(Un);
    vertBi = vertices(Bi);
    vertDi = vertices(Di);
    vertRest = toList(set(vertDi) - (set(vertUn)+set(vertBi))); 
    for i to #vertRest-1 do (
	flag = true;
	for j to #(edges(Di))-1 do (
	    ed = (edges(Di))#j;
	    if (ed#1 == vertRest#i and isSubset(set{ed#0}, set(vertBi))) then (
		print(vertRest#i);
		print(ed);
		print("vertRest#i = ed#1 and ed#0 is in vertBi");
		vertBi = append(vertBi, ed#1);
		flag = false;
		break;
	    );
	);
    	if (flag == true) then (
	    vertUn = append(vertUn, vertRest#i);
	);
    );    
    if (#(set(vertBi)*set(vertUn)) > 0) then (
	print("The bidirected and undirected part cannot have any vertices in common!");
	return FALSE;
    );
    vertUn = sort(vertUn);
    vertBi = sort(vertBi);
    numUn = #vertUn;
    numBi = #vertBi;
    perm = vertUn | vertBi;
    perm1 = apply(#perm, i -> (perm#i-1));
    perm2 = inversePermutation perm1;
    newG = permuteVerticesOfGraph(G, perm2);
    return {newG, perm2};    
);

--------------------------------------------
-- This function rearranges the vertices withing the undirected and bidirected
-- part so that all directed edges go from a smaller vertex to a bigger vertex.
--------------------------------------------
orderingOfMixedGraphDirectedOrdering = method();
orderingOfMixedGraphDirectedOrdering(MixedGraph) := (G) -> (    
    Un = G#graph#Graph;	   
    Bi = G#graph#Bigraph;    
    Di = G#graph#Digraph;
    vertUn = vertices(Un);
    vertBi = vertices(Bi);
    vertUn = sort(vertUn);
    vertBi = sort(vertBi);
    vertRest = sort(toList(set(vertDi) - (set(vertUn)+set(vertBi)))); 
    for i to #vertRest-1 do (
	flag = true;
	for j to #(edges(Di))-1 do (
	    ed = (edges(Di))#j;
	    if (ed#1 == vertRest#i and isSubset(set{ed#0}, set(vertBi))) then (
		print(vertRest#i);
		print(ed);
		print("vertRest#i = ed#1 and ed#0 is in vertBi");
		vertBi = append(vertBi, ed#1);
		flag = false;
		break;
	    );
	);
    	if (flag == true) then (
	    vertUn = append(vertUn, vertRest#i);
	);
    );
    numUn = #vertUn;
    numBi = #vertBi;
    if (numUn + numBi < #(vertices(G))) then (
    	print ("Problem!");
    );
    -- arrange the directed edges going from small to large in the undirected
    -- part of the graph.
    UnM = mutableMatrix(ZZ, 1, numUn);
    for i to (#edges(Di)-1) do (
	ed = (edges(Di))#i;
	print(ed);
	if (ed#0 <= numUn and ed#1 <= numUn) then (
    	    UnM_(0, ed#1-1) = UnM_(0, ed#1-1) + 1;
	);
    );
    S = set{};
    for i to numUn-1 do (
    	S = S + set{{i+1, UnM_(0, i)}};
    );
    perm = {};
    -- this is not as efficient as it could be; for example we can use a heap
    numS = #S;
    for s to numS-1 do (
	el = elements(S);
    	for i to #el-1 do (
	    if UnM_(0,el#i#0-1) == 0 then (
		perm = append(perm, el#i#0-1);
		ch = children(Di, el#i#0);
		elch = elements(ch);
		for j to #elch-1 do (
		    if(elch#j <= numUn) then (
		    	UnM_(0, elch#j-1) = UnM_(0, elch#j-1) - 1;
		    );
	        );
	    	S = S - set{el#i};
		break;
	    );
    	);	
    );
    -- rearrange vertices so that directed edges go from small to large in the bidirected
    -- part
    BiM = mutableMatrix(ZZ, 1, numBi);
    for i to (#edges(Di)-1) do (
	ed = (edges(Di))#i;
	print(ed);
	if (ed#0 <= numBi+numUn and numUn < ed#0 and ed#1 <= numBi+numUn and numUn < ed#1) then (
    	    BiM_(0, ed#1-numUn-1) = BiM_(0, ed#1-numUn-1) + 1;
	);
    );
    S = set{};
    for i to numBi-1 do (
    	S = S + set{{i+numUn+1, BiM_(0, i)}};
    );
    print (S);
    -- this is not as efficient as it could be; for example we can use a heap
    numS = #S;
    for s to numS-1 do (
	el = elements(S);
    	for i to #el-1 do (
	    if BiM_(0,el#i#0-numUn-1) == 0 then (
		perm = append(perm, el#i#0-1);
		ch = children(Di, el#i#0);
		elch = elements(ch);
		for j to #elch-1 do (
		    if(elch#j <= numBi + numUn and numUn < elch#j) then (
		    	BiM_(0, elch#j-numUn-1) = BiM_(0, elch#j-numUn-1) - 1;
		    );
	        );
	    	S = S - set{el#i};
		break;
	    );
    	);	
    );
    perm1 = inversePermutation perm;
    newG = permuteVerticesOfGraph(G, perm1);
    return {newG, perm1};    
);

-------------------------------
-- This function takes a mixed graph and permutes the vertices so that the
-- requirements on p.9 from Sullivant, Talaska and Draisma's paper, called
-- Trek Separation for Gaussian Graphical Models, are satisfied.
------------------------------
orderingOfMixedGraph = method();
orderingOfMixedGraph(MixedGraph) := (G) -> (
    l1 = orderingOfMixedGraphUndirectedBidirectedSeparation(G);
    l2 = orderingOfMixedGraphDirectedOrdering(l1#0);
    perm = {};
    for i to (#l1-1) do (
	perm = append(perm, l2#1#(l1#1#i - 1));
    );
    return {l2#0, perm};
);    

getVariables = method();
getVaiables(MixedGraph) := (G) -> (
    Un = G#graph#Graph;	   
    Bi = G#graph#Bigraph;    
    Di = G#graph#Digraph;
    edgesUn = edges(Un);
    edgesBi = edges(Bi);
    edgesDi = edges(Di);
    numVert = #vertices(G);
    variables = (toList(l_(1,1)..l_(numVert,numVert))) | (toList(p_(1,1)..p_(numVert,numVert)))| (toList(k_(1,1)..k_(numVert,numVert)));
    varListDi = {};
    varListBi = {};
    varListUn = {};
    for i to (#edgesDi - 1) do (
	varListDi = append(varListDi, l_(edgesDi#i#0, edgesDi#i#1));
    )
    for i to (#edgesUn - 1) do (
	
    )
);

covarianceMatrixOfMixedGraph = method();
covarianceMatrixOfMixedGraph(MixedGraph) := (G) -> (
    R := gaussianRing(G);
    L := directedEdgesMatrix R;
    
);

scoreEquationsFromCovarianceMatrix = method();
scoreEquationsFromCovarianceMatrix(MixedGraph,List) := (G, U) -> (
    L1 = reorderingOfMixedGraph(G);
    G = L1#0;
    perm = L1#1;
-- we have to shuffle U according to perm
    variables = getVariables(G);
    S := sampleCovarianceMatrix(U);   
    R := gaussianRing(G); 
--    use R;   
    -- Lambda
    L := directedEdgesMatrix R;
    -- d is equal to the number of vertices in G
    d := numRows L;
    -- Omega
    W := bidirectedEdgesMatrix R;
    -- move to a new ring, lpR, which does not have the s variables
    numSvars:=lift(d*(d+1)/2,ZZ);
    --lp rings is the ring without the s variables
    lpRvarlist:=apply(numgens(R)-numSvars,i->(gens(R))_i);
    KK:=coefficientRing(R);
    lpR:=KK[lpRvarlist];
    lpRTarget:=apply(numgens(R),i-> if i<= numgens(R)-numSvars-1 then (gens(lpR))_i else 0);
    F:=map(lpR,R,lpRTarget);
    L = matRtolpR(L,F);
    W = matRtolpR(W,F);
    FR := frac(lpR);
    K := inverse (id_(lpR^d)-L);
    Sigma := (transpose K) * W * K;
    SigmaInv := inverse substitute(Sigma, FR);    
    C1 := trace(SigmaInv * V)/2;
    C1derivative := JacobianMatrixOfRationalFunction(trace(SigmaInv * V)/2);
    LL := (substitute((jacobian(matrix{{det(S)}})), FR))*matrix{{(-1/(2*det(S)))}} - (transpose C1derivative);
    LL=flatten entries(LL);
    denoms := apply(#LL, i -> lift(denominator(LL_i), lpR));
    prod := product(denoms);
    J:=ideal apply(#LL, i -> lift(numerator(LL_i),lpR));
    J = saturate(J, prod);
    
    return J;
);
