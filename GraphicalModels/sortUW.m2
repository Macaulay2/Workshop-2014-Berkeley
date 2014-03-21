

-----------------------------------
-- This function takes a mixed graph and a permutation and
-- creates the same graph but with permuted vertices.
-----------------------------------

permuteVerticesOfGraph = (G, perm) -> (
    perm1 := apply(#perm, i -> (perm#i +1));
    edgesUn := edges(Un);
    edgesBi := edges(Bi);
    edgesDi := edges(Di);
    newEdgesUn := {};
    newEdgesBi := {};
    newEdgesDi := {};
    ed:={};
    ed1:={};
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
    newG := mixedGraph(graph(newEdgesUn), digraph(newEdgesDi), bigraph(newEdgesBi));
    return newG;
);

------------------------------------------
-- This function separates the vertices into two groups: the undirected
-- and the bidirected subgraphs. If a vertex does not have any adjacent
-- undirected or bidirected edge, then, it is distributed to one of the
-- two sets in a way to not create collisions.
-----------------------------------------

orderingOfMixedGraphUndirectedBidirectedSeparation = (G) -> (
    Un := G#graph#Graph;	   
    Bi := G#graph#Bigraph;    
    Di := G#graph#Digraph;
    vertUn := vertices(Un);
    vertBi := vertices(Bi);
    vertDi := vertices(Di);
    flag:=true;
    ed:={};
    vertRest := toList(set(vertDi) - (set(vertUn)+set(vertBi))); 
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
    numUn := #vertUn;
    numBi := #vertBi;
    perm := vertUn | vertBi;
    perm1 := apply(#perm, i -> (perm#i-1));
    perm2 := inversePermutation perm1;
    newG := permuteVerticesOfGraph(G, perm2);
    return {newG, perm2, numUn, numBi};    
);

--------------------------------------------
-- This function rearranges the vertices within the undirected and bidirected
-- part so that all directed edges go from a smaller vertex to a bigger vertex.
--------------------------------------------

orderingOfMixedGraphDirectedOrdering = (G) -> (    
    Un := G#graph#Graph;	   
    Bi := G#graph#Bigraph;    
    Di := G#graph#Digraph;
    vertUn := vertices(Un);
    vertBi := vertices(Bi);
    vertUn = sort(vertUn);
    vertBi = sort(vertBi);
    vertRest := sort(toList(set(vertDi) - (set(vertUn)+set(vertBi)))); 
    flag:=true;
    ed:={};
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
    numUn := #vertUn;
    numBi := #vertBi;
    if (numUn + numBi < #(vertices(G))) then (
    	print ("Problem!");
    );
    -- arrange the directed edges going from small to large in the undirected
    -- part of the graph.
    UnM := mutableMatrix(ZZ, 1, numUn);
    for i to (#edges(Di)-1) do (
	ed = (edges(Di))#i;
	print(ed);
	if (ed#0 <= numUn and ed#1 <= numUn) then (
    	    UnM_(0, ed#1-1) = UnM_(0, ed#1-1) + 1;
	);
    );
    S := set{};
    for i to numUn-1 do (
    	S = S + set{{i+1, UnM_(0, i)}};
    );
    perm := {};
    -- this is not as efficient as it could be; for example we can use a heap
    numS := #S;
    el:={};
    ch:={};
    elch:={};
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
    BiM := mutableMatrix(ZZ, 1, numBi);
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
    perm1 := inversePermutation perm;
    newG := permuteVerticesOfGraph(G, perm1);
    return {newG, perm1};    
);

-------------------------------
-- This function takes a mixed graph and permutes the vertices so that the
-- requirements on p.9 from Sullivant, Talaska and Draisma's paper, called
-- Trek Separation for Gaussian Graphical Models, are satisfied.
------------------------------

orderingOfMixedGraph = (G) -> (
    l1 := orderingOfMixedGraphUndirectedBidirectedSeparation(G);
    l2 := orderingOfMixedGraphDirectedOrdering(l1#0);
    perm := {};
    for i to (#l1-1) do (
	perm = append(perm, l2#1#(l1#1#i - 1));
    );
    return {l2#0, perm, l1#2, l1#3};
);    
