--Goal: compute the UW decomposition
--and check the Sullivant-Talaska-Draisma conditions
--Do not permute (this can be added later)

UWDecomposition = (G) -> (
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
	error("The bidirected and undirected part cannot have any vertices in common!")
    );
    vertUn = sort(vertUn);
    vertBi = sort(vertBi);
    return {vertUn,vertBi};    
);


{*
From Sullivant-Talaska-Draisma, 
"Trek separation for Gaussian graphical models"
page 9, Section 2.3.
Sullivant-Talaska-Draisma conditions:
1. no directed edge from W to U
2. vertices in U come before vertices in W
3. if i->j is a directed edge then i<j
*}

checkSullivantTalaskaDraismaConditions = (G) -> (
    Un := G#graph#Graph;	   
    Bi := G#graph#Bigraph;    
    Di := G#graph#Digraph;
    vertUn := vertices(Un);
    vertBi := vertices(Bi);
    vertDi := vertices(Di);
    UW:=UWDecomposition(G);
    U:=UW_0;
    W:=UW_1;
    --Check condition 1: no directed edge from W to U
    edgesDi:=edges(Di);
    e:={};
    for i from 0 to #edgesDi-1 do (
	e=edgesDi_i;
	if isSubset(set({e_0}),set(W)) and isSubset(set({e_1}),set(U)) then (
	    print "There cannot be any directed edges from W to U";
	    return false   
	)
    );
    --Check condition 2: vertices in U come before vertices in W 
    if min(W) < max(U) then (
	print "The vertices in U must come before the vertices in W"; 
        return false	
    );
    --Check condition 3: if i->j is a directed edge then i<j
    for i from 0 to #edgesDi-1 do (
	e=edgesDi_i;
	if e_1 < e_0 then (
	    print "A directed edge i->j must satisfy i<j";
	    return false   
	)
    );
    return true
);




{*
Examples for checking: 

This graph should satisfy the Sullivant-Talaska-Draisma conditions:
G = mixedGraph(graph {{1,2}},digraph {{1,{2}},{2,{3,4}}},bigraph {{3,4}});
checkSullivantTalaskaDraismaConditions(G)

This graph should fail condition 1: no directed edge from W to U
G = mixedGraph(graph {{1,2}},digraph {{1,{2}},{3,{2,4}}},bigraph {{3,4}});
checkSullivantTalaskaDraismaConditions(G)

This graph should fail condition 2: vertices in U come before vertices in W
G = mixedGraph(graph {{3,4}},digraph {{3,{4}},{4,{1,2}}},bigraph {{1,2}});
checkSullivantTalaskaDraismaConditions(G)

This graph should fail condition 3: if i->j is a directed edge then i<j 
G = mixedGraph(graph {{1,2}},digraph {{2,{1}},{2,{3,4}}},bigraph {{3,4}});
checkSullivantTalaskaDraismaConditions(G)
*}