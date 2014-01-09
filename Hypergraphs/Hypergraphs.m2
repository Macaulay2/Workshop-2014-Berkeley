restart

loadPackage "Graphs"
Hypergraph = new Type of HashTable

hypergraph = method(Options => {symbol Singletons => null})

hypergraph HashTable
hypergraph (List, List) := Hypergraph => opts -> (V,E) -> (
    E = toList \ E;
    if any(E, e -> not instance(e, List)) then error "Edges must be lists.";
    if any(E, e -> not isSubset(e,V)) then error "Edges must be subsets of the vertices.";
    V = unique join(V, if instance(opts.Singletons, List) then opts.Singletons else {});
    A := if #V == 0 then map(ZZ^0,ZZ^0,0) else matrix apply(#V, i -> apply(#E, j -> if member(V_i,E_j) then 1 else 0));
    H := new Hypergraph from hashTable { symbol edges => E,
	symbol vertices => V,
	symbol adjacencyMatrix => A};
    return H
    )


isSimple = method ()
isSimple Hypergraph := Boolean => H -> (
    if any(0..#H#edges-1, I -> any(0..I-1, J-> isSubset(H#edges#I, H#edges#J) or isSubset(H#edges#J, H#edges#I)))then error "Edges satisfy a inclusion relation, i.e. is not simple";
     
     

V = toList (a..e)
E = {{b,c,d},{c,d,e},{a,b,d,e},{a,c}}
H = hypergraph (V,E)
