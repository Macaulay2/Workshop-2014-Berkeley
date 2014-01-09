newPackage(
	"Hypergraphs",
	Version => "0.0.1",
	Date => "3-5-2012",
	Authors => {
		{  Name => "Augustine O'Keefe", Email => "???" },
		{  Name => "Nicholas Armenoff", Email => "nicholas.armenoff@uky.edu" },
	},
	DebuggingMode => true,
	Reload => true
)

export {
	Hypergraph,
	hypergraph,
	Singletons,
	isHypergraphSimple,
	edges,
	vertices,
	incidenceMatrix,
	vertexContainments,
	neighbors
}

--to do:
-- 1 - add a net function for the Hypergraph class

--the classes defined in this package
Hypergraph = new Type of HashTable;

--constructor for Hypergraph class
hypergraph = method(TypicalValue => Hypergraph, Options => {Singletons => null});
hypergraph(List, List) := Hypergraph => opts -> (V, E) -> (
    E = toList \ E; --force E to be a list of lists
    if any(E, e -> not instance(e, List)) then error "Edges must be lists.";
    if any(E, e -> not isSubset(e, V)) then error "Edges must be subsets of the vertices.";
    
    V = unique join(V, if instance(opts.Singletons, List) then opts.Singletons else {});
    A := if #V == 0 then map(ZZ^0,ZZ^0,0) else matrix apply(#V, i -> apply(#E, j -> if member(V_i, E_j) then 1 else 0));
    vContainments := for v in V list ( v => for i from 0 to #E-1 list ( if member(v, E#i) then i else continue ) );
	nbors := for c in vContainments list ( unique sort delete(c#0, flatten E_(c#1)) );
    
    return new Hypergraph from hashTable {
		edges => E,
		vertices => V,
		incidenceMatrix => A,
		vertexContainments => hashTable vContainments,
		neighbors => nbors
	};
)

isHypergraphSimple = method(TypicalValue => Boolean);
isHypergraphSimple(Hypergraph) := Boolean => H -> (
    if any(0 .. #H#edges-1, I -> any(0 .. I-1, J -> isSubset(H#edges#I, H#edges#J) or isSubset(H#edges#J, H#edges#I))) then return false else return true;
)
