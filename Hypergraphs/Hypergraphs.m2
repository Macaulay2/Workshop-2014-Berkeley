newPackage(
	"Hypergraphs",
	Version => "0.0.1",
	Date => "3-5-2012",
	Authors => {
		{  Name => "Augustine O'Keefe", Email => "augustine.okeefe@gmail.com" },
		{  Name => "Nicholas Armenoff", Email => "nicholas.armenoff@uky.edu" },
	},
	DebuggingMode => true,
	Reload => true
)

export {
	-- Hypergraph class, constructor, options, and properties
	"Hypergraph",
	"hypergraph",
	"Singletons",
	"hypergraphEdges",
	"hypergraphVertices",
	"hypergraphIncidenceMatrix",
	"hypergraphVertexContainments",
	"hypergraphNeighbors",
	
	"vertices",
	"edges",
	"incidenceMatrix",
	"neighbors",

	-- Functions accepting a hypergraph as input
	"isHypergraphSimple",
	"inducedSubhypergraph",
	"hypergraphDual"
}

{*
	To do:
	- add various other constructors (contructors accepting incidence matrix (done), list of edges (done))
	- overload the net function for the Hypergraph class
	- Add accessors:
		- vertices(Hypergraph)
		- edges(Hypergraph)
		- incidenceMatrix(Hypergraph)
		- neighbors(Hypergraph, vertex)

	Functions We Might Want to Include:
	From EdgeIdeals:
	- chromaticNumber
	- complementGraph
	- connectedComponents
	- numberConnectedComponents
	- hypergraphToSimplicialComplex
	- coverIdeal (This should maybe be of class MonomialIdeal?
	- edgeIdeal (Also MonomialIdeal?)
	- deleteEdges
	- independenceComplex
	- lineGraph
	- simplicialComplexToHypergraph
	- vertexCovers
	- vertexCoverNumber
	- isGraph
	- isCM
	- isSCM
	- isConnected
	- isForest
	- isLeaf

	From Nauty:
	- areIsomorphic
	- addEdges
	- generateHypergraphs
	
	Future priorities:
	- Implement weighted hypergraphs
	- Implement directed hypergraphs
*}

--------------------------------------------

--the classes defined in this package
Hypergraph = new Type of HashTable;

--constructors for Hypergraph class
hypergraph = method(TypicalValue => Hypergraph, Options => {Singletons => null});
hypergraph(List, List) := Hypergraph => opts -> (V, E) -> (
    if any(E, e -> not instance(e, BasicList)) then error "edges must be Lists";
    E = toList \ E; --force E to be a List containing Lists
    if any(E, e -> not isSubset(e, V)) then error "edges must be subsets of the vertex set";
    
    V = unique join(V, if instance(opts.Singletons, List) then opts.Singletons else {});
    A := if #V == 0 then map(ZZ^0,ZZ^0,0) else matrix apply(#V, i -> apply(#E, j -> if member(V_i, E_j) then 1 else 0));
    vContainments := for v in V list ( v => for i from 0 to #E-1 list ( if member(v, E#i) then i else continue ) );
	nbors := for c in vContainments list ( c#0 => unique delete(c#0, flatten E_(c#1)) );
    
    new Hypergraph from hashTable {
		hypergraphEdges => E,
		hypergraphVertices => V,
		hypergraphIncidenceMatrix => A,
		hypergraphVertexContainments => hashTable vContainments, --keys are vertices and values are lists of edges numbered 0 through #E-1
		hypergraphNeighbors => hashTable nbors
	}
)

hypergraph(List) := Hypergraph => opts -> E -> hypergraph(unique flatten E, E, opts)

--Output: returns a Hypergraph given an incidence matrix.  The vertices are 0 .. numRows(incMatrix)-1.
hypergraph(Matrix) := Hypergraph => opts -> (incMatrix) -> (
	V := toList(0 .. numRows(incMatrix)-1);
	E := for j from 0 to numColumns(incMatrix)-1 list (
		for i from 0 to numRows(incMatrix)-1 list (
			if incMatrix_(i,j) != 0 then i
			else continue
		)
	);

	hypergraph(V, E, opts)
)

vertices = method(TypicalValue => List);
vertices(Hypergraph) := List => (H) -> H.hypergraphVertices

edges = method(TypicalValue => List);
edges(Hypergraph) := List => (H) -> H.hypergraphEdges

incidenceMatrix = method(TypicalValue => Matrix);
incidenceMatrix(Hypergraph) := Matrix => (H) -> H.hypergraphIncidenceMatrix

neighbors = method(TypicalValue => List);
neighbors(Hypergraph, Thing) := List => (H, vertex) -> (
	if member(vertex, keys H.hypergraphNeighbors) then H.hypergraphNeighbors#vertex
	else error "expected a vertex in the vertex set of H"
)

inducedSubhypergraph = method(TypicalValue => Hypergraph);
inducedSubhypergraph(Hypergraph, List) := Hypergraph => (H, V) -> (
	if any(V, x -> not member(x, H.hypergraphVertices)) then error "expected a subset of the vertex set of H";
    vComplement := select (H.hypergraphVertices, x -> not member(x, V));
    eComplement := unique flatten apply (vComplement, v -> H.hypergraphVertexContainments#v); --returns the indices of the edges to delete
    E := H.hypergraphEdges_(select(toList(0 .. #H.hypergraphEdges-1), e -> not member (e, eComplement)));
    hypergraph(V, E)
)

hypergraphDual = method(TypicalValue => Hypergraph);
hypergraphDual(Hypergraph) := Hypergraph => opts -> (H) -> hypergraph(transpose H.hypergraphIncidenceMatrix, opts)

isHypergraphSimple = method(TypicalValue => Boolean);
isHypergraphSimple(Hypergraph) := Boolean => H -> (
    if any(0 .. #H#hypergraphEdges-1, I -> any(0 .. I-1, J -> isSubset(H#hypergraphEdges#I, H#hypergraphEdges#J) or isSubset(H#hypergraphEdges#J, H#hypergraphEdges#I))) then false else true
)
