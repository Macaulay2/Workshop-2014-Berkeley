newPackage(
	"Hypergraphs",
	Version => "0.0.2",
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
	
	-- Hypergraph accessors
	"vertices",
	"edges",
	"incidenceMatrix",
	"neighbors",

	-- Functions accepting a Hypergraph as input
	"isHypergraphSimple",
	"inducedSubhypergraph",
	"hypergraphDual"
}

{*
	To do:
	- Add various other constructors:
		- list of vertices + list of edges (done)
		- list of edges (done)
		- hash table of form {a => {0,1}, ...} indicating vertex a belongs to edges 0 and 1, and so on
		- incidence matrix + automatically assigned vertex names (done)
		- incidence matrix + user provided vertex names
		- constructors similar to those in EdgeIdeals (for reverse compatibility)
	- Overload the net function for the Hypergraph class
	- Add accessors:
		- vertices(Hypergraph)			- done
		- edges(Hypergraph)				- done
		- incidenceMatrix(Hypergraph)	- done
		- neighbors(Hypergraph, vertex) - done
		
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

--------------------------------------------------

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
    vComplement := select(H.hypergraphVertices, x -> not member(x, V));
    eComplement := unique flatten apply(vComplement, v -> H.hypergraphVertexContainments#v); --returns the indices of the edges to delete
    E := H.hypergraphEdges_(select(toList(0 .. #H.hypergraphEdges-1), e -> not member(e, eComplement)));
    hypergraph(V, E)
)

hypergraphDual = method(TypicalValue => Hypergraph);
hypergraphDual(Hypergraph) := Hypergraph => opts -> (H) -> hypergraph(transpose H.hypergraphIncidenceMatrix, opts)

isHypergraphSimple = method(TypicalValue => Boolean);
isHypergraphSimple(Hypergraph) := Boolean => H -> (
    if any(0 .. #H#hypergraphEdges-1, I -> any(0 .. I-1, J -> isSubset(H#hypergraphEdges#I, H#hypergraphEdges#J) or isSubset(H#hypergraphEdges#J, H#hypergraphEdges#I))) then false else true
)


--------------------------------------------------

beginDocumentation();

--test 0: tests of constructor, correctness of isSimple with a simple hypergraph, and correctness of inducedSubhypergraph
TEST ///
	V = toList(1 .. 5);
	E = {{2,3,4}, {3,4,5}, {1,2,4,5}, {1,3}};
	H = hypergraph(V, E);
	assert(edges H === E);
	assert(vertices H === V);
	assert(entries incidenceMatrix H === {{0,0,1,1}, {1,0,1,0}, {1,1,0,1}, {1,1,1,0}, {0,1,1,0}});
	assert(neighbors(H, V_0) === {2,4,5,3});
	assert(neighbors(H, V_1) === {3,4,1,5});
	assert(neighbors(H, V_2) === {2,4,5,1});
	assert(neighbors(H, V_3) === {2,3,5,1});
	assert(neighbors(H, V_4) === {3,4,1,2});
	P = pairs H.hypergraphVertexContainments;
	P' = {(1,{2,3}), (2,{0,2}), (3,{0,1,3}), (4,{0,1,2}), (5,{1,2})};
	assert(isSubset(P, P') and isSubset(P', P));
	assert(isHypergraphSimple H === true);
	S1 = inducedSubhypergraph(H, V);
	assert(S1 === H);
	V2 = {};
	S2 = inducedSubhypergraph(H, V2);
	assert(vertices S2 === V2);
	assert(edges S2 === {});
	assert(entries incidenceMatrix S2 === {});
	assert(pairs S2.hypergraphNeighbors === {});
	assert(pairs S2.hypergraphVertexContainments === {});
	V3 = {1,3,4,5};
	S3 = inducedSubhypergraph(H, V3);
	assert(vertices S3 === V3);
	assert(edges S3 === {{3,4,5}, {1,3}});
	assert(entries incidenceMatrix S3 === {{0,1}, {1,1}, {1,0}, {1,0}});
	P3 = pairs S3.hypergraphNeighbors;
	P3' = {(1, {3}), (3, {4,5,1}), (4, {3,5}), (5, {3,4})};
	assert(isSubset(P3, P3') and isSubset(P3', P3));
	P3 = pairs S3.hypergraphVertexContainments;
	P3' = {(1, {1}), (3, {0,1}), (4, {0}), (5, {0})};
	assert(isSubset(P3, P3') and isSubset(P3', P3));
///

--test 1: tests constructor with non-numeric vertices and correctness of isHypergraphSimple on a non-simple hypergraph with a simple induced subhypergraph
TEST ///
	V = {a,b,c,d,e};
	E = {{b,d}, {e}, {a,b,d,e}, {a,c}};
	H = hypergraph(V, E);
	assert(edges H === E);
	assert(vertices H === V);
	assert(isHypergraphSimple H === false);
	assert(isHypergraphSimple inducedSubhypergraph(H, {a,b,c,d}) === true);
///

--test 2: tests constructor with vertices of mixed type and a hypergraph with an isolated vertex
TEST ///
	R = ZZ[x];
	V = {x, a, 1, z -> z};
	E = {{V_3}, {V_0, V_2}};
	H = hypergraph(V, E);
	assert(vertices H === V);
	assert(edges H === E);
	assert(entries incidenceMatrix H === {{0,1}, {0,0}, {0,1}, {1,0}});
	P = pairs H.hypergraphVertexContainments;
	P' = {(V_0, {1}), (V_1, {}), (V_2, {1}), (V_3, {0})};
	assert(isSubset(P, P') and isSubset(P', P));
///

end
