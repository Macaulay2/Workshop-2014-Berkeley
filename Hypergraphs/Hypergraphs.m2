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
	-- Hypergraph class and hypergraph properties
	"Hypergraph",
	"hypergraph",
	"Singletons",
	"edges",
	"vertices",
	"incidenceMatrix",
	"vertexContainments",
	"neighbors",

	-- Functions accepting a hypergraph as input
	"isHypergraphSimple",
	"inducedSubhypergraph"
}

{*
	To do:
	- add various other constructors (contructors accepting incidence matrix (done), list of edges (done))
	- overload the net function for the Hypergraph class

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
*}

--------------------------------------------

--the classes defined in this package
Hypergraph = new Type of HashTable;

--constructors for Hypergraph class
hypergraph = method(TypicalValue => Hypergraph, Options => {Singletons => null});
hypergraph(List, List) := Hypergraph => opts -> (V, E) -> (
    E = toList \ E; --force E to be a list of lists
    if any(E, e -> not instance(e, List)) then error "Edges must be lists.";
    if any(E, e -> not isSubset(e, V)) then error "Edges must be subsets of the vertices.";
    
    V = unique join(V, if instance(opts.Singletons, List) then opts.Singletons else {});
    A := if #V == 0 then map(ZZ^0,ZZ^0,0) else matrix apply(#V, i -> apply(#E, j -> if member(V_i, E_j) then 1 else 0));
    vContainments := for v in V list ( v => for i from 0 to #E-1 list ( if member(v, E#i) then i else continue ) );
	nbors := for c in vContainments list ( c#0 => unique delete(c#0, flatten E_(c#1)) );
    
    return new Hypergraph from hashTable {
		edges => E,
		vertices => V,
		incidenceMatrix => A,
		vertexContainments => hashTable vContainments, --keys are vertices and values are lists of edges numbered 0 through #E-1
		neighbors => nbors
	};
)

hypergraph(List) := Hypergraph => opts -> E -> (
	V := unique flatten E;
	return hypergraph(V, E, opts);
)

--Output: returns a Hypergraph given an incidence matrix.  The vertices are 0 .. numRows(incMatrix)-1.
hypergraph(Matrix) := Hypergraph => opts -> (incMatrix) -> (
	V := toList(0 .. numRows(incMatrix)-1);
	E := for j from 0 to numColumns(incMatrix)-1 list (
		for i from 0 to numRows(incMatrix)-1 list (
			if incMatrix_(i,j) != 0 then i
			else continue
		)
	);

	return hypergraph(V, E, opts);
)

<<<<<<< Updated upstream
--------------------------------------------
=======
>>>>>>> Stashed changes

inducedSubhypergraph = method(TypicalValue => Hypergraph);
inducedSubhypergraph(List,Hypergraph) := Hypergraph => (V,H) -> (
    vComplement := select (H.vertices, x -> not member(x,V));
    eComplement := unique flatten apply (vComplement, v -> H.vertexContainments#v); --returns the indices of the edges to delete
    E := H.edges_(select(toList(0 .. #H.edges-1), e -> not member (e, eComplement)));
    return hypergraph(V, E);
)

hypergraphDual = method(TypicalValue => Hypergraph);
hypergraphDual(Hypergraph) := Hypergraph => opts -> (h) -> return hypergraph(transpose h.incidenceMatrix, opts);

isHypergraphSimple = method(TypicalValue => Boolean);
isHypergraphSimple(Hypergraph) := Boolean => H -> (
    if any(0 .. #H#edges-1, I -> any(0 .. I-1, J -> isSubset(H#edges#I, H#edges#J) or isSubset(H#edges#J, H#edges#I))) then return false else return true;
)
