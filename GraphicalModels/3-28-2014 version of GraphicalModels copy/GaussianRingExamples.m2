Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases, PrimaryDecomposition, ReesAlgebra, TangentCone

i1 : loadPackage("GraphicalModels")

o1 = GraphicalModels

o1 : Package

i2 : --If you enter a graph you get the gaussianRing of a graph:
     g=graph {{1,3},{2,4}};

i3 : H=new HashTable from {Graph=>g, Digraph=>digraph {}, Bigraph => bigraph {} };

i4 : gG = new MixedGraph from { {graph,H}, {cache,new CacheTable from {}} }

o4 = MixedGraph{Bigraph => Bigraph{}        }
                Digraph => Digraph{}
                Graph => Graph{1 => set {3}}
                               2 => set {4}
                               3 => set {1}
                               4 => set {2}

o4 : MixedGraph

i5 : gaussianRing(gG)

o5 = QQ[k   , k   , k   , k   , k   , k   , s   , s   , s   , s   , s   , s   , s   , s   , s   , s   ]
         1,1   2,2   3,3   4,4   1,3   2,4   1,1   1,2   1,3   1,4   2,2   2,3   2,4   3,3   3,4   4,4

o5 : PolynomialRing

i6 : gaussianRing(g)

o6 = QQ[k   , k   , k   , k   , k   , k   , s   , s   , s   , s   , s   , s   , s   , s   , s   , s   ]
         1,1   2,2   3,3   4,4   1,3   2,4   1,1   1,2   1,3   1,4   2,2   2,3   2,4   3,3   3,4   4,4

o6 : PolynomialRing

i7 : --If you enter a digraph you get the gaussianRing of a digraph:
     gD = mixedGraph(digraph {{1,{2}},{2,{3,4}}});

i8 : gaussianRing(gD)

o8 = QQ[s   , s   , s   , s   , s   , s   , s   , s   , s   , s   ]
         1,1   1,2   1,3   1,4   2,2   2,3   2,4   3,3   3,4   4,4

o8 : PolynomialRing

i9 : gaussianRing(digraph {{1,{2}},{2,{3,4}}})

o9 = QQ[s   , s   , s   , s   , s   , s   , s   , s   , s   , s   ]
         1,1   1,2   1,3   1,4   2,2   2,3   2,4   3,3   3,4   4,4

o9 : PolynomialRing

i10 : --Example: chain graph (a graph with no bidirected edges)
      gGD = mixedGraph(graph {{1,3},{2,4}} , digraph {{1,{2}},{2,{3,4}}});

i11 : gGD = mixedGraph(graph {{1,3},{2,4}} , digraph {{1,{2}},{2,{3,4}}})

o11 = MixedGraph{Bigraph => Bigraph{}               }
                 Digraph => Digraph{1 => set {2}   }
                                    2 => set {3, 4}
                                    3 => set {}
                                    4 => set {}
                 Graph => Graph{1 => set {3}}
                                2 => set {4}
                                3 => set {1}
                                4 => set {2}

o11 : MixedGraph

i12 : gaussianRing(gGD)

o12 = QQ[l   , l   , l   , k   , k   , k   , k   , k   , k   , s   , s   , s   , s   , s   , s   , s   , s   , s   , s   ]
          1,2   2,4   2,3   1,1   2,2   3,3   4,4   1,3   2,4   1,1   1,2   1,3   1,4   2,2   2,3   2,4   3,3   3,4   4,4

o12 : PolynomialRing

i13 : gDB = mixedGraph(digraph {{1,{2}},{2,{3,4}}}, bigraph {{1,3},{2,4}})

o13 = MixedGraph{Bigraph => Bigraph{1 => set {3}}   }
                                    2 => set {4}
                                    3 => set {1}
                                    4 => set {2}
                 Digraph => Digraph{1 => set {2}   }
                                    2 => set {3, 4}
                                    3 => set {}
                                    4 => set {}
                 Graph => Graph{}

o13 : MixedGraph

i14 : gaussianRing(gDB)

o14 = QQ[l   , l   , l   , p   , p   , p   , p   , p   , p   , s   , s   , s   , s   , s   , s   , s   , s   , s   , s   ]
          1,2   2,4   2,3   1,1   2,2   3,3   4,4   1,3   2,4   1,1   1,2   1,3   1,4   2,2   2,3   2,4   3,3   3,4   4,4

o14 : PolynomialRing
