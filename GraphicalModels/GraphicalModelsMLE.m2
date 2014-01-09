-- -*- coding: utf-8-unix -*-

newPackage(
     "GraphicalModelsMLE",
     Version => "0.1",
     Date => "January 8, 2014",
     Authors => {
	  {Name => "Luis Garcia-Puente",
	   Email => "lgarcia@shsu.edu",
	   HomePage => "http://www.shsu.edu/~ldg005"},
          {Name=> "David Swinarski", 
	   Email=> "dswinarski@fordham.edu",
	   HomePage=>"http://faculty.fordham.edu/dswinarski"}, 
          {Name=> "Elina Robeva", 
	   Email=> "erobeva@gmail.com",
	   HomePage=>"http://math.berkeley.edu/~erobeva"}
	  },
     Headline => "MLE estimates for structural equation models",
     DebuggingMode => true
     )
export {"sampleCovarianceMatrix",
    "JacobianMatrixOfRationalFunction",
    "scoreEquationsFromCovarianceMatrix"
       	} 
     
needsPackage "Graphs"     
needsPackage "GraphicalModels"
needsPackage "Bertini"



--**************************--
--  INTERNAL ROUTINES       --
--**************************--



--*************************************--
--  Functions (local) used throughout  --
--*************************************--




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



--**************************--
--  METHODS 	      	   	  --
--**************************--
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

scoreEquationsFromCovarianceMatrix = method();
scoreEquationsFromCovarianceMatrix(MixedGraph,List) := (G, U) -> (
    V := sampleCovarianceMatrix(U);   
    R := gaussianRing(G); 
    use R;   
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
    S := (transpose K) * W * K;
    Sinv := inverse substitute(S, FR);    
    C1 := trace(Sinv * V)/2;
    C1derivative := JacobianMatrixOfRationalFunction(trace(Sinv * V)/2);
    LL := (substitute((jacobian(matrix{{det(S)}})), FR))*matrix{{(-(#U)/(2*det(S)))}} - (transpose C1derivative);
    LL=flatten entries(LL);
    J:=ideal apply(#LL, i -> lift(numerator(LL_i),lpR));
    J = saturate(J, det(S));
    return J
);

scoreEquationsCovariance1 = method();
scoreEquationsCovariance1(MixedGraph, List) := (G, U) -> (
    R := gaussianRing G;
    V := sampleCovarianceMatrix(U);
    use R;
    L := directedEdgesMatrix R;
    -- d is equal to the number of vertices in G
    d := numRows L;
    -- Omega
    W := bidirectedEdgesMatrix R;
    K := inverse (id_(R^d)-L);
    -- the matrix Sigma expressed in terms of the l and p-variables.
    S := (transpose K) * W * K;
    -- Get the ideal in l's and p's.
    J := scoreEquationsFromCovarianceMatrix(G, U);
    -- same ideal in the ring R.
    J1 := substitute(J, R);
    Sigma := mutableMatrix(R, d, d); 
    for i to d-1 do (
	for j from i to d-1 do (
	    use R;
	    Sigma_(i,j) = s_(i+1,j+1);
	    use R;
 	    Sigma_(j,i) = s_(i+1,j+1);
	);
    );
    J2 := ideal(matrix(Sigma) - S);
    J3 := J1 + J2;
    numSvars := lift(d*(d+1)/2,ZZ);
    -- the l and p-variables
    lpRvarlist := apply(numgens(R)-numSvars,i->(gens(R))_i);
    -- the s-variables
    sRvarlist := apply(numSvars, i->(gens(R))_(i + numgens(R) - numSvars));
    J4 := eliminate(lpRvarlist, J3);
    sR := coefficientRing(R)[sRvarlist];
    J5 := substitute(J4, sR);
    return J5;
);


--******************************************--
-- DOCUMENTATION     	       	    	    -- 
--******************************************--

beginDocumentation()

doc ///
    Key
        GraphicalModelsMLE
    Headline
        a package for MLE estimates of parameters for statistical graphical models 
    Description        
        Text
            Add some text and an example.   
	    
	    In the example below, we create the score equations (defining the critical points of the log likelihood function written in terms of the covariance matrix) associated to the four data vectors $(1,2,1,-1)$, $(2,1,3,0)$, $(-1,0,1,1)$, $(-5,3,4,-6)$ for a graphical model with four vertices, five directed edges, and one bidirected edge.
	    
        Example
	    needsPackage("Graphs");
	    needsPackage("GraphicalModels");
	    needsPackage("GraphicalModelsMLE");
	    G = mixedGraph(digraph {{1,{2,3}},{2,{3}},{3,{4}}},bigraph {{3,4}})
	    U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}
            scoreEquationsFromCovarianceMatrix(G,U)
	    scoreEquationsCovariance1(G, U)
    Caveat
        GraphicalModelsMLE requires Graphs.m2 and GraphicalModels.m2. 
///;

--------------------------------
-- Documentation
--------------------------------

doc /// 
    Key
        sampleCovarianceMatrix
        (sampleCovarianceMatrix,List) 
    Headline
        compute the sample covariance matrix of a list of data vectors
    Usage
        sampleCovarianceMatrix U
    Inputs
        U:List
    Outputs
         :Matrix
    Description 
        Text
	    The sample covariance matrix is $S = \frac{1}{n} \sum_{i=1}^{n} (X^{(i)}-\bar{X}) (X^{(i)}-\bar{X})^T$.  The entries here are not the unbiased estimators of the variance/covariance; that is, the entries here correspond to the outputs of the commands VAR.P and COVARIANCE.P in Excel, not VAR and COVARIANCE in Excel.
	    
	    We assume that the data vectors are entered as a list of row matrices, all with the same width.
        Example
            M = {matrix{{1, 2, 0}}, matrix{{-1, 0, 5/1}}, matrix{{3, 5, 2/1}}, matrix{{-1, -4, 1/1}}};
	    sampleCovarianceMatrix(M)
	    U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}};
	    sampleCovarianceMatrix(U)	    
///

doc /// 
    Key
        JacobianMatrixOfRationalFunction
        (JacobianMatrixOfRationalFunction,RingElement) 
    Headline
        compute the Jacobian matrix of a rational function
    Usage
        JacobianMatrixOfRationalFunction(F)
    Inputs
        F:RingElement
    Outputs
         :Matrix
    Description 
        Text
	    This function computes the Jacobian matrix of a rational function
        Example
	    R=QQ[x,y];
	    FR=frac R;
	    F=1/(x^2+y^2);
            JacobianMatrixOfRationalFunction(F)
	    R=QQ[t_1,t_2,t_3];
	    FR=frac R;
	    JacobianMatrixOfRationalFunction( (t_1^2*t_2)/(t_1+t_2^2+t_3^3) )
///

doc /// 
    Key
        scoreEquationsFromCovarianceMatrix
        (scoreEquationsFromCovarianceMatrix,MixedGraph,List) 
    Headline
        computes the score equations that arise from the log likelihood formula in terms of the covariance matrix Sigma
    Usage
        scoreEquationsFromCovarianceMatrix(G,U)
    Inputs
        G:MixedGraph
	U:List
    Outputs
         :Ideal
    Description 
        Text
	    This function computes the score equations that arise from the log likelihood formula in terms of the covariance matrix $\Sigma$.  
	    Suppose we are given a list of data vectors, each a row matrix.  
	    The log likelihood function we want to maximize is given in Prop. 2.1.12 of Sturmfels's lecture notes (to do: update this reference to the printed version).  
	    
        Example
	    needsPackage("Graphs");
	    needsPackage("GraphicalModels");
	    G = mixedGraph(digraph {{1,{2,3}},{2,{3}},{3,{4}}},bigraph {{3,4}})
	    U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}
            scoreEquationsFromCovarianceMatrix(G,U)
///


doc /// 
    Key
        scoreEquationsCovariance1
        (scoreEquationsCovariance1,MixedGraph,List) 
    Headline
        computes the score equations that arise from the log likelihood formula in terms of the covariance matrix Sigma
    Usage
        scoreEquationsCovariance1(G,U)
    Inputs
        G:MixedGraph
	U:List
    Outputs
         :Ideal
    Description 
        Text
	    This function computes the score equations that arise from the log likelihood formula in terms of the covariance matrix $\Sigma$.  
	    Suppose we are given a list of data vectors, each a row matrix.  
	    The log likelihood function we want to maximize is given in Prop. 2.1.12 of Sturmfels's lecture notes (to do: update this reference to the printed version).  
	    
        Example
	    needsPackage("Graphs");
	    needsPackage("GraphicalModels");
	    needsPackage("GraphicalModelsMLE");
	    G = mixedGraph(digraph {{1,{2,3}},{2,{3}},{3,{4}}},bigraph {{3,4}})
	    U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}
            scoreEquationsCovariance1(G,U)
///


--******************************************--
-- TESTS     	       	    	      	    --
--******************************************--

TEST /// 
R=QQ[x,y];
FR=frac R;
F=1/(x^2+y^2);
M=entries JacobianMatrixOfRationalFunction(F);
N={{-2*x/(x^2 + y^2)^2,-2*y/(x^2 + y^2)^2 }};
assert(M === N)
///

TEST ///
R=QQ[x_1,x_2,x_3];
FR=frac R;
M=entries JacobianMatrixOfRationalFunction( (x_1^2*x_2)/(x_1+x_2^2+x_3^3) );
N={{2*x_1*x_2/(x_2^2 + x_3^3 + x_1) - x_1^2*x_2/(x_2^2 + x_3^3 + x_1)^2, -2*x_1^2*x_2^2/(x_2^2 + x_3^3 + x_1)^2 + x_1^2/(x_2^2 + x_3^3 + x_1) , -3*x_1^2*x_2*x_3^2/(x_2^2 + x_3^3 + x_1)^2 }};
assert(M === N)
/// 

TEST ///
M = {matrix{{1, 2, 0}}, matrix{{-1, 0, 5/1}}, matrix{{3, 5, 2/1}}, matrix{{-1, -4, 1/1}}};
N = sampleCovarianceMatrix(M);
A = matrix {{11/4, 39/8, -1}, {39/8, 171/16, 0}, {-1, 0, 7/2}};
assert(N===A)	
///

TEST ///
X = {matrix {{36, -3, -25, -36}}, matrix {{-10, 11, -29, -20}}, matrix {{-18, 33, -15, -11}}, matrix {{-42, 0, 20, 43}}, matrix {{-30, -26, 32, 2}}, matrix {{2, -38, -24, -43}} };
Y = sampleCovarianceMatrix(X);
B = matrix matrix {{5621/9, -1037/18, -7835/18, -10565/18}, {-1037/18, 19505/36, -4897/36, 5147/36}, {-7835/18, -4897/36, 20465/36, 18941/36}, {-10565/18, 5147/36, 18941/36, 28889/36}};
assert(Y===B)	
///

TEST ///
needsPackage("Graphs");
needsPackage("GraphicalModels");
G = mixedGraph(digraph {{1,{2,3}},{2,{3}},{3,{4}}},bigraph {{3,4}})
U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}
J=scoreEquationsFromCovarianceMatrix(G,U);
I=ideal(80*p_(3,4)+39,200*p_(4,4)-271,1760416*p_(3,3)-742363,920*p_(2,2)-203,64*p_(1,1)-115,5*l_(3,4)+2,110026*l_(2,3)-2575,55013*l_(1,3)-600,115*l_(1,2)+26);
assert(J===I)
///     

TEST ///
needsPackage("Graphs");
needsPackage("GraphicalModels");
G = mixedGraph(digraph {{1,{2,3}},{2,{3}},{3,{4}}},bigraph {{3,4}})
U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}
J=scoreEquationsCovariance1(G,U);
I=ideal(16*s_(4,4)-29,32*s_(3,4)+21,64*s_(3,3)-27,4336*s_(2,4)+5,8672*s_(2,3)-25,16*s_(2,2)-5,8672*s_(1,4)+35,17344*s_(1,3)-175,32*s_(1,2)+13,64*s_(1,1)-115);
assert(J===I)
///     
--------------------------------------
--------------------------------------
end
--------------------------------------
--------------------------------------


--blank documentation node:
doc /// 
    Key
       gaussianMatrix
       (gaussianMatrix,Digraph,Matrix,List) 
    Headline
    Usage
    Inputs
    Outputs
    Description 
        Text
        Example
    SeeAlso
///


uninstallPackage "GraphicalModelsMLE"
restart
--installPackage("Graphs", UserMode=>true)
installPackage ("GraphicalModelsMLE", RemakeAllDocumentation => true, UserMode=>true)
viewHelp GraphicalModelsMLE
installPackage("GraphicalModelsMLE",UserMode=>true,DebuggingMode => true)


----------------------
-- Parameterization -- ????????????????????????????????????????????????????????????????????????
---------------------- 
---- We need this for both directed and undirected graphs:

----  parameterizations and for toric varieties the corresponding matrix. 
----  In the case of toric varieties the matrix is easy.  Here is the code, 
----  commented out to be used later when we are ready. 
---- 
----  toAMatrix = method()
----  toAMatrix List := Matrix => (M) -> (
----      if any(M,isMonomial)
----         then error "this parameterization does not correspond to a toric ideal." 
----         else (
----              Mexp := apply(M, exponents);
----              transpose matrix apply(Mexp, flatten)))
----
---- isMonomial = method()
---- isMonomial RingElement := Boolean => (m) -> (
----      termList := terms m;
----      if #termList == 1 then true else false)

---- isMonomial works well as long as m is actually a polynomial or monomial and not 
---- an element of ZZ, QQ, RR, etc.


end;
restart
installPackage "GraphicalModelsMLE"
check "GraphicalModelsMLE"


