--An example from Luis's website
loadPackage("Graphs");
loadPackage("GraphicalModels");
loadPackage("Bertini");
loadPackage("GraphicalModelsMLE");
G = mixedGraph(digraph {{1,{2}},{2,{3,4}}},bigraph {{3,4}});
U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}};
--dim 0 deg 1

--Another example from Luis's website
loadPackage("Graphs");
loadPackage("GraphicalModels");
loadPackage("Bertini");
loadPackage("GraphicalModelsMLE");
G = mixedGraph(digraph {{1,{2}},{2,{3,4}}},bigraph {{1,3},{2,4}});
R=gaussianRing(G);
I=gaussianVanishingIdeal(R)
U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}};
J=scoreEquationsFromCovarianceMatrix(G,U);
dim J
degree J

--dim 0, deg 2, factors over QQ

ideal(184*p_(2,4)+39,12992*p_(1,3)+3105,2116*p_(4,4)-939,2637376*p_(3,3)-733435,16*p_(2,2)-5,64*p_(1,1)-115,203*l_(2,3)-107,23*l_(2,4)+35,l_(1,2)), 

ideal(52*p_(2,4)-23,12992*p_(1,3)+3105,338*p_(4,4)-313,2637376*p_(3,3)-733435,920*p_(2,2)-203,64*p_(1,1)-115,203*l_(2,3)-107,13*l_(2,4)+47,115*l_(1,2)+26)



--Produce the equations for the covariance matrix in terms of L and W

    R = gaussianRing G;
--    I = gaussianVanishingIdeal R;
    -- Lambda
    L = directedEdgesMatrix R;
    -- d is equal to the number of vertices in G
d=numRows(L);
    -- Omega
    W = bidirectedEdgesMatrix R;
    K = inverse (id_(R^d)-L);
    S = (transpose K) * W * K;


--First component
A = ideal(-p_(1,1)+s_(1,1),-l_(1,2)*p_(1,1)+s_(1,2),-l_(1,2)*l_(2,3)*p_(1,1)-p_(1,3)+s_(1,3),-l_(1,2)*l_(2,4)*p_(1,1)+s_(1,4),-l_(1,2)*p_(1,1)+s_(1,2),-l_(1,2)^2*p_(1,1)-p_(2,2)+s_(2,2),-l_(1,2)^2*l_(2,3)*p_(1,1)-l_(2,3)*p_(2,2)-l_(1,2)*p_(1,3)+s_(2,3),-l_(1,2)^2*l_(2,4)*p_(1,1)-l_(2,4)*p_(2,2)-p_(2,4)+s_(2,4),-l_(1,2)*l_(2,3)*p_(1,1)-p_(1,3)+s_(1,3),-l_(1,2)^2*l_(2,3)*p_(1,1)-l_(2,3)*p_(2,2)-l_(1,2)*p_(1,3)+s_(2,3),-l_(1,2)^2*l_(2,3)^2*p_(1,1)-l_(2,3)^2*p_(2,2)-2*l_(1,2)*l_(2,3)*p_(1,3)-p_(3,3)+s_(3,3),-l_(1,2)^2*l_(2,4)*l_(2,3)*p_(1,1)-l_(2,4)*l_(2,3)*p_(2,2)-l_(1,2)*l_(2,4)*p_(1,3)-l_(2,3)*p_(2,4)+s_(3,4),-l_(1,2)*l_(2,4)*p_(1,1)+s_(1,4),-l_(1,2)^2*l_(2,4)*p_(1,1)-l_(2,4)*p_(2,2)-p_(2,4)+s_(2,4),-l_(1,2)^2*l_(2,4)*l_(2,3)*p_(1,1)-l_(2,4)*l_(2,3)*p_(2,2)-l_(1,2)*l_(2,4)*p_(1,3)-l_(2,3)*p_(2,4)+s_(3,4),-l_(1,2)^2*l_(2,4)^2*p_(1,1)-l_(2,4)^2*p_(2,2)-2*l_(2,4)*p_(2,4)-p_(4,4)+s_(4,4))
B=ideal(184*p_(2,4)+39,12992*p_(1,3)+3105,2116*p_(4,4)-939,2637376*p_(3,3)-733435,16*p_(2,2)-5,64*p_(1,1)-115,203*l_(2,3)-107,23*l_(2,4)+35,l_(1,2));
C=A+B;

eliminate({l_(1,2), l_(2,4), l_(2,3), p_(1,1), p_(2,2), p_(3,3), p_(4,4), p_(1,3), p_(2,4)},C)

ideal(16*s_(4,4)-29,3248*s_(3,4)+1177,2637376*s_(3,3)-962415,16*s_(2,4)+11,3248*s_(2,3)-535,16*s_(2,2)-5,s_(1,4),12992*s_(1,3)+3105,s_(1,2),64*s_(1,1)-115)

16*s_(4,4)-29
3248*s_(3,4)+1177
2637376*s_(3,3)-962415
16*s_(2,4)+11
3248*s_(2,3)-535
16*s_(2,2)-5
s_(1,4)
12992*s_(1,3)+3105
s_(1,2)
64*s_(1,1)-115

matrix {{115/64, 0, -3105/12992, 0},
{0, 5/16, 535/3248, -11/16}, 
{-3105/12992, 535/3248, 962415/2637376, -1177/3248},
{0, -11/16, -1177/3248, 29/16}}

i28 : eigenvalues matrix {{115/64, 0, -3105/12992, 0},
      {0, 5/16, 535/3248, -11/16}, 
      {-3105/12992, 535/3248, 962415/2637376, -1177/3248},
      {0, -11/16, -1177/3248, 29/16}}

o28 = {1.82437 }
      {2.17514 }
      {.244516 }
      {.0427617}



--Second component
A = ideal(-p_(1,1)+s_(1,1),-l_(1,2)*p_(1,1)+s_(1,2),-l_(1,2)*l_(2,3)*p_(1,1)-p_(1,3)+s_(1,3),-l_(1,2)*l_(2,4)*p_(1,1)+s_(1,4),-l_(1,2)*p_(1,1)+s_(1,2),-l_(1,2)^2*p_(1,1)-p_(2,2)+s_(2,2),-l_(1,2)^2*l_(2,3)*p_(1,1)-l_(2,3)*p_(2,2)-l_(1,2)*p_(1,3)+s_(2,3),-l_(1,2)^2*l_(2,4)*p_(1,1)-l_(2,4)*p_(2,2)-p_(2,4)+s_(2,4),-l_(1,2)*l_(2,3)*p_(1,1)-p_(1,3)+s_(1,3),-l_(1,2)^2*l_(2,3)*p_(1,1)-l_(2,3)*p_(2,2)-l_(1,2)*p_(1,3)+s_(2,3),-l_(1,2)^2*l_(2,3)^2*p_(1,1)-l_(2,3)^2*p_(2,2)-2*l_(1,2)*l_(2,3)*p_(1,3)-p_(3,3)+s_(3,3),-l_(1,2)^2*l_(2,4)*l_(2,3)*p_(1,1)-l_(2,4)*l_(2,3)*p_(2,2)-l_(1,2)*l_(2,4)*p_(1,3)-l_(2,3)*p_(2,4)+s_(3,4),-l_(1,2)*l_(2,4)*p_(1,1)+s_(1,4),-l_(1,2)^2*l_(2,4)*p_(1,1)-l_(2,4)*p_(2,2)-p_(2,4)+s_(2,4),-l_(1,2)^2*l_(2,4)*l_(2,3)*p_(1,1)-l_(2,4)*l_(2,3)*p_(2,2)-l_(1,2)*l_(2,4)*p_(1,3)-l_(2,3)*p_(2,4)+s_(3,4),-l_(1,2)^2*l_(2,4)^2*p_(1,1)-l_(2,4)^2*p_(2,2)-2*l_(2,4)*p_(2,4)-p_(4,4)+s_(4,4))
D=ideal(52*p_(2,4)-23,12992*p_(1,3)+3105,338*p_(4,4)-313,2637376*p_(3,3)-733435,920*p_(2,2)-203,64*p_(1,1)-115,203*l_(2,3)-107,13*l_(2,4)+47,115*l_(1,2)+26);
E=A+D;

eliminate({l_(1,2), l_(2,4), l_(2,3), p_(1,1), p_(2,2), p_(3,3), p_(4,4), p_(1,3), p_(2,4)},E)

16*s_(4,4)-29
6496*s_(3,4)+3623
64*s_(3,3)-27
16*s_(2,4)+11
32*s_(2,3)-7
16*s_(2,2)-5
32*s_(1,4)-47
64*s_(1,3)+29
32*s_(1,2)+13
64*s_(1,1)-115

eigenvalues matrix {{115/64, -13/32, -29/64, 47/32},
{-13/32, 5/16, 7/32, -11/16},
{-29/64, 7/32, 27/64, -3623/6496},
{47/32, -11/16, -3623/6496, 29/16}}

o35 = {3.63809  }
      {.470419  }
      {.00998306}
      {.225253  }
      
--To do: construct the Hessian and test whether it is positive definite, negative definite, or neither


--Example 3: same graph as Example 2, but with more generic data
U=apply(4, i -> matrix {apply(4, j -> random(-50,50)/(random(1,50)))})

loadPackage("Graphs");
loadPackage("GraphicalModels");
loadPackage("Bertini");
loadPackage("GraphicalModelsMLE");
G = mixedGraph(digraph {{1,{2}},{2,{3,4}}},bigraph {{1,3},{2,4}});
R=gaussianRing(G);
I=gaussianVanishingIdeal(R)
--U = {matrix{{1.05,2.1,1.0,-1.1}}, matrix{{2.1,1.0,3.2,0.2}}, matrix{{-1.1, 0.1, 1.1, 1.1}}, matrix{{-5.1, 3.1, 4.1, -6.2}}};
U={matrix {{1, 4/5, -5/22, -2/41}}, matrix {{5/9, -1/3, 7/19, -10}}, matrix {{-13/16, 0, -1/13, 1/49}}, matrix {{-3, -32/7, 47/6, 11/36}}};
J=scoreEquationsFromCovarianceMatrix(G,U);
dim J
degree J

--Again get dim 0, degree 2, factors over QQ.  Interesting!
--I is irreducible.  It is determinantal
--Can we solve for the two solutions symbolically?

matrix {{115/64, 0, -3105/12992, 0},
{0, 5/16, 535/3248, -11/16}, 
{-3105/12992, 535/3248, 962415/2637376, -1177/3248},
{0, -11/16, -1177/3248, 29/16}}

matrix {{115/64, -13/32, -29/64, 47/32},
{-13/32, 5/16, 7/32, -11/16},
{-29/64, 7/32, 27/64, -3623/6496},
{47/32, -11/16, -3623/6496, 29/16}}



--Transcripts

Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases, PrimaryDecomposition, ReesAlgebra,
               TangentCone

i1 : needsPackage "Bertini"
--loading configuration for package "Bertini" from file /Users/davids/Library/Application Support/Macaulay2/init-Bertini.m2

o1 = Bertini

o1 : Package

i2 : needsPackage "GraphicalModels"

o2 = GraphicalModels

o2 : Package

i3 : 
     --G = mixedGraph(digraph {{1,{2,3}},{2,{3}},{3,{4}}},bigraph {{1,2},{2,4}});
     G = mixedGraph(digraph {{1,{2,3}},{2,{3}},{3,{4}}},bigraph {{3,4}});

i4 : --G = mixedGraph(digraph {{1, {2}}, {2, {3}}}, bigraph {{2, 3}});
     
     -- list whose rows are the sampels
     --U = {matrix{{1, 2, -10}}, matrix{{-1, -20, 5}}, matrix{{3, -50, 2}}, matrix{{-1, -4, -30}}};
     U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}

o4 = {| 1 2 1 -1 |, | 2 1 3 0 |, | -1 0 1 1 |, | -5 3 4 -6 |}

o4 : List

i5 : 
     matRtolpR = (M) -> (
         E:=entries(M);    
         return matrix apply(#E, i -> apply(#(E_i), j -> F(E_i_j)))    
     );

i6 : 
     matZZtoQQ = (M) -> (
         E:=entries(M);    
         return matrix apply(#E, i -> apply(#(E_i), j -> (1/1)*E_i_j))    
     );

i7 : 
     sampleCovarianceMatrix = (U) -> (
         n = #U;
         U=apply(#U, i -> if ring(U_i)===ZZ then matZZtoQQ(U_i) else U_i);
         Ubar = matrix{{(1/n)}} * sum(U);
         return ((1/n)*(sum apply(n, i -> (transpose (U#i-Ubar))*(U#i-Ubar))));        
     );

i8 : JacobianMatrixOfRationalFunction = (F) -> (
         f:=numerator(F);
         g:=denominator(F);
         R:=ring(f);
         answer:=diff(vars(R), f) * g - diff(vars(R), g)*f;
         answer=substitute(answer, ring(F));
         matrix({{(1/g)^2}})*answer
     );

i9 : MLEmixedGraph = (G, U) -> (
         V = sampleCovarianceMatrix(U);
         R = gaussianRing G;
         I = gaussianVanishingIdeal R;
         -- Lambda
         L = directedEdgesMatrix R;
         -- d is equal to the number of vertices in G
         d = numRows L;
         -- Omega
         W = bidirectedEdgesMatrix R;
         -- move to a new ring, lpR, which does not have the s variables
         numSvars=lift(d*(d+1)/2,ZZ);
         --lp rings is the ring without the s variables
         lpRvarlist=apply(numgens(R)-numSvars,i->(gens(R))_i);
         KK=coefficientRing(R);
         lpR=KK[lpRvarlist];
         lpRTarget=apply(numgens(R),i-> if i<= numgens(R)-numSvars-1 then (gens(lpR))_i else 0);
         F=map(lpR,R,lpRTarget);
         I=F(I);    
         L = matRtolpR L;
         W = matRtolpR W;
         K = inverse (id_(lpR^d)-L);
         S = (transpose K) * W * K;
         FR = frac(lpR);
         Sinv = inverse substitute(S, FR);
         C1 = trace(Sinv * V)/2;
         f = numerator(C1);
         g = denominator(C1);
         CrazyExp = JacobianMatrixOfRationalFunction(trace(Sinv * V)/2);
         LL = (substitute((jacobian(matrix{{det(S)}})), FR))*matrix{{(-(#U)/(2*det(S)))}} - (transpose CrazyExp);
         ELL=flatten entries(LL);
         J=ideal apply(#ELL, i -> lift(numerator(ELL_i),lpR));
         J = saturate(J, det(S));
         print toString(J) << endl;
         print toString(ring(J)) << endl;
         print toString(dim J) << endl;
         degree J;
         if (dim J == 0) then (
             print concatenate("The dimension of the ideal of likelihood equations is 0 and its degree is ", toString (degree J));
         ) else (
             print concatenate("The dimension of the ideal of likelihood equations is", toString (dim J),
         "and its degree is", toString (degree J), ". The model is not identifiable.");
         );
         print toString(J == radical J) << endl;
         print toString(mingens J) << endl;
     );

i10 : MLEmixedGraph(G,U)
ideal(80*p_(3,4)+39,200*p_(4,4)-271,1760416*p_(3,3)-742363,920*p_(2,2)-203,64*p_(1,1)-115,5*l_(3,4)+2,110026*l_(2,3)-2575,55013*l_(1,3)-600,115*l_(1,2)+26)
lpR
0
The dimension of the ideal of likelihood equations is 0 and its degree is 1
true
matrix {{80*p_(3,4)+39, 200*p_(4,4)-271, 1760416*p_(3,3)-742363, 920*p_(2,2)-203, 64*p_(1,1)-115, 5*l_(3,4)+2, 110026*l_(2,3)-2575, 55013*l_(1,3)-600, 115*l_(1,2)+26}}

i11 : 
      scoreEquationsFromCovarianceMatrix = (G, U) -> (
          V := sampleCovarianceMatrix(U);
          R := gaussianRing G;
          I := gaussianVanishingIdeal R;
          -- Lambda
          L := directedEdgesMatrix R;
          -- d is equal to the number of vertices in G
          d := numRows L;
          -- Omega
          W := bidirectedEdgesMatrix R;
          -- move to a new ring, lpR, which does not have the s variables
          numSvars:=lift(d*(d+1)/2,ZZ);
          --lpR is the ring R without the s variables
          lpRvarlist:=apply(numgens(R)-numSvars,i->(gens(R))_i);
          KK:=coefficientRing(R);
          lpR:=KK[lpRvarlist];
          lpRTarget:=apply(numgens(R),i-> if i<= numgens(R)-numSvars-1 then (gens(lpR))_i else 0);
          F:=map(lpR,R,lpRTarget);
          I=F(I);    
          E:=entries(L);    
          L = matrix apply(#E, i -> apply(#(E_i), j -> F(E_i_j)));   
          E=entries(W);
          W = matrix apply(#E, i -> apply(#(E_i), j -> F(E_i_j)));
          K := inverse (id_(lpR^d)-L);
          S := (transpose K) * W * K;
          FR := frac(lpR);
          Sinv := inverse substitute(S, FR);
          C1 := trace(Sinv * V)/2;
          C1derivative := JacobianMatrixOfRationalFunction(trace(Sinv * V)/2);
          LL := (substitute((jacobian(matrix{{det(S)}})), FR))*matrix{{(-(#U)/(2*det(S)))}} - (transpose C1derivative);
          LL=flatten entries(LL);
          J:=ideal apply(#LL, i -> lift(numerator(LL_i),lpR));
          J = saturate(J, det(S));
          return J
      );

i12 : scoreEquationsFromCovarianceMatrix(G,U)
stdio:87:10:(3):[1]: error: no method found for applying promote to:
     argument 1 :  p    (of class FR)
                    1,1
     argument 2 :  R

i13 : MLEmixedGraph = (G, U) -> (
          V := sampleCovarianceMatrix(U);
          R := gaussianRing G;
          I := gaussianVanishingIdeal R;
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
          I=F(I);    
          L = matRtolpR L;
          W = matRtolpR W;
          K := inverse (id_(lpR^d)-L);
          S := (transpose K) * W * K;
          FR := frac(lpR);
          Sinv := inverse substitute(S, FR);
          C1 := trace(Sinv * V)/2;
          CrazyExp := JacobianMatrixOfRationalFunction(trace(Sinv * V)/2);
          LL := (substitute((jacobian(matrix{{det(S)}})), FR))*matrix{{(-(#U)/(2*det(S)))}} - (transpose CrazyExp);
          ELL:=flatten entries(LL);
          J:=ideal apply(#ELL, i -> lift(numerator(ELL_i),lpR));
          J = saturate(J, det(S));
          print toString(J) << endl;
          print toString(ring(J)) << endl;
          print toString(dim J) << endl;
          degree J;
          if (dim J == 0) then (
              print concatenate("The dimension of the ideal of likelihood equations is 0 and its degree is ", toString (degree J));
          ) else (
              print concatenate("The dimension of the ideal of likelihood equations is", toString (dim J),
          "and its degree is", toString (degree J), ". The model is not identifiable.");
          );
          print toString(J == radical J) << endl;
          print toString(mingens J) << endl;
      );

i14 : MLEmixedGraph(G,U)
stdio:123:10:(3):[1]: error: no method found for applying promote to:
     argument 1 :  p    (of class FR)
                    1,1
     argument 2 :  R

i15 : restart
Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases, PrimaryDecomposition, ReesAlgebra,
               TangentCone

i1 : needsPackage "Bertini"
--loading configuration for package "Bertini" from file /Users/davids/Library/Application Support/Macaulay2/init-Bertini.m2

o1 = Bertini

o1 : Package

i2 : needsPackage "GraphicalModels"

o2 = GraphicalModels

o2 : Package

i3 : --G = mixedGraph(digraph {{1,{2,3}},{2,{3}},{3,{4}}},bigraph {{1,2},{2,4}});
     G = mixedGraph(digraph {{1,{2,3}},{2,{3}},{3,{4}}},bigraph {{3,4}});

i4 : --G = mixedGraph(digraph {{1, {2}}, {2, {3}}}, bigraph {{2, 3}});
     
     -- list whose rows are the sampels
     --U = {matrix{{1, 2, -10}}, matrix{{-1, -20, 5}}, matrix{{3, -50, 2}}, matrix{{-1, -4, -30}}};
     U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}

o4 = {| 1 2 1 -1 |, | 2 1 3 0 |, | -1 0 1 1 |, | -5 3 4 -6 |}

o4 : List

i5 : 
     matRtolpR = (M) -> (
         E:=entries(M);    
         return matrix apply(#E, i -> apply(#(E_i), j -> F(E_i_j)))    
     );

i6 : 
     matZZtoQQ = (M) -> (
         E:=entries(M);    
         return matrix apply(#E, i -> apply(#(E_i), j -> (1/1)*E_i_j))    
     );

i7 : 
     sampleCovarianceMatrix = (U) -> (
         n = #U;
         U=apply(#U, i -> if ring(U_i)===ZZ then matZZtoQQ(U_i) else U_i);
         Ubar = matrix{{(1/n)}} * sum(U);
         return ((1/n)*(sum apply(n, i -> (transpose (U#i-Ubar))*(U#i-Ubar))));        
     );

i8 : 
     JacobianMatrixOfRationalFunction = (F) -> (
         f:=numerator(F);
         g:=denominator(F);
         R:=ring(f);
         answer:=diff(vars(R), f) * g - diff(vars(R), g)*f;
         answer=substitute(answer, FR);
         matrix({{(1/g)^2}})*answer
     );

i9 : 
     MLEmixedGraph = (G, U) -> (
         V := sampleCovarianceMatrix(U);
         R := gaussianRing G;
         I := gaussianVanishingIdeal R;
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
         FR := frac(lpR);
         lpRTarget:=apply(numgens(R),i-> if i<= numgens(R)-numSvars-1 then (gens(lpR))_i else 0);
         F:=map(lpR,R,lpRTarget);
         I=F(I);    
         L = matRtolpR L;
         W = matRtolpR W;
         K := inverse (id_(lpR^d)-L);
         S := (transpose K) * W * K;
         Sinv := inverse substitute(S, FR);
         C1 := trace(Sinv * V)/2;
         CrazyExp := JacobianMatrixOfRationalFunction(trace(Sinv * V)/2);
         LL := (substitute((jacobian(matrix{{det(S)}})), FR))*matrix{{(-(#U)/(2*det(S)))}} - (transpose CrazyExp);
         ELL:=flatten entries(LL);
         J:=ideal apply(#ELL, i -> lift(numerator(ELL_i),lpR));
         J = saturate(J, det(S));
         print toString(J) << endl;
         print toString(ring(J)) << endl;
         print toString(dim J) << endl;
         degree J;
         if (dim J == 0) then (
             print concatenate("The dimension of the ideal of likelihood equations is 0 and its degree is ", toString (degree J));
         ) else (
             print concatenate("The dimension of the ideal of likelihood equations is", toString (dim J),
         "and its degree is", toString (degree J), ". The model is not identifiable.");
         );
         print toString(J == radical J) << endl;
         print toString(mingens J) << endl;
     );

i10 : MLEmixedGraph(G,U)
stdio:40:10:(3):[1]: error: no method found for applying promote to:
     argument 1 :  s    (of class IndexedVariable)
                    1,1
     argument 2 :  QQ[l   , l   , l   , l   , p   , p   , p   , p   , p   , s .
                       1,2   1,3   2,3   3,4   1,1   2,2   3,3   4,4   3,4   1.

i11 : matRtolpR = (M,F) -> (
          E:=entries(M);    
          return matrix apply(#E, i -> apply(#(E_i), j -> F(E_i_j)))    
      );

i12 : MLEmixedGraph = (G, U) -> (
          V := sampleCovarianceMatrix(U);
          R := gaussianRing G;
          I := gaussianVanishingIdeal R;
          -- Lambda
          L := directedEdgesMatrix R;
          -- d is equal to the number of vertices in G
          d := numRows L;
          -- Omega
          W := bidirectedEdgesMatrix R;
          print toString "I got to 1" << endl;
          -- move to a new ring, lpR, which does not have the s variables
          numSvars:=lift(d*(d+1)/2,ZZ);
          --lp rings is the ring without the s variables
          lpRvarlist:=apply(numgens(R)-numSvars,i->(gens(R))_i);
          KK:=coefficientRing(R);
          lpR:=KK[lpRvarlist];
          FR := frac(lpR);
          lpRTarget:=apply(numgens(R),i-> if i<= numgens(R)-numSvars-1 then (gens(lpR))_i else 0);
          F:=map(lpR,R,lpRTarget);
          I=F(I);    
          L = matRtolpR(L,F);
          W = matRtolpR(W,F);
          print toString "I got to 2" << endl;
          K := inverse (id_(lpR^d)-L);
          S := (transpose K) * W * K;
          Sinv := inverse substitute(S, FR);
          C1 := trace(Sinv * V)/2;
          CrazyExp := JacobianMatrixOfRationalFunction(trace(Sinv * V)/2);
          LL := (substitute((jacobian(matrix{{det(S)}})), FR))*matrix{{(-(#U)/(2*det(S)))}} - (transpose CrazyExp);
          ELL:=flatten entries(LL);
          J:=ideal apply(#ELL, i -> lift(numerator(ELL_i),lpR));
          J = saturate(J, det(S));
          print toString(J) << endl;
          print toString(ring(J)) << endl;
          print toString(dim J) << endl;
          degree J;
          if (dim J == 0) then (
              print concatenate("The dimension of the ideal of likelihood equations is 0 and its degree is ", toString (degree J));
          ) else (
              print concatenate("The dimension of the ideal of likelihood equations is", toString (dim J),
          "and its degree is", toString (degree J), ". The model is not identifiable.");
          );
          print toString(J == radical J) << endl;
          print toString(mingens J) << endl;
      );

i13 : MLEmixedGraph(G,U)
stdio:89:10:(3):[1]: error: no method found for applying promote to:
     argument 1 :  s    (of class IndexedVariable)
                    1,1
     argument 2 :  QQ[l   , l   , l   , l   , p   , p   , p   , p   , p   , s .
                       1,2   1,3   2,3   3,4   1,1   2,2   3,3   4,4   3,4   1.

i14 : MLEmixedGraph = (G, U) -> (
          V := sampleCovarianceMatrix(U);
          print toString "I got to 0" << endl;    
          R := gaussianRing G;
          I := gaussianVanishingIdeal R;
          -- Lambda
          L := directedEdgesMatrix R;
          -- d is equal to the number of vertices in G
          d := numRows L;
          -- Omega
          W := bidirectedEdgesMatrix R;
          print toString "I got to 1" << endl;
          -- move to a new ring, lpR, which does not have the s variables
          numSvars:=lift(d*(d+1)/2,ZZ);
          --lp rings is the ring without the s variables
          lpRvarlist:=apply(numgens(R)-numSvars,i->(gens(R))_i);
          KK:=coefficientRing(R);
          lpR:=KK[lpRvarlist];
          FR := frac(lpR);
          lpRTarget:=apply(numgens(R),i-> if i<= numgens(R)-numSvars-1 then (gens(lpR))_i else 0);
          F:=map(lpR,R,lpRTarget);
          I=F(I);    
          L = matRtolpR(L,F);
          W = matRtolpR(W,F);
          print toString "I got to 2" << endl;
          K := inverse (id_(lpR^d)-L);
          S := (transpose K) * W * K;
          Sinv := inverse substitute(S, FR);
          C1 := trace(Sinv * V)/2;
          CrazyExp := JacobianMatrixOfRationalFunction(trace(Sinv * V)/2);
          LL := (substitute((jacobian(matrix{{det(S)}})), FR))*matrix{{(-(#U)/(2*det(S)))}} - (transpose CrazyExp);
          ELL:=flatten entries(LL);
          J:=ideal apply(#ELL, i -> lift(numerator(ELL_i),lpR));
          J = saturate(J, det(S));
          print toString(J) << endl;
          print toString(ring(J)) << endl;
          print toString(dim J) << endl;
          degree J;
          if (dim J == 0) then (
              print concatenate("The dimension of the ideal of likelihood equations is 0 and its degree is ", toString (degree J));
          ) else (
              print concatenate("The dimension of the ideal of likelihood equations is", toString (dim J),
          "and its degree is", toString (degree J), ". The model is not identifiable.");
          );
          print toString(J == radical J) << endl;
          print toString(mingens J) << endl;
      );

i15 : MLEmixedGraph(G,U)
I got to 0
stdio:137:10:(3):[1]: error: no method found for applying promote to:
     argument 1 :  s    (of class IndexedVariable)
                    1,1
     argument 2 :  QQ[l   , l   , l   , l   , p   , p   , p   , p   , p   , s .
                       1,2   1,3   2,3   3,4   1,1   2,2   3,3   4,4   3,4   1.

i16 : restart
Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases, PrimaryDecomposition, ReesAlgebra,
               TangentCone

i1 : needsPackage "Bertini"
--loading configuration for package "Bertini" from file /Users/davids/Library/Application Support/Macaulay2/init-Bertini.m2

o1 = Bertini

o1 : Package

i2 : needsPackage "GraphicalModels"

o2 = GraphicalModels

o2 : Package

i3 : --G = mixedGraph(digraph {{1,{2,3}},{2,{3}},{3,{4}}},bigraph {{1,2},{2,4}});
     G = mixedGraph(digraph {{1,{2,3}},{2,{3}},{3,{4}}},bigraph {{3,4}});

i4 : --G = mixedGraph(digraph {{1, {2}}, {2, {3}}}, bigraph {{2, 3}});
     
     -- list whose rows are the sampels
     --U = {matrix{{1, 2, -10}}, matrix{{-1, -20, 5}}, matrix{{3, -50, 2}}, matrix{{-1, -4, -30}}};
     U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}

o4 = {| 1 2 1 -1 |, | 2 1 3 0 |, | -1 0 1 1 |, | -5 3 4 -6 |}

o4 : List

i5 : 
     matRtolpR = (M,F) -> (
         E:=entries(M);    
         return matrix apply(#E, i -> apply(#(E_i), j -> F(E_i_j)))    
     );

i6 : 
     matZZtoQQ = (M) -> (
         E:=entries(M);    
         return matrix apply(#E, i -> apply(#(E_i), j -> (1/1)*E_i_j))    
     );

i7 : 
     sampleCovarianceMatrix = (U) -> (
         n = #U;
         U=apply(#U, i -> if ring(U_i)===ZZ then matZZtoQQ(U_i) else U_i);
         Ubar = matrix{{(1/n)}} * sum(U);
         return ((1/n)*(sum apply(n, i -> (transpose (U#i-Ubar))*(U#i-Ubar))));        
     );

i8 : 
     JacobianMatrixOfRationalFunction = (F) -> (
         f:=numerator(F);
         g:=denominator(F);
         R:=ring(f);
         answer:=diff(vars(R), f) * g - diff(vars(R), g)*f;
         answer=substitute(answer, FR);
         matrix({{(1/g)^2}})*answer
     );

i9 : MLEmixedGraph = (G, U) -> (
         V := sampleCovarianceMatrix(U);
         print toString "I got to 0" << endl;    
         R := gaussianRing G;
         I := gaussianVanishingIdeal R;
         -- Lambda
         L := directedEdgesMatrix R;
         -- d is equal to the number of vertices in G
         d := numRows L;
         -- Omega
         W := bidirectedEdgesMatrix R;
         print toString "I got to 1" << endl;
         -- move to a new ring, lpR, which does not have the s variables
         numSvars:=lift(d*(d+1)/2,ZZ);
         --lp rings is the ring without the s variables
         lpRvarlist:=apply(numgens(R)-numSvars,i->(gens(R))_i);
         KK:=coefficientRing(R);
         lpR:=KK[lpRvarlist];
         FR := frac(lpR);
         lpRTarget:=apply(numgens(R),i-> if i<= numgens(R)-numSvars-1 then (gens(lpR))_i else 0);
         F:=map(lpR,R,lpRTarget);
         I=F(I);    
         L = matRtolpR(L,F);
         W = matRtolpR(W,F);
         print toString "I got to 2" << endl;
         K := inverse (id_(lpR^d)-L);
         S := (transpose K) * W * K;
         Sinv := inverse substitute(S, FR);
         C1 := trace(Sinv * V)/2;
         CrazyExp := JacobianMatrixOfRationalFunction(trace(Sinv * V)/2);
         LL := (substitute((jacobian(matrix{{det(S)}})), FR))*matrix{{(-(#U)/(2*det(S)))}} - (transpose CrazyExp);
         ELL:=flatten entries(LL);
         J:=ideal apply(#ELL, i -> lift(numerator(ELL_i),lpR));
         J = saturate(J, det(S));
         print toString(J) << endl;
         print toString(ring(J)) << endl;
         print toString(dim J) << endl;
         degree J;
         if (dim J == 0) then (
             print concatenate("The dimension of the ideal of likelihood equations is 0 and its degree is ", toString (degree J));
         ) else (
             print concatenate("The dimension of the ideal of likelihood equations is", toString (dim J),
         "and its degree is", toString (degree J), ". The model is not identifiable.");
         );
         print toString(J == radical J) << endl;
         print toString(mingens J) << endl;
     );

i10 : MLEmixedGraph(G,U)
I got to 0
stdio:40:10:(3):[1]: error: no method found for applying promote to:
     argument 1 :  s    (of class IndexedVariable)
                    1,1
     argument 2 :  QQ[l   , l   , l   , l   , p   , p   , p   , p   , p   , s .
                       1,2   1,3   2,3   3,4   1,1   2,2   3,3   4,4   3,4   1.

i11 : gaussianRing(G)

o11 = QQ[l   , l   , l   , l   , p   , p   , p   , p   , p   , s   , s   , s   , s   , s   , s   , s   , s   , s   , s   ]
          1,2   1,3   2,3   3,4   1,1   2,2   3,3   4,4   3,4   1,1   1,2   1,3   1,4   2,2   2,3   2,4   3,3   3,4   4,4

o11 : PolynomialRing

i12 : gaussianVanishingIdeal R
stdio:85:1:(3): error: no method found for applying gaussianVanishingIdeal to:
     argument   :  R (of class Symbol)

i13 : gaussianVanishingIdeal o11
stdio:86:1:(3): error: no method found for applying promote to:
     argument 1 :  s    (of class IndexedVariable)
                    1,1
     argument 2 :  QQ[l   , l   , l   , l   , p   , p   , p   , p   , p   , s .
                       1,2   1,3   2,3   3,4   1,1   2,2   3,3   4,4   3,4   1.

i14 : restart
Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases, PrimaryDecomposition, ReesAlgebra,
               TangentCone

i1 : needsPackage "Bertini"
--loading configuration for package "Bertini" from file /Users/davids/Library/Application Support/Macaulay2/init-Bertini.m2

o1 = Bertini

o1 : Package

i2 : needsPackage "GraphicalModels"

o2 = GraphicalModels

o2 : Package

i3 : G = mixedGraph(digraph {{1,{2,3}},{2,{3}},{3,{4}}},bigraph {{3,4}});

i4 : R=gaussianRing(G)

o4 = R

o4 : PolynomialRing

i5 : I=gaussianVanishingIdeal R

o5 = ideal(s   s    - s   s   )
            1,4 2,3    1,3 2,4

o5 : Ideal of R

i6 : restart
Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases, PrimaryDecomposition, ReesAlgebra,
               TangentCone

i1 : needsPackage "Bertini"
--loading configuration for package "Bertini" from file /Users/davids/Library/Application Support/Macaulay2/init-Bertini.m2

o1 = Bertini

o1 : Package

i2 : needsPackage "GraphicalModels"

o2 = GraphicalModels

o2 : Package

i3 : 
     G = mixedGraph(digraph {{1,{2,3}},{2,{3}},{3,{4}}},bigraph {{3,4}});

i4 : 
     U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}

o4 = {| 1 2 1 -1 |, | 2 1 3 0 |, | -1 0 1 1 |, | -5 3 4 -6 |}

o4 : List

i5 : 
     matRtolpR = (M) -> (
         E:=entries(M);    
         return matrix apply(#E, i -> apply(#(E_i), j -> F(E_i_j)))    
     );

i6 : 
     matZZtoQQ = (M) -> (
         E:=entries(M);    
         return matrix apply(#E, i -> apply(#(E_i), j -> (1/1)*E_i_j))    
     );

i7 : 
     sampleCovarianceMatrix = (U) -> (
         n:= #U;
         U=apply(#U, i -> if ring(U_i)===ZZ then matZZtoQQ(U_i) else U_i);
         Ubar := matrix{{(1/n)}} * sum(U);
         return ((1/n)*(sum apply(n, i -> (transpose (U#i-Ubar))*(U#i-Ubar))));        
     );

i8 : JacobianMatrixOfRationalFunction = (F) -> (
         f:=numerator(F);
         g:=denominator(F);
         R:=ring(f);
         answer:=diff(vars(R), f) * g - diff(vars(R), g)*f;
         answer=substitute(answer, ring(F));
         matrix({{(1/g)^2}})*answer
     );

i9 : MLEmixedGraph = (G, U) -> (
         V := sampleCovarianceMatrix(U);
         print toString "I got to 0" << endl;    
         R := gaussianRing(G);
         I := gaussianVanishingIdeal(R);
         print toString "I got to 1" << endl;    
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
         FR := frac(lpR);
         lpRTarget:=apply(numgens(R),i-> if i<= numgens(R)-numSvars-1 then (gens(lpR))_i else 0);
         F:=map(lpR,R,lpRTarget);
         I=F(I);    
         L = matRtolpR(L,F);
         W = matRtolpR(W,F);
         print toString "I got to 2" << endl;
         K := inverse (id_(lpR^d)-L);
         S := (transpose K) * W * K;
         Sinv := inverse substitute(S, FR);
         C1 := trace(Sinv * V)/2;
         CrazyExp := JacobianMatrixOfRationalFunction(trace(Sinv * V)/2);
         LL := (substitute((jacobian(matrix{{det(S)}})), FR))*matrix{{(-(#U)/(2*det(S)))}} - (transpose CrazyExp);
         ELL:=flatten entries(LL);
         J:=ideal apply(#ELL, i -> lift(numerator(ELL_i),lpR));
         J = saturate(J, det(S));
         print toString(J) << endl;
         print toString(ring(J)) << endl;
         print toString(dim J) << endl;
         degree J;
         if (dim J == 0) then (
             print concatenate("The dimension of the ideal of likelihood equations is 0 and its degree is ", toString (degree J));
         ) else (
             print concatenate("The dimension of the ideal of likelihood equations is", toString (dim J),
         "and its degree is", toString (degree J), ". The model is not identifiable.");
         );
         print toString(J == radical J) << endl;
         print toString(mingens J) << endl;
     );

i10 : MLEmixedGraph(G,U)
I got to 0
stdio:36:10:(3):[1]: error: no method found for applying promote to:
     argument 1 :  s    (of class IndexedVariable)
                    1,1
     argument 2 :  QQ[l   , l   , l   , l   , p   , p   , p   , p   , p   , s .
                       1,2   1,3   2,3   3,4   1,1   2,2   3,3   4,4   3,4   1.

i11 : R=gaussianRing(G)

o11 = R

o11 : PolynomialRing

i12 : I=gaussianVanishingIdeal(R)

o12 = ideal(s   s    - s   s   )
             1,4 2,3    1,3 2,4

o12 : Ideal of R

i13 : MLEmixedGraph = (G, U) -> (
          V := sampleCovarianceMatrix(U);
          print toString "I got to 0" << endl;    
          R := gaussianRing(G);
      --    I := gaussianVanishingIdeal(R);
          print toString "I got to 1" << endl;    
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
          FR := frac(lpR);
          lpRTarget:=apply(numgens(R),i-> if i<= numgens(R)-numSvars-1 then (gens(lpR))_i else 0);
          F:=map(lpR,R,lpRTarget);
      --    I=F(I);    
          L = matRtolpR(L,F);
          W = matRtolpR(W,F);
          print toString "I got to 2" << endl;
          K := inverse (id_(lpR^d)-L);
          S := (transpose K) * W * K;
          Sinv := inverse substitute(S, FR);
          C1 := trace(Sinv * V)/2;
          CrazyExp := JacobianMatrixOfRationalFunction(trace(Sinv * V)/2);
          LL := (substitute((jacobian(matrix{{det(S)}})), FR))*matrix{{(-(#U)/(2*det(S)))}} - (transpose CrazyExp);
          ELL:=flatten entries(LL);
          J:=ideal apply(#ELL, i -> lift(numerator(ELL_i),lpR));
          J = saturate(J, det(S));
          print toString(J) << endl;
          print toString(ring(J)) << endl;
          print toString(dim J) << endl;
          degree J;
          if (dim J == 0) then (
              print concatenate("The dimension of the ideal of likelihood equations is 0 and its degree is ", toString (degree J));
          ) else (
              print concatenate("The dimension of the ideal of likelihood equations is", toString (dim J),
          "and its degree is", toString (degree J), ". The model is not identifiable.");
          );
          print toString(J == radical J) << endl;
          print toString(mingens J) << endl;
      );

i14 : MLEmixedGraph(G,U)
I got to 0
I got to 1
stdio:8:17:(3):[1]: error: expected 1 argument, but got 2

i15 : matRtolpR = (M,F) -> (
          E:=entries(M);    
          return matrix apply(#E, i -> apply(#(E_i), j -> F(E_i_j)))    
      );

i16 : MLEmixedGraph(G,U)
I got to 0
I got to 1
stdio:89:10:(3):[1]: error: no method found for applying promote to:
     argument 1 :  l    (of class QQ[l   , l   , l   , l   , p .)
                    1,2               1,2   1,3   2,3   3,4   1.
     argument 2 :  R

i17 : restart
Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases, PrimaryDecomposition, ReesAlgebra,
               TangentCone

i1 : needsPackage "Bertini"
--loading configuration for package "Bertini" from file /Users/davids/Library/Application Support/Macaulay2/init-Bertini.m2

o1 = Bertini

o1 : Package

i2 : needsPackage "GraphicalModels"

o2 = GraphicalModels

o2 : Package

i3 : 
     G = mixedGraph(digraph {{1,{2,3}},{2,{3}},{3,{4}}},bigraph {{3,4}});

i4 : 
     U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}

o4 = {| 1 2 1 -1 |, | 2 1 3 0 |, | -1 0 1 1 |, | -5 3 4 -6 |}

o4 : List

i5 : 
     matRtolpR = (M,F) -> (
         E:=entries(M);    
         return matrix apply(#E, i -> apply(#(E_i), j -> F(E_i_j)))    
     );

i6 : 
     matZZtoQQ = (M) -> (
         E:=entries(M);    
         return matrix apply(#E, i -> apply(#(E_i), j -> (1/1)*E_i_j))    
     );

i7 : 
     sampleCovarianceMatrix = (U) -> (
         n:= #U;
         U=apply(#U, i -> if ring(U_i)===ZZ then matZZtoQQ(U_i) else U_i);
         Ubar := matrix{{(1/n)}} * sum(U);
         return ((1/n)*(sum apply(n, i -> (transpose (U#i-Ubar))*(U#i-Ubar))));        
     );

i8 : 
     JacobianMatrixOfRationalFunction = (F) -> (
         f:=numerator(F);
         g:=denominator(F);
         R:=ring(f);
         answer:=diff(vars(R), f) * g - diff(vars(R), g)*f;
         answer=substitute(answer, ring(F));
         matrix({{(1/g)^2}})*answer
     );

i9 : V = sampleCovarianceMatrix(U);

              4        4
o9 : Matrix QQ  <--- QQ

i10 : R = gaussianRing(G);

i11 : L = directedEdgesMatrix R;

              4       4
o11 : Matrix R  <--- R

i12 : d = numRows L;

i13 : W = bidirectedEdgesMatrix R;

              4       4
o13 : Matrix R  <--- R

i14 : numSvars=lift(d*(d+1)/2,ZZ);

i15 : lpRvarlist=apply(numgens(R)-numSvars,i->(gens(R))_i);

i16 : KK=coefficientRing(R);

i17 : lpR=KK[lpRvarlist];

i18 : lpRTarget=apply(numgens(R),i-> if i<= numgens(R)-numSvars-1 then (gens(lpR))_i else 0);

i19 : F=map(lpR,R,lpRTarget);

o19 : RingMap lpR <--- R

i20 : L = matRtolpR(L,F);

                4         4
o20 : Matrix lpR  <--- lpR

i21 : W = matRtolpR(W,F);

                4         4
o21 : Matrix lpR  <--- lpR

i22 : FR = frac(lpR);

i23 : restart
Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases, PrimaryDecomposition, ReesAlgebra,
               TangentCone

i1 : needsPackage "Bertini"
--loading configuration for package "Bertini" from file /Users/davids/Library/Application Support/Macaulay2/init-Bertini.m2

o1 = Bertini

o1 : Package

i2 : needsPackage "GraphicalModels"

o2 = GraphicalModels

o2 : Package

i3 : 
     G = mixedGraph(digraph {{1,{2,3}},{2,{3}},{3,{4}}},bigraph {{3,4}});

i4 : 
     U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}

o4 = {| 1 2 1 -1 |, | 2 1 3 0 |, | -1 0 1 1 |, | -5 3 4 -6 |}

o4 : List

i5 : 
     matRtolpR = (M,F) -> (
         E:=entries(M);    
         return matrix apply(#E, i -> apply(#(E_i), j -> F(E_i_j)))    
     );

i6 : 
     matZZtoQQ = (M) -> (
         E:=entries(M);    
         return matrix apply(#E, i -> apply(#(E_i), j -> (1/1)*E_i_j))    
     );

i7 : 
     sampleCovarianceMatrix = (U) -> (
         n:= #U;
         U=apply(#U, i -> if ring(U_i)===ZZ then matZZtoQQ(U_i) else U_i);
         Ubar := matrix{{(1/n)}} * sum(U);
         return ((1/n)*(sum apply(n, i -> (transpose (U#i-Ubar))*(U#i-Ubar))));        
     );

i8 : 
     JacobianMatrixOfRationalFunction = (F) -> (
         f:=numerator(F);
         g:=denominator(F);
         R:=ring(f);
         answer:=diff(vars(R), f) * g - diff(vars(R), g)*f;
         answer=substitute(answer, ring(F));
         matrix({{(1/g)^2}})*answer
     );

i9 : MLEmixedGraph = (G, U) -> (
         V := sampleCovarianceMatrix(U);
         print toString "I got to 0" << endl;    
         R := gaussianRing(G);
     --    I := gaussianVanishingIdeal(R);
         print toString "I got to 1" << endl;    
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
     --    I=F(I);    
         L = matRtolpR(L,F);
         W = matRtolpR(W,F);
         print toString "I got to 2" << endl;
         FR := frac(lpR);
         K := inverse (id_(lpR^d)-L);
         S := (transpose K) * W * K;
         Sinv := inverse substitute(S, FR);
         print toString "I got to 3" << endl;    
         C1 := trace(Sinv * V)/2;
         CrazyExp := JacobianMatrixOfRationalFunction(trace(Sinv * V)/2);
         LL := (substitute((jacobian(matrix{{det(S)}})), FR))*matrix{{(-(#U)/(2*det(S)))}} - (transpose CrazyExp);
         ELL:=flatten entries(LL);
         J:=ideal apply(#ELL, i -> lift(numerator(ELL_i),lpR));
         J = saturate(J, det(S));
         print toString(J) << endl;
         print toString(ring(J)) << endl;
         print toString(dim J) << endl;
         degree J;
         if (dim J == 0) then (
             print concatenate("The dimension of the ideal of likelihood equations is 0 and its degree is ", toString (degree J));
         ) else (
             print concatenate("The dimension of the ideal of likelihood equations is", toString (dim J),
         "and its degree is", toString (degree J), ". The model is not identifiable.");
         );
         print toString(J == radical J) << endl;
         print toString(mingens J) << endl;
     );

i10 : MLEmixedGraph(G,U)
I got to 0
I got to 1
stdio:40:10:(3):[1]: error: no method found for applying promote to:
     argument 1 :  l    (of class IndexedVariable)
                    1,2
     argument 2 :  QQ[l   , l   , l   , l   , p   , p   , p   , p   , p   , s .
                       1,2   1,3   2,3   3,4   1,1   2,2   3,3   4,4   3,4   1.

i11 : MLEmixedGraph = (G, U) -> (
          V := sampleCovarianceMatrix(U);
          print toString "I got to 0" << endl;    
          R := gaussianRing(G);
      --    I := gaussianVanishingIdeal(R);
          print toString "I got to 1" << endl;    
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
      print toString "I got to 2" << endl;    
          KK:=coefficientRing(R);
          lpR:=KK[lpRvarlist];
          lpRTarget:=apply(numgens(R),i-> if i<= numgens(R)-numSvars-1 then (gens(lpR))_i else 0);
          F:=map(lpR,R,lpRTarget);
      --    I=F(I);    
          L = matRtolpR(L,F);
          W = matRtolpR(W,F);
          FR := frac(lpR);
          K := inverse (id_(lpR^d)-L);
          S := (transpose K) * W * K;
          Sinv := inverse substitute(S, FR);
          print toString "I got to 3" << endl;    
          C1 := trace(Sinv * V)/2;
          CrazyExp := JacobianMatrixOfRationalFunction(trace(Sinv * V)/2);
          LL := (substitute((jacobian(matrix{{det(S)}})), FR))*matrix{{(-(#U)/(2*det(S)))}} - (transpose CrazyExp);
          ELL:=flatten entries(LL);
          J:=ideal apply(#ELL, i -> lift(numerator(ELL_i),lpR));
          J = saturate(J, det(S));
          print toString(J) << endl;
          print toString(ring(J)) << endl;
          print toString(dim J) << endl;
          degree J;
          if (dim J == 0) then (
              print concatenate("The dimension of the ideal of likelihood equations is 0 and its degree is ", toString (degree J));
          ) else (
              print concatenate("The dimension of the ideal of likelihood equations is", toString (dim J),
          "and its degree is", toString (degree J), ". The model is not identifiable.");
          );
          print toString(J == radical J) << endl;
          print toString(mingens J) << endl;
      );

i12 : MLEmixedGraph(G,U)
I got to 0
I got to 1
stdio:89:10:(3):[1]: error: no method found for applying promote to:
     argument 1 :  l    (of class IndexedVariable)
                    1,2
     argument 2 :  QQ[l   , l   , l   , l   , p   , p   , p   , p   , p   , s .
                       1,2   1,3   2,3   3,4   1,1   2,2   3,3   4,4   3,4   1.

i13 : restart
Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases, PrimaryDecomposition, ReesAlgebra,
               TangentCone

i1 : needsPackage "Bertini"
--loading configuration for package "Bertini" from file /Users/davids/Library/Application Support/Macaulay2/init-Bertini.m2

o1 = Bertini

o1 : Package

i2 : needsPackage "GraphicalModels"

o2 = GraphicalModels

o2 : Package

i3 : 
     G = mixedGraph(digraph {{1,{2,3}},{2,{3}},{3,{4}}},bigraph {{3,4}});

i4 : 
     U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}

o4 = {| 1 2 1 -1 |, | 2 1 3 0 |, | -1 0 1 1 |, | -5 3 4 -6 |}

o4 : List

i5 : 
     matRtolpR = (M,F) -> (
         E:=entries(M);    
         return matrix apply(#E, i -> apply(#(E_i), j -> F(E_i_j)))    
     );

i6 : 
     matZZtoQQ = (M) -> (
         E:=entries(M);    
         return matrix apply(#E, i -> apply(#(E_i), j -> (1/1)*E_i_j))    
     );

i7 : 
     sampleCovarianceMatrix = (U) -> (
         n:= #U;
         U=apply(#U, i -> if ring(U_i)===ZZ then matZZtoQQ(U_i) else U_i);
         Ubar := matrix{{(1/n)}} * sum(U);
         return ((1/n)*(sum apply(n, i -> (transpose (U#i-Ubar))*(U#i-Ubar))));        
     );

i8 : 
     JacobianMatrixOfRationalFunction = (F) -> (
         f:=numerator(F);
         g:=denominator(F);
         R:=ring(f);
         answer:=diff(vars(R), f) * g - diff(vars(R), g)*f;
         answer=substitute(answer, ring(F));
         matrix({{(1/g)^2}})*answer
     );

i9 : MLEmixedGraph = (G, U) -> (
         V := sampleCovarianceMatrix(U);
         print toString "I got to 0" << endl;    
         R := gaussianRing(G);
     --    I := gaussianVanishingIdeal(R);
         print toString "I got to 1" << endl;    
         -- Lambda
         L := directedEdgesMatrix R;
         print toString "I got to 2" << endl;   
         -- d is equal to the number of vertices in G
         d := numRows L;
         -- Omega
         W := bidirectedEdgesMatrix R; return {L,W})

o9 = MLEmixedGraph

o9 : FunctionClosure

i10 : MLEmixedGraph(G,U)
I got to 0
I got to 1
stdio:40:10:(3):[1]: error: no method found for applying promote to:
     argument 1 :  l    (of class IndexedVariable)
                    1,2
     argument 2 :  QQ[l   , l   , l   , l   , p   , p   , p   , p   , p   , s .
                       1,2   1,3   2,3   3,4   1,1   2,2   3,3   4,4   3,4   1.

i11 : restart
Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases, PrimaryDecomposition, ReesAlgebra,
               TangentCone

i1 : needsPackage "Bertini"
--loading configuration for package "Bertini" from file /Users/davids/Library/Application Support/Macaulay2/init-Bertini.m2

o1 = Bertini

o1 : Package

i2 : needsPackage "GraphicalModels"

o2 = GraphicalModels

o2 : Package

i3 : 
     G = mixedGraph(digraph {{1,{2,3}},{2,{3}},{3,{4}}},bigraph {{3,4}});

i4 : 
     U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}

o4 = {| 1 2 1 -1 |, | 2 1 3 0 |, | -1 0 1 1 |, | -5 3 4 -6 |}

o4 : List

i5 : 
     matRtolpR = (M,F) -> (
         E:=entries(M);    
         return matrix apply(#E, i -> apply(#(E_i), j -> F(E_i_j)))    
     );

i6 : 
     matZZtoQQ = (M) -> (
         E:=entries(M);    
         return matrix apply(#E, i -> apply(#(E_i), j -> (1/1)*E_i_j))    
     );

i7 : 
     sampleCovarianceMatrix = (U) -> (
         n:= #U;
         U=apply(#U, i -> if ring(U_i)===ZZ then matZZtoQQ(U_i) else U_i);
         Ubar := matrix{{(1/n)}} * sum(U);
         return ((1/n)*(sum apply(n, i -> (transpose (U#i-Ubar))*(U#i-Ubar))));        
     );

i8 : 
     JacobianMatrixOfRationalFunction = (F) -> (
         f:=numerator(F);
         g:=denominator(F);
         R:=ring(f);
         answer:=diff(vars(R), f) * g - diff(vars(R), g)*f;
         answer=substitute(answer, ring(F));
         matrix({{(1/g)^2}})*answer
     );

i9 : MLEmixedGraph = (G, U) -> (
         V := sampleCovarianceMatrix(U);   
         R := gaussianRing(G); 
         use R;   
         -- Lambda
         L := directedEdgesMatrix R;
         print toString "I got to 2" << endl;   
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
     --    I=F(I);    
         L = matRtolpR(L,F);
         W = matRtolpR(W,F);
         FR := frac(lpR);
         K := inverse (id_(lpR^d)-L);
         S := (transpose K) * W * K;
         Sinv := inverse substitute(S, FR);
         print toString "I got to 3" << endl;    
         C1 := trace(Sinv * V)/2;
         CrazyExp := JacobianMatrixOfRationalFunction(trace(Sinv * V)/2);
         LL := (substitute((jacobian(matrix{{det(S)}})), FR))*matrix{{(-(#U)/(2*det(S)))}} - (transpose CrazyExp);
         ELL:=flatten entries(LL);
         J:=ideal apply(#ELL, i -> lift(numerator(ELL_i),lpR));
         J = saturate(J, det(S));
         print toString(J) << endl;
         print toString(ring(J)) << endl;
         print toString(dim J) << endl;
         degree J;
         if (dim J == 0) then (
             print concatenate("The dimension of the ideal of likelihood equations is 0 and its degree is ", toString (degree J));
         ) else (
             print concatenate("The dimension of the ideal of likelihood equations is", toString (dim J),
         "and its degree is", toString (degree J), ". The model is not identifiable.");
         );
         print toString(J == radical J) << endl;
         print toString(mingens J) << endl;
     );

i10 : 
      MLEmixedGraph(G,U)
I got to 2
I got to 3
ideal(80*p_(3,4)+39,200*p_(4,4)-271,1760416*p_(3,3)-742363,920*p_(2,2)-203,64*p_(1,1)-115,5*l_(3,4)+2,110026*l_(2,3)-2575,55013*l_(1,3)-600,115*l_(1,2)+26)
QQ[l_(1,2), l_(1,3), l_(2,3), l_(3,4), p_(1,1), p_(2,2), p_(3,3), p_(4,4), p_(3,4)]
0
The dimension of the ideal of likelihood equations is 0 and its degree is 1
true
matrix {{80*p_(3,4)+39, 200*p_(4,4)-271, 1760416*p_(3,3)-742363, 920*p_(2,2)-203, 64*p_(1,1)-115, 5*l_(3,4)+2, 110026*l_(2,3)-2575, 55013*l_(1,3)-600, 115*l_(1,2)+26}}

i11 : exit
Fordham-David-Swinarski:~ davids$ pwd
/Users/davids
Fordham-David-Swinarski:~ davids$ cd Desktop/Dropbox/AlgStat/
Fordham-David-Swinarski:AlgStat davids$ M2
Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases, PrimaryDecomposition, ReesAlgebra,
               TangentCone

i1 : installPackage("GraphicalModelsMLE")
--loading configuration for package "Bertini" from file /Users/davids/Library/Application Support/Macaulay2/init-Bertini.m2
GraphicalModelsMLE.m2:196:1:(3):[10]: error: documentation key for 'scoreEquationsFromCovarianceMatrix(MixedGraph,List)' encountered, but no method installed
GraphicalModelsMLE.m2:196:1:(3):[10]: --entering debugger (type help to see debugger commands)
GraphicalModelsMLE.m2:196:1-196:1: --source code:
doc /// 

ii2 : restart
Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases, PrimaryDecomposition, ReesAlgebra,
               TangentCone

i1 : installPackage("GraphicalModelsMLE")
--loading configuration for package "Bertini" from file /Users/davids/Library/Application Support/Macaulay2/init-Bertini.m2
warning: 'JacobianMatrixOfRationalFunction' redefined
: here is the first use of 'JacobianMatrixOfRationalFunction'
       : debug with expression   debug 4103   or with command line option   --debug 4103
GraphicalModelsMLE.m2:87:53:(3):[10]: error: expected left hand parameter to be a function, type, or a hash table
GraphicalModelsMLE.m2:87:53:(3):[10]: --entering debugger (type help to see debugger commands)
GraphicalModelsMLE.m2:87:1-118:12: --source code:
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

ii2 : restart
Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases, PrimaryDecomposition, ReesAlgebra,
               TangentCone

i1 : installPackage("GraphicalModelsMLE")
--loading configuration for package "Bertini" from file /Users/davids/Library/Application Support/Macaulay2/init-Bertini.m2
--installing package GraphicalModelsMLE in ../../../Library/Application Support/Macaulay2/local/
--using package sources found in /Users/davids/Desktop/Dropbox/AlgStat/
--storing raw documentation in ../../../Library/Application Support/Macaulay2/local/lib/Macaulay2/x86_64-MacOS-10.8/GraphicalModelsMLE/cache/rawdocumentation-dcba-8.db
--running tests that are functions
--making example result files in ../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/
--making example results for scoreEquationsFromCovarianceMatrix in file ../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/_score__Equations__From__Covariance__Matrix.out
ulimit -t 700; ulimit -m 850000; ulimit -v 850000; ulimit -s 8192; cd /var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/1-rundir/; M2 --silent --print-width 77 --stop --int --no-readline -e 'loadPackage("GraphicalModelsMLE", Reload => true, FileName => "/Users/davids/Desktop/Dropbox/AlgStat/GraphicalModelsMLE.m2")' <"/var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/0_score__Equations__From__Covariance__Matrix.m2" >>"/Users/davids/Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/_score__Equations__From__Covariance__Matrix.errors" 2>&1
../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/_score__Equations__From__Covariance__Matrix.errors:0:1: (output file) error: Macaulay2 exited with return code 256
/var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/0_score__Equations__From__Covariance__Matrix.m2:0:1: (input file)
M2: *** [check] Error 1
stdio:1:1:(3): error: 1 error(s) occurred running examples for package GraphicalModelsMLE

i2 : uninstallPackage("GraphicalModelsMLE")

i3 : restart
Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases, PrimaryDecomposition, ReesAlgebra,
               TangentCone

i1 : installPackage("GraphicalModelsMLE")
--loading configuration for package "Bertini" from file /Users/davids/Library/Application Support/Macaulay2/init-Bertini.m2
--installing package GraphicalModelsMLE in ../../../Library/Application Support/Macaulay2/local/
--using package sources found in /Users/davids/Desktop/Dropbox/AlgStat/
--storing raw documentation in ../../../Library/Application Support/Macaulay2/local/lib/Macaulay2/x86_64-MacOS-10.8/GraphicalModelsMLE/cache/rawdocumentation-dcba-8.db
--running tests that are functions
--making example result files in ../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/
--making example results for sampleCovarianceMatrix in file ../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/_sample__Covariance__Matrix.out
ulimit -t 700; ulimit -m 850000; ulimit -v 850000; ulimit -s 8192; cd /var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/1-rundir/; M2 --silent --print-width 77 --stop --int --no-readline -e 'loadPackage("GraphicalModelsMLE", Reload => true, FileName => "/Users/davids/Desktop/Dropbox/AlgStat/GraphicalModelsMLE.m2")' <"/var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/0_sample__Covariance__Matrix.m2" >>"/Users/davids/Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/_sample__Covariance__Matrix.errors" 2>&1
--making example results for JacobianMatrixOfRationalFunction in file ../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/___Jacobian__Matrix__Of__Rational__Function.out
ulimit -t 700; ulimit -m 850000; ulimit -v 850000; ulimit -s 8192; cd /var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/2-rundir/; M2 --silent --print-width 77 --stop --int --no-readline -e 'loadPackage("GraphicalModelsMLE", Reload => true, FileName => "/Users/davids/Desktop/Dropbox/AlgStat/GraphicalModelsMLE.m2")' <"/var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/0___Jacobian__Matrix__Of__Rational__Function.m2" >>"/Users/davids/Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/___Jacobian__Matrix__Of__Rational__Function.errors" 2>&1
--making example results for GraphicalModelsMLE in file ../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/___Graphical__Models__M__L__E.out
ulimit -t 700; ulimit -m 850000; ulimit -v 850000; ulimit -s 8192; cd /var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/3-rundir/; M2 --silent --print-width 77 --stop --int --no-readline -e 'loadPackage("GraphicalModelsMLE", Reload => true, FileName => "/Users/davids/Desktop/Dropbox/AlgStat/GraphicalModelsMLE.m2")' <"/var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/0___Graphical__Models__M__L__E.m2" >>"/Users/davids/Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/___Graphical__Models__M__L__E.errors" 2>&1
--processing documentation nodes...
--assembling table of contents
--making info file in ../../../Library/Application Support/Macaulay2/local/share/info/GraphicalModelsMLE.info.tmp
--completed info file moved to ../../../Library/Application Support/Macaulay2/local/share/info/GraphicalModelsMLE.info
--making empty html pages in ../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/html/
--making html pages in ../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/html/
--making 'GraphicalModelsMLE : Index' in ../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/html/master.html
--making  GraphicalModelsMLE : Table of Contents' in ../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/html/toc.html
--file created: ../../../Library/Application Support/Macaulay2/local/lib/Macaulay2/x86_64-MacOS-10.8/GraphicalModelsMLE/.installed
--installed package GraphicalModelsMLE in ../../../Library/Application Support/Macaulay2/local/

i2 : G = mixedGraph(digraph {{1,{2,3}},{2,{3}},{3,{4}}},bigraph {{3,4}})
stdio:2:16:(3): error: no method for adjacent objects:
--            digraph (of class Symbol)
--    SPACE   {{1, {2, 3}}, {2, {3}}, {3, {4}}} (of class List)

i3 : restart
Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases, PrimaryDecomposition, ReesAlgebra,
               TangentCone

i1 : loadPackage("GraphicalModelsMLE")
--loading configuration for package "Bertini" from file /Users/davids/Library/Application Support/Macaulay2/init-Bertini.m2

o1 = GraphicalModelsMLE

o1 : Package

i2 : G = mixedGraph(digraph {{1,{2,3}},{2,{3}},{3,{4}}},bigraph {{3,4}})
stdio:2:16:(3): error: no method for adjacent objects:
--            digraph (of class Symbol)
--    SPACE   {{1, {2, 3}}, {2, {3}}, {3, {4}}} (of class List)

i3 : uninstallPackage("GraphicalModelsMLE")

i4 : restart
Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases, PrimaryDecomposition, ReesAlgebra,
               TangentCone

i1 : installPackage("GraphicalModelsMLE")
--loading configuration for package "Bertini" from file /Users/davids/Library/Application Support/Macaulay2/init-Bertini.m2
--installing package GraphicalModelsMLE in ../../../Library/Application Support/Macaulay2/local/
--using package sources found in /Users/davids/Desktop/Dropbox/AlgStat/
--storing raw documentation in ../../../Library/Application Support/Macaulay2/local/lib/Macaulay2/x86_64-MacOS-10.8/GraphicalModelsMLE/cache/rawdocumentation-dcba-8.db
--running tests that are functions
--making example result files in ../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/
--making example results for sampleCovarianceMatrix in file ../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/_sample__Covariance__Matrix.out
ulimit -t 700; ulimit -m 850000; ulimit -v 850000; ulimit -s 8192; cd /var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/1-rundir/; M2 --silent --print-width 77 --stop --int --no-readline -e 'loadPackage("GraphicalModelsMLE", Reload => true, FileName => "/Users/davids/Desktop/Dropbox/AlgStat/GraphicalModelsMLE.m2")' <"/var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/0_sample__Covariance__Matrix.m2" >>"/Users/davids/Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/_sample__Covariance__Matrix.errors" 2>&1
--making example results for GraphicalModelsMLE in file ../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/___Graphical__Models__M__L__E.out
ulimit -t 700; ulimit -m 850000; ulimit -v 850000; ulimit -s 8192; cd /var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/2-rundir/; M2 --silent --print-width 77 --stop --int --no-readline -e 'loadPackage("GraphicalModelsMLE", Reload => true, FileName => "/Users/davids/Desktop/Dropbox/AlgStat/GraphicalModelsMLE.m2")' <"/var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/0___Graphical__Models__M__L__E.m2" >>"/Users/davids/Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/___Graphical__Models__M__L__E.errors" 2>&1
--making example results for scoreEquationsFromCovarianceMatrix in file ../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/_score__Equations__From__Covariance__Matrix.out
ulimit -t 700; ulimit -m 850000; ulimit -v 850000; ulimit -s 8192; cd /var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/3-rundir/; M2 --silent --print-width 77 --stop --int --no-readline -e 'loadPackage("GraphicalModelsMLE", Reload => true, FileName => "/Users/davids/Desktop/Dropbox/AlgStat/GraphicalModelsMLE.m2")' <"/var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/0_score__Equations__From__Covariance__Matrix.m2" >>"/Users/davids/Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/_score__Equations__From__Covariance__Matrix.errors" 2>&1
../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/_score__Equations__From__Covariance__Matrix.errors:0:1: (output file) error: Macaulay2 exited with return code 256
/var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/0_score__Equations__From__Covariance__Matrix.m2:0:1: (input file)
M2: *** [check] Error 1
--making example results for JacobianMatrixOfRationalFunction in file ../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/___Jacobian__Matrix__Of__Rational__Function.out
ulimit -t 700; ulimit -m 850000; ulimit -v 850000; ulimit -s 8192; cd /var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/4-rundir/; M2 --silent --print-width 77 --stop --int --no-readline -e 'loadPackage("GraphicalModelsMLE", Reload => true, FileName => "/Users/davids/Desktop/Dropbox/AlgStat/GraphicalModelsMLE.m2")' <"/var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/0___Jacobian__Matrix__Of__Rational__Function.m2" >>"/Users/davids/Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/___Jacobian__Matrix__Of__Rational__Function.errors" 2>&1
stdio:1:1:(3): error: 1 error(s) occurred running examples for package GraphicalModelsMLE

i2 : uninstallPackage("GraphicalModelsMLE")

i3 : restart
Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases, PrimaryDecomposition, ReesAlgebra,
               TangentCone

i1 : installPackage("GraphicalModelsMLE")
--loading configuration for package "Bertini" from file /Users/davids/Library/Application Support/Macaulay2/init-Bertini.m2
--installing package GraphicalModelsMLE in ../../../Library/Application Support/Macaulay2/local/
--using package sources found in /Users/davids/Desktop/Dropbox/AlgStat/
--storing raw documentation in ../../../Library/Application Support/Macaulay2/local/lib/Macaulay2/x86_64-MacOS-10.8/GraphicalModelsMLE/cache/rawdocumentation-dcba-8.db
--running tests that are functions
--making example result files in ../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/
--making example results for sampleCovarianceMatrix in file ../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/_sample__Covariance__Matrix.out
ulimit -t 700; ulimit -m 850000; ulimit -v 850000; ulimit -s 8192; cd /var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/1-rundir/; M2 --silent --print-width 77 --stop --int --no-readline -e 'loadPackage("GraphicalModelsMLE", Reload => true, FileName => "/Users/davids/Desktop/Dropbox/AlgStat/GraphicalModelsMLE.m2")' <"/var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/0_sample__Covariance__Matrix.m2" >>"/Users/davids/Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/_sample__Covariance__Matrix.errors" 2>&1
--making example results for GraphicalModelsMLE in file ../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/___Graphical__Models__M__L__E.out
ulimit -t 700; ulimit -m 850000; ulimit -v 850000; ulimit -s 8192; cd /var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/2-rundir/; M2 --silent --print-width 77 --stop --int --no-readline -e 'loadPackage("GraphicalModelsMLE", Reload => true, FileName => "/Users/davids/Desktop/Dropbox/AlgStat/GraphicalModelsMLE.m2")' <"/var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/0___Graphical__Models__M__L__E.m2" >>"/Users/davids/Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/___Graphical__Models__M__L__E.errors" 2>&1
--making example results for scoreEquationsFromCovarianceMatrix in file ../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/_score__Equations__From__Covariance__Matrix.out
ulimit -t 700; ulimit -m 850000; ulimit -v 850000; ulimit -s 8192; cd /var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/3-rundir/; M2 --silent --print-width 77 --stop --int --no-readline -e 'loadPackage("GraphicalModelsMLE", Reload => true, FileName => "/Users/davids/Desktop/Dropbox/AlgStat/GraphicalModelsMLE.m2")' <"/var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/0_score__Equations__From__Covariance__Matrix.m2" >>"/Users/davids/Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/_score__Equations__From__Covariance__Matrix.errors" 2>&1
--making example results for JacobianMatrixOfRationalFunction in file ../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/___Jacobian__Matrix__Of__Rational__Function.out
ulimit -t 700; ulimit -m 850000; ulimit -v 850000; ulimit -s 8192; cd /var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/4-rundir/; M2 --silent --print-width 77 --stop --int --no-readline -e 'loadPackage("GraphicalModelsMLE", Reload => true, FileName => "/Users/davids/Desktop/Dropbox/AlgStat/GraphicalModelsMLE.m2")' <"/var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/0___Jacobian__Matrix__Of__Rational__Function.m2" >>"/Users/davids/Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/___Jacobian__Matrix__Of__Rational__Function.errors" 2>&1
--processing documentation nodes...
--assembling table of contents
--making info file in ../../../Library/Application Support/Macaulay2/local/share/info/GraphicalModelsMLE.info.tmp
--completed info file moved to ../../../Library/Application Support/Macaulay2/local/share/info/GraphicalModelsMLE.info
--making empty html pages in ../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/html/
--making html pages in ../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/html/
--making 'GraphicalModelsMLE : Index' in ../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/html/master.html
--making  GraphicalModelsMLE : Table of Contents' in ../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/html/toc.html
--file created: ../../../Library/Application Support/Macaulay2/local/lib/Macaulay2/x86_64-MacOS-10.8/GraphicalModelsMLE/.installed
--installed package GraphicalModelsMLE in ../../../Library/Application Support/Macaulay2/local/

i2 : restart
Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases, PrimaryDecomposition, ReesAlgebra,
               TangentCone

i1 : loadPackage("GraphicalModelsMLE")
--loading configuration for package "Bertini" from file /Users/davids/Library/Application Support/Macaulay2/init-Bertini.m2

o1 = GraphicalModelsMLE

o1 : Package

i2 : viewHelp "GraphicalModelsMLE"

i3 : G = mixedGraph(digraph {{1,{2,3}},{2,{3}},{3,{4}}},bigraph {{3,4}})
stdio:3:16:(3): error: no method for adjacent objects:
--            digraph (of class Symbol)
--    SPACE   {{1, {2, 3}}, {2, {3}}, {3, {4}}} (of class List)

i4 : loadPackage("Graphs")
/Applications/Macaulay2-1.6/share/Macaulay2/Graphs.m2:11:1:(3):[7]: error: package Graphs not reloaded; try Reload => true
/Applications/Macaulay2-1.6/share/Macaulay2/Graphs.m2:11:1:(3):[7]: --entering debugger (type help to see debugger commands)
/Applications/Macaulay2-1.6/share/Macaulay2/Graphs.m2:11:1-21:14: --source code:
newPackage("Graphs",
     Authors => {
          {Name => "Amelia Taylor", Email => "originalbrickhouse@gmail.com"},
          {Name => "Augustine O'Keefe", Email => "aokeefe@tulane.edu"}
          },
     ---- Also Doug Torrance.  --- clearly a current author.  Current role of Amelia and Tina?
     ---- Shaowei Lin and Alex Diaz contributed mixedGraph
     DebuggingMode => false,
     Headline => "Data types, visualization, and basic functions for graphs",
     Version => "0.1"
     )

ii5 : break
/Applications/Macaulay2-1.6/share/Macaulay2/Graphs.m2:1021:1:(3):[7]: error: item to be documented comes from another package: Graphs :: Graphs
/Applications/Macaulay2-1.6/share/Macaulay2/Graphs.m2:1021:1:(3):[7]: --entering debugger (type help to see debugger commands)
/Applications/Macaulay2-1.6/share/Macaulay2/Graphs.m2:1021:1-1021:1: --source code:
doc ///

ii6 : restart
Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases, PrimaryDecomposition, ReesAlgebra,
               TangentCone

i1 : loadPackage("Graphs")

o1 = Graphs

o1 : Package

i2 : loadPackage("GraphicalModels")

o2 = GraphicalModels

o2 : Package

i3 : loadPackage("GraphicalModelsMLE")
--loading configuration for package "Bertini" from file /Users/davids/Library/Application Support/Macaulay2/init-Bertini.m2

o3 = GraphicalModelsMLE

o3 : Package

i4 : G = mixedGraph(digraph {{1,{2,3}},{2,{3}},{3,{4}}},bigraph {{3,4}})

o4 = MixedGraph{Bigraph => Bigraph{3 => set {4}}   }
                                   4 => set {3}
                Digraph => Digraph{1 => set {2, 3}}
                                   2 => set {3}
                                   3 => set {4}
                                   4 => set {}
                Graph => Graph{}

o4 : MixedGraph

i5 : U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}

o5 = {| 1 2 1 -1 |, | 2 1 3 0 |, | -1 0 1 1 |, | -5 3 4 -6 |}

o5 : List

i6 : 
     scoreEquationsFromCovarianceMatrix(G,U)

o6 = ideal (80p    + 39, 200p    - 271, 1760416p    - 742363, 920p    - 203, 64p    - 115, 5l    + 2, 110026l    -
               3,4           4,4                3,3               2,2           1,1          3,4             2,3  
     ---------------------------------------------------------------------------------------------------------------
     2575, 55013l    - 600, 115l    + 26)
                 1,3            1,2

o6 : Ideal of QQ[l   , l   , l   , l   , p   , p   , p   , p   , p   ]
                  1,2   1,3   2,3   3,4   1,1   2,2   3,3   4,4   3,4

i7 : print toString o6
ideal(80*p_(3,4)+39,200*p_(4,4)-271,1760416*p_(3,3)-742363,920*p_(2,2)-203,64*p_(1,1)-115,5*l_(3,4)+2,110026*l_(2,3)-2575,55013*l_(1,3)-600,115*l_(1,2)+26)

i8 : restart
Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases, PrimaryDecomposition, ReesAlgebra,
               TangentCone

i1 : uninstallPackage("GraphicalModelsMLE")

i2 : restart
Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases, PrimaryDecomposition, ReesAlgebra,
               TangentCone

i1 : installPackage("GraphicalModelsMLE")
--loading configuration for package "Bertini" from file /Users/davids/Library/Application Support/Macaulay2/init-Bertini.m2
--installing package GraphicalModelsMLE in ../../../Library/Application Support/Macaulay2/local/
--using package sources found in /Users/davids/Desktop/Dropbox/AlgStat/
--storing raw documentation in ../../../Library/Application Support/Macaulay2/local/lib/Macaulay2/x86_64-MacOS-10.8/GraphicalModelsMLE/cache/rawdocumentation-dcba-8.db
--running tests that are functions
--making example result files in ../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/
--making example results for sampleCovarianceMatrix in file ../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/_sample__Covariance__Matrix.out
ulimit -t 700; ulimit -m 850000; ulimit -v 850000; ulimit -s 8192; cd /var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/1-rundir/; M2 --silent --print-width 77 --stop --int --no-readline -e 'loadPackage("GraphicalModelsMLE", Reload => true, FileName => "/Users/davids/Desktop/Dropbox/AlgStat/GraphicalModelsMLE.m2")' <"/var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/0_sample__Covariance__Matrix.m2" >>"/Users/davids/Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/_sample__Covariance__Matrix.errors" 2>&1
--making example results for GraphicalModelsMLE in file ../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/___Graphical__Models__M__L__E.out
ulimit -t 700; ulimit -m 850000; ulimit -v 850000; ulimit -s 8192; cd /var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/2-rundir/; M2 --silent --print-width 77 --stop --int --no-readline -e 'loadPackage("GraphicalModelsMLE", Reload => true, FileName => "/Users/davids/Desktop/Dropbox/AlgStat/GraphicalModelsMLE.m2")' <"/var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/0___Graphical__Models__M__L__E.m2" >>"/Users/davids/Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/___Graphical__Models__M__L__E.errors" 2>&1
--making example results for scoreEquationsFromCovarianceMatrix in file ../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/_score__Equations__From__Covariance__Matrix.out
ulimit -t 700; ulimit -m 850000; ulimit -v 850000; ulimit -s 8192; cd /var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/3-rundir/; M2 --silent --print-width 77 --stop --int --no-readline -e 'loadPackage("GraphicalModelsMLE", Reload => true, FileName => "/Users/davids/Desktop/Dropbox/AlgStat/GraphicalModelsMLE.m2")' <"/var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/0_score__Equations__From__Covariance__Matrix.m2" >>"/Users/davids/Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/_score__Equations__From__Covariance__Matrix.errors" 2>&1
--making example results for JacobianMatrixOfRationalFunction in file ../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/___Jacobian__Matrix__Of__Rational__Function.out
ulimit -t 700; ulimit -m 850000; ulimit -v 850000; ulimit -s 8192; cd /var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/4-rundir/; M2 --silent --print-width 77 --stop --int --no-readline -e 'loadPackage("GraphicalModelsMLE", Reload => true, FileName => "/Users/davids/Desktop/Dropbox/AlgStat/GraphicalModelsMLE.m2")' <"/var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/0___Jacobian__Matrix__Of__Rational__Function.m2" >>"/Users/davids/Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/example-output/___Jacobian__Matrix__Of__Rational__Function.errors" 2>&1
--processing documentation nodes...
--assembling table of contents
--making info file in ../../../Library/Application Support/Macaulay2/local/share/info/GraphicalModelsMLE.info.tmp
--completed info file moved to ../../../Library/Application Support/Macaulay2/local/share/info/GraphicalModelsMLE.info
--making empty html pages in ../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/html/
--making html pages in ../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/html/
--making 'GraphicalModelsMLE : Index' in ../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/html/master.html
--making  GraphicalModelsMLE : Table of Contents' in ../../../Library/Application Support/Macaulay2/local/share/doc/Macaulay2/GraphicalModelsMLE/html/toc.html
--file created: ../../../Library/Application Support/Macaulay2/local/lib/Macaulay2/x86_64-MacOS-10.8/GraphicalModelsMLE/.installed
--installed package GraphicalModelsMLE in ../../../Library/Application Support/Macaulay2/local/

i2 : check "GraphicalModelsMLE"
--running test 0 of package GraphicalModelsMLE on line 243 in file ./GraphicalModelsMLE.m2
--    rerun with: check_0 "GraphicalModelsMLE"
--making test results in file /var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/5.out
ulimit -t 700; ulimit -m 850000; ulimit -v 850000; ulimit -s 8192; cd /var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/6-rundir/; M2 --silent --print-width 77 --stop --int --no-readline -e 'loadPackage("GraphicalModelsMLE", Reload => true, FileName => "/Users/davids/Desktop/Dropbox/AlgStat/GraphicalModelsMLE.m2")' <"/var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/5.m2" >>"/var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/5.tmp" 2>&1
--running test 1 of package GraphicalModelsMLE on line 251 in file ./GraphicalModelsMLE.m2
--    rerun with: check_1 "GraphicalModelsMLE"
--making test results in file /var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/7.out
ulimit -t 700; ulimit -m 850000; ulimit -v 850000; ulimit -s 8192; cd /var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/8-rundir/; M2 --silent --print-width 77 --stop --int --no-readline -e 'loadPackage("GraphicalModelsMLE", Reload => true, FileName => "/Users/davids/Desktop/Dropbox/AlgStat/GraphicalModelsMLE.m2")' <"/var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/7.m2" >>"/var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/7.tmp" 2>&1
--running test 2 of package GraphicalModelsMLE on line 258 in file ./GraphicalModelsMLE.m2
--    rerun with: check_2 "GraphicalModelsMLE"
--making test results in file /var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/9.out
ulimit -t 700; ulimit -m 850000; ulimit -v 850000; ulimit -s 8192; cd /var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/10-rundir/; M2 --silent --print-width 77 --stop --int --no-readline -e 'loadPackage("GraphicalModelsMLE", Reload => true, FileName => "/Users/davids/Desktop/Dropbox/AlgStat/GraphicalModelsMLE.m2")' <"/var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/9.m2" >>"/var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/9.tmp" 2>&1
--running test 3 of package GraphicalModelsMLE on line 265 in file ./GraphicalModelsMLE.m2
--    rerun with: check_3 "GraphicalModelsMLE"
--making test results in file /var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/11.out
ulimit -t 700; ulimit -m 850000; ulimit -v 850000; ulimit -s 8192; cd /var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/12-rundir/; M2 --silent --print-width 77 --stop --int --no-readline -e 'loadPackage("GraphicalModelsMLE", Reload => true, FileName => "/Users/davids/Desktop/Dropbox/AlgStat/GraphicalModelsMLE.m2")' <"/var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/11.m2" >>"/var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/11.tmp" 2>&1
--running test 4 of package GraphicalModelsMLE on line 275 in file ./GraphicalModelsMLE.m2
--    rerun with: check_4 "GraphicalModelsMLE"
--making test results in file /var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/13.out
ulimit -t 700; ulimit -m 850000; ulimit -v 850000; ulimit -s 8192; cd /var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/14-rundir/; M2 --silent --print-width 77 --stop --int --no-readline -e 'loadPackage("GraphicalModelsMLE", Reload => true, FileName => "/Users/davids/Desktop/Dropbox/AlgStat/GraphicalModelsMLE.m2")' <"/var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/13.m2" >>"/var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/13.tmp" 2>&1

i3 : restart
Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases, PrimaryDecomposition, ReesAlgebra,
               TangentCone

i1 : loadPackage("Graphs")

o1 = Graphs

o1 : Package

i2 : loadPackage("GraphicalModels")

o2 = GraphicalModels

o2 : Package

i3 : loadPackage("GraphicalModelsMLE")
--loading configuration for package "Bertini" from file /Users/davids/Library/Application Support/Macaulay2/init-Bertini.m2

o3 = GraphicalModelsMLE

o3 : Package

i4 : viewHelp "GraphicalModelsMLE"

i5 : G = mixedGraph(digraph {{1,{2,3}},{2,{3}},{3,{4}}},bigraph {{3,4}});

i6 : 
     U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}}

o6 = {| 1 2 1 -1 |, | 2 1 3 0 |, | -1 0 1 1 |, | -5 3 4 -6 |}

o6 : List

i7 : J=scoreEquationsFromCovarianceMatrix(G,U)

o7 = ideal (80p    + 39, 200p    - 271, 1760416p    - 742363, 920p    - 203, 64p    - 115, 5l    + 2, 110026l    -
               3,4           4,4                3,3               2,2           1,1          3,4             2,3  
     ---------------------------------------------------------------------------------------------------------------
     2575, 55013l    - 600, 115l    + 26)
                 1,3            1,2

o7 : Ideal of QQ[l   , l   , l   , l   , p   , p   , p   , p   , p   ]
                  1,2   1,3   2,3   3,4   1,1   2,2   3,3   4,4   3,4

i8 : loadPackage("Bertini")
Bertini.m2:2:1:(3):[7]: error: package Bertini not reloaded; try Reload => true
Bertini.m2:2:1:(3):[7]: --entering debugger (type help to see debugger commands)
Bertini.m2:2:1-26:25: --source code:
newPackage(
  "Bertini",
  Version => "0.2", 
  Date => "May 23, 2013",
  Authors => {
    {Name => "Dan Bates",
     Email => "bates@math.colostate.edu",
     HomePage => "http://www.math.colostate.edu/~bates"}, 
    {Name => "Elizabeth Gross"},
    {Name => "Jose Israel Rodriguez",
     Email => "jo.ro@berkeley.edu",
     HomePage => "http://math.berkeley.edu/~jrodrig/"},
    {Name => "Anton Leykin",
     HomePage => "http://www.math.gatech.edu/~leykin"}
  },
  Headline => "Interface to Bertini",
  Configuration => { 
    "path" => "",
    "BERTINIexe"=>"bertini", 
    "keep files" => true
  },
  DebuggingMode => true,
  AuxiliaryFiles => true,
  CacheExampleOutput => true
)

ii9 : break
Bertini/doc.m2:1:1:(3):[11]: error: item to be documented comes from another package: Bertini :: Bertini
Bertini/doc.m2:1:1:(3):[11]: --entering debugger (type help to see debugger commands)
Bertini/doc.m2:1:1-1:1: --source code:
doc ///

ii10 : restart
Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases, PrimaryDecomposition, ReesAlgebra,
               TangentCone

i1 : loadPackage("Graphs")

o1 = Graphs

o1 : Package

i2 : loadPackage("GraphicalModels")

o2 = GraphicalModels

o2 : Package

i3 : loadPackage("Bertini")
--loading configuration for package "Bertini" from file /Users/davids/Library/Application Support/Macaulay2/init-Bertini.m2

o3 = Bertini

o3 : Package

i4 : loadPackage("GraphicalModelsMLE")

o4 = GraphicalModelsMLE

o4 : Package

i5 : restart
Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases, PrimaryDecomposition, ReesAlgebra,
               TangentCone

i1 : loadPackage("Graphs");

i2 : loadPackage("GraphicalModels");

i3 : loadPackage("Bertini");
--loading configuration for package "Bertini" from file /Users/davids/Library/Application Support/Macaulay2/init-Bertini.m2

i4 : loadPackage("GraphicalModelsMLE");

i5 : G = mixedGraph(digraph {{1,{2,3}},{2,{3}},{3,{4}}},bigraph {{3,4}});

i6 : U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}};

i7 : J=scoreEquationsFromCovarianceMatrix(G,U);

o7 : Ideal of QQ[l   , l   , l   , l   , p   , p   , p   , p   , p   ]
                  1,2   1,3   2,3   3,4   1,1   2,2   3,3   4,4   3,4

i8 : dim J

o8 = 0

i9 : degree J

o9 = 1

i10 : J

o10 = ideal (80p    + 39, 200p    - 271, 1760416p    - 742363, 920p    - 203, 64p    - 115, 5l    + 2, 110026l    -
                3,4           4,4                3,3               2,2           1,1          3,4             2,3  
      --------------------------------------------------------------------------------------------------------------
      2575, 55013l    - 600, 115l    + 26)
                  1,3            1,2

o10 : Ideal of QQ[l   , l   , l   , l   , p   , p   , p   , p   , p   ]
                   1,2   1,3   2,3   3,4   1,1   2,2   3,3   4,4   3,4

i11 : mingens J

o11 = | 80p_(3,4)+39 200p_(4,4)-271 1760416p_(3,3)-742363 920p_(2,2)-203 64p_(1,1)-115 5l_(3,4)+2 110026l_(2,3)-2575
      --------------------------------------------------------------------------------------------------------------
      55013l_(1,3)-600 115l_(1,2)+26 |

                                                                       1                                                                9
o11 : Matrix (QQ[l   , l   , l   , l   , p   , p   , p   , p   , p   ])  <--- (QQ[l   , l   , l   , l   , p   , p   , p   , p   , p   ])
                  1,2   1,3   2,3   3,4   1,1   2,2   3,3   4,4   3,4              1,2   1,3   2,3   3,4   1,1   2,2   3,3   4,4   3,4

i12 : entries oo

o12 = {{80p    + 39, 200p    - 271, 1760416p    - 742363, 920p    - 203, 64p    - 115, 5l    + 2, 110026l    - 2575,
           3,4           4,4                3,3               2,2           1,1          3,4             2,3        
      --------------------------------------------------------------------------------------------------------------
      55013l    - 600, 115l    + 26}}
            1,3            1,2

o12 : List

i13 : viewHelp "Bertini"

i14 : gens J

o14 = | 80p_(3,4)+39 200p_(4,4)-271 1760416p_(3,3)-742363 920p_(2,2)-203 64p_(1,1)-115 5l_(3,4)+2 110026l_(2,3)-2575
      --------------------------------------------------------------------------------------------------------------
      55013l_(1,3)-600 115l_(1,2)+26 |

                                                                       1                                                                9
o14 : Matrix (QQ[l   , l   , l   , l   , p   , p   , p   , p   , p   ])  <--- (QQ[l   , l   , l   , l   , p   , p   , p   , p   , p   ])
                  1,2   1,3   2,3   3,4   1,1   2,2   3,3   4,4   3,4              1,2   1,3   2,3   3,4   1,1   2,2   3,3   4,4   3,4

i15 : flatten entries gens J

o15 = {80p    + 39, 200p    - 271, 1760416p    - 742363, 920p    - 203, 64p    - 115, 5l    + 2, 110026l    - 2575,
          3,4           4,4                3,3               2,2           1,1          3,4             2,3        
      --------------------------------------------------------------------------------------------------------------
      55013l    - 600, 115l    + 26}
            1,3            1,2

o15 : List

i16 : restart
Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases, PrimaryDecomposition, ReesAlgebra,
               TangentCone

i1 : loadPackage("Graphs");

i2 : loadPackage("GraphicalModels");

i3 : loadPackage("Bertini");
--loading configuration for package "Bertini" from file /Users/davids/Library/Application Support/Macaulay2/init-Bertini.m2

i4 : loadPackage("GraphicalModelsMLE");

i5 : G = mixedGraph(digraph {{1,{2,3}},{2,{3}},{3,{4}}},bigraph {{3,4}});

i6 : U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}};

i7 : J=scoreEquationsFromCovarianceMatrix(G,U);

o7 : Ideal of QQ[l   , l   , l   , l   , p   , p   , p   , p   , p   ]
                  1,2   1,3   2,3   3,4   1,1   2,2   3,3   4,4   3,4

i8 : dim J

o8 = 0

i9 : degree J

o9 = 1

i10 : L=flatten entries gens J

o10 = {80p    + 39, 200p    - 271, 1760416p    - 742363, 920p    - 203, 64p    - 115, 5l    + 2, 110026l    - 2575,
          3,4           4,4                3,3               2,2           1,1          3,4             2,3        
      --------------------------------------------------------------------------------------------------------------
      55013l    - 600, 115l    + 26}
            1,3            1,2

o10 : List

i11 : bertiniZeroDimSolve(L)
Temporary directory for input and output files:/var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/0

The version of Bertini you have installed on your computer was used for this run. 
Bertini is under ongoing development by D. Bates, J. Hauenstein, A. Sommese, and C. Wampler.

Bertini.m2:596:19:(3):[10]: error: opening input file "/var/folders/93/0ltyt89x74g8qb5qd_68ppnr0000gp/T/M2-20374-0/0/raw_data" failed: No such file or directory
Bertini.m2:596:19:(3):[10]: --entering debugger (type help to see debugger commands)
Bertini.m2:596:19-596:27: --source code:
       l := lines get (dir|"/raw_data"); -- grabs all lines of the file

ii12 : break

i13 : L

o13 = {80p    + 39, 200p    - 271, 1760416p    - 742363, 920p    - 203, 64p    - 115, 5l    + 2, 110026l    - 2575,
          3,4           4,4                3,3               2,2           1,1          3,4             2,3        
      --------------------------------------------------------------------------------------------------------------
      55013l    - 600, 115l    + 26}
            1,3            1,2

o13 : List

i14 : restart
Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases, PrimaryDecomposition, ReesAlgebra,
               TangentCone

i1 : loadPackage("Graphs");

i2 : loadPackage("GraphicalModels");

i3 : loadPackage("Bertini");
--loading configuration for package "Bertini" from file /Users/davids/Library/Application Support/Macaulay2/init-Bertini.m2

i4 : loadPackage("GraphicalModelsMLE");

i5 : G = mixedGraph(digraph {{1,{2}},{2,{4}},{3,{2}}},bigraph {{3,4}});

i6 : U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}};

i7 : R=gaussianRing(G)

o7 = R

o7 : PolynomialRing

i8 : I=gaussianVanishingIdeal(R)

                       2
o8 = ideal (s   , s   s    - s   s   s    + s   s   s    - s   s   s   )
             1,3   1,4 2,3    1,4 2,2 3,3    1,2 2,4 3,3    1,2 2,3 3,4

o8 : Ideal of R

i9 : restart
Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases, PrimaryDecomposition, ReesAlgebra,
               TangentCone

i1 : loadPackage("Graphs");

i2 : loadPackage("GraphicalModels");

i3 : loadPackage("Bertini");
--loading configuration for package "Bertini" from file /Users/davids/Library/Application Support/Macaulay2/init-Bertini.m2

i4 : loadPackage("GraphicalModelsMLE");

i5 : G = mixedGraph(digraph {{1,{2}},{2,{3,4}}},bigraph {{3,4}});

i6 : U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}};

i7 : R=gaussianRing(G)

o7 = R

o7 : PolynomialRing

i8 : I=gaussianVanishingIdeal(R)

o8 = ideal (s   s    - s   s   , s   s    - s   s   , s   s    - s   s   )
             1,4 2,3    1,3 2,4   1,4 2,2    1,2 2,4   1,3 2,2    1,2 2,3

o8 : Ideal of R

i9 : J=scoreEquationsFromCovarianceMatrix(G,U)

o9 = ideal (40p    + 7, 10p    - 3, 160p    - 43, 920p    - 203, 64p    - 115, 10l    - 7, 5l    + 11, 115l    + 26)
               3,4         4,4          3,3           2,2           1,1           2,3        2,4           1,2

o9 : Ideal of QQ[l   , l   , l   , p   , p   , p   , p   , p   ]
                  1,2   2,4   2,3   1,1   2,2   3,3   4,4   3,4

i10 : dim J

o10 = 0

i11 : degree J

o11 = 1

i12 : restart
Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases, PrimaryDecomposition, ReesAlgebra,
               TangentCone

i1 : loadPackage("Graphs");

i2 : loadPackage("GraphicalModels");

i3 : loadPackage("Bertini");
--loading configuration for package "Bertini" from file /Users/davids/Library/Application Support/Macaulay2/init-Bertini.m2

i4 : loadPackage("GraphicalModelsMLE");

i5 : G = mixedGraph(digraph {{1,{2}},{2,{3,4}}},bigraph {{1,3},{2,4}});

i6 : R=gaussianRing(G);

i7 : I=gaussianVanishingIdeal(R)

                                                                        2
o7 = ideal(s   s   s    - s   s   s    - s   s   s    + s   s   s    + s   s    - s   s   s   )
            1,3 1,4 2,2    1,2 1,4 2,3    1,2 1,3 2,4    1,1 2,3 2,4    1,2 3,4    1,1 2,2 3,4

o7 : Ideal of R

i8 : U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}, matrix{{-5, 3, 4, -6}}};

i9 : J=scoreEquationsFromCovarianceMatrix(G,U);

o9 : Ideal of QQ[l   , l   , l   , p   , p   , p   , p   , p   , p   ]
                  1,2   2,4   2,3   1,1   2,2   3,3   4,4   1,3   2,4

i10 : J

                                                                                                                   
o10 = ideal (12992p    + 3105, 1495p    - 1102p    - 897, 126063p    + 146687p   , 359950p    + 50531p    - 101774,
                   1,3              4,4        2,4               3,3          1,3         2,2         2,4          
      --------------------------------------------------------------------------------------------------------------
                                                                                                 2               
      27p    + 203p   , 203l    - 107, 5l    + 16p    + 11, 179975l    + 62192p    + 13182, 9568p    - 2204p    -
         1,1       1,3      2,3          2,4      2,4              1,2         2,4               2,4        2,4  
      --------------------------------------------------------------------------------------------------------------
                                               2
      897, 12992p   p    + 3105p   , 168792064p    - 9641025)
                 1,3 2,4        2,4            1,3

o10 : Ideal of QQ[l   , l   , l   , p   , p   , p   , p   , p   , p   ]
                   1,2   2,4   2,3   1,1   2,2   3,3   4,4   1,3   2,4

i11 : dim J

o11 = 0

i12 : degree J

o12 = 2

i13 : decompose J

o13 = {ideal (184p    + 39, 12992p    + 3105, 2116p    - 939, 2637376p    - 733435, 16p    - 5, 64p    - 115,
                  2,4             1,3              4,4                3,3              2,2         1,1       
      --------------------------------------------------------------------------------------------------------------
      203l    - 107, 23l    + 35, l   ), ideal (52p    - 23, 12992p    + 3105, 338p    - 313, 2637376p    - 733435,
          2,3           2,4        1,2             2,4             1,3             4,4                3,3          
      --------------------------------------------------------------------------------------------------------------
      920p    - 203, 64p    - 115, 203l    - 107, 13l    + 47, 115l    + 26)}
          2,2           1,1            2,3           2,4           1,2

o13 : List

i14 : print toString J
ideal(12992*p_(1,3)+3105,1495*p_(4,4)-1102*p_(2,4)-897,126063*p_(3,3)+146687*p_(1,3),359950*p_(2,2)+50531*p_(2,4)-101774,27*p_(1,1)+203*p_(1,3),203*l_(2,3)-107,5*l_(2,4)+16*p_(2,4)+11,179975*l_(1,2)+62192*p_(2,4)+13182,9568*p_(2,4)^2-2204*p_(2,4)-897,12992*p_(1,3)*p_(2,4)+3105*p_(2,4),168792064*p_(1,3)^2-9641025)

i15 : print toString o13
{ideal(184*p_(2,4)+39,12992*p_(1,3)+3105,2116*p_(4,4)-939,2637376*p_(3,3)-733435,16*p_(2,2)-5,64*p_(1,1)-115,203*l_(2,3)-107,23*l_(2,4)+35,l_(1,2)), ideal(52*p_(2,4)-23,12992*p_(1,3)+3105,338*p_(4,4)-313,2637376*p_(3,3)-733435,920*p_(2,2)-203,64*p_(1,1)-115,203*l_(2,3)-107,13*l_(2,4)+47,115*l_(1,2)+26)}

i16 : debug GraphicalModelsMLE

i17 : V = sampleCovarianceMatrix(U);

               4        4
o17 : Matrix QQ  <--- QQ

i18 :     R = gaussianRing G;

i19 : --    I = gaussianVanishingIdeal R;
          -- Lambda
          L = directedEdgesMatrix R;

              4       4
o19 : Matrix R  <--- R

i20 :     -- d is equal to the number of vertices in G
          d = numRows L;

i21 :     -- Omega
          W = bidirectedEdgesMatrix R;

              4       4
o21 : Matrix R  <--- R

i22 :     -- move to a new ring, lpR, which does not have the s variables
          numSvars=lift(d*(d+1)/2,ZZ);

i23 :     --lp rings is the ring without the s variables
          lpRvarlist=apply(numgens(R)-numSvars,i->(gens(R))_i);

i24 :     KK=coefficientRing(R);

i25 :     lpR=KK[lpRvarlist];

i26 :     lpRTarget=apply(numgens(R),i-> if i<= numgens(R)-numSvars-1 then (gens(lpR))_i else 0);

i27 :     F=map(lpR,R,lpRTarget);

o27 : RingMap lpR <--- R

i28 : --    I=F(I);    
          L = matRtolpR L;
GraphicalModelsMLE.m2:58:19:(3): error: expected 2 arguments but got 1

i29 :     L = matRtolpR(L,F);

                4         4
o29 : Matrix lpR  <--- lpR

i30 :     W = matRtolpR(W,F);

                4         4
o30 : Matrix lpR  <--- lpR

i31 :     L = directedEdgesMatrix R;
stdio:38:9:(3): error: no method found for applying promote to:
     argument 1 :  l    (of class lpR)
                    1,2
     argument 2 :  R

i32 :     -- d is equal to the number of vertices in G
          d = numRows L;

i33 :     -- Omega
          W = bidirectedEdgesMatrix R;
stdio:42:9:(3): error: no method found for applying promote to:
     argument 1 :  p    (of class lpR)
                    1,1
     argument 2 :  R

i34 : decompose J

o34 = {ideal (184p    + 39, 12992p    + 3105, 2116p    - 939, 2637376p    - 733435, 16p    - 5, 64p    - 115,
                  2,4             1,3              4,4                3,3              2,2         1,1       
      --------------------------------------------------------------------------------------------------------------
      203l    - 107, 23l    + 35, l   ), ideal (52p    - 23, 12992p    + 3105, 338p    - 313, 2637376p    - 733435,
          2,3           2,4        1,2             2,4             1,3             4,4                3,3          
      --------------------------------------------------------------------------------------------------------------
      920p    - 203, 64p    - 115, 203l    - 107, 13l    + 47, 115l    + 26)}
          2,2           1,1            2,3           2,4           1,2

o34 : List

i35 : print toString(o34_1)
ideal(52*p_(2,4)-23,12992*p_(1,3)+3105,338*p_(4,4)-313,2637376*p_(3,3)-733435,920*p_(2,2)-203,64*p_(1,1)-115,203*l_(2,3)-107,13*l_(2,4)+47,115*l_(1,2)+26)


Fordham-David-Swinarski:~ davids$ M2
Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases, PrimaryDecomposition, ReesAlgebra,
               TangentCone

i1 : loadPackage("Graphs");

i2 : loadPackage("GraphicalModels");

i3 : loadPackage("Bertini");
--loading configuration for package "Bertini" from file /Users/davids/Library/Application Support/Macaulay2/init-Bertini.m2

i4 : loadPackage("GraphicalModelsMLE");

i5 : G = mixedGraph(digraph {{1,{2}},{2,{3,4}}},bigraph {{1,3},{2,4}});

i6 :     R = gaussianRing G;

i7 : --    I = gaussianVanishingIdeal R;
         -- Lambda
         L = directedEdgesMatrix R;

             4       4
o7 : Matrix R  <--- R

i8 :     -- d is equal to the number of vertices in G
     
         -- Omega
         W = bidirectedEdgesMatrix R;

             4       4
o8 : Matrix R  <--- R

i9 :     K = inverse (id_(R^d)-L);
stdio:14:23:(3): error: no method for binary operator ^ applied to objects:
--            R (of class PolynomialRing)
--      ^     d (of class Symbol)

i10 : restart
Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases, PrimaryDecomposition, ReesAlgebra,
               TangentCone

i1 : loadPackage("Graphs");

i2 : loadPackage("GraphicalModels");

i3 : loadPackage("Bertini");
--loading configuration for package "Bertini" from file /Users/davids/Library/Application Support/Macaulay2/init-Bertini.m2

i4 : loadPackage("GraphicalModelsMLE");

i5 : G = mixedGraph(digraph {{1,{2}},{2,{3,4}}},bigraph {{1,3},{2,4}});

i6 : R=gaussianRing(G);

i7 :    L = directedEdgesMatrix R;

             4       4
o7 : Matrix R  <--- R

i8 :     -- d is equal to the number of vertices in G
     d=numRows(L);

i9 :     -- Omega
         W = bidirectedEdgesMatrix R;

             4       4
o9 : Matrix R  <--- R

i10 :     K = inverse (id_(R^d)-L);

              4       4
o10 : Matrix R  <--- R

i11 :     S = (transpose K) * W * K;

              4       4
o11 : Matrix R  <--- R

i12 : S

o12 = | p_(1,1)                       l_(1,2)p_(1,1)                                       
      | l_(1,2)p_(1,1)                l_(1,2)^2p_(1,1)+p_(2,2)                             
      | l_(1,2)l_(2,3)p_(1,1)+p_(1,3) l_(1,2)^2l_(2,3)p_(1,1)+l_(2,3)p_(2,2)+l_(1,2)p_(1,3)
      | l_(1,2)l_(2,4)p_(1,1)         l_(1,2)^2l_(2,4)p_(1,1)+l_(2,4)p_(2,2)+p_(2,4)       
      --------------------------------------------------------------------------------------------------------------
      l_(1,2)l_(2,3)p_(1,1)+p_(1,3)                                                            
      l_(1,2)^2l_(2,3)p_(1,1)+l_(2,3)p_(2,2)+l_(1,2)p_(1,3)                                    
      l_(1,2)^2l_(2,3)^2p_(1,1)+l_(2,3)^2p_(2,2)+2l_(1,2)l_(2,3)p_(1,3)+p_(3,3)                
      l_(1,2)^2l_(2,4)l_(2,3)p_(1,1)+l_(2,4)l_(2,3)p_(2,2)+l_(1,2)l_(2,4)p_(1,3)+l_(2,3)p_(2,4)
      --------------------------------------------------------------------------------------------------------------
      l_(1,2)l_(2,4)p_(1,1)                                                                     |
      l_(1,2)^2l_(2,4)p_(1,1)+l_(2,4)p_(2,2)+p_(2,4)                                            |
      l_(1,2)^2l_(2,4)l_(2,3)p_(1,1)+l_(2,4)l_(2,3)p_(2,2)+l_(1,2)l_(2,4)p_(1,3)+l_(2,3)p_(2,4) |
      l_(1,2)^2l_(2,4)^2p_(1,1)+l_(2,4)^2p_(2,2)+2l_(2,4)p_(2,4)+p_(4,4)                        |

              4       4
o12 : Matrix R  <--- R

i13 : Sigma=covarianceMatrix(R)

o13 = | s_(1,1) s_(1,2) s_(1,3) s_(1,4) |
      | s_(1,2) s_(2,2) s_(2,3) s_(2,4) |
      | s_(1,3) s_(2,3) s_(3,3) s_(3,4) |
      | s_(1,4) s_(2,4) s_(3,4) s_(4,4) |

              4       4
o13 : Matrix R  <--- R

i14 : A=ideal(Sigma-S)

                                                                                                                
o14 = ideal (- p    + s   , - l   p    + s   , - l   l   p    - p    + s   , - l   l   p    + s   , - l   p    +
                1,1    1,1     1,2 1,1    1,2     1,2 2,3 1,1    1,3    1,3     1,2 2,4 1,1    1,4     1,2 1,1  
      --------------------------------------------------------------------------------------------------------------
               2                         2                                            2                            
      s   , - l   p    - p    + s   , - l   l   p    - l   p    - l   p    + s   , - l   l   p    - l   p    - p   
       1,2     1,2 1,1    2,2    2,2     1,2 2,3 1,1    2,3 2,2    1,2 1,3    2,3     1,2 2,4 1,1    2,4 2,2    2,4
      --------------------------------------------------------------------------------------------------------------
                                               2                                            2   2          2        
      + s   , - l   l   p    - p    + s   , - l   l   p    - l   p    - l   p    + s   , - l   l   p    - l   p    -
         2,4     1,2 2,3 1,1    1,3    1,3     1,2 2,3 1,1    2,3 2,2    1,2 1,3    2,3     1,2 2,3 1,1    2,3 2,2  
      --------------------------------------------------------------------------------------------------------------
                                      2                                                                 
      2l   l   p    - p    + s   , - l   l   l   p    - l   l   p    - l   l   p    - l   p    + s   , -
        1,2 2,3 1,3    3,3    3,3     1,2 2,4 2,3 1,1    2,4 2,3 2,2    1,2 2,4 1,3    2,3 2,4    3,4   
      --------------------------------------------------------------------------------------------------------------
                              2                                        2                                            
      l   l   p    + s   , - l   l   p    - l   p    - p    + s   , - l   l   l   p    - l   l   p    - l   l   p   
       1,2 2,4 1,1    1,4     1,2 2,4 1,1    2,4 2,2    2,4    2,4     1,2 2,4 2,3 1,1    2,4 2,3 2,2    1,2 2,4 1,3
      --------------------------------------------------------------------------------------------------------------
                            2   2          2
      - l   p    + s   , - l   l   p    - l   p    - 2l   p    - p    + s   )
         2,3 2,4    3,4     1,2 2,4 1,1    2,4 2,2     2,4 2,4    4,4    4,4

o14 : Ideal of R

i15 : print toString A
ideal(-p_(1,1)+s_(1,1),-l_(1,2)*p_(1,1)+s_(1,2),-l_(1,2)*l_(2,3)*p_(1,1)-p_(1,3)+s_(1,3),-l_(1,2)*l_(2,4)*p_(1,1)+s_(1,4),-l_(1,2)*p_(1,1)+s_(1,2),-l_(1,2)^2*p_(1,1)-p_(2,2)+s_(2,2),-l_(1,2)^2*l_(2,3)*p_(1,1)-l_(2,3)*p_(2,2)-l_(1,2)*p_(1,3)+s_(2,3),-l_(1,2)^2*l_(2,4)*p_(1,1)-l_(2,4)*p_(2,2)-p_(2,4)+s_(2,4),-l_(1,2)*l_(2,3)*p_(1,1)-p_(1,3)+s_(1,3),-l_(1,2)^2*l_(2,3)*p_(1,1)-l_(2,3)*p_(2,2)-l_(1,2)*p_(1,3)+s_(2,3),-l_(1,2)^2*l_(2,3)^2*p_(1,1)-l_(2,3)^2*p_(2,2)-2*l_(1,2)*l_(2,3)*p_(1,3)-p_(3,3)+s_(3,3),-l_(1,2)^2*l_(2,4)*l_(2,3)*p_(1,1)-l_(2,4)*l_(2,3)*p_(2,2)-l_(1,2)*l_(2,4)*p_(1,3)-l_(2,3)*p_(2,4)+s_(3,4),-l_(1,2)*l_(2,4)*p_(1,1)+s_(1,4),-l_(1,2)^2*l_(2,4)*p_(1,1)-l_(2,4)*p_(2,2)-p_(2,4)+s_(2,4),-l_(1,2)^2*l_(2,4)*l_(2,3)*p_(1,1)-l_(2,4)*l_(2,3)*p_(2,2)-l_(1,2)*l_(2,4)*p_(1,3)-l_(2,3)*p_(2,4)+s_(3,4),-l_(1,2)^2*l_(2,4)^2*p_(1,1)-l_(2,4)^2*p_(2,2)-2*l_(2,4)*p_(2,4)-p_(4,4)+s_(4,4))

i16 : A = ideal(-p_(1,1)+s_(1,1),-l_(1,2)*p_(1,1)+s_(1,2),-l_(1,2)*l_(2,3)*p_(1,1)-p_(1,3)+s_(1,3),-l_(1,2)*l_(2,4)*p_(1,1)+s_(1,4),-l_(1,2)*p_(1,1)+s_(1,2),-l_(1,2)^2*p_(1,1)-p_(2,2)+s_(2,2),-l_(1,2)^2*l_(2,3)*p_(1,1)-l_(2,3)*p_(2,2)-l_(1,2)*p_(1,3)+s_(2,3),-l_(1,2)^2*l_(2,4)*p_(1,1)-l_(2,4)*p_(2,2)-p_(2,4)+s_(2,4),-l_(1,2)*l_(2,3)*p_(1,1)-p_(1,3)+s_(1,3),-l_(1,2)^2*l_(2,3)*p_(1,1)-l_(2,3)*p_(2,2)-l_(1,2)*p_(1,3)+s_(2,3),-l_(1,2)^2*l_(2,3)^2*p_(1,1)-l_(2,3)^2*p_(2,2)-2*l_(1,2)*l_(2,3)*p_(1,3)-p_(3,3)+s_(3,3),-l_(1,2)^2*l_(2,4)*l_(2,3)*p_(1,1)-l_(2,4)*l_(2,3)*p_(2,2)-l_(1,2)*l_(2,4)*p_(1,3)-l_(2,3)*p_(2,4)+s_(3,4),-l_(1,2)*l_(2,4)*p_(1,1)+s_(1,4),-l_(1,2)^2*l_(2,4)*p_(1,1)-l_(2,4)*p_(2,2)-p_(2,4)+s_(2,4),-l_(1,2)^2*l_(2,4)*l_(2,3)*p_(1,1)-l_(2,4)*l_(2,3)*p_(2,2)-l_(1,2)*l_(2,4)*p_(1,3)-l_(2,3)*p_(2,4)+s_(3,4),-l_(1,2)^2*l_(2,4)^2*p_(1,1)-l_(2,4)^2*p_(2,2)-2*l_(2,4)*p_(2,4)-p_(4,4)+s_(4,4))

                                                                                                                
o16 = ideal (- p    + s   , - l   p    + s   , - l   l   p    - p    + s   , - l   l   p    + s   , - l   p    +
                1,1    1,1     1,2 1,1    1,2     1,2 2,3 1,1    1,3    1,3     1,2 2,4 1,1    1,4     1,2 1,1  
      --------------------------------------------------------------------------------------------------------------
               2                         2                                            2                            
      s   , - l   p    - p    + s   , - l   l   p    - l   p    - l   p    + s   , - l   l   p    - l   p    - p   
       1,2     1,2 1,1    2,2    2,2     1,2 2,3 1,1    2,3 2,2    1,2 1,3    2,3     1,2 2,4 1,1    2,4 2,2    2,4
      --------------------------------------------------------------------------------------------------------------
                                               2                                            2   2          2        
      + s   , - l   l   p    - p    + s   , - l   l   p    - l   p    - l   p    + s   , - l   l   p    - l   p    -
         2,4     1,2 2,3 1,1    1,3    1,3     1,2 2,3 1,1    2,3 2,2    1,2 1,3    2,3     1,2 2,3 1,1    2,3 2,2  
      --------------------------------------------------------------------------------------------------------------
                                      2                                                                 
      2l   l   p    - p    + s   , - l   l   l   p    - l   l   p    - l   l   p    - l   p    + s   , -
        1,2 2,3 1,3    3,3    3,3     1,2 2,4 2,3 1,1    2,4 2,3 2,2    1,2 2,4 1,3    2,3 2,4    3,4   
      --------------------------------------------------------------------------------------------------------------
                              2                                        2                                            
      l   l   p    + s   , - l   l   p    - l   p    - p    + s   , - l   l   l   p    - l   l   p    - l   l   p   
       1,2 2,4 1,1    1,4     1,2 2,4 1,1    2,4 2,2    2,4    2,4     1,2 2,4 2,3 1,1    2,4 2,3 2,2    1,2 2,4 1,3
      --------------------------------------------------------------------------------------------------------------
                            2   2          2
      - l   p    + s   , - l   l   p    - l   p    - 2l   p    - p    + s   )
         2,3 2,4    3,4     1,2 2,4 1,1    2,4 2,2     2,4 2,4    4,4    4,4

o16 : Ideal of R

i17 : B=ideal(184*p_(2,4)+39,12992*p_(1,3)+3105,2116*p_(4,4)-939,2637376*p_(3,3)-733435,16*p_(2,2)-5,64*p_(1,1)-115,203*l_(2,3)-107,23*l_(2,4)+35,l_(1,2));

o17 : Ideal of R

i18 : C=A+B;

o18 : Ideal of R

i19 : viewHelp "eliminate"

i20 : gens R

o20 = {l   , l   , l   , p   , p   , p   , p   , p   , p   , s   , s   , s   , s   , s   , s   , s   , s   , s   ,
        1,2   2,4   2,3   1,1   2,2   3,3   4,4   1,3   2,4   1,1   1,2   1,3   1,4   2,2   2,3   2,4   3,3   3,4 
      --------------------------------------------------------------------------------------------------------------
      s   }
       4,4

o20 : List

i21 : print toString oo
{l_(1,2), l_(2,4), l_(2,3), p_(1,1), p_(2,2), p_(3,3), p_(4,4), p_(1,3), p_(2,4), s_(1,1), s_(1,2), s_(1,3), s_(1,4), s_(2,2), s_(2,3), s_(2,4), s_(3,3), s_(3,4), s_(4,4)}

i22 : eliminate({l_(1,2), l_(2,4), l_(2,3), p_(1,1), p_(2,2), p_(3,3), p_(4,4), p_(1,3), p_(2,4)},C)

o22 = ideal (16s    - 29, 3248s    + 1177, 2637376s    - 962415, 16s    + 11, 3248s    - 535, 16s    - 5, s   ,
                4,4            3,4                 3,3              2,4            2,3           2,2       1,4 
      --------------------------------------------------------------------------------------------------------------
      12992s    + 3105, s   , 64s    - 115)
            1,3          1,2     1,1

o22 : Ideal of R

i23 : print toString oo
ideal(16*s_(4,4)-29,3248*s_(3,4)+1177,2637376*s_(3,3)-962415,16*s_(2,4)+11,3248*s_(2,3)-535,16*s_(2,2)-5,s_(1,4),12992*s_(1,3)+3105,s_(1,2),64*s_(1,1)-115)

i24 : flatten entries gens o22

o24 = {16s    - 29, 3248s    + 1177, 2637376s    - 962415, 16s    + 11, 3248s    - 535, 16s    - 5, s   , 12992s   
          4,4            3,4                 3,3              2,4            2,3           2,2       1,4        1,3
      --------------------------------------------------------------------------------------------------------------
      + 3105, s   , 64s    - 115}
               1,2     1,1

o24 : List

i25 : #oo

o25 = 10

i26 : for i from 0 to #o24-1 do print toString(o24_i) << endl
16*s_(4,4)-29
3248*s_(3,4)+1177
2637376*s_(3,3)-962415
16*s_(2,4)+11
3248*s_(2,3)-535
16*s_(2,2)-5
s_(1,4)
12992*s_(1,3)+3105
s_(1,2)
64*s_(1,1)-115

i27 : eig matrix {{115/64, 0, -3105/12992, 0},
      {0, 5/16, 535/3248, -11/16}, 
      {-3105/12992, 535/3248, 962415/2637376, -1177/3248},
      {0, -11/16, -1177/3248, 29/16}}
stdio:29:1:(3): error: no method for adjacent objects:
--            eig (of class Symbol)
--    SPACE   | 115/64      0        -3105/12992    0          | (of class Matrix)
--            | 0           5/16     535/3248       -11/16     |
--            | -3105/12992 535/3248 962415/2637376 -1177/3248 |
--            | 0           -11/16   -1177/3248     29/16      |

i28 : eigenvalues matrix {{115/64, 0, -3105/12992, 0},
      {0, 5/16, 535/3248, -11/16}, 
      {-3105/12992, 535/3248, 962415/2637376, -1177/3248},
      {0, -11/16, -1177/3248, 29/16}}

o28 = {1.82437 }
      {2.17514 }
      {.244516 }
      {.0427617}

o28 : VerticalList

i29 : D=ideal(52*p_(2,4)-23,12992*p_(1,3)+3105,338*p_(4,4)-313,2637376*p_(3,3)-733435,920*p_(2,2)-203,64*p_(1,1)-115,203*l_(2,3)-107,13*l_(2,4)+47,115*l_(1,2)+26);

o29 : Ideal of R

i30 : E=A+D;

o30 : Ideal of R

i31 : eliminate({l_(1,2), l_(2,4), l_(2,3), p_(1,1), p_(2,2), p_(3,3), p_(4,4), p_(1,3), p_(2,4)},E)

o31 = ideal (16s    - 29, 6496s    + 3623, 64s    - 27, 16s    + 11, 32s    - 7, 16s    - 5, 32s    - 47, 64s    +
                4,4            3,4            3,3          2,4          2,3         2,2         1,4          1,3  
      --------------------------------------------------------------------------------------------------------------
      29, 32s    + 13, 64s    - 115)
             1,2          1,1

o31 : Ideal of R

i32 : for i from 0 to #o31-1 do (print toString(o31_i) << endl)
16*s_(4,4)-29
6496*s_(3,4)+3623
64*s_(3,3)-27

i33 : L=flatten entries gens o31

o33 = {16s    - 29, 6496s    + 3623, 64s    - 27, 16s    + 11, 32s    - 7, 16s    - 5, 32s    - 47, 64s    + 29,
          4,4            3,4            3,3          2,4          2,3         2,2         1,4          1,3      
      --------------------------------------------------------------------------------------------------------------
      32s    + 13, 64s    - 115}
         1,2          1,1

o33 : List

i34 : for i from 0 to #o33-1 do (print toString(o33_i) << endl)
16*s_(4,4)-29
6496*s_(3,4)+3623
64*s_(3,3)-27
16*s_(2,4)+11
32*s_(2,3)-7
16*s_(2,2)-5
32*s_(1,4)-47
64*s_(1,3)+29
32*s_(1,2)+13
64*s_(1,1)-115

i35 : eigenvalues matrix {{115/64, -13/32, -29/64, 47/32},
      {-13/32, 5/16, 7/32, -11/16},
      {-29/64, 7/32, 27/64, -3623/6496},
      {47/32, -11/16, -3623/6496, 29/16}}

o35 = {3.63809  }
      {.470419  }
      {.00998306}
      {.225253  }

o35 : VerticalList

Last login: Wed Jan  8 11:15:08 on ttys003
Fordham-David-Swinarski:~ davids$ M2
Macaulay2, version 1.6
with packages: ConwayPolynomials, Elimination, IntegralClosure, LLLBases, PrimaryDecomposition, ReesAlgebra,
               TangentCone

i1 : loadPackage("Graphs");

i2 : loadPackage("GraphicalModels");

i3 : loadPackage("Bertini");
--loading configuration for package "Bertini" from file /Users/davids/Library/Application Support/Macaulay2/init-Bertini.m2

i4 : loadPackage("GraphicalModelsMLE");

i5 : G = mixedGraph(digraph {{1,{2}},{2,{3,4}}},bigraph {{1,3},{2,4}});

i6 : U={matrix {{1, 4/5, -5/22, -2/41}}, matrix {{5/9, -1/3, 7/19, -10}}, matrix {{-13/16, 0, -1/13, 1/49}}, matrix {{-3, -32/7, 47/6, 11/36}}};

i7 : J=scoreEquationsFromCovarianceMatrix(G,U);

o7 : Ideal of QQ[l   , l   , l   , p   , p   , p   , p   , p   , p   ]
                  1,2   2,4   2,3   1,1   2,2   3,3   4,4   1,3   2,4

i8 : dim J

o8 = 0

i9 : degree J

o9 = 2

i10 : mingens J

o10 = | 65778937915392p_(1,3)-8237297521993
      --------------------------------------------------------------------------------------------------------------
      1492493108720121470007687801362112p_(4,4)+10778189343184084993552404883334496p_(2,4)-
      --------------------------------------------------------------------------------------------------------------
      13566522828716826793628943722343233 10671016041188896431946479p_(3,3)-7037216301628946369868232p_(1,3)
      --------------------------------------------------------------------------------------------------------------
      96377412412190476121003629907047848800p_(2,2)+9916360626515293055871246684584203200p_(2,4)-
      --------------------------------------------------------------------------------------------------------------
      21180802427550955603964691063839056711 737511576p_(1,1)-3568735781p_(1,3) 235536561546l_(2,3)+411891124855
      --------------------------------------------------------------------------------------------------------------
      2648308572l_(2,4)+2430086400p_(2,4)+1227190355
      --------------------------------------------------------------------------------------------------------------
      86051261082312925108038955274149865l_(1,2)-11792371307797840476379097810995200p_(2,4)-
      --------------------------------------------------------------------------------------------------------------
      99714594281479456065358281100079568
      --------------------------------------------------------------------------------------------------------------
      2739021610956359952916914110668800p_(2,4)^2+21556378686368169987104809766668992p_(2,4)-
      --------------------------------------------------------------------------------------------------------------
      13566522828716826793628943722343233 65778937915392p_(1,3)p_(2,4)-8237297521993p_(2,4)
      --------------------------------------------------------------------------------------------------------------
      4326868673276995234550513664p_(1,3)^2-67853070465832018318692049 |

                                                                       1                                                                11
o10 : Matrix (QQ[l   , l   , l   , p   , p   , p   , p   , p   , p   ])  <--- (QQ[l   , l   , l   , p   , p   , p   , p   , p   , p   ])
                  1,2   2,4   2,3   1,1   2,2   3,3   4,4   1,3   2,4              1,2   2,4   2,3   1,1   2,2   3,3   4,4   1,3   2,4

i11 : decompose J

o11 = {ideal (4188635785188832521600p    + 35418501250618079506319, 65778937915392p    - 8237297521993,
                                     2,4                                           1,3                 
      --------------------------------------------------------------------------------------------------------------
      74594684272808110287953863236453292800p    - 5233164403743350431937741986520702387491,
                                             4,4                                            
      --------------------------------------------------------------------------------------------------------------
      10651674590383228014742272p    - 879652037703618296233529, 235200p    - 256321, 442368p    - 268057,
                                 3,3                                    2,2                  1,1          
      --------------------------------------------------------------------------------------------------------------
      235536561546l    + 411891124855, 35617651234598916l    - 259855773249341113, l   ), ideal (653917349568p    -
                   2,3                                   2,4                        1,2                       2,4  
      --------------------------------------------------------------------------------------------------------------
      383034921007, 65778937915392p    - 8237297521993, 7363146238382030334336p    - 35783091041145029884987,
                                   1,3                                         4,4                           
      --------------------------------------------------------------------------------------------------------------
      10651674590383228014742272p    - 879652037703618296233529, 2626958600p    - 419001367, 442368p    - 268057,
                                 3,3                                        2,2                     1,1          
      --------------------------------------------------------------------------------------------------------------
      235536561546l    + 411891124855, 2502234756l    + 2504422465, 9381995l    - 11624784)}
                   2,3                            2,4                       1,2

o11 : List

i12 : R=gaussianRing(G);

i13 : I=gaussianVanishingIdeal(R)

                                                                         2
o13 = ideal(s   s   s    - s   s   s    - s   s   s    + s   s   s    + s   s    - s   s   s   )
             1,3 1,4 2,2    1,2 1,4 2,3    1,2 1,3 2,4    1,1 2,3 2,4    1,2 3,4    1,1 2,2 3,4

o13 : Ideal of R

i14 : decompose I

                                                                            2
o14 = {ideal(- s   s   s    + s   s   s    + s   s   s    - s   s   s    - s   s    + s   s   s   )}
                1,3 1,4 2,2    1,2 1,4 2,3    1,2 1,3 2,4    1,1 2,3 2,4    1,2 3,4    1,1 2,2 3,4

o14 : List

i15 : matrix {{115/64, 0, -3105/12992, 0},
      {0, 5/16, 535/3248, -11/16}, 
      {-3105/12992, 535/3248, 962415/2637376, -1177/3248},
      {0, -11/16, -1177/3248, 29/16}}

o15 = | 115/64      0        -3105/12992    0          |
      | 0           5/16     535/3248       -11/16     |
      | -3105/12992 535/3248 962415/2637376 -1177/3248 |
      | 0           -11/16   -1177/3248     29/16      |

               4        4
o15 : Matrix QQ  <--- QQ

i16 : matrix {{115/64, -13/32, -29/64, 47/32},
      {-13/32, 5/16, 7/32, -11/16},
      {-29/64, 7/32, 27/64, -3623/6496},
      {47/32, -11/16, -3623/6496, 29/16}}

o16 = | 115/64 -13/32 -29/64     47/32      |
      | -13/32 5/16   7/32       -11/16     |
      | -29/64 7/32   27/64      -3623/6496 |
      | 47/32  -11/16 -3623/6496 29/16      |

               4        4
o16 : Matrix QQ  <--- QQ