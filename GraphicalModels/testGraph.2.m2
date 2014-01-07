needsPackage "GraphicalModels"

--G = mixedGraph(digraph {{1,{2,3}},{2,{3}},{3,{4}}},bigraph {{1,2},{2,4}});
G = mixedGraph(digraph {{1, {2}}, {2, {3}}}, bigraph {{1, 2}});
-- number of vertices in G
-- d = 3;
-- list whose rows are the sampels
--U = {matrix{{1,2,1,-1}}, matrix{{2,1,3,0}}, matrix{{-1, 0, 1, 1}}};
U = {matrix{{1, 2, 0}}, matrix{{-1, 0, 5/1}}, matrix{{3, 5, 2/1}}, matrix{{-1, -4, 1/1}}};
--U = matrix{{1, 2, 0}, {-1, 0, 5}, {3, 5, 2},{-1, -4, 1.0}};
	
R = gaussianRing G
I = gaussianVanishingIdeal R;
-- Lambda
L = directedEdgesMatrix R
-- Omega
W = bidirectedEdgesMatrix R
d=numRows(L);
numSvars=lift(d*(d+1)/2,ZZ);
lpRvarlist=apply(numgens(R)-numSvars,i->(gens(R))_i);
lpRTarget=apply(numgens(R),i-> if i<= numgens(R)-numSvars-1 then (gens(R))_i else 0);
KK=coefficientRing(R);
lpR=KK[lpRvarlist];
F=map(lpR,R,lpRTarget)



K = inverse (id_(R^d)-L)
--the matrix S is the matrix Sigma in Sturmfels's notes
S = (transpose K) * W * K
FR = frac(R)
Sinv = inverse substitute(S, FR)

{*
-- V = sum of U_i^T * U_i
V = (transpose U#0)*(U#0);
for i from 1 to #U-1 do (V = V + (transpose U#i)*(U#i));
*}

--Take a matrix over ZZ and make it over QQ instead
matZZtoQQ = (M) -> (
E:=entries(M);    
return matrix apply(#E, i -> apply(#(E_i), j -> (1/1)*E_i_j))    
);


--Compute the sample covariance matrix
--See Bernd's notes, p. 44
--Assume data is entered as a list of matrices (row vectors)
sampleCovarianceMatrix = (L) -> (
n:=#L;
L=apply(#L, i -> if ring(L_i)===ZZ then matZZtoQQ(L_i) else L_i);
Xbar:=(1/n)*(sum L);
return (1/n)*(sum apply(n, i -> (transpose (L_i-Xbar))*(L_i-Xbar)) );    
);

V=sampleCovarianceMatrix(U);


JacobianMatrixOfRationalFunction = (F) -> (
f:=numerator(F);
g:=denominator(F);
R:=ring(f);
answer:=diff(vars(R), f) * g - diff(vars(R), g)*f;
answer=substitute(answer, FR);
matrix({{(1/g)^2}})*answer
);


C1 = trace(Sinv * V)
C1derivative=JacobianMatrixOfRationalFunction(C1);

{*
Use the log likelihood from Bernd's Prop. 2.1.12, p. 45
*}

LL = (substitute((jacobian(matrix{{det(S)}})), FR))*matrix{{(-(#U)/(2*det(S)))}} - (transpose C1derivative);
print LL
I = ideal(0);
for i to (numgens target LL)-1 do (I = I + ideal(numerator(LL_(i, 0))));
use R
ELL=flatten entries(LL);
J=ideal apply(#ELL, i -> lift(numerator(ELL_i),R));
J = saturate(J, det(S));





Sigma = matrix{{s_(1,1), s_(1,2), s_(1,3)}, {s_(1,2), s_(2,2), s_(2,3)}, {s_(1,3), s_(2,3), s_(3,3)}}
S
J = ideal(Sigma - S)
I = I + J
I
dim I
degree I
I == radical I
mingens I
