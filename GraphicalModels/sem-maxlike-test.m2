
restart
--We input a graphical model G, and output its likelihood equations for data U. 
needsPackage "GraphicalModels";
G = mixedGraph(digraph {{1,{2,3}},{2,{3}},{3,{4}}},bigraph {{1,2},{2,4}});
--R = gaussianRing (G,sVariableName=>t)
R = gaussianRing (G)

L=directedEdgesMatrix R
P=bidirectedEdgesMatrix R
S=covarianceMatrix R 
K=1-L

paramS=(transpose inverse K)*P*inverse K

allI=ideal matrix for i from 0 to 4-1 list for j from 0 to 4-1 list ((S_(i,j)-paramS_(i,j)))
idealS=eliminate({l_(1,2), l_(1,3), l_(2,3), l_(3,4), p_(1,1), p_(2,2), p_(3,3), p_(4,4),p_(1,2), p_(2,4)},allI)

R2=QQ[ s_(1,1), s_(1,2), s_(1,3), s_(1,4), s_(2,2), s_(2,3),
      s_(2,4), s_(3,3), s_(3,4), s_(4,4),k_(1,1), k_(1,2), k_(1,3), k_(1,4),
 k_(2,2), k_(2,3), k_(2,4), k_(3,3), k_(3,4), k_(4,4)]

genK=genericSymmetricMatrix(R2,k_(1,1),4)---inverse of S

Ib=ideal(    sub(S,R2)*genK-1)+sub(idealS,R2)
Iwin=eliminate({s_(1,1), s_(1,2), s_(1,3), s_(1,4), s_(2,2), s_(2,3),
      s_(2,4), s_(3,3), s_(3,4), s_(4,4)},Ib)
detGenK=det genK
myF=(sub(idealS,R2))_0

U=matrix {{2,3,5,7},
    {3,11,13,17},
    {5,13,19,23},
    {7,17,23,29_R2}}


--matrix for i from 1 to 4 list for j from 1 to 4 list (diff( k_(i,j),detGenK )/detGenK-diff(k_(i,j),trace(U*genK)))
mA=matrix{flatten  for i from 1 to 4 list for j from i to 4 list (diff( k_(i,j),detGenK )-detGenK*diff(k_(i,j),trace(U*genK)))}
mB=matrix{flatten  for i from 1 to 4 list for j from i to 4 list diff(k_(i,j),Iwin_0)}
J=flatten mA||flatten mB
likeEq=minors(2,J)+Iwin
--maybe we need to saturate singular locus. 
degree likeEq
