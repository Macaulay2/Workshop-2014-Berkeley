-- Add here various tests to check whether procedures aren't broken



restart
loadPackage "PosChar"
R=ZZ/2[x,y,z,t]

u=random(4,R);
I=ideal(random(2,R),random(4,R),random(4,R));
ethRoot (I, 1)
minimalCompatible(I,u,1)
minimalCompatible(I,u,1,1)


A=matrix{{x^2,y^2},{x*y*z,1}};
ethRoot (A, 1)
U=A*A;
minimalCompatible(A,U,1)




--------------------------------------------------

restart
loadPackage "PosChar"

p=2;
R0=ZZ/p[u,v,x_1,x_2,x_3,x_4, MonomialOrder=>{2,4}];
I0=ideal(x_1-u^4, x_2-u^3*v, x_3-u*v^3, x_4-v^4);
G=gens gb I0;
I=selectInSubring(1,G);

R=ZZ/p[x_1,x_2,x_3,x_4];
I=substitute(I,R);
generatingMorphisms=findGeneratingMorphisms (I) 


A=relations ((generatingMorphisms#2)#0)
U=(generatingMorphisms#2)#1
findHSLloci(A,U)


nonFInjectiveLocus (I)

