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
