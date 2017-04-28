newPackage(
        "UnitTestsNCA",
        Version => "0.1", 
        Date => "",
        Authors => {{Name => "", 
                  Email => "", 
                  HomePage => ""}},
        Headline => "Short tests for basic correctness of NCAlgebra.m2",
        DebuggingMode => true
        )

--- check basic ring operations code
TEST ///
needsPackage "NCAlgebra"
A = QQ{x,y,z}
f = y*z + z*y - x^2
g = x*z + z*x - y^2
h = z^2 - x*y - y*x
assert(f*g == z*y*z*x-z*y^3+z*y*x*z+y*z^2*x-y*z*y^2+y*z*x*z-x^2*z*x+x^2*y^2-x^3*z)
assert(f^2 == z*y*z*y+z*y^2*z-z*y*x^2+y*z^2*y+y*z*y*z-y*z*x^2-x^2*z*y-x^2*y*z+x^4)
assert(f-g == z*y-z*x+y*z+y^2-x*z-x^2)
assert(3*g == 3*z*x-3*y^2+3*x*z)
assert(f+g == z*y+z*x+y*z-y^2+x*z-x^2)
B = A/ncIdeal{f,g,h}
j = -y^3-x*y*z+y*x*z+x^3
k = x^2 + y^2 + z^2
j*k
k^3
assert(k^3 == y^6+3*y*x*y^4+3*y*x*y*x*y^2+y*x*y*x*y*x+3*x*y^5+3*x*y*x*y^3+x*y*x*y*x*y+9*x^2*y^4+9*x^2*y*x*y^2+3*x^2*y*x*y*x+9*x^3*y^3+3*x^3*y*x*y+9*x^4*y^2+3*x^4*y*x+3*x^5*y+x^6)
assert(j*k == -y^5+y*x*y^2*z-y*x*y^3+y*x*y*x*z-x*y^3*z-x*y^4-x*y*x*y*z-x^2*y^3+x^2*y*x*z-x^3*y*z+x^3*y^2+x^3*y*x+x^4*y+x^5)
assert(j-k == -y^3-y^2+y*x*z-y*x-x*y*z-x*y+x^3-x^2)
assert(j+k == -y^3+y^2+y*x*z+y*x-x*y*z+x*y+x^3+x^2)
assert(3*k == 3*y^2+3*y*x+3*x*y+3*x^2)
///

--- checking central elements (skylanin)
TEST ///
needsPackage "NCAlgebra"
A = QQ{x,y,z}
I = ncIdeal { y*z + z*y - x^2,x*z + z*x - y^2,z^2 - x*y - y*x}
B = A/I
g = -y^3-x*y*z+y*x*z+x^3
h = x^2 + y^2 + z^2
assert (isCentral h)
assert (isCentral g)
assert(centralElements(B,2) == ncMatrix {{x^2, y*x+x*y, y^2}})
assert(centralElements(B,3) == ncMatrix {{y^3-y*x*z+x*y*z-x^3}})
///

--- checking normal form code with non-field base
TEST ///
needsPackage "NCAlgebra"
R = QQ[a,b,c,d]/ideal{a*b+c*d}
A = R {x,y,z}
I = ncIdeal {a*x*y,b*z^2}
Igb = ncGroebnerBasis(I, InstallGB=>true)
assert(c*z^2 % Igb == c*z^2)
assert(b*z^2 % Igb == 0)
///

--- checking normal form bergman code
TEST ///
needsPackage "NCAlgebra"
A = QQ{x,y,z}
f = y*z + z*y - x^2
g = x*z + z*x - y^2
h = z^2 - x*y - y*x
I = ncIdeal {f,g,h}
Igb = ncGroebnerBasis I
red = y*x*y*x*y*x*y*x*y*x*y*x*y*x*y*x*z+x*y*x*y*x*y*x*y*x*y*x*y*x*y*x*y*z+8*x^2*y*x*y*x*y*x*y*x*y*x*y*x*y^2*z+8*x^3*y*x*y*x*y*x*y*x*y*x*y^3*z+28*x^4*y*x*y*x*y*x*y*x*y^4*z+28*x^5*y*x*y*x*y*x*y^5*z+56*x^6*y*x*y*x*y^6*z+56*x^7*y*x*y^7*z+70*x^8*y^8*z
assert(normalFormBergman(z^17,Igb) == red)
///

--- checking basis of algebra (sklyanin)
TEST ///
needsPackage "NCAlgebra"
A = QQ{x,y,z}
f = y*z + z*y - x^2
g = x*z + z*x - y^2
h = z^2 - x*y - y*x
I = ncIdeal {f,g,h}
B = A/I
basis1 = ncMatrix {{x^4, x^2*y*x, y*x*y*x, x^3*y, x*y*x*y, x^2*y^2, y*x*y^2, x*y^3, y^4, x^3*z, x*y*x*z, x^2*y*z, y*x*y*z, x*y^2*z, y^3*z}}
basis2 = basis(4,B)
assert(basis1 == basis2)
///

--- checking basis of algebra (quantum polynomial ring)
TEST ///
needsPackage "NCAlgebra"
B = skewPolynomialRing(QQ,(-1)_QQ,{x,y,z,w})
basis1 = ncMatrix {{x^3, x^2*y, x*y^2, y^3, x^2*z, x*y*z, y^2*z, x*z^2, y*z^2, z^3, x^2*w, x*y*w, y^2*w, x*z*w, y*z*w, z^2*w, x*w^2, y*w^2, z*w^2, w^3}}
basis2 = basis(3,B)
assert(basis1 == basis2)
///

--- checking ore extensions
TEST ///
needsPackage "NCAlgebra"
B = skewPolynomialRing(QQ,(-1)_QQ,{x,y,z,w})
sigma = ncMap(B,B,{y,z,w,x})
C = oreExtension(B,sigma,a)
sigmaC = ncMap(C,C,{y,z,w,x,a})
assert(normalElements(sigmaC,1) == ncMatrix {{a}})
assert(normalElements(sigmaC,2) == 0)
assert(normalElements(sigmaC @@ sigmaC,2) == ncMatrix {{a^2}})
assert(matrix normalAutomorphism a == ncMatrix {{y,z,w,x,a}})
assert(matrix normalAutomorphism a^2 == ncMatrix {{z,w,x,y,a}})
///

--- checking endomorphism ring code
TEST ///
debug needsPackage "NCAlgebra"
Q = QQ[a,b,c]
R = Q/ideal{a*b-c^2}
kRes = res(coker vars R, LengthLimit=>7);
M = coker kRes.dd_5
B = endomorphismRing(M,X);
gensI = gens ideal B;
gensIMin = minimizeRelations(gensI, Verbosity=>1);
checkHomRelations(gensIMin,B.cache.endomorphismRingGens)
///

--- Another example testing out endomorphism code
TEST ///
debug needsPackage "NCAlgebra"
Q = QQ[a,b,c,d]
R = Q/ideal{a*b+c*d}
kRes = res(coker vars R, LengthLimit=>7);
M = coker kRes.dd_5
time B = endomorphismRing(M,X);
gensI = gens ideal B;
time gensIMin = minimizeRelations(gensI, Verbosity=>1);
checkHomRelations(gensIMin,B.cache.endomorphismRingGens);
///

--- Test creation of matrices, as matrices of matrices
TEST ///
needsPackage "NCAlgebra"
A = QQ{a,b,c,d}
M = ncMatrix {{a,b,c,d}}
N = ncMatrix {{M,2*M,3*M},{4*M,5*M,6*M}}
assert (N == ncMatrix {{a, b, c, d, 2*a, 2*b, 2*c, 2*d, 3*a, 3*b, 3*c, 3*d}, {4*a, 4*b, 4*c, 4*d, 5*a, 5*b, 5*c, 5*d, 6*a, 6*b, 6*c, 6*d}})
///

--- Test NCMatrix commands
TEST ///
needsPackage "NCAlgebra"
A = QQ{x,y,z}
f = y*z + z*y - x^2
g = x*z + z*x - y^2
h = z^2 - x*y - y*x
I = ncIdeal {f,g,h}
Igb = ncGroebnerBasis I
M = ncMatrix {{x, y, z}}
sigma = ncMap(A,A,{y,z,x})
N = ncMatrix {{M},{sigma M}, {sigma sigma M}}
Nred = N^3 % Igb
B = A/I
phi = ncMap(B,A,gens B)
NB = phi N
N3B = NB^3
X = NB + 3*NB
Y = NB | 2*NB
Z = X || NB
Xans = ncMatrix {{4*x, 4*y, 4*z}, {4*y, 4*z, 4*x}, {4*z, 4*x, 4*y}}
Yans = ncMatrix {{x, y, z, 2*x, 2*y, 2*z}, {y, z, x, 2*y, 2*z, 2*x}, {z, x, y, 2*z, 2*x, 2*y}}
Zans = ncMatrix {{4*x, 4*y, 4*z}, {4*y, 4*z, 4*x}, {4*z, 4*x, 4*y}, {x, y, z}, {y, z, x}, {z, x, y}}
answer = ncMatrix {{-y^2*z+y^3+y*x*z-y*x*y+x*y*z+x*y^2+2*x*y*x+x^2*z+3*x^2*y, y^2*z+y*x*z+2*y*x*y+x*y*z+3*x*y^2-x*y*x-x^2*z+x^2*y+x^3, 2*y^2*z+y^3+y*x*y+x*y*x+2*x^2*z+x^3}, {y^2*z+y*x*z+2*y*x*y+x*y*z+3*x*y^2-x*y*x-x^2*z+x^2*y+x^3, 2*y^2*z+y^3+y*x*y+x*y*x+2*x^2*z+x^3, -y^2*z+y^3+y*x*z-y*x*y+x*y*z+x*y^2+2*x*y*x+x^2*z+3*x^2*y}, {2*y^2*z+y^3+y*x*y+x*y*x+2*x^2*z+x^3, -y^2*z+y^3+y*x*z-y*x*y+x*y*z+x*y^2+2*x*y*x+x^2*z+3*x^2*y, y^2*z+y*x*z+2*y*x*y+x*y*z+3*x*y^2-x*y*x-x^2*z+x^2*y+x^3}}
assert(X == Xans)
assert(Y == Yans)
assert(Z == Zans)
assert(N3B == phi Nred)
assert(N3B == answer)
///

--- Check rightKernelBergman
TEST ///
needsPackage "NCAlgebra"
A = QQ{x,y,z}
f1 = y*z + z*y - x^2
f2 = x*z + z*x - y^2
f3 = z^2 - x*y - y*x
g = -y^3-x*y*z+y*x*z+x^3
I = ncIdeal {f1,f2,f3,g}
B = A/I
M3 = ncMatrix {{x,y,z,0},
               {-y*z-2*x^2,-y*x,z*x-x*z,x},
               {x*y-2*y*x,x*z,-x^2,y},
               {-y^2-z*x,x^2,-x*y,z}}
assignDegrees(M3,{1,0,0,0},{2,2,2,1})
assert isHomogeneous M3
ker1M3 = rightKernelBergman(M3)
assert isHomogeneous ker1M3
assert(M3*ker1M3 == 0)
ker2M3 = rightKernelBergman(ker1M3)
assert isHomogeneous ker2M3
assert(ker1M3*ker2M3 == 0)
ker3M3 = rightKernelBergman(ker2M3)
assert isHomogeneous ker3M3  
assert(ker2M3*ker3M3 == 0)
///

--- Check mingens of a NCMatrix (code not yet finished)
TEST ///
///

--- Check NCIdeal operations
--- can't really do much with these yet, except define quotient rings
TEST ///
needsPackage "NCAlgebra"
A = QQ{x,y,z}
f1 = y*z + z*y - x^2
f2 = x*z + z*x - y^2
f3 = z^2 - x*y - y*x
g = -y^3-x*y*z+y*x*z+x^3
I = ncIdeal {f1,f2,f3}
J = ncIdeal g
I + J
IJgb = ncGroebnerBasis (I+J)
-- need to make some tests here
///

{*
--- Check basis of ideal code
--- was this code ever added?
TEST ///
///
*}

--- Check NCRightIdeal operations
--- can't really do much with these yet.
TEST ///
needsPackage "NCAlgebra"
A = QQ{x,y,z}
f1 = y*z + z*y - x^2
f2 = x*z + z*x - y^2
f3 = z^2 - x*y - y*x
g = -y^3-x*y*z+y*x*z+x^3
I = ncRightIdeal {f1,f2,f3}
J = ncRightIdeal g
I + J
///

--- Check NCLeftIdeal operations
--- can't really do much with these yet.
TEST ///
needsPackage "NCAlgebra"
A = QQ{x,y,z}
f1 = y*z + z*y - x^2
f2 = x*z + z*x - y^2
f3 = z^2 - x*y - y*x
g = -y^3-x*y*z+y*x*z+x^3
I = ncLeftIdeal {f1,f2,f3}
J = ncLeftIdeal g
I + J
///

{*
--- Test gbFromOutputFile
TEST ///
-- this test works as of 7/19/2012
-- however, it does not work in 'check UnitTestsNCA', because
-- I don't know how to get the checker to see the gb file.  How can I do this?
restart
needsPackage "NCAlgebra"
A=QQ{a, b, c, d, e, f, g, h}
-- this got very slow all of a sudden... speed this up?  Not really sure how...
I = time gbFromOutputFile(A,"UghABCgb6.txt", ReturnIdeal=>true);
B=A/I
IdegListAns = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}
assert(IdegListAns == (I.generators / degree))
F = a^7+b^7+c^7+d^7+e^7+f^7+g^7+h^7;
--- sometimes bergman calls are great!
bas=basis(2,B);
X = flatten entries (F*bas);
-- change back to avoid calls to bergman engine for 'value' and equality test.
XA = apply(X, x -> promote(x,A));
use A
firstFewAns = {e^7*a*b+e^7*a^2+7*e^6*a*b*h-7*e^6*a*b*e+3*e^5*b*e^2*b-3*e^5*a*e^2*b+21*e^5*a*b*h^2-42*e^5*a*b*e*h+21*e^5*a*b*e^2+9*e^4*b*e^3*b-10*e^4*b*e^3*a+28*e^4*b*e^2*b*h-57*e^4*b*e^2*b*e+e^4*a*e^3*b-28*e^4*a*e^2*b*h+57*e^4*a*e^2*b*e+35*e^4*a*b*h^3-105*e^4*a*b*e*h^2+105*e^4*a*b*e^2*h-35*e^4*a*b*e^3-20*e^3*b*e^4*b+19*e^3*b*e^4*a+28*e^3*b*e^3*b*h+44*e^3*b*e^3*b*e-49*e^3*b*e^3*a*h+34*e^3*b*e^3*a*e+105*e^3*b*e^2*b*h^2-357*e^3*b*e^2*b*e*h+231*e^3*b*e^2*b*e^2+e^3*a*e^4*b+21*e^3*a*e^3*b*h-78*e^3*a*e^3*b*e-105*e^3*a*e^2*b*h^2+357*e^3*a*e^2*b*e*h-231*e^3*a*e^2*b*e^2+35*e^3*a*b*h^4-140*e^3*a*b*e*h^3+210*e^3*a*b*e^2*h^2-140*e^3*a*b*e^3*h+35*e^3*a*b*e^4+15*e^2*b*e^5*b-15*e^2*b*e^5*a-70*e^2*b*e^4*b*h+70*e^2*b*e^4*a*h-27*e^2*b*e^4*a*e+280*e^2*b*e^3*b*e*h-198*e^2*b*e^3*b*e^2-105*e^2*b*e^3*a*h^2+224*e^2*b*e^3*a*e*h-189*e^2*b*e^3*a*e^2+210*e^2*b*e^2*b*h^3-945*e^2*b*e^2*b*e*h^2+1281*e^2*b*e^2*b*e^2*h-615*e^2*b*e^2*b*e^3+27*e^2*a*e^4*b*e+105*e^2*a*e^3*b*h^2-504*e^2*a*e^3*b*e*h+387*e^2*a*e^3*b*e^2-210*e^2*a*e^2*b*h^3+945*e^2*a*e^2*b*e*h^2-1281*e^2*a*e^2*b*e^2*h+615*e^2*a*e^2*b*e^3+21*e^2*a*b*h^5-105*e^2*a*b*e*h^4+210*e^2*a*b*e^2*h^3-210*e^2*a*b*e^3*h^2+105*e^2*a*b*e^4*h-21*e^2*a*b*e^5-7*e*b*e^6*b+6*e*b*e^6*a+28*e*b*e^5*b*h+6*e*b*e^5*b*e-35*e*b*e^5*a*h-105*e*b*e^4*b*h^2+35*e*b*e^4*b*e*h-15*e*b*e^4*b*e^2+105*e*b*e^4*a*h^2-182*e*b*e^4*a*e*h+171*e*b*e^4*a*e^2-140*e*b*e^3*b*h^3+840*e*b*e^3*b*e*h^2-1204*e*b*e^3*b*e^2*h+628*e*b*e^3*b*e^3-140*e*b*e^3*a*h^3+630*e*b*e^3*a*e*h^2-1085*e*b*e^3*a*e^2*h+600*e*b*e^3*a*e^3+210*e*b*e^2*b*h^4-1260*e*b*e^2*b*e*h^3+2835*e*b*e^2*b*e^2*h^2-2975*e*b*e^2*b*e^3*h+1194*e*b*e^2*b*e^4-105*e*b*e^2*a*h^4+840*e*b*e^2*a*e*h^3-2520*e*b*e^2*a*e^2*h^2+3360*e*b*e^2*a*e^3*h-1680*e*b*e^2*a*e^4+42*e*b*e*b*h^5-420*e*b*e*b*e*h^4+1680*e*b*e*b*e^2*h^3-3360*e*b*e*b*e^3*h^2+3360*e*b*e*b*e^4*h-1344*e*b*e*b*e^5+e*a*e^6*b+7*e*a*e^5*b*h-6*e*a*e^5*b*e+147*e*a*e^4*b*e*h-156*e*a*e^4*b*e^2+280*e*a*e^3*b*h^3-1470*e*a*e^3*b*e*h^2+2289*e*a*e^3*b*e^2*h-1228*e*a*e^3*b*e^3-105*e*a*e^2*b*h^4+420*e*a*e^2*b*e*h^3-315*e*a*e^2*b*e^2*h^2-385*e*a*e^2*b*e^3*h+486*e*a*e^2*b*e^4-42*e*a*e*b*h^5+420*e*a*e*b*e*h^4-1680*e*a*e*b*e^2*h^3+3360*e*a*e*b*e^3*h^2-3360*e*a*e*b*e^4*h+1344*e*a*e*b*e^5+7*e*a*b*h^6-42*e*a*b*e*h^5+105*e*a*b*e^2*h^4-140*e*a*b*e^3*h^3+105*e*a*b*e^4*h^2-42*e*a*b*e^5*h+7*e*a*b*e^6+d^7*a^2+c^7*a^2-b*e^7*a-14*b*e^6*b*h+7*b*e^6*b*e+7*b*e^6*a*h+42*b*e^5*b*e*h-21*b*e^5*b*e^2-21*b*e^5*a*h^2-70*b*e^4*b*h^3+105*b*e^4*b*e*h^2-105*b*e^4*b*e^2*h+35*b*e^4*b*e^3+105*b*e^4*a*h^3-525*b*e^4*a*e*h^2+1008*b*e^4*a*e^2*h-623*b*e^4*a*e^3-210*b*e^3*b*h^4+1260*b*e^3*b*e*h^3-2940*b*e^3*b*e^2*h^2+3304*b*e^3*b*e^3*h-1449*b*e^3*b*e^4+70*b*e^3*a*h^4-700*b*e^3*a*e*h^3+2415*b*e^3*a*e^2*h^2-3598*b*e^3*a*e^3*h+1953*b*e^3*a*e^4-42*b*e^2*b*h^5+420*b*e^2*b*e*h^4-1890*b*e^2*b*e^2*h^3+4305*b*e^2*b*e^3*h^2-4851*b*e^2*b*e^4*h+2142*b*e^2*b*e^5-21*b*e^2*a*h^5+105*b*e^2*a*e*h^4-840*b*e^2*a*e^3*h^2+1680*b*e^2*a*e^4*h-1008*b*e^2*a*e^5-42*b*e*b*e*h^5+210*b*e*b*e^2*h^4-280*b*e*b*e^3*h^3-210*b*e*b*e^4*h^2+756*b*e*b*e^5*h-462*b*e*b*e^6-7*b*e*a*h^6+84*b*e*a*e*h^5-420*b*e*a*e^2*h^4+1120*b*e*a*e^3*h^3-1680*b*e*a*e^4*h^2+1344*b*e*a*e^5*h-448*b*e*a*e^6+42*b^2*e^2*h^5-210*b^2*e^3*h^4+490*b^2*e^4*h^3-630*b^2*e^5*h^2+434*b^2*e^6*h-126*b^2*e^7+a*e^7*b+7*a*e^6*b*h-7*a*e^6*b*e+21*a*e^5*b*h^2-42*a*e^5*b*e*h+21*a*e^5*b*e^2-35*a*e^4*b*h^3+420*a*e^4*b*e*h^2-903*a*e^4*b*e^2*h+588*a*e^4*b*e^3+140*a*e^3*b*h^4-560*a*e^3*b*e*h^3+525*a*e^3*b*e^2*h^2+294*a*e^3*b*e^3*h-504*a*e^3*b*e^4+63*a*e^2*b*h^5-525*a*e^2*b*e*h^4+1890*a*e^2*b*e^2*h^3-3465*a*e^2*b*e^3*h^2+3171*a*e^2*b*e^4*h-1134*a*e^2*b*e^5+7*a*e*b*h^6-42*a*e*b*e*h^5+210*a*e*b*e^2*h^4-840*a*e*b*e^3*h^3+1890*a*e*b*e^4*h^2-2100*a*e*b*e^5*h+910*a*e*b*e^6-7*a*b*e*h^6-21*a*b*e^2*h^5+175*a*b*e^3*h^4-455*a*b*e^4*h^3+609*a*b*e^5*h^2-427*a*b*e^6*h+125*a*b*e^7+a^2*h^7+a^2*g^7+a^2*f^7+a^2*b^7+a^9, 2*e^7*a*b+7*e^6*a*b*h-7*e^6*a*b*e+6*e^5*b*e^2*b-6*e^5*a*e^2*b+21*e^5*a*b*h^2-42*e^5*a*b*e*h+21*e^5*a*b*e^2-4*e^4*b*e^3*b-4*e^4*b*e^3*a+42*e^4*b*e^2*b*h-54*e^4*b*e^2*b*e+8*e^4*a*e^3*b-42*e^4*a*e^2*b*h+54*e^4*a*e^2*b*e+35*e^4*a*b*h^3-105*e^4*a*b*e*h^2+105*e^4*a*b*e^2*h-35*e^4*a*b*e^3+3*e^3*b*e^4*a-28*e^3*b*e^3*b*h+40*e^3*b*e^3*b*e-28*e^3*b*e^3*a*h+40*e^3*b*e^3*a*e+126*e^3*b*e^2*b*h^2-336*e^3*b*e^2*b*e*h+228*e^3*b*e^2*b*e^2-3*e^3*a*e^4*b+56*e^3*a*e^3*b*h-80*e^3*a*e^3*b*e-126*e^3*a*e^2*b*h^2+336*e^3*a*e^2*b*e*h-228*e^3*a*e^2*b*e^2+35*e^3*a*b*h^4-140*e^3*a*b*e*h^3+210*e^3*a*b*e^2*h^2-140*e^3*a*b*e^3*h+35*e^3*a*b*e^4+21*e^2*b*e^4*a*h-33*e^2*b*e^4*a*e-84*e^2*b*e^3*b*h^2+252*e^2*b*e^3*b*e*h-192*e^2*b*e^3*b*e^2-84*e^2*b*e^3*a*h^2+252*e^2*b*e^3*a*e*h-192*e^2*b*e^3*a*e^2+210*e^2*b*e^2*b*h^3-882*e^2*b*e^2*b*e*h^2+1260*e^2*b*e^2*b*e^2*h-612*e^2*b*e^2*b*e^3-21*e^2*a*e^4*b*h+33*e^2*a*e^4*b*e+168*e^2*a*e^3*b*h^2-504*e^2*a*e^3*b*e*h+384*e^2*a*e^3*b*e^2-210*e^2*a*e^2*b*h^3+882*e^2*a*e^2*b*e*h^2-1260*e^2*a*e^2*b*e^2*h+612*e^2*a*e^2*b*e^3+21*e^2*a*b*h^5-105*e^2*a*b*e*h^4+210*e^2*a*b*e^2*h^3-210*e^2*a*b*e^3*h^2+105*e^2*a*b*e^4*h-21*e^2*a*b*e^5+63*e*b*e^4*a*h^2-210*e*b*e^4*a*e*h+177*e*b*e^4*a*e^2-140*e*b*e^3*b*h^3+672*e*b*e^3*b*e*h^2-1092*e*b*e^3*b*e^2*h+600*e*b*e^3*b*e^3-140*e*b*e^3*a*h^3+672*e*b*e^3*a*e*h^2-1092*e*b*e^3*a*e^2*h+600*e*b*e^3*a*e^3+210*e*b*e^2*b*h^4-1260*e*b*e^2*b*e*h^3+2898*e*b*e^2*b*e^2*h^2-3024*e*b*e^2*b*e^3*h+1206*e*b*e^2*b*e^4-105*e*b*e^2*a*h^4+840*e*b*e^2*a*e*h^3-2520*e*b*e^2*a*e^2*h^2+3360*e*b*e^2*a*e^3*h-1680*e*b*e^2*a*e^4+42*e*b*e*b*h^5-420*e*b*e*b*e*h^4+1680*e*b*e*b*e^2*h^3-3360*e*b*e*b*e^3*h^2+3360*e*b*e*b*e^4*h-1344*e*b*e*b*e^5-63*e*a*e^4*b*h^2+210*e*a*e^4*b*e*h-177*e*a*e^4*b*e^2+280*e*a*e^3*b*h^3-1344*e*a*e^3*b*e*h^2+2184*e*a*e^3*b*e^2*h-1200*e*a*e^3*b*e^3-105*e*a*e^2*b*h^4+420*e*a*e^2*b*e*h^3-378*e*a*e^2*b*e^2*h^2-336*e*a*e^2*b*e^3*h+474*e*a*e^2*b*e^4-42*e*a*e*b*h^5+420*e*a*e*b*e*h^4-1680*e*a*e*b*e^2*h^3+3360*e*a*e*b*e^3*h^2-3360*e*a*e*b*e^4*h+1344*e*a*e*b*e^5+7*e*a*b*h^6-42*e*a*b*e*h^5+105*e*a*b*e^2*h^4-140*e*a*b*e^3*h^3+105*e*a*b*e^4*h^2-42*e*a*b*e^5*h+7*e*a*b*e^6+6*d^5*a*d^2*b+4*d^4*b*d^3*a-8*d^4*a*d^3*b-12*d^4*a*d^2*b*d-3*d^3*b*d^4*a-12*d^3*b*d^3*a*d+3*d^3*a*d^4*b+24*d^3*a*d^3*b*d+18*d^3*a*d^2*b*d^2+12*d^2*b*d^4*a*d+24*d^2*b*d^3*a*d^2-12*d^2*a*d^4*b*d-48*d^2*a*d^3*b*d^2-24*d^2*a*d^2*b*d^3-30*d*b*d^4*a*d^2-40*d*b*d^3*a*d^3+105*d*b*d^2*a*d^4+30*d*a*d^4*b*d^2+80*d*a*d^3*b*d^3-75*d*a*d^2*b*d^4-42*d*a*d*b*d^5+c^7*a^2+6*c^5*a*c^2*b-6*c^5*a*c^2*a+4*c^4*b*c^3*a-8*c^4*a*c^3*b+4*c^4*a*c^3*a-12*c^4*a*c^2*b*c+12*c^4*a*c^2*a*c-3*c^3*b*c^4*a-12*c^3*b*c^3*a*c+3*c^3*a*c^4*b+24*c^3*a*c^3*b*c-12*c^3*a*c^3*a*c+18*c^3*a*c^2*b*c^2-18*c^3*a*c^2*a*c^2+12*c^2*b*c^4*a*c+24*c^2*b*c^3*a*c^2-12*c^2*a*c^4*b*c-48*c^2*a*c^3*b*c^2+24*c^2*a*c^3*a*c^2-24*c^2*a*c^2*b*c^3+24*c^2*a*c^2*a*c^3-30*c*b*c^4*a*c^2-40*c*b*c^3*a*c^3+105*c*b*c^2*a*c^4+30*c*a*c^4*b*c^2+80*c*a*c^3*b*c^3-40*c*a*c^3*a*c^3-75*c*a*c^2*b*c^4-30*c*a*c^2*a*c^4-42*c*a*c*b*c^5+42*c*a*c*a*c^5-b*e^7*b-7*b*e^6*b*h+7*b*e^6*b*e-21*b*e^5*b*h^2+42*b*e^5*b*e*h-21*b*e^5*b*e^2-35*b*e^4*b*h^3+105*b*e^4*b*e*h^2-105*b*e^4*b*e^2*h+35*b*e^4*b*e^3+105*b*e^4*a*h^3-567*b*e^4*a*e*h^2+1029*b*e^4*a*e^2*h-627*b*e^4*a*e^3-175*b*e^3*b*h^4+1120*b*e^3*b*e*h^3-2814*b*e^3*b*e^2*h^2+3248*b*e^3*b*e^3*h-1439*b*e^3*b*e^4+70*b*e^3*a*h^4-700*b*e^3*a*e*h^3+2436*b*e^3*a*e^2*h^2-3612*b*e^3*a*e^3*h+1956*b*e^3*a*e^4-21*b*e^2*b*h^5+315*b*e^2*b*e*h^4-1680*b*e^2*b*e^2*h^3+4116*b*e^2*b*e^3*h^2-4767*b*e^2*b*e^4*h+2127*b*e^2*b*e^5-21*b*e^2*a*h^5+105*b*e^2*a*e*h^4-840*b*e^2*a*e^3*h^2+1680*b*e^2*a*e^4*h-1008*b*e^2*a*e^5+7*b*e*b*h^6-84*b*e*b*e*h^5+315*b*e*b*e^2*h^4-420*b*e*b*e^3*h^3-105*b*e*b*e^4*h^2+714*b*e*b*e^5*h-455*b*e*b*e^6-7*b*e*a*h^6+84*b*e*a*e*h^5-420*b*e*a*e^2*h^4+1120*b*e*a*e^3*h^3-1680*b*e*a*e^4*h^2+1344*b*e*a*e^5*h-448*b*e*a*e^6+60*b*d^4*a*d^3-150*b*d^3*a*d^4+84*b*d^2*a*d^5+7*b*d*a*d^6+60*b*c^4*a*c^3-150*b*c^3*a*c^4+84*b*c^2*a*c^5+7*b*c*a*c^6-7*b^2*e*h^6+63*b^2*e^2*h^5-245*b^2*e^3*h^4+525*b^2*e^4*h^3-651*b^2*e^5*h^2+441*b^2*e^6*h-127*b^2*e^7+a*e^7*b+7*a*e^6*b*h-7*a*e^6*b*e+21*a*e^5*b*h^2-42*a*e^5*b*e*h+21*a*e^5*b*e^2-70*a*e^4*b*h^3+462*a*e^4*b*e*h^2-924*a*e^4*b*e^2*h+592*a*e^4*b*e^3+105*a*e^3*b*h^4-420*a*e^3*b*e*h^3+378*a*e^3*b*e^2*h^2+364*a*e^3*b*e^3*h-517*a*e^3*b*e^4+42*a*e^2*b*h^5-420*a*e^2*b*e*h^4+1680*a*e^2*b*e^2*h^3-3276*a*e^2*b*e^3*h^2+3087*a*e^2*b*e^4*h-1119*a*e^2*b*e^5+105*a*e*b*e^2*h^4-700*a*e*b*e^3*h^3+1785*a*e*b*e^4*h^2-2058*a*e*b*e^5*h+903*a*e*b*e^6+a*d^7*b+7*a*d^6*b*g-7*a*d^6*b*d+21*a*d^5*b*g^2-42*a*d^5*b*d*g+21*a*d^5*b*d^2+35*a*d^4*b*g^3-105*a*d^4*b*d*g^2+105*a*d^4*b*d^2*g-95*a*d^4*b*d^3+35*a*d^3*b*g^4-140*a*d^3*b*d*g^3+210*a*d^3*b*d^2*g^2-140*a*d^3*b*d^3*g+125*a*d^3*b*d^4+21*a*d^2*b*g^5-105*a*d^2*b*d*g^4+210*a*d^2*b*d^2*g^3-210*a*d^2*b*d^3*g^2+105*a*d^2*b*d^4*g-15*a*d^2*b*d^5+7*a*d*b*g^6-42*a*d*b*d*g^5+105*a*d*b*d^2*g^4-140*a*d*b*d^3*g^3+105*a*d*b*d^4*g^2-42*a*d*b*d^5*g-28*a*d*b*d^6+a*c^7*b-a*c^7*a+7*a*c^6*b*f-7*a*c^6*b*c-7*a*c^6*a*f+7*a*c^6*a*c+21*a*c^5*b*f^2-42*a*c^5*b*c*f+21*a*c^5*b*c^2-21*a*c^5*a*f^2+42*a*c^5*a*c*f-21*a*c^5*a*c^2+35*a*c^4*b*f^3-105*a*c^4*b*c*f^2+105*a*c^4*b*c^2*f-95*a*c^4*b*c^3-35*a*c^4*a*f^3+105*a*c^4*a*c*f^2-105*a*c^4*a*c^2*f+35*a*c^4*a*c^3+35*a*c^3*b*f^4-140*a*c^3*b*c*f^3+210*a*c^3*b*c^2*f^2-140*a*c^3*b*c^3*f+125*a*c^3*b*c^4-35*a*c^3*a*f^4+140*a*c^3*a*c*f^3-210*a*c^3*a*c^2*f^2+140*a*c^3*a*c^3*f+25*a*c^3*a*c^4+21*a*c^2*b*f^5-105*a*c^2*b*c*f^4+210*a*c^2*b*c^2*f^3-210*a*c^2*b*c^3*f^2+105*a*c^2*b*c^4*f-15*a*c^2*b*c^5-21*a*c^2*a*f^5+105*a*c^2*a*c*f^4-210*a*c^2*a*c^2*f^3+210*a*c^2*a*c^3*f^2-105*a*c^2*a*c^4*f-69*a*c^2*a*c^5+7*a*c*b*f^6-42*a*c*b*c*f^5+105*a*c*b*c^2*f^4-140*a*c*b*c^3*f^3+105*a*c*b*c^4*f^2-42*a*c*b*c^5*f-28*a*c*b*c^6-7*a*c*a*f^6+42*a*c*a*c*f^5-105*a*c*a*c^2*f^4+140*a*c*a*c^3*f^3-105*a*c*a*c^4*f^2+42*a*c*a*c^5*f+21*a*c*a*c^6+a*b*h^7+a*b*g^7+a*b*f^7-42*a*b*e^2*h^5+210*a*b*e^3*h^4-490*a*b*e^4*h^3+630*a*b*e^5*h^2-434*a*b*e^6*h+126*a*b*e^7-7*a*b*d*g^6+21*a*b*d^2*g^5-35*a*b*d^3*g^4+35*a*b*d^4*g^3-21*a*b*d^5*g^2+7*a*b*d^6*g-2*a*b*d^7-7*a*b*c*f^6+21*a*b*c^2*f^5-35*a*b*c^3*f^4+35*a*b*c^4*f^3-21*a*b*c^5*f^2+7*a*b*c^6*f-2*a*b*c^7+a*b^8+7*a^2*c*f^6-21*a^2*c^2*f^5+35*a^2*c^3*f^4-35*a^2*c^4*f^3+21*a^2*c^5*f^2-7*a^2*c^6*f+2*a^2*c^7+a^8*b, e^7*b*c+e^7*a*c+7*e^6*b*c*h-7*e^6*b*c*e+21*e^5*b*c*h^2-42*e^5*b*c*e*h+21*e^5*b*c*e^2+35*e^4*b*c*h^3-105*e^4*b*c*e*h^2+105*e^4*b*c*e^2*h-35*e^4*b*c*e^3+35*e^3*b*c*h^4-140*e^3*b*c*e*h^3+210*e^3*b*c*e^2*h^2-140*e^3*b*c*e^3*h+35*e^3*b*c*e^4+21*e^2*b*c*h^5-105*e^2*b*c*e*h^4+210*e^2*b*c*e^2*h^3-210*e^2*b*c*e^3*h^2+105*e^2*b*c*e^4*h-21*e^2*b*c*e^5+7*e*b*c*h^6-42*e*b*c*e*h^5+105*e*b*c*e^2*h^4-140*e*b*c*e^3*h^3+105*e*b*c*e^4*h^2-42*e*b*c*e^5*h+7*e*b*c*e^6+d^7*a*c+c^7*a*c-b*e^7*c-7*b*e^6*c*h+7*b*e^6*c*e-21*b*e^5*c*h^2+42*b*e^5*c*e*h-21*b*e^5*c*e^2-35*b*e^4*c*h^3+105*b*e^4*c*e*h^2-105*b*e^4*c*e^2*h+35*b*e^4*c*e^3-35*b*e^3*c*h^4+140*b*e^3*c*e*h^3-210*b*e^3*c*e^2*h^2+140*b*e^3*c*e^3*h-35*b*e^3*c*e^4-21*b*e^2*c*h^5+105*b*e^2*c*e*h^4-210*b*e^2*c*e^2*h^3+210*b*e^2*c*e^3*h^2-105*b*e^2*c*e^4*h+21*b*e^2*c*e^5-7*b*e*c*h^6+42*b*e*c*e*h^5-105*b*e*c*e^2*h^4+140*b*e*c*e^3*h^3-105*b*e*c*e^4*h^2+42*b*e*c*e^5*h-7*b*e*c*e^6+a*e^7*c+7*a*e^6*c*h-7*a*e^6*c*e+21*a*e^5*c*h^2-42*a*e^5*c*e*h+21*a*e^5*c*e^2+35*a*e^4*c*h^3-105*a*e^4*c*e*h^2+105*a*e^4*c*e^2*h-35*a*e^4*c*e^3+35*a*e^3*c*h^4-140*a*e^3*c*e*h^3+210*a*e^3*c*e^2*h^2-140*a*e^3*c*e^3*h+35*a*e^3*c*e^4+21*a*e^2*c*h^5-105*a*e^2*c*e*h^4+210*a*e^2*c*e^2*h^3-210*a*e^2*c*e^3*h^2+105*a*e^2*c*e^4*h-21*a*e^2*c*e^5+7*a*e*c*h^6-42*a*e*c*e*h^5+105*a*e*c*e^2*h^4-140*a*e*c*e^3*h^3+105*a*e*c*e^4*h^2-42*a*e*c*e^5*h+7*a*e*c*e^6+a*d^7*c+7*a*d^6*c*g-7*a*d^6*c*d+21*a*d^5*c*g^2-42*a*d^5*c*d*g+21*a*d^5*c*d^2+35*a*d^4*c*g^3-105*a*d^4*c*d*g^2+105*a*d^4*c*d^2*g-35*a*d^4*c*d^3+35*a*d^3*c*g^4-140*a*d^3*c*d*g^3+210*a*d^3*c*d^2*g^2-140*a*d^3*c*d^3*g+35*a*d^3*c*d^4+21*a*d^2*c*g^5-105*a*d^2*c*d*g^4+210*a*d^2*c*d^2*g^3-210*a*d^2*c*d^3*g^2+105*a*d^2*c*d^4*g-21*a*d^2*c*d^5+7*a*d*c*g^6-42*a*d*c*d*g^5+105*a*d*c*d^2*g^4-140*a*d*c*d^3*g^3+105*a*d*c*d^4*g^2-42*a*d*c*d^5*g+7*a*d*c*d^6+a*c*h^7+a*c*g^7+a*c*f^7-7*a*c*e*h^6+21*a*c*e^2*h^5-35*a*c*e^3*h^4+35*a*c*e^4*h^3-21*a*c*e^5*h^2+7*a*c*e^6*h-a*c*e^7-7*a*c*d*g^6+21*a*c*d^2*g^5-35*a*c*d^3*g^4+35*a*c*d^4*g^3-21*a*c*d^5*g^2+7*a*c*d^6*g-a*c*d^7+a*b^7*c+a^8*c, e^7*b*d+e^7*a*d-7*e^6*b*e*d+7*e^6*b*d*h+21*e^5*b*e^2*d-42*e^5*b*e*d*h+21*e^5*b*d*h^2-35*e^4*b*e^3*d+105*e^4*b*e^2*d*h-105*e^4*b*e*d*h^2+35*e^4*b*d*h^3+35*e^3*b*e^4*d-140*e^3*b*e^3*d*h+210*e^3*b*e^2*d*h^2-140*e^3*b*e*d*h^3+35*e^3*b*d*h^4-21*e^2*b*e^5*d+105*e^2*b*e^4*d*h-210*e^2*b*e^3*d*h^2+210*e^2*b*e^2*d*h^3-105*e^2*b*e*d*h^4+21*e^2*b*d*h^5+7*e*b*e^6*d-42*e*b*e^5*d*h+105*e*b*e^4*d*h^2-140*e*b*e^3*d*h^3+105*e*b*e^2*d*h^4-42*e*b*e*d*h^5+7*e*b*d*h^6+d^7*a*d+c^7*a*d-b*e^7*d+7*b*e^6*d*h-21*b*e^5*d*h^2+35*b*e^4*d*h^3-35*b*e^3*d*h^4+21*b*e^2*d*h^5-7*b*e*d*h^6+a*d*h^7+a*d*g^7+a*d*f^7-7*a*d*c*f^6+21*a*d*c^2*f^5-35*a*d*c^3*f^4+35*a*d*c^4*f^3-21*a*d*c^5*f^2+7*a*d*c^6*f-a*d*c^7+7*a*c*d*f^6-42*a*c*d*c*f^5+105*a*c*d*c^2*f^4-140*a*c*d*c^3*f^3+105*a*c*d*c^4*f^2-42*a*c*d*c^5*f+7*a*c*d*c^6+21*a*c^2*d*f^5-105*a*c^2*d*c*f^4+210*a*c^2*d*c^2*f^3-210*a*c^2*d*c^3*f^2+105*a*c^2*d*c^4*f-21*a*c^2*d*c^5+35*a*c^3*d*f^4-140*a*c^3*d*c*f^3+210*a*c^3*d*c^2*f^2-140*a*c^3*d*c^3*f+35*a*c^3*d*c^4+35*a*c^4*d*f^3-105*a*c^4*d*c*f^2+105*a*c^4*d*c^2*f-35*a*c^4*d*c^3+21*a*c^5*d*f^2-42*a*c^5*d*c*f+21*a*c^5*d*c^2+7*a*c^6*d*f-7*a*c^6*d*c+a*c^7*d+a*b^7*d+a^8*d, e^7*b*e+e^7*a*e+7*e^6*b*e*h-7*e^6*b*e^2+21*e^5*b*e*h^2-42*e^5*b*e^2*h+21*e^5*b*e^3+35*e^4*b*e*h^3-105*e^4*b*e^2*h^2+105*e^4*b*e^3*h-35*e^4*b*e^4+35*e^3*b*e*h^4-140*e^3*b*e^2*h^3+210*e^3*b*e^3*h^2-140*e^3*b*e^4*h+35*e^3*b*e^5+21*e^2*b*e*h^5-105*e^2*b*e^2*h^4+210*e^2*b*e^3*h^3-210*e^2*b*e^4*h^2+105*e^2*b*e^5*h-21*e^2*b*e^6+7*e*b*e*h^6-42*e*b*e^2*h^5+105*e*b*e^3*h^4-140*e*b*e^4*h^3+105*e*b*e^5*h^2-42*e*b*e^6*h+7*e*b*e^7+d^7*a*e+c^7*a*e-7*b*e^2*h^6+21*b*e^3*h^5-35*b*e^4*h^4+35*b*e^5*h^3-21*b*e^6*h^2+7*b*e^7*h-b*e^8+a*e*h^7+a*e*g^7+a*e*f^7+a*b^7*e+a^8*e};
assert(XA_{0,1,2,3,4} == firstFewAns)
///
*}

{*
--- Check hilbert series code
TEST ///
restart
///
*}

--- Test coefficients
--- coefficients should be fixed so that the return value
--- can recreate the input in some way.
TEST ///
restart
needsPackage "NCAlgebra"
A=QQ{a, b, c, d, e, f, g, h}
F = a^6+b^6+c^6+d^6+e^6+f^6+g^6+h^6;
bas = time flatten entries basis(6,A);
#bas
coeffs = time coefficients(F,Monomials=>bas);
time sparseCoeffs({F},Monomials=>bas);
time sparseCoeffs(toList (10:F),Monomials=>bas);
///

--- Check variables with non-standard weights
TEST ///
restart
///

--- Factoring code with skew polynomial ring example/check
TEST ///
restart
needsPackage "NCAlgebra"
-- n is the number of variables in skew poly ring
n = 3
rk = 2^(n-1)
B = skewPolynomialRing(QQ,-1_QQ,{x_1..x_n})
assert(x_1*x_2 + x_2*x_1 == 0)
assert(x_1*x_3 + x_3*x_1 == 0)
assert(x_2*x_3 + x_3*x_2 == 0)
assert(ring x_1 === B)
A = ambient B
assert(ring x_1 === A)
-- tau is cyclic ring automorphism
tau = ncMap(B,B,drop(gens B,1) | {(gens B)#0})
-- create an Ore extension with a new variable w and tau
oreIdeal(B,tau,w)
C = oreExtension(B,tau,w)
assert(w*x_1-x_2*w == 0)
assert(w*x_2-x_3*w == 0)
assert(w*x_3-x_1*w == 0)
assert(ring x_1 === C)
A' = ambient C
assert(ring x_1 === A')
use A'
J = ideal C
J = J + ncIdeal {w^2}
D = A'/J
assert(ring x_1 === D)
DtoC = ncMap(C,D,gens C)
-- need inverse maps, since \varphi^\sigma is obtained by applying \sigma^{-1}
-- to all the matrix entries
tauCInv = ncMap(C,C,{(gens C)#(n-1)} | drop(drop(gens C,-1),-1) | {DtoC w})
-- have to compose twice, because the normalizing automorphism of w^2 is tau^2
sigmaInv = tauCInv @@ tauCInv
-- build the presentation of the example module
M1 = ncMatrix {drop(gens D,-2) | {w}}
assignDegrees M1
-- take a sufficiently high syzygy over D to make sure we have a MCM module
-- Note that since the algebra is Koszul, we may use the rightKernel code
-- that does not make a call to Bergman.
hiSyz = M1
scan(n, i -> (<< "Computing " << i+2 << "th syzygy" << endl; hiSyz = rightKernel(hiSyz,1)))
-- lift the map to the ambient ring
hiSyzC = DtoC hiSyz
assignDegrees hiSyzC;
assert(isHomogeneous hiSyzC)
use C
-- build f*id map
-- (this needs to be easier to make)
fId = w^2*(ncMatrix applyTable(entries id_(ZZ^rk), i -> promote(i,C)));
assignDegrees fId;
assert(isHomogeneous fId)
-- now factor through w^2*identity
hiSyzC' = fId // hiSyzC
--- check that these satisfy the twisted matrix factorization condition
assert(hiSyzC*hiSyzC' == fId)
assert(hiSyzC'*(sigmaInv hiSyzC) == fId)
///

--- Factoring map code over Sklyanin example/check
TEST ///
needsPackage "NCAlgebra"
A = QQ{x,y,z}
f1 = y*z + z*y - x^2
f2 = x*z + z*x - y^2
f3 = z^2 - x*y - y*x
I = ncIdeal {f1,f2,f3}
-- note that this is not a 'generic' sklyanin, and has a finite GB
Igb = ncGroebnerBasis I
g = 2*(-y^3-x*y*z+y*x*z+x^3)
J = I + ncIdeal {g}
B = A/I -- this is the sklyanin
B' = A/J -- factor of sklyanin
BprimeToB = ncMap(B,B',gens B)
k = ncMatrix {{x,y,z}}
assignDegrees k
M = BprimeToB rightKernelBergman rightKernelBergman k
fId = g*(ncMatrix applyTable(entries id_(ZZ^4), i -> promote(i,B)));
assignDegrees(fId,{2,2,2,3},{5,5,5,6});
assert(isHomogeneous fId)
-- now factor through g*id
M' = fId // M
assert(M*M' == fId)
assignDegrees(fId,{3,4,4,4},{6,7,7,7});
assert(M'*(M[-3]) == fId)
///

{*
--- NEED TO FIX THIS
--- Skew polynomial ring, abelianization, and isCommutative
TEST ///
needsPackage "NCAlgebra"
R = QQ[q]/ideal{q^4+q^3+q^2+q+1}
B = skewPolynomialRing(R,q,{x,y,z,w})
assert(x*y == q*y*x)
C = skewPolynomialRing(QQ,1_QQ, {x,y,z,w})
assert(isCommutative C)
assert(not isCommutative B)
abC = abelianization C
abC' = abelianization ambient C
assert(ideal abC == 0)
assert(ideal abC' == 0)
Bop = oppositeRing B
assert(y*x == q*x*y)
///
*}

end

--- Some examples are below and could be turned into part of the documentation

--- matrix factorizations over sklyanin algebra
restart
debug needsPackage "NCAlgebra"
A = QQ{x,y,z}
f1 = y*z + z*y - x^2
f2 = x*z + z*x - y^2
f3 = z^2 - x*y - y*x
I = ncIdeal {f1,f2,f3}
Igb = ncGroebnerBasis I
z*x^3 % Igb
g = -y^3-x*y*z+y*x*z+x^3
J = I + ncIdeal {g}
g % Igb
normalFormBergman(g,Igb)
B = A/I
use A
B' = A/J
use B
leftMultiplicationMap(x,3)
leftMultiplicationMap(z,3)
centralElements(B,3)
basis(3,B)
g = -y^3-x*y*z+y*x*z+x^3
---- broken?
isLeftRegular(g,6)
M=ncMatrix{{z,-x,-y},{-y,z,-x},{x,y,z}}
rightKernel(M,1)
--- again, waaay too many calls to bergman!
rightKernel(basis(1,B),10)
--- skip the next line if you want to work in the tensor algebra
h = x^2 + y^2 + z^2
isCentral h
isCentral g
M3 = ncMatrix {{x,y,z,0},
               {-y*z-2*x^2,-y*x,z*x-x*z,x},
               {x*y-2*y*x,x*z,-x^2,y},
               {-y^2-z*x,x^2,-x*y,z}}
M2 = ncMatrix {{-z*y,-x,z,y},
               {z*x-x*z,z,-y,x},
               {x*y,y,x,-z},
               {2*x*y*z-4*x^3,-2*x^2,2*y^2,2*x*y-2*y*x}}
M3*M2
M2*M3
--- checking homogeneity
assignDegrees M3
isHomogeneous 3M
assignDegrees(M3,{1,0,0,0},{2,2,2,1})
isHomogeneous M3
---
M3' = M3^2
--- can now work in quotient ring!
M3*M4
M4*M3
M3*(x*y) - (x*y)*M3
M1 = ncMatrix {{x}}
M2 = ncMatrix {{y}}
M1*M2
--- apparently it is very important to reduce your entries along the way.
wallTiming (() -> M3^6)
--- still much slower!  It seems that reducing all along the way is much more efficient.
wallTiming (() -> M3'^5)

--- or can work over free algebra and reduce later
M4*M3 % ncgb
M3*M4 % ncgb
M3' = M3 % ncgb'
M4' = M4 % ncgb'
M3'*M4' % ncgb
M4'*M3' % ncgb

M3+M4
2*M3
(g*M3 - M3*g) % ncgb
M3^4 % ncgb
wallTiming (() -> M3^4 % Igb)
---------------------------------------------

---- working in quantum polynomial ring -----
restart
needsPackage "NCAlgebra"
R = toField(QQ[q]/ideal{q^4+q^3+q^2+q+1})
A = R{x,y,z}
--- this is a gb of the poly ring skewed by a fifth root of unity.
I = ncIdeal {y*x - q*x*y, z*y - q*y*z, z*x - q*x*z}
ncgb = ncGroebnerBasis(I,InstallGB=>true)
B = A / I
-- get a basis of the degree n piece of A over the base ring
time bas = basis(10,B);
coefficients(x*y+q^2*x*z)
bas2 = flatten entries basis(2,B)
coeffs = coefficients(x*y+q^2*x*z, Monomials => bas2)
first flatten entries ((ncMatrix{bas2})*coeffs)
-- yay!
matrix {{coeffs, coeffs},{coeffs,coeffs}}
basis(2,B)
basis(3,B)
leftMultiplicationMap(x,2)
rightMultiplicationMap(x,2)
centralElements(B,4)
centralElements(B,5)
--- we can verify that f is central in this ring, for example
f = x^5 + y^5 + z^5
g = x^4 + y^4 + z^4
isCentral f
isCentral g
-- example computation
-- takes waaay too long!  Need to try and speed up reduction code.
time h = f^2
time h = f^3
time h = f^4
------------------------------------------------------

--- testing out Bergman interface
restart
debug needsPackage "NCAlgebra"
A = QQ{x,y,z}
f1 = y*z + z*y - x^2
f2 = x*z + z*x - y^2
f3 = z^2 - x*y - y*x
I = ncIdeal {f1,f2,f3}
B = A / I
hs = hilbertBergman(B,DegreeLimit=>25)

restart
debug needsPackage "NCAlgebra"
A = QQ{x,y,z}
setWeights(A,{2,2,2})
f1 = y*z + z*y - x^2
f2 = x*z + z*x - y^2
f3 = z^2 - x*y - y*x
I = ncIdeal {f1,f2,f3}
B = A/I
hs = hilbertBergman(B,DegreeLimit=>20)
--------------------------------------------

-----------
-- this doesn't work since it is not homogeneous unless you use degree q = 0, which is not allowed.
restart
needsPackage "NCAlgebra"
A = QQ{q,x,y,z}
I = ncIdeal {q^4+q^3+q^2+q+1,q*x-x*q,q*y-y*q,q*z-z*q,y*x-q*x*y,z*y-q*y*z,z*x-q*x*z}
Igb = twoSidedNCGroebnerBasisBergman I

---- ore extensions
quit
restart
needsPackage "NCAlgebra"
A = QQ{x,y,w}
I = ncIdeal {x*y-y*x, x*w-w*x-w*y,y*w-w*x,w*w}
B=A/I
M1=ncMatrix{{x,w}}
assignDegrees M1
M2=rightKernelBergman(M1)
M3=rightKernelBergman(M2)
M4=rightKernelBergman(M3,DegreeLimit=>8)
M5=rightKernelBergman(M4,DegreeLimit=>8)
quit
restart
needsPackage "NCAlgebra"
A = QQ{x,y,z,w}
I = ncIdeal {x*y+y*x,x*z+z*x,y*z+z*y,x*w-w*y,y*w-w*z,z*w-w*x,w^2}
J = ncIdeal {x*y-y*x,x*z-z*x,y*z-z*y,x*w-w*y-w*z,y*w-w*z-w*x,z*w-w*x-w*y,w^2}
B = A/I
C=A/J
M1 = ncMatrix {{x,y,w}}
M2 = rightKernel(M1,1)
M3 = rightKernel(M2,1)
M4 = rightKernel(M3,1)
M5 = rightKernel(M4,1)
M6 = rightKernel(M5,1)
M7 = rightKernel(M6,1)
M8 = rightKernel(M7,1)
M9 = rightKernel(M8,1)
Minfty = ncMatrix {{-1/2*w, 0, 0, 0}, {-z-2*x,w,-3/2*w,0},{y-x,0,-1/2*w,0},{0,y-x,z+2*x,w}}
M0 = rightKernel(Minfty,1)
N = matrix {{1,1,0},{1,0,1},{0,1,1}}
N^5
M4A = promote(M4,A)
M5A = promote(M5,A)
M6A = promote(M6,A)
M4A
M6A

---- ore extensions
restart
needsPackage "NCAlgebra"
A = QQ{x,y,z,w}
f1 = y*z + z*y - x^2
f2 = x*z + z*x - y^2
f3 = z^2 - x*y - y*x
I = ncIdeal {f1,f2,f3,x*w-w*y,y*w-w*z,z*w-w*x,w^2}
B = A/I
M1 = ncMatrix {{x,y,z,w}}
M2 = rightKernel(M1,1)
M3 = rightKernel(M2,1)
M4 = rightKernel(M3,1)
M5 = rightKernel(M4,1)
M6 = rightKernel(M5,1)
M4A = promote(M4,A)
M5A = promote(M5,A)
M6A = promote(M6,A)
M4A
M6A

---- ore extension of skew ring
restart
debug needsPackage "NCAlgebra"
A = QQ{x,y,z,w}
I = ncIdeal {x*y+y*x,x*z+z*x,y*z+z*y,x*w-w*y,y*w-w*z,z*w-w*x,w^2}
B = A/I
M1 = ncMatrix {{x,y,w}}
assignDegrees M1
-- another bug!
rightKernelBergman M1
M2 = rightKernel(M1,1)
M3 = rightKernel(M2,1)
M4 = rightKernel(M3,1)
M5 = rightKernel(M4,1)
M6 = rightKernel(M5,1)
M7 = rightKernel(M6,1)
M8 = rightKernel(M7,1)
M9 = rightKernel(M8,1)
M10 = rightKernel(M9,1)
M3A = promote(M3,A)
M4A = promote(M4,A)
M5A = promote(M5,A)
M6A = promote(M6,A)
M7A = promote(M7,A)
M8A = promote(M8,A)
M9A = promote(M9,A)
M10A = promote(M10,A)

--- make sure ores with subscripts work
restart
needsPackage "NCAlgebra"
n=4
A = QQ {x_1..x_n}
gensA = gens A;
I = ncIdeal apply(subsets(toList(0..(n-1)),2), p -> gensA#(p#0)*gensA#(p#1)+gensA#(p#1)*gensA#(p#0))
B = A/I
sigma = ncMap(B,B,drop(gens B,1) | {(gens B)#0})
matrix sigma
C = oreExtension(B,sigma,w)

---- ore extension of skew ring
restart
debug needsPackage "NCAlgebra"
n = 6
rk = 2^(n-1)
A = QQ {x_1..x_n}
gensA = gens A;
I = ncIdeal apply(subsets(toList(0..(n-1)),2), p -> gensA#(p#0)*gensA#(p#1)+gensA#(p#1)*gensA#(p#0))
B = A/I
sigma = ncMap(B,B,drop(gens B,1) | {(gens B)#0})
matrix sigma
C = oreExtension(B,sigma,w)
A' = ambient C
use A'
J = ideal C
J = J + ncIdeal {w^2}
D = A'/J
DtoC = ncMap(C,D,gens C)
sigmaCInv = ncMap(C,C,{(gens C)#(n-1)} | drop(drop(gens C,-1),-1) | {DtoC w})
M1 = ncMatrix {drop(gens D,-2) | {w}}
assignDegrees M1
hiSyz = M1
scan(n, i -> (<< "Computing " << i+2 << "th syzygy" << endl; hiSyz = rightKernel(hiSyz,1)));
hiSyzC = DtoC hiSyz
assignDegrees hiSyzC;
isHomogeneous hiSyzC
----

--- try factoring code
use C
fId = w^2*(ncMatrix applyTable(entries id_(ZZ^rk), i -> promote(i,C)));
assignDegrees fId;
isHomogeneous fId
hiSyzC' = fId // hiSyzC
hiSyzC*hiSyzC'
hiSyzC'*(sigmaCInv sigmaCInv hiSyzC)
----
M5 = rightKernel(M4,1)
M5C = DtoC M5
M4C*M5C   -- huzzah!

---- ore extension of sklyanin
restart
needsPackage "NCAlgebra"
A = QQ{x,y,z,w}
f1 = y*z + z*y - x^2
f2 = x*z + z*x - y^2
f3 = z^2 - x*y - y*x
I = ncIdeal {f1,f2,f3,x*w-w*y,y*w-w*z,z*w-w*x,w^2}
B = A/I
M1 = ncMatrix {{x,y,w}}
M2 = rightKernel(M1,1)
M2 = rightKernel(M1,2)
M2 = rightKernel(M1,3)
M2 = rightKernel(M1,4)
M3 = rightKernel(M2,1)
M4 = rightKernel(M3,1)
M5 = rightKernel(M4,1)
M6 = rightKernel(M5,1)
M7 = rightKernel(M6,1)
M8 = rightKernel(M7,1)
M9 = rightKernel(M8,1)
M10 = rightKernel(M9,1)
M3A = promote(M3,A)
M4A = promote(M4,A)
M5A = promote(M5,A)
M6A = promote(M6,A)
M7A = promote(M7,A)
M8A = promote(M8,A)
M9A = promote(M9,A)
M10A = promote(M10,A)

---- ore extension of skew ring with infinite automorphism
restart
debug needsPackage "NCAlgebra"
A = QQ{x,y,z,w}
I = ncIdeal {x*y+y*x,x*z+z*x,y*z+z*y,x*w-w*(y+z),y*w-w*(z+x),z*w-w*(x+y),w^2}
B = A/I
M1 = ncMatrix {{x,y,w}}
M2 = rightKernel(M1,1)
M3 = rightKernel(M2,1)
M4 = rightKernel(M3,1)
M5 = rightKernel(M4,1)
M6 = rightKernel(M5,1)
M7 = rightKernel(M6,1)
M8 = rightKernel(M7,1)
M9 = rightKernel(M8,1)
M10 = rightKernel(M9,1)
M11 = rightKernel(M10,1)
M12 = rightKernel(M11,1)
M13 = rightKernel(M12,1)
M14 = rightKernel(M13,1)
M15 = rightKernel(M14,1)

assignDegrees M1
M2 = rightKernelBergman M1
M3 = rightKernelBergman M2
M4 = rightKernelBergman M3


--- test for speed of reduction code
restart
needsPackage "NCAlgebra"
A = QQ{x,y,z,w}
I = ncIdeal {x*y+y*x,x*z+z*x,y*z+z*y,x*w-2*w*y,y*w-2*w*z,z*w-2*w*x,w^2}
B = A/I
M1 = ncMatrix {{x,y,w}}
time M2 = rightKernel(M1,7,Verbosity=>1);
time M3 = rightKernel(M2,3,Verbosity=>1);

---- andy's example
restart
debug needsPackage "NCAlgebra"
A=QQ{a, b, c, d, e, f, g, h}
-- this got very slow all of a sudden!
I = gbFromOutputFile(A,"NCAlgebra/UghABCgb6.txt", ReturnIdeal=>true);
Igb = ncGroebnerBasis I;
B=A/I;
M1 = ncMatrix {gens B};
M2 = rightKernel(M1,1)
--- call that takes too long!!!  Shouldn't call Bergman again for NF.
M3 = rightKernel(M2,4)

M2=ncMatrix{{-b,-f,-c,0,0,0,-g,0,0,0,0,-h,0,0,0,0},
    {a,0,c-f,0,0,0,0,d-g,0,0,0,e,e-h,0,0,0},
    {0,0,a-b,-f,-d,0,0,0,d-g,0,0,0,0,e-h,0,0},
    {0,0,0,0,c-f,0,0,-b,-c,-g,0,0,0,0,-h,0},
    {0,0,0,0,0,-f,0,0,0,0,-g,-b,-b,-c,0,-h},
    {0,a,b,c,d,e,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,a,b,c,d,e,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,a,b,c,d,e}};
M=basis(1,B);

bas=basis(2,B)
use A
f = a^7+b^7+c^7+d^7+e^7+f^7+g^7+h^7
f = promote(f,B);
--- sometimes bergman calls are great!
time X = flatten entries (f*bas);
netList X

-------------------------------
-- Testing NCRingMap code
restart
needsPackage "NCAlgebra"
A = QQ{x,y,z}
I = ncIdeal {x*y-y*x,x*z-z*x,y*z-z*y}
B = A/I
sigma = ncMap(B,B,{y,z,x})
delta = ncMap(B,B,apply(gens B, x -> promote(1,B)))
isWellDefined sigma
oreExtension(B,sigma,w)
-- doesn't really matter that we can handle derivations yet, since bergman doesn't do inhomogeneous
-- gbs very well.
oreIdeal(B,sigma,delta,w)
oreExtension(B,sigma,delta,w)
-------------------------------
--- Testing multiple rings code
restart
needsPackage "NCAlgebra"
A = QQ{x,y,z}
I = ncIdeal {x*y+y*x,x*z+z*x,y*z+z*y}
C = QQ{x,y,z}
B = A/I

-----------------------------------
--- Testing out normal element code
restart
debug needsPackage "NCAlgebra"
A = QQ{a,b,c}
I = ncIdeal {a*b+b*a,a*c+c*a,b*c+c*b}
B = A/I
sigma = ncMap(B,B,{b,c,a})
sigma_2   -- testing restriction of NCMap to degree code
isWellDefined sigma
C = oreExtension(B,sigma,w)
normalElements(C,1,x,y)
tau = ncMap(C,C,{b,c,a,w})
normalElements(tau,3)
normalElements(tau,1)
normalElements(tau @@ tau,1)
normalElements(tau @@ tau,2)
normalElements(tau @@ tau,3)
normalElements(tau @@ tau,4)
findNormalComplement(w,x)
isNormal w
phi = normalAutomorphism w
matrix phi
phi2 = normalAutomorphism w^2
matrix phi2

-- another normal element test
restart
debug needsPackage "NCAlgebra"
A = QQ{a,b,c}
I = ncIdeal {a*b+b*a,a*c+c*a,b*c+c*b}
B = A/I
normalElements(B,1,x,y)

-- another test
restart
debug needsPackage "NCAlgebra"
A = QQ{x,y,z}
I = ncIdeal {y*z + z*y - x^2,x*z + z*x - y^2,z^2 - x*y - y*x}
B = A/I
basis(2,B)
normalElements(B,2,a,b)
f = x^2+2*x*y+2*y*x+3*y^2
isNormal f
isCentral f
findNormalComplement(f,x)
findNormalComplement(f,y)
findNormalComplement(f,z)

--- one more test...
restart
debug needsPackage "NCAlgebra"
A = (QQ[a_0..a_5]) {x,y,z};
gensI = {y^2*x-x*y^2, y*x^2-x^2*y, z*x-y^2+x*z, z*y+y*z-x^2, z^2-y*x-x*y}
I = ncIdeal gensI
Igb = ncGroebnerBasis(I, InstallGB=>true)
B = A/I
f = a_0*x^2 + a_1*x*y + a_1*y*x + a_4*y^2
centralElements(B,3)
normalElements(B,3,b,c)
basis(3,B)

------------------------------------
--- Testing out endomorphism code
restart
debug needsPackage "NCAlgebra"
Q = QQ[a,b,c]
R = Q/ideal{a*b-c^2}
kRes = res(coker vars R, LengthLimit=>7);
M = coker kRes.dd_5
B = endomorphismRing(M,X);
gensI = gens ideal B;
gensIMin = minimizeRelations(gensI, Verbosity=>1)
checkHomRelations(gensIMin,B.cache.endomorphismRingGens);
------------------------------------
--- Testing out endomorphism code
restart
debug needsPackage "NCAlgebra"
Q = QQ[a,b,c,d]
R = Q/ideal{a*b+c*d}
kRes = res(coker vars R, LengthLimit=>7);
M = coker kRes.dd_5
time B = endomorphismRing(M,X);
gensI = gens ideal B;
time gensIMin = minimizeRelations(gensI, Verbosity=>1);   -- new bug here!!!
checkHomRelations(gensIMin,B.cache.endomorphismRingGens);
--------------------------------------
--- Skew group ring example?
restart
debug needsPackage "NCAlgebra"
S = QQ[x,y,z]
Q = QQ[w_1..w_6,Degrees=>{2,2,2,2,2,2}]
phi = map(S,Q,matrix{{x^2,x*y,x*z,y^2,y*z,z^2}})
I = ker phi
R = Q/I
phi = map(S,R,matrix{{x^2,x*y,x*z,y^2,y*z,z^2}})
M = pushForward(phi,S^1)
B = endomorphismRing(M,X)
gensI = gens ideal B;
gensIMin = minimizeRelations(gensI, Verbosity=>1)
checkHomRelations(gensIMin,B.cache.endomorphismRingGens);

first gensI
first newGensI

netList B.cache.endomorphismRingGens
f = sum take(flatten entries basis(1,B),10)
f^2
HomM = Hom(M,M)
map1 = homomorphism(HomM_{0})
map2 = homomorphism(HomM_{1})
gensHomM = gens HomM
elt = transpose flatten matrix (map2*map1)
gensHomM // elt

restart
debug needsPackage "NCAlgebra"
A = QQ{x,y,z}
mon = first first pairs (x*y*z).terms
substrings(mon,3)

QQ toList vars(0..14)
QQ{x_1,x_2}
QQ{a..d}

restart
debug needsPackage "NCAlgebra"
R = ZZ/101[a,b,c,d]/ideal{b^5,c^5}
cSubstrings(first exponents (a*b^2*c^3*d^4),0,4)
cSubstrings(first exponents (a*b^2*c^3*d^4),0,0)


-------------------------------------
--- Testing lead terms
restart
debug needsPackage "NCAlgebra"
A = QQ{a,b,c,d}
leadTerm(d*b+d*a*c)

-- Andy test
restart
debug needsPackage "NCAlgebra"
A = QQ {a,b,c,d,e,f,g,h}
J = ncIdeal  {b*a-a*b, f*a-a*f, f*b-c*b+c*a-b*f+b*c-a*c, f*c-c*f, f*d-d*f+d*c-c*d,
     f*e-e*f, g*a-a*g, g*b-d*b-b*g+b*d, g*c-d*c-c*g+c*d, g*d-d*g, g*e-e*g, h*a-e*b+b*e-a*h,
     h*b-e*b-b*h+b*e, h*c-e*c-c*h+c*e, h*d-d*h, h*e-e*h}
Jgb = ncGroebnerBasis(J,InstallGB=>true)
time basis(5,J)
time basis(6,J)

--- nonmonomial rels ker test -- not a good MF example, but a good kernel bug.
restart
debug needsPackage "NCAlgebra"
A = QQ{x,y,z,w}
I = ncIdeal {x*y+y*x,x*z+z*x,y*z+z*y,x*w-w*(y+z),y*w-w*(z+x),z*w-w*(x+y),w^2}
B = A/I
M1 = ncMatrix {{x,y,w}}
assignDegrees M1
M2 = rightKernelBergman M1
M3 = rightKernelBergman M2
M4 = rightKernelBergman M3
rightKernelBergman (M2_{0,1,2,3,5,7})
M2_{0,1,2,3,5,7}
rightMingens M2
isHomogeneous M2
---

--- matrix factoring test
M2a = M2_{0,1,2,3}
isHomogeneous M2a
M2b = M2_{4}
isHomogeneous M2b
M2b = M2_{2}*x + M2_{1}*y
isHomogeneous M2b
liftb = M2b // M2a
M2b - M2a*liftb
errTerm = transpose ncMatrix {{x^2,0,0}}
assignDegrees(errTerm,{1,1,1},{3})
M2c = M2_{2}*x + M2_{1}*y - errTerm
isHomogeneous M2c
liftc = M2c // M2a
M2c - M2a*liftc
M2d = M2_{4,5,6,7,8,9,10,11,12,13}
liftd = M2d // M2a
M2d - M2a*liftd

--- submatrices
M2' = M2_{0,1,2,3,4}
M2' = M2_{0,1,2,3,5}
M2' = M2_{0,1,2,3,6}
M2' = M2_{0,1,2,3,7}
M2' = M2_{0,1,2,3,8}
assignDegrees(M2', {0,0,0},{1,1,1,1,2})
isHomogeneous M2'
M3' = rightKernelBergman M2'  --- bug!?!
M4 = rightKernelBergman M3'

-- debugging
g = newGBGens#3
minKerGens2 = minKerGens | {g}
nonminKerGens = flatten apply(take(gens ring first minKerGens2,4), v -> minKerGens2 * v)
tempGb =  twoSidedNCGroebnerBasisBergman(gbIdealB | nonminKerGens | {g},
                                         "NumModuleVars" => numgens C - numgens B,
                                         DegreeLimit => opts#DegreeLimit,
                                         ClearDenominators=>true,
                                         CacheBergmanGB=>false);
apply(newGBGens, f -> f % tempGb)

--- skew poly ring
restart
debug needsPackage "NCAlgebra"
R = QQ[q]/ideal{q^4+q^3+q^2+q+1}
B = skewPolynomialRing(R,q,{x,y,z,w})
C = skewPolynomialRing(QQ,1_QQ, {x,y,z,w})
isCommutative C
abC = abelianization C
abC' = abelianization ambient C

--- opposite ring
restart
debug needsPackage "NCAlgebra"
R = QQ[q]/ideal{q^4+q^3+q^2+q+1}
B = skewPolynomialRing(R,q,{x,y,z,w})
oppositeRing B

--- check NCAlgebra
restart
needsPackage "NCAlgebra"
check "UnitTestsNCA"

