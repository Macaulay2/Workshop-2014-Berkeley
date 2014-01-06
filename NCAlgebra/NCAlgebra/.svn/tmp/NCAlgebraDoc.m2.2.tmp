undocumented {(net,NCGroebnerBasis),
              (net,NCIdeal),
	      (net,NCLeftIdeal),
	      (net,NCRightIdeal),
	      (net,NCRing),
	      (net,NCRingElement),
	      (net,NCMatrix),
	      (net,NCRingMap),
	      (expression, NCMatrix),
	      (net,NCQuotientRing),
	      functionHash,
	      (NewFromMethod,NCPolynomialRing,List),
	      (NewFromMethod,NCQuotientRing,List)}

beginDocumentation()

-------------------------
----- Types
-------------------------
    
doc ///
  Key
    NCAlgebra
  Description
    Text
      This package is used to define and manipulate noncommutative algebras.  Many of the
      commands contain calls to the existing noncommutative algebra package Bergman.
  Subnodes
    "Basic operations on noncommutative algebras"
    "Using the Bergman interface"
///

doc ///
   Key
      NCRing
   Headline
      Type of a noncommutative ring
   Description
      Text
         All noncommutative rings have this as an ancestor type.  It is the parent of the
	 types @ TO NCPolynomialRing @ and @ TO NCQuotientRing @. 
      Text
         In addition to defining a ring as a quotient of a @ TO NCPolynomialRing @, some common ways to create
	 NCRings include @ TO skewPolynomialRing @, @ TO endomorphismRing @, and @ TO oreExtension @.      
      
         Let's consider a three dimensional Sklyanin algebra.  We first define the tensor algebra:
      Example
         A = QQ{x,y,z}
      Text
         Then input the defining relations, and put them in an ideal:
      Example
	 f = y*z + z*y - x^2
	 g = x*z + z*x - y^2
	 h = z^2 - x*y - y*x
     	 I=ncIdeal{f,g,h}
      Text
         Next, define the quotient ring (as well as try a few functions on this new ring).  Note that
	 when the quotient ring is defined, a call is made to Bergman to compute the Groebner basis
	 of I (out to a certain degree, should the Groebner basis be infinite).
      Example
	 B=A/I
	 generators B
	 numgens B
	 isCommutative B
	 coefficientRing B
      Text
	 As we can see, x is an element of B.
      Example
         x
      Text
         If we define a new ring containing x, x is now part of that new ring
      Example
      	 C = skewPolynomialRing(QQ,(-1)_QQ,{x,y,z,w}) 
         x
      Text
         We can 'go back' to B using the command @ TO (use, NCRing) @.
      Example
	 use B
	 x
	 use C
      Text
         We can also create an Ore extension.  First define a @ TO NCRingMap @ with @ TO ncMap @.
      Example
	 sigma = ncMap(C,C,{y,z,w,x})
      Text
         Then call the command @ TO oreExtension @.
      Example
	 D = oreExtension(C,sigma,a)
	 generators D
	 numgens D
   SeeAlso
      "Basic operations on noncommutative algebras"
///

doc ///
  Key
    (generators, NCRing)
  Headline
    The list of algebra generators of an NCRing
  Usage
    gensA = generators A
  Inputs
    A : NCRing
  Outputs
    gensA : List
  Description
    Text
       This function returns the generators of an NCRing as a list.  As usual,
       gens is a synonym for generators.
    Example
       A = QQ{x,y,z}
       generators A
       gens A
///

doc ///
  Key
    (numgens, NCRing)
  Headline
    The number of algebra generators of an NCRing
  Usage
    numgensA = numgens A
  Inputs
    A : NCRing
  Outputs
    numgensA : ZZ
  Description
    Text
       This function returns the number of generators of an NCRing.
    Example
       A = QQ{x,y,z}
       numgens A
///

doc ///
  Key
    (isCommutative, NCRing)
    isExterior
    (isExterior, NCRing)
--    (isExterior, Ring)
  Headline
    Returns whether an NCRing is commutative
  Usage
    isComm = isCommutative A
  Inputs
    A : NCRing
  Outputs
    isComm : Boolean
  Description
    Text
       This function returns whether an NCRing is commutative
    Example
       A = QQ{x,y,z}
       isCommutative A
       B = skewPolynomialRing(QQ,(-1)_QQ,{x,y,z})
       isCommutative B
       C = skewPolynomialRing(QQ,1_QQ,{x,y,z})
       isCommutative C
///

doc ///
  Key
    (coefficientRing, NCRing)
  Headline
    Returns the base ring of an NCRing
  Usage
    k = coefficientRing NCRing
  Inputs
    A : NCRing
  Outputs
    k : Ring
  Description
    Text
       This function returns the base ring of an NCRing
    Example
       A = QQ{x,y,z}
       coefficientRing A
       R = ZZ/101[a,b,c,d]/(ideal(a^2-b^2))
       B = R{x,y,z}
       coefficientRing B
///

doc ///
  Key
     (use, NCRing)
  Headline
     Brings the variables of a particular NCRing in scope
  Usage
    use A
  Inputs
    A : NCRing
  Description
    Text
       This function brings the variables of a particular NCRing in scope.
       For an illustration:
    Example
       A = QQ{x,y,z}
       coefficientRing A
       B = skewPolynomialRing(QQ,(-1)_QQ,{x,y,z})
       x
    Text
       As you can see, at this point the interpreter treats x,y and z as elements of B.  To go back to
       A, we run the command use A:
    Example
       use A
       x
///

doc ///
   Key
      NCPolynomialRing
   Headline
      Type of a noncommutative polynomial ring
   Usage
      A = QQ{x,y}
   Description
      Text
         This is the type of a noncommutative polynomial ring over a commutative
	 ring R (i.e. a tensor algebra over R).  It has parent type @ TO NCRing @.
      Example
         A = QQ{x,y}
///

doc ///
   Key
      (ideal, NCPolynomialRing)
   Headline
      The defining ideal of an NCPolynomialRing
   Usage
      I = ideal A
   Inputs
      A : NCPolynomialRing
   Outputs
      I : NCIdeal
   Description
      Text
         This returns the defining ideal of an NCPolynomialRing, which 
	 will be the zero ideal in the noncommutative polynomial ring.
      Example
         A = QQ{x,y}
	 ideal A
///

doc ///
   Key
      NCQuotientRing
   Headline
      Type of a noncommutative ring
   Description
      Text
         This is the type of a quotient of a tensor algebra by a two-sided ideal.
    
         At this point, one cannot define quotients of quotients.
///

doc ///
   Key
     (symbol /, NCPolynomialRing, NCIdeal)
   Headline
     Construct a NCQuotientRing
   Usage
     B = A/I
   Inputs
     A : NCPolynomialRing
     I : NCIdeal
   Outputs
     B : NCQuotientRing
   Description
      Text
         This is one way to create a quotient of the tensor algebra modulo some relations.
    
         At this point, one cannot define quotients of quotients.
	 
	 If the base ring is QQ or a finite field of order p, then Bergman is called to compute a
	 Groebner basis.
	 
	 If not, then one has a couple of options.  The first is to take the defining ideal of the algebra, and provide a
	 Groebner Basis by calling @ TO ncGroebnerBasis @ with the InstallGB flag set to true.  Of course, if this generating
	 set is not a Groebner basis, then you will get incorrect results upon calls to functions like @ TO (basis, ZZ, NCRing) @.
	 
	 The alternative is to use the built in commands @ TO skewPolynomialRing @ and @ TO oreExtension @ which
	 has the same effect as above occuring behind the scenes.  Just be careful using these commands to create your
	 ring if your base ring is not a field Bergman can work with, as the generating sets created may not be a Groebner
	 basis for the defining ideal (this is more often a problem for @ TO oreExtension @ than @ TO skewPolynomialRing @).
      Example
         A = QQ{x,y,z}
         f = y*z + z*y - x^2
	 g = x*z + z*x - y^2
	 h = z^2 - x*y - y*x
     	 I=ncIdeal{f,g,h}
    	 B = A/I
         z^2
	 R = toField(QQ[a]/ideal(a^4+a^3+a^2+a+1))
	 C = skewPolynomialRing(R,a,{x,y,z})
	 y*x
///

doc ///
   Key
     (ambient, NCQuotientRing)
   Headline
     Ambient ring of an NCQuotientRing
   Usage
     A = ambient B 
   Inputs
     B : NCQuotientRing
   Outputs
     A : NCPolynomialRing
   Description
      Text
         Returns the ambient ring of an @ TO NCQuotientRing @.  
	 
	 As quotients of NCQuotientRings are added, this will return the top-level ambient ring.
	 
      Example
         B = skewPolynomialRing(QQ,(-1)_QQ,{x,y,z})
	 A = ambient B
///

doc ///
   Key
     (ideal, NCQuotientRing)
   Headline
     Defining ideal of an NCQuotientRing in its ambient ring
   Usage
     I = ideal B
   Inputs
     B : NCQuotientRing
   Outputs
     I : NCIdeal
   Description
      Text
         This returns the defining ideal of an NCQuotientRing in its ambient ring.  As of now,
	 this is always an ideal in an NCPolynomialRing, but when quotients of @ TO NCQuotientRing @s
	 are added, this will no longer be the case.
      Example
         B = skewPolynomialRing(QQ,(-1)_QQ,{x,y,z})
	 A = ambient B
	 I = ideal B
	 ring I === A
///

doc ///
   Key
      NCMatrix
   Headline
      Type of a matrix over a noncommutative ring
   Description
      Text
         This is the type of a matrix over a noncommutative ring.  These represent homomorphisms between two free modules in
	 chosen bases (whether you think of it as a map of left or right modules is up you).  Modules themselves are not
	 implemented yet in the @ TO NCAlgebra @ package, but are slated for a later release.
      Text
         Common ways to make (and use) a matrix include
      Code
         UL {TO (ncMatrix, List),
	     TO (basis, ZZ, NCRing),
	     TO (rightKernel, NCMatrix, ZZ),
	     TO (rightKernelBergman, NCMatrix)}
      Text
         Common ways to get information about matrices
      Code
         UL {TO (ring, NCMatrix),
	     TO (entries, NCMatrix)}
      Text
         Common operations on matrices:
      Code
         UL {TO (symbol +, NCMatrix,NCMatrix),
	     TO (symbol -, NCMatrix,NCMatrix),
	     TO (symbol %, NCMatrix,NCGroebnerBasis),
             TO (symbol *, NCMatrix,NCMatrix),
	     TO (symbol *, NCMatrix,NCRingElement),
	     TO (symbol *, NCMatrix,RingElement),
	     TO (symbol //, NCMatrix,NCMatrix),
	     TO (symbol _, NCMatrix,List),
	     TO (symbol ==, NCMatrix, NCMatrix),
	     TO (symbol |, NCMatrix,NCMatrix),
	     TO (symbol ||, NCMatrix,NCMatrix),
	     TO (symbol ^, NCMatrix,List),
	     TO (symbol ^, NCMatrix,ZZ),
	     }
      Text
         This is the type of a matrix with entries in an NCRing.  Many of the basic operations
	 one can perform on a @ TO Matrix @ are also allowed with an @ TO NCMatrix @, and
	 the behavior of the functions should be similar to the corresponding 'usual' command.
	 Some examples of creating and using NCMatrices are given below.
      Example
         A = QQ{a,b,c,d}
	 M = ncMatrix {{a,b,c,d}}
	 N = ncMatrix {{M,2*M,3*M},{4*M,5*M,6*M}}

         B = QQ{x,y,z}
	 f = y*z + z*y - x^2
	 g = x*z + z*x - y^2
	 h = z^2 - x*y - y*x
	 I = ncIdeal {f,g,h}
	 Igb = ncGroebnerBasis I
	 M = ncMatrix {{x, y, z}}
	 sigma = ncMap(B,B,{y,z,x})
	 N = ncMatrix {{M},{sigma M}, {sigma sigma M}}
	 Nred = N^3 % Igb
	 C = B/I
	 phi = ncMap(C,B,gens C)
	 NC = phi N
	 N3C = NC^3
	 X = NC + 3*NC
	 Y = NC | 2*NC
	 Z = X || NC
///

doc ///
   Key
      ncMatrix
      (ncMatrix,List)
   Headline
      Create an NCMatrix
   Usage
      M = ncMatrix entriesList
   Inputs
      entriesList : List
   Outputs
      M : NCMatrix
   Description
      Text
         This command creates an NCMatrix.  As with the @ TO matrix @ command, the user
	 may provide this matrix as a doubly nested list of NCRingElements, or as a
	 doubly nested list of NCMatrices.
      Example
         A = QQ{a,b,c,d}
	 M = ncMatrix {{a,b,c,d}}
	 N = ncMatrix {{M,2*M,3*M},{4*M,5*M,6*M}}
///

doc ///
   Key
      (symbol -, NCMatrix, NCMatrix)
   Headline
      Subtract NCMatrices
   Usage
      L = M - N
   Inputs
     M : NCMatrix
     N : NCMatrix
   Outputs
     L : NCMatrix
   Description
      Text
         This subtracts NCMatrices.
      Example
         A = QQ{x,y,z}
	 M = ncMatrix {{x, y, z}}
	 sigma = ncMap(A,A,{y,z,x})
	 N = ncMatrix {{M},{sigma M}, {sigma sigma M}}
	 N' = ncMatrix {{sigma sigma M}, {sigma M}, {M}}
	 N - N'
///

doc ///
   Key
      (symbol -, NCMatrix)
   Headline
      Negates NCMatrices
   Usage
     L = -M
   Inputs
     M : NCMatrix
   Outputs
     L : NCMatrix
   Description
      Text
         This negates NCMatrices.
      Example
         A = QQ{x,y,z}
	 M = ncMatrix {{x, y, z}}
	 sigma = ncMap(A,A,{y,z,x})
	 N = ncMatrix {{M},{sigma M}, {sigma sigma M}}
	 -N
///

doc ///
   Key
      (symbol +, NCMatrix, NCMatrix)
   Headline
      Add NCMatrices
   Usage
      L = M + N
   Inputs
     M : NCMatrix
     N : NCMatrix
   Outputs
     L : NCMatrix
   Description
      Text
         This adds NCMatrices.
      Example
         A = QQ{x,y,z}
	 M = ncMatrix {{x, y, z}}
	 sigma = ncMap(A,A,{y,z,x})
	 N = ncMatrix {{M},{sigma M}, {sigma sigma M}}
	 N' = ncMatrix {{sigma sigma M}, {sigma M}, {M}}
	 N + N'
///

doc ///
   Key
      (symbol %, NCMatrix, NCGroebnerBasis)
   Headline
      Reduces the entries of an NCMatrix with respect to an NCGroebnerBasis
   Usage
      L = M % Igb
   Inputs
     M : NCMatrix
     Igb : NCGroebnerBasis
   Outputs
     L : NCMatrix
   Description
      Text
         This command reduces the entries of an NCMatrix with respect to an NCGroebnerBasis.
      Example
         A = QQ{x,y,z}
	 f = y*z + z*y - x^2
	 g = x*z + z*x - y^2
	 h = z^2 - x*y - y*x
	 I = ncIdeal {f,g,h}
	 Igb = ncGroebnerBasis I
	 M = ncMatrix {{x, y, z}}
	 sigma = ncMap(A,A,{y,z,x})
	 N = ncMatrix {{M},{sigma M}, {sigma sigma M}}
	 N3 = N^3
	 N3red = N3 % Igb
///

doc ///
   Key
      (symbol *, NCMatrix, NCMatrix)
   Headline
      Product of NCMatrices
   Usage
      L = M*N
   Inputs
     M : NCMatrix
     N : NCMatrix
   Outputs
     L : NCMatrix
   Description
      Text
         This command allows for the product of composable NCMatrices (or ordinary matrices over the base).
      Example
         A = QQ{x,y,z}
	 B = A/ncIdeal{y*z + z*y - x^2, x*z + z*x - y^2, z^2 - x*y - y*x}
	 M = ncMatrix {{x, y, z}}
	 sigma = ncMap(B,B,{y,z,x})
	 N = ncMatrix {{M},{sigma M}, {sigma sigma M}}
	 N' = ncMatrix {{sigma sigma M}, {sigma M}, {M}}
	 N*N'
	 N'*N
///

doc ///
   Key
      (symbol *, NCMatrix, Matrix)
   Headline
      Product of NCMatrices
   Usage
      L = M*N
   Inputs
     M : NCMatrix
     N : Matrix
   Outputs
     L : NCMatrix
   Description
      Text
         This command allows for the product of composable NCMatrices (or ordinary matrices over the base).
      Example
         A = QQ{x,y,z}
	 B = A/ncIdeal{y*z + z*y - x^2, x*z + z*x - y^2, z^2 - x*y - y*x}
	 M = ncMatrix {{x, y, z}}
	 sigma = ncMap(B,B,{y,z,x})
	 N = ncMatrix {{M},{sigma M}, {sigma sigma M}}
         L = map(QQ^3,QQ^3,{{2,0,0},{1,2,0},{1,2,3}})
	 N*L
///

doc ///
   Key
      (symbol *, NCMatrix, NCRingElement)
   Headline
      Product of NCMatrices
   Usage
      L = M*f
   Inputs
     M : NCMatrix
     f : NCRingElement
   Outputs
     L : NCMatrix
   Description
      Text
         This command allows for the scalar multiplication of an NCMatrix by an NCRingElement on the right.
      Example
         A = QQ{x,y,z}
	 B = A/ncIdeal{y*z + z*y - x^2, x*z + z*x - y^2, z^2 - x*y - y*x}
	 M = ncMatrix {{x, y, z}}
	 sigma = ncMap(B,B,{y,z,x})
	 N = ncMatrix {{M},{sigma M}, {sigma sigma M}}
	 N*x^2
///

doc ///
   Key
      (symbol *, NCMatrix, RingElement)
   Headline
      Product of NCMatrices
   Usage
      L = M*f
   Inputs
     M : NCMatrix
     N : RingElement
   Outputs
     L : NCMatrix
   Description
      Text
         This command allows for the scaling of an NCMatrix by an element in the base ring.
      Example
         R = frac(QQ[a])
	 B = skewPolynomialRing(R,a,{x,y,z})
	 M = ncMatrix {{x, y, z}}
	 sigma = ncMap(B,B,{y,z,x})
	 N = ncMatrix {{M},{sigma M}, {sigma sigma M}}
    	 N*a
///

doc ///
   Key
      (symbol *, NCMatrix, QQ)
   Headline
      Product of NCMatrices
   Usage
      L = M*a
   Inputs
     M : NCMatrix
     a : QQ
   Outputs
     L : NCMatrix
   Description
      Text
         This command allows for the scaling of an @ TO NCMatrix @ by an element in @ TO QQ @.
      Example
         A = QQ{x,y,z}
	 B = A/ncIdeal{y*z + z*y - x^2, x*z + z*x - y^2, z^2 - x*y - y*x}
	 M = ncMatrix {{x, y, z}}
	 sigma = ncMap(B,B,{y,z,x})
	 N = ncMatrix {{M},{sigma M}, {sigma sigma M}}
	 N*(1/2)
///

doc ///
   Key
      (symbol *, NCMatrix, ZZ)
   Headline
      Product of NCMatrices
   Usage
      L = M*a
   Inputs
     M : NCMatrix
     a : ZZ
   Outputs
     L : NCMatrix
   Description
      Text
         This command allows for the scaling of an @ TO NCMatrix @ by an element in @ TO ZZ @.
      Example
         A = QQ{x,y,z}
	 B = A/ncIdeal{y*z + z*y - x^2, x*z + z*x - y^2, z^2 - x*y - y*x}
	 M = ncMatrix {{x, y, z}}
	 sigma = ncMap(B,B,{y,z,x})
	 N = ncMatrix {{M},{sigma M}, {sigma sigma M}}
	 N*3
///

doc ///
   Key
      (symbol *, Matrix, NCMatrix)
   Headline
      Product of NCMatrices
   Usage
      L = N*M
   Inputs
     N : Matrix
     M : NCMatrix
   Outputs
     L : NCMatrix
   Description
      Text
         This command allows for the product of composable NCMatrices (or ordinary matrices over the base).
      Example
         A = QQ{x,y,z}
	 B = A/ncIdeal{y*z + z*y - x^2, x*z + z*x - y^2, z^2 - x*y - y*x}
	 M = ncMatrix {{x, y, z}}
	 sigma = ncMap(B,B,{y,z,x})
	 N = ncMatrix {{M},{sigma M}, {sigma sigma M}}
         L = map(QQ^3,QQ^3,{{2,0,0},{1,2,0},{1,2,3}})
	 L*N
///

doc ///
   Key
      (symbol *, NCRingElement, NCMatrix)
   Headline
      Product of NCMatrices
   Usage
      L = f*M
   Inputs
     f : NCRingElement
     M : NCMatrix
   Outputs
     L : NCMatrix
   Description
      Text
         This command allows for the scalar multiplication of an NCMatrix by an NCRingElement on the left.
      Example
         A = QQ{x,y,z}
	 B = A/ncIdeal{y*z + z*y - x^2, x*z + z*x - y^2, z^2 - x*y - y*x}
	 M = ncMatrix {{x, y, z}}
	 sigma = ncMap(B,B,{y,z,x})
	 N = ncMatrix {{M},{sigma M}, {sigma sigma M}}
	 x^2*N
///

doc ///
   Key
      (symbol *, RingElement, NCMatrix)
   Headline
      Product of NCMatrices
   Usage
      L = f*M
   Inputs
     f : RingElement
     M : NCMatrix
   Outputs
     L : NCMatrix
   Description
      Text
         This command allows for the scaling of an NCMatrix by an element in the base ring.
      Example
         R = frac(QQ[a])
	 B = skewPolynomialRing(R,a,{x,y,z})
	 M = ncMatrix {{x, y, z}}
	 sigma = ncMap(B,B,{y,z,x})
	 N = ncMatrix {{M},{sigma M}, {sigma sigma M}}
    	 a*N
///

doc ///
   Key
      (symbol *, QQ, NCMatrix)
   Headline
      Product of NCMatrices
   Usage
      L = a*M
   Inputs
     a : QQ
     M : NCMatrix
   Outputs
     L : NCMatrix
   Description
      Text
         This command allows for the scaling of an @ TO NCMatrix @ by an element in @ TO QQ @.
      Example
         A = QQ{x,y,z}
	 B = A/ncIdeal{y*z + z*y - x^2, x*z + z*x - y^2, z^2 - x*y - y*x}
	 M = ncMatrix {{x, y, z}}
	 sigma = ncMap(B,B,{y,z,x})
	 N = ncMatrix {{M},{sigma M}, {sigma sigma M}}
	 (1/2)*N
///

doc ///
   Key
      (symbol *, ZZ, NCMatrix)
   Headline
      Product of NCMatrices
   Usage
      L = a*M
   Inputs
     a : ZZ
     M : NCMatrix
   Outputs
     L : NCMatrix
   Description
      Text
         This command allows for the scaling of an @ TO NCMatrix @ by an element in @ TO ZZ @.
      Example
         A = QQ{x,y,z}
	 B = A/ncIdeal{y*z + z*y - x^2, x*z + z*x - y^2, z^2 - x*y - y*x}
	 M = ncMatrix {{x, y, z}}
	 sigma = ncMap(B,B,{y,z,x})
	 N = ncMatrix {{M},{sigma M}, {sigma sigma M}}
	 3*N
///

doc ///
   Key
      (symbol ==, NCMatrix, NCMatrix)
      (symbol ==, NCMatrix, ZZ)
      (symbol ==, ZZ, NCMatrix)
   Headline
      Test equality of matrices
   Usage
      isEqual = M == N
   Inputs
      M : NCMatrix
      N : NCMatrix
   Outputs
      isEqual : Boolean
   Description
      Text
         This command tests equality for matrices.  If one of the inputs is an integer, then the test
	 only will work if the integer is zero.  Below, we test the well-definedness of the exponentiation
	 operation using Groebner bases.
      Example
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
	 (phi Nred) == N3B
///

doc ///
   Key
      (symbol |, NCMatrix, NCMatrix)
   Headline
      Join NCMatrices horizontally
   Usage
      L = M | N
   Inputs
      M : NCMatrix
      N : NCMatrix
   Outputs
      L : NCMatrix
   Description
      Text
         This command joins NCMatrices horizontally.
      Example
         A = QQ{a,b,c,d}
	 M = ncMatrix {{a,b,c,d},{b,c,d,a}}
	 N = M | 2*M | -3*M
///

doc ///
   Key
      (symbol ||, NCMatrix, NCMatrix)
   Headline
      Join NCMatrices vertically
   Usage
      L = M || N
   Inputs
      M : NCMatrix
      N : NCMatrix
   Outputs
      L : NCMatrix
   Description
      Text
         This command joins NCMatrices vertically.
      Example
         A = QQ{a,b,c,d}
	 M = ncMatrix {{a,b,c,d},{b,c,d,a}}
	 N = M || 2*M || -3*M
///

doc ///
   Key
      (symbol ^, NCMatrix, List)
   Headline
      Select some rows of an NCMatrix
   Usage
      L = M^rows
   Inputs
      M : NCMatrix
      rows : List
   Outputs
      L : NCMatrix
   Description
      Text
         This command selects some rows of an NCMatrix.
      Example
         A = QQ{a,b,c,d}
	 M = ncMatrix {{a,b,c,d},{b,c,d,a}}
	 N = M || 2*M || -3*M
	 N^{0,3,4}
///

doc ///
   Key
      (symbol _, NCMatrix, List)
   Headline
      Select some columns of an NCMatrix
   Usage
      L = M_cols
   Inputs
      M : NCMatrix
      cols : List
   Outputs
      L : NCMatrix
   Description
      Text
         This command selects some columns of an NCMatrix.
      Example
         A = QQ{a,b,c,d}
	 M = ncMatrix {{a,b,c,d},{b,c,d,a}}
	 N = M || 2*M || -3*M
	 N_{0,2}
///

doc ///
   Key
      (symbol ^, NCMatrix, ZZ)
   Headline
      Exponentiate an NCMatrix
   Usage
      L = M^n
   Inputs
      M : NCMatrix
      n : ZZ
   Outputs
      L : NCMatrix
   Description
      Text
         This exponentiates an NCMatrix.  It should be remarked that the matrix is reduced
	 with the GB of the ring it is over on each iteration of the product.  If your algebra
	 is significantly smaller than the tensor algebra, this is a large savings.
	 The input is assumed to be a nonnegative integer at this time.
      Example
         A = QQ{x,y,z}
	 M = ncMatrix {{x, y, z}}
	 sigma = ncMap(A,A,{y,z,x})
	 N = ncMatrix {{M},{sigma M}, {sigma sigma M}}
	 N^3
	 B = A/ncIdeal{y*z + z*y - x^2, x*z + z*x - y^2, z^2 - x*y - y*x}
	 NB = promote(N,B)
	 NB^3
///

doc ///
   Key
      (symbol //, NCMatrix, NCMatrix)
   Headline
      Factor one map through another
   Usage
      L = M // N
   Inputs
      M : NCMatrix
      N : NCMatrix
   Outputs
      L : NCMatrix
   Description
      Text
         This command factors one map through another.  One nice application
	 of this is to compute twisted matrix factorizations.  
	 If the maps input are homogeneous, then the degrees must match up for the command to work.

         If M does not factor through N, then the return value L is such that M - N*L is the reduction
	 of M modulo a Groebner basis for the image of N.

	 Here is an example of doing so over a PI Sklyanin algebra.
      Example
         A = QQ{x,y,z}
	 g = 2*(-y^3-x*y*z+y*x*z+x^3)
	 B = A/ncIdeal{y*z + z*y - x^2, x*z + z*x - y^2, z^2 - x*y - y*x} -- this is the Sklyanin
	 J = (ideal B) + ncIdeal {g}
	 B' = A/J -- Factor of sklyanin
	 BprimeToB = ncMap(B,B',gens B) -- way to lift back from B' to B
	 k = ncMatrix {{x,y,z}}
	 assignDegrees k
	 M = BprimeToB rightKernelBergman rightKernelBergman k  -- second syzygy of k over B
      Text
         At this point, M is maximal Cohen-Macaulay B'-module,
	 and hence the projective dimension of M as a B-module
	 is 1.  Since M is a B' module, multiplication by g on the
	 complex that gives the resolution over B is null homotopic.  This means
	 we may factor the map f through f times the identity.  We do so below.
      Example
	 gId = g*(ncMatrix applyTable(entries id_(ZZ^4), i -> promote(i,B)))
	 assignDegrees(gId,{2,2,2,3},{5,5,5,6});
	 -- now factor through g*id
	 M' = gId // M
	 M*M' == gId
///

doc ///
   Key
      (transpose, NCMatrix)
   Headline
      Transposes an NCMatrix
   Usage
      L = transpose M
   Inputs
      M : NCMatrix
   Outputs
      L : NCMatrix
   Description
      Text
         This command transposes an NCMatrix
      Example
         A = QQ{a,b,c,d}
	 M = ncMatrix {{a,b,c,d},{b,c,d,a}}
	 N = M || 2*M || -3*M
	 transpose N
///

doc ///
   Key
      (ring, NCMatrix)
   Headline
      Gives the ring of the NCMatrix
   Usage
      L = ring M
   Inputs
      M : NCMatrix
   Outputs
      L : NCRing
   Description
      Text
         This command returns the ring over which the NCMatrix is defined.
      Example
         A = QQ{a,b,c,d}
	 M = ncMatrix {{a,b,c,d},{b,c,d,a}}
    	 ring M
///

doc ///
   Key
      (entries, NCMatrix)
   Headline
      Returns the entries of the NCMatrix
   Usage
      L = entries M
   Inputs
      M : NCMatrix
   Outputs
      L : List
   Description
      Text
         Returns the entries of the NCMatrix as a doubly nested list.
      Example
         A = QQ{a,b,c,d}
	 M = ncMatrix {{a,b,c,d},{b,c,d,a}}
	 N = M || 2*M || -3*M
	 entries N
///

doc ///
   Key
      (lift, NCMatrix)
   Headline
      Lifts an NCMatrix
   Usage
      L = lift M
   Inputs
      M : NCMatrix
   Outputs
      L : NCMatrix
   Description
      Text
         This command lifts an NCMatrix to a matrix over its @ TO ambient @ NCRing.
      Example
         A = QQ{x,y,z,w}
	 B = A/ncIdeal{y*z + z*y - x^2, x*z + z*x - y^2, z^2 - x*y - y*x}
	 M = ncMatrix {{x,y,z,w},{y,z,w,x}}
	 N = M || 2*M || -3*M
	 ring N
	 ring lift N
///

doc ///
   Key
      NCRingElement
   Headline
      Type of an element in a noncommutative ring
   --Usage
   --Inputs
   --Outputs
   Description
      Text
        This is the type of an element in a noncommutative graded ring.  One can deal with these elements
	in much the same way as in the commutative case.  See @ TO RingElement @ for details.
///

doc ///
   Key
      (degree, NCRingElement)
   Headline
      Returns the degree of an NCRingElement
   Usage
     d = degree f
   Inputs
     f : NCRingElement
   Outputs
     d : ZZ
   Description
      Text
        Returns the degree of an NCRingElement.  At the moment, multigraded NCRings are not supported.
      Example
        A = QQ{x,y,z,w}
        f = x^2+y^2+z^3
	degree f
        isHomogeneous f
	setWeights(A,{3,3,2,1})
	degree f
	isHomogeneous f
///

doc ///
   Key
      (ring, NCRingElement)
   Headline
      Returns the NCRing of an NCRingElement
   Usage
     A = ring f
   Inputs
     f : NCRingElement
   Outputs
     A : NCRing
   Description
      Text
        Returns the ring of an NCRingElement
      Example
        A = QQ{x,y,z,w}
        f = x^2+y^2+z^3
        ring f
///

doc ///
   Key
      (terms, NCRingElement)
   Headline
      Returns the terms of an NCRingElement
   Usage
     t = terms f
   Inputs
     f : NCRingElement
   Outputs
     t : List
   Description
      Text
        Returns the list of terms that make up the NCRingElement.  It is a list of
	NCRingElements.
      Example
        A = QQ{x,y,z,w}
        f = 2*x^2+y^2+z^3
        t = terms f
        first t
///

doc ///
   Key
      (size, NCRingElement)
   Headline
      Returns the number of terms of an NCRingElement
   Usage
     n = size f
   Inputs
     f : NCRingElement
   Outputs
     n : ZZ
   Description
      Text
        Returns the number of terms of an NCRingElement.
      Example
        A = QQ{x,y,z,w}
        f = x^2+y^2+z^3
        size f
///

doc ///
   Key
      (support, NCRingElement)
   Headline
      Returns the variables appearing in the NCRingElement
   Usage
     sup = support f
   Inputs
     f : NCRingElement
   Outputs
     sup : List
   Description
      Text
        Returns the variables appearing in f (as elements of the ring of f).
      Example
        A = QQ{x,y,z,w}
        f = x^2+y^2+z^3
        sup = support f
	first sup
///

doc ///
   Key
      (monomials, NCRingElement)
   Headline
      Returns the monomials appearing in the NCRingElement
   Usage
     mons = support f
   Inputs
     f : NCRingElement
   Outputs
     mons : NCMatrix
   Description
      Text
        Returns the monomials appearing in NCRingElement as an NCMatrix.
      Example
        A = QQ{x,y,z,w}
        f = x^2+y^2+z^3
        mons = monomials f
///

doc ///
   Key
      (leadMonomial, NCRingElement)
   Headline
      Returns the lead monomial of an NCRingElement
   Usage
     mon = leadMonomial f
   Inputs
     f : NCRingElement
   Outputs
     mon : NCRingElement
   Description
      Text
        Returns the lead monomial of an NCRingElement (as an NCRingElement).
      Example
        A = QQ{x,y,z,w}
        f = x^2+y^2+z^3
        mon = leadMonomial f
///

doc ///
   Key
      (leadCoefficient, NCRingElement)
   Headline
      Returns the lead monomial of an NCRingElement
   Usage
     coeff = leadCoefficient f
   Inputs
     f : NCRingElement
   Outputs
     coeff : RingElement
   Description
      Text
        Returns the lead coefficient of an NCRingElement (as an element of the base).
      Example
        A = QQ{x,y,z,w}
        f = x^2+y^2+z^3
        coeff = leadCoefficient f
///

doc ///
   Key
      (leadTerm, NCRingElement)
   Headline
      Returns the lead term of an NCRingElement
   Usage
     coeff = leadTerm f
   Inputs
     f : NCRingElement
   Outputs
     coeff : NCRingElement
   Description
      Text
        Returns the lead term of an NCRingElement (as an NCRingElement).
      Example
        A = QQ{x,y,z,w}
        f = x^2+y^2+2*z^3
        coeff = leadTerm f
///

doc ///
   Key
     (isConstant, NCRingElement)
   Headline
     Returns whether the NCRingElement is constant
   Usage
     t = isConstant f
   Inputs
     f : NCRingElement
   Outputs
     t : Boolean
   Description
      Text
        Returns whether the NCRingElement is constant.
      Example
        A = QQ{x,y,z,w}
        f = x^2+y^2+2*z^3
        isConstant f
	g = promote(1,A)
	isConstant g
///

doc ///
   Key
     (baseName, NCRingElement)
   Headline
     Returns the base name of a generator of an NCRing
   Usage
     name = baseName f
   Inputs
     f : NCRingElement
   Description
      Text
        Returns the base name of a generator of an NCRing.  This could be something of type
	@ TO IndexedVariable @ or a @ TO Symbol @.
      Example
        A = QQ{x,y,z,w}
        baseName x
	B = QQ{p_1..p_6}
	baseName p_1
///

doc ///
   Key
      (toString, NCRingElement)
   Headline
      Converts an NCRingElement to a string
   Usage
     str = toString f
   Inputs
     f : NCRingElement
   Outputs
     str : String
   Description
      Text
        Converts an NCRingElement to a string.  This should be readable by both Macaulay2 as well as Bergman.
      Example
        A = QQ{x,y,z,w}
        f = x^2+y^2+z^3
        toString f
///

doc ///
   Key
      (symbol *, List, NCRingElement)
   Headline
      Scales a list by an NCRingElement on the right
   Usage
     xsf = xs*f
   Inputs
     xs : List
     f : NCRingElement
   Outputs
     xsf : List
   Description
      Text
        Scales a list by an NCRingElement on the right.
      Example
        A = QQ{x,y}
        f = x^2+y^2
        bas = flatten entries basis(3,A)
	bas*f
///

doc ///
   Key
      (symbol *, NCRingElement, List)
   Headline
      Scales a list by an NCRingElement on the left
   Usage
     fxs = f*xs
   Inputs
     f : NCRingElement
     xs : List
   Outputs
     xsf : List
   Description
      Text
        Scales a list by an NCRingElement on the right.
      Example
        A = QQ{x,y}
        f = x^2+y^2
        bas = flatten entries basis(3,A)
	f*bas
///

doc ///
   Key
      NCGroebnerBasis
   Headline
      Type of a Groebner basis for an NCIdeal in an NCRing.
   Description
     Text
       This is the type for a Groebner basis of an ideal in the tensor algebra.
       One can provide one using the @ TO InstallGB @ option of @ TO ncGroebnerBasis @
       if you happen to know it.
       
       One also canhave Macaulay2 call Bergman and have it computed via the function
       @ TO twoSidedNCGroebnerBasisBergman @.  This command is automatically called when defining
       a quotient ring, if the defining ideal does not yet have a cached Groebner basis.
       
       You can also install one from a Bergman output file if you have that handy; see
       @ TO gbFromOutputFile @.
       
       Below are a couple of examples.
     Example
       R = QQ[a,b,c,d]/ideal{a*b+c*d}
       A = R {x,y,z}
       I = ncIdeal {a*x*y,b*z^2}
       Igb = ncGroebnerBasis(I, InstallGB=>true)
     Text
       Note that after the InstallGB flag is set, no checking is done
       to ensure that the input is in fact a Groebner basis.
     Example
       c*z^2 % Igb 
       b*z^2 % Igb

       A = QQ{x,y,z}
       p = y*z + z*y - x^2
       q = x*z + z*x - y^2
       r = z^2 - x*y - y*x
       I = ncIdeal {p,q,r}
       Igb = ncGroebnerBasis I
       normalFormBergman(z^17,Igb)

     Text
       stuff
///

doc ///
   Key
      [ncGroebnerBasis,InstallGB]
      InstallGB
   Headline
      Install a NCGroebnerBasis (without verifying that it is one).
   Usage
     Igb = ncGroebnerBasis(I, InstallGB=>true)
   Inputs
     I : NCIdeal
     InstallGB : Boolean
   Outputs
     Igb : NCGroebnerBasis
///

doc ///
   Key
      [ncGroebnerBasis,DegreeLimit]
   Headline
      Degree limit on generators of a Groebner basis.
   Usage
     Igb = ncGroebnerBasis(I, DegreeLimit=>n)
   Inputs
     I : NCIdeal
     DegreeLimit : ZZ
   Outputs
     Igb : NCGroebnerBasis
///

doc ///
   Key
      ncGroebnerBasis
      (ncGroebnerBasis,List)
      (ncGroebnerBasis,NCIdeal)
   Headline
      Compute a noncommutative Groebner basis.
   Usage
     Igb = ncGroebnerBasis I
   Inputs
     I : NCIdeal
     DegreeLimit : ZZ
     InstallGB : Boolean
   Outputs
     Igb : NCGroebnerBasis
   Description
     Example
       R = QQ[a,b,c,d]/ideal{a*b+c*d}
       A = R {x,y,z}
       I = ncIdeal {a*x*y,b*z^2}
       Igb = ncGroebnerBasis(I, InstallGB=>true)
       c*z^2 % Igb 
       b*z^2 % Igb
     Text
       Note that after the InstallGB flag is set, no checking is done
       to ensure that the input is in fact a Groebner basis.
     Example
       A = QQ{x,y,z}
       p = y*z + z*y - x^2
       q = x*z + z*x - y^2
       r = z^2 - x*y - y*x
       I = ncIdeal {p,q,r}
       Igb = ncGroebnerBasis I
     Text
       If the InstallGB flag is not set, then a call to Bergman is made, if the base ring is @ TO QQ @ or a finite field
       of characteristic p.  Otherwise, an error is raised.
       
       Now we can do things with an ncgb, like compute normal forms (using the Bergman interface).
     Example
       normalFormBergman(z^17,Igb)
     Text
       Or using the built in reduction code:
     Example
       z^17 % Igb
     Text
       Calls to Bergman are usually faster, except for when the polynomial is small.  See the documentation
       for @ TO (symbol %, NCRingElement, NCGroebnerBasis) @ for details on controlling when Bergman is called.
///

doc ///
   Key
      CacheBergmanGB
   Headline
      Whether or not to cache the gb from Bergman to a file for later use
   Usage
      Igb = gbFromOutputFile(A,fileName, CacheBergmanGB=>true)
   Inputs
      A : NCPolynomialRing
      fileName : String
      CacheBergmanGB : Boolean
   Outputs
      Igb : NCGroebnerBasis
///

doc ///
   Key
      MakeMonic
   Headline
      An option that specifies Bergman output should be made monic
///

doc ///
   Key
      ReturnIdeal
   Headline
      An option that specifies certain NCAlgebra functions should return an ideal rather than a Groebner basis.
///

doc ///
   Key
      gbFromOutputFile
      (gbFromOutputFile,NCPolynomialRing,String)
      [gbFromOutputFile,CacheBergmanGB]
      [gbFromOutputFile,ReturnIdeal]
      [gbFromOutputFile,MakeMonic]
   Headline
      Read in a NCGroebnerBasis from a Bergman output file.
   Usage
      Igb = gbFromOutputFile(A,fileName)
   Inputs
      A : NCPolynomialRing
      fileName : String
      CacheBergmanGB : Boolean
      MakeMonic : Boolean
      ReturnIdeal : Boolean
   Outputs
      Igb : NCGroebnerBasis
   Description
      Text
        This command reads in an Groebner basis from a Bergman output file.
	It can be useful if you have performed a lengthy computation before,
	and wish to load in a previously computed result to do some computations.
	
	The below code will work in your Macaulay2 session, as long as the file
	mentioned below is on your path.
	A=QQ{a, b, c, d, e, f, g, h}
        --- This doesn't work in help generator because it can't find the file.
	-- I = gbFromOutputFile(A,"NCAlgebra/UghABCgb6.txt", ReturnIdeal=>true)
	--- This is not the 'real' ideal, but just to prevent an error upon running examples
	I = ncIdeal {a^9,b^9,c^9,d^9,e^9,f^9,g^9,h^9}
	B=A/I
	F = a^7+b^7+c^7+d^7+e^7+f^7+g^7+h^7
	bas=basis(2,B)
	X = flatten entries (F*bas)
	XA = apply(X, x -> promote(x,A))
	use A
	XA_{0,1,2,3,4}
///

doc ///
  Key
    (generators, NCGroebnerBasis)
  Headline
    The list of algebra generators of an NCGroebnerBasis
  Usage
    gensIgb = generators Igb
  Inputs
    Igb : NCGroebnerBasis
  Outputs
    gensIgb : List
  Description
    Text
       This function returns the generators of an NCGroebnerBasis as a list.  As usual,
       gens is a synonym for generators.
    Example
       A = QQ{x,y,z}
       p = y*z + z*y - x^2
       q = x*z + z*x - y^2
       r = z^2 - x*y - y*x
       I = ncIdeal {p,q,r}
       Igb = ncGroebnerBasis I
       gens Igb
///

doc ///
   Key
      NumModuleVars
   Headline
      An option specifying the number of module variables in the ring of the Groebner basis.
///

doc ///
   Key
      (symbol %, NCRingElement, NCGroebnerBasis)
      (symbol %, QQ, NCGroebnerBasis)
      (symbol %, ZZ, NCGroebnerBasis)
   Headline
      Reduces a NCRingElement by a NCGroebnerBasis
   Usage
     fred = f % Igb
   Inputs
     f : NCRingElement
     Igb : NCGroebnerBasis
   Outputs
     fred : NCRingElement
   Description
     Text
       This command reduces the input modulo a noncommutative Groebner basis.
       It will either reduce it using top-level Macaulay code, or via a call to
       Bergman, depending on the size and degree of the input element.
     Example
       A = QQ{x,y,z}
       p = y*z + z*y - x^2
       q = x*z + z*x - y^2
       r = z^2 - x*y - y*x
       I = ncIdeal {p,q,r}
       Igb = ncGroebnerBasis I
       z^6 % Igb
///

doc ///
   Key
      twoSidedNCGroebnerBasisBergman
      (twoSidedNCGroebnerBasisBergman,List)
      (twoSidedNCGroebnerBasisBergman,NCIdeal)
      [twoSidedNCGroebnerBasisBergman,DegreeLimit]
      [twoSidedNCGroebnerBasisBergman,NumModuleVars]
      [twoSidedNCGroebnerBasisBergman,CacheBergmanGB]
      [twoSidedNCGroebnerBasisBergman,MakeMonic]
   Headline
      Calls Bergman to compute a two sided noncommutative Groebner Basis.
   Usage
      Igb = twoSidedNCGroebnerBasisBergman I
   Inputs
      I : NCIdeal
      DegreeLimit : ZZ
      NumModuleVars : ZZ
      CacheBergmanGB : Boolean
      MakeMonic : Boolean
   Outputs
      Igb : NCGroebnerBasis
   Description
     Text
        This command calls the computer algebra system Bergman to
	compute a noncommutative Groebner basis.
	
	Since Groebner bases in the tensor algebra need not be
	finitely generated, one should specify a degree limit on the
	computation unless one has a reason to believe the Groebner
	basis is finite.
     Example
       A = QQ{x,y,z}
       p = y*z + z*y - x^2
       q = x*z + z*x - y^2
       r = z^2 - x*y - y*x
       I = ncIdeal {p,q,r}
       Igb = twoSidedNCGroebnerBasisBergman I
///

doc ///
   Key
      NCLeftIdeal
   Headline
      Type of a left ideal in a noncommutative ring
   Description
      Text
         This defines a left ideal in a noncommutative ring.  Not much
	 can be done with these objects at this point (as one can tell
	 by the dearth of operations that take an NCLeftIdeal as
	 input), , but eventually it will be a 'fully featured'
	 object.
///

doc ///
   Key
      ncLeftIdeal
      (ncLeftIdeal, List)
      (ncLeftIdeal, NCRingElement)
   Headline
      Define a left ideal in a noncommutative ring
   Usage
      I = ncLeftIdeal fs
   Inputs
      fs : List
   Outputs
      I : NCLeftIdeal      
   Description
      Text
         This defines a left ideal in a noncommutative ring.  Not much
	 can be done with these objects at this point (as one can tell
	 by the dearth of operations that take an NCLeftIdeal as
	 input), but eventually it will be a 'fully featured'
	 object.
      Example
         A = QQ{x,y,z}
         p = y*z + z*y - x^2
         q = x*z + z*x - y^2
         r = z^2 - x*y - y*x
         I = ncLeftIdeal{p,q,r}
///

doc ///
   Key
      (generators, NCLeftIdeal)
   Headline
      Returns the generators of an NCLeftIdeal
   Usage
      gensI = generators I
   Inputs
      I : NCLeftIdeal      
   Outputs
      gensI : List
   Description
      Text
         Returns the generators of an NCLeftIdeal.
      Example
         A = QQ{x,y,z}
         p = y*z + z*y - x^2
         q = x*z + z*x - y^2
         r = z^2 - x*y - y*x
         I = ncLeftIdeal{p,q,r}
	 gens I
///

doc ///
   Key
      (ring, NCLeftIdeal)
   Headline
      Returns the ring of an NCLeftIdeal
   Usage
      A = ring I
   Inputs
      I : NCLeftIdeal      
   Outputs
      A : NCRing
   Description
      Text
         Returns the ring of an NCLeftIdeal.
      Example
         A = QQ{x,y,z}
         p = y*z + z*y - x^2
         q = x*z + z*x - y^2
         r = z^2 - x*y - y*x
         I = ncLeftIdeal{p,q,r}
         ring I
///

doc ///
   Key
      (symbol +, NCLeftIdeal, NCLeftIdeal)
   Headline
      Sum of NCLeftIdeals
   Usage
      K = I + J
   Inputs
      I : NCLeftIdeal      
      J : NCLeftIdeal      
   Outputs
      K : NCLeftIdeal
   Description
      Text
         This command sums two NCLeftIdeals.
      Example
         A = QQ{x,y,z}
         p = y*z + z*y - x^2
         q = x*z + z*x - y^2
         r = z^2 - x*y - y*x
         I = ncLeftIdeal{p,q}
         J = ncLeftIdeal r
	 I + J
///

doc ///
   Key
      NCRightIdeal
   Headline
      Type of a right ideal in a noncommutative ring
   Description
      Text
         This defines a right ideal in a noncommutative ring.  Not much
	 can be done with these objects at this point (as one can tell
	 by the dearth of operations that take an NCRightIdeal as
	 input), , but eventually it will be a 'fully featured'
	 object.
///

doc ///
   Key
      ncRightIdeal
      (ncRightIdeal, List)
      (ncRightIdeal, NCRingElement)
   Headline
      Define a right ideal in a noncommutative ring
   Usage
      I = ncRightIdeal fs
   Inputs
      fs : List
   Outputs
      I : NCRightIdeal      
   Description
      Text
         This defines a right ideal in a noncommutative ring.  Not much
	 can be done with these objects at this point (as one can tell
	 by the dearth of operations that take an NCRightIdeal as
	 input), but eventually it will be a 'fully featured'
	 object.
      Example
         A = QQ{x,y,z}
         p = y*z + z*y - x^2
         q = x*z + z*x - y^2
         r = z^2 - x*y - y*x
         I = ncRightIdeal{p,q,r}
///

doc ///
   Key
      (generators, NCRightIdeal)
   Headline
      Returns the generators of an NCRightIdeal
   Usage
      gensI = generators I
   Inputs
      I : NCRightIdeal      
   Outputs
      gensI : List
   Description
      Text
         Returns the generators of an NCRightIdeal.
      Example
         A = QQ{x,y,z}
         p = y*z + z*y - x^2
         q = x*z + z*x - y^2
         r = z^2 - x*y - y*x
         I = ncRightIdeal{p,q,r}
	 gens I
///

doc ///
   Key
      (ring, NCRightIdeal)
   Headline
      Returns the ring of an NCRightIdeal
   Usage
      A = ring I
   Inputs
      I : NCRightIdeal      
   Outputs
      A : NCRing
   Description
      Text
         Returns the ring of an NCRightIdeal.
      Example
         A = QQ{x,y,z}
         p = y*z + z*y - x^2
         q = x*z + z*x - y^2
         r = z^2 - x*y - y*x
         I = ncRightIdeal{p,q,r}
         ring I
///

doc ///
   Key
      (symbol +, NCRightIdeal, NCRightIdeal)
   Headline
      Sum of NCRightIdeals
   Usage
      K = I + J
   Inputs
      I : NCRightIdeal      
      J : NCRightIdeal      
   Outputs
      K : NCRightIdeal
   Description
      Text
         This command sums two NCRightIdeals.
      Example
         A = QQ{x,y,z}
         p = y*z + z*y - x^2
         q = x*z + z*x - y^2
         r = z^2 - x*y - y*x
         I = ncRightIdeal{p,q}
         J = ncRightIdeal r
	 I + J
///

doc ///
   Key
      NCIdeal
   Headline
      Type of a two-sided ideal in a noncommutative ring
   Description
      Text
         This defines a right ideal in a noncommutative ring.
///

doc ///
   Key
      ncIdeal
      (ncIdeal, List)
      (ncIdeal, NCRingElement)
   Headline
      Define a two-sided ideal in a noncommutative ring
   Usage
      I = ncIdeal fs
   Inputs
      fs : List
   Outputs
      I : NCIdeal      
   Description
      Text
         This defines a two-sided ideal in a noncommutative ring.
      Example
         A = QQ{x,y,z}
         p = y*z + z*y - x^2
         q = x*z + z*x - y^2
         r = z^2 - x*y - y*x
         I = ncIdeal{p,q,r}
///

doc ///
   Key
      (generators, NCIdeal)
   Headline
      Returns the generators of an NCIdeal
   Usage
      gensI = generators I
   Inputs
      I : NCIdeal      
   Outputs
      gensI : List
   Description
      Text
         Returns the generators of an NCIdeal.
      Example
         A = QQ{x,y,z}
         p = y*z + z*y - x^2
         q = x*z + z*x - y^2
         r = z^2 - x*y - y*x
         I = ncIdeal{p,q,r}
	 gens I
///

doc ///
   Key
      (ring, NCIdeal)
   Headline
      Returns the ring of an NCIdeal
   Usage
      A = ring I
   Inputs
      I : NCIdeal      
   Outputs
      A : NCRing
   Description
      Text
         Returns the ring of an NCIdeal.
      Example
         A = QQ{x,y,z}
         p = y*z + z*y - x^2
         q = x*z + z*x - y^2
         r = z^2 - x*y - y*x
         I = ncIdeal{p,q,r}
         ring I
///

doc ///
   Key
      (symbol +, NCIdeal, NCIdeal)
   Headline
      Sum of NCIdeals
   Usage
      K = I + J
   Inputs
      I : NCIdeal      
      J : NCIdeal      
   Outputs
      K : NCIdeal
   Description
      Text
         This command sums two NCIdeals.
      Example
         A = QQ{x,y,z}
         p = y*z + z*y - x^2
         q = x*z + z*x - y^2
         r = z^2 - x*y - y*x
         I = ncIdeal{p,q}
         J = ncIdeal r
	 I + J
///

doc ///
   Key
      (basis, ZZ, NCIdeal)
   Headline
      Returns a basis of an NCIdeal in a particular degree.
   Usage
      bas = basis(n,I)
   Inputs
      n : ZZ
      I : NCIdeal
   Outputs
      bas : NCMatrix
   Description
      Text
         This command returns a basis (or minimal generating set, if
	 the ground ring is not a field), of a homogeneous two-sided
	 ideal in a noncommutative ring.
      Example
         A = QQ{x,y,z}
         p = y*z + z*y - x^2
         q = x*z + z*x - y^2
         r = z^2 - x*y - y*x
         I = ncIdeal{p,q,r}
	 bas = basis(3,I)
///

doc ///
   Key
      (basis, ZZ, NCLeftIdeal)
   Headline
      Returns a basis of an NCLeftIdeal in a particular degree.
   Usage
      bas = basis(n,I)
   Inputs
      n : ZZ
      I : NCLeftIdeal
   Outputs
      bas : NCMatrix
   Description
      Text
         This command returns a basis (or minimal generating set, if
	 the ground ring is not a field), of a homogeneous left
	 ideal in a noncommutative ring.
      Example
         A = QQ{x,y,z}
         p = y*z + z*y - x^2
         q = x*z + z*x - y^2
         r = z^2 - x*y - y*x
         I = ncLeftIdeal{p,q,r}
	 bas = basis(3,I)
///

doc ///
   Key
      (basis, ZZ, NCRightIdeal)
   Headline
      Returns a basis of an NCRightIdeal in a particular degree.
   Usage
      bas = basis(n,I)
   Inputs
      n : ZZ
      I : NCIdeal
   Outputs
      bas : NCMatrix
   Description
      Text
         This command returns a basis (or minimal generating set, if
	 the ground ring is not a field), of a homogeneous right
	 ideal in a noncommutative ring.
      Example
         A = QQ{x,y,z}
         p = y*z + z*y - x^2
         q = x*z + z*x - y^2
         r = z^2 - x*y - y*x
         I = ncRightIdeal{p,q,r}
	 bas = basis(3,I)
///

doc ///
   Key
      (basis, ZZ, NCRing)
   Headline
      Returns a basis of an NCRing in a particular degree.
   Usage
      bas = basis(n,B)
   Inputs
      n : ZZ
      I : NCRing
   Outputs
      bas : NCMatrix
   Description
      Text
         This command returns a basis (or minimal generating set, if
	 the ground ring is not a field), of a graded noncommutative
         ring.
      Example
         A = QQ{x,y,z}
         p = y*z + z*y - x^2
         q = x*z + z*x - y^2
         r = z^2 - x*y - y*x
         I = ncIdeal{p,q,r}
	 B = A/I
	 bas = basis(4,B)
///

doc ///
   Key
     setWeights
     (setWeights,NCRing,List)
   Headline
      Set a nonstandard grading for a NCRing.
   --Usage
   --Inputs
   --Outputs
   Description
      Example
      -- need to finish unit tests
      Text
        stuff
///

doc ///
   Key
      (isHomogeneous, NCIdeal)
      (isHomogeneous, NCRightIdeal)
      (isHomogeneous, NCLeftIdeal)
      (isHomogeneous, NCRing)
      (isHomogeneous, NCMatrix)
      (isHomogeneous, NCRingElement)
   Headline
      Determines whether the input defines a homogeneous object
   --Usage
   --Inputs
   --Outputs
   Description
      Example
      -- need to finish unit tests
      Text
        stuff
///


doc ///
   Key
      rightKernelBergman
      (rightKernelBergman,NCMatrix)
      assignDegrees
      (assignDegrees,NCMatrix)
      (assignDegrees,NCMatrix,List,List)
      [rightKernelBergman,DegreeLimit]
   Headline
      Methods for computing kernels of matrices over noncommutative rings using Bergman
   --Usage
   --Inputs
   --Outputs
   Description
      Example
         A = QQ{x,y,z}
         f1 = y*z + z*y - x^2
         f2 = x*z + z*x - y^2
         f3 = z^2 - x*y - y*x
         g = -y^3-x*y*z+y*x*z+x^3
         I = ncIdeal {f1,f2,f3,g}
         B = A/I
         M3 = ncMatrix {{x,y,z,0}, {-y*z-2*x^2,-y*x,z*x-x*z,x},{x*y-2*y*x,x*z,-x^2,y}, {-y^2-z*x,x^2,-x*y,z}}
         assignDegrees(M3,{1,0,0,0},{2,2,2,1})
         ker1M3 = rightKernelBergman(M3)
         M3*ker1M3 == 0
         ker2M3 = rightKernelBergman(ker1M3)
         ker1M3*ker2M3 == 0
         ker3M3 = rightKernelBergman(ker2M3)
         ker2M3*ker3M3 == 0
      Text
         stuff
///

doc ///
   Key
      isLeftRegular
      (isLeftRegular,NCRingElement,ZZ)
      isRightRegular
      (isRightRegular,NCRingElement,ZZ)
   Headline
      Determines if a given (homogeneous) element is regular in a given degree
   --Usage
   --Inputs
   --Outputs
   Description
      Text
         stuff
///

doc ///
   Key
      isCentral
      (isCentral,NCRingElement)
      (isCentral,NCRingElement,NCGroebnerBasis)
      centralElements
      (centralElements, NCRing, ZZ)
   Headline
      Methods for finding/checking central elements
   --Usage
   --Inputs
   --Outputs
   Description
      Example
        A = QQ{x,y,z}
        I = ncIdeal { y*z + z*y - x^2,x*z + z*x - y^2,z^2 - x*y - y*x}
        B = A/I
        g = -y^3-x*y*z+y*x*z+x^3
        h = x^2 + y^2 + z^2
        isCentral h
        isCentral g
        centralElements(B,2)
        centralElements(B,3)

      Text
         We have not yet implemented the check in a fixed degree.
///




doc ///
   Key
      normalElements
 --     (normalElements, NCRingMap, ZZ) -- does this key exist?
      (normalElements, NCQuotientRing, ZZ, Symbol, Symbol)
      (normalElements, NCRingMap, ZZ)
      (isNormal, NCRingElement)
   Headline
      Computes normal monomials and components of the variety of normal elements in a given degree
   --Usage
   --Inputs
   --Outputs
   Description
      Text
         stuff
///


doc ///
   Key
      normalAutomorphism
      (normalAutomorphism,NCRingElement)
   Headline
      Computes the automorphism determined by a normal homogeneous element
   --Usage
   --Inputs
   --Outputs
   Description
      Text
        This is the type of a matrix over a noncommutative graded ring.
///

doc ///
   Key
      leftMultiplicationMap
      (leftMultiplicationMap,NCRingElement,ZZ)
      (leftMultiplicationMap,NCRingElement,List,List)
      rightMultiplicationMap
      (rightMultiplicationMap,NCRingElement,ZZ)
      (rightMultiplicationMap,NCRingElement,List,List)
   Headline
      Computes a matrix for left or right multiplication by a homogeneous element
   --Usage
   --Inputs
   --Outputs
   Description
      Text
         stuff
///

doc ///
   Key
      rightKernel
      (rightKernel,NCMatrix,ZZ)
      NumberOfBins
      [rightKernel,NumberOfBins]
      [rightKernel,Verbosity]
   Headline
      Method for computing kernels of matrices over noncommutative rings in a given degree without using Bergman
   --Usage
   --Inputs
   --Outputs
   Description
      Text
         stuff
///


doc ///
   Key
      quadraticClosure
      (quadraticClosure,NCIdeal)
      (quadraticClosure,NCQuotientRing)
   Headline
      Creates the subideal generated by quadratic elements of a given ideal
   --Usage
   --Inputs
   --Outputs
   Description
      Text
         stuff
///

doc ///
   Key
      homogDual
      (homogDual,NCIdeal)
      (homogDual,NCQuotientRing)
   Headline
      Computes the dual of a pure homogeneous ideal
   --Usage
   --Inputs
   --Outputs
   Description
      Text
         stuff
///

doc ///
   Key
      sparseCoeffs
      (sparseCoeffs,List)
      (sparseCoeffs,NCRingElement)
      [sparseCoeffs,Monomials]
   Headline
      Converts ring elements into vectors over the coefficient ring
   --Usage
   --Inputs
   --Outputs
   Description
      Example
         A=QQ{a, b, c, d, e, f, g, h}
	 F = a^2+b^2+c^2+d^2+e^2+f^2+g^2+h^2;
	 bas = flatten entries basis(2,A);
	 #bas
	 sparseCoeffs(F,Monomials=>bas)
	 sparseCoeffs(toList (10:F),Monomials=>bas)
      Text
         stuff
///


doc ///
   Key
      NCRingMap
      ncMap
      (ncMap,NCRing,NCRing,List)
      (ncMap,Ring,NCRing,List)
      (ncMap,NCRing,Ring,List)
      (ambient, NCRingMap)
      (isHomogeneous, NCRingMap)
      (isWellDefined, NCRingMap)
      (symbol /, List, NCRingMap)
      (matrix, NCRingMap)
      (symbol @@, NCRingMap, NCRingMap)
      (symbol SPACE, NCRingMap, NCRingElement)
      (symbol SPACE, NCRingMap, NCMatrix)
      (symbol _, NCRingMap, ZZ)
      (source, NCRingMap)
      (target, NCRingMap)
   Headline
      Creates a map from a non-commutative ring
   --Usage
   --Inputs
   --Outputs
   Description
      Text
         stuff
///

doc ///
   Key
      oreExtension
      (oreExtension,NCRing,NCRingMap,NCRingMap,NCRingElement)
      (oreExtension,NCRing,NCRingMap,NCRingMap,Symbol)
      (oreExtension,NCRing,NCRingMap,NCRingElement)
      (oreExtension,NCRing,NCRingMap,Symbol)
   Headline
      Creates an Ore extension of a noncommutative ring
   --Usage
   --Inputs
   --Outputs
   Description
      Example
         B = skewPolynomialRing(QQ,(-1)_QQ,{x,y,z,w})
	 sigma = ncMap(B,B,{y,z,w,x})
	 C = oreExtension(B,sigma,a)
      Text
         stuff
///

doc ///
   Key
      oreIdeal
      (oreIdeal,NCRing,NCRingMap,NCRingMap,NCRingElement)
      (oreIdeal,NCRing,NCRingMap,NCRingMap,Symbol)
      (oreIdeal,NCRing,NCRingMap,NCRingElement)
      (oreIdeal,NCRing,NCRingMap,Symbol)
   Headline
      Creates the defining ideal of an Ore extension of a noncommutative ring
   --Usage
   --Inputs
   --Outputs
   Description
      Example
         B = skewPolynomialRing(QQ,(-1)_QQ,{x,y,z,w})
	 sigma = ncMap(B,B,{y,z,w,x})
	 C = oreIdeal(B,sigma,a)
      Text
         stuff
///

doc ///
   Key
      endomorphismRing
      (endomorphismRing,Module,Symbol)
      endomorphismRingGens
       minimizeRelations
      (minimizeRelations,List)
   Headline
      Methods for creating endomorphism rings of modules over a commutative ring
   --Usage
   --Inputs
   --Outputs
   Description
      Example
         Q = QQ[a,b,c]
         R = Q/ideal{a*b-c^2}
         kRes = res(coker vars R, LengthLimit=>7);
         M = coker kRes.dd_5
         B = endomorphismRing(M,X);
         gensI = gens ideal B;
         gensIMin = minimizeRelations(gensI);

         Q = QQ[a,b,c,d]
         R = Q/ideal{a*b+c*d}
         kRes = res(coker vars R, LengthLimit=>7);
         M = coker kRes.dd_5
         B = endomorphismRing(M,Y);
         gensI = gens ideal B;
         gensIMin = minimizeRelations(gensI);
      Text
         stuff
///



doc ///
   Key
      skewPolynomialRing
      (skewPolynomialRing,Ring,Matrix,List)
   Headline
      Defines a skew polynomial ring via a skewing matrix
   --Usage
   --Inputs
   --Outputs
   Description
      Example
         R = QQ[q]/ideal{q^4+q^3+q^2+q+1}
         B = skewPolynomialRing(R,q,{x,y,z,w})
         x*y == q*y*x
         C = skewPolynomialRing(QQ,1_QQ, {x,y,z,w})
         isCommutative C
         isCommutative B
         --abC = abelianization C
         --abC' = abelianization ambient C
         abC = toM2Ring C
	 abC' = toM2Ring ambient C
	 ideal abC == 0
         ideal abC' == 0
         Bop = oppositeRing B
         y*x == q*x*y
      Text
         Link to oppositeRing.
///


doc ///
   Key
      skewPolynomialRing
      (skewPolynomialRing,Ring,QQ,List)
      (skewPolynomialRing,Ring,ZZ,List)
      (skewPolynomialRing,Ring,RingElement,List)
   Headline
      Methods for working defining a skew polynomial ring where the skewing factor
      is the same for each pair of elements.
   --Usage
   --Inputs
   --Outputs
   Description
      Example
         R = QQ[q]/ideal{q^4+q^3+q^2+q+1}
         B = skewPolynomialRing(R,q,{x,y,z,w})
         x*y == q*y*x
         C = skewPolynomialRing(QQ,1_QQ, {x,y,z,w})
         isCommutative C
         isCommutative B
         --abC = abelianization C
         --abC' = abelianization ambient C
         abC = toM2Ring C
	 abC' = toM2Ring ambient C
	 ideal abC == 0
         ideal abC' == 0
         Bop = oppositeRing B
         y*x == q*x*y
      Text
         Link to oppositeRing.
///

doc ///
   Key
      toM2Ring
      (toM2Ring,NCRing)
   Headline
      Compute the abelianization of an NCRing and returns a PolynomialRing. 
   Usage
   --Inputs
   --Outputs
   Description
      Example

      Text
         Link to toNCRing.
///

doc ///
   Key
      toNCRing
      (toNCRing,Ring)
   Headline
      Converts a PolynomialRing to an NCRing
   Usage
      Anc = toNCRing A 
   Inputs
      A : PolynomialRing
   Outputs
      Anc : NCRing
   Description
      Text
         This function converts commutative rings and quotients of 
	 exterior algebras created by Macaulay2 to the type of an NCRing. An
	 error is returned if the input ring has some commutative and some
	 skew-commutative generators.
      Example
         
      Text
         Link to toM2Ring.
///


doc ///
   Key
      oppositeRing
      (oppositeRing,NCRing)
   Headline
      Creates the opposite ring of a noncommutative ring
   Usage                    
      Aop = oppositeRing A  
   Inputs
      A : NCRing	
   Outputs        
      Aop : NCRing 
   Description
      Text 
         Given an NCRing A, this creates an NCRing whose defining NCIdeal is generated by 
	 the "opposites" - elements whose noncommutative monomial terms have been reversed - 
	 of the generators of the defining NCIdeal of A. If the coefficient ring of A is a
	 bergman ring, an NCGroebnerBasis is computed for Aop.
      Text 
         Link to skewPolynomialRing  
      Example
          R = QQ[q]/ideal{q^4+q^3+q^2+q+1}
          A = skewPolynomialRing(R,q,{x,y,z,w}) 
	  x*y == q*y*x
          Aop = oppositeRing A
	  y*x == q*x*y 		
///

doc ///
   Key
      normalFormBergman
      (normalFormBergman,List,NCGroebnerBasis)
      (normalFormBergman,NCRingElement,NCGroebnerBasis)
      [normalFormBergman,CacheBergmanGB]
      [normalFormBergman,MakeMonic]
      [normalFormBergman,NumModuleVars]
      [normalFormBergman,DegreeLimit]
   Headline
      Calls Bergman for a normal form calculation
   --Usage
   --Inputs
   --Outputs
   Description
      Text
         stuff
///


doc ///
   Key
      isReduced
   Headline
      Determines if a given element is in normal form with respect to a Groebner basis
   Usage
      b = w.isReduced
   Inputs
      w : NCRingElement
   Outputs
      b : Boolean
   Description
      Text
         This flag is set ot true if a given element of an NCRing (considered 
	 as an element of an appropriate NCPolynomialRing) is in normal form 
         relative to a Groebner basis.
      Example
         A = QQ{x,y}
	 I = ncIdeal{y^2-x^2}
	 Igb = ncGroebnerBasis(I)
	 w = y^3+x*y^2
	 w.isReduced
	 w'=w % Igb
	 w'.isReduced
///


doc ///
   Key
      hilbertBergman
      (hilbertBergman, NCQuotientRing)
      [hilbertBergman,DegreeLimit]
      DegreeVariable
   Headline
      Calls Bergman to compute the Hilbert series of an NCQuotientRing
   Usage
     hseries = hilbertBergman(A,DegreeLimit=>d)
   Inputs
     A : NCQuotientRing
     d : ZZ
   Outputs
     hseries : ZZ[T]
   Description
      Text
         This method calls the bergman function ncpbhgroebner to compute the Hilbert
	 series of an NCQuotientRing. The input ring must be a ring over QQ or ZZ/p.
	 At this time, the output is correct only for NCRings with a standard grading -
	 all generators have degree 1. The output is returned as a polynomial in ZZ[T].
      Example
///

doc ///
   Key
      "Basic operations on noncommutative algebras"
   Description
      Example
         A = QQ{x,y,z}
	 f = y*z + z*y - x^2
	 g = x*z + z*x - y^2
	 h = z^2 - x*y - y*x
	 f*g
	 f^2
	 f-g 
         3*g
         f+g
	 B = A/ncIdeal{f,g,h}
	 j = -y^3-x*y*z+y*x*z+x^3
	 k = x^2 + y^2 + z^2
	 j*k
	 k^3
      Text
         Here will go an extended example
///

doc ///
   Key
      "Using the Bergman interface"
   Description
      Text
         Here will go an extended example
///
