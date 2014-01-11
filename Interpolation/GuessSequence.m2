newPackage(
             "GuessSequence",
             Version => "0.1", 
             Date => "10 January 2014",
             Authors => {{Name => "David Cook II", Email => "dcook8@nd.edu", HomePage => "http://www.nd.edu/~dcook8/"},
			{Name => "Caroline Jansen", Email => "cjansen@nd.edu"}},
             Headline => "Package for guessing sequence functions through interpolation and other methods",
             DebuggingMode => false
             )

export {
	"polynomialInterpolation",
	"cubicSplines",
	"evaluateCubicSpline",
	"rationalInterpolation",
	"floaterHormann",
	"guessRational",
	"guessLinearRecurrence",
	"applyLinearRecurrence"
}

---------------------------
--polynomial tools---------
---------------------------

--finds a polynomial function fitting a sequence through Lagrange interpolation
polynomialInterpolation = method ()
polynomialInterpolation (List,PolynomialRing) := RingElement => (L,R) -> (
	if #gens R != 1 then error "must be a single variable polynomial ring";
	n := #L;
	t := first gens R;
	p := for j from 0 to n-1 list (
		L_j * product(0..(j-1) | (j+1) .. (n-1), k -> (t - (k+1)) / ((j+1) - (k+1)))
		);
	sum(p)
)
polynomialInterpolation (List) := RingElement => L -> (
	t := symbol t;
	polynomialInterpolation(L, QQ[t])
)

cubicSplines = method()
cubicSplines (List,PolynomialRing) := List => (L,R) -> (
	cubicCurve := (p0, p1, m0, m1, t) -> (2*t^3 - 3*t^2 + 1) * p0 + (t^3 - 2*t^2 + t) * m0 + (-2*t^3 + 3*t^2) * p1 + (t^3 - t^2) * m1;
	t := first gens R;
	n := #L;
	firstTerm := cubicCurve(L_0, L_1, (L_1 - L_0)/2, (L_2 - L_0)/2, t);
	P := for i from 1 to n-3 list cubicCurve(L_i, L_(i+1), (L_(i+1) - L_(i-1))/2, (L_(i+2) - L_i)/2, t);
	lastTerm := cubicCurve(L_(n-2), L_(n-1), (L_(n-1) - L_(n-3))/2, (L_(n-1) - L_(n-2))/2, t);
	{firstTerm} | P | {lastTerm}
)
cubicSplines List := RingElement => L -> (
	t := symbol t;
	cubicSplines(L,QQ[t])
)

evaluateCubicSpline = method()
evaluateCubicSpline (List, QQ) := QQ => (P, v) -> (
	if v < 1 or v > #P+1 then error "The value v must be in the interval [1, n].";
	t := first gens ring P_0;
	i := floor v;
	j := v - i;
	if i == #P+1 then sub(P_(-1), t => 1) else sub(P_(i-1), t => j)
)
evaluateCubicSpline (List, ZZ) := QQ => (P, v) -> evaluateCubicSpline(P, sub(v, QQ))
evaluateCubicSpline (List, RR) := QQ => (P, v) -> evaluateCubicSpline(P, lift(v, QQ))

-------------------------
-----rational tools------
-------------------------

--finds a rational function fitting the sequence through Chebyshev interpolation
rationalInterpolation = method ()
rationalInterpolation (List, PolynomialRing) := RingElement => (L,R) -> (
	if #gens R != 1 then error "must be a single variable polynomial ring";
	t := first gens R;
	n := #L;
	(sum(1..(n-2), j -> ((-1)^j * L_j)/(t-(j+1))) + 1/2 * (L_0)/(t-1) + 1/2 * ((-1)^(n-1) * L_(n-1))/(t-n)) / 
		(sum(1..(n-2), j -> ((-1)^j)/(t-(j+1))) + 1/2 * ((-1)^0)/(t-1) + 1/2 * ((-1)^(n-1))/(t-n))
)
rationalInterpolation (List) := RingElement => L -> (
	t := symbol t;
	rationalInterpolation(L, QQ[t])
)

--interpolates a rational function using the floater hormann algorithm
--user can choose a d st 0 <= d <= n-1 (where n is the size of list
--to compute a family of interpolations
floaterHormann = method ()
floaterHormann (List,PolynomialRing,ZZ) := RingElement => (L,R,d) -> (
	n := #L;
	t := first gens R;
	if d > n-1 or d < 0 then error "d must be <= n-1 and >= 0";
	p := for i from 0 to n-1-d list (
		polynomialInterpolation(drop(L,i),R)
	);
	lambda := for i from 0 to n-1-d list (
		(-1)^i/product(toList (0..d), j -> t - L_(j+i))
	);
	sum(toList(0..(n-1-d)), j -> product({lambda_j , p_j}))/sum(toList(0..(n-1-d)), j -> lambda_j)
)

floaterHormann(List,ZZ) := RingElement => (L,d) -> (
	t := symbol t;
	floaterHormann(L,QQ[t],d)
)

--finds a rational function that fits the sequence via polynomialInterpolation of the numerator and the denominator
guessRational = method ()
guessRational (List,PolynomialRing) := RingElement => (L,R) -> (
	n := #L;
	L = apply(L, i -> i/1);
	a := apply(L, i -> numerator i);
	b := apply(L, i -> denominator i);
	p := polynomialInterpolation(a,R);
	q := polynomialInterpolation(b,R);
	p/q
)
guessRational (List) := RingElement => L -> (
	t := symbol t;
	guessRational(L, QQ[t])
)

-------------------------
---linear recurrence-----
-------------------------

--guesses linear recurrence for a sequence
guessLinearRecurrence = method()
guessLinearRecurrence List := List => L -> (
    -- Prepare the matrix of terms
    m := #L//2;
    M := matrix apply(m, i -> L_(toList(i..i+m-1)));
    -- Determine the order of the sequence
    r := rank M;
    N := M_{0..r-1}^{0..r-1};
    -- M2 inverts matrices in the ring they live in, so we may need to change it.
    if all(L, l -> instance(l, ZZ)) then N = sub(N, QQ);
    -- Solve for the coefficients
    first entries(matrix {L_(toList(r..2*r-1))} * inverse N)
    )

applyLinearRecurrence = method()
applyLinearRecurrence (List, List, ZZ) := List => (C, L, n) -> (
    if #C > #L then error "Not enough initial data is known.";
    for i to n-1 do L = append(L, sum(#C, i -> C_i * L_(#L-#C+i)));
    L
    )

-----------------------------
-------documentation---------
-----------------------------

beginDocumentation()

--polynomialInterpolation
doc ///
	Key
		polynomialInterpolation
		(polynomialInterpolation, List, PolynomialRing)
		(polynomialInterpolation, List)
	Headline
		Interpolates a list of points using Lagrange interpolations
	Usage
		P = polynomialInterpolation(L,R)
		P = polynomialInterpolation(L)
	Inputs
		L:List
		R:PolynomialRing
	Outputs
		P:RingElement
			The unique polynomial interpolated from the list of points supplied, in the ring specified.
	Description
		Text
			Lagrange interpolation finds a unique polynomial describing a given sequence of points, e.g., the interpolation of the list {1,4,9,16} would be t^2.  The polynomial has degree less than or equal to the length of the list, and is unique.  The list the user inputs is assumed to be indexed from 1.  The method will produce a polynomial in the polynomial ring specified by the user; if the user does not specify a ring, the method assumes QQ[t].
		Example
			L = {1,4,9,16}
			P = polynomialInterpolation(L, QQ[x])
			L = {2, 10, 30, 64, 130}
			P = polynomialInterpolation L
	SeeAlso
		cubicSplines
///

     --cubicSplines
doc ///
	Key
		cubicSplines
		(cubicSplines, List, PolynomialRing)
		(cubicSplines, List)
	Headline
		Calculates cubic splines based on a sequence
	Usage
		P = cubicSplines(L,R)
		P = cubicSplines(L)
	Inputs
		L:List
		R:PolynomialRing
	Outputs
		P:List
			The list of cubic functions over intervals of length 1 approximating the function described by the sequence.
	Description
		Text
			The user inputs their desired sequence as a list indexed from 1, e.g. a_1, a_2, ... , a_n.  The method will return a list of polynomials defined on the interval [0,1], i.e. the Cubic Hermite Splines.  If the user wants to find a_i using P, the cubic splines calculated by the method, the user enters: sub(P_(i-1), t => 0) or sub(P_(i-2), t => 1) where t is the polynomial ring generator.  To estimate a value in between i and i + 1, say i + j with j in (0,1), then the user would enter:  sub(P_(i-1), t => j).  This estimation is best for values in between 1 and n, where n is the length of the sequence you record.  The polynomials are in the ring specified by the user; if there is no ring specified, the default is QQ[t].
		Example
			L = {1,4,9,16}
			P = cubicSplines L
			R = QQ[x]
			P = cubicSplines(L,R)
	SeeAlso
		polynomialInterpolation
		evaluateCubicSpline
///

--evaluateCubicSpline
doc ///
	Key
		evaluateCubicSpline
		(evaluateCubicSpline, List, QQ)
		(evaluateCubicSpline, List, ZZ)
		(evaluateCubicSpline, List, RR)
	Headline
		Evaluates a cubic spline at a specified number
	Usage
		E = evaluateCubicSpline(S,q)
		E = evaluateCubicSpline(S,z)
		E = evaluateCubicSpline(S,r)
	Inputs
		S:List
		q:QQ
		z:ZZ
		r:RR
	Outputs
		E:QQ
			The rational number approximating the spline at the supplied point.
	Description
		Text
			The user should input a list of polynomials obtained from generating cubic splines using the cubicSplines method, and a number they wish the spline to be evaluated at, which must be between 1 and the length of their sequence.  The method will evaluate the spline at the appropriate interval using the appropriate polynomial, and return a rational number.
		Example
			L = {1,4,9,16}
			R = QQ[x]
			P = cubicSplines(L,R)
			evaluateCubicSpline(P,2)
			evaluateCubicSpline(P,1.5)
			evaluateCubicSpline(P, sqrt(2))
	SeeAlso
		polynomialInterpolation
		cubicSplines
///

--rationalInterpolation
doc ///
	Key
		rationalInterpolation
		(rationalInterpolation, List, PolynomialRing)
		(rationalInterpolation, List)
	Headline
		Interpolates a list of points using Chebyshev's algorithm for rational interpolation
	Usage
		P = rationalInterpolation(L,R)
		P = rationalInterpolation(L)
	Inputs
		L:List
		R:PolynomialRing
	Outputs
		P:RingElement
			The rational function interpolated from the list of points supplied, in the fraction field of the ring specified.
	Description
		Text
			This method uses rational interpolation to guess a rational function that matches the specified points.  It will return a rational expression in the fraction field of the ring specified; if no ring is specified, the fraction returned will be in the fraction field of QQ[t].  The method assumes the user is indexing 1 to n.
		Example
			L = {1,1/4,1/9}
			P = rationalInterpolation(L, QQ[x])
			L = {1,2,3,4,5}
			P = polynomialInterpolation L
	SeeAlso
			floaterHormann
			guessRational
///

--floaterHormann
doc ///
	Key
		floaterHormann
		(floaterHormann, List, PolynomialRing, ZZ)
		(floaterHormann, List, ZZ)
	Headline
		Interpolates a list of points using the Floater-Hormann algorithm for interpolation.
	Usage
		P = floaterHormann(L,R,d)
		P = floaterHormann(L,d)
	Inputs
		L:List
		R:PolynomialRing
		d:ZZ
	Outputs
		P:RingElement
			The rational function interpolated from the list of points supplied, in the fraction field of the ring specified.
	Description
		Text
			This method uses the Floater-Hormann algorithm for rational interpolation to guess a rational function that matches the specified points.  It will return a rational expression in the fraction field of the ring specified; if no ring is specified, the fraction returned will be in the fraction field of QQ[t].  The method assumes the user is indexing 1 to n.  The value d can be greater than or equal to 0 and less than or equal to n-1.  By varying d, the user can obtain a family of interpolated rational functions matching the sequence given.
		Example
			L = {1,1/4,1/9,1/16,1/25}
			P = floaterHormann(L, QQ[x], 3)
			L = {1,2,3,4,5}
			P = floaterHormann (L,4)
	SeeAlso
			rationalInterpolation
			guessRational
///

--guessRational
doc ///
	Key
		guessRational
		(guessRational, List, PolynomialRing)
		(guessRational, List)
	Headline
		Interpolates a list of points as a rational function by using Lagrange interpolation on the numerator and the denominator.
	Usage
		P = guessRational(L,R)
		P = guessRational(L)
	Inputs
		L:List
		R:PolynomialRing
	Outputs
		P:RingElement
			The rational function interpolated from the list of points supplied, in the fraction field of the ring specified.
	Description
		Text
			This method uses the Lagrange interpolation on the sequence of the numerators of a rational sequence and also interpolates the sequence of denominators of the rational sequence.  It then takes their fraction, in the fraction field of the ring specified; if no ring is specified by the user, the method will return a fraction in the fraction field of QQ[t].   The method assumes the user is indexing 1 to n.
		Example
			L = {1,1/4,1/9,1/16,1/25}
			P = guessRational(L, QQ[x])
			L = {1,2,3,4,5}
			P = guessRational (L)
	SeeAlso
			rationalInterpolation
			floaterHormann
///

--guessLinearRecurrence
doc ///
	Key
		guessLinearRecurrence
		(guessLinearRecurrence, List)
	Headline
		Guesses a rule for a sequence based on linear recurrence.
	Usage
		P = guessLinearRecurrence(L)
	Inputs
		L:List
	Outputs
		P:List
	Description
		Text
			Guesses a homogenous linear recurrence with constant coefficients.  The user inputs a sequence a_1, ... , a_n, and the method finds a linear recurrence defined as:  a_t = c_1*a_(t-k) + ... + c_k*a_(t-1).  The output is the list of coefficients c_1, ..., c_k.
		Example
			-- a_n = a_(n-3) + a_(n-1)
			guessLinearRecurrence {1, 2, 3, 4, 6, 9, 13, 19, 28, 41, 60, 88, 129}
			-- Fibonacci Sequence
			guessLinearRecurrence {1, 1, 2, 3, 5, 8, 13, 21}
			-- Perrin Sequence
			guessLinearRecurrence {3, 0, 2, 3, 2, 5, 5, 7, 10, 12, 17}
			-- Pell numbers
			guessLinearRecurrence {0, 1, 2, 5, 12, 29, 70, 169, 408, 985}
			-- Padovan Sequence
			guessLinearRecurrence {1, 0, 0, 1, 0, 1, 1, 1, 2, 2, 3, 4, 5}
	SeeAlso
			applyLinearRecurrence
///

--applyLinearRecurrence
doc ///
	Key
		applyLinearRecurrence
		(applyLinearRecurrence, List, List, ZZ)
	Headline
		Calculates further terms in a sequence based on linear recurrence.
	Usage
		P = applyLinearRecurrence(C,L,n)
	Inputs
		C:List
			list of coefficients used in recurrence
		L:List
			beginning of sequence to be extended
		n:ZZ
			number of further terms of the sequence to calculate
	Outputs
		P:List
			the extended sequence
	Description
		Text
			The user enters a part of the sequence {a_n} as the second list in the constructor and coefficients as the first list in the constructor; the method will calculate the next n terms in the sequence using the formula for linear recurrence: a_t = c_1*a_(t-k) + ... + c_k*a_(t-1), where c_1, ... , c_k are the coefficients (indexed from 1).  The output will be a list consisting of the sequence inputted appended to the next n terms of the sequence calculated by the method.
		Example
	SeeAlso
			guessLinearRecurrence
///

-----------
--TESTS----
-----------

--polynomialInterpolation
TEST ///

	L = {1,4,9,16,25,36};
	P = polynomialInterpolation L;
	assert(first degree P == 2)

	L = {2, 10, 30, 68};
	R = QQ[y];
	P = polynomialInterpolation(L, R);
	assert(P == (y^3 + y)_R)

///

--cubicSplines
TEST ///

	L = {1,4,9};
	P = cubicSplines L;
	assert(sub(P_(0), first(gens(ring P_(0))) => 0) == 1 and sub(P_(1), first(gens(ring P_(1))) => 0) == 4 and sub(P_(2), first(gens(ring P_(2))) => 0) == 9)
	R = QQ[x];
	assert(sub(P_(0), first(gens(ring P_(0))) => 0) == 1 and sub(P_(1), first(gens(ring P_(1))) => 0) == 4 and sub(P_(2), first(gens(ring P_(2))) => 0) == 9)

///

--rationalInterpolation
TEST ///

	L = {1/3,1/9,1/27};
	P = rationalInterpolation (L);
	assert(sub(P, first(gens(ring P)) => 1) == 1/3 and sub(P, first(gens(ring P)) => 2) == 1/9 and sub(P, first(gens(ring P)) => 3) == 1/27)

	L = {1/3,1/9,1/27};
	R = QQ[x];
	P = rationalInterpolation (L,R);
	assert(sub(P, first(gens(ring P)) => 1) == 1/3 and sub(P, first(gens(ring P)) => 2) == 1/9 and sub(P, first(gens(ring P)) => 3) == 1/27) 
	

///

--floaterHormann
TEST ///

	L = {1/3,1/9,1/27};
	P = floaterHormann (L,1);
	assert(sub(P, first(gens(ring P)) => 1) == 1/3 and sub(P, first(gens(ring P)) => 2) == 1/9 and sub(P, first(gens(ring P)) => 3) == 1/27)

	L = {1/3,1/9,1/27};
	R = QQ[x];
	P = floaterHorman (L,R);
	assert(sub(P, first(gens(ring P)) => 1) == 1/3 and sub(P, first(gens(ring P)) => 2) == 1/9 and sub(P, first(gens(ring P)) => 3) == 1/27)

///

--guessRational
TEST ///

	L = {1/1, 1/2, 1/3, 1/4};
	P = guessRational L;
	assert(first degree P == -1)

	L = {1/1, 1/2, 1/3, 1/4};
	R = QQ[x];
	P = guessRational (L,R);
	assert(P == 1/(x_R))	

///

guessLinearRecurrence
TEST ///

	L = {1,1,2,3,5,8};
	R = guessLinearRecurrence L;
	assert(R == {1,1})

///

applyLinearRecurrence
TEST ///

	L = {1,1,2,3,5};
	R = applyLinearRecurrence({1,1},L, 3);
	assert(R == {1,1,2,3,5,8,13,21})

///

end;
