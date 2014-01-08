--finds a polynomial function fitting a sequence through Lagrange interpolation
polynomialInterpolation = method ()
polynomialInterpolation (List,PolynomialRing) := RingElement => (L,R) -> (
	if #gens R != 1 then error "must be a single variable polynomial ring";
	n := #L;
	p := for j from 0 to n-1 list (
		L_j * product(0..(j-1) | (j+1) .. (n-1), k -> (t - (k+1)) / ((j+1) - (k+1)))
		);
	sum(p)
)
polynomialInterpolation (List) := RingElements => L -> polynomialInterpolation(L, QQ[t])

--finds a rational function fitting the sequence through Chebyshev interpolation
rationalInterpolation = method ()
rationalInterpolation (List, PolynomialRing) := RingElement => (L,R) -> (
	if #gens R != 1 then error "must be a single variable polynomial ring";
	QQ[t];
	n := #L;
	(sum(1..(n-2), j -> ((-1)^j * L_j)/(t-(j+1))) + 1/2 * (L_0)/(t-1) + 1/2 * ((-1)^(n-1) * L_(n-1))/(t-n)) / 
		(sum(1..(n-2), j -> ((-1)^j)/(t-(j+1))) + 1/2 * ((-1)^0)/(t-1) + 1/2 * ((-1)^(n-1))/(t-n))
)
rationalInterpolation (List) := RingElement => L -> rationalInterpolation(L, QQ[t])

--finds a rational function that fits the sequence via polynomialInterpolation of the numerator and the denominator
guessRational = method ()
guessRational (List,PolynomialRing) := RingElement => (L,R) -> (
	n := #L;
	a := apply(L, i -> numerator i);
	b := apply(L, i -> denominator i);
	p = polynomialInterpolation(a,R);
	q = polynomialInterpolation(b,R);
	p/q
)
guessRational (List) := RingElement => L -> guessRational(L, QQ[t])
