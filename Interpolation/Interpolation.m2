--finds a polynomial function fitting a sequence through Lagrange interpolation
polynomialInterpolation = method ()
polynomialInterpolation (List,PolynomialRing) := RingElement => (L,R) -> (
	if #gens R != 1 then error "must be a single variable polynomial ring";
	n := #L;
	y := for k from 0 to n-1 list L_k;
	x := toList (1..n);
	p := for j from 0 to n-1 list (
		y_j * product(0..(j-1) | (j+1) .. (n-1), k -> (t - (k+1)) / ((j+1) - (k+1)))
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
	y := for k from 0 to n-1 list L_k;
	x := toList (1..n);
	(sum(1..(n-2), j -> ((-1)^j * y_j)/(t-(j+1))) + 1/2 * (y_0)/(t-1) + 1/2 * ((-1)^(n-1) * y_(n-1))/(t-n)) / 
		(sum(1..(n-2), j -> ((-1)^j)/(t-(j+1))) + 1/2 * ((-1)^0)/(t-1) + 1/2 * ((-1)^(n-1))/(t-n))
)
rationalInterpolation (List) := RingElement => L -> rationalInterpolation(L, QQ[t])

--finds a rational function that fits the sequence via polynomialInterpolation of the numerator and the denominator
guessRational = method ()
guessRational (List,PolynomialRing) := RingElement => (L,R) -> (
	n := #L;
	y := for k from 0 to n-1 list L_k;
	a := apply(y, i -> numerator i);
	b := apply(y, i -> denominator i);
	p = polynomialInterpolation(a,R);
	q = polynomialInterpolation(b,R);
	p/q
)
guessRational (List) := RingElement => L -> guessRational(L, QQ[t])
