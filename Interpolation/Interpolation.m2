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
polynomialInterpolation (List) := RingElements => L -> polynomialInterpolation(L, QQ[t])

--finds a rational function fitting the sequence through Chebyshev interpolation
rationalInterpolation = method ()
rationalInterpolation (List, PolynomialRing) := RingElement => (L,R) -> (
	if #gens R != 1 then error "must be a single variable polynomial ring";
	t := first gens R;
	n := #L;
	(sum(1..(n-2), j -> ((-1)^j * L_j)/(t-(j+1))) + 1/2 * (L_0)/(t-1) + 1/2 * ((-1)^(n-1) * L_(n-1))/(t-n)) / 
		(sum(1..(n-2), j -> ((-1)^j)/(t-(j+1))) + 1/2 * ((-1)^0)/(t-1) + 1/2 * ((-1)^(n-1))/(t-n))
)
rationalInterpolation (List) := RingElement => L -> rationalInterpolation(L, QQ[t])

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

floaterHormann(List,ZZ) := RingElement => (L,R,d) -> floaterHormann(L,QQ[t],d)

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
guessRational (List) := RingElement => L -> guessRational(L, QQ[t])

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
    for i to n do L = append(L, sum(#C, i -> C_i * L_(#L-#C+i)));
    L
    )

-- Examples
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
