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

cubicSplines = method ()
cubicSplines (List,PolynomialRing) := RingElement => (L,R) -> (
	t := first gens R;
	n := #L;
	P := for i from 1 to n-3 list (
		(2*t^3 - 3*t^2 + 1)*(L_i) + (t^3 - 2*t^2 + t)*(1/2)*(L_(i+1) - L_(i-1)) + (-2*t^3 + 3*t^2)*(L_(i+1)) + (t^3 - t^2)*(1/2)*(L_(i+2) - L_i)
	);
	firstTerm := (2*t^3 - 3*t^2 + 1)*(L_0) + (t^3 - 2*t^2 + t)*(1/2)*(L_1 - L_0) + (-2*t^3 + 3*t^2)*(L_1) + (t^3 - t^2)*(1/2)*(L_2 - L_0);
	finalTerm := (2*t^3 - 3*t^2 + 1)*(L_(n-2)) + (t^3 - 2*t^2 + t)*(1/2)*(L_(n-1) - L_(n-3)) + (-2*t^3 + 3*t^2)*(L_(n-1)) + (t^3 - t^2)*(1/2)*(L_(n-1) - L_(n-2));
	{firstTerm} | P | {finalTerm}
)

--to evaluate at the point you want, say you have a1, a2, ... , an, and u want ai by typing in "i"
--then u want to type in sub(P_(i-1), t => 0) or sub(P_(i-2), t => 1) where t is ur ring generator
--if you want a value in between i and i + 1, say i + j with j in (0,1), then type:
--sub(P_(i-1), t => j)
--this is only a good estimation for values in between 1 and n, where n is the length of the sequence you
--record.
