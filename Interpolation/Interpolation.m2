guessPolynomial = method ()
guessPolynomial (List,PolynomialRing) := RingElement => (L,R) -> (
	if #gens R != 1 then error "must be a single variable polynomial ring";
	n := #L;
	y := for k from 0 to n-1 list (L_k)_1;
	x := for k from 0 to n-1 list (L_k)_0;
	p := for j from 0 to n-1 list (
		y_j * product(0..(j-1) | (j+1) .. (n-1), k -> (t - x_k) / (x_j - x_k))
		);
	sum(p)
)
guessPolynomial (List) := RingElements => L -> guessPolynomial(L, QQ[t])
