--This function computes the parameter test module of a ring, it returns it as a submodule of a canonical ideal.
--this is a slightly modified function originally written by Moty Katzman for "Parameter test ideals of Cohen Macaulay rings"
--it returns the lift of the canonical module to the ambient ring
canonicalIdeal = (R1) -> (
	S1 := ambient R1;
	I1 := ideal(R1);
	d1 := (dim S1) - (dim R1);
	
	canModuleMatrix := relations(prune( Ext^d1(S1^1/I1, S^1)));
	
	answer:=0;
	s1:=syz transpose substitute(canModuleMatrix,R1);
	s2:=entries transpose s1;
	use S1;
	apply(s2, t->
	{
		s3:=substitute(syz gens ideal t,S1);
---		print(s3%canModuleMatrix);
		if ((s3%canModuleMatrix)==0) then
		{
			answer=substitute(mingens ideal t,S1);
			break;
		};
	});
ideal answer
)

--the following function computes the u of a canonical ideal in a polynomial ring
--it uses previous work of Katzman
finduOfIdeal = (canIdeal, defIdeal) -> (
	Ip := frobeniusPower(defIdeal, 1);
	tempIdeal := intersect( (frobeniusPower(canIdeal, 1)) : canIdeal, Ip : defIdeal );
	
	M1 := compress ((gens tempIdeal)%(gens Ip));
	first first entries M1
)

--computes the parameter test submodule of a given ring.  It 
paraTestModuleAmbient = (R1) -> (
	S1 := ambient R1;
	I1 := ideal(R1);
	
	canIdeal := canonicalIdeal(R1);
	
	J1 := findTestElementAmbient(R1);
	tau0 := J1*canIdeal; --this is the starting test element times the ideal
	
	u1 = finduOfIdeal(canIdeal, I1); --this is the multiplying object that gives us (u*omega)^{[1/p]} \subseteq omega.
	
	tauOut = ascendIdeal(tau0, u1, 1);
	
	(tauOut, canIdeal)
)

paraTestIdealAmbient = (R1) -> (
	tempList := paraTestModuleAmbient(R1);
	(tempList#0) : (tempList#1)
)
