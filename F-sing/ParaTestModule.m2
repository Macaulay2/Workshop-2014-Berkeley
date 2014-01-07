--This function computes the parameter test module of a ring, it returns it as a submodule of a canonical ideal.
--this is a slightly modified function originally written by Moty Katzman for "Parameter test ideals of Cohen Macaulay rings"
--it returns the lift of the canonical module to the ambient ring
canonicalIdeal = (R1) -> (
	S1 := ambient R1;
	I1 := ideal(R1);
	d1 := (dim S1) - (dim R1);
	
	canModuleMatrix := relations(prune( Ext^d1(S1^1/I1, S1^1)));
	
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

--computes the parameter test submodule of a given ring.  It outputs the parameter test module (as an ideal), it then outputs the canonical module (as an ideal), and finally it outputs the term u used as the action on the ideal
paraTestModuleAmbient = (R1) -> (
	S1 := ambient R1;
	I1 := ideal(R1);
	
	canIdeal := canonicalIdeal(R1);
	
	J1 := findTestElementAmbient(R1);
	tau0 := J1*canIdeal; --this is the starting test element times the ideal
	
	u1 = finduOfIdeal(canIdeal, I1); --this is the multiplying object that gives us (u*omega)^{[1/p]} \subseteq omega.
	
	tauOut = ascendIdeal(tau0, u1, 1);
	
	(sub(tauOut, R1), sub(canIdeal, R1), u1)
)

--computes the parameter test ideal of an ambient ring
paraTestIdealAmbient = (R1) -> (
	tempList := paraTestModuleAmbient(R1);
	(tempList#0) : (tempList#1)
)

--this computes the parameter test module \tau(R, f^t).  It does not assume that R is a polynomial ring.
paraTestModule = (fk, t1) -> (
	R1 := ring fk;
	S1 := ambient R1;
	f1 := sub(fk, S1);
	I1 := ideal R1;
	pp := char R1;
	funList := divideFraction(t1, pp);
	
	aa := funList#0;
	bb := funList#1;
	cc := funList#2;
	
--	tempList := paraTestModuleAmbient(R1);
--	tauAmb := sub(tempList#0, S1);
--	omegaAmb := sub(tempList#1, S1);
--	u1 := tempList#2;

	omegaAmb := canonicalIdeal(R1);
	J1 := findTestElementAmbient(R1)*omegaAmb;
	u1 := finduOfIdeal(omegaAmb, I1);

	
	uPower := 1;
	if (cc != 0) then
		uPower = floor((pp^cc-1)/(pp-1));
	firstTau := J1;
	
	
	if (cc != 0) then	
		firstTau = ascendIdeal(J1*ideal(f1^aa), f1^aa*u1^(uPower), cc)
		--I should write an ascendIdealSafe that works for multiple elements raised to powers...	
	else 
		firstTau = ascendIdeal(J1, u1^(uPower), 1)*ideal(f1^aa);
	
	secondTau := firstTau;
	if (bb != 0) then
		secondTau = ethRoot(u1, firstTau, uPower, bb);

	(sub(secondTau, R1), omegaAmb, u1)
)
