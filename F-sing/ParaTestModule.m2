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

--computes the parameter test submodule of a given ring.  It 
paraTestModuleAmbient = (R1) -> (
	S1 := ambient R1;
	I1 := ideal(R1);
	
	canIdeal = canonicalIdeal(R1);
	
	
)
