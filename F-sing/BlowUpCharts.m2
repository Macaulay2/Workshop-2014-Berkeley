needsPackage "BGG"; --we'll be pushing forward...

getChartOfBlowup = (I1, m1) -> ( --this will return a chart of the blowup based on the m1'th generator of I1.  
--It assumes that Rm is a polynomial ring
	Rm := ring I1;
	n:=rank source vars(Rm); --number of variables
     vv:=first entries vars(Rm); --the variables
     genList := (entries gens I1)#0;
     k1 := coefficientRing(Rm);
     l := #genList;
	YYY := local YYY;
	
	myMon := monoid[ (vv | toList(YYY_1..YYY_l)), MonomialSize=>64];	
	
	R2 := k1 myMon;--we have some relations
	
	newGenList := apply(genList, z-> sub(z, R2) );--the list of my generators in the new ring
	myYYYvarList := drop((first entries vars(R2)), {0, #vv-1});
	
	
	relList1 := apply(0..(m1-1), i -> (newGenList#i - (myYYYvarList#m1)*(myYYYvarList#i)) );
	relList2 := apply((m1+1)..(l-1), i -> (newGenList#i - (myYYYvarList#m1)*(myYYYvarList#i)) );
	
	J1 := ideal(relList1) + ideal(relList2) + ideal(newGenList#m1 - myYYYvarList#m1);--the naive blowup strict transform
	J2 := saturate(J1, ideal(myYYYvarList#m1));--should be the actual strict transform 
	
	R2/J2
)

flattenedReesAlgebra = (I1) -> (--takes an ideal, forms the rees algebra, and returns the rees algebra in two ways, first with flattened variables and the second without
	S1 := reesAlgebra I1;
	J1 := ideal S1;
	tempMonoid := S1.FlatMonoid;
	k1 := coefficientRing (ring I1);
	S2 := k1 tempMonoid;
	
	J2 := sub(J1, S2);
	
	(S2/J2, S1)
)

ascendBlowupSafe = (J1, irrI, u1, ak, ek) -> (--you pass it the ideal you want to ascend (in the ring with flatened variables), and the irrelevant ideal.  It ascends the ideal but it ignores things up to the irrelevant ideal
	S1 := ring J1;
	pp := char S1;
	IN := S1;
	irrI1 := sub(irrI, S1);
	IP := ideal(0_S1);
	u2 := sub(u1, S1);
	
	while (isSubset(IN, IP) == false) do(
     	IP = saturate(IN, irrI1);
--	  error "help";
	  	IN = ethRootSafe(u1, IP, ak, ek)+IP
     );
	IP	
)




--computes the ethRoot of the pushForward of an appropriate twisted parameter
ethImageOfBlowup = (I1, a1, e1) -> (
	if ( not(codim(I1) > 1)) then error "We can only handle ideals of codimension > 1.";

	reesList := flattenedReesAlgebra I1;
	A1 := reesList#0; --this one has flattened variables
	A2 := reesList#1;
 	irrIdeal := sub(ideal(first entries vars A1), A1);
 	singLocus := ideal singularLocus (A1);
 	
 	IRees := sub(I1, A2);
 	
 	
 	canList := canonicalIdeal(A1, FullMap=>true);
 	canIdeal := canList#0;
 	canMap := canList#1;
 	
 	paraTest := paraTestModuleAmbient(A1); --this is dumb, we are doing work twice (computing the Ext)...  It makes mer nervous on the off chance that the presentation of the Ext group is different for two different calls.
 	
 	newMap := map(A1^1/(paraTest#0), source(canMap), matrix(canMap));
 	newKer := (ker newMap)**A2;
 	
 	myDirectImage := HH_0(directImageComplex(IRees^a1*newKer));
 	
 	assert false;
 	
 	directIdeal := moduleToIdeal(myDirectImage, ring I1);
 	print "we got somewhere";
 	if ( codim(directIdeal)==1) then error "This needs to be written in a better way so this error doesn't happen.";
 	print "we got somewhere 2";
 	ethRoot(directIdeal, e1)
)
	
