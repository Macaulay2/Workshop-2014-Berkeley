newPackage( "PosChar",
Version => "0.2a", Date => "May 30th, 2015", Authors => {
     {Name => "Erin Bela",
     Email=> "ebela@nd.edu",
     },
     {Name => "DJ Bruce",
     Email=> "djbruce@math.wic.edu",
     HomePage => "http://www.math.wisc.edu/~djbruce/"
     },
     {Name => "Daniel Hernandez",
     Email=> "dhernan@math.utah.edu",
     HomePage=> "http://www.math.utah.edu/~dhernan/"
     },
     {Name => "Zhibek Kadyrsizova",
     Email=> "zhikadyr@umich.edu",
     },
     {Name => "Mordechai Katzman",
     Email=> "m.katzman@sheffield.ac.uk",
     HomePage=> "http://www.katzman.staff.shef.ac.uk/"
     },
     {Name => "Sara Malec",
     Email=> "smalec@gsu.edu"
     },
     {Name => "Karl Schwede",
     Email => "schwede@math.psu.edu",
     HomePage => "http://math.utah.edu/~schwede/"
     },
     {Name => "Pedro Teixeira",
     Email => "pteixeir@knox.edu",
     HomePage => "http://www.knox.edu/academics/faculty/teixeira-pedro.html"
     },
     {Name=> "Emily Witt",
     Email=> "ewitt@umn.edu",
     HomePage => "http://math.umn.edu/~ewitt/"
     }
},
Headline => "A package for calculations of singularities in positive characteristic", DebuggingMode => true, Reload => true )
export{
    "aPower",
    "ascendIdeal", 
    "ascendIdealSafe",
    "ascendIdealSafeList",
    "AscentCount",
    "basePExp",
    "basePExpMaxE",
    "BinomialCheck",
    "binomialFPT",
    "canonicalIdeal",
    "canVector",
    "carryTest",
    "digit", 	 
    "denom",
    "DiagonalCheck", 
    "diagonalFPT",
    "divideFraction",
    "estFPT",
    "ethRoot",
    "ethRootSafe", 		--MK
    "factorList",
    "fancyEthRoot",		--MK
    "fastExp",
    "findCPBelow",
    "findGeneratingMorphisms",     --MK
    "findHSLloci",                 --MK
    "findTestElementAmbient",
    "FinalCheck",
    "findAllCompatibleIdeals", 	--- MK
    "findQGorGen",
    "finduOfIdeal",
    "firstCarry", 
    "FPTApproxList",     
    "FPT2VarHomog",     
    "FPT2VarHomogInternal",
    "fracPart",
    "frobenius",
    "frobeniusPower",
    "fSig",
    "FTApproxList",
    "FTHatApproxList", 
    "FullMap",--specifies whether the full data should be returned
    "getNumAndDenom",
    "genFrobeniusPower",
    "guessFPT",
    "HSL",
    "imageOfRelativeCanonical",
    "imageOfTrace", --doesn't work!
    "isBinomial",
    "isCP",
    "isDiagonal",
    "isFJumpingNumberPoly",
    "isFPTPoly",
    "isFPure",
    "isFRegularPoly",
    "isFRegularQGor",
    "isInLowerRegion",
    "isInUpperRegion",
    "isJToAInIToPe",
    "isSharplyFPurePoly",
    "isMapSplit",
    "MaxExp",
    "minimalCompatible",		--- MK
---    "Mstar",			--- MK
    "multOrder",
    "MultiThread",
    "nonFInjectiveLocus",   --MK
    "Nontrivial",
    "nu",
    "NuCheck",
    "nuHat",
    "nuHatList",
    "nuList",
    "num",
    "Origin",
    "OutputRange",
    "paraTestModule",
    "paraTestModuleAmbient",
    "PrintCP",
    "setFTData",
    "sigmaAOverPEMinus1Poly", 
    "sigmaQGorAmb", --needs optimization
    "sigmaAOverPEMinus1QGor",      --needs optimization
    "splittingField",
    "tauPoly",
    "tauNonPrincipalAOverPEPoly",
    "tauAOverPEMinus1Poly",
    "tauGor",--needs optimization
    "tauGorAmb",--needs optimization
    "tauQGor",--needs optimization
    "tauQGorAmb",--needs optimization
    "taxicabNorm",
    "truncation",
    "truncationBaseP"
}
--This file has "finished" functions from the Macaulay2 workshop at Wake 
--Forest in August 2012.  Sara Malec, Karl Schwede and Emily Witt contributed
--to it.  Some functions, are based on code written by Eric Canton and Moty
-- Katzman
--
--UPDATE January 2014 at Macaulay2 workshop at MSRI:  Daniel Hernandez, Moty 
--Katzman, Karl Schwede, Pedro Teixeira, Emily Witt added more functionality.
--
--UPDATE May 2015 at Macaulay2 workshop atBoise State:  Erin Bela, DJ Bruce,
-- Daniel Hernandez, Zhibek Kadyrsizova, and Emily Witt improved/fixed/added 
--functionality.

----------------------------------------------------------------
--************************************************************--
--Functions for doing particular factorizations with numbers. --
--************************************************************--
----------------------------------------------------------------

--the functions denom and num are needed because M2 complains about 
--computing numerators and denominators of integers (we need this)

denom = method(); --Finds the denominator of a rational number or integer
denom QQ := x -> denominator x;
denom ZZ := x -> 1;

num = method(); --Finds the numerator of a rational number or integer
num QQ := x -> numerator x;
num ZZ := x -> x;

fracPart = (x) -> (x - floor(x)) --Finds the fractional part of a number

aPower = (x,p) -> --find the largest power of p dividing the denominator of x
(
    den:=denom(x);
    a:=1;
    while den%p^a==0 do a=a+1;
    a-1
)
     
-- This function takes in a fraction t and a prime p and spits out a list
-- {a,b,c}, where t = (a/p^b)(1/(p^c-1))
-- if c = 0, then this means that t = (a/p^b)
divideFraction = (t1,pp) -> (
     a := num t1; -- finding a is easy, for now
     b := aPower(t1,pp); -- finding b is easy based upon aPower (written by Emily)
     temp := denom(t1*pp^b); --find the p^c-1 part of the denominator
     pow := 0; --we will look around looking for the power of pp that is 1 mod temp. 
     done := false; --when we found the power, this is set to true.
     if (temp == 1) then done = true; --if there is nothing to do, do nothing.
     while (done==false)  do (
          pow = pow + 1;
	  if (pp^pow % temp == 1) then done = true
     );
     c := pow; --we found c, now we return the list
     if (c > 0) then a = lift(a*(pp^c-1)/temp, ZZ); --after we fix a
     {a,b,c}
)

--Finds the a/pp^e1 nearest t1 from above
findNearPthPowerAbove = (t1, pp, e1) -> (
     ceiling(t1*pp^e1)/pp^e1
)

--Finds the a/pp^e1 nearest t1 from below
findNearPthPowerBelow = (t1, pp, e1) -> (
     floor(t1*pp^e1)/pp^e1
)

--Returns the digits in nn which are nonzero in binary 
--for example, 5 in binary is 101, so this would return {0,2}
--the second term tells me where to start the count, so passing
--5,0 gives {0,2} but 5,1 is sent to {1,3}.  i should be
--used only for recursive purposes
getNonzeroBinaryDigits = (nn, i) -> (
--    error "breakme";
    halfsies := nn//2;
    val1 := nn%2;
    val2 := false; 
    if (halfsies > 0) then val2 = (getNonzeroBinaryDigits(nn//2,i+1));
    if ( (val1 != 0) and (not (val2 === false))) then (
	 flatten {i, val2}
    )
    else if (val1 != 0) then (
	 {i}
    )
    else if ( not (val2 === false)) then (
	 flatten {val2}
    )
    else(
	 false
    )
)

--Returns the entries of myList specified by entryList
--For example, ( {1,2,3}, {0,2}) is sent to {1,3}
getSublistOfList = (myList, entryList) -> (
     --error "help";
     apply( #entryList, i->myList#(entryList#i) )
)

--Returns the power set of a given list, except it leaves out
--the emptyset.  
--For example {2,3,4} becomes { (2),(3),(4),(2,3),(2,4),(3,4),(2,3,4) }
nontrivialPowerSet = (myList) ->(
     apply(2^(#myList)-1, i-> getSublistOfList(myList, getNonzeroBinaryDigits(i+1,0) ) )
)

--This turns a number into a list of factors with repeats
--For example, 12 becomes (2,2,3)
numberToPrimeFactorList = (nn)->(
     prod := factor nn;
     flatten (apply(#prod, i -> toList(((prod#(i))#1):((prod#(i))#0)) ))
)

--Returns a list of all proper factors of nn, for use with sieving...
getFactorList = (nn) ->(
     if (nn < 1) then error "getFactorList: expected an integer greater than 1.";
     powSet := nontrivialPowerSet(numberToPrimeFactorList(nn)); 
     toList ( set apply(#powSet, i->product(powSet#i)) )
)

--This function finds rational numbers in the range of the interval
--with the given denominator
findNumberBetweenWithDenom = (myInterv, myDenom)->(
     upperBound := floor((myInterv#1)*myDenom)/myDenom; 
          --finds the number with denominator myDenom less than the upper 
	  --bound of myInterv
     lowerBound := ceiling((myInterv#0)*myDenom)/myDenom; 
          --finds the number with denominator myDenom greater than the lower
	  -- bound of myInterv
     if (upperBound >= lowerBound) then (
	  --first we check whether there is anything to search for
	  apply( 1+numerator((upperBound-lowerBound)*myDenom), i-> lowerBound+(i/myDenom) )
     )
     else(
	  {}
     )
)

--This function finds rational numbers in the range of 
--the interval; the max denominator allowed is listed. 
findNumberBetween = (myInterv, maxDenom)->(
     divisionChecks :=  new MutableList from maxDenom:true; 
         -- creates a list with maxDenom elements all set to true.
     outList := {};
     i := maxDenom;
     while (i > 0) do (
	  if ((divisionChecks#(i-1)) == true) then --if we need to do a computation..
	      outList = join(outList,findNumberBetweenWithDenom(myInterv, i));
	  factorList := getFactorList(i);
     	  apply(#factorList, j-> (divisionChecks#( (factorList#j)-1) = false) );
	  i = i - 1;
     );
     sort(toList set outList)
)


--Computes the non-terminating base p expansion of an integer
basePExp = (N,p) ->
(
    if N==0 then return {0};
    e:= floor(log_p N);
    E:=new MutableList;
    scan(0..e,i-> 
    	(
     	    a := N//p^(e-i);
     	    E#(e-i) = a;
     	    N = N - (a)*p^(e-i);
    	)
    );
    new List from E
)

--Computes the non-terminating base p expansion of an integer 
--from digits zero to e-1 (little-endian first)
basePExpMaxE = (N,p,e1) ->
(
    e:=e1-1;
    E:=new MutableList;
    scan(0..e,i-> 
    	(
     	    a := N//p^(e-i);
     	    E#(e-i) = a;
     	    N = N - (a)*p^(e-i);
    	)
    );
    new List from E
)


---------------------------------------------------------------
--***********************************************************--
--Basic functions for computing powers of elements in        --
--characteristic p>0.                                        --
--***********************************************************--
---------------------------------------------------------------


--Computes powers of elements in char p>0, using that Frobenius is an endomorphism
-- If N = N_0 + N_1 p + ... + N_e p^e, then this computes f^N as f^N = f^(N_0) f^(N_1)^p ... (f^(N_e))^(p^e)

fastExp = (f,N) ->
(
     p:=char ring f;
     E:=basePExp(N,p);
     product(#E, e -> (sum(terms f^(E#e), g->g^(p^e))))
)

-- Old version of fastExp. 
-- If N = N_0 + N_1 p + ... + N_e p^e, then this computes f^N as f^N = f^(N_0) (f^p)^(N_1) ... (f^(p^e))^(N_e)
fastExpOld = (f,N) ->
(
     p:=char ring f;
     E:=basePExp(N,p);
     product(#E, e -> (sum(terms f, g->g^(p^e)))^(E#e) )
)


---------------------------------------------------------------
--***********************************************************--
--Functions for computing \nu_I(p^e), \nu_f(p^e), and using  --
--these to compute estimates of FPTs.                        --
--***********************************************************--
---------------------------------------------------------------

-- If I is contained in Rad(J) then this finds the minimal N 
-- such that I^n is conatined in J

effRad = (I1,J1) ->(
       d1 := 1;
       if isSubset(I1, radical(J1))==false then (print "Error: I Not Contained in Rad(J)")
       else(
       while isSubset(I1^d1,J1) == false do (
	   d1 = d1+1
   	   );
        d1)
)


nuList = method()

nuList(Ideal, Ideal,  ZZ) := (I1, J1, e1) -> ( --this is a faster nuList computation, it tries to do a smart nu list computation
	d1 := 0;
	p1 := char ring I1;
	local top;--for the binary search
	local bottom;--for the binary search
 	local middle;--for the binary search
	local answer; --a boolean for determining if we go up, or down
	
	if isSubset(I1, radical(J1))==false then (print "Error: Nu Undefined")
	else(
	myList := new MutableList;
	nuPrev := effRad(I1,J1);
	N := numgens(trim(J1));
	top = nuPrev*(N*p1-1);
	bottom = 0;
	
	for d1 from 1 to e1 do (
		while (top - 1 > bottom) do (--the bottom value is always not in m, the top is always in m
			middle := floor((top + bottom)/2);
			answer = isSubset(I1^middle, frobeniusPower(J1, d1));
			if (answer == false) then bottom = middle else top = middle;
		);
		nuPrev = bottom;
		myList#(d1-1) = nuPrev;
		top = (nuPrev+1)*(N*p1-1);
		bottom = p1*nuPrev;
	);
	toList myList)
)

nuList(Ideal, ZZ) := (I1, e1) -> (
    nuList(I1,ideal(first entries vars ring I1), e1)
    )

nuList(RingElement,ZZ) := (f,e) -> nuList(ideal(f),e)

nu = method()

nu(Ideal, Ideal, ZZ) := (I1, J1, e1) -> ( --this does a fast nu computation
	p1 := char ring I1;
	local top;--for the binary search
	local bottom;--for the binary search
	local middle;--for the binary search
	local answer; --a boolean for determining if we go up, or down 
	if isSubset(I1, radical(J1))==false then (print "Error: Nu Undefined")
	else(
	N := 0;
	myList := new MutableList;
	nuPrev := effRad(I1,J1);
	N = numgens(trim(J1));
	top = nuPrev*N*p1^e1-1;
	bottom = 0;
			
	while (top - 1 > bottom) do (--the bottom value is always not in m, the top is always in m
		middle = floor((top + bottom)/2);
		answer = isSubset(I1^middle, frobeniusPower(J1, e1));
		if (answer == false) then bottom = middle else top = middle;
	);
	bottom)
)

nu(Ideal, ZZ) := (I1, e1) -> (
    nu(I1,ideal(first entries vars ring I1), e1)
    )

nu(RingElement, ZZ) := (f,e) -> nu(ideal(f),e)

--Approximates the F-pure Threshold
--Gives a list of nu_I(p^d)/p^d for d=1,...,e
FPTApproxList = method();
FPTApproxList (Ideal,ZZ) := (I,e) ->
(
     p := char ring I;
     apply(#nuList(I,e), i->((nuList(I,e))#i)/p^(i+1)) 
)
FPTApproxList (RingElement,ZZ) := (f,e) -> FPTApproxList(ideal(f),e)

--Approximates the F-Threshold with respect to an ideal J
--More specifically, this gives a list of nu_I^J(p^d)/p^d for d=1,...,e

FTApproxList = method();

FTApproxList(Ideal,Ideal,ZZ) := (I1,J1,e1) ->
(
    if isSubset(I1, radical(J1))==false then (print "Error: F-Threshold Undefined")
    else(
     p1 := char ring I1;
     apply(#nuList(I1,J1,e1), i->((nuList(I1,J1,e1))#i)/p1^(i+1)))
)

FTApproxList (RingElement,Ideal,ZZ) := (f1,J1,e1) -> FTApproxList(ideal(f1),J1,e1)

---------------------------------------------------------------
--***********************************************************--
--This seems wrong becuase the top bound is incorrect for    --
--ideals that are not principal. (5/29/15)-Boise
--***********************************************************--
--Functions for computing \nu_I(p^e), \nu_f(p^e), and using  --
--these to compute estimates of FPTs.                        --
--***********************************************************--
---------------------------------------------------------------

--isJToAInIToPe = (J1, a1, I1, e1) -> (--checks whether or not f1^a1 is in I1^(p^e1).  It seems to be much faster than raising f1 to a power
--	root := ethRoot(J1, a1, e1); 
--	
--	isSubset(root, I1)
--)
--
--nuList = method()
--
--nuList(Ideal, Ideal,  ZZ) := (I1, J1, e1) -> ( --this is a faster nuList computation, it tries to do a smart nu list computation
--	d1 := 0;
--	p1 := char ring I1;
--	local top;--for the binary search
--	local bottom;--for the binary search
--	local middle;--for the binary search
--	local answer; --a boolean for determining if we go up, or down
--	N := 0;
--	myList := new MutableList;
--	curPower := 0;
	
--	for d1 from 1 to e1 do (
--		if (curPower == 0) then curPower = 1 else curPower = p1^(d1-1)*curPower;
--		top = p1;
--		bottom = 0;
--		
--		while (top - 1 > bottom) do (--the bottom value is always not in m, the top is always in m
--			middle := floor((top + bottom)/2);
--			answer = isJToAInIToPe(I1, curPower + middle, J1, d1);
--			print "Here we are";
--			print (bottom, top);
--			if (answer == false) then bottom = middle else top = middle
--			--print "Here we are"
--		);
--	--	print (bottom, top, curPower);
--		curPower = curPower + bottom;
--		myList#(d1-1) = curPower;
--		curPower = curPower*p1;
--	);
--	toList myList
--)
--
--nuList(Ideal, ZZ) := (I1, e1) -> (
--    nuList(I1,ideal(first entries vars ring I1), e1)
--    )
--
--nuList(RingElement,ZZ) := (f,e) -> nuList(ideal(f),e)
--
--nu = method()
--
--nu(Ideal, Ideal, ZZ) := (I1, J1, e1) -> ( --this does a fast nu computation
--	p1 := char ring I1;
--	local top;--for the binary search
--	local bottom1;--for the binary search
--	local middle;--for the binary search
--	local answer; --a boolean for determining if we go up, or down 
--	N := 0;
--	myList := new MutableList;
--	curPower := 0;
--	
--	bottom1 = 0;
--	top = p1^e1;		
--	while (top - 1 > bottom1) do (--the bottom value is always not in m, the top is always in m
--		middle = floor((top + bottom1)/2);
--		answer = isJToAInIToPe(I1, middle, J1, e1);
--		if (answer == false) then bottom1 = middle else top = middle;
--	);
--	bottom1
--)
--
--nu(Ideal, ZZ) := (I1, e1) -> (
--    nu(I1,ideal(first entries vars ring I1), e1)
--    )
--
--nu(RingElement, ZZ) := (f,e) -> nu(ideal(f),e)
--
--Approximates the F-pure Threshold
--Gives a list of nu_I(p^d)/p^d for d=1,...,e
--FPTApproxList = method();
--FPTApproxList (Ideal,ZZ) := (I,e) ->
--(
--     p := char ring I;
--     apply(#nuList(I,e), i->((nuList(I,e))#i)/p^(i+1)) 
--)
--FPTApproxList (RingElement,ZZ) := (f,e) -> FPTApproxList(ideal(f),e)
--
--Approximates the F-Threshold with respect to an ideal J
--More specifically, this gives a list of nu_I^J(p^d)/p^d for d=1,...,e
--
--FTApproxList = method();
--
--FTApproxList(Ideal,Ideal,ZZ) := (I1,J1,e1) ->
--(
--    if isSubset(I1, radical(J1))==false then (print "Error: F-Threshold Undefined")
--    else(
--     p1 := char ring I1;
--     apply(#nuList(I1,J1,e1), i->((nuList(I1,J1,e1))#i)/p1^(i+1)))
--)
--
--FTApproxList (RingElement,Ideal,ZZ) := (f1,J1,e1) -> FTApproxList(ideal(f1),J1,e1)
--

---------------------------------------------------------------
--***********************************************************--
--Functions for computing \nuHat_I(p^e), \nHatu_f(p^e), and  --
--using these to compute estimates of FThat's.                        --
--***********************************************************--
---------------------------------------------------------------


nuHatList = method()

nuHatList(Ideal, Ideal,  ZZ) := (I1, J1, e1) -> ( --this is a faster nuList computation, it tries to do a smart nu list computation
	d1 := 0;
	p1 := char ring I1;
	local top;--for the binary search
	local bottom;--for the binary search
 	local middle;--for the binary search
	local answer; --a boolean for determining if we go up, or down
    	if isSubset(I1, radical(J1))==false then (print "Error: NuHat Undefined")
	else(
	myList := new MutableList;
	nuPrev := effRad(I1,J1);
--	N = numgens(trim(J1));
	top =nuPrev*p1;
	bottom = 0;
	
	
	for d1 from 1 to e1 do (
		while (top - 1 > bottom) do (--the bottom value is always not in m, the top is always in m
			middle := floor((top + bottom)/2);
			answer = isSubset(genFrobeniusPower(I1,middle), frobeniusPower(J1,d1));
			if (answer == false) then bottom = middle else top = middle;
		);
		nuPrev = bottom;
		myList#(d1-1) = bottom;
		top = p1*(nuPrev+1);
		bottom = p1*nuPrev;
	);
	toList myList)
)

nuHatList (Ideal, ZZ) := (I1, e1) -> (
    nuHatList(I1,ideal(first entries vars ring I1), e1)
    )

nuHat = method()

nuHat (Ideal, Ideal, ZZ) := (I1, J1, e1) -> ( --this does a fast nu computation
	p1 := char ring I1;
	local top;--for the binary search
	local bottom;--for the binary search
	local middle;--for the binary search
	local answer; --a boolean for determining if we go up, or down
    	if isSubset(I1, radical(J1))==false then (print "Error: NuHat Undefined")
	else(
	myList := new MutableList;
	nuPrev := effRad(I1,J1);
--	N = numgens(trim(J1));
	top = nuPrev*p1^e1;
	bottom = 0;
			
	while (top - 1 > bottom) do (--the bottom value is always not in m, the top is always in m
		middle = floor((top + bottom)/2);
		answer = isSubset(genFrobeniusPower(I1, middle), frobeniusPower(J1,e1));
		if (answer == false) then bottom = middle else top = middle;
	);
	bottom)
)

nuHat (Ideal, ZZ) := (I1, e1) -> (
    nuHat(I1,ideal(first entries vars ring I1), e1)
    )


--Aproximates the F-Threshold with respects to an ideal J

FTHatApproxList = method();

FTHatApproxList(Ideal,Ideal,ZZ) := (I1,J1,e1) ->
(
    if isSubset(I1, radical(J1))==false then (print "Error: F-Threshold Undefined")
    else(
     p1 := char ring I1;
     apply(#nuHatList(I1,J1,e1), i->((nuHatList(I1,J1,e1))#i)/p1^(i+1)))
)

FTHatApproxList (RingElement,Ideal,ZZ) := (f1,J1,e1) -> FTHatApproxList(ideal(f1),J1,e1)

---------------------------------------------------------------
--***********************************************************--
--Basic functions for Frobenius powers of ideals and related --
--constructions (colons).                                    --
--***********************************************************--
---------------------------------------------------------------
 
 
--The following raises an ideal to a Frobenius power; it was written by Moty Katzman
frobeniusPower=method()

frobeniusPower(Ideal,ZZ) := (I1,e1) ->(
     R1:=ring I1;
     p1:=char R1;
     local answer;
     G1:=first entries gens I1;
     if (#G1==0) then answer=ideal(0_R1) else answer=ideal(apply(G1, j->fastExp(j, (p1^e1))));
     answer
);

frobeniusPower(Matrix,ZZ) := (M,e) ->
(
    p:=char ring M;
    matrix apply(entries M,u->apply(u,j->j^(p^e)))
)

-- The following raises an ideal to a generalized Frobenius power i.e. if N=n_0+n_1P+...+n_eP^e then
-- I^N = I^n_0*(I^n_1)^[P]*...*(I^n_e)^[P^e].

genFrobeniusPower = (I1,e1) ->(
     R1:=ring I1;
     p1:=char R1;
     E1:=basePExp(e1,p1);
     local answer;
     answer = product(#E1, q -> (frobeniusPower(I1^(E1#q),q)));
     answer
)

-- This function computes the element in the ambient ring S of R=S/I such that
-- I^{[p^e]}:I = (f) + I^{[p^e]}
-- If there is no such unique element, the function returns zero

findQGorGen=method();
findQGorGen (Ring,ZZ) := (Rk,ek) -> (
     Sk := ambient Rk; -- the ambient ring
     Ik := ideal Rk; -- the defining ideal
     pp := char Sk; --the characteristic
     Ikpp := frobeniusPower(Ik,ek);
     
     J1 := trim (Ikpp : Ik); --compute the colon
     Tk := Sk/Ikpp; --determine the ideal in 
     
     J2 := trim sub(J1, Tk);
     
     Lk := first entries gens J2;
     
     nk := #Lk;
     val := 0_Sk;
     
     if (nk != 1) then (
	  error "findGorGen: this ring does not appear to be (Q-)Gorenstein, or
	   you might need to work on a smaller chart.  Or the index may not divide p^e-1
	   for the e you have selected.";
     )
     else (
	  val = lift(Lk#0, Sk);
     );    
     val 
)
findQGorGen(Ring) := (R2) -> ( findQGorGen(R2, 1) )


---------------------------------------------------------------------
--*****************************************************************--
--Functions for computing the F-pure threshold of a diagonal       --
--or binomial hypersurface using the algorithms of Daniel Hernandez--
--(written by E. Witt)                                             --                      
--*****************************************************************--
---------------------------------------------------------------------

--Gives the e-th digit of the non-terminating base p expansion of x in [0,1] 
digit = method()

digit (ZZ,QQ,ZZ) := (e, x, p) -> 
(
     y := 0;
     if fracPart(p^e*x) != 0 then y = floor(p^e*x) - p*floor(p^(e-1)*x);
     if fracPart(p^e*x) == 0 then y = floor(p^e*x) - p*floor(p^(e-1)*x) - 1;
     if fracPart(p^(e-1)*x) == 0 then y = p-1;
     y
)

--digit (ZZ,List,ZZ) threads over lists.
digit (ZZ,List,ZZ) := (e,u,p) -> apply(u,x->digit(e,x,p))

--Gives the e-th truncation of the non-terminating base p expansion of a nonnegative 
--rational x as a fraction
truncation = method()

truncation (ZZ,QQ,ZZ) := (e,x,p) -> (ceiling(p^e*x)-1)/p^e

--truncation (ZZ,List,ZZ) threads over lists.
truncation (ZZ,List,ZZ) := (e,u,p) -> apply(u,x->truncation(e,x,p))

--Gives the first e digits of the non-terminating base p expansion of x in [0,1]
--as a list
truncationBaseP = (e,x,p) -> 
(
     L := new MutableList;
     for i from 0 to e-1 do L#i = digit(i+1,x,p);
     L
)

--Given a rational number x, if a is the power of p dividing its denomiator, 
--finds an integer b so that p^a(p^b-1)*a is an integer
bPower = (n,p) ->
(
     if aPower(n,p)>0 then n = n*p^(aPower(n,p));
     denom(n)
)

--Given a vector w={x,y}, x and y rational in [0,1], returns a number of digits 
--such that it suffices to check to see if x and y add without carrying in base p
carryTest = (w,p) ->
(
     c := 0; for i from 0 to #w-1 do c = max(c, aPower(w#i, p));
     d := 1; for j from 0 to #w-1 do if bPower(w#j, p)!=0 then d = lcm(d, bPower(w#j, p));
     c+d+1
)

--Given a vector w={x,y} of rational integers in [0,1], returns the first spot 
--e where the x and y carry in base p; i.e., 
--(the e-th digit of x)+(the e-th digit of y) >= p
firstCarry = (w,p) ->
(     
    i:=0;
    d:=0;
    carry:=0;
    zeroTest := false;
    for j from 0 to #w-1 do if w#j == 0 then zeroTest=true;
    if zeroTest == true then carry = -1 else
     (
	       i = 0; while d < p and i < carryTest(w,p) do 
	       (
	       	    i = i + 1;
	       	    d = 0; for j from 0 to #w-1 do  d = d + digit(i,w#j,p);
	   	);
      	       if i == carryTest(w,p) then i = -1;
      	       carry = i;
      );
      carry
)

--Given a vector w, returns a vector of the reciprocals of the entries of w
reciprocal = w ->
(
     v := new MutableList from w;
     for c from 0 to #w-1 do v#c = 1/w#c;
     v
)

--Computes the F-pure threshold of a diagonal hypersurface 
--x_1^(a_1) + ... +x_n^(a_n) using Daniel Hernandez' algorithm
diagonalFPT = f ->
(
     p := char ring f;
     w := apply(terms f, g->first degree(g));
     y := 0; if firstCarry(reciprocal(w),p)==-1 then for i from 0 to #w-1 do y = y + 1/w#i else
     (
	  x := 0; for c from 0 to #w-1 do x = x + truncation(firstCarry(reciprocal(w),p)-1, 1/w#c, p); 
	  y = x+1/p^(firstCarry(reciprocal(w),p)-1);
     );
     y
)

--Given a polynomial f, outputs a list of multi-degrees (under the usual grading)
--of f as lists
multiDegree = f ->
(
     variables := first entries vars ring f;
     apply(terms f, g -> apply(#variables, i ->degree(variables#i,g)))
)

--Determines whether a polynomial f is diagonal; i.e., of the form 
--x_1^(a_1)+...+x_n^(a_n) (up to renumbering variables)
isDiagonal = f ->
(
     d := multiDegree(f);
     alert1 := true;
     alert2 := true;
     for i from 0 to #d-1 do
     (
	  for j from 0 to #(d#0)-1 do
	  (
	       if (d#i)#j!=0 and alert1==false then alert2=false;
	       if (d#i)#j!=0 and alert1==true then alert1=false;
	  );
     alert1=true;
     );
     for j from 0 to #(d#0)-1 do
     (
	  for i from 0 to #d-1 do 
	  (
     	       if alert1==false and (d#i)#j!=0 then alert2=false;
     	       if alert1==true and (d#i)#j!=0 then alert1=false;
	  );
     alert1=true;
     );
     alert2
)

--Given input vectors v={a_1,...,a_n} and w={b_1,...,b_n}, gives the
--corresponding vectors that omit all a_i and b_i such that a_i=b_i
factorOutMonomial = (v,w) ->
(
     v1 := new MutableList;
     w1 := new MutableList;
     c := 0; i := 0; for i from 0 to #v-1 do (if v#i != w#i then (v1#c = v#i; w1#c = w#i; c = c+1; ); );
     (v1,w1)
)

--Given input vectors v={a_1,...,a_n} and w={b_1,...,b_n}, gives the
--vector of the a_i for which a_i=b_i
monomialFactor = (v,w) ->
(
     a := new MutableList;
     c := 0; i := 0; for i from 0 to #v-1 do (if v#i == w#i then (a#c = v#i; c = c+1; ); );
     a
)

--Given two vectors v={v0,v1} and w={w0,w1} in the real plane, finds 
--the intersection of the associated lines v0*x+w0*y=1 and v1*x+w1*y=1
twoIntersection = (v,w) ->
(
     if v#0*w#1-v#1*w#0 != 0 then 
     (
	  x := (w#1-w#0)/(v#0*w#1-v#1*w#0);
	  y := (v#0 - v#1)/(v#0*w#1-v#1*w#0);
	  P := {x,y};
     ) else P = {0,0};
P
)

--Given two vectors v={v0,...vn} and w={w0,...,wn}, list the 
--intersections of all lines vi*x+wi*y=1 and vj*x+wj*y=1
allIntersections = (v,w) ->
(
     L := new MutableList;
     c := 0;
     for i from 0 to #v-1 do
     (
	  for j from i+1 to #v-1 do 
     	  (
     	       if twoIntersection({v#i,v#j}, {w#i,w#j}) != {0,0} then 
     	       (
	  	    L#c = twoIntersection({v#i,v#j}, {w#i,w#j});
	  	    c = c+1;
     	       );
	  );
     );
     for i from 0 to #v-1 do
     (
	  if v#i !=0 then  
	  (
	       L#c = {1/(v#i), 0};
	       c = c + 1;
	  );
     );
     for i from 0 to #v-1 do
     (
	  if w#i !=0 then  
	  (
	       L#c = {0, 1/(w#i)};
	       c = c + 1;
	  );
     ); 
     K := new MutableList;
     c = 0; for i from 0 to #L-1 do
     (
	  if (L#i)#0 >= 0 and (L#i)#1 >=0 then (K#c = {(L#i)#0, (L#i)#1}; c = c+1);
     );
     K
)

--Given a point a=(x,y) in the real plane and two vectors v={v0,...,vn} and w={w0,...,wn}, checks whether a is in the polytope defined by the equations vi*x+wi*y<=1
isInPolytope = (a,v,w) ->
(
     alert := true;
     for i from 0 to #v-1 do
     (
	  if v#i*a#0 + w#i*a#1 > 1 then alert = false;
     );
     alert
)


--Given a point a=(x,y) in the real plane and two vectors
--v={v0,...,vn} and w={w0,...,wn}, checks whether a is in
--the polytope defined by the equations vi*x+wi*y<=1
isInInteriorPolytope = (a,v,w) ->
(
     alert := true;
     for i from 0 to #v-1 do
     (
	  if v#i*a#0 + w#i*a#1 >= 1 then alert = false;
     );
     alert
)

--Given two vectors v and w of the same length, outputs 
--a list of the defining vectors of the polytope as in isInPolytope
polytopeDefiningPoints = (v,w) ->
(
     L := allIntersections(v,w);
     K := new MutableList;
     c := 0;
     for j from 0 to #L-1 do
     (
	  if isInPolytope(L#j,v,w) == true then (K#c = {(L#j)#0, (L#j)#1}; c = c+1;)
     );
     K
)

--Given a list of coordinates in the real plane, 
--outputs the one with the largest coordinate sum
maxCoordinateSum = L ->
(
     K := new MutableList from {0,0};
     for i from 0 to #L-1 do if (L#i)#0 + (L#i)#1 > K#0 + K#1 then K = {(L#i)#0, (L#i)#1};
     K
)

--Finds the "delta" in Daniel Hernandez's algorithm
--for F-pure thresholds of binomials
dCalculation = (w,N,p) ->
(
     d := 0; for j from 0 to #w-1 do  d = d + digit(N+1,w#j,p);
     i := N; while d > p-2 do 
     (
	  d = 0; for j from 0 to #w-1 do  d = d + digit(i,w#j,p);
	  i = i - 1;
     );
     i + 1
)

--Given the "truncation" point (P1,P2) and two vectors 
--defining the binomial v and w, outputs the "epsilon" in 
--Daniel Hernandez's algorithm for F-thresholds of binomials
calculateEpsilon = (P1,P2,v,w) ->
(
     X := new MutableList;
     Y := new MutableList;
     c:=0; d := 0; for i from 0 to #v-1 do 
     (
	  if w#i != 0 then 
     	  (
	       X#c = (1 - (v#i)*(P2#0) - (w#i)*(P2#1))/(w#i);
	       c = c+1;
	  );
          if v#i != 0 then 
	  (
	       Y#d = (1 - (v#i)*(P1#0) - (w#i)*(P1#1))/(v#i);
	       d = d+1;
	  );
     );
     i:=0;
     epsilon:=0;
     if isInInteriorPolytope(P1,v,w)==false and isInInteriorPolytope(P2,v,w)==false then epsilon = -1 else
     (
	  if isInInteriorPolytope(P1,v,w)==false then for i from 0 to #v-1 do X#1 = 0;
	  if isInInteriorPolytope(P2,v,w)==false then for i from 0 to #v-1 do Y#1 = 0;
	  M := X#0; 
	  N := Y#0;
	  for i from 1 to #X-1 do M = min(M, X#i);
	  for j from 1 to #Y-1 do N = min(N, Y#j);
	  epsilon = max(M,N); 
     );
     epsilon
)

--Computes the FPT of a binomial, based on the work of Daniel Hernandez 
--(implemented by Emily Witt)
binomialFPT = g ->
(
     p := char ring g;
     v := (multiDegree(g))#0;
     w := (multiDegree(g))#1;
     FPT := 0;
     f := monomialFactor(v,w);
     x := factorOutMonomial(v,w);
     v = x#0;
     w = x#1;
     Q := maxCoordinateSum(polytopeDefiningPoints(v,w));
     if Q#0+Q#1 > 1 then FPT = 1 else
     (
	  L :=  firstCarry(Q,p);
	  if L == -1 then FPT = Q#0+Q#1 else
     	  (
     	       d := dCalculation(Q,L-1,p);
     	       P := (truncation(d,Q#0,p),  truncation(d,Q#1,p));
     	       P1 := {P#0, P#1+1/p^d};
     	       P2 := {P#0+1/p^d,P#1};
     	       FPT = truncation(L-1,Q#0+Q#1,p);
     	       if calculateEpsilon(P1,P2, v, w) != -1 then FPT = FPT +  calculateEpsilon(P1, P2, v, w);
     	  );
     );
     monFPT := infinity;
     for i from 0 to #f-1 do (if f#i!=0 then monFPT = min(monFPT, 1/(f#i)););
     if monFPT != 0 then FPT = min(FPT, monFPT);
     FPT
)

--Returns true if the polynomial is binomial.
isBinomial = f ->
(
     alert := true;
     if #(terms f)>2 then alert = false;
     alert
)

---------------------------------------------------------------------
--*****************************************************************--
--Functions for computing F-thresholds of forms in two variables   --
--over finite fields. Based on the work of Hernandez and Teixeira. -- 	                                           --                      
--*****************************************************************--
---------------------------------------------------------------------

{*
    Remark: At this point, only commands for computations of F-pure thresholds are
    implemented. Eventually computations of F-thresholds with respect to more general
    ideals will be implemented, and perhaps of more general polynomials. Some structures 
    and functions below are already designed to handle such greater generality. 
*}
    
{*
    Types and auxiliary commands
*}

--FTData is a HashTable that stores the data necessary in F-threshold calculations
--(for conveniently passing those data from one function to another).
--It contains the following keys:
--    "ring": the ring of the polynomial in question;
--    "char": the characteristic of ring;
--    "ideal": the ideal with respect to which we want to compute the F-threshold;
--    "gens": the generators of the ideal;
--    "polylist": a list of the (non-associated) factors of the polynomial in question;
--    "numpolys": the number of factors.
FTData = new Type of HashTable

--setFTData takes a list of generators of the ideal or the ideal itself and a list
--    of polynomials, and builds an FTData from them.
setFTData = method()

setFTData (List,List) := (gen,polylist) -> 
(
    	A:=ring gen_0;
    	p:= char A;	
	new FTData from {"char"=>p,"ring"=>A, "ideal"=>ideal gen, "gens" => gen,
	    "numpolys"=>#polylist,"polylist"=>polylist}
)

setFTData (Ideal,List) := (I,polylist) -> setFTData(I_*,polylist)

{*
    Tests and auxiliary functions
*}

--isInUpperRegion(a,q,S)/isInUpperRegion(u,S) test if the point u=a/q is in the
--upper region attached to S. Suppose I is the ideal of the FTData S under consideration 
--and L={L_1,...,L_n} is the "polylist". Then a point a/q (where a=(a_1,...,a_n) is a 
--nonnegative integer vector and q a power of "char") is in the "upper region" if 
--L_1^(a_1)...L_n^(a_n) is in I^[q]; otherwise it is in the lower region.
isInUpperRegion = method()

isInUpperRegion (List,ZZ,FTData) := (a,q,S) -> 
(
--    p:=S#"char";
--    try e:=lift(log_p q,ZZ) else error "isInUpperRegion: expected second entry to be a power of the characteristic.";
--    frob:=frobeniusPower(S#"ideal",e);
-- **** Had to modify this, because M2 thinks, e.g., that log_5 125 is not an integer.
    frob:=ideal apply(S#"gens",f->f^q);
    F:=product(S#"polylist",a,(f,i)->fastExp(f,i));
    (F % frob) == 0
)

isInUpperRegion (List,FTData) := (u,S) ->
    isInUpperRegion append(getNumAndDenom(u),S)

--isInLoweRegion(a,q,S)/isInLoweRegion(u,S) test if the point u=a/q is in the
--lower region attached to S.
isInLowerRegion = method()

isInLowerRegion (List,ZZ,FTData) := (a,q,S) -> not isInUpperRegion(a,q,S)

isInLowerRegion (List,FTData) := (u,S) -> not isInUpperRegion(u,S)

--neighborInUpperRegion(a,q,S)/neighborInUpperRegion(u,S): auxiliary commands that, 
--given a point u=a/q in the upper region, try to find a "neighbor" of the form 
--(a-e_i)/q that also lies in the upper region. If the search is successful, they return
--the first such neighbor found; otherwise they return nothing.
neighborInUpperRegion = method()

neighborInUpperRegion (List,ZZ,FTData) := (a,q,S) ->
(
    if isInLowerRegion(a,q,S) then (error "Expected point in the upper region.");
    n := S#"numpolys";
    posEntries := positions(a,k->(k>0));
    found := false;
    i:=0;
    local candidate;
    local neighbor;
    while ((not found) and (i<#posEntries)) do 
    (
	candidate=a-canVector(posEntries_i,n);
	if isInUpperRegion(candidate,q,S) then (found=true; neighbor=candidate);
	i=i+1;
    );
    if (not found) then null else (neighbor,q)
)

neighborInUpperRegion (List,FTData) := (u,S) -> 
(
    nbr:=neighborInUpperRegion append(getNumAndDenom(u),S);
    if nbr===null then nbr else (nbr_0)/(nbr_1)
)

--isCP(a,q,S)/isCP(u,S) test if u=a/q is a critical point, that is, if u is in the
--upper region but each neighbor (a-e_i)/q (where a_i>0) is not.
isCP = method()

isCP (List,ZZ,FTData) := (a,q,S) -> 
(
    if isInLowerRegion(a,q,S) then return false;
    neighborInUpperRegion(a,q,S)===null
)

isCP (List,FTData) := (u,S) -> isCP append(getNumAndDenom(u),S)

--findCPBelow(u,S) takes a point u in the upper region attached to S and finds a 
--critical point <= u with the same denominator.
findCPBelow = method()

findCPBelow (List,FTData) := (pt,S) ->
(
    if isInLowerRegion(pt,S) then (error "The point must be in the upper region.");
    nbr:=neighborInUpperRegion(pt,S);
    if nbr===null then return pt else findCPBelow(nbr,S)
)

{*
    Computation of FPTs
*}

--FPT2VarHomogInternal({a1,...an},S): if S#"polylist={L1,...,Ln} is a list of linear
--forms, FPT2VarHomogInternal({a1,...an},S) finds the FPT of the polynomial
--F=L1^(a1)...Ln^(an)
FPT2VarHomogInternal = method(Options => {MaxExp => infinity, PrintCP => false, Nontrivial => false})

FPT2VarHomogInternal (List,FTData) := opt -> (a,S) ->
(
    deg:=taxicabNorm(a);
    pos:=positions(a,k->(k>=deg/2));
    if (pos!={}) then return(1/a_(pos_0)); 
       -- if some multiplicity a_i is "too big", return 1/a_i
    p:=S#"char";
    den:=denom(2/deg);
    local mult;
    if (opt.Nontrivial) then mult = infinity
    else
    ( 
    	if gcd(S#"char",den)==1 then mult = multOrder(p,den)
	else
	(
	    F:=product(S#"polylist",a,(f,i)->f^i);
	    if isFPTPoly(F,2/deg) then (return (2/deg))
	    else mult = infinity
	)
    );    
    rng:=S#"ring";
    polys:=S#"polylist";
    I:=S#"ideal";
    ideals:={I};
    e:=0;
    dgt:=0;
    u:=2*a/deg;
    while (I != ideal(1_rng) and e < (opt.MaxExp) and e < mult) do 
    (
	e=e+1;
	dgt=digit(e,u,p);
	I=frobeniusPower(I,1):product(polys,dgt,(f,k)->f^k);
	ideals=append(ideals,I)
    );
    if I!=ideal(1_rng) then 
    (
	if e == mult then (return (2/deg)) 
	else error "Reached MaxExp."
    );    
    e0:=e-1;
    S1:=setFTData(ideals_e0,polys);
    cp:=findCPBelow(dgt/p,S1); 
    	--if some coordinate of cp is 0, its magnification may not be a CP
    while product(cp)==0 do 
    (
	e0=e0-1;
        -- zoom out one step and look for CP again
    	S1=setFTData(ideals_e0,polys);
	cp=findCPBelow(cp/p+digit(e0+1,u,p)/p,S1) 
    );
    cp=cp/p^e0+truncation(e0,u,p); -- "zoom out"
    if opt.PrintCP then print(toString cp);
    max apply(cp,a,(c,k)->c/k)
)

-----------------------
FPT2VarHomog = method(Options => {MaxExp => infinity, PrintCP => false})

--FPT2VarHomog(RingElement)
--FPT(F) computes the F-pure threshold of a form F in two variables. 
--KNOWN ISSUE: if the splitting field of F is too big, factor will not work.
FPT2VarHomog (RingElement) :=  opt ->  F ->
(    
   if not isNonConstantBinaryForm(F) then (
	error "FPT2VarHomog expects a nonconstant homogeneous polynomial in 2 variables."
    );
    -- because factoring is the weakness of this algorithm, we try to avoid it
    -- by first checking if fpt=lct
    deg:=(degree F)_0;
    if isFPTPoly(F,2/deg) then return 2/deg;
    R:=ring F;
    vv:=R_*;
    kk:=splittingField(F);
    a:= symbol a;
    b:= symbol b;
    S:=kk[a,b];
    G:=sub(F,{(vv#0)=>a,(vv#1)=>b});
    (L,m):=toSequence transpose factorList(G);
    FPT2VarHomogInternal(m,setFTData(S_*,L),MaxExp=>(opt.MaxExp),PrintCP=>(opt.PrintCP),Nontrivial=>true)
)

--FPT2VarHomog(List,List)
--Given a list L={L_1,...,L_n} of linear forms in 2 variables and a list m={m_1,...,m_n}
--of multiplicities, FPT2VarHomog(L,m) returns the F-pure threshold of the polynomial 
--L_1^(m_1)*...*L_n^(m_n). 
FPT2VarHomog (List,List) :=  opt -> (L,m) -> 
    FPT2VarHomogInternal(m,setFTData(gens ring L_0,L),MaxExp=>(opt.MaxExp),PrintCP=>(opt.PrintCP))


{*
    Miscellaneous.
*}

-- Some commands for dealing with vectors --

--canVector(i,n) returns the i-th canonical basis vector in dimension n
--Warning: for convenience, this uses Macaulay2's convention of indexing lists starting 
--with 0; so, for example, {1,0,0,0} is canVector(0,4), not canVector(1,4).
canVector = method()

canVector (ZZ,ZZ) := (i,n) -> 
(
    if ((i<0) or (i>=n)) 
        then (error "canVector(i,n) expects integers i and n with 0<=i<n.");   
    apply(n,j->if i==j then 1 else 0)
)
 
-- getNumAndDenom(u) takes a rational vector u and returns a pair (a,q), where a 
--is an integer vector and q an integer such that u=a/q.
getNumAndDenom = method()

getNumAndDenom (List) := u -> 
(
    den := lcm apply(u,n->denom n);
    a := apply(u,n->lift(n*den,ZZ));
    (a,den)        
)

--Computes the taxicab norm of a vector.
taxicabNorm = method()

taxicabNorm (List) := u -> sum(u,x->abs(x))

-- multOrder(a,b) finds the multiplicative order of a modulo b
multOrder = method()

multOrder (ZZ,ZZ) := (a,b) ->
(
    if gcd(a,b) != 1 then (error "Expected arguments to be relatively prime.");
    n := 1;
    x := a % b;
    while ( x != 1 ) do 
    (
	n = n+1;
	x = (x*a) % b
    );
    n	      
)     

-- Factorization of polynomials and splitting fields --

--factorList(F) factors the RingElement F and returns a list of pairs of the form
--{factor,multiplicity}.
factorList = method()

factorList (RingElement) := F ->  
(
    prod := factor F;
    apply(#prod, i -> {(prod#i)#0,(prod#i)#1}) 
)

--splittingField returns the splittingField of a polynomial over a finite field
splittingField = method()

splittingField (RingElement) := F -> 
(
    if not isPolynomialOverFiniteField(F) 
        then (error "splittingField expects a polynomial over a finite field");
    p:=char ring F;
    ord:=(coefficientRing(ring F))#order;
    factors:=first transpose factorList(F);
    deg:=lcm select(flatten apply(factors,degree),i->i>0);
    GF(p,deg*floor(log_p ord))
)

-- Some tests

--isBinaryForm(F) checks if F is a homogeneous polynomial in two variables.
--WARNING: what we are really testing is if the *ring* of F is a polynomial ring in two 
--variables, and not whether F explicitly involves two variables. (For example, if F=x+y 
--is an element of QQ[x,y,z], this test will return "false"; if G=x is an element of 
--QQ[x,y], this test will return "true".)
isBinaryForm = method()

isBinaryForm (RingElement) := F ->
(
    R:=ring F;
    isPolynomialRing(R) and numgens(R)==2 and isHomogeneous(F)
)

--isNonconstantBinaryForm(F) checks if F is a nonconstant homogeneous polynomial in two 
--variables. See warning under "isBinaryForm".
isNonConstantBinaryForm = method()

isNonConstantBinaryForm (RingElement) := F -> (isBinaryForm(F) and (degree(F))_0>0)

--isLinearBinaryForm(F) checks if F is a linearform in two variables. See warning 
--under "isBinaryForm".
isLinearBinaryForm = method()

isLinearBinaryForm (RingElement) := F -> (isBinaryForm(F) and (degree(F))_0==1)

--isPolynomialOverFiniteField(F) checks if F is a polynomial over a finite field.
isPolynomialOverFiniteField = method()

isPolynomialOverFiniteField (RingElement) := F ->
(
    R:=ring F;
    kk:=coefficientRing(R);
    try kk#order then (isPolynomialRing(R) and isField(kk))
    	else false   
)    

----------------------------------------------------------------
--************************************************************--
--Functions for computing eth roots                           --
--************************************************************--
----------------------------------------------------------------

ethRoot = method(); --- MK


--Computes I^{[1/p^e]}, we must be over a perfect field. and working with a polynomial ring
--This is a slightly stripped down function due to Moty Katzman, with some changes to avoid the
--use(Rm) which is commented out below
--The real meat of the function is ethRootInternal, this function just looks for a certain error and calls 
--the other function depending on that error.
ethRoot(Ideal,ZZ) := (Im,e) -> (
     J := Im;
     success := false;
     count := 0;
     try J = ethRootInternal(J,e) then success = true else (
--     print "blew a buffer";
	 while(count < e) do (	 	
	      J = ethRootInternal(J, 1);
	      count = count + 1
	 )
     );
     J
)

--This tries to compute (f^a*I)^{[1/p^e]} in such a way that we don't blow exponent buffers.  It can be much faster as well.
--We should probably just use it.  It relies on the fact that (f^(ap+b))^{[1/p^2]} = (f^a(f^b)^{[1/p]})^{[1/p]}.
ethRootSafe = (f, I, a, e) -> (
	R1 := ring I;
	p1 := char R1;
	
	aRem := a%(p1^e);
	aQuot := floor(a/p1^e);
	
	expOfA := basePExpMaxE(aRem,p1,e); --this gives "a base p", with the left-most term the smallest "endian".
	
	IN1 := I;
	
	if (e > 0) then (
		IN1 = IN1*ideal(f^(expOfA#0));
		IN1 = ethRoot(IN1, 1);
		i := 1;
	
		while(i < #expOfA) do (
			IN1 = ethRoot( IN1*ideal(f^(expOfA#i)), 1);
			i = i + 1;
		)
	);
	IN1*ideal(f^(aQuot))
)

--This tries to compute (f1^a1*f2^a2*...fk^ak*I)^{[1/p^e]} in such a way that we don't blow exponent buffers.  It can be much faster as well.
ethRootSafeList = (elmtList, I1, aList, e1) -> (
	   R1 := ring I1;
        p1 := char R1;
        
        aListRem := apply(aList, z1 -> z1%(p1^e1) );
        aListQuot := apply(aList, z1 -> floor(z1/p1^e1) );
        
        expOfaList := apply(aListRem, z1-> basePExpMaxE(z1, p1, e1) );
        
        aPowerList := apply(elmtList, expOfaList, (f1, z1) -> f1^(z1#0));
        
        IN1 := I1*ideal(fold(times, aPowerList));
        if (e1 > 0) then (
                IN1 = ethRoot(IN1, 1);
                i := 1;
                while(i < e1) do (
                        aPowerList = apply(elmtList, expOfaList, (f1, z1) -> f1^(z1#i));
                        IN1 = ethRoot( IN1*ideal(fold(times, aPowerList)), 1);
                        i = i + 1;
                )
        );
        aPowerList = apply(elmtList, aListQuot, (f1, z1) -> f1^z1);
        IN1*ideal(fold(times, aPowerList))
)

ethRoot(RingElement, Ideal, ZZ, ZZ) := (f, I, a, e) -> ethRootSafe (f, I, a, e) ---MK

ethRootInternalOld = (Im,e) -> (
     if (isIdeal(Im) != true) then (
     	  error "ethRoot: Expted a nonnegative integer."; 
     );
     if (not (e >= 0)) then (error "ethRoot: Expected a nonnegative integer.");
     Rm:=ring(Im); --Ambient ring
     if (not (class Rm === PolynomialRing)) then (error "ethRoot: Expected an ideal in a PolynomialRing.");
     pp:=char(Rm); --characteristic
     Sm:=coefficientRing(Rm); --base field
     n:=rank source vars(Rm); --number of variables
     vv:=first entries vars(Rm); --the variables
     YY:=local YY; -- this is an attempt to avoid the ring overwriting
                         -- the ring in the users terminal
			 -- MonomialOrder=>ProductOrder{n,n}
     myMon := monoid[ (vv | toList(YY_1..YY_n)), MonomialOrder=>ProductOrder{n,n},MonomialSize=>64];
     R1:=Sm myMon; -- a new ring with new variables
     vv2 := first entries vars R1;
     J0:=apply(1..n, i->vv2#(n+i-1)-vv2#(i-1)^(pp^e)); -- 
     --print J0;
     M:=toList apply(1..n, i->vv2#(n+i-1)=>substitute(vv#(i-1),R1));

     G:=first entries compress( (gens substitute(Im,R1))%gens(ideal(J0)) );

     L:=ideal 0_R1;
     apply(G, t-> --this appears to just be a for loop
	  {
    	       L=L+ideal((coefficients(t,Variables=>vv))#1);
	  });
     L2:=mingens L;
     L3:=first entries L2;
     L4:=apply(L3, t->substitute(t,M));
     --use(Rm);
     substitute(ideal L4,Rm)
)

ethRootInternal = (I,e) -> (
     if (not isIdeal(I)) then (error "ethRoot: Expected first argument to be an ideal.");
     if (not e >= 0) then (error "ethRoot: Expected second argument to be a nonnegative integer.");
     R:=ring(I); --Ambient ring
     if (class R =!= PolynomialRing) then (error "ethRoot: Expected an ideal in a PolynomialRing.");
     p:=char(R); --characteristic
     kk:=coefficientRing(R); --base field
     if ((kk =!= ZZ/p) and (class(kk) =!= GaloisField)) then (error "ethRoot: Expected the coefficient field to be ZZ/p or a GaloisField.");
     n:=rank source vars(R); --number of variables
     var:=first entries vars(R); --the variables (henceforth denoted X_i)
     Y:=local Y;
     S:=kk(monoid[(var | toList(Y_1..Y_n)), MonomialOrder=>ProductOrder{n,n},MonomialSize=>64]);
         -- brand new ring, with a variable Y_i for each X_i
     newvar := first entries vars(S);
     J:=matrix {toList apply(0..(n-1), i->newvar#(n+i)-newvar#(i)^(p^e))}; 
         -- J = (Y_i-x_i^(p^e)) 
     rules:=toList apply(0..(n-1), i->newvar#(n+i)=>substitute(var#(i),S)); 
         -- {Y_i =>X_i} 
     G:=first entries compress((gens substitute(I,S)) % J);
     	 -- replaces X_i^(p^e) with Y_i 
     L:=sum(G,t->ideal((coefficients(t,Variables=>var))#1));	 
     L=first entries mingens L;
     L=apply(L, t->substitute(t,rules));
     q:=kk#order;
     if (q > p) then 
     (
	 a:=(gens kk)#0;
	 e0:=floor(log(p,q)); 
	 root:=a^(p^(e0-(e%e0)));
     	 L=apply(L,t->substitute(t,a=>root));
	     -- substitute generator of kk with its p^e-th root
     );
     substitute(ideal L,R)
)


--A short version of ethRoot
eR = (I1,e1)-> (ethRoot(I1,e1) )

---------------------------------------------------------------------------------------
--- The following code was written in order to more quickly compute eth roots of (f^n*I)
--- It is used in fancyEthRoot
----------------------------------------------------------------------------------------
--- Find all ORDERED partitions of n with k parts
allPartitions = (n,k)->
(
	PP0:=matrix{ toList(1..k) };
	PP:=mutableMatrix PP0;
	allPartitionsInnards (n,k,PP,{})
)

allPartitionsInnards = (n,k,PP,answer)->
(
	local i;
	if (k==1) then 
	{
		PP_(0,k-1)=n;
		answer=append(answer,first entries (PP));
	}
	else
	{
		for i from 1 to n-(k-1) do
		{
			PP_(0,k-1)=i;
			answer=allPartitionsInnards (n-i,k-1,PP,answer)	;	
		};
	};
	answer
)




--- write n=a*p^e+a_{e-1} p^{e-1} + \dots + a_0 where 0\leq e_j <p 
baseP1 = (n,p,e)->
(
	a:=n//(p^e);
	answer:=1:a;
	m:=n-a*(p^e);
	f:=e-1; 
	while (f>=0) do
	{
		d:=m//(p^f);
		answer=append(answer,d);
		m=m-d*(p^f);
		f=f-1;
	};
	answer
)	


fancyEthRoot = (I,m,e) ->
(
	G:=first entries mingens I;
	k:=#G;
	P:=allPartitions(m,k);
	R:=ring(I);
	p:=char(R);
	answer:=ideal(0_R);
	apply(P, u->
	{
	---print("Partition: ",u);
		a:=ideal(1_R);
		U:=apply(u, v->baseP1(v,p,e));
		for i from 0 to e do
		{
			j:=e-i;
			g:=1_R;
			for l from 0 to k-1 do g=g*(G#l)^((U#l)#j); 
			a=ideal(g)*a;
			if (i<e) then a=ethRoot(a ,1);
---print(g,answer);
		};
		answer=answer+a;
	});
	ideal(mingens(answer))
)

ethRoot (Ideal, ZZ, ZZ) := (I,m,e) -> fancyEthRoot (I,m,e)  --- MK



--Computes I^{[1/p^e]}, we must be over a perfect field. and working with a polynomial ring
--This is a slightly stripped down function due to Moty Katzman, with some changes to avoid the
--use(Rm) which is commented out below
--The real meat of the function is ethRootInternal, this function just looks for a certain error and calls 
--the other function depending on that error.
ethRoot(Ideal,ZZ) := (Im,e) -> (
     J := Im;
     success := false;
     count := 0;
     try J = ethRootInternal(J,e) then success = true else (
--     print "blew a buffer";
	 while(count < e) do (	 	
	      J = ethRootInternal(J, 1);
	      count = count + 1
	 )
     );
     J
)

----------------------------------------------------------------
--************************************************************--
--Functions for computing compatibly split ideals             --
--************************************************************--
----------------------------------------------------------------

-----------------------------------------------------------------------


--- Start of MK ---------------------------------------------------------------------------------------------------

-- FIND IDEALS COMPATIBLE WITH A GIVEN NEAR-SPLITTING
-- This is an implementation of the algorithm described in
-- Moty Katzman and Karl Schwede's paper 
-- "An algorithm for computing compatibly Frobenius split subvarieties"
-- J. Symbolic Comput. 47 (2012), no. 8, 996–1008. 

----------------------------------------------------------------------------------------


--- Input:
---   	an element u of the polynomial ring R OVER A PRIME FIELD.
--- Output:
---	A list of all prime ideals P such that
---	(a) u P \subseteq P^{[p]}, and
---	(b) the action of uT on the the annihilator of P on the injective hull of the residue field of R 
---	is not the zero Frobenius map.


findAllCompatibleIdeals = (u) ->(
	L:={}; R:=ring u; p:=char R;
	P:=ideal(0_R);
	J:=ethRoot(ideal(u),1);
	t:=1_R % (gens J);
	if (t != 0_R) then print("*** WARNING *** Frobenius action has nilpotent elements");
	findAllCompatibleIdealsInnards (u,L,P)
)



findAllCompatibleIdealsInnards = (u,L,P) ->(
	R:=ring u;
	p:=char R;
	local tau;
	local Plist;
	P1:=frobeniusPower(P,1);
	C1:=ideal((singularLocus(P)).relations);
	---tau=ideal mingens star(C1,u,1) ; ---OLD VERSION
	tau=ideal mingens ascendIdeal (C1, u, 1);
	Plist=minimalPrimes tau;
	local Q;
	local T;
	apply(Plist, Q->
	{
		f:= any(L,T -> T == Q);
---print(L,Q,f);
		if (not f) then
		{
			L=append(L,Q);
			L=unique(L | findAllCompatibleIdealsInnards(u,L,Q));
		};
	});
---
	C2:=(P1+ideal(u)):(P1:P);
---	JB:=C1*C2; ---MK
---print(mingens P, mingens JB);
---tau=ideal mingens star(C2,u,1) ;  --- OLD VERSION
	tau=ideal mingens ascendIdeal  (C2, u, 1);
	Plist=minimalPrimes tau;
	local Q;
	local T;
	apply(Plist, Q->
	{
		f:= any(L,T -> T == Q);
	---print(L,Q,f);
		if (not f) then
		{
			L=append(L,Q);
			L=unique(L | findAllCompatibleIdealsInnards(u,L,Q));
		};
	});
	---
	L
)



-----------------------------------------------------------------------------
--- Extend the Frobenius p^e th roots and star operations to submodules of
--- free modules (over polynomial rings with *prime* coeeficient field)
--- This implements the methods described in 
--- Moty Katzman and Wenliang Zhang's paper
--- "Annihilators of Artinian modules compatible with a Frobenius map"
--- Journal of Symbolic computation, 2014

-----------------------------------------------------------------------------


mEthRootOfOneElement = (v,e) ->(
	local i;
	local d;
	local w;
	local m;
	R:=ring(v); p:=char R;
	F:=coefficientRing(R);
	n:=rank source vars(R);
	V:=ideal vars(R);
	vv:=first entries vars(R);
	YY:=local YY;
	R1:=F[vv, YY_1..YY_n, MonomialOrder=>ProductOrder{n,n},MonomialSize=>16];
	V=substitute(V,R1);
	---------------------------
	M0:=set {1_R1};
	apply(vv, w->
	{
		ww:=substitute(w,R1);
		M1:=set toList apply(0..p-1, i-> ww^i);
		M0=M0**M1
	});
	M:=toList apply(elements(M0), w-> product toList deepSplice(w));
---------------------------
	J0:=gens ideal apply(1..n, i->YY_i-substitute(vv#(i-1)^(p^e),R1));
	S:=toList apply(1..n, i->YY_i=>substitute(vv#(i-1),R1));
	Ie:=transpose matrix{{(rank target v):0_R1}}; 
	ev:=entries substitute(v,R1);
	apply(M, m->
	{
		L:={};
		apply(ev, t->
		{
			tt:=((t#0)%J0);
			q1:=coefficients( tt , Variables=>(first entries gens V), Monomials=>{m});
			q2:=q1#1;
			q3:=first first entries q2;
			q3=substitute(q3,S);
			L=append(L,q3);
	---		print(m,tt,q3);
		});
	---	print(ev,L,m);
		Ie=Ie | (transpose matrix {L});
	});
---	use R;
	compress(substitute(Ie,R))
)




mEthRoot = (A,e) ->(
	local i;
	local answer;
	answer1:=apply(1..(rank source A), i->mEthRootOfOneElement (A_{i-1},e));
	if (#answer1==0) then 
	{
		answer=A;
	}	
	else
	{
		answer=answer1#0;
		apply(2..(#answer1), i->answer=answer | answer1#(i-1));
	};
	mingens( image answer )
)	


ethRoot (Matrix, ZZ) := (A,e) -> mEthRoot (A,e)  --- MK






--- Mstar is the implementaion of the star closure operation desribed in 
--- M Katzman's "Parameter test ideals of Cohen Macaulay rings" 
--- Input:
---    ideals I and element u of the same polynomial ring R OVER A PRIME FIELD.
---    a positive integer e
---    a prime p which is the characteristic of R
--- Output:
---    the smallest ideal J of R containing I with the property that u^(1+p+...+p^(e-1)) J is in J^{[p^e]}
Mstar = (A,U,e) ->(
	local answer;
	R:=ring(A); p:=char R;
	if (A==0) then
	{
		answer=A;
	}
	else
	{
		f:=true;
		Ne:=sum toList(apply(0..(e-1), i->p^i));
		lastA:= A;
		while (f) do
		{
			f=false;
			A1:=mEthRoot(mingens image ((U^Ne)*lastA),e);
			A1=A1 | lastA;
			t1:=compress ((A1))%((lastA));
			if (t1!=0) then 
			{
				f=true;
				lastA=mingens image A1;
			};
		};
		answer=mingens (image A1);
	};
	answer
)


--- end of MK ---------------------------------------------------------------------------------------------------


----------------------------------------------------------------
--************************************************************--
--Functions for computing test ideals, and related objects.   --
--************************************************************--
----------------------------------------------------------------


--Finds the smallest phi-stable ideal containing the given ideal Jk
--in a polynomial ring Sk
--Jk is the given ideal, ek is the power of Frobenius to use, hk is the function to multiply 
--trace by to give phi:  phi(_) = Tr^(ek)(hk._)
--This is based on ideas of Moty Katzman, and his star closure
ascendIdeal = (Jk, hk, ek) -> (
     Sk := ring Jk;
     pp := char Sk;
     IN := Jk;
     IP := ideal(0_Sk);
     --we want to make the largest ideal that is phi-stable, following Moty Katzman's idea
     --we do the following
     while (isSubset(IN, IP) == false) do(
     	  IP = IN;
--	  error "help";
	  IN = eR(ideal(hk)*IP, ek)+IP
     );

     --trim the output
     trim IP
)

--Works like ascendIdeal but tries to minimize the exponents elements are taken to
ascendIdealSafe = (Jk, hk, ak, ek) -> (
	Sk := ring Jk;
     pp := char Sk;
     IN := Jk;
     IP := ideal(0_Sk);
     --we want to make the largest ideal that is phi-stable, following Moty Katzman's idea
     --we do the following
     while (isSubset(IN, IP) == false) do(
     	  IP = IN;
--	  error "help";
	  	IN = ethRootSafe(hk, IP, ak, ek)+IP
     );

     --trim the output
     trim IP
)




--works just like ascendIdealSafe but also handles lists of hk to powers...
ascendIdealSafeList ={AscentCount=>false} >> o ->  (Jk, hkList, akList, ek) -> (
	Sk := ring Jk;
	pp := char Sk;
	IN := Jk;
	IP := ideal(0_Sk);
	i1 := 0;
	--we ascend the ideal as above
	while (isSubset(IN, IP) == false) do(
		i1 = i1 + 1; --
		IP = IN;
		IN = ethRootSafeList( hkList, IP, akList, ek) + IP
	);
	
	--trim the output
	if (o.AscentCount == false) then 
		trim IP
	else (trim IP, i1)
)

--MKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMK
-- minimalCompatible is a method which is implemented as:
-- (1) the finding of the smallest ideal J which satisfies uJ\subset J^{[p^e]} 
---    containg a given ideal for a given ring element u,
-- (2) the finding of the smallest submodule V of a free module which satisfies UV\subset V^{[p^e]} 
--     containg a given submodule for a given matrix U.
minimalCompatible = method();
minimalCompatible(Ideal,RingElement,ZZ) :=  (Jk, hk, ek) -> ascendIdeal (Jk, hk, ek)
minimalCompatible(Ideal,RingElement,ZZ,ZZ) :=  (Jk, hk, ak, ek) -> ascendIdealSafe (Jk, hk, ak, ek)
minimalCompatible(Matrix,Matrix,ZZ) := (A,U,e) -> Mstar (A,U,e)

--MKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMK


--Finds a test element of a ring R = k[x, y, ...]/I (or at least an ideal 
--containing a nonzero test element).  It views it as an element of the ambient ring
--of R.  It returns an ideal with some of these elements in it.
--One could make this faster by not computing the entire Jacobian / singular locus
--instead, if we just find one element of the Jacobian not in I, then that would also work
--and perhaps be substantially faster
findTestElementAmbient = Rk -> (
     --Sk := ambient Rk;
     -- Ik := ideal Sk;
     
     Jk := ideal singularLocus(Rk);
     if (isSubset(Jk, ideal Rk) == true) then 
          error "findTestElementAmbient: No test elements found, is the ring non-reduced?";
	  
     
     Jk          
)


--Outputs the test ideal of a (Q-)Gorenstein ring (with no pair or exponent)
--ek is the number such that the index divides (p^ek - 1)
--It actually spits out the appropriate stable/fixed ideal inside the ambient ring
tauQGorAmb = (Rk, ek) -> (
     Jk := findTestElementAmbient(Rk);
     hk := findQGorGen(Rk, ek);

     sub(ascendIdeal(Jk,hk,ek),Rk)
)

--Computes the test ideal of an ambient Gorenstein ring
tauGorAmb = (Rk) -> (tauQGorAmb(Rk, 1))

--Computes the test ideal of (R, f^(a/(p^e - 1)))
--when R is a polynomial ring.  This is based upon ideas of Moty Katzman.
tauAOverPEMinus1Poly = (fm, a1, e1) -> (
     Rm := ring fm;
     pp := char Rm;
     a2 := a1 % (pp^e1 - 1);
     k2 := a1 // (pp^e1 - 1); --it seems faster to use the fact that tau(f^(1+k)) = f*tau(f^k) 
     --this should be placed inside a try, and then if it fails we should be smarter...
     --fpow := fastExp(fm,a2);
     --IN := eR(ideal(fpow*fm),e1);  --the idea contained inside the test ideal.
     IN := ethRootSafe(fm, ideal(fm), a2, e1);
     
     IN = ascendIdealSafe(IN, fm, a2, e1);
     -- this is going to be the new value.  The *fm is a test element
     --5/0;
     --return the final ideal
     IN*ideal(fm^k2)
)

--Computes the test ideal of (R, f^t) when R 
--is a polynomial ring over a perfect field.
tauPolyOld = (fm, t1) -> (
     Rm := ring fm; 
     pp := char Rm;
     L1 := divideFraction(t1,pp); --this breaks up t1 into the pieces we need
     local I1;
     --first we compute tau(fm^{a/(p^c-1)})
     if (L1#2 != 0) then 
          I1 = tauAOverPEMinus1Poly(fm,L1#0,L1#2) else I1 = ideal(fm^(L1#0));     
	  
       
     
     --now we compute the test ideal using the fact that 
     --tau(fm^t)^{[1/p^a]} = tau(fm^(t/p^a))
     if (L1#1 != 0) then 
          ethRoot(I1, L1#1) else I1
)

--a slightly faster tauPoly
tauPoly = (fm, t1) -> (
     Rm := ring fm; 
     pp := char Rm;
     L1 := divideFraction(t1,pp); --this breaks up t1 into the pieces we need
     local I1;
     --first we compute tau(fm^{a/(p^c-1)})
     if (L1#2 != 0) then (
     	I1 = tauAOverPEMinus1Poly(fm,L1#0,L1#2);
     	if (L1#1 != 0) then
     		I1 = ethRoot(I1, L1#1)
     	)
     else (
     	if (L1#1 != 0) then
     		I1 = ethRootSafe(fm, ideal( sub(1, Rm)), L1#0, L1#1 )
     	else
 	    		I1 = ideal(fm^(L1#0))
 	    	);
     I1
)

--This is an internal function
--It is used to compute the test ideals of pairs (R, fm^(a1/p^e1-1)) where
--R = Sk/Ik.
--Inputs are Jk, a nonzero ideal contained in the test ideal
--hk, the multiple used to form phi of the ambient ring.  ek is the power associated with hk
--a1 and e1 and fm are as above
tauAOverPEMinus1QGorAmbOld = (Sk, Jk, hk, ek, fm, a1, e1) -> (
     pp := char Sk;
     et := lcm(ek, e1);
     hk1 := (hk)^(numerator ((pp^et - 1)/(pp^ek - 1)));  
       --hk for higher powers are simply appropriate powers of hk for lower powers, 
       --may as well take advantage of that
     a3 := numerator (a1*(pp^et - 1)/(pp^e1 - 1)); --we need to use a common e for both the 
                                               --index of R and of our divisor.
     
     a2 := a3 % (pp^et - 1);
     k2 := a3 // (pp^et - 1); --it seems faster to use the fact 
                              --that tau(f^(1+k)) = f*tau(f^k) 
     fpow := fm^a2; 
     
     Iasc := ascendIdeal(Jk*ideal(fm), fpow*hk1, et);
    
     Iasc*ideal(fm^k2)
)

tauAOverPEMinus1QGorAmb = (Sk, Jk, hk, ek, fm, a1, e1) -> (
     pp := char Sk;
     et := lcm(ek, e1);
     
     ak1 := numerator ((pp^et - 1)/(pp^ek - 1)); --an exponent for hk
     a3 := numerator (a1*(pp^et - 1)/(pp^e1 - 1)); --we need to use a common e for both the 
                                               --index of R and of our divisor.
                                               
	a2 := a3 % (pp^et - 1);
     k2 := a3 // (pp^et - 1); --it seems faster to use the fact that we can do simple Skoda for tau
     
     Jl := ascendIdealSafe(Jk, hk, 1, ek);
                      
        --          assert false;                             
     Iasc := ascendIdealSafeList(Jk*ideal(fm)^(ceiling(a3/(pp^et - 1))), (fm, hk), (a2, numerator ((pp^et - 1)/(pp^ek - 1))), et);
     
--     assert false;
     
     Iasc*ideal(fm^k2)
)


--Computes the test ideal of (Rk, fk^t1).  Here we assume the index of the canonical
--divides (p^ek - 1)
tauQGor = (Rk, ek, fk, t1) -> (
     Sk := ambient Rk;
     pp := char Sk;
     L1 := divideFraction(t1,pp); --this breaks up t1 into the pieces we need
     hk := findQGorGen(Rk, ek); --the term in the test ideal
     Jk := findTestElementAmbient(Rk); --this finds some test elements (lifted on the ambient
                                       --ring).  Right now it is slow because of a call to 
				       --singularLocus (ie, computing a Jacobian).
     I1 := ideal(0_Sk); I2 := ideal(0_Sk);
     fm := lift(fk, Sk); --we lift our f to the ambient polynomial ring
     a1 := L1#0; e1 := L1#2; pPow := L1#1; --t1 = a1 / (pp^pPow * (pp^e1 - 1))
     
     --before continuing, we need to rewrite a1 / (pp^pPow * (pp^e1 - 1)) as 
     --                                      a3 / (pp^(n1*ek) * (pp^e1 - 1))
     --the point is that ek needs to divide pPow
     remain := pPow % ek;
     dualRemain := ek - remain;
     
     pPow = pPow + dualRemain; --note in the end here, ek divides pPow
     a3 := a1*pp^(dualRemain);
     
     if (e1 != 0) then assert (t1 == a3/((pp^e1-1)*pp^pPow) ) else assert (t1 == a3/(pp^pPow) );
     
     d1 := pp^(pPow); if (e1 != 0) then d1 = d1*(pp^e1-1); --this is our denominator, used
                                                           --for simplifying computations
     a2 := a3 % d1;
     k2 := a3 // d1; --it seems faster to use the fact 
                              --that tau(f^(k+t)) = f^k*tau(f^t).  We'll toss on the multiple 
			      --f^k at the end
	     			  
     local I2;
     --first we compute tau(fk^{a2/(p^e1-1)})
     if (e1 != 0) then (
          I1 = tauAOverPEMinus1QGorAmb(Sk,Jk,hk,ek,fm,a2,e1);
          if (pPow != 0) then (
          	I2 = ethRootSafe(hk, I1, numerator((pp^pPow - 1)/(pp^ek - 1)), pPow)
		)
		else I2 = I1
     )
     else (
	  	I1 = ascendIdeal(Jk, hk, ek);
	  	if (pPow != 0) then (
	  		I2 = ethRootSafeList( (hk, fm), I1, (numerator((pp^pPow - 1)/(pp^ek - 1)), a2), pPow)
	  	)
	  	else I2 = I1
     );

	  
     sub(I2, Rk)*ideal(fk^k2)
)

--Computes tau(Rk,fk^tk), assuming Gorenstein rings
tauGor = (Rg,fg,tg) -> tauQGor (Rg,1,fg,tg)


----------------------------------------------------------------
--************************************************************--
--Test ideals for non-principal ideals                        --
--************************************************************--
----------------------------------------------------------------

flattenedReesAlgebra = (I1) -> (--takes an ideal, forms the rees algebra, and returns the rees algebra in two ways, first with flattened variables and the second without
	S1 := reesAlgebra I1;
	J1 := ideal S1;
	tempMonoid := S1.FlatMonoid;
	k1 := coefficientRing (ring I1);
	S2 := k1 tempMonoid;
	
	J2 := sub(J1, S2);
	
	(S2/J2, S1)
)

needsPackage "BGG"; --we'll be pushing forward...

needsPackage "Divisor";

tauNonPrincipalAOverPEPoly = {Verbose=> false}>> o -> (I1, a1, e1) -> ( -- computes \tau(I^{a/p^e}) for I an ideal in a polynomial ring
	if ( not(codim(I1) > 1)) then error "We can only handle ideals of codimension > 1 at this time.";
	
	--this function currently doesn't take advantage of Skoda's theorem, this will need to be done

	reesList := flattenedReesAlgebra I1;
	A1 := reesList#0; --this one has flattened variables
	A2 := reesList#1;
 	irrIdeal := sub(ideal(first entries vars A1), A1);
 	singLocus := ideal singularLocus (A1);
 	
 	IRees := sub(I1, A2);
 	
 	canList := canonicalIdeal(A1, FullMap=>true);
 	canIdeal := canList#0;
 	canMap := canList#1;
 	
 	paraTest := paraTestModuleAmbient(A1, canIdeal); 
 		
 	newMap := map(A1^1/(paraTest#0), canList#2, matrix(canMap));
 	newKer := (ker newMap)**A2; --this is the parameter test submodule of the canonical module  

	flag := false;
	i1 := e1;
	R1 := ring I1;
	p1 := char R1;
	ascend := I1; --dummy variables for checking whether we are done
	descend := ideal(sub(1, R1)); --dummy variables for checking whether we are done
	
	while (flag == false) do (
		ascend = fancyEthRoot(I1, a1*p1^(i1-e1), i1);
		if (o.Verbose == true) then (print  "Ascending ideal"; print ascend);
		
		flag = isSubset(descend, ascend);
		if (o.Verbose == true) then (print "flag"; print flag);
		if (flag == false) then (
			
			myDirectImage := HH_0(directImageComplex(IRees^(a1*p1^(i1-e1))*newKer, Regularity=>(10+a1))); 	
 	
		 	directIdeal := module2Ideal(myDirectImage, R1);
 			if ( codim(directIdeal)==1) then error "This function produced a codimension 1 ideal.";
 	
 			descend = ethRoot(directIdeal, i1);
 			if (o.Verbose == true) then (print  "Descending ideal"; print descend)
		);
		

		flag = isSubset(descend, ascend);
		
		--the following should be removed eventually, it is only here for debug purposes
		if ((flag == true) and (isSubset(ascend, descend)==false)) then error "Major error detected";
		i1 = i1+1;
		if (o.Verbose==true) then (print "Loop complete, continue?"; print (not flag) );
	);
	
	ascend
)

----------------------------------------------------------------
--************************************************************--
--Functions for computing sigma                               --
--************************************************************--
----------------------------------------------------------------


--Computes Non-Sharply-F-Pure ideals over polynomial rings for (R, fm^{a/(p^{e1}-1)}), 
--at least defined as in Fujino-Schwede-Takagi.
sigmaAOverPEMinus1Poly ={HSL=> false}>> o -> (fm, a1, e1) -> ( 
     Rm := ring fm;
     pp := char Rm;
     m1 := 0;
	e2 := e1;
	a2 := a1;
	--if e1 = 0, we treat (p^e-1) as 1.  
     if (e2 == 0) then (e2 = 1; a2 = a1*(pp-1));
     if (a2 > pp^e2-1) then (m1 = floor((a2-1)/(pp^e2-1)); a2=((a2-1)%(pp^e2-1)) + 1 );
     --fpow := fm^a2;
     IN := eR(ideal(1_Rm),e2); -- this is going to be the new value.
     -- the previous commands should use the fast power raising when Emily finishes it
     IP := ideal(0_Rm); -- this is going to be the old value.
     count := 0;
     
     --our initial value is something containing sigma.  This stops after finitely many steps.  
     while (IN != IP) do(
		IP = IN;
	  	IN = ethRootSafe(fm,IP,a2,e2); -- eR(ideal(fpow)*IP,e2);
	  	count = count + 1
     );

     --return the final ideal and the HSL number of this function
     if (o.HSL == true) then {IP*ideal(fm^m1),count} else IP*ideal(fm^m1)
)

--Computes Non-Sharply-F-pure ideals for non-polynomial rings with respect to no pair.
sigmaQGor ={HSL=> false}>> o -> (Rm, gg) -> (
     Sm := ambient Rm; --the polynomial ring that Rm is a quotient of
     hk := findQGorGen(Rm, gg);
     
     IN := ideal(1_Sm);
     count := 0;
     IP := ideal(0_Sm);
     
     while (IN != IP) do(
     	IP = IN;
     	IN = eR(ideal(hk)*IP,gg);
     	count = count + 1
     );
     
     --return the ideal and HSL
     if (o.HSL == true) then {sub(IP,Rm), count} else sub(IP, Rm)
)

--Computes Non-Sharply-F-Pure ideals for non-polynomial rings
--gg is the Gorenstein index
sigmaAOverPEMinus1QGor  ={HSL=> false}>> o -> (fk, a1, e1, gg) -> (
     Rm := ring fk;
     Sm := ambient Rm; --the polynomial ring that Rm is a quotient of
     pp := char Rm;
     ek := lcm(e1,gg); --we need to raise things to common powers
     hk := findQGorGen(Rm, gg); --it will be faster to find the Q-Gorenstein generator
     hl := hk^(sub((pp^ek - 1)/(pp^gg - 1), ZZ) ); --
	fm := lift(fk, Sm); --lift fk to the ambient ring
	fpow := fm^(a1*sub( (pp^ek - 1)/(pp^e1-1), ZZ) );


	IN := sigmaAOverPEMinus1Poly(hk,1,gg);
	count := 0;
	IP := ideal(0_Sm);

	while (IN != IP) do(
		IP = IN;
		IN = eR(ideal(fpow*hl)*IP, e1);
		count = count + 1
	);
	
     --return the final ideal
     if (o.HSL == true) then {sub(IP,Rm), count} else sub(IP,Rm)
	
)

----------------------------------------------------------------
--************************************************************--
--Functions for computing parameter test modules and ideals   --
--************************************************************--
----------------------------------------------------------------


--This function computes the parameter test module of a ring, it returns it as a submodule of a canonical ideal.
--this is a slightly modified function originally written by Moty Katzman for "Parameter test ideals of Cohen Macaulay rings"
--it returns the lift of the canonical module to the ambient ring

canonicalIdeal ={FullMap=> false} >> o -> (R1) -> (
	S1 := ambient R1;
	I1 := ideal(R1);
	d1 := (dim S1) - (dim R1);
	local answer2;
	
	degShift := sum degrees S1;
	myExt := prune( Ext^d1(S1^1/I1, S1^{-degShift}));
	canModuleMatrix := relations(myExt);
	
	answer:=0;
	s1:=syz transpose substitute(canModuleMatrix,R1);
	s2:=entries transpose s1;
	--use S1;
	apply(s2, t->
	{
		s3:=substitute(syz gens ideal t,S1);
---		print(s3%canModuleMatrix);
		if ((s3%canModuleMatrix)==0) then
		{
			answer2 = t;
			answer=substitute(mingens ideal t,S1);
			break;
		};
	});
	
	
	
	if (o.FullMap == true) then (ideal answer, map(R1^1, myExt**R1, matrix {answer2}), (myExt**S1^{-degShift})**R1) else ideal answer
)

--moduleToIdeal = (M1, R1) -> (--turns a module to an ideal of a ring, it returns the lift of the ideal to the ambient ring
--	S1 := ambient R1;
---	myMatrix := substitute(relations prune M1, S1);
--	
--	answer:=0;
--	s1:=syz transpose substitute(myMatrix,R1);
--	s2:=entries transpose s1;
--	
--	apply(s2, t->
--	{
--		s3:=substitute(syz gens ideal t,S1);
---		print(s3%canModuleMatrix);
--		if ((s3%myMatrix)==0) then
--		{
--			answer=substitute(mingens ideal t,S1);
--			break;
--		};
--	});
--	ideal answer	
--)

--the following function computes the u of a canonical ideal in a polynomial ring
--it uses previous work of Katzman
finduOfIdeal = (canIdeal, defIdeal) -> (
	Ip := frobeniusPower(defIdeal, 1);
	tempIdeal := intersect( (frobeniusPower(canIdeal, 1)) : canIdeal, Ip : defIdeal );
	
	M1 := compress ((gens tempIdeal)%(gens Ip));
	first first entries M1
)

--computes the parameter test submodule of a given ring.  It outputs the parameter test module (as an ideal), it then outputs the canonical module (as an ideal), and finally it outputs the term u used as the action on the ideal
paraTestModuleAmbient = method();

paraTestModuleAmbient (Ring) := (R1) -> (
	S1 := ambient R1;
	I1 := ideal(R1);
	
	canIdeal := canonicalIdeal(R1);
	
	J1 := findTestElementAmbient(R1);
	tau0 := J1*canIdeal; --this is the starting test element times the ideal
	
	u1 := finduOfIdeal(canIdeal, I1); --this is the multiplying object that gives us (u*omega)^{[1/p]} \subseteq omega.
	
	tauOut := ascendIdeal(tau0, u1, 1);
	
	(sub(tauOut, R1), sub(canIdeal, R1), u1)
)

paraTestModuleAmbient (Ring, Ideal) := (R1, canIdeal) -> (
	S1 := ambient R1;
	I1 := ideal(R1);
	
	J1 := findTestElementAmbient(R1);
	tau0 := J1*canIdeal; --this is the starting test element times the ideal
	
	u1 := finduOfIdeal(canIdeal, I1); --this is the multiplying object that gives us (u*omega)^{[1/p]} \subseteq omega.
	
	tauOut := ascendIdeal(tau0, u1, 1);
	
	(sub(tauOut, R1), sub(canIdeal, R1), u1)
)

--computes the parameter test ideal of an ambient ring
paraTestIdealAmbient = (R1) -> (
	tempList := paraTestModuleAmbient(R1);
	(tempList#0) : (tempList#1)
)

--this computes the parameter test module \tau(R, f^t).  It does not assume that R is a polynomial ring.
paraTestModule ={AscentCount=>false} >> o -> (fk, t1) -> ( --maintained by Karl
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
	local tempList;
	ascendingCount := 0;
--	assert false;
	if (cc != 0) then	
		if (o.AscentCount == false) then (firstTau = ascendIdealSafeList( J1*ideal(f1^(pp^bb*ceiling(t1))), (f1, u1), (aa, uPower), cc))
		else (tempList = ascendIdealSafeList( J1*ideal(f1^(pp^bb*ceiling(t1))), (f1, u1), (aa, uPower), cc, AscentCount=>true);
			firstTau = tempList#0;
			ascendingCount = tempList#1;
		)
--		firstTau = ascendIdeal(J1*ideal(f1^(aa)), f1^aa*u1^(uPower), cc)
		--I should write an ascendIdealSafe that works for multiple elements raised to powers...	
	else 
--		firstTau = ascendIdeal(J1, u1^(uPower), 1)*ideal(f1^aa);
		firstTau = ascendIdealSafe(J1, u1, uPower, 1);
			
	secondTau := firstTau;
	if (bb != 0) then
		secondTau = ethRootSafe(u1, firstTau, floor((pp^bb-1)/(pp-1)) , bb);

	if (o.AscentCount == false) then (sub(secondTau, R1), omegaAmb, u1) else (sub(secondTau, R1), omegaAmb, u1, ascendingCount)
)



----------------------------------------------------------------
--************************************************************--
--Functions for checking whether a ring/pair is F-pure/regular--
--************************************************************--
----------------------------------------------------------------

-- Given an ideal I of polynomial ring R
-- this uses Fedder's Criterion to check if R/I is F-pure
-- Recall that this involves checking if I^[p]:I is in m^[p]
-- Note:  We first check if I is generated by a regular sequence.

isFPure = I1->(
    maxIdeal:= monomialIdeal(first entries vars ring I1);
    local answer;
    local cond;
    p1:=char ring I1;
    if codim(I1)==numgens(I1) then(
	L:=flatten entries gens I1;
	cond = isSubset(ideal(product(#L, l-> fastExp(L#l, p1-1))),frobeniusPower(maxIdeal,1));
	if(cond==false) then answer=true else answer=false;
	)
    else(
	cond = isSubset((frobeniusPower(I1,1)):I1,frobeniusPower(maxIdeal,1));
	if(cond==false) then answer=true else answer=false;
	);
    answer
)

isFRegularPoly = method();

--Determines if a pair (R, f^t) is F-regular at a prime
--ideal Q in R, R is a polynomial ring  
isFRegularPoly (RingElement, QQ, Ideal) := (f1, t1, Q1) -> (
     not isSubset(tauPoly(f1,t1), Q1)
)

--Determines if a pair (R, f^t) is F-regular, R a polynomial ring
isFRegularPoly (RingElement, QQ) := (f1, t1) -> (
     isSubset(ideal(1_(ring f1)), tauPoly(f1,t1))
)

--Checks whether (R, f1^(a1/(p^e1-1)) is sharply F-pure at the prime ideal m1
isSharplyFPurePoly = (f1, a1, e1,m1) -> (
     if (isPrime m1 == false) then error "isSharplyFPurePoly: expected a prime ideal.";
     not (isSubset(ideal(f1^a1), frobeniusPower(m1,e1)))
)

--Checks whether a Q-Gorenstein pair is strongly F-regular 
isFRegularQGor = method();

--Computes whether (R, f1^t1) is F-regular, assuming the index of R divides p^e1-1
isFRegularQGor (ZZ, RingElement, QQ) := (e1,f1, t1) ->(
     R := ring f1;
     isSubset(ideal(1_R), tauQGor((ring f1),e1,f1,t1))
)

--Computes whether (R, f1^t1) is F-regular at Q1, assuming the index of R divides p^e1-1
isFRegularQGor (ZZ, RingElement, QQ, Ideal) := (e1,f1, t1, Q1) ->(
     not isSubset(tauQGor((ring f1),e1,f1,t1), Q1)
)

--Assuming no pair
isFRegularQGor (Ring,ZZ) := (R,e1) ->(
     isSubset(ideal(1_R), tauQGor(R,e1,1_R,1/1 ) )
)

--Assuming no pair checking at Q1
isFRegularQGor (Ring,ZZ,Ideal) := (R,e1,Q1) ->(
     not isSubset(tauQGor(R,e1,1_R,1/1 ),Q1 )
)


----------------------------------------------------------------
--************************************************************--
--Auxiliary functions for F-signature and Fpt computations.   --
--************************************************************--
----------------------------------------------------------------

--Finds the x-intercept of a line passing through two points
xInt = (x1, y1, x2, y2) ->  x1-(y1/((y1-y2)/(x1-x2)))
 
 
--- Computes the F-signature for a specific value a1/p^e1
--- Input:
---	e1 - some positive integer
---	a1 - some positive integer between 0 and p^e
---	f1 - some polynomial in two or three variables in a ring R of PRIME characteristic
--- Output:
---	returns value of the F-signature of the pair (R, f1^{a1/p^e1})
--- Code is based on work of Eric Canton
fSig = (f1, a1, e1) -> (
     R1:=ring f1;
     pp:= char ring f1;     
     1-(1/pp^(dim(R1)*e1))*
          degree( (ideal(apply(first entries vars R1, i->i^(pp^e1)))+ideal(fastExp(f1,a1) ))) 
)  

--Calculates the x-int of the secant line between two guesses for the fpt
--Input:
--     t - some positive rational number
--     b - the f-signature of (R,f^{t/p^e})
--     e - some positive integer
--     t1- another rational number > t
--     f - some polynomial in two or three variables in a ring of PRIME characteristic
--
-- Output:
--	fSig applied to (f,t1,e)
--	x-intercept of the line passing through (t,b) and (t1,fSig(f,t1,e))
threshInt = (f,e,t,b,t1)-> (
     b1:=fSig(f,t1,e);
{b1,xInt(t,b,t1/(char ring f)^e,b1)}
)
 
--Guesses the FPT of ff.  It returns a list of all numbers in 
--the range suggested by nu(ff,e1) with maxDenom as the maximum denominator
guessFPT ={OutputRange=>false}>>o -> (ff, e1, maxDenom) ->(
     nn := nu(ff, e1);
     pp := char ring ff;
     if (o.OutputRange == false) then 
          findNumberBetween({nn/(pp^e1-1), (nn+1)/(pp^e1)}, maxDenom)
     else
          {{ nn/(pp^e1-1), (nn+1)/(pp^e1)}, findNumberBetween({nn/(pp^e1-1), (nn+1)/(pp^e1)}, maxDenom)}
)

--F-pure threshold estimation, at the origin
--e is the max depth to search in
--FinalCheck is whether the last isFRegularPoly is run (it is possibly very slow) 
--If MultiThread is set to true, it will compute the two F-signatures simultaneously
estFPT={FinalCheck=> true, Verbose=> false, MultiThread=>false, DiagonalCheck=>true, BinomialCheck=>true, NuCheck=>true} >> o -> (ff,ee)->(
     print "starting estFPT";
     
     maxIdeal := ideal( first entries vars( ring ff) );   --the maximal ideal we are computing the fpt at  

     foundAnswer := false; --this will be set to true as soon as we found the answer.  Setting it to true will stop further tests from being run
     answer := null; --this stores the answer until it can be returned.
     
     --first check if it is diagonal:
     if ( (o.DiagonalCheck==true) and (foundAnswer == false) ) then (
	   if (isDiagonal(ff)==true) then ( 
		if (o.Verbose==true) then print "Polynomial is diagonal."; 
		answer = diagonalFPT(ff); 
		foundAnswer = true
	   )
     );

     --now check if it is binomial:
     if ( (o.BinomialCheck==true) and (foundAnswer == false) ) then (
	  if (isBinomial(ff)==true) then ( 
	       if  (o.Verbose==true) then print "Polynomial is binomial.";
	       answer = binomialFPT(ff);
	       foundAnswer = true
	  )
     );
     
     --compute nu's
     if (foundAnswer == false) then (
     	  pp:=char ring ff;
     	  nn:=nu(ff,ee);
	  if  (o.Verbose==true) then print "nu's have been computed";

     	  --if our nu's aren't fine enough, we just spit back some information
       	  if nn==0 then (
	       answer = {0,1/pp};
	       foundAnswer = true
	   )
      );
 
      --check to see if nu/(p^e-1) is the fpt
      if ((o.NuCheck==true) and (foundAnswer == false)) then (
	   if (isFRegularPoly(ff,(nn/(pp^ee-1)),maxIdeal)==false) then ( 
		if  (o.Verbose==true) then print "Found answer via nu/(p^e-1)."; 
		answer = nn/(pp^ee-1);
		foundAnswer = true
	   ) 
      	   else (
	   	if  (o.Verbose==true) then print "nu/(p^e - 1) is not the fpt.";
	   )
      );
	 
	--check to see if (nu+1)/p^e is the FPT
	if ((o.NuCheck==true) and (foundAnswer == false)) then(
		if (isFPTPoly(ff, (nn+1)/pp^ee,Origin=>true) == true) then (
			answer = (nn+1)/pp^ee;
			foundAnswer = true
		)
	);

     --do the F-signature computation
     if (foundAnswer == false) then (
	   ak := 0;
	   if (o.MultiThread==false ) then (ak=threshInt(ff,ee,(nn-1)/pp^ee,fSig(ff,nn-1,ee),nn) ) else(
		if (o.Verbose==true) then print "Beginning multithreaded F-signature";
		allowableThreads = 4;
		numVars := rank source vars (ring ff);
		YY := local YY;
		myMon := monoid[  toList(YY_1..YY_numVars), MonomialOrder=>RevLex,Global=>false];
		--print myMon;
     		R1:=(coefficientRing ring ff) myMon;
		rMap := map(R1, ring ff, vars myMon);
		gg := rMap(ff);
		
		
		H := (fff,aaa,eee) -> () -> fSig(fff,aaa,eee);
		newSig1 := H(gg,nn-1,ee);
		t1 := schedule newSig1;
	     	s2 := fSig(ff,nn,ee);	
		if (o.Verbose==true) then print "One signature down";
		while ((not isReady t1)) do sleep 1;
		s1 := taskResult t1;
     	      --  print s1; print s2;
		ak = {s2,xInt( (nn-1)/pp^ee, s1, (nn)/pp^ee,s2)};
		--print nn;		
	   );
	   if  (o.Verbose==true) then print "Computed F-signatures.";
	   --now check to see if we cross at (nu+1)/p^e, if that's the case, then that's the fpt.
	   if ( (nn+1)/pp^ee == (ak#1) ) then (
		if  (o.Verbose==true) then print "F-signature line crosses at (nu+1)/p^e."; 
		answer = ak#1;
		foundAnswer = true
	   )
      );	  
      	
      --if we run the final check, do the following
      if ( (foundAnswer == false) and (o.FinalCheck == true)) then ( 
	  if  (o.Verbose==true) then print "Starting FinalCheck."; 
          	if ((isFRegularPoly(ff,(ak#1),maxIdeal)) ==false ) then (	
	      		if  (o.Verbose==true) then print "FinalCheck successful"; 
	      		answer = (ak#1);
	      		foundAnswer = true 
      	  	)
	  		else ( 
	      		if  (o.Verbose==true) then print "FinalCheck didn't find the fpt."; 
	      		answer = {(ak#1),(nn+1)/pp^ee};
	      		foundAnswer = true
	  		)
       );
       
       --if we don't run the final check, do the following
       if ((foundAnswer == false) and (o.FinalCheck == false) ) then (
	  if  (o.Verbose==true) then print "FinalCheck not run.";
	  answer = {(ak#1),(nn+1)/pp^ee};
      	  foundAnswer = true
       );
     
     --return the answer
     answer
)

--isFPTPoly, determines if a given rational number is the FPT of a pair in a polynomial ring. 
--if Origin is specified, it only checks at the origin. 

isFPTPoly ={Verbose=> false,Origin=>false}>> o -> (f1, t1) -> (
	pp := char ring f1;
	if (o.Origin == true) then org := ideal(vars (ring f1));
	funList := divideFraction(t1, pp);
	--this writes t1 = a/(p^b(p^c-1))
	aa := funList#0;
	bb := funList#1;
	cc := funList#2;
	mySigma := ideal(f1);
	myTau := ideal(sub(1, ring f1));
	myA := aa;
	myA2 := 0;
	
	if (cc != 0) then (
		myA = floor(aa / (pp^cc - 1));
		myTau = tauPoly( f1, (aa%(pp^cc-1))/(pp^cc-1) )
	);
	
	if (o.Verbose==true) then print "higher tau Computed";

	--first we check whether this is even a jumping number.
	if (cc == 0) then (
		myA2 = aa-1;
		mySigma = sigmaAOverPEMinus1Poly(f1, (pp-1), 1)
	)
	else (
		myA2 = floor((aa-1)/(pp^cc-1));
		mySigma = (sigmaAOverPEMinus1Poly(f1, ((aa-1)%(pp^cc-1))+1, cc))
	);
	if (o.Verbose==true) then print "higher sigma Computed";

	returnValue := false;
	
	if ( isSubset(ideal(sub(1, ring f1)), ethRootSafe(f1, mySigma, myA2, bb) )) then (
		if (o.Verbose==true) then print "we know t1 <= FPT";
		if (not isSubset(ideal(sub(1, ring f1)), ethRootSafe(f1, myTau, myA, bb) ))  then returnValue = true 
	);
		
	returnValue
)


--********************************************
--Some functions for the purpose of checking whether a map of rings is a splitting.  It also computes images of (field) trace.
--********************************************

needsPackage "PushForward"; 

--checks whether f1 : R1 -> S1 splits as a map of R1-modules
isMapSplit = (f1) -> (
	J1 := imageOfRelativeCanonical(f1);
	val := false;
	if (1 % J1 == 0) then val = true;
	
	val
)

--computes the image of Hom_R1(S1, R1) -> R1.
imageOfRelativeCanonical = (f1) -> (
	outList := pushFwd(f1);
--	myGenerators := first entries (outList#1);	
	target1 := (outList#2)(sub(1, target f1));
	
	h1 := map(outList#0, (source f1)^1, target1);
	
	d1 := Hom(h1, (source f1)^1);
	
	trim ideal( first entries gens image d1)
)

--computes the image of trace : S \to R if S is a finite R-module.
imageOfTrace = (f1) -> (
	print "Warning, this only works right now if S is a free module.  We should try to fix it...";
	outList := pushFwd(f1);
	myGenerators := first entries (outList#1);	
	i := 0;
	traceList := {};
	newMap := null;
	newMatrix := null;
	S1 := target f1;
	
	while (i < #myGenerators) do (
		newMap = map(S1^1, S1^1, {{myGenerators#i}});
		newMatrix = pushFwd(newMap, f1);
		traceList = append(traceList, trace newMatrix);
		i = i+1;
	);
	
	trim ideal traceList
)

--computes the relative e-iterated Frobenius over the base ring (the absolute Frobenius in the case 
frobenius = (R1, e1) -> (
	p1 := char R1;
	genList := first entries gens R1;
	fPowerList := apply(genList, z->z^p1 );
	map(R1, R1, fPowerList);
)

--isFJumpingNumberPoly determines if a given rational number is an F-jumping number
--***************************************************************************
--This needs to be speeded up, like the above function
--***************************************************************************
isFJumpingNumberPoly ={Verbose=> false}>> o -> (f1, t1) -> (
	pp := char ring f1;
	funList := divideFraction(t1, pp);
	--this writes t1 = a/(p^b(p^c-1))
	aa := funList#0;
	bb := funList#1;
	cc := funList#2;
	mySigma := ideal(f1);
	myTau := ethRoot(tauPoly(f1, t1*pp^bb), bb);
	if (o.Verbose==true) then print "higher tau Computed";

	--first we check whether this is even a jumping number.
	if (cc == 0) then
		mySigma = ethRoot((ideal(f1^(aa-1)))*((sigmaAOverPEMinus1Poly(f1, (pp-1), 1))), bb)
	else 
		mySigma = ethRoot((sigmaAOverPEMinus1Poly(f1, aa, cc)),bb);
	if (o.Verbose==true) then print "sigma Computed";

	not (isSubset(mySigma, myTau))
)



--MKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMK
-- F-loci (HSL, nonFInjective, simple)

-- Produce a sequence of maps Ext^i(R/I,R) ->  Ext^i(R/I^{[p]},R) induced
-- by the surjections R/I^{[p]} -> R/I
-- for i=1..pdim(coker I)
-- The output consists of a sequence of pairs (A,B) where the induced maps are
-- B: coker A -> coker A^{[p]}
findGeneratingMorphisms = (I) ->
(
	local i;
	Ip:=frobeniusPower(I,1);
	M:=coker I;
	Mp:=coker Ip;
	resM:=res M;
	resMp:=res Mp;
	f:=inducedMap(M,Mp);
	resf:=res f;
	resLength:=length(resM);
	answer:=();
	apply(1..resLength, i->
	{
		G:=resf#i; G=transpose(G);
		F1:=(resM.dd)_(i+1); F1=transpose(F1);
		F0:=(resM.dd)_(i); F0=transpose(F0);
		K:=ker F1;
		C:=subquotient(gens K,F0);
		C1:=prune(C);
		h:=C1.cache.pruningMap;
--
		generatingMorphism0:=G*gens(K)*matrix(entries h);
		F1p:=(resMp.dd)_(i+1); F1p=transpose(F1p);
		F0p:=(resMp.dd)_(i); F0p=transpose(F0p);
		Kp:=ker F1p; 
		Cp:=subquotient(gens Kp,F0p);
		C1p:=prune(Cp);
		hp:=C1p.cache.pruningMap;
--
		A0:=gens(Kp)*matrix(entries hp); 
		A:=A0| F0p;
		gbA:=gb(A, ChangeMatrix => true) ;
		B:=generatingMorphism0// A;
--- Now generatingMorphism0=A*B
		k:=rank source A0;
		generatingMorphism:=submatrix(B,toList(0..(k-1)),);
		answer=append(answer, (C1,generatingMorphism));
---		print(generatingMorphism);
	});
answer
)


----------------------------------------------------------------------------------------
--- Given an Artinian module with Frobenius action F whose Delta functor (=the Matlis dual
--- which keeps track of the given Frobenius action) U: coker A -> coker A^{[p]}
--- produce a sequence (I_0, I_1, ..., I_h) so that
--- the locus of primes on which HSL(F)>h is V(I_h).
--- I_h is always the unit ideal the "global HSL number" is h-1.
--- Note that the "non-F-injective" locus is V(I_0), 
----------------------------------------------------------------------------------------
findHSLloci = (A,U0) ->
(
U:=U0;
local M1;
local M2;
M2=id_(target A); M2=matrix entries M2;
answer:=();
f:=true;
e:=1;
while (f) do
{
	M1=M2;
	M2=ethRoot(U,e) | A;
	W:=subquotient(M1,M2);
	locus:=annihilator W;
	answer=append(answer,locus);
	if (locus==1) then f=false;  --- is this Kosher?
	e=e+1;
	U=U*frobeniusPower(U,1);
};
answer
)


nonFInjectiveLocus = (I) ->
(
R:=ring(I);
answer:=ideal(1_R);
local t;
generatingMorphisms:=findGeneratingMorphisms (I); 
apply (generatingMorphisms, t->
{
	A:=relations ((t)#0);
	if (A!=0) then
	{
		U:=(t)#1;
		HSLLoci:=findHSLloci(A,U);
		answer=intersect(answer,HSLLoci#0);
	};	
});
answer
)


--MKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMKMK

--****************************************************--
--*****************Documentation**********************--
--****************************************************--

beginDocumentation()
doc ///
   Key
      PosChar 
   Headline
      A package for calculations in positive characteristic 
   Description
      Text    
         This will do a lot of cool stuff someday. 
///

doc ///
     Key
     	aPower
     Headline
        Finds the largest power of p dividing x.
     Usage
     	 aPower(x,p)
     Inputs 
		x:ZZ
		p:ZZ
     Outputs
         :ZZ
     Description
	Text
	    Returns the largest exponent e such that p^e divides x.
///

doc ///
     Key
     	ascendIdeal
     Headline
        Finds the smallest phi-stable ideal containing a given ideal in a polynomial ring.
     Usage
     	 ascendIdeal(J, h, e)
     Inputs
     	 J:Ideal 
	h:RingElement
	e:ZZ
     Outputs
         :Ideal
     Description
	Text
	     Let phi be the p^(-e) linear map obtained by multiplying e-th Frobenius trace by h.  Then this function finds the smallest phi-stable ideal containing J.  The idea is to consider the ascending chain J, J+phi(J), J+phi(J)+phi^2(J), etc.  We return the stable value.  For instance, this can be used to compute the test ideal.  This method appared first in the work of Mordechai Katzman on star closure.
///

doc ///
     Key
     	ascendIdealSafe
     Headline
        Finds the smallest phi-stable ideal containing a given ideal in a polynomial ring.
     Usage
     	 ascendIdealSafe(J, h, a, e)
     Inputs
     	 J:Ideal 
	h:RingElement
	a:ZZ
	e:ZZ
     Outputs
         :Ideal
     Description
	Text
	     Let phi be the p^(-e) linear map obtained by multiplying e-th Frobenius trace by h^a.  Then this function finds the smallest phi-stable ideal containing J.  The idea is to consider the ascending chain J, J+phi(J), J+phi(J)+phi^2(J), etc.  We return the stable value.  For instance, this can be used to compute the test ideal.  This method appared first in the work of Mordechai Katzman on star closure.  It differs from ascendIdeal in that it minimizes the exponents that h is raised to, this can make it faster or slower depending on the circumstances.
///

doc ///
     Key
     	basePExp 
     Headline
        Base p Expansion of an integer N
     Usage
     	  basePExp(N,p) 
     Inputs
         N:ZZ
	 p:ZZ
     Outputs
        :List
     Description
     	Text
	     Given an integer N and a prime p, outputs the digits of the base p expansion of N in base p.
///

doc ///
     Key
     	basePExpMaxE
     Headline
        Computes non-terminating base-p expansion of N from digits zero to e-1.
     Usage
     	 basePExpMaxE(N,p,e)
     Inputs
     	 N:ZZ
	 	 p:ZZ
	 	 e:ZZ
     Outputs
         :List
     Description
	Text
	     This computes the base p expansion of N, from digits 0 to e-1.  The digits are given in a list, and come with leading zeros.  If fewer than e digits are required, the list is padded with zeros.  If more digits are required, the final digit lists them.  Little endian is first.  For example, if p=5 and N = 16, the basePExpMaxE(16,5,4) will return {1,3,0,0} (1 one, 3 fives, 0 twentyfives, 0 onehundred twentyfives).
///

doc ///
     Key
     	binomialFPT
     Headline
        Computes the F-pure threshold of a binomial polynomial.
     Usage
     	 binomialFPT(f)
     Inputs 
		f:RingElement
     Outputs
         :QQ
     Description
	Text
	    Returns the F-pure threshold of a binomial in a polynomial ring.  This is based on the work of Daniel Hernandez.
///

doc ///
     Key
     	carryTest
     Headline
        Finds the number of digits we must check to see whether x and y add without carrying.
     Usage
     	 carryTest(w,p)
     Inputs 
		w:List
	p:ZZ
     Outputs
         :ZZ
     Description
	Text
	     Set w = {x,y} a list of rational numbers in [0,1].  This function finds the number of digit places we must check to see if x and y add without carrying.
///

doc ///
     Key
     	denom
     	(denom,ZZ)
     	(denom,QQ)
     Headline
        Returns the denominator of a rational number.
     Usage
     	 denom(x)
     	 denom(y)
     Inputs 
		x:QQ
		y:ZZ
     Outputs
         :ZZ
     Description
	Text
	    Returns the denominator of a rational number or integer (in the latter case it returns 1).
///

doc ///
     Key
     	diagonalFPT
     Headline
        Computes the F-pure threshold of a diagonal polynomial.
     Usage
     	 diagonalFPT(f)
     Inputs 
		f:RingElement
     Outputs
         :QQ
     Description
	Text
	    Returns the F-pure threshold of a diagonal hypersurface in a polynomial ring.  This is based on the work of Daniel Hernandez.
///

doc ///
     Key
     	digit
	(digit,ZZ,QQ,ZZ)
	(digit,ZZ,List,ZZ)
     Headline
        Gives the e-th digit of the base p expansion 
     Usage
     	 d=digit(e,x,p), D=digit(e,X,p)
     Inputs
	e:ZZ
	x:QQ
	p:ZZ
	X:List
	   consisting of rational numbers
     Outputs
        d:ZZ
	    which is the e-th digit of the non-terminating base p expansion of x
	D:List
	    which contains the e-th digits of the entries of the list X
     Description
	Text
	     Gives the e-th digit, to the right of the decimal point, of the non-terminating base p expansion of x in [0,1]; threads over lists of rational numbers. 
///

doc ///
     Key
     	divideFraction
     Headline
        Converts a rational number into something of the form (a/(p^b p^(c-1)).
     Usage
     	 divideFraction(t, p)
     Inputs 
		t:QQ
	p:ZZ
     Outputs
         :List
     Description
	Text
	     Given a rational number t and prime p, this function finds a list of integers {a,b,c} such that t= (a/(p^b p^(c-1)).
///

doc ///
     Key
     	 estFPT
     Headline
         Atempts to compute the F-pure threshold, where e is the max depth to search in.  
     Usage
     	  estFPT(f,e,finalCheck=>V,Verbose=>W)
     Inputs
     	 f:RingElement
         e:ZZ
	 V:Boolean
	 W:Boolean
     Outputs
        L:List
	Q:QQ
     Description
     	  Text 
	      This tries to find an exact value for the fpt.  If it can, it returns that value.  Otherwise it should return a range of possible values (eventually).  It first checks to see if the ring is binonmial or diagonal.  In either case it uses methods of D. Hernandez.  Next it tries to estimate the range of the FPT using nu's.  Finally, it tries to use this to deduce the actual FPT via taking advantage of convexity of the F-signature function and a secant line argument.  finalCheck is a Boolean with default value True that determines whether the last isFRegularPoly is run (it is possibly very slow).  If FinalCheck is false, then a last time consuming check won't be tried.  If it is true, it will be.  Verbose set to true displays verbose output.
///

doc ///
     Key
     	 ethRoot
     Headline
        Computes $I^{[1/p^e]}$ in a polynomial ring over a perfect field
     Usage
     	  ethRoot(I,e) 
     Inputs
     	 I:Ideal
         e:ZZ
     Outputs
        :Ideal
     Description
	Text
	     In a polynomial ring k[x1, ..., xn], I^{[1/p^e]} is the smallest ideal J such that J^{[p^e]} = FrobeniusPower(J,e) \supseteq I.  This function computes it.
///

doc ///
     Key
     	ethRootSafe
     Headline
        Computes (f^a*I)^{[1/p^e]} in such a way that we don not blow exponent buffers.
     Usage
     	 ethRootSafe(f, I, a, e)
     Inputs
     	 f:RingElement
	 I:Ideal
	 a:ZZ
	 e:ZZ
     Outputs
         :Ideal
     Description
	Text
	     Computes the 1/p^e-th root of (f^a*I).  It does it while trying to minimize the power that f gets raised to (in case a is a large number).  This can either be faster or slower than ethRoot.
///

doc ///
     Key
     	fastExp 
     Headline
        Computes powers of elements in rings of characteristic p>0 quickly.
     Usage
     	  fastExp(f,N) 
     Inputs
     	 f:RingElement
         N:ZZ
     Outputs
        :RingElement
     Description
	Text
	     In prime characteristic p > 0, raising a sum (a+b) to a power is more quickly done by simply computing a^p and b^p and adding them.  The basic strategy is to break up the exponent into its base p expansion, and then use the exponent rules.  For example, (x+y)^(3*p^2 + 5*p+2) = ((x+y)^3)^(p^2)*((x+y)^5)^p*(x+y)^2.
///

doc ///
     Key
     	findQGorGen
     Headline
        If R = S/I where S is a polynomial ring, returns the ring element with I^{[p^e]} : I = (f) + I^{[p^e]}.
     Usage
     	 findQGorGen(R, e)
     Inputs
     	 R:Ring
     Outputs
         :RingElement
     Description
	Text
	     If R is Q-Gorenstein with index not divisible by p, then I^{[p^e]} : I = (f) + I^{[p^e]}.  For some e.  This function tries to find the f.  If the argument e is left out then e is assumed to be 1.
///

doc ///
     Key
     	firstCarry
     Headline
        Finds the first spot where (the eth digit of x) + (the eth digit of y) >= p.
     Usage
     	 firstCarry(w,p)
     Inputs 
		w:List
	p:ZZ
     Outputs
         :ZZ
     Description
	Text
	     Set w = {x,y} a list of rational numbers in [0,1].  Finds the first place where (the eth digit of x) + (the eth digit of y) >= p, in other words where the numbers add with carrying.
///

doc ///
     Key
     	 FPTApproxList
	 (FPTApproxList,Ideal,ZZ)
	 (FPTApproxList,RingElement,ZZ)
     Headline
        Gives a list of nu_I(p^d)/p^d for d=1,...,e.
     Usage
     	  FPTApproxList(I,e)
	  FPTApproxList(f,e) 
     Inputs
     	 I:Ideal
	 f:RingElement
         e:ZZ
     Outputs
         :List
     Description
	Text 
 	     This returns a list of nu_I(p^d)/p^d for d = 1, ..., e.  The {nu_I(p^d)/p^d} converge to the F-pure threshold.	     
///

doc ///
     Key
     	FPT2VarHomog
	(FPT2VarHomog,RingElement)
	(FPT2VarHomog,List,List)
     Headline
        F-pure threshold of a form in two variables
     Usage
     	  fpt=FPT2VarHomog(G), fpt=FPT2VarHomog(factors,multiplicities)
     Inputs 
	factors:List
	    which contains the linear factors of a form G in two variables 
	multiplicities:List
	    which contains the multiplicities of those linear factors in G
	G:RingElement
	    a form in two variables
     Outputs
        fpt:QQ
     Description
	Text
	    FPT2VarHomog computes the F-pure threshold of a homogeneous polynomial G
	    	in two variables. 
	    The polynomial G can be entered directly, or if the user knows a factorization
	    	G=L1^(a1)...Ln^(an) into linear forms, that can be used for improved 
		performance: FPT2VarHomog({L1,...,Ln},{a1,...,an}).
///

doc ///
     Key
     	 frobeniusPower
     Headline
        The following raises an ideal to the $p^e$th power.
     Usage
     	  frobeniusPower(I,e) 
     Inputs
     	 I:Ideal
         e:ZZ
     Outputs
        :Ideal
     Description
	Text
	     If I = ideal(x1, ..., xn), then frobeniusPower(I,e) outputs ideal(x1^(p^e), ..., xn^(p^e)) where p is the characteristic of the ring.
///

doc ///

     Key
     	 fSig
     Headline
        Computes the F-signature for a specific value $a/p^e$.
     Usage
     	  fSig(f,a,e)
     Inputs
     	 f:RingElement
	 a:ZZ
         e:ZZ
     Outputs
        :QQ
     Description
	Text
	     This computes the F-signature $s(R, f^{a/p^e})$ if R is a polynomial ring over a perfect field.
///

doc ///
     Key
     	 FTApproxList
	 (FTApproxList,Ideal,Ideal, ZZ)
	 (FTApproxList,RingElement,Ideal,ZZ)
     Headline
        Gives a list of nu_I^J(p^d)/p^d for d=1,...,e.
     Usage
     	  FPTApproxList(I,J,e)
	  FPTApproxList(f,J,e) 
     Inputs
     	 I:Ideal
	 J:Ideal
	 f:RingElement
         e:ZZ
     Outputs
         :List
     Description
	Text 
 	     This returns a list of nu_I^J(p^d)/p^d for d = 1, ..., e.  The {nu_I^J(p^d)/p^d} converge to the F-threshold.	     
///

doc ///
     Key
     	genFrobeniusPower 
     Headline
        Computes the generalized Frobenius power of an ideal
     Usage
     	  frobeniusPower(I1,e1)  
     Inputs
         	I1:Ideal
	 	e1:ZZ
     Outputs
        :Ideal
     Description
     	Text
	     Computes I^[N] for an ideal I and an integer N, where I^[N] is defined as follows. If N's base P-expansion is N=n_0+n_1P+...+n_eP^e then I^[N]=I^(n_0)*(I^(n_1))^[P]*...*(I^(n_e))^[P^e]. When P is prime I^[P^e] is the usual Frobenius power.
 ///
 
doc ///
     Key
     	guessFPT 
     Headline
        Tries to guess the FPT in a really naive way (this should be improved).
     Usage
     	  guessFPT(f,e,d) 
     Inputs
     	 f:RingElement
         e:ZZ
	 d:ZZ
     Outputs
        :List
     Description
	Text
	     This tries to guess the FPT.  In particular, it computes the number nu such that nu/(p^e - 1) <= FPT < (nu+1)/p^e.  It then outputs a list of all rational numbers with denominators less than or equal to d, which lie in that range.  WARNING:  There are several improvements which should be made to this function to rule out many of the possibilies.
///

doc ///
     Key
     	isBinomial 
     Headline
        Checks whether a polynomial is binomial.
     Usage
     	 isBinomial(f)
     Inputs 
		f:RingElement
     Outputs
         :Boolean
     Description
	Text
	    Returns true if f is a binomial, otherwise returns false.
///

doc ///
     Key
     	isDiagonal 
     Headline
        Checks whether a polynomial is diagonal.
     Usage
     	 isDiagonal(f)
     Inputs 
		f:RingElement
     Outputs
         :Boolean
     Description
	Text
	    Returns true if f is a diagonal, otherwise returns false.  Recall f is called diagonal if it is of the form x_1^(a_1)+...+x_n^(a_n) up to renumbering of the variables.
///

doc ///
     Key
     	isFJumpingNumberPoly 
     Headline
        Checks whether a given number is the FPT
     Usage
     	  isFJumpingNumberPoly(f,t,Verbose=>V)  
     Inputs
         	f:RingElement
	 	t:ZZ
		W:Boolean
     Outputs
        :Boolean
     Description
     	Text
	     Returns true if t is an F-jumping number, otherwise it returns false.
///

doc ///
     Key
     	isFPTPoly 
     Headline
        Checks whether a given number is the FPT
     Usage
     	  isFPTPoly(f,t,Verbose=>V,Origin=>W)  
     Inputs
         	f:RingElement
	 	t:ZZ
		W:Boolean
		W:Origin
     Outputs
        :Boolean
     Description
     	Text
	     Returns true if t is the FPT, otherwise it returns false.  If Origin is true, it only checks it at ideal(vars ring f).
///

doc ///
     Key
     	 isFPure 
     Headline
         Tests for a given ideal I, if R/I is F-pure. 
     Usage
     	 isFPure(I)
     Inputs
     	 I:Ideal
     Outputs
         :Boolean
     Description
	Text 
	    In the case where I is a complete intersection,this function applies Fedder's Criterion.
	    Otherwise, checks if I^[p]:I is contained in m^[p]. 
 	   
///

doc ///
     Key
     	 isFRegularPoly
     Headline
        Determines if a pair $(R, f^t)$ is F-regular when R is a polynomial ring. 
     Usage
     	  isFRegularPoly
     Inputs
     	 f:RingElement
         t:QQ
     Outputs
        :Boolean
     Description
	Text
	     This computes the test ideal.  The ring is F-regular if the test ideal is the whole ring, in which case this function returns true.  Otherwise, this function returns false.

///

doc ///
     Key
     	isFRegularQGor
	 (isFRegularQGor, ZZ, RingElement, QQ)
	 (isFRegularQGor, ZZ, RingElement, QQ, Ideal)
	 (isFRegularQGor, Ring, ZZ)
	 (isFRegularQGor, Ring, ZZ, Ideal)
     Headline
        Checks whether a ring or a pair is Q-Gorenstein.
     Usage
     	 isFRegularQGor(e,f,t)
     	 isFRegularQGor(e,f,t,Q)
     	 isFRegularQGor(R,e)
     	 isFRegularQGor(R,e,Q)
     Inputs 
		R:Ring
	 f:RingElement
	 e:ZZ
	 t:QQ
	 Q:Ideal
     Outputs
         :Boolean
     Description
	Text
	     Checks whether R, or the pair (R, f^t),  is strongly F-regular at Q (respectively the origin).  It assumes the Q-Gorenstein index divides (p^e - 1).
///

doc ///
     Key
     	 isSharplyFPurePoly
     Headline
        Checks whether (R, f^(a/(p^e - 1))) is F-pure at the prime ideal m.
     Usage
     	 isSharplyFPurePoly(f,a,e,m)
     Inputs
     	 f:RingElement
	 a:ZZ
         e:ZZ
	 m:Ideal
     Outputs
         :Boolean
     Description
	Text
	     This checks whether (R, f^(a/(p^e-1))) is F-pure at the prime ideal m at least in the case that R is a polynomial ring.
///

doc ///
     Key
     	 nu
	 (nu,Ideal,ZZ)
	 (nu,RingElement,ZZ)
     Headline
        Gives $\nu_I(p^e)$.
     Usage
     	  nu(I,e)
	  nu(f,e) 
     Inputs
     	 I:Ideal
	 f:RingElement
         e:ZZ
     Outputs
        :ZZ
     Description
	Text
	    Given an ideal I in a polynomial ring k[x1, ..., xn], this function outputs the smallest integer nu such that I^nu is not in ideal(x1^(p^e), ..., xn^(p^e) ).  If a RingElement is passed, it computes nu of the principal ideal generated by this element. This is used frequently to compute the F-pure threshold.
///

doc ///
     Key
     	 nuList
	 (nuList,Ideal,ZZ)
	 (nuList,RingElement,ZZ)
     Headline
        Lists $\nu_I(p^d)$ for d = 1,...,e.
     Usage
     	  nuList(I,e)
	  nuList(f,e) 
     Inputs
     	 I:Ideal
	 f:RingElement
         e:ZZ
     Outputs
        :List
     Description
	Text
	     Given an ideal I in a polynomial ring k[x1,...,xn], this function computes nu(I,d) for d = 1,...,e.
///

doc ///
     Key
     	sigmaAOverPEMinus1Poly
     Headline
        Computes the non-sharply F-pure ideal of (R, f^{a/(p^e-1)}) when R is a polynomial ring.
     Usage
     	 sigmaAOverPEMinus1Poly (f, a, e, HSL=>W)
     Inputs 
		f:RingElement
	a:ZZ
	e:ZZ
	W:Boolean
     Outputs
         :Ideal
     Description
	Text
	     Let phi be the p^(-e) linear map obtained by multiplying e-th Frobenius trace by f^a.  This computes \phi^n(R) for large n.  This stabilizes by Hartshorne-Speiser-Lyubeznik-Gabber.  If HSL is true, then the function returns a list where the first entry is sigma and the second entry is the HSL number.
///

doc ///
     Key
     	sigmaAOverPEMinus1QGor
     Headline
        Computes the non-sharply F-pure ideal of (R, f^{a/(p^e-1)}).
     Usage
     	 sigmaAOverPEMinus1QGor(f, a, e, g,HSL=>W)
     Inputs 
		f:RingElement
		a:ZZ
		e:ZZ
		g:ZZ
		W:Boolean
     Outputs
         :Ideal
     Description
	Text
	     Let phi be the p^(-e) linear map obtained by multiplying e-th Frobenius trace of R by f^a (we assume that the Q-Gorenstein index of R divides p^g-1).  This computes \phi^n(R) for large n.  This stabilizes by Hartshorne-Speiser-Lyubeznik-Gabber.  If HSL is true, then the function returns a list where the first entry is sigma and the second entry is the HSL number of sigma(ring f) relative to f.
///

doc ///
     Key
     	sigmaQGor
     Headline
        Computes the non-sharply F-pure ideal of R, where R is Q-Gorenstein with index dividing (p^g-1).
     Usage
     	 sigmaQGor(R, g,HSL=>W)
     Inputs 
		R:Ring
		g:ZZ
		W:Boolean
     Outputs
         :Ideal
     Description
	Text
	     Let phi be the  g-th Frobenius trace of R (we assume that g is the Q-Gorenstein index of R).  This computes \phi^n(R) for large n.  This stabilizes by Hartshorne-Speiser-Lyubeznik-Gabber.  If HSL is true, then the function returns a list where the first entry is sigma and the second entry is the HSL number.
///

doc ///
     Key
     	tauAOverPEMinus1Poly
     Headline
        Computes the test ideal of f^(a/(p^e-1)) if f is in a polynomial ring.
     Usage
     	 tauAOverPEMinus1Poly(f, a, e)
     Inputs
     	 f:RingElement
	 a:ZZ
	 e:ZZ
     Outputs
         :Ideal
     Description
	Text
	     Computes the test ideal tau(f^(a/(p^e-1)) ).  The basic idea first appeared in a paper of Mordechai Katzman.
///

doc ///
     Key
     	tauGor
     Headline
        Computes tau(R,f^t) for a Gorenstein ring such that the index divides p^e-1.
     Usage
     	 tauGor(R,f,t)
     Inputs
     	 R:Ring
	 f:RingElement
	 t:QQ
     Outputs
         :Ideal
     Description
	Text
	     This computes the test ideal tau(R, f^t) for a Gorenstein ring.  First the test ideal of the ambient space is computed (and computed on a polynomial ring S of which R is a quotient).  Then writing t = a/(p^b-1)p^c we compute tau(R, f^{a/(p^b-1)}), or rather a preimage of it on S, by summing the images of the map induced by f^{a/(p^b-1)}.  We then compute tau(R, f^t) by multiplying by the output of a findQGorGen on S, and taking [1/p^e]th roots on S.
///

doc ///
     Key
     	tauGorAmb
     Headline
        Computes tau(R) for a Gorenstein ring.
     Usage
     	 tauGorAmb(R)
     Inputs
     	 R:Ring
     Outputs
         :Ideal
     Description
	Text
	     This computes the test ideal tau(R) for a quasi-Gorenstein ring R.  It uses the fact that if R is a quotient of a polynomial ring S, then tau(R) can be computed as a sort of test/adjoint ideal on S.  The function findQGorGen is used to find the map to use on S. 
///

doc ///
     Key
     	 tauPoly
     Headline
        Computes the test ideal of $(R, f^t)$.
     Usage
     	  tauPoly(f,t) 
     Inputs
     	 f:RingElement
         t:QQ
     Outputs
        :Ideal
     Description
	Text
	     This computes the test ideal of (R, f^t) when R is a polynomial ring over a perfect field.  It is done as follows.  If t = a/(p^e - 1) then tau(R, f^t) is computed as a sum of (f^{\lceil t \rceil}*f^{\lceil t(p^e-1) \rceil})^{[1/p^e]} until the sum stabilizes.  For the more general case, we use the formula tau(R, f^t)^{[1/p^d]} = tau(R, f^{t/p^d}).
///

doc ///
     Key
     	tauQGorAmb
     Headline
        Computes tau(R) for a Q-Gorenstein ring with index not dividing p^e - 1.
     Usage
     	 tauQGorAmb(R,e)
     Inputs
     	 R:Ring
	 e:ZZ
     Outputs
         :Ideal
     Description
	Text
	     This computes the test ideal tau(R) for a Q-Gorenstein ring R with index dividing p^e - 1.  It uses the fact that if R is a quotient of a polynomial ring S, then tau(R) can be computed as a sort of test/adjoint ideal on S.  The function findQGorGen is used to find the map to use on S.  e is the index of the canonical divisor on R.
///

doc ///
     Key
     	tauQGor
     Headline
        Computes tau(R,f^t) for a Q-Gorenstein ring such that the index divides p^e-1.
     Usage
     	 tauQGor(R,e,f,t)
     Inputs
     	 R:Ring
	 e:ZZ
	 f:RingElement
	 t:QQ
     Outputs
         :Ideal
     Description
	Text
	     This computes the test ideal tau(R, f^t) for a Q-Gorenstein ring such that the index divides p^e -1.  First the test ideal of the ambient space is computed (and computed on a polynomial ring S of which R is a quotient).  Then writing t = a/(p^b-1)p^c we compute tau(R, f^{a/(p^b-1)}), or rather a preimage of it on S, by summing images of the map induced by f^{a/(p^b-1)}.  We then compute tau(R, f^t) by multiplying by the output of a findQGorGen on S, and taking [1/p^e]th roots on S.
///

doc ///
     Key
     	truncation
	(truncation,ZZ,QQ,ZZ)
	(truncation,ZZ,List,ZZ)
     Headline
        Truncations of base p expansions of rational numbers
     Usage
     	 t=truncation(e,x,p), T=truncation(e,X,p)
     Inputs 
	e:ZZ
	x:QQ
	p:ZZ
	X:List
	   which contains rational numbers
     Outputs
        t:QQ
	    which is the e-th truncation of the non-terminating base p expansion of x
	T:List
	    which contains the e-th truncations of the entries of the list X
     Description
	Text
	     Gives the first e digits of the non-terminating base p expansion of a nonnegative rational number x, as a fraction; threads over lists of rational numbers.
///

doc ///
     Key
     	truncationBaseP
     Headline
        Gives the first e digits of the non-terminating base p expansion of x.
     Usage
     	 truncationBaseP(e,x,p)
     Inputs 
		e:ZZ
	x:QQ
	p:ZZ
     Outputs
         :List
     Description
	Text
	     Gives the first e digits of the non-terminating base p expansion of x in [0,1], as a list.
///

end

--**********************************
--Changes in 0.2a
----Fixed some typos in documentation and comments
----Commented out moduleToIdeal, replaced with needsPackage "Divisor" which has a better version of moduleToIdeal
--- Fixed things

---Zhibek was here 