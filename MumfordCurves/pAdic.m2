PAdicField = new Type of InexactField
PAdicFieldElement = new Type of HashTable

new PAdicField from List := (PAdicField, inits) -> new PAdicField of PAdicFieldElement from new HashTable from inits

valuation = method()
relativePrecision = method()

makePAdicField:=(R,p)->(
   A := new PAdicField from {(symbol baseRing) => R,
                            (symbol uniformizingParameter) => p};
   precision A := a->a#"precision";
   valuation A := a->(if #a#"expansion">0 then return min a#"expansion"_0;
	infinity);
   relativePrecision A:= a -> (precision a)-(valuation a);
   net A := a->(expans:=a#"expansion";
	keylist:=expans_0;
       	((concatenate apply(#keylist,i->
		  toString(expans_1_i)|"*"|toString p|"^"
		  |toString keylist_i|"+"))
		|"O("|toString p|"^"|toString(precision a)|")"));
   computeCarryingOver := (aKeys,aValues,prec) -> (
	newKeys := ();
	newValues := ();
	carryingOver := 0_R;
	aPointer := 0;
	while (aPointer<#aKeys and aKeys#aPointer<=prec) do (
	     currentKey := aKeys#aPointer;
	     currentValue := aValues#aPointer+carryingOver;
	     carryingOver = 0_R;
	     while currentValue!=0 do (
	     	  q := currentValue%p;
		  currentValue = (currentValue-q)//p;
		  if q!=0 then (
		       newKeys = (newKeys,currentKey);
		       newValues = (newValues,q);
		       );
		  currentKey = currentKey+1;
		  if (currentKey>=prec or 
		       ((aPointer+1<#aKeys) and 
			    (currentKey>=aKeys#(aPointer+1)))) then (
		       carryingOver = currentValue;
		       break;
		       );
	     	  );
	     aPointer = aPointer+1;
	     );
	new A from {"precision"=>prec,
	     "expansion"=>{toList deepSplice newKeys,
		  toList deepSplice newValues}}
	);
   A + A := (a,b) -> (
	newPrecision := min(a#"precision",b#"precision");
	aKeys := a#"expansion"_0;
	aValues := a#"expansion"_1;
	bKeys := b#"expansion"_0;
	bValues := b#"expansion"_1;
	newKeys := ();
	newValues := ();
	aPointer := 0;
	bPointer := 0;
	while (aPointer<#aKeys) or (bPointer<#bKeys) do (
	     newKey := Nothing;
	     newValue := Nothing;
	     if (bPointer==#bKeys) then (
		  newKey = aKeys#aPointer;
		  newValue = aValues#aPointer;
		  aPointer = aPointer + 1;
		  ) else if (aPointer==#aKeys) then (
		  newKey = bKeys#bPointer;
		  newValue = bValues#bPointer;
		  bPointer = bPointer + 1;
		  ) else if (aKeys#aPointer<bKeys#bPointer) then (
		  newKey = aKeys#aPointer;
		  newValue = aValues#aPointer;
		  aPointer = aPointer + 1;
		  ) else if (bKeys#bPointer<aKeys#aPointer) then (
		  newKey = bKeys#bPointer;
		  newValue = bValues#bPointer;
		  bPointer = bPointer + 1;
		  ) else (
		  newKey = bKeys#bPointer;
		  newValue = aValues#aPointer + bValues#bPointer;
		  aPointer = aPointer + 1;
		  bPointer = bPointer + 1;
		  );
	     if (newKey>=newPrecision) then (
		  break;
		  );
	     newKeys = (newKeys,newKey);
	     newValues = (newValues,newValue);
	     );
	newKeys = toList deepSplice newKeys;
	newValues = toList deepSplice newValues;
	computeCarryingOver(newKeys,newValues,newPrecision)
	);
   A * A := (a,b)->(
	newPrecision := min(precision a+min(precision b,valuation b),
	     precision b+min(precision a,valuation a));
	aKeys := a#"expansion"_0;
	aValues := a#"expansion"_1;
	bKeys := b#"expansion"_0;
	bValues := b#"expansion"_1;
	prod := new MutableHashTable;
	for i in 0..#aKeys-1 do (
	     for j in 0..#bKeys-1 do (
		  newKey := aKeys#i+bKeys#j;
		  if newKey<newPrecision then (
		       newValue := aValues#i*bValues#j;
		       if prod#?newKey then (
			    prod#newKey = prod#newKey + newValue;
			    ) else (
			    prod#newKey = newValue;
			    );
		       );
		  );
	     );
	newKeys := sort keys prod;
	newValues := for i in newKeys list prod#i;
	computeCarryingOver(newKeys,newValues,newPrecision)
	);
     A 
)

pAdicField = method()
pAdicField(Ring,RingElement):=(R,p)->makePAdicField(R,p)
pAdicField(Ring,ZZ):=(R,p)->makePAdicField(R,p)
pAdicField(ZZ):=(p)->makePAdicField(ZZ,p)

-- PAdicField Elements are hashtables with following keys:
-- precision (value ZZ)
-- expansion (hashtable, two entries: exponents, coefficients)




toPAdicFieldElement = method()

toPAdicFieldElement (List,PAdicField) := (L,S) -> (
   n:=#L;
   local expans;
   if all(L,i->i==0) then expans={{},{}}    
   else expans=entries transpose matrix select(apply(n,i->{i,L_i}),j->not j_1==0);
   new S from {"precision"=>n,"expansion"=>expans}
   )

end
----------------------------
--Qingchun's testing area
----------------------------
restart
load "/Users/qingchun/Desktop/M2Berkeley/Workshop-2014-Berkeley/MumfordCurves/pAdic.m2"
Q3 = pAdicField(3)
x = toPAdicFieldElement({1,2,0,1,0},Q3);
y = new Q3 from {"precision"=>3,"expansion"=>{{},{}}};
print(x+x)
print(x*x)
print(x+y)
print(x*y)
print(y*y)
end

----------------------------
--Nathan's testing area
----------------------------
restart
load "pAdic.m2"
load "~/Workshop-2014-Berkeley/MumfordCurves/pAdic.m2"
Q3=pAdicField(3)

x=toPAdicFieldElement({0,1,0,0,1,0},Q3)

precision x
relativePrecision x
valuation x
f=x


ZZ[y]
f=y^2+y^6
jacobian f
help jacobian
diff(y,f)
help diff
