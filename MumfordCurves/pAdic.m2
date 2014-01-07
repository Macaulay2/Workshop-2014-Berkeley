PAdicField = new Type of InexactField
PAdicFieldElement = new Type of HashTable

new PAdicField from List := (PAdicField, inits) -> new PAdicField of PAdicFieldElement from new HashTable from inits


valuation = method()
relativePrecision = method()





makePAdicField:=(R,p)->(
   A := new PAdicField from {(symbol baseRing) => R,
                            (symbol uniformizingParameter) => p};
   precision A := a->a#"precision";
   valuation A := a->min a#"expansion"_0;
   relativePrecision A:= a -> (precision a)-(valuation a);
   net A := a->(expans:=a#"expansion";
	keylist:=expans_0;
	((concatenate apply(#keylist,i->
		  toString(expans_1_i)|"*"|toString p|"^"
		  |toString keylist_i|"+"))
	|"O("|toString p|"^"|toString(precision a)|")"));
   A + A := (a,b) -> 3;
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
   expans:=select(apply(n,i->(i,L_i)),j->not j_1==0);
   new S from {"precision"=>n,"expansion"=>expans}
   )

end

restart
load "~/git/Workshop-2014-Berkeley/MumfordCurves/pAdic.m2"
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
