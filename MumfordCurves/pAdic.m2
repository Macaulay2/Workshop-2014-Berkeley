PAdicField = new Type of InexactField
PAdicFieldElement = new Type of HashTable

new PAdicField from List := (PAdicField, inits) -> new PAdicField of PAdicFieldElement from new HashTable from inits


valuation = method()
relativePrecision = method()





makePAdicField:=(R,p)->(
   A := new PAdicField from {(symbol baseRing) => R,
                            (symbol uniformizingParameter) => p};

     computeCarryOver := (a,prec) -> (
	  
	  );

     A + A := (a,b) -> (
	  newPrecision := min(a#"precision",b#"precision");
	  s := new MutableHashTable from a#"expansion";
	  for i in keys a#"expansion" do (
	       if i<newPrecision then (
		    s#i := a#"expansion"#i;
	       );
	  );
     	  for i in keys b#"expansion" do (
	       if i<newPrecision then (
	       	    if s#?i then (
		    	 s#i := s#i + b#"expansion"#i;
		    ) else (
		    	 s#i := b#"expansion"#i;
		    );
	       );
	  computeCarryOver(s,newPrecision);
	  new A from {"precision"=>newPrecision,"expansion"=>s};
	  );

   precision A := a->a#"precision";
   valuation A := a->min keys a#"expansion";
   relativePrecision A:= a -> (precision a)-(valuation a);
   net A := a->(expans:=a#"expansion";
	keylist:=sort keys expans;
	((concatenate apply(keylist,i->
		  toString(expans#i)|"*"|toString p|"^"
		  |toString i|"+"))
	|"O("|toString p|"^"|toString(precision a)|")"));
     A 
)

pAdicField = method()
pAdicField(Ring,RingElement):=(R,p)->makePAdicField(R,p)
pAdicField(Ring,ZZ):=(R,p)->makePAdicField(R,p)
pAdicField(ZZ):=(p)->makePAdicField(ZZ,p)

-- PAdicField Elements are hashtables with following keys:
-- precision (value ZZ)
-- expansion (hashtable, keys are integers, values elements of baseRing)




toPAdicFieldElement = method()

toPAdicFieldElement (List,PAdicField) := (L,S) -> (
   n:=#L;
   expans:=new HashTable from select(apply(n,i->(i,L_i)),j->not j_1==0);
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

