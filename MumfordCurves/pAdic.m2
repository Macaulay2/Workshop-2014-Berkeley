PAdicField = new Type of InexactField
PAdicFieldElement = new Type of HashTable

new PAdicField from List := (PAdicField, inits) -> new PAdicField of PAdicFieldElement from new HashTable from inits

pAdicField = method()

makePAdicField:=(R,p)->(
   A := new PAdicField from {(symbol baseRing) => R,
                            (symbol uniformizingParameter) => p};
 --  net A := a->(expans:=a#"expansion";
--	d:=min keys expans;
--	n:=a#"precision";
--	if d<0 then 
--		  if 
   A + A := (a,b) -> 3;
     A 
)
-3..2
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

4p500===4p300

restart
load "~/git/Workshop-2014-Berkeley/MumfordCurves/pAdic.m2"
Q3=pAdicField(3)

x=toPAdicFieldElement({1,0,0,1},Q3)

kZ3 = pAdicIntegers(3)
x = toPAdicElement(Z3,{(0,1)})
x+x
class 3
ZZ[x]
x
class x
new RingElement from 3
help RingElement
RingElement
