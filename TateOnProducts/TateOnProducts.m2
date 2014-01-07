newPackage(
	"TateOnProducts",
    	Version => "0.1", 
    	Date => "October 18, 2013",
    	Authors => {{Name => "Frank-Olaf Schreyer", 
		  Email => "schreyer@math.uni-sb.de", 
		  HomePage => "http://www.math.uni-sb.de/ag/schreyer/"}},
--	             {Name => "David Eisenbud", 
--		  Email => "de@msri.org", 
--		  HomePage => "http://www.msri.org/~de/"}
--                   },
    	Headline => "Tate resolutions on products of projective spaces",
    	DebuggingMode => true
    	)

export {
    setupRings,
    numFactors,
    symExt,
    cohomologyTable,
    tallyDegrees,
    upperCorner,
    lowerCorner,
    beilinsonWindow,
    tateExtensionOld,
    truncateInE,
    tateExtension,
    firstQuadrantComplex,
    lastQuadrantComplex,
    cornerComplex,
    isChainComplex,
    regionComplex,
    strand
    }

setupRings=method()
setupRings(Ring,List) := (kk,n) -> (
     t:= #n;
     x:= symbol x;
     xx:=flatten apply(t,i->apply(n_i+1,j->x_(i,j)));
     degs:=flatten apply(t,i->apply(n_i+1,k->apply(t,j->if i==j then 1 else 0)));
     Sloc:=kk[xx,Degrees=>degs];
     e:= symbol e;
     ee:=flatten apply(t,i->apply(n_i+1,j->e_(i,j)));
     Eloc:=kk[ee,Degrees=>degs,SkewCommutative=>true];
--     h:=symbol h;
--     H:=ZZ[h];
     return(Sloc,Eloc))
TEST ///
n={1,2}
kk=ZZ/101
(S,E)=setupRings(kk,n)
scan(#n,i->scan(n_i+1,j->x_(i,j)=S_(sum(i,k->n_k+1)+j)))
scan(#n,i->scan(n_i+1,j->e_(i,j)=E_(sum(i,k->n_k+1)+j)))

m=matrix{{x_(0,0),x_(1,0)},
       {x_(0,1),0},
       {0,x_(1,1)},
       {0,x_(1,2)}}
betti m
mE=symExt(m,E)
betti (T= res coker mE)
tallyDegrees T
cohomologyTable( T, -{2,2},{6,6})
///

numFactors=method(TypicalValue=>List)
    -- given the symmetric or exterior Cox ring E of the product pf projective spaces
    -- compute the number t of factors and return the List {0,...,t-1} 
numFactors(Ring) := E-> (
    t:= #unique degrees E;
    return toList(0..(t-1)))




cohomologyTable=method()

cohomologyTable(ChainComplex,List,List) := (F,da,db) -> (
     E:= ring F; 
     if not #unique degrees E==2 then error "works only for two factors";
     L:=flatten apply(toList(min F..max F),k->apply(degrees F_k,deg->sum deg-k));
     minL:=min L;maxL := max L;
     h:=symbol h;
     H:=ZZ[h];
     C:=matrix apply(toList(da_0..db_0),i->apply(toList(da_1..db_1),j->sum(0..(maxL- minL) ,p-> h^p*(tally degrees F_(i+j-p-minL))_{i,j})));
     C_(reverse apply(rank source C,i->i))
     )
TEST ///
n={1,2}
(S,E)=setupRings(ZZ/101,n)
use S
vars S
m=map(S^4,S^{{ -1,0},{0,-1}}, transpose matrix{{S_0,S_1,0,0},{S_2,0,S_3,S_4}})
mE=symExt(m,E)
time fB=dual res(coker transpose mE,LengthLimit=>10)
cohomologyTable(fB,{ -5,-5},{0,0})
tallyDegrees fB
///

tallyDegrees=method()
tallyDegrees(ChainComplex) := C -> (
    apply(toList(min C..max C),k->tally degrees C_k))



boxDegrees=method()
boxDegrees(Ring) := E -> (
     degs:= unique degrees E;
     t:=#degs;
     n:=apply(t,k->#select(degrees E,d->d==degs_k)-1);
     deg:=0*n;
     box:={deg};
     scan(#n,k->
	  box=flatten apply(box,deg-> apply(n_k+1,i->deg+i*degs_k))
	  );
     box)	    
TEST ///
n={1,2,2}
(S,E)=setupRings(ZZ/101, n);
boxDegrees E
///


upperCorner=method()
upperCorner(ChainComplex,List) := (F,deg) ->(
     E:=ring F;
     degsE:= unique degrees E;
     n:=apply(degsE,dege->#select(degrees E,d->d==dege)-1);
     assert(
	  #(degrees E)_0 == # deg
	  );
     k:=sum deg;
     degsa:=degrees F_(k-1);
     L1:=select(#degsa,i->#select(#deg,j->degsa_i_j <= deg_j and degsa_i_j >= deg_j-n_j-1 )==#deg);
     degsb:=degrees F_k;
     L2:=select(#degsb,i->#select(#deg,j->degsb_i_j==deg_j)==#deg);
     ((F.dd_(k))^L1)_L2     
     )
TEST ///
n={1,2}
kk=ZZ/101
(S,E)=setupRings(kk,n)
F=dual res((ker transpose vars E)**E^{{ 2,3}},LengthLimit=>4)
betti F
tallyDegrees F
deg={ -2,-1} 
m=upperCorner(F,deg);
tally degrees source m, tally degrees target m
Fm=(res(coker m,LengthLimit=>4))[-sum deg-#n]
betti Fm
cohomologyTable(Fm,deg-{1,1},deg+{5,5})
///

lowerCorner=method()
lowerCorner(ChainComplex,List) := (F,deg) ->(
     E:=ring F;
     degsE:= unique degrees E;
     n:=apply(degsE,dege->#select(degrees E,d->d==dege)-1);
     assert(
	  #(degrees E)_0 == # deg
	  );
     k:= sum deg;
     degsa:=degrees F_k;     
     L1:=select(#degsa,i->#select(#deg,j->degsa_i_j==deg_j)==#deg);
     degsb:=degrees F_(k+1);
     L2:=select(#degsb,i->#select(#deg,j->degsb_i_j<=deg_j+n_j+1 and degsb_i_j>=deg_j)==#deg);
     ((F.dd_(k+1))^L1)_L2
     )
TEST///
n={1,1}
(S,E)=setupRings(ZZ/101,n)

time fB=dual res(coker random(E^7,E^{13:{ -1,0},11:{0,-1}}),LengthLimit=>10);	 	  
cohomologyTable(fB,{ -5,-5},{0,0})
deg={ -3,-3}
m= lowerCorner(fB,deg);
f= res( ker  m,LengthLimit=> 16)[-sum deg-2]
betti m, tally degrees target m, tally degrees source m
m= lowerCorner(fB,-6,{ -4,-4});
betti m, tally degrees target m, tally degrees source m
deg={3,3}
m=lowerCorner(f,deg);
betti m, tally degrees target m, tally degrees source m
m=lowerCorner(f,6,deg);
betti m, tally degrees target m, tally degrees source m
///
quadrant=method()
quadrant(ChainComplex,List,List) := (F,deg1,deg2) ->(
     m:=upperCorner(F,deg1);
     n:=apply(unique degrees ring F,e->#select(degrees ring F,d->d==e)-1);
     injRes:=dual res( coker transpose m,LengthLimit=> (sum deg1-sum deg2))[-sum deg1];
     cohomologyTable(injRes,deg2,deg1+n+{1,1})
     )



quadrant(ChainComplex,Sequence) := (F,L) ->(
     (k,deg1,deg2):=L;
     m:=upperCorner(F,k,deg1);
     n:=apply(unique degrees ring F,e->#select(degrees ring F,d->d==e)-1);
     injRes:=dual res( coker transpose m,LengthLimit=> (sum deg1-sum deg2))[-sum deg1];
     cohomologyTable(injRes,deg2,deg1+n+{1,1})
     )


tateCut=method()
tateCut(ChainComplex,Sequence):= (F,L)->(
     (k,deg1,deg2,deg3):=L;
     m:=lowerCorner(F,deg1);
     n:=apply(unique degrees ring F,e->#select(degrees ring F,d->d==e)-1);
     injRes:=dual res( coker transpose m,LengthLimit=> (sum deg1-sum deg2))[-sum deg1];
     projRes:=res(ker m,LengthLimit => (sum deg3-sum deg1))[-sum deg1-2];
     (cohomologyTable(injRes,deg2,deg1+n+{1,1}),cohomologyTable(projRes,deg1+((k-sum deg1)//2)*{1,1}-n,deg3),betti injRes,betti projRes)
     )     



symExt=method()
symExt(Matrix,Ring) := (m,E) -> (
     ev := map(E,ring m,vars E);
     mt := transpose jacobian m;
     jn := ev(syz mt);
     a:=(vars E**E^(rank target m));
--betti a,tally degrees source a, isHomogeneous a     
--betti jn, tally degrees target jn, tally degrees source jn, isHomogeneous jn
--tally(     degrees target jn+degrees source a)
     b:=a*jn;
--betti b, tally degrees target b, tally degrees source b
     c:=map(target b,E^(degrees source jn),b);
     transpose c)

    
TEST ///       
n={1,2}
(S,E)=setupRings(ZZ/101,n)
use S
vars S
m=map(S^4,S^{{ -1,0},{0,-1}}, transpose matrix{{S_0,S_1,0,0},{S_2,0,S_3,S_4}})
mE=symExt(m,E)

///
beilinsonWindow=method()
beilinsonWindow(ChainComplex) := C-> (
     E:= ring C;
     length C;
--    C':=C[min C];
--   T:=#unique degrees E;
     n:=apply(unique degrees E,d-> (#select( degrees  E, e-> e==d)-1)); -- dimension vector; could use numFactors
     Ck:=0;sourceCK:=0; targetCK:=0;d:=0;
     W:=chainComplex apply(min C+1..max C,k-> (Ck=(-1)^(min C)*C.dd_k;
	  --source indices and target rows and columns in the Beilison window
	  sourceCK = select(rank source Ck,i-> (d=degree (source Ck)_i;#select(#d,i->d_i>=0 and d_i<=n_i)==#n));	
          targetCK =  select(rank target Ck,i-> (d=degree (target Ck)_i;#select(#d,i->d_i>=0 and d_i<=n_i)==#n));	
     	  (Ck^targetCK)_sourceCK));
     --W':=W[-min C'];
     return W[-min C]) 

TEST ///
n={4}
(S,E)=setupRings(ZZ/101,n)
C=res ideal vars E
C1=C**E^{{ +1}}[0]
W=beilinsonWindow C1
apply(min W+1 ..max W,k->(W.dd_k==C1.dd_k,betti W.dd_k))
W
///

isChainComplex=method()
isChainComplex(ChainComplex) := W -> (
     lengthW:= max W- min W;
     #select(min W+1..max W-1,i->( if (source W.dd_i==0 or W.dd_(i+1)==0) then  true else W.dd_i*W.dd_(i+1)==0)) ==lengthW-1)


outsideBeilinsonRange=method()
outsideBeilinsonRange(Matrix) :=  m -> (
     E:= ring m;
     t:=#unique degrees E;
     n:=apply(unique degrees E,d-> (#select( degrees  E, e-> e==d)-1));
     d:=0;
	  --source indices not in the Beilison window
     sourcem := select(rank source m,i-> (d=degree (source m)_i;#select(#d,i->d_i<0 or d_i>n_i)>0));
     m_sourcem)	
          
outsideBeilinsonSource=method()
outsideBeilinsonSource(Matrix) :=  m -> (
     E:= ring m;
     t:=#unique degrees E;
     n:=apply(unique degrees E,d-> (#select( degrees  E, e-> e==d)-1));
	  --source indices not in the Beilison window
     d:=0;
     sourcem := select(rank source m,i-> (d:=degree (source m)_i;#select(#d,i->d_i<0 or d_i>n_i)>0));
     m_sourcem)	

outsideBeilinsonTarget=method()
outsideBeilinsonTarget(Matrix) :=  m -> (
     E:= ring m;
     t:=#unique degrees E;
     n:=apply(unique degrees E,d-> (#select( degrees  E, e-> e==d)-1));
	  --source indices not in the Beilison window
     d:=0;
     targetm := select(rank target m,i-> (d=degree (target m)_i;#select(#d,i->d_i<0 or d_i>n_i)>0));
     transpose (transpose m)_targetm)          
     

tateExtensionOld=method()
tateExtensionOld(ChainComplex) := C -> (
     W':=beilinsonWindow C;
--W', min W', max W' -- problem: min W'and max W' is not necesaarily correctly calculated
     minW':= min W';
     while W'_minW' ==0 do minW'=minW'+1;
     W:=W'[minW']; --now min W==0
--betti W, minW
     E:=ring C;
     t:=#unique degrees E;
     v:=apply(unique degrees E,d-> (#select( degrees  E, e-> e==d)));
     minW := 0;
     maxW := max W;
     while W_maxW == 0 do maxW=maxW-1;
--betti W, minW, maxW
     k:=minW;
     A:=W.dd_(k-1);
     B:=W.dd_(k); D:=0;
     L:=apply(toList(minW..maxW+1),k->(     
     	     D=W.dd_(k+1);
     	     A=B|mingens image(gens truncateInE(v,image syz A) %(gb image B)); 
     	     B=D||map(source outsideBeilinsonSource A,source D,0);
	     A));  -- the List of matrices from .dd_1 on   
     W1':=dual ((chainComplex L[1]));
--betti W, betti dual W1' -- both start in position the same position
     W1:=W1'[min W1'-2]; -- now min W1==2                                              min W1' shift
--betti W1
          --sloppy version: W a subcomplex is only true up to isomorphism
          --print betti(    dual (res(coker W1.dd_(3),LengthLimit=>length W'+2)[-min W1'+minW']));
--betti W'
--W1
     minW1:= min W1;    
     maxW1 := max W1;
--betti W1, min W1, max W1     
     k=minW1;
     A=W1.dd_(k+1);
     B=W1.dd_(k+2);D1:=0;
     L=apply(toList(minW1+2..maxW1),k->(     
     	     D=W1.dd_(k+1);D1=mingens image (syz A %gb image B); 
     	     A=B|D1; 
     	     B=D||map(source D1,source D,0);
	     A));
    W2':=dual( (chainComplex L)[-3]); -- the fact that we once started with an odd index the other time by an even index in the L's gives a sign difference, shifts cancel the sign
--betti W2', betti beilinsonWindow dual W1 -- beilinsonWindow dual W1 is a visible subquotient complex up, W2' should have one more term in front. 
--betti(b1=beilinsonWindow W2'),betti(b2=beilinsonWindow dual W1)
--b1,b2
--apply(min b1-1..max b1+3,k-> b1.dd_k==-b2.dd_k) --consistently opposit signs
    W2:= W2'[min W2'];
--betti W2
    maxW2:=max W2;minW2:= min W2;
    W3:=(chainComplex append(apply(toList(min W2+1..maxW2),k->W2.dd_k),syz W2.dd_maxW2))[-min W2' ];-- a sign changebetti W3,betti W2', betti dual W1 -- a visible complex only  different range
--betti(b1=beilinsonWindow W3),betti(b2=beilinsonWindow dual W1)
--b1,b2
--apply(min b1-1..max b1+3,k-> b1.dd_k==b2.dd_k) --consistently the same sign
--betti W3[min W1'-2], betti W
    W4:=W3[min W1'-2-minW' ];
--betti(b1=beilinsonWindow W4),betti(b2=beilinsonWindow C);
--b1,b2
--print apply(min b1-1..max b1+3,k-> b1.dd_k==b2.dd_k);
    return (W4)
     )

TEST ///
n={4}
(S,E)=setupRings(ZZ/101,n)
vars E
vars S
betti (fE=res ideal vars E)
A=transpose gens truncate({1},image transpose fE.dd_1)
C= chainComplex {A,fE.dd_2,map(source fE.dd_2,E^0,0) }**E^{{ -1}}[3]
isChainComplex C
betti C
W=beilinsonWindow C
betti W
T=tateExtensionOld(C)
betti T
WT=beilinsonWindow(T)
betti W, betti WT
WT, W
--the following two assertion should now hold
 apply(toList(min WT+1..max WT+1),k->WT.dd_k==W.dd_k)
betti res(coker T.dd_(min T+1),LengthLimit=> max T-min T)[-min T] == betti T
Cdual=dual C**E^{ -n}
--C=Cdual
betti beilinsonWindow Cdual, betti Cdual, betti W
Tdual=tateExtensionOld Cdual
betti Cdual
betti Tdual, betti T, betti (b1=beilinsonWindow Tdual), betti (b2=beilinsonWindow Cdual)
apply(min b1-1..max b1+3,k-> b1.dd_k==b2.dd_k) 
///

truncateInE=method()
truncateInE(List,Module):= (d,M) -> (
    base:=basis(M);
    degs:=degrees source base;
    m:=base_(select(#degs,k->#select(#d,i->degs_k_i >= d_i)==#d));
    image m 	
    )    

-- add
-- truncateInE(ZZ,Module):= (d,M) -> (
-- in case of ZZ-grading
  
truncateAboveBeilinsonWindow=method() -- maybe below is the better name ?
truncateAboveBeilinsonWindow(Matrix):= m -> (
    E:= ring m;
    t:=#unique degrees E;
    n:=apply(unique degrees E,d-> (#select( degrees  E, e-> e==d))-1);
    M:=image m;
    m1:=basis M;
    m1':= map(source m,source m1,m1);
    degs := degrees source m1';
    positionLists:=apply(t,k->select(#degs,j->degs_j_k>n_k));
    matrices:=apply(positionLists,L->m1'_L);
    m2:= matrices_0;
    scan(1..t-1,i-> m2=m2|matrices_i);
    mingens image(m*m2))

    

TEST ///
viewHelp
n={1,1}
(S,E)=setupRings(ZZ/101,n)
m=random(E^{2:{0,1}},E^3);
M= coker m
d={2,1}
truncateInE(d,M)==truncate(d,M)
///


tateExtension=method()
tateExtension(ChainComplex) := C -> (
         W':=beilinsonWindow C;
--W', min W', max W' -- problem: min W'and max W' is not necesaarily correctly calculated
     minW':= min W';
     while W'_minW' ==0 do minW'=minW'+1;
     W:=W'[minW']; --now min W==0
--betti W, minW
     E:=ring C;
     t:=#unique degrees E;
     v:=apply(unique degrees E,d-> (#select( degrees  E, e-> e==d)));
     minW := 0;
     maxW := max W;
     while W_maxW == 0 do maxW=maxW-1;
--betti W, minW, maxW
     k:=minW;
     A:=W.dd_(k-1);
     B:=W.dd_(k); D:=0;
     L:=apply(toList(minW..maxW+1),k->(     
     	     D=W.dd_(k+1);
     	     A=B|mingens image(gens truncateInE(v,image syz A) %(gb image B)); 
     	     B=D||map(source outsideBeilinsonSource A,source D,0);
	     A));  -- the List of matrices from .dd_1 on   
     W1':=dual ((chainComplex L[1]));
--betti W, betti dual W1' -- both start in position the same position
     W1:=W1'[min W1'-2]; -- now min W1==2                                              min W1' shift
--betti W1
          --sloppy version: W a subcomplex is only true up to isomorphism
          --print betti(    dual (res(coker W1.dd_(3),LengthLimit=>length W'+2)[-min W1'+minW']));
--betti W'
--W1
     minW1:= min W1;    
     maxW1 := max W1;
--betti W1, min W1, max W1     
     k=minW1;
     A=W1.dd_(k+1);
     B=W1.dd_(k+2);D1:=0;
     L=apply(toList(minW1+2..maxW1),k->(     
     	     D=W1.dd_(k+1);D1=mingens image (syz A %gb image B); 
     	     A=B|D1; 
     	     B=D||map(source D1,source D,0);
	     A));
    W2':=dual( (chainComplex L)[-3]); -- the fact that we once started with an odd index the other time by an even index in the L's gives a sign difference, shifts cancel the sign
--betti W2', betti beilinsonWindow dual W1 -- beilinsonWindow dual W1 is a visible subquotient complex up, W2' should have one more term in front. 
--betti(b1=beilinsonWindow W2'),betti(b2=beilinsonWindow dual W1)
--b1,b2
--apply(min b1-1..max b1+3,k-> b1.dd_k==-b2.dd_k) --consistently opposit signs
    W2:= W2'[min W2'];
--betti W2
    maxW2:=max W2;minW2:= min W2;
    W3:=(chainComplex append(apply(toList(min W2+1..maxW2),k->W2.dd_k),syz W2.dd_maxW2))[-min W2' ];-- a sign changebetti W3,betti W2', betti dual W1 -- a visible complex only  different range
--betti(b1=beilinsonWindow W3),betti(b2=beilinsonWindow dual W1)
--b1,b2
--apply(min b1-1..max b1+3,k-> b1.dd_k==b2.dd_k) --consistently the same sign
--betti W3[min W1'-2], betti W
    W4:=W3[min W1'-2-minW' ];
--betti(b1=beilinsonWindow W4),betti(b2=beilinsonWindow C);
--b1,b2
--print apply(min b1-1..max b1+3,k-> b1.dd_k==b2.dd_k);
    return (W4)
     )
    
    
--Start here for slow subcomplex
TEST ///
n={1,2}
(S,E)=setupRings(ZZ/101,n)
betti (fE= res(ideal vars E,LengthLimit=>8))

A=transpose mingens truncateInE({1,0},image transpose fE.dd_1)
C=chainComplex {A,fE.dd_2,fE.dd_3,map(source fE.dd_3,E^0,0) }**E^{{ 0,0}}[1]
isChainComplex C
betti C, betti dual C
Cdual= dual C ** E^{{0,0}-n}[-sum n+1]
betti Cdual
betti C
time W=beilinsonWindow C, betti W
T=tateExtension(W)
betti W, betti T ,betti (WT=beilinsonWindow T)
--the following two assertion should now hold
apply(toList(min WT..max WT+1),k->WT.dd_k==W.dd_k)
betti res(coker T.dd_(min T+1),LengthLimit=> max T-min T)[-min T] == betti T
---

m1=min T
betti T
betti (T1=(dual res(coker transpose  T.dd_(m1+1), LengthLimit=>12))[-m1-1])

tallyDegrees T1
qT=lastQuadrantComplex(T1,{ -6,-6});betti qT
betti (qT1=res(coker qT.dd_(-8),LengthLimit=>18)[8+1])
time bT= beilinsonWindow  qT1 --Mike this is slow
tallyDegrees bT

T2=dual qT1**E^{{ -1,0}}
W2=beilinsonWindow T2
tallyDegrees W2
T2=tateExtension W2
betti T2
m1=min T2
betti T2
betti (T21=(dual res(coker transpose  T2.dd_(m1+1), LengthLimit=>10))[-m1-1])

qT2=lastQuadrantComplex(T21,{ -8,-8});betti qT2
betti (qT21=res(coker qT2.dd_(-9),LengthLimit=>18)[5+1])
bT2= beilinsonWindow  qT21
tallyDegrees bT2
cohomologyTable(qT21,-{4,4},{4,4})
cohomologyTable(qT1,-{4,4},{4,4})




--dual approach
Cdual =dual C ** E^{{0,0}-n}[-sum n+2]
Wdual=beilinsonWindow Cdual
cohomologyTable(W,{0,0},n), cohomologyTable(Wdual,{0,0},n)
betti Cdual, betti beilinsonWindow Cdual
Tdual=tateExtension(Cdual)
betti Tdual
betti beilinsonWindow(Tdual)
tallyDegrees beilinsonWindow(Tdual)
betti Tdual
c=min Tdual
betti (Tdual1=(dual res(coker transpose Tdual.dd_(c+1),LengthLimit=>12))[-c-1])
qTdual=lastQuadrantComplex(Tdual1,{ -6,-6});betti qTdual
betti (qTdual1=res(coker qTdual.dd_(-8),LengthLimit=>18)[8+1])
bTdual= beilinsonWindow  qTdual1
tallyDegrees bTdual




///
-----------------
--  The corner Complex
-------------------------


quadrantMap=method()
quadrantMap(Matrix,List,List,List) := (M,c,I,J) -> (
    degSource:=degrees source M;
    degTarget:=degrees target M;
    t:= numFactors ring M; --the list {0,...,t-1}
    I':=select(t,j-> not member(j,I));
    J':=select(t,j-> not member(j,J));  
    goodColumns:=select(#degSource,k -> (#select(I,i-> degSource_k_i >= c_i)==#I and #select(I',i->degSource_k_i < c_i)==#I')); 
    goodRows:=select(#degTarget,k -> (#select(J,i-> degTarget_k_i >= c_i)==#J and #select(J',i->degTarget_k_i < c_i)==#J')); 	
    return ((M^goodRows)_goodColumns)) 


firstQuadrantComplex=method()
firstQuadrantComplex(ChainComplex,List) := (C,c) -> (
    s:=min C;
    C':=C[s];
    -- I:= numFactors ring C;
    Cge:=chainComplex apply(max C'-1,d -> quadrantMap(C'.dd_(d+1),c,{},{}));
    return Cge[-s])


lastQuadrantComplex=method()
lastQuadrantComplex(ChainComplex,List) := (C,c) -> (
    s:=min C;
    C':=C[s];
    I:= numFactors ring C;
    Cge:=chainComplex apply(max C'-1,d -> quadrantMap(C'.dd_(d+1),c,I,I));
    return Cge[-s])

cornerMap=method()
cornerMap(ChainComplex,List,ZZ) := (C,c,d) -> (
    E := ring C;
    t := numFactors E;
    Is:=reverse apply(t,i->select(t,j->j<i));
    M:= quadrantMap(C.dd_d,c,t,Is_0);
    Ms:=apply(#t-1,j->quadrantMap(C.dd_(d-j-1),c,Is_j,Is_(j+1))); 
    -- multiplication of empty matrices some times does not work! so there is a work around.
    scan(Ms, N-> if source N == E^0 then M=map(target N, source M,0) else M=N*M); 
    return M)

regionComplex=method()
regionComplex(ChainComplex,List,Sequence) := (T,c,IJK) -> (
    Ls:=apply(toList(min T..max T),k->goodColumns(T_k,c,IJK));
    rT:=chainComplex apply(max T- min T-1,k-> (-1)^(min T)*((T.dd_(min T + k +1))^(Ls_k))_(Ls_(k+1)));	
    rT[min T])

goodColumns=method()
goodColumns(Module,List,Sequence) := (F,c,IJK) -> (
    degF:=degrees F;
    select(#degF,g-> goodDegree(degF_g,c,IJK))
    )

goodDegree=method()
goodDegree(List,List,Sequence) := (d,c,IJK) -> (
    (I,J,K) := IJK;
    #select(I,i-> d_i<c_i)==#I and #select(J,j->d_j==c_j)==#J and #select(K,k->d_k>=c_k)==#K
    )	

TEST ///
        n={1,2,2}
        (S,E)=setupRings(ZZ/101,n)
        betti (fE= res(ideal vars E,LengthLimit=>8))
        T= fE[2]
        IJK=({0},{},{})
	c={2,3,2}
	rT=regionComplex(T,c,IJK)
	betti T, betti rT 
    	rT=regionComplex(T,c,({0},{1},{}))
    	tallyDegrees rT
       	rT=regionComplex(T,c,({0},{1},{2}))
    	tallyDegrees rT 
///

strand=method()
strand(ChainComplex,List,List) := (T,c,I) -> (
    regionComplex(T,c,({},I,{})))



cornerComplex=method()
cornerComplex(ChainComplex,List) := (C,c) -> (
    t:= numFactors ring C;
    if max C -min C < #t then error " need a complex of length at least t";
    C':= C[min C+1]; 
    Cge := firstQuadrantComplex(C'[-#t+1],c);
    Cle := lastQuadrantComplex(C',c);
--    <<(betti Cge, betti Cle) <<endl;
    A:=0;B:=0;AB:=0;d:=0;
    Ccorner:= chainComplex apply(max C- min C - #t,e-> (d:=e+#t-1; A=Cge.dd_(d);B= Cle.dd_(d); AB = cornerMap(C',c,d);
--	   print((betti A,betti AB,betti B));
	    (A|AB)||(map(target B, source A,0)|B)));
    return Ccorner[-min C-1])

TEST ///
n={1,2}
(S,E)=setupRings(ZZ/101,n)
betti (fE= res(ideal vars E,LengthLimit=>8))
T= fE[2]
betti T
betti (cT=cornerComplex(T,{2,3}))
tallyDegrees cT
tallyDegrees T 
tallyDegrees (cTextended=dual res(coker transpose  cT.dd_4,LengthLimit=>10))
betti cTextended
cohomologyTable(dual cTextended,{ -6,-6},{8,8})
///


beginDocumentation()

document { 
  Key => TateOnProducts,
  Headline => "Computation of parts of the Tate resoultion on products",
  "This package contains implementations of the algorithm from our paper ",
   HREF("http://arxiv.org/abs/","Tate Resolutions on Products of Projective Spaces") ,
  ". It allows to compute the direct image complexes of coherent sheaves on various factors.
     In the moment the function tateExtension is not completed, a version for single projective space
     however works nicely. The documentation and comments in the code are in a bad shape.
     Some function have to be removed, other wait for their implementations ",


   PARA{},
   SUBSECTION "From graded modules to Tate resolutions",  
   UL{   TO setupRings,
         TO numFactors,
	 TO symExt,
         TO upperCorner,
	 TO lowerCorner
      },
   SUBSECTION "Numerical Information",
   UL{ 
      TO cohomologyTable,
      TO tallyDegrees
     },
    SUBSECTION "Subcomplexes",
    UL{
       TO regionComplex,
       TO strand,
       TO firstQuadrantComplex,
       TO lastQuadrantComplex,
       TO cornerComplex,    
       TO isChainComplex
      },    
    SUBSECTION "Beilinson monads",
    UL{ 
	TO beilinsonWindow,
	TO tateExtensionOld,
--	TO truncateInE,
	TO tateExtension
        },
   }

doc ///
  Key
    setupRings
    (setupRings,Ring,List)
  Headline
    create the Cox ring and the corresponding exterior algebra for a product of projective spaces
  Usage
    (S,E)=setupRings(kk,n)
  Inputs
   kk: Ring
       the ground field, typically QQ or a finite field.
    n: List
       of dimensions of the factors
  Outputs
     : Sequence
       (S,E) of the Cox ring and the corresponding exterior algebra
  Description
     Text
     Example
        n={1,2}
        (S,E)=setupRings(ZZ/101,n)
	vars S
	vars E
///



doc ///
  Key
    numFactors
    (numFactors,Ring)
  Headline
    compute the number of factors of the Cox ring S or exterior algebra E
  Usage
     numFactors(S)
  Inputs
   S: Ring
      the Cow ring or the corresponding exterior algebra
  Outputs
     : List
       the List \{0,...,t-1\}, where t is the number of factors
  Description
     Text
     Example
        n={1,2,1,3}
        (S,E)=setupRings(ZZ/101,n)
	t=numFactors E
	#t==#n
///

doc ///
  Key
    symExt
    (symExt,Matrix,Ring)
  Headline
    from a linear presentation matrix over S to a linear presentation matrix over E and conversely 
  Usage
    symExt(m,E)
  Inputs
    m: Matrix
       a linear presentation matrix over S 
    E: Ring
       the Koszul dual ring of S
  Outputs
     : Matrix
       the corresponding linear presentation matrix over E
  Description
     Text
       Same method as in the single factor case
     Example
        n={1,2}
        (S,E)=setupRings(ZZ/101,n)
	vars S, vars E
        m=map(S^4,S^{{ -1,0},{0,-1}}, transpose matrix{{S_0,S_1,0,0},{S_2,0,S_3,S_4}})
        mE=symExt(m,E)

///

doc ///
  Key
    upperCorner
    (upperCorner,ChainComplex,List)
  Headline
    compute the upper corner 
  Usage
     m=upperCorner(F,d)
  Inputs
   F: ChainComplex
      over the exterior algebra
   d: List
      a (multi)-degree 
  Outputs
     : Matrix
       a submatrix of the differential $F_k -> F_{k-1}$
  Description
     Text
       Let $k = |d|$ be the total degree and $G \subset F_k$ the summand spanned by the generators of $F_k$ in degree d, 
       $H \subset F_{k-1}$ the summand spanned by generators of degree d' with $0 \le d'-d \le n$. The function returns
       the corresponding submatrix $m: G -> H$ of the differential.
     Example
        n={1,2}; kk=ZZ/101; (S,E)=setupRings(kk,n);
        F=dual res((ker transpose vars E)**E^{{ 2,3}},LengthLimit=>4)
        betti F
	tallyDegrees F
        deg={ -2,-1} 
        m=upperCorner(F,deg);
        tally degrees target m, tally degrees source m
        Fm=(res(coker m,LengthLimit=>4))[-sum deg-#n]
        betti Fm
        cohomologyTable(Fm,deg-{1,1},deg+{5,5})
///








doc ///
  Key
    lowerCorner
    (lowerCorner,ChainComplex,List)
  Headline
    compute the lower corner 
  Usage
     m=lowerCorner(F,d)
  Inputs
   F: ChainComplex
      over the exterior algebra
   d: List
      a (multi)-degree 
  Outputs
     : Matrix
       a submatrix of the differential $F_{k+1} -> F_{k}$
  Description
     Text
       Let $k = |deg|$ be the total degree and $G \subset F_k$ the summand spanned by the generators of $F_k$ in degree d, 
       $H \subset F_{k+1}$ the summand spanned by generators of degree d' with $0 \le d-d' \le n$. The function returns
       the corresponding submatrix $m: H -> G$ of the differential.
     Example
        n={1,2}; kk=ZZ/101; (S,E)=setupRings(kk,n);
        F=dual res((ker transpose vars E)**E^{{ 2,3}},LengthLimit=>4)
        betti F
	tallyDegrees F
        deg={ -2,-1} 
        m=lowerCorner(F,deg);
        tally degrees target m, tally degrees source m
        Fm=(res(coker m,LengthLimit=>7))[-sum deg]
        betti Fm
        cohomologyTable(Fm,deg,deg+{5,5})
///



doc ///
  Key
    cohomologyTable
    (cohomologyTable,ChainComplex,List,List)
  Headline
    cohomology Table
  Usage
    cohomologyTable(T,d1,d2)
  Inputs
    T: ChainComplex
       over the exterior algebra  
    d1: List
    d2: List
        two bi-degrees
  Outputs
     : Matrix
       incoding the dimension of the cohomology groups for twist $d$ with $d1 \le d \le d2$ 
  Description
     Text
        If $T$ is a (part) of the Tate resolution of a sheaf $F$ then the matrix of cohomology groups in the specified range will be returned
	encode by polynomials in ZZ[h]. The coefficient of $h^0,\ldots,h^N$ corresponds to $h^0(F(d),\ldots,h^N(F(d))$
     Example
        n={1,2}
        (S,E)=setupRings(ZZ/101,n)
        m=map(S^4,S^{{ -1,0},{0,-1}}, transpose matrix{{S_0,S_1,0,0},{S_2,0,S_3,S_4}})
        mE=symExt(m,E);
        fB=dual res(coker transpose mE,LengthLimit=>3)
	d={ -2,-1}
    	m=lowerCorner(fB, { -2,-1});
    	fm =res(coker m,LengthLimit=>10)[-sum d]
        cohomologyTable(fm,{ -2,-1},{4,5})       
    	tallyDegrees chainComplex lowerCorner(fm,{ -2,0})
	tallyDegrees chainComplex upperCorner(fm[-3],{ 2,4})
///
--      tallyDegrees fm
--	betti upperCorner(fm,{ -2,0})
--	tallyDegrees chainComplex lowerCorner(fm,{ -2,0})



doc ///
  Key
    tallyDegrees
    (tallyDegrees,ChainComplex)
  Headline
    number of generators in a multigarded complex
  Usage
    tallyDegrees T
  Inputs
    T: ChainComplex
       over the exterior algebra  
  Outputs
     : List
       a list of tally of the degree of the generators
  Description
     Text
     Example
        n={1,2}
        (S,E)=setupRings(ZZ/101,n)
        fm=dual res(coker vars E,LengthLimit=>3)
    	tallyDegrees fm
	fn=res(coker upperCorner(fm,{ -1,-1}),LengthLimit=>5)
	betti fn
        tallyDegrees fn
///




doc ///
  Key
    regionComplex
    (regionComplex,ChainComplex,List,Sequence)
  Headline
    region complex
  Usage
    regionComplex(T,c,IJK)
  Inputs
    T: ChainComplex
       over the exterior algebra
    c: List
       a (multi) degree 
    IJK: Sequence
       a sequence (I,J,K) of disjoint subsets of \{0..t-1\}       
  Outputs
     : ChainComplex
       a region complex of T
  Description
     Text
        We compute the region complex of T as defined in,
        HREF("http://arxiv.org/abs/","Tate Resolutions on Products of Projective Spaces"),
	Thm x.y 
     Example
        n={1,2,2}
        (S,E)=setupRings(ZZ/101,n)
        betti (fE= res(ideal vars E,LengthLimit=>8))
        T= fE[2]
        IJK=({0},{},{})
	c={2,3,2}
	rT=regionComplex(T,c,IJK)
	betti T, betti rT 
    	rT=regionComplex(T,c,({0},{1},{}))
    	tallyDegrees rT
       	rT=regionComplex(T,c,({0},{1},{2}))
    	tallyDegrees rT 
///


doc ///
  Key
    strand
    (strand,ChainComplex,List,List)
  Headline
    strand
  Usage
    strand(T,c,I)
  Inputs
    T: ChainComplex
       over the exterior algebra
    c: List
       a (multi) degree 
    I: List
       a sublist of \{0..t-1\} , where t denotes the number of factors      
  Outputs
     : ChainComplex
       a strand of T
  Description
     Text
        We compute the strand of T as defined in,
        HREF("http://arxiv.org/abs/","Tate Resolutions on Products of Projective Spaces"),
	Thm x.y. If T is (part of) the Tate resolution of a sheaf $F$, then the strand of $T$
	correponds to the Tate resolution $R\pi_*^{I'}(F(c))$ where $\pi^{I'}$ denotes the projection onto
	the factors corresponding to $I' =\{0,\ldots,t-1\} - I$ the complementary indices. 
     Example
        n={1,2}
        (S,E)=setupRings(ZZ/101,n)
        betti (T= res(ideal vars E,LengthLimit=>8))
	cohomologyTable(T,{0,0},{4,4})
	c={2,3}
	sT0=strand(T,c,{0}) 
    	sT1=strand(T,c,{1})
	sT01=strand(T,c,{0,1})  
      	tallyDegrees sT0
	tallyDegrees sT1
	
///


doc ///
  Key
     firstQuadrantComplex
     lastQuadrantComplex
    (firstQuadrantComplex,ChainComplex,List)
    (lastQuadrantComplex,ChainComplex,List)
  Headline
    first/last quadrant complex
  Usage
    lastQuadrantComplex(T,c)
    firstQuadrantComplex(T,c)
  Inputs
    T: ChainComplex
       over the exterior algebra
    c: List
       a (multi) degree     
  Outputs
     : ChainComplex
       a quadrant complex of T
  Description
      Text
         We compute the last quadrant complex of T as defined in 
         HREF("http://arxiv.org/abs/","Tate Resolutions on Products of Projective Spaces"),
         Thm x.y. To call, UL{ TO regionComplex}, would be an alternative. 
     Example
        n={1,2}
        (S,E)=setupRings(ZZ/101,n)
        betti (T= res(ideal vars E,LengthLimit=>8))
	cohomologyTable(T,{0,0},{4,5})
	c={2,3}
	lqT=lastQuadrantComplex(T,c)
	fqT=firstQuadrantComplex(T,c)
	tallyDegrees lqT
	tallyDegrees fqT
	cutT=res(coker fqT.dd_3,LengthLimit=>4)[-2] 
        cornerComplex(T,c)
///

doc ///
  Key
     cornerComplex
    (cornerComplex,ChainComplex,List)
  Headline
    corner complex
  Usage
    cornerComplex(T,c)
  Inputs
    T: ChainComplex
       over the exterior algebra
    c: List
       a (multi) degree     
  Outputs
     : ChainComplex
       the corner complex of T
  Description
      Text
         We compute the corner complex of T as defined in 
         HREF("http://arxiv.org/abs/","Tate Resolutions on Products of Projective Spaces"),
         Thm x.y. A corner complex is the mapping cone between a lastQuadrantComplex and a firstQuadrantComplex
     Example
        n={1,2}
        (S,E)=setupRings(ZZ/101,n)
	m=mingens (ideal vars E)^2
        betti (T1= (dual res( coker m,LengthLimit=>4))**E^{{0,-1}})
	n=lowerCorner(T1,-{3,1});
	T=res(coker n,LengthLimit=>9)[4]
        tallyDegrees T
	c={2,2}
	lqT=lastQuadrantComplex(T,c)
	fqT=firstQuadrantComplex(T,c)
        cT=cornerComplex(T,c)
	isChainComplex cT
	betti cT,betti fqT,betti lqT
	cohomologyTable(T,{ -3,-1},{4,5})	
///

doc ///
  Key
    isChainComplex
    (isChainComplex,ChainComplex)
  Headline
    is the input a chain complex?
  Usage
    isChainComplex C
  Inputs
    C: ChainComplex
  Outputs
     : Boolean 
  Description
     Text
        Checks whether in a chain complex the differentials really compose to zero
     Example
        C=chainComplex(matrix{{1}},matrix{{1}})
	isChainComplex C
	a=(matrix{{1,1}})
	b=matrix{{1},-{1}}
	C=chainComplex(a,b)
	isChainComplex C
///

doc ///
  Key
    beilinsonWindow
    (beilinsonWindow,ChainComplex)
  Headline
    Beilinson window
  Usage
    beilinsoWindow
  Inputs
    T: ChainComplex
  Outputs
     : ChainComplex
       a subquotient complex of T 
  Description
     Text
        The Beilinson monad ot a Tate resolution depends only on a finite free subquotient complex, which we call the Beilinson window.
     Example
        n={4}
        (S,E)=setupRings(ZZ/101,n)
        C=res ideal vars E
        C1=C**E^{{ +1}}[0]
        W=beilinsonWindow C1
        betti W, betti C	
///

doc ///
  Key
    tateExtensionOld
    (tateExtensionOld,ChainComplex)
  Headline
    Tate extension
  Usage
    eT=tateExtensionOld(T)
  Inputs
    T: ChainComplex
  Outputs
     : ChainComplex
       the Tate extension 
  Description
     Text
        The Beilinson monad representing an arbitrary object $\mathcal F \in D^b(\mathbb P)$ in the derived category can be 
	obtaine from a complex T with T=beilinsonWindow T.
	The function computes a larger complex eT which is finite piece of the Tate resolution, such that T = beilinsonWindow eT.
	So far this works only for a single projective space.
     Example
        n={4}
	(S,E)=setupRings(ZZ/101,n)
	betti (fE=res ideal vars E)
	A=transpose gens truncate({1},image transpose fE.dd_1)
	C= chainComplex {A,fE.dd_2,map(source fE.dd_2,E^0,0) }**E^{{ -1}}[3]
	isChainComplex C
	betti C
	W=beilinsonWindow C
	betti W
	T=tateExtensionOld(C)
	betti T
	WT=beilinsonWindow(T)
	betti W, betti WT
	--the following two assertion should now hold
 	apply(toList(min WT+1..max WT+1),k->WT.dd_k==W.dd_k)
        betti res(coker T.dd_(min T+1),LengthLimit=> max T-min T)[-min T] == betti T
        Cdual=dual C**E^{ -n}
	betti beilinsonWindow Cdual, betti Cdual, betti W
	Tdual=tateExtensionOld Cdual
	betti Cdual
        betti Tdual, betti (b1=beilinsonWindow Tdual), betti (b2=beilinsonWindow Cdual)
        apply(min b1-1..max b1+3,k-> b1.dd_k==b2.dd_k) 
	
///

doc ///
  Key
    tateExtension
    (tateExtension,ChainComplex)
  Headline
    Tate extension
  Usage
    T=tateExtension(C)
  Inputs
    C: ChainComplex
  Outputs
     : ChainComplex
       the Tate extension 
  Description
     Text
        The Beilinson monad representing an arbitrary object $\mathcal F \in D^b(\mathbb P)$ in the derived category can be 
	obtaine from a complex C with C=beilinsonWindow C.
	The function computes a larger complex T which is a finite piece of the Tate resolution, such that C = beilinsonWindow C.
	
	There is no correct implementation in the several factor case.
	
     Example
        n={1,2}
	(S,E)=setupRings(ZZ/101,n)
	betti (fE= res(ideal vars E,LengthLimit=>8))
	A=transpose mingens truncateInE({1,0},image transpose fE.dd_1)
	C=chainComplex {A,fE.dd_2,fE.dd_3,map(source fE.dd_3,E^0,0) }**E^{{ 0,1}}[1]
	isChainComplex C
	betti C, betti dual C
	W=beilinsonWindow C, betti W
	T=tateExtension(C)
	betti W, betti T ,betti (WT=beilinsonWindow T)
	--the following two assertion should now hold
	apply(toList(min WT..max WT+1),k->WT.dd_k==W.dd_k)
	betti res(coker T.dd_(min T+1),LengthLimit=> max T-min T)[-min T] == betti T
	--dual approach
	Cdual =dual C ** E^{{0,0}-n}[-sum n+2]
	Wdual=beilinsonWindow Cdual
	cohomologyTable(W,{0,0},n), cohomologyTable(Wdual,{0,0},n)
	betti Cdual, betti beilinsonWindow Cdual
	Tdual=tateExtension(Cdual)
	WTdual=beilinsonWindow(Tdual)
	apply(toList(min WTdual..max WTdual+1),k->WTdual.dd_k==Wdual.dd_k)
	betti WT, betti WTdual
	betti T, betti Tdual
	cohomologyTable(T,-n,2*n), cohomologyTable(Tdual,-n,2*n)

	
///






doc ///
  Key
    truncateInE
    (truncateInE,List,Module)
  Headline
    truncate in modules over an exterior algebra
  Usage
    truncateInE(d,M)
  Inputs
    d: List
       a (multi-) degree
  Outputs
     : Module
       the submodule of elemnts of degree $\ge d$
  Description
     Text
        a variant of the truncation command which might be needed 
     Example
        n={1,2}
        (S,E)=setupRings(ZZ/101,n)
        M=image matrix{{E_0+E_1,E_2+E_4}}
	d={2,1}
        betti truncate(d,M)
        betti truncateInE(d,M)	
///


end


--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end

restart
uninstallPackage"TateOnProducts"
installPackage"TateOnProducts"

loadPackage("TateOnProducts")
viewHelp TateOnProducts

 -- for talk in Hanoi
kk=ZZ/101
n={1,1}
(S,E) = setupRings(kk,n)
m=mingens (ideal vars E)^2
betti res coker transpose m
betti (fm=(res coker m)[-1])
tallyDegrees fm
n=lowerCorner(dual fm, -{3,3})
betti (fn=res(coker n,LengthLimit=> 10))
tex cohomologyTable(fn,-{3,3},{3,3})
betti dual fm[-1]
