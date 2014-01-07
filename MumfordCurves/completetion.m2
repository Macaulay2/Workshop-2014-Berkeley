makeCompleteRing=(R,t,n)->(
     OR=R[x];
     NR=new PolynomialRing from OR;
     NR+NR:=(a,b)->(
	  g:=0_OR;
	  co:=0_OR;
     	  scan(n,i->(
		    ca:=coefficient(x^i,a);
		    cb:=coefficient(x^i,b);
		    s:=ca+cb+co;
		    q:=s%t;
		    co=(s-q)//t;
		    g=g+q*x^i));
	 sub(g,NR));
    NR
    )
