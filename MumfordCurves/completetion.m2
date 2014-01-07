makeCompleteRing = (R,t,n)->(
     NR = new Ring from List;
     NR+NR := (a,b)->(
	  s := new MutableList from apply(n,i->a_i+b_i);
	  co := 0_R;
     	  for i from 0 to n-1 do (
		q := (s#i+co)%t;
	        co := (s#i+co-q)//t;
		s#i := q));
      	  new NR from s
     	  );
     NR
     );

PP = makeCompleteRing(ZZ,3,5);
a = new PP from {1,2,0,0,0};
b = new PP from {2,2,2,0,2};
a+b		