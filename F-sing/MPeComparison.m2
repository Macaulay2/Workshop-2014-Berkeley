functionA = (f, a, e) -> (
	fa := fastExp(f, a);
	varList := first entries vars ring f;
	m := ideal(varList);
	me := frobeniusPower(m, e);
	print "did powers";
	isSubset(ideal(fa), me)
)

functionB = (f, a, e) -> (
	R1 := ring f;
	varList := first entries vars R1;
	m := ideal(varList);
	my1 := sub(1, R1);
	root := ethRootSafe(f, ideal(my1), a, e); 
	
	isSubset(root, m)
)

functionC = (f, a, e) -> (
	R1 := ring f;
	varList := first entries vars R1;
	m := ideal(varList);
	my1 := sub(1, R1);
	root := ethRoot(ideal(f), a, e); 
	
	isSubset(root, m)
)
