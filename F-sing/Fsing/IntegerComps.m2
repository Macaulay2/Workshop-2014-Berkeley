--*************************************************
--*************************************************
--This file is used for doing computations with 
--integers that support other functions in the Fsing
--package.  
--*************************************************
--*************************************************

-- floorlog(b,x) computes floor(log_b x), correcting problems due to rounding
floorlog = ( b, x ) -> 
(
    flog := floor( log_b x ); -- first approximation (assumed to be <= correct value)
    while b^flog <= x do flog = flog + 1;
    flog - 1       
)

-- multOrder(a,b) finds the multiplicative order of a modulo b
multOrder = method()

multOrder( ZZ, ZZ ) := ( a, b ) ->
(
    if gcd( a, b ) != 1 then error "Expected numbers to be relatively prime.";
    n := 1;
    x := 1;
    while  (x = (x*a) % b) != 1  do n = n+1;
    n	      
)     


